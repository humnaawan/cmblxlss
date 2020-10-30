import numpy as np
import os
import healpy as hp
import pyccl as ccl
import time
from pixell import reproject,enmap

from helpers import get_dndz, generate_correlated_alm
from utils_plot import plot_power_spectrum, plot_skymap
from add_fakelss import add_fakelss
from settings import *

outdir = '/global/cscratch1/sd/awan/soxlsst/output/'

one_sim_only = False
grecon = True
cmb_data_type = 'TT' # must be one of the following: 'TT', 'mv', 'mvpol'
cmb_data_path = '/global/cscratch1/sd/msyriac/data/depot/solsst/'
lensing_mask_path = '%s/v1_mask.fits' % cmb_data_path

fakelss_map_path = '/global/homes/a/awan/so/cmblxlss/data/fakelss_baseline-v1.4_nside1024.npz'

nside = 1024
lmax = 4000 #3 * nside

# set up for CCL
cparams = {'Omega_c' : 0.27, 'Omega_b' : 0.045, 'Omega_k' : 0.0, 'Omega_nu' :0.0,
           'h0' : 0.67, 'sigma_8' : 0.8, 'n_s' : 0.96, 'w': -1, 'wa': 0 ,
           'transfer_function' : 'eisenstein_hu', 'matter_power_spectrum' : 'linear',
           'has_rsd' : False, 'has_magnification' : False}

zedges = np.array( [ (0.1, 1.2) ] ) #[(0.6, 0.8), (0.8, 1.0), (1.0, 1.2)] )
nzbins = len(zedges)

# LSST-related params
bias_func = lambda z: 1 + 0.84 * z
zmin, zmax = 1e-2, 3 # range over which the dndz, etc will be considered
ilim_ext = 25.3   # extended source limiting magnitude

start_time = time.time()

print('## setting up things for CCL ...')
# set up cosmology for CCL
cosmo = ccl.Cosmology(Omega_c=cparams['Omega_c'],
                          Omega_b=cparams['Omega_b'],
                          h=cparams['h0'],
                          sigma8=cparams['sigma_8'],
                          n_s=cparams['n_s'],
                          w0=cparams['w'],
                          wa=cparams.get('wa', 0.0),
                          transfer_function=cparams['transfer_function'],
                          matter_power_spectrum=cparams['matter_power_spectrum']
                         )

z_arr, dndz = get_dndz(outdir=outdir, zedges=zedges, zmin_overall=zmin,
                       zmax_overall=zmax, ilim=ilim_ext)
bias = bias_func(z_arr)

print('## getting gg spectra ...')
# get the gg cells
gal_counts = ([ccl.NumberCountsTracer(cosmo, has_rsd=False,
                                      dndz=(z_arr, dndz[zi]),
                                      bias=(z_arr, bias)
                                     )
                   for zi in range(0, nzbins)])
n_tracer = len(gal_counts)

lmin_, lmax_ = 0, lmax
ells = np.arange(lmin_, lmax_, 1.0)
c_gg = ([[ccl.angular_cl(cosmo, gal_counts[ni], gal_counts[nj], ells)
          for ni in range(0, n_tracer)] for nj in range(0, n_tracer)])

print('## getting kk spectra ...')
# get the kk cells
lensing = ccl.CMBLensingTracer(cosmo=cosmo, z_source=1100)
c_kk = ccl.angular_cl(cosmo, lensing, lensing, ells)

print('## getting kg spectra ...')
# get the kg cells
c_kg = ([ccl.angular_cl(cosmo, gal_counts[ni], lensing, ells)
            for ni in range(0, n_tracer)])

# plot things out
plt.clf()
for i in range(nzbins):
    for j in range(i, nzbins):
        plt.plot(ells, c_gg[i][j], '+:', linewidth=2, label='gg %s-%s' % (i, j))

for i in range(nzbins):
    plt.plot(ells, c_kg[i], 'x--', linewidth=2, label=r'g%s - $\kappa$' % i)

plt.plot(ells, c_kk, '.-', linewidth=2, label=r'$\kappa\kappa$')
    
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('$\ell$')
ax.set_ylabel('$C_\ell$')
ax.legend(bbox_to_anchor=(1,1))
fname = 'plot_spectra-ccl_%sbins.png' % nzbins
plt.savefig('%s/%s' % (outdir, fname), format='png', bbox_inches='tight')
plt.close()
print('## saved spectra plots in %s' % fname)

print('## reading in the CMB alms ...')

cmb_data_tag = ''
if grecon:
    fnames_alms = [ f for f in os.listdir(cmb_data_path) if \
                   f.__contains__('grecon') ]
    cmb_data_tag = 'grecon'
else:
    fnames_alms = [ f for f in os.listdir(cmb_data_path) if \
                  not f.__contains__('grecon') ]
    cmb_data_tag = 'recon'
# now select the specific kind of data
fnames_alms = [ f for f in fnames_alms if \
               f.__contains__(cmb_data_type) ][0:100]
cmb_data_tag += '-%s' % cmb_data_type

if one_sim_only:
    fnames_alms = fnames_alms[0:1]

lsst_mask = None
lensing_mask = None
joint_mask = None

outdir_plots = '%s/interm-plots_%s' % (outdir, cmb_data_tag)
os.makedirs(outdir_plots, exist_ok=True)

outdir_data = '%s/interm-data_%s' % (outdir, cmb_data_tag)
os.makedirs(outdir_data, exist_ok=True)

print('## %s sims to consider.' % len(fnames_alms))
for i, fname_alms in enumerate( fnames_alms ):
    start_time0 = time.time()
    fname_tag = fname_alms.split('_')[-1].split('.fits')[0]
    print('\n## running things for sim %s: %sth' % (fname_tag, i+1))
    # okay now read in the CMB alms
    kappa_alms = hp.read_alm('%s/%s' % (cmb_data_path, fname_alms))
    kappa_alms[~np.isfinite(kappa_alms)] = 0
    for nzbin in range(nzbins):
        print('## generating correlated alms ...')
        # generate galaxy alm that are correlated
        # try for just one bin
        gal_density_alm = generate_correlated_alm(input_alm_f1=kappa_alms,
                                                  Clf1f1=c_gg[nzbin][nzbin],
                                                  Clf2f2=c_kk,
                                                  Clf1f2=c_kg[nzbin], seed=nzbin)

        # construct the map from the alms
        gal_density_map = hp.alm2map(gal_density_alm, nside=nside)

        if i == 0: save_plot = True
        else: save_plot = False

        plot_skymap(map_in=gal_density_map, title='theory, correlated g',
                    data_label='density', show_plot=False, save_plot=save_plot,
                    outdir=outdir_plots, file_tag='')

        print('## adding fakelss ...')
        # add fake lss
        gal_density_wfakelss_map = add_fakelss(theory_density_map=gal_density_map,
                                               fakelss_map_path=fakelss_map_path,
                                               show_plot=False, save_plot=save_plot,
                                               outdir=outdir_plots)

        if joint_mask is None:
            # set up lsst mask
            if lsst_mask is None:
                lsst_mask = 1 - gal_density_wfakelss_map.mask.astype(int)
            # set up the lensing mask
            if lensing_mask is None:
                lensing_mask = enmap.read_map(lensing_mask_path)
                lensing_mask = reproject.healpix_from_enmap(lensing_mask, lmax=6000, nside=nside)
                lensing_mask[lensing_mask<0.99] = 0
            # set up the joint mask
            joint_mask = lensing_mask * lsst_mask
            
        # update the mask
        gal_density_wfakelss_map.mask = joint_mask.copy()
        plot_skymap(map_in=gal_density_wfakelss_map,
                    title='correlated g + fakelss; joint mask',
                    data_label='density', show_plot=False, save_plot=save_plot,
                    outdir=outdir_plots, file_tag='')
        
        gal_density_wfakelss_alms = hp.map2alm(gal_density_wfakelss_map, lmax=lmax)
        
        fname = 'alms_correlated-g-wfakelss-jmasked_%s-%s.fits' % (cmb_data_tag, fname_tag)
        hp.fitsfunc.write_alm(filename='%s/%s' % (outdir_data, fname),
                              alms=gal_density_wfakelss_alms, overwrite=True)
        print('## saved alms in %s' % fname)
        
    # mask the kappa map too    
    kappa_map = hp.alm2map(kappa_alms, nside=nside) * joint_mask
    plot_skymap(map_in=kappa_map, title='theory, kappa', data_label='kappa',
                show_plot=False, save_plot=save_plot,
                outdir=outdir_plots, file_tag='')
    kappa_alms = hp.map2alm(kappa_map, lmax=lmax)
    fname = 'alms_kappa-jmasked_%s-%s.fits' % (cmb_data_tag, fname_tag)
    hp.fitsfunc.write_alm(filename='%s/%s' % (outdir_data, fname), alms=kappa_alms,
                          overwrite=True)
    print('## saved alms in %s' % fname)
    print('## time taken: %.2f min' % ((time.time() - start_time0)/60.))

print('## time taken: %.2f min' % ((time.time() - start_time)/60.))