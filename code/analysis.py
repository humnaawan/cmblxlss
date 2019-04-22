import matplotlib
matplotlib.use('Agg')
import numpy as np
import healpy as hp
import pandas as pd
import os
from collections import OrderedDict
import time
from optparse import OptionParser
from orphics import cosmology
# ------------------------------------------------------------------------------
# imports from this folder
from utils_plot import plot_cls_dict, plot_mollview
from utils_analysis import *
#-------------------------------------------------------------------------------
parser = OptionParser()
# required params
parser.add_option('--outdir', dest='main_outdir',
                  help='Output directory')
parser.add_option('--lensed-cmb-path', dest='lensed_cmb_map_path',
                  help='Path to the file containing lensed cmb map')
parser.add_option('--cosmology-path', dest='cosmo_path',
                  help='Path to orphics cosmology')
parser.add_option('--kappa-norm-path', dest='kappa_norm_path',
                  help='Path to the file containing the normalization of the kappa alms')
parser.add_option('--kappa-alm-theory_path', dest='kappa_alm_theory_path',
                  help='Path to the file containing the theory kappa alms')
parser.add_option('--gcals-path', dest='gcls_path',
                  help='Path to the file containing the theory gk and kk cls')
parser.add_option('--fg-map-path', dest='fg_map_path',
                  help='Path to the file containing the foregrounds map')
parser.add_option('--lsst-path', dest='lsst_path',
                  help='Path to the npz file containing the lsst opsim output')
parser.add_option('--lsstdata-tag', dest='lsst_data_tag',
                  help='Tag to identify lsst data')
parser.add_option('--lmax', dest='lmax', type='int',
                  help='Maximum multipole to consider')
# optional
parser.add_option('--completeness-thresh', dest='completeness_threshold',
                  help='Threshold on deltaNByN to define the binary (0/1) completeness map',
                  default=-0.2)
parser.add_option('--smoothing-fwhm', dest='smoothing_fwhm',
                  help='FWHM (in radians) of the gaussian beam used to smooth the halo density map',
                  default=np.radians(1.))
#--------------------------------------------------------------------------------------
# parse the options
(options, args) = parser.parse_args()
main_outdir = options.main_outdir
lensed_cmb_map_path = options.lensed_cmb_map_path
cosmo_path = options.cosmo_path
kappa_norm_path = options.kappa_norm_path
kappa_alm_theory_path = options.kappa_alm_theory_path
gcls_path = options.gcls_path
fg_map_path = options.fg_map_path
lsst_path = options.lsst_path
lsst_data_tag = options.lsst_data_tag
lmax = options.lmax
completeness_threshold = options.completeness_threshold
smoothing_fwhm = options.smoothing_fwhm

# ------------------------------------------------------------------------------
time0 = time.time()
# set up the readme
readme = ''
readme_tag = 'lmax%s'%(lmax)
# set up the outdir
outdir = '%s/output_%s_lsst-%s'%(main_outdir, readme_tag, lsst_data_tag)
if not os.path.exists(outdir): os.makedirs(outdir)
# update readme
readme = print_update(update='\nOutput directory: %s\n'%outdir.split(main_outdir)[1],
                      readme=readme)
# set up dictionary that will contain all the spectra we care about
c_ells = OrderedDict()
# set up the ells
ells = np.arange(0, lmax)
mlmax = lmax + 1000
# -----------------------------------------------
# read in the lense cmb map
readme = print_update(update='\n## Reading in lensed cmb map ... \n',
                      readme=readme)
lensed_cmb_map = hp.read_map(lensed_cmb_map_path)
# get the nside
nside_target = hp.npix2nside(len(lensed_cmb_map))
# plot
filename = plot_mollview(map_in=lensed_cmb_map,
                         title='lensed_cmb_map',
                         data_label='',
                         outdir=outdir,
                         file_tag='lensed_cmb_map',
                         save_plot=True, show_plot=False)
readme = print_update(update='Saved lensed cmb map in %s\n'%(filename),
                      readme=readme)
# -----------------------------------------------
# read in the fg map
readme = print_update(update='\n## Reading in foregrounds map ... \n',
                      readme=readme)
fg_map = hp.read_map(fg_map_path)
# save the fg map
filename = plot_mollview(map_in=fg_map,
                         title='fg_map',
                         data_label='',
                         outdir=outdir,
                         file_tag='fg_map',
                         save_plot=True, show_plot=False)
readme = print_update(update='Saved fg map in %s\n'%(filename),
                    readme=readme)
# rotate the map
readme = print_update(update='Rotating the map from galactic to equatorial coords ... ',
                      readme=readme)
r = hp.rotator.Rotator(coord=['G','C']) # rotator object: G=galactic --> C=equatorial
fg_map = rotate_map(r, fg_map)  # rotate the map
# save the fg map
filename = plot_mollview(map_in=fg_map,
                         title='fg_map_rotated',
                         data_label='',
                         outdir=outdir,
                         file_tag='fg_map_rotated',
                         save_plot=True, show_plot=False)
readme = print_update(update='Saved rotated fg map in %s\n'%(filename),
                      readme=readme)
# -----------------------------------------------
# get some theory stuff from orphics
theory = cosmology.loadTheorySpectraFromCAMB(cosmo_path,
                                             get_dimensionless=False)
cl_tt_theory = theory.lCl(XYType='TT', ell=ells)
c_ells[r'$\kappa\kappa$ theory'] = theory.gCl('kk', ells)

# -----------------------------------------------
# read in theory spectra for gg and kg
gcls = np.genfromtxt(gcls_path)
#ells_gcls = gcls[:, 0]
c_ells[r'$gg$ theory'] = gcls[:, 1][0: lmax]
c_ells[r'$\kappa g$ theory'] = gcls[:, 2][0: lmax]
# -----------------------------------------------
# read in theory alms for phi (related to kappa)
readme = print_update(update='\n## Reading in theory phi alms ... \n',
                      readme=readme)
kappa_theory_phi_alms = hp.read_alm(kappa_alm_theory_path)
kappa_theory_phi_alms = kappa_theory_phi_alms.astype(np.complex128)
kappa_theory_map = hp.alm2map(kappa_theory_phi_alms, nside=nside_target)
kappa_theory_phi_alms_red = hp.map2alm(kappa_theory_map, lmax=mlmax)
kappa_theory_filt = hp.almxfl(alm=kappa_theory_phi_alms_red, fl=ells*(ells+1)/2)
# -----------------------------------------------
# read in the normalization for reconstructed alms
readme = print_update(update='\b## Reading in kappa normalization ... \b',
                      readme=readme)
kappa_norm = np.genfromtxt(kappa_norm_path)[:, 1]
# -----------------------------------------------
# get the lsst dn/n map and the completeness mask (binary; apodized)
out = get_lsst_maps(data_file=lsst_path, data_tag=lsst_data_tag,
                    data_label='dNbyN', nside_out=nside_target,
                    completeness_threshold=completeness_threshold,
                    smoothing_fwhm=smoothing_fwhm, outdir=outdir,
                    plot_interm=True)
lsst_mask_smoothed, lsst_data_map, lsst_fsky, readme = out
# -----------------------------------------------
# construct correlated galaxy density alms
readme = print_update(update='\n## Generating correlated density field ... \n',
                      readme=readme)
gal_density_alm = generate_correlated_alm(input_alm_f1=kappa_theory_filt,
                                          Clf1f1=c_ells[r'$\kappa\kappa$ theory'],
                                          Clf2f2=c_ells[r'$gg$ theory'],
                                          Clf1f2=c_ells[r'$\kappa g$ theory'],
                                          seed=1)
# construct the map from the alms
gal_density_map = hp.alm2map(gal_density_alm, nside=nside_target)
# plot
filename = plot_mollview(map_in=gal_density_map,
                         title='correlated galaxy density',
                         data_label='',
                         outdir=outdir,
                         file_tag='galdensity',
                         save_plot=True, show_plot=False)
readme = print_update(update='Saved the correlated galaxy density map in %s\n'%(filename),
                      readme=readme)
# modulate the galaxy density map with fake LSS
gal_density_map_mod = (gal_density_map + 1) * (lsst_data_map + 1) - 1
# now multiply by the apodized mask
gal_density_map_mod *= lsst_mask_smoothed
# plot
filename = plot_mollview(map_in=gal_density_map_mod,
                         title='modulated correlated galaxy density',
                         data_label='',
                         outdir=outdir,
                         file_tag='modulated-galdensity',
                         save_plot=True, show_plot=False)
readme = print_update(update='Saved the moduldated correlated galaxy density map in %s\n'%(filename),
                      readme=readme)
gal_density_alm_mod = hp.map2alm(gal_density_map_mod, lmax=mlmax)
# -----------------------------------------------
# set up the filter for the cmb reconstruction
kappa_filter = 1/cl_tt_theory
kappa_filter[ells < 300] = 0
# -----------------------------------------------
# baseline reconstruction
readme = print_update(update='\n## Reconstructing kappa (no fg or mask) ... \n',
                      readme=readme)
out = get_reconstructed_kappa_alm(lensed_cmb_map=lensed_cmb_map,
                                  kappa_filter=kappa_filter, kappa_norm=kappa_norm,
                                  lmax=lmax, mlmax=lmax+1000, outdir=outdir,
                                  lsst_mask_map=None, fg_map=None)
kappa_alms_normed, readme = out
# -----------------------------------------------
# reconstruction with foregrounds
readme = print_update(update='\n## Reconstructing kappa with fg ... \n',
                      readme=readme)
out = get_reconstructed_kappa_alm(lensed_cmb_map=lensed_cmb_map,
                                  kappa_filter=kappa_filter, kappa_norm=kappa_norm,
                                  lmax=lmax, mlmax=lmax+1000, outdir=outdir,
                                  lsst_mask_map=None, fg_map=fg_map)
kappa_alms_normed_fg, readme = out
# -----------------------------------------------
# reconstruction with foregrounds + lsst mask
readme = print_update(update='\n## Reconstructing kappa with fg & lsst mask ... \n',
                      readme=readme)
out = get_reconstructed_kappa_alm(lensed_cmb_map=lensed_cmb_map,
                                  kappa_filter=kappa_filter, kappa_norm=kappa_norm,
                                  lmax=lmax, mlmax=mlmax, outdir=outdir,
                                  lsst_mask_map=lsst_mask_smoothed,
                                  fg_map=fg_map)
kappa_alms_normed_fg_masked, readme = out
# -----------------------------------------------
# calculate correlations
readme = print_update(update='\n## Calculating correlations ... \n',
                      readme=readme)
# autos
c_ells[r'$\kappa\kappa$ baseline'] = hp.alm2cl(kappa_alms_normed, lmax=lmax)[0:lmax]
c_ells[r'$\kappa\kappa$ w/ fg'] = hp.alm2cl(kappa_alms_normed_fg, lmax=lmax)[0:lmax]
c_ells[r'$\kappa\kappa$ w/ fg + lsst mask'] = hp.alm2cl(kappa_alms_normed_fg_masked, lmax=lmax)[0:lmax]/lsst_fsky
# cross correlation
c_ells[r'$\kappa$ baseline x correlated $g$'] = hp.alm2cl(kappa_alms_normed, gal_density_alm, lmax=lmax)[0:lmax]
c_ells[r'$\kappa$ baseline x lsst modulated correlated $g$'] = hp.alm2cl(kappa_alms_normed, gal_density_alm_mod, lmax=lmax)[0:lmax]
c_ells[r'$\kappa$ w/ fg x correlated $g$'] = hp.alm2cl(kappa_alms_normed_fg, gal_density_alm, lmax=lmax)[0:lmax]
c_ells[r'$\kappa$ w/ fg + lsst mask x correlated $g$'] = hp.alm2cl(kappa_alms_normed_fg_masked, gal_density_alm, lmax=lmax)[0:lmax]/lsst_fsky
c_ells[r'$\kappa$ w/ fg + lsst mask x lsst modulated correlated $g$'] = hp.alm2cl(kappa_alms_normed_fg_masked, gal_density_alm_mod, lmax=lmax)[0:lmax]

# -----------------------------------------------
# plot the cross-correlations in sets
readme = print_update(update='\n## Plotting spectra ... \n',
                      readme=readme)
# first plot the auto spectra for checks
c_ells_to_plot = {}
for key in [f for f in c_ells.keys() if f.__contains__(r'\kappa\kappa')]:
    c_ells_to_plot[key] = c_ells[key]
# plot
filename = plot_cls_dict(cls_in=c_ells_to_plot, outdir=outdir, file_tag='kk-only',
                         save_plot=True, show_plot=False,
                         cross_convention=True,
                         sci_yticks=True,
                         binned=False, bin_width=20, lmax=lmax)
readme = print_update(update='Saved %s'%filename,
                      readme=readme)
# plot binned
filename = plot_cls_dict(cls_in=c_ells_to_plot, outdir=outdir, file_tag='kk-only',
                         save_plot=True, show_plot=False,
                         cross_convention=True,
                         sci_yticks=True,
                         binned=True, bin_width=50, lmax=lmax)
readme = print_update(update='Saved %s'%filename,
                      readme=readme)
# -----------------------------------------------
# now plot the cross correlations
c_ells_to_plot = {}
c_ells_to_plot[r'$\kappa g$ theory'] = c_ells[r'$\kappa g$ theory']
for key in [f for f in c_ells.keys() if f.__contains__(r'correlated $g$')]:
    c_ells_to_plot[key] = c_ells[key]
# plot
filename = plot_cls_dict(cls_in=c_ells_to_plot, outdir=outdir, file_tag='kg-only',
                         save_plot=True, show_plot=False,
                         cross_convention=True,
                         sci_yticks=True,
                         binned=False, bin_width=20, lmax=lmax)
readme = print_update(update='Saved %s'%filename,
                      readme=readme)
# binned
filename = plot_cls_dict(cls_in=c_ells_to_plot, outdir=outdir, file_tag='kg-only',
                         save_plot=True, show_plot=False,
                         cross_convention=True,
                         sci_yticks=True,
                         binned=True, bin_width=50, lmax=lmax)
readme = print_update(update='Saved %s'%filename,
                      readme=readme)
# -----------------------------------------------
# save the corrs
filename = 'cross_specs.csv'
pd.DataFrame(c_ells).to_csv('%s/%s'%(outdir, filename), index=False)
# update readme
readme = print_update(update='Saved cross-spectra in %s.\n' % filename,
                      readme=readme)
readme = print_update(update='\nTime taken for the whole thing: %.3f min\n' % ((time.time()-time0)/60.),
                      readme=readme)

readme_file = open('%s/readme_%s.txt' % (outdir, readme_tag), 'a')
readme_file.write(readme)
readme_file.close()
