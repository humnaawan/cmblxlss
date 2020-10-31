import numpy as np
import os
import healpy as hp
import time

from settings import *
# ------------------------------------------------------------------------------
from optparse import OptionParser
parser = OptionParser()
parser.add_option('--cmb_data_type',
                  dest='cmb_data_type',
                  help='Must be TT, mv, mvpol.')
parser.add_option('--grecon',
                  action='store_true', dest='grecon', default=False,
                  help='Using grecon.')
# ------------------------------------------------------------------------------
(options, args) = parser.parse_args()
print('\n## inputs: %s' % options)
# read in the inputs
grecon = options.grecon
cmb_data_type = options.cmb_data_type

outdir = '/global/cscratch1/sd/awan/soxlsst/output/'

cmb_data_tag = ''
if grecon:
    cmb_data_tag = 'grecon'
else: 
    cmb_data_tag = 'recon'
cmb_data_tag += '-%s' % cmb_data_type

outdir_data = '%s/interm-data_%s' % (outdir, cmb_data_tag)
if grecon:
    alms_files = [ f for f in os.listdir(outdir_data) if \
                  f.endswith('fits') and f.__contains__('grecon')]
else:
    alms_files = [ f for f in os.listdir(outdir_data) if \
                  f.endswith('fits') and not f.__contains__('grecon')]

kappa_files = [ f for f in alms_files if f.__contains__('kappa') ]

nside = 256
lmax = 4000

for i, kappa_file in enumerate( kappa_files ):
    fname_tag = kappa_file.split('_')[-1].split('.npz')[0]
    
    print('## reading in %sth CMB alms ...' % (i+1))
    # okay now read in the CMB alms
    kappa_alms = hp.read_alm('%s/%s' % (outdir_data, kappa_file))

    print('## reading in %sth galaxy alms ...' % (i+1))
    # okay now read in the galaxy alms
    galaxy_file = [ f for f in alms_files if f.__contains__('correlated-g') \
                   and f.__contains__(fname_tag) ][0]
    galaxy_alms = hp.read_alm('%s/%s' % (outdir_data, galaxy_file))
    
    # spectra
    cls_kk = hp.alm2cl(kappa_alms)[0:lmax]
    cls_gg = hp.alm2cl(galaxy_alms)[0:lmax]
    cls_kg = hp.alm2cl(kappa_alms, galaxy_alms)[0:lmax]

    if i == 0:
        cls_kk_stack = [ cls_kk.copy() ]
        cls_kg_stack = [ cls_kg.copy() ]
        cls_gg_stack = [ cls_gg.copy() ] 
    else:
        cls_kk_stack = np.vstack( [cls_kk_stack, cls_kk.copy()] )
        cls_kg_stack = np.vstack( [cls_kg_stack, cls_kg.copy()] )
        cls_gg_stack = np.vstack( [cls_gg_stack, cls_gg.copy()] )
        
# now plot
ells = np.arange(lmax)
# calculate the mean and stddev
def deal_with_stack(stack):
    mean = np.mean(stack, axis=0)
    stddev = np.std(stack, axis=0)
    y = mean #(mean * ells * (ells+1)) / 2.0 / np.pi
    yerr = stddev #(stddev * ells * (ells+1)) / 2.0 / np.pi
    return y, yerr

c_ells, c_ells_err = {}, {}
key = r'$\kappa\kappa$'
c_ells[key], c_ells_err[key] = deal_with_stack(stack=cls_kk_stack)

key = r'$\kappa g$'
c_ells[key], c_ells_err[key] = deal_with_stack(stack=cls_kg_stack)

key = r'$gg$'
c_ells[key], c_ells_err[key] = deal_with_stack(stack=cls_gg_stack)

# plot things out
alpha = 0.25

plt.clf()
for key in c_ells:
    p = plt.plot(ells, c_ells[key], '.-', linewidth=2, label=key)
    plt.fill_between(ells, c_ells[key], c_ells[key] + c_ells_err[key],
                     alpha=alpha, facecolor=p[0].get_color())
    plt.fill_between(ells, c_ells[key], c_ells[key] - c_ells_err[key],
                     alpha=alpha, facecolor=p[0].get_color())

ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('$\ell$')
ax.set_ylabel('$C_\ell$')
ax.set_ylim([1e-13, 1e-5])
#plt.ylabel(r'$\ell(\ell+1)C_\ell/(2\pi)$')
ax.legend(bbox_to_anchor=(1,1))
fname = 'plot_spectra_%s-sims_%s.png' % (len(kappa_files), cmb_data_tag)
plt.savefig('%s/%s' % (outdir, fname), format='png', bbox_inches='tight')
plt.close()
print('## saved spectra plots in %s' % fname)

fname = 'data_cells-mean-stddev_%s-sims_%s.npz' % (len(kappa_files), cmb_data_tag)
np.savez('%s/%s' % (outdir, fname), cells=c_ells, c_ells_err=c_ells_err)
print('## saved data in %s' % fname)
