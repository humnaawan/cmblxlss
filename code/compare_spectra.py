import numpy as np
import os
import healpy as hp
import time
import pickle

from settings import *
# ------------------------------------------------------------------------------
from optparse import OptionParser
parser = OptionParser()
parser.add_option('--cmb_data_type',
                  dest='cmb_data_type',
                  help='Must be TT, mv, mvpol.')
parser.add_option('--nsims',
                  dest='nsims',
                  help='Number of sims considered.', default=100)
# ------------------------------------------------------------------------------
(options, args) = parser.parse_args()
print('\n## inputs: %s' % options)
# read in the inputs
cmb_data_type = options.cmb_data_type
nsims = options.nsims

# outdir
outdir = '/global/cscratch1/sd/awan/soxlsst/output/'

fname = 'data_cells-mean-stddev_%s-sims_grecon-%s.pickle' % (nsims, cmb_data_type)
with open('%s/%s'%(outdir, fname), 'rb') as handle:
    with_fgs = pickle.load(file=handle)

fname = 'data_cells-mean-stddev_%s-sims_recon-%s.pickle' % (nsims, cmb_data_type)
with open('%s/%s'%(outdir, fname), 'rb') as handle:
    wo_fgs = pickle.load(file=handle)

alpha = 0.5
lmax = 4000
ells = np.arange(lmax)
for key in with_fgs['cells'].keys():
    grid = plt.GridSpec(nrows=3, ncols=1, wspace=0.4, hspace=0.3)
    ax1 = plt.subplot(grid[0:2, 0])
    ax2 = plt.subplot(grid[2, 0])

    # plot the actual spectra first
    ax = ax1
    
    label = 'with fgs'
    y, yerr = with_fgs['cells'][key], with_fgs['cells_err'][key]
    p = ax.plot(ells, y, '--', linewidth=2, label=label)
    ax.fill_between(ells, y, y + yerr,
                     alpha=alpha, facecolor=p[0].get_color())
    ax.fill_between(ells, y, y - yerr,
                     alpha=alpha, facecolor=p[0].get_color())
    
    label = 'w/o fgs'
    y, yerr = wo_fgs['cells'][key], wo_fgs['cells_err'][key]
    ax.plot(ells, y, '.-', linewidth=2, label=label) #, color=p[0].get_color())
    ax.fill_between(ells, y, y + yerr,
                     alpha=alpha) #, facecolor=p[0].get_color())
    ax.fill_between(ells, y, y - yerr,
                     alpha=alpha) #, facecolor=p[0].get_color())

    # now plot the differences
    ax = ax2
    label = key
    y = with_fgs['cells'][key] - wo_fgs['cells'][key] 
    yerr = np.sqrt( with_fgs['cells_err'][key]**2 + wo_fgs['cells_err'][key]**2 )
    
    ax.plot(ells, y, '-', linewidth=2, label=label, color=p[0].get_color())
    ax.fill_between(ells, y, y + yerr,
                     alpha=alpha, facecolor=p[0].get_color())
    ax.fill_between(ells, y, y - yerr,
                     alpha=alpha, facecolor=p[0].get_color())

    for ax in [ax1, ax2]:
        ax.set_xscale('log')
        ax.set_yscale('log')

    ax1.set_ylabel('$C_\ell$')
    ax2.set_ylabel('$\Delta C_\ell$')
    ax2.set_xlabel('$\ell$')
    
    ax1.set_title(key)
    ax1.legend(bbox_to_anchor=(1,1))
    
    fname = 'plot_compare-spectra-%s_%s-sims_%s.png' % (key, nsims, cmb_data_type)
    plt.savefig('%s/%s' % (outdir, fname), format='png', bbox_inches='tight')
    plt.close()
    print('## saved spectra plots in %s' % fname)