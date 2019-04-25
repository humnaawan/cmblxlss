################################################################################
# This script reads in the saved alms for kappa and galaxy density (w different
# variations, e.g. with lsst mask and without), calculates the various spectra,
# plots them and saves them.
#
# ** This script needs the alms saved to outdir/alms_dir ** (acheived by running
# save_needed_alms.py script.)
################################################################################
import matplotlib
matplotlib.use('Agg')
import numpy as np
import healpy as hp
import pandas as pd
import os
from collections import OrderedDict
import time, datetime
from optparse import OptionParser
# ------------------------------------------------------------------------------
# imports from this folder
from utils_plot import plot_cls_dict
from utils_files import read_pickle
from utils_analysis import print_update
#-------------------------------------------------------------------------------
parser = OptionParser()
# required params
parser.add_option('--outdir', dest='outdir',
                  help='Output directory')
#-------------------------------------------------------------------------------
# parse the options
(options, args) = parser.parse_args()
outdir = options.outdir
# ------------------------------------------------------------------------------
time0 = time.time()
readme = '\n------------------------------------------------------------\n'
readme += '%s\n'%datetime.datetime.now()
readme += '## Running get_spectra.'
# ------------------------------------------------------------------------------
# read in the alms
alms_dir = '%s/alms_dir/'%outdir

# baseline kappa alms
filename = 'kappa_alms_normed.pickle'
kappa_alms_normed = read_pickle(filename='%s/%s'%(alms_dir, filename))
readme = print_update(update='\nReading in %s'%filename,
                      readme=readme)
# kappa alms with fg
filename = 'kappa_alms_normed_fg.pickle'
kappa_alms_normed_fg = read_pickle(filename='%s/%s'%(alms_dir, filename))
readme = print_update(update='\nReading in %s'%filename,
                      readme=readme)
# kappa alms with lsst mask
filename = 'kappa_alms_normed_masked.pickle'
kappa_alms_normed_masked = read_pickle(filename='%s/%s'%(alms_dir, filename))
readme = print_update(update='\nReading in %s'%filename,
                      readme=readme)
# kappa alms with fg + lsst mask
filename = 'kappa_alms_normed_fg_masked.pickle'
kappa_alms_normed_fg_masked = read_pickle(filename='%s/%s'%(alms_dir, filename))
readme = print_update(update='\nReading in %s'%filename,
                      readme=readme)
# correlated g field
filename = 'gal_density_alm.pickle'
gal_density_alm = read_pickle(filename='%s/%s'%(alms_dir, filename))
readme = print_update(update='\nReading in %s'%filename,
                      readme=readme)
# lsst-modulated correlated g field
filename = 'gal_density_alm_mod.pickle'
gal_density_alm_mod = read_pickle(filename='%s/%s'%(alms_dir, filename))
readme = print_update(update='\nReading in %s'%filename,
                      readme=readme)
# lsst-modulated x mask correlated g field
filename = 'gal_density_alm_mod_xmask.pickle'
gal_density_alm_mod_masked = read_pickle(filename='%s/%s'%(alms_dir, filename))
readme = print_update(update='\nReading in %s'%filename,
                      readme=readme)
# misc
filename = 'misc_info.pickle'
misc_info = read_pickle(filename='%s/%s'%(alms_dir, filename))
readme = print_update(update='\nReading in %s'%filename,
                      readme=readme)
lsst_fsky = misc_info['lsst_fsky']
lmax = misc_info['lmax']

# ------------------------------------------------------------------------------
# set up the readme
readme_tag = 'lmax%s_get-spectra'%(lmax)
# update readme
readme = print_update(update='\nOutput directory: %s\n'%outdir,
                      readme=readme)

# set up dictionary that will contain all the spectra we care about
c_ells = OrderedDict()
# set up the ells
ells = np.arange(0, lmax)
# -----------------------------------------------
# calculate correlations
readme = print_update(update='## Calculating correlations ... \n',
                      readme=readme)
# autos
#c_ells[r'$\kappa\kappa$ baseline'] = hp.alm2cl(kappa_alms_normed, lmax=lmax)[0:lmax]
c_ells[r'$\kappa\kappa$ w/ fg'] = hp.alm2cl(kappa_alms_normed_fg, lmax=lmax)[0:lmax]
c_ells[r'$\kappa\kappa$ w/ lsst mask'] = hp.alm2cl(kappa_alms_normed_masked, lmax=lmax)[0:lmax]/lsst_fsky
c_ells[r'$\kappa\kappa$ w/ fg + lsst mask'] = hp.alm2cl(kappa_alms_normed_fg_masked, lmax=lmax)[0:lmax]/lsst_fsky
# cross correlation
#c_ells[r'$\kappa$ baseline x correlated $g$'] = hp.alm2cl(kappa_alms_normed, gal_density_alm, lmax=lmax)[0:lmax]
#c_ells[r'$\kappa$ baseline x lsst modulated correlated $g$'] = hp.alm2cl(kappa_alms_normed, gal_density_alm_mod, lmax=lmax)[0:lmax]
#c_ells[r'$\kappa$ w/ lsst mask x correlated $g$'] = hp.alm2cl(kappa_alms_normed_masked, gal_density_alm, lmax=lmax)[0:lmax]/lsst_fsky
c_ells[r'$\kappa$ w/ lsst mask x $g$ w/ lsst mask'] = hp.alm2cl(kappa_alms_normed_masked, gal_density_alm_mod, lmax=lmax)[0:lmax]
c_ells[r'$\kappa$ w/ lsst mask x $g$ w/ lsst mask + modulation'] = hp.alm2cl(kappa_alms_normed_masked, gal_density_alm_mod_masked, lmax=lmax)[0:lmax]
#c_ells[r'$\kappa$ w/ fg x correlated $g$'] = hp.alm2cl(kappa_alms_normed_fg, gal_density_alm, lmax=lmax)[0:lmax]
#c_ells[r'$\kappa$ w/ fg + lsst mask x correlated $g$'] = hp.alm2cl(kappa_alms_normed_fg_masked, gal_density_alm, lmax=lmax)[0:lmax]/lsst_fsky
c_ells[r'$\kappa$ w/ lsst mask + fg x $g$ w/ lsst mask + modulation'] = hp.alm2cl(kappa_alms_normed_fg_masked, gal_density_alm_mod_masked, lmax=lmax)[0:lmax]
# -----------------------------------------------
cls_dir = '%s/cls_dir/'%outdir
if not os.path.exists(cls_dir): os.makedirs(cls_dir)

# plot the cross-correlations in sets
readme = print_update(update='## Plotting spectra + saving them in %s ... \n' % (cls_dir.split(outdir)[-1]),
                      readme=readme)
bin_width = 50
# first plot the auto spectra for checks
c_ells_to_plot = {}
for key in [f for f in c_ells.keys() if f.__contains__(r'\kappa\kappa')]:
    c_ells_to_plot[key] = c_ells[key]
# plot
filename = plot_cls_dict(cls_in=c_ells_to_plot, outdir=cls_dir, file_tag='kk-only',
                         save_plot=True, show_plot=False,
                         cross_convention=True,
                         sci_yticks=True,
                         binned=False, bin_width=bin_width, lmax=lmax)
readme = print_update(update='Saved %s'%filename,
                      readme=readme)
# plot binned
filename = plot_cls_dict(cls_in=c_ells_to_plot, outdir=cls_dir, file_tag='kk-only',
                         save_plot=True, show_plot=False,
                         cross_convention=True,
                         sci_yticks=True,
                         binned=True, bin_width=bin_width, lmax=lmax)
readme = print_update(update='Saved %s'%filename,
                      readme=readme)
# -----------------------------------------------
# now plot the cross correlations
c_ells_to_plot = {}
for key in [f for f in c_ells.keys() if f.__contains__(' x ')]:
    c_ells_to_plot[key] = c_ells[key]

# set up colors, markers
markers = ['+', '1', '.', 'x']
colors = ['b', 'm', 'orangered', 'k']
# plot
filename = plot_cls_dict(cls_in=c_ells_to_plot, outdir=cls_dir, file_tag='kg-only',
                         save_plot=True, show_plot=False,
                         cross_convention=True,
                         sci_yticks=True, colors=colors, markers=markers,
                         binned=False, bin_width=bin_width, lmax=lmax)
readme = print_update(update='Saved %s'%filename,
                      readme=readme)
# binned
filename = plot_cls_dict(cls_in=c_ells_to_plot, outdir=cls_dir, file_tag='kg-only',
                         save_plot=True, show_plot=False,
                         cross_convention=True,
                         sci_yticks=True, colors=colors, markers=markers,
                         binned=True, bin_width=bin_width, lmax=lmax)
readme = print_update(update='Saved %s'%filename,
                      readme=readme)
# binned: residuals
filename = plot_cls_dict(cls_in=c_ells_to_plot, outdir=cls_dir, file_tag='kg-only',
                         save_plot=True, show_plot=False,
                         cross_convention=True, colors=colors, markers=markers,
                         residuals=True, baseline_key=r'$\kappa$ w/ lsst mask x $g$ w/ lsst mask',
                         sci_yticks=True, loglog=False,
                         binned=True, bin_width=bin_width, lmax=lmax)
readme = print_update(update='Saved %s'%filename,
                      readme=readme)
# binned: residuals: zoomed
filename = plot_cls_dict(cls_in=c_ells_to_plot, outdir=cls_dir, file_tag='kg-only',
                         save_plot=True, show_plot=False,
                         cross_convention=True, colors=colors, markers=markers,
                         residuals=True, baseline_key=r'$\kappa$ w/ lsst mask x $g$ w/ lsst mask',
                         sci_yticks=True, loglog=True,
                         binned=True, bin_width=20, lmax=lmax, lmin=10)
readme = print_update(update='Saved %s'%filename,
                      readme=readme)
# -----------------------------------------------
# save the corrs
filename = 'all_spectra.csv'
pd.DataFrame(c_ells).to_csv('%s/%s'%(cls_dir, filename), index=False)
# update readme
readme = print_update(update='Saved cross-spectra in %s.\n' % filename,
                      readme=readme)
readme = print_update(update='\nTime taken for the whole thing: %.3f min\n' % ((time.time()-time0)/60.),
                      readme=readme)
# write the readme
readme_file = open('%s/readme_%s.txt' % (outdir, readme_tag), 'a')
readme_file.write(readme)
readme_file.close()
