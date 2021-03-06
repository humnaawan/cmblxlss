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

if outdir.__contains__('-nodust'):
    outdir = outdir.split('-nodust')[0]

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
# g field
filename = 'gal_density_alm.pickle'
gal_density_alm = read_pickle(filename='%s/%s'%(alms_dir, filename))
readme = print_update(update='\nReading in %s'%filename,
                      readme=readme)
# masked g field
filename = 'gal_density_alm_xmask.pickle'
gal_density_alm_masked = read_pickle(filename='%s/%s'%(alms_dir, filename))
readme = print_update(update='\nReading in %s'%filename,
                      readme=readme)
# lsst-modulated g field
filename = 'gal_density_alm_mod.pickle'
gal_density_alm_mod_wdust = read_pickle(filename='%s/%s'%(alms_dir, filename))
readme = print_update(update='\nReading in %s'%filename,
                      readme=readme)
# lsst-modulated x mask g field
filename = 'gal_density_alm_mod_xmask.pickle'
gal_density_alm_mod_masked_wdust = read_pickle(filename='%s/%s'%(alms_dir, filename))
readme = print_update(update='\nReading in %s'%filename,
                      readme=readme)
# lsst-modulated g field no dust
filename = 'gal_density_alm_mod-nodust.pickle'
gal_density_alm_mod_nodust = read_pickle(filename='%s/%s'%(alms_dir, filename))
readme = print_update(update='\nReading in %s'%filename,
                      readme=readme)
# lsst-modulated x mask g field w/ dust
filename = 'gal_density_alm_mod-nodust_xmask.pickle'
gal_density_alm_mod_masked_nodust = read_pickle(filename='%s/%s'%(alms_dir, filename))
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
# kappa kappa
c_ells[r'$\kappa\kappa$'] = hp.alm2cl(kappa_alms_normed)[0:lmax]
c_ells[r'$\kappa\kappa$ w/ fg'] = hp.alm2cl(kappa_alms_normed_fg)[0:lmax]
c_ells[r'$\kappa\kappa$ w/ lsst mask'] = hp.alm2cl(kappa_alms_normed_masked)[0:lmax]
c_ells[r'$\kappa\kappa$ w/ fg + lsst mask'] = hp.alm2cl(kappa_alms_normed_fg_masked)[0:lmax]

# gg
c_ells[r'$gg$'] = hp.alm2cl(gal_density_alm)[0:lmax]
c_ells[r'$gg$ w/ lsst mask'] = hp.alm2cl(gal_density_alm_masked)[0:lmax]
c_ells[r'$gg$ w/ modulated w/ dust'] = hp.alm2cl(gal_density_alm_mod_wdust)[0:lmax]
c_ells[r'$gg$ w/ modulated w/ dust + lsst mask'] = hp.alm2cl(gal_density_alm_mod_masked_wdust)[0:lmax]
c_ells[r'$gg$ w/ modulated w/o dust'] = hp.alm2cl(gal_density_alm_mod_nodust)[0:lmax]
c_ells[r'$gg$ w/ modulated w/o dust + lsst mask'] = hp.alm2cl(gal_density_alm_mod_masked_nodust)[0:lmax]

# cross correlation: correlate stuff s.t. both maps have the lsst mask
c_ells[r'$\kappa$ w/ lsst mask x $g$ w/ lsst mask'] = hp.alm2cl(kappa_alms_normed_masked, gal_density_alm_masked)[0:lmax]
c_ells[r'$\kappa$ w/ lsst mask x $g$ w/ lsst mask + modulation w/o dust'] = hp.alm2cl(kappa_alms_normed_masked, gal_density_alm_mod_masked_nodust)[0:lmax]
c_ells[r'$\kappa$ w/ lsst mask x $g$ w/ lsst mask + modulation w/ dust'] = hp.alm2cl(kappa_alms_normed_masked, gal_density_alm_mod_masked_wdust)[0:lmax]
c_ells[r'$\kappa$ w/ lsst mask + fg x $g$ w/ lsst mask + modulation w/o dust'] = hp.alm2cl(kappa_alms_normed_fg_masked, gal_density_alm_mod_masked_nodust)[0:lmax]
c_ells[r'$\kappa$ w/ lsst mask + fg x $g$ w/ lsst mask + modulation w/ dust'] = hp.alm2cl(kappa_alms_normed_fg_masked, gal_density_alm_mod_masked_wdust)[0:lmax]
# -----------------------------------------------
cls_dir = '%s/cls_dir/'%outdir
if not os.path.exists(cls_dir): os.makedirs(cls_dir)
# nominal bin width
bin_width = 50
# plot the cross-correlations in sets
# -----------------------------------------------
readme = print_update(update='## Plotting spectra + saving them in %s ... \n' % (cls_dir.split(outdir)[-1]),
                      readme=readme)
# -----------------------------------------------
# first plot all the kappa kappa auto corrs
c_ells_to_plot = {}
for key in [f for f in c_ells.keys() if f.__contains__(r'\kappa\kappa')]:
    c_ells_to_plot[key] = c_ells[key]
# plot
filename = plot_cls_dict(cls_in=c_ells_to_plot, outdir=cls_dir, file_tag='kk-only',
                         save_plot=True, show_plot=False,
                         ylog=True, xlog=False, lmin=None, lmax=lmax,
                         binned=False, bin_width=bin_width)
readme = print_update(update='Saved %s'%filename,
                      readme=readme)
# plot binned
filename = plot_cls_dict(cls_in=c_ells_to_plot, outdir=cls_dir, file_tag='kk-only',
                         save_plot=True, show_plot=False,
                         ylog=True, xlog=False, lmin=None, lmax=lmax,
                         binned=True, bin_width=bin_width)
readme = print_update(update='Saved %s'%filename,
                      readme=readme)
# -----------------------------------------------
# now plot all the gg auto corrs
c_ells_to_plot = {}
for key in [f for f in c_ells.keys() if f.__contains__(r'gg')]:
    c_ells_to_plot[key] = c_ells[key]
# plot
filename = plot_cls_dict(cls_in=c_ells_to_plot, outdir=cls_dir, file_tag='gg-only',
                         save_plot=True, show_plot=False,
                         ylog=True, xlog=False, lmin=None, lmax=lmax,
                         binned=False, bin_width=bin_width)
readme = print_update(update='Saved %s'%filename,
                      readme=readme)
# plot binned
filename = plot_cls_dict(cls_in=c_ells_to_plot, outdir=cls_dir, file_tag='gg-only',
                         save_plot=True, show_plot=False,
                         ylog=True, xlog=False, lmin=None, lmax=lmax,
                         binned=True, bin_width=bin_width)
readme = print_update(update='Saved %s'%filename,
                      readme=readme)
# -----------------------------------------------
# and now plot the cross correlations
c_ells_to_plot = {}
for key in [f for f in c_ells.keys() if f.__contains__(' x ')]:
    c_ells_to_plot[key] = c_ells[key]

# set up colors, markers
markers = ['+', 'x', '.', '3', 'x']
colors = ['b', 'orangered', 'm', 'darkcyan', 'k']
linestyles = ['-', '-', '--', '--', '-']
# plot
filename = plot_cls_dict(cls_in=c_ells_to_plot, outdir=cls_dir, file_tag='kg-only',
                         save_plot=True, show_plot=False,
                         ylog=True, xlog=False, lmin=None, lmax=lmax,
                         colors=colors, markers=markers, linestyles=linestyles,
                         binned=False, bin_width=bin_width)
readme = print_update(update='Saved %s'%filename,
                      readme=readme)
# binned
filename = plot_cls_dict(cls_in=c_ells_to_plot, outdir=cls_dir, file_tag='kg-only',
                         save_plot=True, show_plot=False,
                         ylog=True, xlog=False, lmin=None, lmax=lmax,
                         colors=colors, markers=markers, linestyles=linestyles,
                         binned=True, bin_width=bin_width)
readme = print_update(update='Saved %s'%filename,
                      readme=readme)
# binned: residuals
filename = plot_cls_dict(cls_in=c_ells_to_plot, outdir=cls_dir, file_tag='kg-only',
                         save_plot=True, show_plot=False,
                         colors=colors, markers=markers, linestyles=linestyles,
                         residuals=True, baseline_key=r'$\kappa$ w/ lsst mask x $g$ w/ lsst mask',
                         ylog=False, xlog=True, lmin=None, lmax=lmax,
                         binned=True, bin_width=bin_width)
readme = print_update(update='Saved %s'%filename,
                      readme=readme)
# binned: residuals: zoomed: all
filename = plot_cls_dict(cls_in=c_ells_to_plot, outdir=cls_dir, file_tag='kg-only',
                         save_plot=True, show_plot=False,
                         colors=colors, markers=markers, linestyles=linestyles,
                         residuals=True, baseline_key=r'$\kappa$ w/ lsst mask x $g$ w/ lsst mask',
                         ylog=False, xlog=True, lmin=10, lmax=lmax,
                         binned=True, bin_width=30)
readme = print_update(update='Saved %s'%filename,
                      readme=readme)
# ----------------------------------
# plot a subset
# set up colors, markers
markers = ['+', '+', 'x']
colors = ['b', 'orangered',  'k']
linestyles = ['-', '-', '-']

c_ells_to_plot = {}
keys = [r'$\kappa$ w/ lsst mask x $g$ w/ lsst mask',
        r'$\kappa$ w/ lsst mask x $g$ w/ lsst mask + modulation w/o dust',
        r'$\kappa$ w/ lsst mask + fg x $g$ w/ lsst mask + modulation w/ dust']
for key in keys:
    c_ells_to_plot[key] = c_ells[key]

# binned
filename = plot_cls_dict(cls_in=c_ells_to_plot, outdir=cls_dir, file_tag='kg-only',
                         save_plot=True, show_plot=False,
                         ylog=True, xlog=True, lmin=None, lmax=lmax,
                         colors=colors, markers=markers, linestyles=linestyles,
                         binned=True, bin_width=bin_width)
readme = print_update(update='Saved %s'%filename,
                      readme=readme)

# binned: residuals: zoomed: some
filename = plot_cls_dict(cls_in=c_ells_to_plot, outdir=cls_dir, file_tag='kg-only',
                         save_plot=True, show_plot=False,
                         colors=colors, markers=markers, linestyles=linestyles,
                         residuals=True, baseline_key=r'$\kappa$ w/ lsst mask x $g$ w/ lsst mask',
                         ylog=False, xlog=True, lmin=10, lmax=lmax,
                         binned=True, bin_width=30)
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
