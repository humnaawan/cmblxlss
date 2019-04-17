import matplotlib
matplotlib.use('Agg')
import numpy as np
import healpy as hp
import pandas as pd
import os
from collections import OrderedDict
import time
from optparse import OptionParser

#--------------------------------------------------------------------------------------
# imports from this folder
from utils_plot import plot_cls_dict
from utils_analysis import *
#--------------------------------------------------------------------------------------
parser = OptionParser()
# required params
parser.add_option('--outdir', dest='main_outdir',
                  help='Output directory')
parser.add_option('--halocat-path', dest='halocat_filepath',
                  help='Path to the halo catalog')
parser.add_option('--convergence-path', dest='convergence_filepath',
                  help='Path to the file containing convergence map')
parser.add_option('--cmb-sim-tag', dest='simtag',
                  help='Tag to identify simulation realization for the halo catalog and the convergence, etc.')
parser.add_option('--lsst-mask-path', dest='lsst_completeness_path',
                  help='Path to the opsim output with lsst completeness map')
parser.add_option('--lsstdata-tag', dest='lsstdata_tag',
                  help='Tag to identify lsst map')
# optional
parser.add_option('-q', '--quiet',
                  action='store_false', dest='verbose', default=True,
                  help='Do not print stuff out')
parser.add_option('--nside', dest='nside', type='int',
                  help='HEALPix resolution parameter.', default=256)
parser.add_option('--lmax', dest='lmax', type='int',
                  help='Maximum multipole to consider', default=1000)
parser.add_option('--cross_corr_rec',
                  action='store_false', dest='cross_corr_rec', default=False,
                  help='Cross-correlate things with reconstructions. Else with convergence map. ')
parser.add_option('--completeness-thresh', dest='completeness_threshold',
                  help='Threshold on deltaNByN to define the binary (0/1) completeness map',
                  default=-0.2)
parser.add_option('--smoothing-fwhm', dest='smoothing_fwhm',
                  help='FWHM (in radians) of the gaussian beam used to smooth the halo density map',
                  default=np.radians(2/60.))
parser.add_option('--plot-diagnostic',
                  action='store_false', dest='plot_data_diagnostics', default=False,
                  help='Plot out diagnostics on the halo catalog data.')
#--------------------------------------------------------------------------------------
# parse the options
(options, args) = parser.parse_args()
verbose = options.verbose
nside = options.nside
lmax = options.lmax
main_outdir = options.main_outdir
simtag = options.simtag
halocat_filepath = options.halocat_filepath
cross_corr_convergence = not options.cross_corr_rec
convergence_filepath = options.convergence_filepath
lsst_completeness_path = options.lsst_completeness_path
lsstdata_tag = options.lsstdata_tag
completeness_threshold = options.completeness_threshold
smoothing_fwhm = options.smoothing_fwhm
plot_data_diagnostics = options.plot_data_diagnostics

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
time0 = time.time()

readme = ''
readme_tag = 'nside%s_lmax%s_takahashi-%s'%(nside, lmax, simtag)

halocat_tag = halocat_filepath.split('_')[-1].split('.csv')[0]
outdir = '%s/output_%s_halocat-%s_lsst-%s'%(main_outdir, readme_tag,
                                            halocat_tag,
                                            lsstdata_tag)
if not os.path.exists(outdir): os.makedirs(outdir)

update = '\nOutput directory: %s\n'%outdir.split(main_outdir)[1]
readme += update
if verbose: print(update)

#--------------------------------------------------------------------------------------
# dictionary that will contain all the maps we care about
maps = OrderedDict()
c_ells = OrderedDict()
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
# get galaxy density map based on the halo catalogs
update = 'Creating the halo galaxy density map ... \n'
readme += update
if verbose: print(update)

out = get_gal_density_map(halocat_filepath=halocat_filepath,
                          halocat_tag=halocat_tag, nside=nside,
                          plot_data_diagnostics=plot_data_diagnostics,
                          outdir=outdir, verbose=verbose)
gal_density_key = list(out[0].keys())[0]
maps.update(out[0])
readme += out[1]

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
conv_recon_tag = ''
if cross_corr_convergence:
	# read in the Takahashi convergence map
	update = 'Reading in the convergence map ... \n'
	readme += update
	if verbose: print(update)
	conv_recon_tag = 'convergence'
	out = get_convergence_map(convergence_filepath=convergence_filepath,
                              wldata_tag=simtag, nside=nside,
                              outdir=outdir)
	completeness_key = list(out[0].keys())[0]
	maps.update(out[0])
	readme += out[1]
else:
	# read in reconstruction
	raise ValueError('reconstruction read-in not developed rn.')
	conv_recon_tag = 'reconstruction'

cmb_key = list(out[0].keys())[0]
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
# create lsst completeness map
update = 'Reading in the lsst completeness map ... \n'
readme += update
if verbose: print(update)
out = get_lsst_completeness_map(lsst_completeness_path=lsst_completeness_path,
                                lsstdata_tag=lsstdata_tag,
                                completeness_threshold=completeness_threshold,
                                nside=nside, outdir=outdir)
completeness_key = list(out[0].keys())[0]
maps.update(out[0])
readme += out[1]
fsky = out[2]

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
# get a smoothed galaxy density map
update = 'Finding a smoothed galaxy density map ... \n'
readme += update
if verbose: print(update)
out = get_smoothed_gal_density_map(map_in=maps[gal_density_key],
                                   fwhm=smoothing_fwhm,
                                   outdir=outdir, data_tag=gal_density_key)

smoothed_gal_density_key = '%s_smoothed' % gal_density_key
maps[smoothed_gal_density_key] = out[0]
readme += out[1]

#--------------------------------------------------------------------------------------
update = 'Maps so far: %s\n'%(maps.keys())
readme += update
if verbose: print(update)
non_comp_keys = [f for f in maps.keys() if f != completeness_key]
#--------------------------------------------------------------------------------------
# modulate the original/smoothed halo density and convergence/recon. map
update = 'Modulating the galaxy density and %s maps ... \n' % conv_recon_tag
readme += update
if verbose: print(update)

for map_key in non_comp_keys:
    out = modulate_map(data_map=maps[map_key],
                       completeness_map=maps[completeness_key],
                       data_tag=map_key, outdir=outdir)
    maps['%s_modulated'%map_key] = out[0]
    readme += out[1]

update = '\nMaps so far now, after modulation: %s\n'%(maps.keys())
readme += update
if verbose: print(update)
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
# calculate the cross-correlations
update = '\nCalculating the cross power spectra ... \n'
readme += update
if verbose: print(update)
key1, key2 = gal_density_key, cmb_key
c_ells['True: %s x %s'%(key1, key2)] = hp.anafast(maps[key1], maps[key2], lmax=lmax)

key1, key2= '%s_modulated'%gal_density_key, '%s_modulated'%cmb_key
c_ells['%s x %s'%(key1, key2)] = hp.anafast(maps[key1], maps[key2], lmax=lmax)/fsky

key1, key2= smoothed_gal_density_key, cmb_key
c_ells['%s x %s'%(key1, key2)] = hp.anafast(maps[key1], maps[key2], lmax=lmax)

key1, key2= '%s_modulated'%smoothed_gal_density_key, '%s_modulated'%cmb_key
c_ells['%s x %s'%(key1, key2)] = hp.anafast(maps[key1], maps[key2], lmax=lmax)/fsky

# plot the cross-correlations
filename = plot_cls_dict(cls_in=c_ells, outdir=outdir, file_tag='',
                         save_plot=True, show_plot=False,
                         cross_convention=True,
                         sci_yticks=True,
                         binned=False, bin_width=20, lmax=lmax)
update = 'Saved cross-spectra plot in %s.\n'%filename
readme += update
if verbose: print(update)

# plot the binned cross-correlations
filename = plot_cls_dict(cls_in=c_ells, outdir=outdir, file_tag='',
                         save_plot=True, show_plot=False,
                         cross_convention=True,
                         sci_yticks=True,
                         binned=True, bin_width=50, lmax=lmax)
update = 'Saved cross-spectra plot in %s.\n'%filename
readme += update
if verbose: print(update)

# save the corrs
filename = 'cross_specs.csv'
pd.DataFrame(c_ells).to_csv('%s/%s'%(outdir, filename), index=False)
update = 'Saved cross-spectra in %s.\n' % filename
readme += update
if verbose: print(update)

update = '\nTime taken for the whole thing: %.3f min\n' % ((time.time()-time0)/60.)
readme += update
if verbose: print(update)

readme_file = open('%s/readme_%s.txt' % (outdir, readme_tag), 'a')
readme_file.write(readme)
readme_file.close()
