import matplotlib as mpl
import matplotlib.pyplot as plt
import healpy as hp
import numpy as np
import pandas as pd
import time
from orphics.catalogs import CatMapper
# ------------------------------------------------------------------------------
# imports from this folder
from utils_plot import plot_seaborn_distplot, plot_seaborn_jointplot, plot_mollview, plot_power_spectrum
from settings import rcparams
# set up the plotting params
for key in rcparams: mpl.rcParams[key] = rcparams[key]
# ------------------------------------------------------------------------------

__all__= ['get_gal_density_map', 'get_convergence_map',
          'get_lsst_completeness_map', 'modulate_map', 'get_smoothed_gal_density_map']
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
def get_gal_density_map(halocat_filepath, halocat_tag, nside,
                        plot_data_diagnostics, outdir, verbose=True):
	time0 = time.time()
	readme = '\n*** get_gal_density_map function ***\n'
	map_dict = {}
	# Read in the halo catalog. Plots some diagnostic plots.
	halocat_data = pd.read_csv(halocat_filepath)
	readme += 'Read in halo cat. data. Shape: %s\n'%(np.shape(halocat_data),)

	# add ra, dec columns based on theta, phi columns on the image plane
	halocat_data['ra'] = np.degrees(halocat_data['phi_i'])
	halocat_data['dec'] = np.degrees(np.pi/2.0 - halocat_data['theta_i'])
	readme += 'Added ra, dec cols.\n'

	if plot_data_diagnostics:
		filenames = plot_seaborn_distplot(data_dict=halocat_data,
                                          keys_to_plot=list(halocat_data.keys()),
                                          outdir=outdir,
                                          file_tag='%s'%halocat_tag,
                                          save_plot=True, show_plot=False)
		filename = plot_seaborn_jointplot(x_in=halocat_data['ra'],
                                          y_in=halocat_data['dec'],
                                          outdir=outdir,
                                          file_tag='%s'%halocat_tag,
                                          save_plot=True, show_plot=False)
		readme += '\nSaved diagnostic plots:\n%s\n\n'%(filenames + [filename])

	#--------------------------------------------------------------------------------------
	# Populate the galaxies on HP map
	halocat_catmap = CatMapper(ras_deg=halocat_data['ra'],
                               decs_deg=halocat_data['dec'], nside=nside, verbose=verbose)
	readme += 'Populated the halo catalog galaxies on an hp map with nside %s\n'%nside

	filename = plot_mollview(map_in=halocat_catmap.counts,
                             title='halo catalog: galaxy counts',
                             data_label='counts',
                             outdir=outdir,
                             file_tag='%s'%halocat_tag,
                             save_plot=True,show_plot=False)
	readme += 'Saved mollview plot for galaxy counts in %s\n'%(filename)

	# get galaxy density
	map_dict['halocat_delta'] = ( halocat_catmap.counts / halocat_catmap.counts.mean() )-1.

	# plot skymap
	filename = plot_mollview(map_in=map_dict['halocat_delta'],
                             title='halocat_delta_map',
                             data_label='dNbyN',
                             outdir=outdir,
                             file_tag='%s'%halocat_tag,
                             save_plot=True, show_plot=False)
	readme += 'Saved mollview plot for galaxy density in %s\n'%(filename)

	# power spectrum of the overdensity
	filename = plot_power_spectrum(map_in=map_dict['halocat_delta'],
                                   lmax=4*nside, title='halocat_delta',
                                   outdir=outdir, file_tag='%s'%halocat_tag,
                                   return_cls=False,
                                   save_plot=True, show_plot=False)
	readme += 'Saved powerspectrum for galaxy density in %s\n'%(filename)

	readme += '\nTime taken: %.3f min\n\n'%((time.time()-time0)/60.)
	return [map_dict, readme]

# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
def get_convergence_map(convergence_filepath, wldata_tag, nside, outdir):
	time0 = time.time()
	readme = '\n*** get_convergence_map function ***\n'
	# read in the data
	datain = open(convergence_filepath,'rb')
	nres = np.fromfile(datain, count=1, dtype=np.int32)[0] # nside = 2^nres
	wl_data = np.fromfile(datain, dtype=np.float32)
	datain.close()

	# order of things based on the provided Fortran code that's supposed to help read in the file
	# http://cosmo.phys.hirosaki-u.ac.jp/takahasi/allsky_raytracing/dat2fits.f90
	npix = hp.nside2npix(nside=2**nres)

	readme += 'WL data read in: npix= %s ; shape(arr): %s\n'%(npix, np.shape(wl_data),)

	if (len(wl_data)!= 4.*npix):
		readme += '!!!!! there are more entries in the .dat file than for 4 healpix maps!?!?!?!\n'

	maps = {}
	maps['kmap'] = wl_data[0:npix]
	maps['shear1'] = wl_data[npix:2*npix]
	maps['shear2'] = wl_data[2*npix:3*npix]
	maps['rotation'] = wl_data[3*npix:4*npix]

	# check to make sure the completenes map is the same nside as wanted
	wl_nside = 2**nres
	if (wl_nside != nside):
		readme+='Converting the read-in map from nside %s to nside %s\n'%(wl_nside, nside)
		for key in maps:
			maps[key] = hp.ud_grade(map_in=maps[key], nside_out=nside)
			readme +='Conversion of the %s map done.\n'%key

	for key in maps:
		# save skymap
		filename = plot_mollview(map_in=maps[key], title=key, data_label='',
                                 outdir=outdir, file_tag='%s_%s'%(key, wldata_tag),
                                 save_plot=True, show_plot=False)
		readme += 'Saved mollview plot for %s in %s\n'%(key, filename)
	readme += '\nTime taken: %.3f min\n\n'%((time.time()-time0)/60.)
	return [{'convergence': maps['kmap']}, readme]

# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
def get_lsst_completeness_map(lsst_completeness_path, lsstdata_tag,
                              completeness_threshold, nside, outdir):
	time0 = time.time()
	readme = '\n*** get_lsst_completeness_map function ***\n'

	lsst_completeness = np.load(lsst_completeness_path)

	# plot skymap
	filename = plot_mollview(map_in=lsst_completeness['metricValues'],
                             title=lsstdata_tag,
                             data_label='',
                             outdir=outdir,
                             file_tag='%s'%lsstdata_tag,
                             save_plot=True, show_plot=False)
	readme += 'Saved mollview lsst galaxy density map in %s\n'%(filename)

	# Construct a completeness map based on the lsst galaxy density map
	map_dict = {}
	map_dict['lsst_completeness'] = np.zeros(len(lsst_completeness['metricValues']))
	ind= np.where((lsst_completeness['metricValues'] > completeness_threshold) & \
	              (lsst_completeness['mask'] == False))[0]
	map_dict['lsst_completeness'][ind] = 1.

	# plot skymap of the completeness map
	filename = plot_mollview(map_in=map_dict['lsst_completeness'],
                             title='completeness map: threshold %s'%completeness_threshold,
                             data_label='',
                             outdir=outdir,
                             file_tag='completness_%s'%lsstdata_tag,
                             save_plot=True, show_plot=False)
	readme += 'Saved mollview lsst completeness map in %s\n'%(filename)

	# check to make sure the completenes map is the same nside as wanted
	lsst_nside = hp.npix2nside(lsst_completeness['slicerShape'])
	if (lsst_nside != nside):
		readme += 'Converting the completness map from nside %s to nside %s\n'%(lsst_nside, nside)
		map_dict['lsst_completeness'] = hp.ud_grade(map_in=map_dict['lsst_completeness'],
                                                    nside_out=nside)
		readme +='Conversion successful.\n'

	fsky = float(len(lsst_completeness['mask'][lsst_completeness['mask'] == False]))/len(lsst_completeness['metricValues'])
	readme += 'fsky: %s\n'%(fsky)

	readme += '\nTime taken: %.3f min\n\n'%((time.time()-time0)/60.)
	return [map_dict, readme, fsky]

# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
def modulate_map(data_map, completeness_map, data_tag, outdir):
    time0 = time.time()
    readme = '\n*** modulate_map function for %s ***\n'%data_tag
    
    if (len(data_map) != len(completeness_map)):
        err = 'Input data_map and completeness_map must be of the same length.'
        err += ' Not: %s vs. %s'%(len(data_map), len(completeness_map))
        raise ValueError(err)
    mod_map = data_map*completeness_map
    # plot skymap of the modulated map
    filename = plot_mollview(map_in=mod_map,
                             title=data_tag,
                             data_label='',
                             outdir=outdir,
                             file_tag='modulated_%s'%data_tag,
                             save_plot=True, show_plot=False)
    readme += 'Saved mollview plot of the modulated map in %s\n'%(filename)
    ### need to make these masked arrays ..
    return [mod_map, readme]

# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
def get_smoothed_gal_density_map(map_in, fwhm, outdir, data_tag):
	# fwhm in radian
	time0 = time.time()
	readme = '\n*** get_smoothed_gal_density_map function ***\n'
	map_smoothed = hp.smoothing(map_in, fwhm=fwhm)

	# plot skymap of the modulated map
	filename = plot_mollview(map_in=map_smoothed,
                             title='smoothed galaxy density',
                             data_label='',
                             outdir=outdir,
                             file_tag='smoothed_%s'%data_tag,
                             save_plot=True, show_plot=False)
	readme += 'Saved mollview plot of the smoothed galaxy density map in %s\n'%(filename)
	readme += '\nTime taken: %.3f min\n\n'%((time.time()-time0)/60.)

	return [map_smoothed, readme]
