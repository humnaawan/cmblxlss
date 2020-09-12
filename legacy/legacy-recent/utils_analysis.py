import matplotlib as mpl
import matplotlib.pyplot as plt
import healpy as hp
import numpy as np
import pandas as pd
import time
from healpy import pixelfunc
from falafel import qe
from pixell import enmap
# ------------------------------------------------------------------------------
# imports from this folder
from utils_plot import plot_mollview
from settings import rcparams
# set up the plotting params
for key in rcparams: mpl.rcParams[key] = rcparams[key]
# ------------------------------------------------------------------------------
__all__= ['print_update', 'get_lsst_maps', 'get_reconstructed_kappa_alm',
          'generate_correlated_alm', 'rotate_map']
# ------------------------------------------------------------------------------
def print_update(update, readme):
    print(update)
    return '%s\n%s'%(readme, update)
# ------------------------------------------------------------------------------
def get_lsst_maps(data_file, data_tag, data_label, nside_out,
                  completeness_threshold, smoothing_fwhm, outdir,
                  plot_interm=True):
    time0 = time.time()
    readme = '\n*** Running get_lsst_maps function ***\n'
    # read in the data
    data_in = np.load(data_file)
    lsst_nside = hp.npix2nside(data_in['slicerShape'])
    # see if need to plot intermediate plots
    if plot_interm:
        # plot skymap of the input data
        filename = plot_mollview(map_in=data_in['metricValues'].copy(),
                                 title='%s : read in'%(data_tag),
                                 data_label=data_label,
                                 outdir=outdir,
                                 file_tag='%s_%s'%(data_tag, data_label),
                                 save_plot=True, show_plot=False)
        # update readme
        readme = print_update(update='Saved mollview map of read-in data in %s\n'%(filename),
                              readme=readme)
    # -----------------------------------------------
    # Construct a completeness map based on the input map
    mask = np.zeros(len(data_in['metricValues']))
    insurvey = np.where((data_in['mask'] == False))[0]
    ind = np.where(data_in['metricValues'][insurvey] > completeness_threshold)[0]
    mask[insurvey[ind]] = 1.
    # just to be sure: zero out the data in all the masked pixels
    data_in['metricValues'][mask == 1] = 0.

    # see if need to plot intermediate plots
    if plot_interm:
        # plot skymap of the completeness map
        filename = plot_mollview(map_in=mask,
                                 title='completeness map: threshold %s'%completeness_threshold,
                                 data_label='mask',
                                 outdir=outdir,
                                 file_tag='lsstmask_%s'%(data_tag),
                                 save_plot=True, show_plot=False)
        readme = print_update(update='Saved mollview lsst completeness map in %s\n'%(filename),
                              readme=readme)
    # -----------------------------------------------
    # check to make sure the completenes map is the same nside as wanted
    if (lsst_nside != nside_out):
        readme = print_update(update='Converting lsst maps from nside %s to nside %s\n'%(lsst_nside, nside_out),
                              readme=readme)
        mask = hp.ud_grade(map_in=mask, nside_out=nside_out)
        data_map = hp.ud_grade(map_in=data_in['metricValues'], nside_out=nside_out)
    # -----------------------------------------------
    # smooth the mask
    readme = print_update(update='Smoothing the mask with fwhm %s radians\n'%(smoothing_fwhm),
                          readme=readme)
    mask_smoothed = hp.smoothing(mask, fwhm=smoothing_fwhm)
    # -----------------------------------------------
    # now plot the to-output maps
    # mask
    filename = plot_mollview(map_in=mask_smoothed.copy(),
                             title='apodized completeness mask',
                             data_label='',
                             outdir=outdir,
                             file_tag='apodized_lsstmask',
                             save_plot=True, show_plot=False)
    readme = print_update(update='Saved the binary apodized mask in %s\n'%(filename),
                          readme=readme)
    # data map
    filename = plot_mollview(map_in=data_map.copy(),
                             title=data_tag,
                             data_label='',
                             outdir=outdir,
                             file_tag='%s_%s'%(data_tag, data_label),
                             save_plot=True, show_plot=False)
    readme = print_update(update='Saved the %s map in %s\n'%(data_label, filename),
                          readme=readme)
    # -----------------------------------------------
    # now calculate the fsky
    fsky = float(len(insurvey))/len(data_in['metricValues'])
    readme = print_update(update='fsky: %s\n'%(fsky),
                          readme=readme)
    # update readme
    readme += '\nTime taken: %.3f min\n\n'%((time.time()-time0)/60.)
    return [mask_smoothed, data_map, fsky, readme]

# ------------------------------------------------------------------------------
def get_reconstructed_kappa_alm(lensed_cmb_map, kappa_filter, kappa_norm,
                                lmax, mlmax, outdir,
                                lsst_mask_map=None, fg_map=None):
    time0 = time.time()
    readme = '\n*** Running get_reconstructed_kappa_alm function ***\n'
    file_tag = 'lensed_cmb_map'
    # -----------------------------------------------
    # see if need to incorporate the foregrounds map to the lensed map
    if fg_map is not None:
        # add the fg map to lensed cmb
        lensed_cmb_map += fg_map
        # update file tag
        file_tag += '_fg'
        # create the plot
        filename = plot_mollview(map_in=lensed_cmb_map,
                                 title=file_tag,
                                 data_label='',
                                 outdir=outdir,
                                 file_tag=file_tag,
                                 save_plot=True, show_plot=False)
        readme = print_update(update='Saved lensed cmb map with fg in %s\n'%(filename),
                              readme=readme)
    # -----------------------------------------------
    # see if need to incorporate the binary mask
    if lsst_mask_map is not None:
        # add the fg map to lensed cmb
        lensed_cmb_map *= lsst_mask_map
        # update file tag
        file_tag += '_xmask'
        # create the plot
        filename = plot_mollview(map_in=lensed_cmb_map,
                                 title=file_tag,
                                 data_label='',
                                 outdir=outdir,
                                 file_tag=file_tag,
                                 save_plot=True, show_plot=False)
        readme = print_update(update='Saved lensed cmb map with lsst mask in %s\n'%(filename),
                              readme=readme)
    # -----------------------------------------------
    # now convert to alms
    lensed_cmb_alms = hp.map2alm(lensed_cmb_map, lmax=lmax)
    # -----------------------------------------------
    # filter the alms
    readme = print_update(update='Filtering lensed_cmb_map alms ... ',
                          readme=readme)
    lensed_cmb_alms_filtered = hp.almxfl(alm=lensed_cmb_alms, fl=kappa_filter)
    # -----------------------------------------------
    # set up the geometry for full fsky reconstruction
    shape, wcs = enmap.fullsky_geometry(res=np.deg2rad(1.7/60.))
    # -----------------------------------------------
    # do the reconstruction
    readme = print_update(update='Reconstructing kappa ... ', readme=readme)
    kappa_alms, _ = qe.qe_tt(shape=shape, wcs=wcs,
                             Xalm=lensed_cmb_alms, Yalm=lensed_cmb_alms_filtered, do_curl=False,
                             mlmax=mlmax, lmax_x=lmax, lmax_y=lmax)
    # -----------------------------------------------
    # normalize with input normalization
    kappa_alms_normed = hp.almxfl(alm=kappa_alms, fl=kappa_norm)
    # update readme
    readme += '\nTime taken: %.3f min\n\n'%((time.time()-time0)/60.)

    return [kappa_alms_normed, readme]

# ------------------------------------------------------------------------------
def generate_correlated_alm(input_alm_f1, Clf1f1, Clf2f2, Clf1f2, seed=None):
    correlated = hp.almxfl(input_alm_f1, Clf1f2/Clf1f1)
    ps_noise = Clf2f2 - np.nan_to_num(Clf1f2**2/Clf1f1)

    assert np.all(ps_noise[2:] >= 0)

    if seed is not None: np.random.seed(seed)
    noise = hp.synalm(ps_noise, lmax=hp.Alm.getlmax(input_alm_f1.size))

    total = correlated + noise
    total[~np.isfinite(total)] = 0

    return total

# ------------------------------------------------------------------------------
# method from here: https://github.com/healpy/healpy/blob/master/healpy/rotator.py
def rotate_map(hp_rotator, m):
        """Rotate a HEALPix map to a new reference frame
        This function first rotates the pixels centers of the new reference
        frame to the original reference frame, then uses hp.get_interp_val
        to interpolate bilinearly the pixel values, finally fixes Q and U
        polarization by the modification to the psi angle caused by
        the Rotator using Rotator.angle_ref.
        Parameters
        ----------
        m : np.ndarray
            Input map, 1 map is considered I, 2 maps:[Q,U], 3 maps:[I,Q,U]
        Returns
        -------
        m_rotated : np.ndarray
            Map in the new reference frame
        """
        if pixelfunc.maptype(m) == 0:  # a single map is converted to a list
            m = [m]
        npix = len(m[0])
        nside = pixelfunc.npix2nside(npix)
        theta_pix_center, phi_pix_center = pixelfunc.pix2ang(nside=nside,
                                                             ipix=np.arange(npix)
                                                             )

        # Rotate the pixels center of the new reference frame to the original frame
        theta_pix_center_rot, phi_pix_center_rot = hp_rotator.I(theta_pix_center,
                                                                phi_pix_center
                                                                )

        # Interpolate the original map to the pixels centers in the new ref frame
        m_rotated = [pixelfunc.get_interp_val(each, theta_pix_center_rot,
                                              phi_pix_center_rot)
                     for each in m]

        # Rotate polarization
        if len(m_rotated) > 1:
            # Create a complex map from QU  and apply the rotation in psi due to the rotation
            # Slice from the end of the array so that it works both for QU and IQU
            L_map = (m_rotated[-2] + m_rotated[-1] * 1j) * \
                    np.exp(1j * 2 * hp_rotator.angle_ref(theta_pix_center,
                                                         phi_pix_center)
                           )

            # Overwrite the Q and U maps with the correct values
            m_rotated[-2] = np.real(L_map)
            m_rotated[-1] = np.imag(L_map)
        else:
            m_rotated = m_rotated[0]

        return m_rotated
