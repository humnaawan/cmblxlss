import matplotlib.pyplot as plt
import healpy as hp
import numpy as np
from matplotlib import cm
import numpy.ma as ma
import copy

# ------------------------------------------------------------------------------
# plot style imports
from settings import *
# ------------------------------------------------------------------------------

__all__= ['plot_power_spectrum', 'plot_skymap']

# ------------------------------------------------------------------------------
def plot_power_spectrum(map_in, lmax, title,
                        return_cls=False, show_plot=True,
                        save_plot=False, outdir=None, file_tag=None):
    """

    This function plots the power spectrum of the input HEALPix smap.

    Required Inputs
    ---------------
    * map_in: array: array containing the HEALPix map.
    * lmax: int/float: maximum ell to consider, e.g., 2 * nside.
    * title: str: title for the plot.

    Optional Inputs
    ---------------
    * return_cls: bool: set to True to get the cells. Default: False
    * show_plot: bool: set to False to not show plot. Default: True
    * save_plot: bool: set to True to save plot. Default: False
    * outdir: str: path to the output directory; must be specified if save_plot
                   is True. Default: None.
    * file_tag: str: tag to include in the output filename. Default: None.

    Returns
    -------
    * [ells, c_ells] IF return_cls is True.

    """
    # figure out the file tag
    if save_plot:
        if outdir is None:
            raise ValueError('outdir cannot be None if plot(s) are to be saved.')
        if file_tag is None:
            file_tag = ''
        if not (file_tag == ''):
            file_tag = '_%s' % file_tag
    # calculate the cls
    cls = hp.anafast(map1=map_in, lmax=lmax)
    # set up the ells
    ells = np.arange(np.size(cls))
    # plot
    plt.clf()
    plt.plot(ells, (cls * ells * (ells+1)) / 2.0 / np.pi)
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$\ell(\ell+1)C_\ell/(2\pi)$')
    plt.gca().set_yscale('log')
    plt.title(title)
    if save_plot:
        filename = 'plot_power-spec%s.png' % file_tag
        plt.savefig('%s/%s' % (outdir, filename),
                    format='png', bbox_inches='tight')
        print('## saved %s' % filename)
    if show_plot: plt.show()
    else: plt.close('all')

    if return_cls:
        return [ells, c_ells]

# ------------------------------------------------------------------------------
def plot_skymap(map_in, title, data_label,
                show_plot=True, save_plot=False, outdir=None, file_tag=None):
    """

    This function plots the skymap of the input map.

    Required Inputs
    ---------------
    * map_in: array: array containing the HEALPix map.
    * title: str: title for the plot.
    * data_label: str: label to designate the contet of that map (e.g. density).

    Optional Inputs
    ---------------
    * show_plot: bool: set to False to not show plot. Default: True
    * save_plot: bool: set to True to save plot. Default: False
    * outdir: str: path to the output directory; must be specified if save_plot
                   is True. Default: None.
    * file_tag: str: tag to include in the output filename. Default: None.

    """
    # figure out the file tag
    if save_plot:
        if outdir is None:
            raise ValueError('outdir cannot be None if plot(s) are to be saved.')
        if file_tag is None:
            file_tag = ''
        if not (file_tag == ''):
            file_tag = '_%s' % file_tag
    # check if input map is masked
    if ma.is_masked(map_in):
        # get the nside for this map
        nside = hp.npix2nside(len(map_in.data))
        # set up the colorbar
        in_survey = np.where(map_in.mask == False)[0]
        if len(in_survey) == 0:
            raise ValueError('everything is masked?!?!?!')
        data_to_consider = map_in.data[in_survey]
    else:
        nside = hp.npix2nside(len(map_in))
        data_to_consider = map_in

    if len(np.unique(data_to_consider)) == 2:
        # i.e. a binary map
        stddev = 0
    else:
        # set up the colorbar
        median = np.nanmedian(data_to_consider)
        stddev = np.nanstd(data_to_consider)

    cmap = copy.copy(mpl.cm.get_cmap("viridis"))
    # figure out the ticks for the colorbar
    nticks = 5
    if stddev == 0:
        color_min, color_max = min(data_to_consider), max(data_to_consider)
    else:
        color_min = median - 1.5 * stddev
        color_max = median + 1.5 * stddev

    increment = (color_max - color_min) / float(nticks)
    ticks = np.arange(color_min + increment, color_max, increment)

    # plot
    if ma.is_masked(map_in):
        to_plot = map_in.filled(map_in.fill_value)
    else:
        to_plot = map_in

    hp.mollview(to_plot, flip='astro', rot=(0,0,0), cmap=cmap,
                min=color_min, max=color_max, title='', cbar=False)
    hp.graticule(dpar=20, dmer=20, verbose=False)
    plt.title('%s \nnside %s; min %.2e; max %.2e' % (title, nside,
                                                     min(data_to_consider),
                                                     max(data_to_consider)
                                                     ))
    im = plt.gca().get_images()[0]
    fig = plt.gcf()
    cbaxes = fig.add_axes([0.1, 0.03, 0.8, 0.04]) # [left, bottom, width, height]
    cb = plt.colorbar(im,  orientation='horizontal',
                      ticks=ticks, format='%.4f', cax=cbaxes)
    cb.set_label(data_label)
    # save fig
    if save_plot:
        if not (data_label == '' or data_label is None):
            data_label = '_%s' % data_label
        filename = 'plot_skymap%s%s_nside%s.png' % (data_label, file_tag, nside)
        plt.savefig('%s/%s' % (outdir, filename), format='png',
                    bbox_inches='tight', transparent=True)
        print('## saved %s' % filename)
    if show_plot: plt.show()
    else: plt.close('all')
