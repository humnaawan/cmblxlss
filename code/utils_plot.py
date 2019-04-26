import matplotlib as mpl
import matplotlib.pyplot as plt
import healpy as hp
import numpy as np
from matplotlib import cm
import numpy.ma as ma
from orphics.stats import bin1D

# ------------------------------------------------------------------------------
# imports from this folder
from settings import rcparams
# set up the plotting params
for key in rcparams: mpl.rcParams[key] = rcparams[key]
# ------------------------------------------------------------------------------

__all__= ['plot_power_spectrum', 'plot_cls_dict',
          'plot_seaborn_distplot', 'plot_seaborn_jointplot']

#--------------------------------------------------------------------------------------
def plot_power_spectrum(map_in, lmax, outdir, file_tag, title,
                        return_cls=True, save_plot=True, show_plot=False):
    # figure out the file tag
    if not (file_tag == '' or file_tag is None):
        file_tag = '_%s'%file_tag
    # calculate the cls
    cl = hp.anafast(map1=map_in, lmax=lmax)
    # set up the ells
    ell = np.arange(np.size(cl))
    # plot
    plt.clf()
    plt.plot(ell, (cl*ell*(ell+1))/2.0/np.pi)
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$\ell(\ell+1)C_\ell/(2\pi)$')
    plt.gca().set_yscale('log')
    plt.title(title)
    if save_plot:
        filename = 'power_spectrum%s.png'%file_tag
        plt.savefig('%s/%s'%(outdir, filename),
                    format='png', bbox_inches='tight')
    if show_plot: plt.show()
    else: plt.close('all')

    if return_cls:
        if save_plot: return [ell, cl, filename]
        else: return [ell, cl]
    else:
        return filename

#--------------------------------------------------------------------------------------
def plot_cls_dict(cls_in, outdir, file_tag, residuals=False, baseline_key=None,
                  save_plot=True, show_plot=False, loglog=False,
                  cross_convention=True, sci_yticks=True,
                  binned=False, bin_width=20, lmax=1000,
                  markers=None, colors=None, linestyles=None, lmin=None):
    #
    if residuals and baseline_key is None and not binned:
        raise ValueError('Need baseline_key + binned == True if want to plot residuals')
    nkeys = len(list(cls_in.keys()))
    # see if need to set up the binner
    if binned:
        if loglog:
            wanted_bin_edges = np.geomspace(10, lmax, bin_width)
        else:
            wanted_bin_edges = np.arange(0, lmax, bin_width)
        binner1d = bin1D(wanted_bin_edges)
    if residuals:
        _ , cl_baseline = binner1d.binned(np.arange(len(cls_in[baseline_key])), cls_in[baseline_key])
    # plot
    if markers is None:
        markers = ['P', 'x', '.', 'X', '+', '1', '3']
    if colors is None:
        colors = ['orangered', 'b', 'darkviolet', 'k', 'olive', 'darkred', 'c']
    if linestyles is None:
        linestyles = ['-']*nkeys
    ntot = len(cls_in.keys())
    plt.clf()
    for i, key in enumerate(cls_in):
        if binned:
            ell_toplot, cl_toplot = binner1d.binned(np.arange(len(cls_in[key])), cls_in[key])
        else:
            ell_toplot, cl_toplot = np.arange(np.size(cls_in[key])), cls_in[key]
        if residuals:
            cl_toplot = (cl_toplot - cl_baseline)/cl_baseline
        plt.plot(ell_toplot, cl_toplot, linestyles[i], label=key, color=colors[i%ntot])
        plt.scatter(ell_toplot, cl_toplot, marker=markers[i%ntot], color=colors[i%ntot])
    # plot details
    plt.xlabel(r'$\ell$')
    title = ''
    if residuals:
        plt.ylabel(r'($C_\ell$-$C_{\ell, 0}$)/$C_{\ell, 0}$')
        title = '$C_{\ell, 0}$ = %s'%baseline_key
    else:
        plt.ylabel(r'$C_\ell$')
        plt.gca().set_yscale('log')
    if loglog:
        plt.gca().set_xscale('log')
    if lmin is not None:
        xmin, xmax = plt.gca().get_xlim()
        plt.xlim(lmin, xmax)
    if binned:
        plt.title('Binned: bin_width: %s\n%s'%(bin_width, title))

    plt.legend(bbox_to_anchor=(1, 1))
    if save_plot:
        if not (file_tag == '' or file_tag is None):
            file_tag = '_%s'%file_tag
        if loglog:
            tag = '_loglog'
        else:
            tag = ''
        if residuals:
            tag += '_residuals'
        if lmin is not None:
            tag += '_lmin%s'%lmin
        if binned:
            filename = 'power_spectra_binned_%sbinwidth%s%s.png'%(bin_width, file_tag, tag)
        else:
            filename = 'power_spectra%s%s.png'%(file_tag, tag)
        plt.savefig('%s/%s'%(outdir, filename), format='png', bbox_inches='tight')
    if show_plot: plt.show()
    else: plt.close('all')

    if save_plot: return filename

#--------------------------------------------------------------------------------------
def plot_mollview(map_in, title, data_label, outdir, file_tag,
                  save_plot=True, show_plot=False):
    # set up the file tag
    if not (file_tag=='' or file_tag is None):
        file_tag = '_%s'%file_tag
    # check if input map is masked
    if not ma.is_masked(map_in): # construct a masked array
        map_in = map_in.view(ma.MaskedArray)
        map_in.mask = [False] * len(map_in.data)
        map_in.mask[abs(map_in.data) < 1e-5] = True
        map_in.fill_value = np.nan
    # get the nside for this map
    nside = hp.npix2nside(len(map_in.data))
    # set up the cmap and the background (dont way anything 'behind' the skymap)
    cmap = cm.viridis
    cmap.set_under('w')
    # set up the colorbar
    in_survey = np.where(map_in.mask == False)[0]
    median = np.nanmedian(map_in.data[in_survey])
    stddev = np.nanstd(map_in.data[in_survey])
    # figure out the ticks for the colorbar
    nticks = 5
    if stddev == 0:
        color_min, color_max = 0, 1
    else:
        color_min = median - 1.5 * stddev
        color_max = median + 1.5 * stddev
    increment = (color_max-color_min)/float(nticks)
    ticks = np.arange(color_min+increment, color_max, increment)
    # plot
    hp.mollview(map_in.filled(map_in.fill_value), flip='astro', rot=(0,0,0), cmap=cmap,
                min=color_min, max=color_max, title='', cbar=False)
    hp.graticule(dpar=20, dmer=20, verbose=False)
    plt.title('%s \nnside %s; min %.2f; max %.2f'%(title, nside,
                                                   min(map_in.data),
                                                   max(map_in.data)
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
            data_label = '_%s'%data_label
        filename = 'mollview%s%s_nside%s.png'%(data_label, file_tag, nside)
        plt.savefig('%s/%s'%(outdir, filename), format='png',
                    bbox_inches='tight', transparent=True)
    if show_plot: plt.show()
    else: plt.close('all')

    if save_plot: return filename
