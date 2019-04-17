import matplotlib as mpl
import matplotlib.pyplot as plt
import healpy as hp
import numpy as np
from matplotlib import cm
import seaborn as sns
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
	if not (file_tag == '' or file_tag is None):
		file_tag = '_%s'%file_tag
	cl = hp.anafast(map1=map_in, lmax=lmax)
	ell = np.arange(np.size(cl))
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
def plot_cls_dict(cls_in, outdir, file_tag,
                  save_plot=True, show_plot=False,
                  cross_convention=True, sci_yticks=True,
                  binned=False, bin_width=20, lmax=1000):
    if binned:        
        wanted_bin_edges = np.arange(0, lmax, bin_width)
        binner1d = bin1D(wanted_bin_edges)
        
    markers = ['d', 'x', '.', '+', ',', '*']
    plt.clf()
    for i, key in enumerate(cls_in):
        if binned:
            ell_toplot, cl_toplot = binner1d.binned(np.arange(len(cls_in[key])), cls_in[key])
        else:
            ell_toplot, cl_toplot = np.arange(np.size(cls_in[key])), cls_in[key]
		# decide on the convention
        if cross_convention:
            plt.plot(ell_toplot, cl_toplot * ell_toplot, label=key)
            plt.scatter(ell_toplot, cl_toplot * ell_toplot,
                        marker=markers[i])
        else:
            plt.plot(ell_toplot, (cl_toplot * ell_toplot * (ell_toplot+1))/2.0/np.pi,
                     'o-', label=key)
    plt.xlabel(r'$\ell$')
    if cross_convention:
        plt.ylabel(r'$\ell C_\ell$')
    else:
        plt.ylabel(r'$\ell(\ell+1)C_\ell/(2\pi)$')

    if sci_yticks: plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.xlim(-1, lmax+bin_width)
    if binned: plt.title('Binned: bin_width: %s'%bin_width)
    plt.legend(bbox_to_anchor=(0.5, 1))
    if save_plot:
        if not (file_tag=='' or file_tag is None):
            file_tag = '_%s'%file_tag
        if binned:
            filename = 'power_spectra_binned_%sbinwidth%s.png'%(bin_width, file_tag)
        else:
            filename = 'power_spectra%s.png'%file_tag
        plt.savefig('%s/%s'%(outdir, filename), format='png', bbox_inches='tight')
    if show_plot: plt.show()
    else: plt.close('all')
        
    if save_plot: return filename

#--------------------------------------------------------------------------------------
def plot_seaborn_distplot(data_dict, keys_to_plot, outdir, file_tag,
                          save_plot=True, show_plot=False):
    filenames = []
    for key in keys_to_plot:
        if key not in data_dict.keys():
            raise ValueError('%s not in the input dictionary.'%key)
        plt.clf()
        sns.distplot(data_dict[key])
        plt.xlabel(key)
        plt.ylabel('Counts')
        if save_plot:
            filename='diagnostic_distplot_%s_%s.png'%(key, file_tag)
            plt.savefig('%s/%s'%(outdir, filename), format='png', bbox_inches='tight')
            filenames.append(filename)
        if show_plot: plt.show()
        else: plt.close('all')
    if save_plot: return filenames

#--------------------------------------------------------------------------------------
def plot_seaborn_jointplot(x_in, y_in, outdir, file_tag,
                           save_plot=True, show_plot=False):
	plt.clf()
	sns.jointplot(x=x_in, y=y_in, kind='hex')
	if save_plot:
		filename='diagnostic_jointplot_%s.png'%(file_tag)
		plt.savefig('%s/%s'%(outdir, filename), format='png', bbox_inches='tight')
	if show_plot: plt.show()
	else: plt.close('all')

	if save_plot: return filename

#--------------------------------------------------------------------------------------
def plot_mollview(map_in, title, data_label, outdir, file_tag,
                  save_plot=True, show_plot=False):
	#
	if not (file_tag=='' or file_tag is None):
		file_tag = '_%s'%file_tag
	# check if input map is masked
	if not ma.is_masked(map_in): # construct a masked array
		map_in = map_in.view(ma.MaskedArray)
		map_in.mask = [False] * len(map_in.data)

	nside = hp.npix2nside(len(map_in.data))

	# set up the colorbar
	in_survey = np.where(map_in.mask == False)[0]
	median = np.median(map_in.data[in_survey])
	stddev = np.std(map_in.data[in_survey])

	nticks= 5
	color_min = median - 1.5 * stddev
	color_max = median + 1.5 * stddev
	increment = (color_max-color_min)/float(nticks)
	ticks = np.arange(color_min+increment, color_max, increment)

	# plot
	hp.mollview(map_in, flip='astro', rot=(0,0,0), cmap=cm.viridis,
                min=color_min, max=color_max, title='', cbar=False)
	hp.graticule(dpar=20, dmer=20, verbose=False)
	plt.title('%s ; nside %s'%(title, hp.npix2nside(len(map_in.data))))
	im = plt.gca().get_images()[0]
	fig = plt.gcf()
	cbaxes = fig.add_axes([0.1, 0.03, 0.8, 0.04]) # [left, bottom, width, height]
	cb = plt.colorbar(im,  orientation='horizontal',
                      ticks=ticks, format='%.2f', cax=cbaxes)
	cb.set_label(data_label)

	# save fig
	if save_plot:
		if not (data_label=='' or data_label is None):
			data_label = '_%s'%data_label
		filename = 'mollview%s%s.png'%(data_label, file_tag)
		plt.savefig('%s/%s'%(outdir, filename), format='png',
                    bbox_inches='tight', transparent=True)
	if show_plot: plt.show()
	else: plt.close('all')

	if save_plot: return filename
