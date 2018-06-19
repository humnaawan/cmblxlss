import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

orphics_path= '/Users/Humna/repos/CAMB/orphics/'
__all__= ['plot_power_spectrum', 'plot_cls_dict']


def plot_power_spectrum(mapIn, lmax, returnCls= True):
    cl= hp.anafast(map1= mapIn, lmax= lmax)
    ell = np.arange(np.size(cl))
    plt.plot(ell, (cl*ell*(ell+1))/2.0/np.pi)
    plt.xlabel(r'$\ell$', fontsize=18)
    plt.ylabel(r'$\ell(\ell+1)C_\ell/(2\pi)$', fontsize=16)
    plt.tick_params(axis='x', labelsize=16)
    plt.tick_params(axis='y', labelsize=16)   
    plt.gca().set_yscale('log')
    plt.show()

    if returnCls: return [ell, cl]
    
def plot_cls_dict(cls, cross_convention= True,
                  binned= False, bin_width= 20, lmax= 1000):
    if binned:
        import sys
        sys.path.insert(0, orphics_path)
        from orphics.stats import bin1D
        
        wanted_bin_edges= np.arange(0, lmax, bin_width)
        binner1d= bin1D(wanted_bin_edges)
    
    plt.clf()
    for key in cls:
        if binned:
            ell_toplot, cl_toplot= binner1d.binned(np.arange(np.size(cls[key])), cls[key])
        else:
            ell_toplot, cl_toplot= np.arange(np.size(cls[key])), cls[key]
        # decide on the convention
        if cross_convention:
            plt.plot(ell_toplot, cl_toplot*ell_toplot, 'o-', label= key)
        else:
            plt.plot(ell_toplot, (cl_toplot*ell_toplot*(ell_toplot+1))/2.0/np.pi, 'o-', label= key)
    plt.xlabel(r'$\ell$', fontsize=18)
    if cross_convention:
        plt.ylabel(r'$\ell C_\ell$', fontsize=16)
    else:
        plt.ylabel(r'$\ell(\ell+1)C_\ell/(2\pi)$', fontsize=16)
    plt.tick_params(axis='x', labelsize=16)
    plt.tick_params(axis='y', labelsize=16)
    #plt.gca().set_yscale('log')
    if binned: plt.title('Binned: bin_width: %s'%bin_width, fontsize= 16)
    plt.gcf().set_size_inches(12, 6)
    plt.legend(fontsize=16)
    plt.show()