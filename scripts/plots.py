import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

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
    
def plot_cls_dict(cls):
    plt.clf()
    for key in cls:
        ell = np.arange(np.size(cls[key]))
        plt.plot(ell, (cls[key]*ell*(ell+1))/2.0/np.pi, 'o-', label= key)
    plt.xlabel(r'$\ell$', fontsize=18)
    plt.ylabel(r'$\ell(\ell+1)C_\ell/(2\pi)$', fontsize=16)
    plt.tick_params(axis='x', labelsize=16)
    plt.tick_params(axis='y', labelsize=16)
    #plt.gca().set_yscale('log')
    plt.gcf().set_size_inches(12, 6)
    plt.legend(fontsize=16)
    plt.show()