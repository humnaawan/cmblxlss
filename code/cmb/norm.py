from __future__ import print_function
from orphics import maps,io,cosmology,stats
from pixell import enmap
import numpy as np
import os,sys
import symlens

"""
This script saves a symlens normalization of the lensing map
"""

def get_norm(uctt,tctt,lmin=300,lmax=3000,plot=False):
    shape,wcs = maps.rect_geometry(width_deg=80.,px_res_arcmin=2.0*3000./lmax)
    emin = maps.minimum_ell(shape,wcs)
    modlmap = enmap.modlmap(shape,wcs)
    tctt = maps.interp(range(len(tctt)),tctt)(modlmap)
    uctt = maps.interp(range(len(uctt)),uctt)(modlmap)
    tmask = maps.mask_kspace(shape,wcs,lmin=lmin,lmax=lmax)
    Al = symlens.A_l(shape, wcs, feed_dict={'uC_T_T':uctt,'tC_T_T':tctt}, estimator="hdv", XY="TT", xmask=tmask, ymask=tmask)
    bin_edges = np.arange(3*emin,lmax,2*emin)
    binner = stats.bin2D(modlmap,bin_edges)
    cents,Al1d = binner.bin(Al)
    ls = np.arange(0,cents.max(),1)
    Als = np.interp(ls,cents,Al1d*cents**2.)/ls**2.
    Als[ls<1] = 0
    if plot:
        pl = io.Plotter(xyscale='loglog')
        pl.add(cents,Al1d*cents**2.)
        pl.add(ls,Als*ls**2.,ls="--")
        pl.done()
    return Als

theory = cosmology.loadTheorySpectraFromCAMB('/global/cscratch1/sd/msyriac/so_lsst/cosmo2017_10K_acc3') 
ells = np.arange(0,4000)
uctt = theory.lCl('TT',ells)
tctt = theory.lCl('TT',ells)
Als = get_norm(uctt,tctt)
print(Als)

io.save_cols('/global/cscratch1/sd/msyriac/so_lsst/norm.txt',(np.arange(Als.size),Als))
