from __future__ import print_function
from orphics import maps,io,cosmology,stats
from pixell import enmap,curvedsky as cs
import numpy as np
import os,sys
import healpy as hp
import mapsims
from collections import namedtuple
import healpy as hp
from falafel import qe
import symlens


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


def get_sim(i,band=145,cmb_dir="/home/msyriac/data/sims/alex/v0.4/",mlmax=3000,verify=False):
    res_3000 = 2.0
    shape,wcs = enmap.fullsky_geometry(res=np.deg2rad(res_3000*3000./mlmax/60.))
    with io.nostdout():
        simgen = mapsims.SOStandalonePrecomputedCMB(
            iteration_num=i,
            shape=shape,
            wcs=wcs,
            lensed=True,
            aberrated=False,
            has_polarization=False,
            cmb_set=0,
            input_units="uK_CMB",
            cmb_dir=cmb_dir)

        ch = namedtuple("Channel", ["telescope", "band"])
        ch.telescope = 'LA'
        ch.band = band
        imap = simgen.simulate(ch, output_units="uK_CMB").astype(np.float32)
        almphi = simgen.get_phi_alm()

    #This is unnecessarily large. Instead just trim phi to lmax=3000
    omap = enmap.zeros(shape[-2:],wcs)
    almphi = cs.map2alm(cs.alm2map(almphi,omap),lmax=mlmax)

    ells = np.arange(0,hp.Alm.getlmax(almphi.size))
    almk = hp.almxfl(almphi,ells*(ells+1.)/2.)

    if verify:
        alms = cs.map2alm(imap,lmax=mlmax)

        cls = hp.alm2cl(alms.astype(np.complex128))
        ells = np.arange(len(cls))
        lcltt = theory.lCl('TT',ells)
        lbeam = maps.gauss_beam(ells,1.4)
        pl = io.Plotter(xyscale='linlog',scalefn=lambda x: x**2./2./np.pi)
        pl.add(ells,cls,ls="--")
        pl.add(ells,lcltt*lbeam**2.)
        pl.done()
        pl = io.Plotter()
        pl.add(ells,lcltt*lbeam**2./cls,ls="--")
        pl.hline(y=1)
        pl.done()
        
    
    return imap,almk,shape,wcs
    

def recon(shape,wcs,dalms,norm,lcltt,nltt_deconvolved=None,lmin=300,lmax=3000,mlmax=None):
    """
    Lensing reconstruction. Possibly deprecated due to API change in falafel.
    """
    nltt_deconvolved = 0*lcltt if nltt_deconvolved is None else nltt_deconvolved
    mlmax = lmax if mlmax is None else mlmax
    ukappa,_ = qe.qe_tt_simple(Xalm=dalms,lcltt=lcltt,ucltt=lcltt, \
                               nltt_deconvolved=nltt_deconvolved,
                               lmin=lmin,lmax=lmax,shape=shape,wcs=wcs,mlmax=mlmax)
    ls = range(len(norm))
    ells = np.arange(0,hp.Alm.getlmax(ukappa.size))
    Als = maps.interp(ls,norm)(ells)
    almrecon = hp.almxfl(ukappa,Als/2.)
    return almrecon

