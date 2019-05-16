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
from so_lsst_gal import get_norm,get_sim,recon


"""
This script uses functions in so_lsst_gal.py to get CMB sims, reconstruct lensing
and verify reconstruction. It might not work anymore due to API changes in falafel.
"""

# Location where the lensed sims are
cmb_dir = "/home/msyriac/data/sims/alex/v0.4/"

nsims = 10

# CMB multipoles to include in reconstruction
lmin = 300
lmax = 3000

# Beam of simulation. This is taken to be exactly the value
# in the 145 band of SO, which is what is being simulated here.

fwhm = 1.4

# Load true CMB theory
theory = cosmology.loadTheorySpectraFromCAMB("cosmo2017_10K_acc3",lpad=9000,get_dimensionless=False)

# stats collector
s = stats.Stats()

for i in range(nsims):
    
    # Get a simulated map, alms of kappa and the geometry specification shape,wcs of the sim
    imap,almk,shape,wcs = get_sim(i,cmb_dir=cmb_dir)
    # Get alms of map
    alms = cs.map2alm(imap,lmax=lmax)

    if i==0:
        # Get beam as function of ells
        ells = np.arange(0,hp.Alm.getlmax(alms.size)+1)
        lbeam = maps.gauss_beam(ells,fwhm)
        # Get normalization for estimator
        lcltt = theory.lCl('TT',ells)
        norm = get_norm(lcltt,lcltt,lmin,lmax,plot=False)
        
    # Deconvolve beam
    dalms = hp.almxfl(alms,1/lbeam)
    # Get reconstruction
    almrecon = recon(shape,wcs,dalms,norm,lcltt,nltt_deconvolved=None,lmin=lmin,lmax=lmax,mlmax=lmax)

    # Get input kappa autopower
    cl_ii = hp.alm2cl(almk,almk)
    # Get recon x input power
    cl_ir = hp.alm2cl(almk,almrecon)
    # Add these to stats collector
    s.add_to_stats('cl_ii',cl_ii)
    s.add_to_stats('cl_ir',cl_ir)
    print(i)

s.get_stats()

# Get means
cl_ii = s.stats['cl_ii']['mean']
cl_ir = s.stats['cl_ir']['mean']

# Plot
pl = io.Plotter(xyscale='loglog')
pl.add(ells,theory.gCl('kk',ells),color='k',lw=3,alpha=0.5)
pl.add(ells,cl_ii,color='k')
pl.add(ells,cl_ir)
pl._ax.set_ylim(1e-9,1e-6)
pl.done("clkk2.png")

pl = io.Plotter(xyscale='loglin')
pl.add(ells,cl_ir/cl_ii)
pl.hline(y=1)
pl._ax.set_ylim(0.5,1.5)
pl.done("clkk3.png")

pl = io.Plotter(xyscale='linlin')
pl.add(ells,cl_ir/cl_ii)
pl.hline(y=1)
pl._ax.set_ylim(0.8,1.2)
pl.done("clkk4.png")

