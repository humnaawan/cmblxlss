from __future__ import print_function
from orphics import maps,io,cosmology
from pixell import enmap,reproject
import numpy as np
import os,sys
import healpy as hp


"""
Reprojects an SO noise map from healpix to CAR to check its coordinate system.
This script was to re-discover the fact that the noise and Galaxy sims have different coordinate systems.
"""


#path = "/project/projectdirs/sobs/v4_sims/mbs/201903_highres_foregrounds/512/dust/0000/simonsobs_dust_uKCMB_la145_nside512_0000.fits"
path = "/project/projectdirs/sobs/v4_sims/mbs/201901_gaussian_fg_lensed_cmb_realistic_noise/4096/noise/0000/simonsobs_noise_uKCMB_la145_nside4096_0000.fits"
hmap = hp.ud_grade(hp.read_map(path),nside_out=512)
shape,wcs = enmap.fullsky_geometry(res=np.deg2rad(6./60.))
imap = reproject.enmap_from_healpix(hmap, shape, wcs, ncomp=1, unit=1, lmax=0,rot=None)
#io.hplot(imap,"reprojected_noise")
io.plot_img(imap,"reprojected_noise_0.png",lim=300)
    
