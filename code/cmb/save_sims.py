from __future__ import print_function
from orphics import maps,io,cosmology
from pixell import enmap
import numpy as np
import os,sys
import mapsims
import healpy as hp

"""
This script will generate and save sims of the SO sky that include
lensing and Galactic foregrounds.
"""


# Lensed sky
simgen = mapsims.SOStandalonePrecomputedCMB(0, nside=2048, lensed=True, aberrated=False, has_polarization=False,input_units='uK_CMB')
#hmap = simgen.simulate(ch=mapsims.Channel(telescope='LA',band=145))[0]
#io.mollview(hmap,'hmap.png')
#hp.write_map("/global/cscratch1/sd/msyriac/so_lsst/lensed_cmb_nside_2048_uK.fits",hmap,overwrite=True)
phi_alm = simgen.get_phi_alm()
hp.write_alm("/global/cscratch1/sd/msyriac/so_lsst/phi_alm.fits",phi_alm,overwrite=True)

