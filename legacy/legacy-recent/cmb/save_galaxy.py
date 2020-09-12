from __future__ import print_function
from orphics import maps,io,cosmology
from pixell import enmap
import numpy as np
import os,sys
import healpy as hp

"""
This script saves the sum of all Galactic foregrounds at 145 GHz.
"""

components = ['ame','dust','freefree','synchrotron']

tmap = 0
for comp in components:
    imap = hp.read_map('/project/projectdirs/sobs/v4_sims/mbs/201903_highres_foregrounds/4096/%s/0000/simonsobs_%s_uKCMB_la145_nside4096_0000.fits' % (comp,comp))

    tmap = tmap + hp.ud_grade(imap,nside_out = 2048)
    print(comp)

hp.write_map('/global/cscratch1/sd/msyriac/so_lsst/galaxy_145.fits',tmap,overwrite=True)
