from __future__ import print_function
from orphics import maps,io,cosmology
from pixell import enmap,reproject
import numpy as np
import os,sys
import healpy as hp

"""
This script explores the SO hit maps
"""

#sname = "total_hits_LA_opportunistic"
sname = "total_hits_LA_classical"
hmap = hp.read_map(os.environ['HOME']+"/repos/mapsims/mapsims/data/%s.fits.gz" % sname)
io.mollview(hmap,sname+".png")

#shape,wcs = enmap.fullsky_geometry(res=np.deg2rad(0.5/60))
shape,wcs = enmap.band_geometry(np.deg2rad((-80,40)),res=np.deg2rad(0.5/60))
#shape,wcs = enmap.band_geometry(np.deg2rad((40,-80)),res=np.deg2rad(0.5/60))




scale = 8
oshape,owcs = enmap.scale_geometry(shape, wcs, 1./scale)
imap = reproject.enmap_from_healpix_interp(hmap, oshape, owcs, interpolate=False,rot=None)
imap[imap<=hmap[hmap!=0].min()] = 0

hmap /= hmap.max()
imap /= imap.max()
print((hmap != 0).sum()/len(hmap))
print((hmap != 0).sum()/hmap.size)
print((imap != 0).sum()/imap.size)
sys.exit()

fname = "reprojected_%s_act_scaled_%d" % (sname,scale)
enmap.write_map(fname+".fits",imap)
io.hplot(imap,fname,min=0,max=350)
io.plot_img(imap,fname+"_lowres.png")
print(hmap.min(),hmap.max())
print(imap.min(),imap.max())
sys.exit()
scales = [1,4,16]

for scale in scales:
    oshape,owcs = enmap.scale_geometry(shape, wcs, 1./scale)
    imap2 = enmap.project(imap, oshape, owcs, order=0)
    fname = "reprojected_%s_act_projected_%d" % (sname,scale)
    enmap.write_map(fname+".fits",imap2)
    io.hplot(imap2,fname)
    io.plot_img(imap2,fname+"_lowres.png")
    print(imap2.min(),imap2.max())

    

