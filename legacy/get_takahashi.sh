#!/bin/bash

cd /global/cscratch1/sd/awan/soxlsst/takahashi_data

echo Getting wl maps
wget http://cosmo.phys.hirosaki-u.ac.jp/takahasi/allsky_raytracing/sub1/nres12/allskymap_nres12r000.zs66.mag.dat

echo Getting lensed CMB map
wget http://cosmo.phys.hirosaki-u.ac.jp/takahasi/allsky_raytracing/sub1/nres12/cmbmaps_betazero/lensed_cmbmap_betazero_nres12r000.fits

echo getting halo catalog
wget http://cosmo.phys.hirosaki-u.ac.jp/takahasi/allsky_raytracing/sub1/nres12/skyhalo_data/skyhalo_nres12r000.halo.gz
gunzip skyhalo_nres12r000.halo.gz