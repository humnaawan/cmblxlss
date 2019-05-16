from __future__ import print_function
from orphics import maps,io,cosmology
from pixell import enmap
import numpy as np
import os,sys


lc = cosmology.LimberCosmology(skipCls=True)
lc.addStepNz("lsst",0.2,1.5,bias=1.6)
ellrange = np.arange(0,3000,1)
lc.generateCls(ellrange)
clgg = lc.getCl("lsst","lsst")
clkg = lc.getCl("cmb","lsst")

#io.save_cols('/global/cscratch1/sd/msyriac/so_lsst/gcls.txt',(ellrange,clgg,clkg))

theory = cosmology.loadTheorySpectraFromCAMB('/global/cscratch1/sd/msyriac/so_lsst/cosmo2017_10K_acc3') 
clkk = theory.gCl('kk',ellrange)

ps_noise = clgg - np.nan_to_num(clkg**2/clkk)
print(ps_noise)
assert np.all(ps_noise>=0)

