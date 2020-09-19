# the goal here is run the add-fakelss code with multiple realizations of the
# theory density map; this will allow getting a handle on the uncertainties
# in the power spectrum and checking whether the final spectrum outputted
# from the code is correct.

import numpy as np
import healpy as hp
import os
import time

import sys
sys.path.append('../code/')

from settings import *
from add_fakelss import add_fakelss
from utils_plot import plot_skymap

time0 = time.time()
# set up the outdir
outdir = '../outputs/'
os.makedirs(outdir, exist_ok=True)

nside = 1024
lmax = int(nside * 2)

# path for the fakelss map
fakelss_map_path = '../data/fakelss_baseline-v1.4_nside1024.npz'

# set up the theory galaxy specturm
gcls = np.genfromtxt('../data/gcls.txt')
# extract the ells
ells = gcls[:, 0]
# extract the cells (i think this is right)
cls = gcls[:, 1]

# set up for the power spectra calculation
# lets read in the fakelss map
fakelss = np.load(fakelss_map_path)
cls_fakelss = hp.anafast(fakelss['data'], lmax=lmax)
ells = np.arange(np.size(cls_fakelss))
fsky = len( np.where( fakelss['mask'] == False )[0] ) / len(fakelss['mask'])

#alpha = 0.25
#plt.clf()
# different realizations
for i, seed in enumerate( range(100) ):
    np.random.seed(seed)
    # now realize a density map using the spectrum
    theory_map = hp.synfast(cls, nside=nside)
    # now run the function
    new_map = add_fakelss(theory_density_map=theory_map,
                          fakelss_map_path=fakelss_map_path,
                          show_plot=False, save_plot=False)

    cls_theory = hp.anafast(theory_map, lmax=lmax)
    cls_newmap = hp.anafast(new_map, lmax=lmax)

    if i == 0:
        cls_theory_stack = cls_theory.copy()
        cls_newmap_stack = cls_newmap.copy()
    else:
        cls_theory_stack = np.vstack( [cls_theory_stack, cls_theory.copy()] )
        cls_newmap_stack = np.vstack( [cls_newmap_stack, cls_newmap.copy()] )

    if False:
        # add each curve to plot
        color = '#1f77b4'
        plt.plot(ells, (cls_theory * ells * (ells+1)) / 2.0 / np.pi,
                 color=color, alpha=alpha)
        color = 'k'
        plt.plot(ells, (cls_newmap * ells * (ells+1)) / 2.0 / np.pi,
                 color=color, alpha=alpha)
        color = '#d62728'
        plt.plot(ells, ((cls_theory + cls_fakelss/fsky ) * ells * (ells+1)) / 2.0 / np.pi,
                 color=color, alpha=alpha)

plt.clf()
alpha = 0.25
# lets create a log-log plot to see the scales a little more clearly
# theory spectrum
color = '#1f77b4'
mean_theory = np.mean(cls_theory_stack, axis=0)
stddev_theory = np.std(cls_theory_stack, axis=0)
x, y = ells, (mean_theory * ells * (ells+1)) / 2.0 / np.pi
yerr = (stddev_theory * ells * (ells+1)) / 2.0 / np.pi
# now plot
plt.plot(x, y, label='theory full sky')
plt.fill_between(x, y, y + yerr, alpha=alpha, facecolor=color)
plt.fill_between(x, y, y - yerr, alpha=alpha, facecolor=color)

# fakelss spectrum
color = '#ff7f0e'
plt.plot(ells, ( (cls_fakelss / fsky) * ells * (ells+1)) / 2.0 / np.pi,
         label='fakelss / fsky', color=color)

# new spectrum
color = 'k'
mean = np.mean(cls_newmap_stack, axis=0)
stddev = np.std(cls_newmap_stack, axis=0)
x, y = ells, ((mean / fsky) * ells * (ells+1)) / 2.0 / np.pi
yerr = ( (stddev/fsky) * ells * (ells+1)) / 2.0 / np.pi
# now plot
plt.plot(x, y, color=color, label='theory+fakelss / fsky (from code)')
plt.fill_between(x, y, y + yerr, alpha=alpha, facecolor=color)
plt.fill_between(x, y, y - yerr, alpha=alpha, facecolor=color)

# plot details
plt.legend(bbox_to_anchor=(1,1))
plt.xlabel(r'$\ell$')
plt.ylabel(r'$\ell(\ell+1)C_\ell/(2\pi)$')
plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
filename = 'plot_theory+fakelss_%srealizations_nside%s.png' % (i, nside)
plt.savefig('%s/%s' % (outdir, filename), format='png',
            bbox_inches='tight', transparent=True)

# add  by hand calculation as a check; this should match the "new curve" above
color = '#d62728'
x, y = ells, ((mean_theory + cls_fakelss/fsky ) * ells * (ells+1)) / 2.0 / np.pi
# now plot
plt.plot(x, y, color=color, label='theory+fakelss / fsky (by hand)')

# plot details
plt.legend(bbox_to_anchor=(1,1))
plt.xlabel(r'$\ell$')
plt.ylabel(r'$\ell(\ell+1)C_\ell/(2\pi)$')
plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
filename = 'plot-check_theory+fakelss_%srealizations_nside%s.png' % (i+1, nside)
plt.savefig('%s/%s' % (outdir, filename), format='png',
            bbox_inches='tight', transparent=True)
plt.close('all')

print('## Time taken: %.2f min' % ( (time.time() - time0 ) / 60.) )
