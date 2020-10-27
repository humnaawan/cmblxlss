import numpy as np
import healpy as hp
import scipy.integrate
import pickle
import os
from settings import *

__all__ = ['get_dndz', 'generate_correlated_alm']

# ------------------------------------------------------------------------------
def get_dndz(outdir, zedges, zmin_overall, zmax_overall, ilim, redo_calc=False, show_plot=False):
    """

    This function generates the dndz, plots it and saves it.

    Required Inputs
    ---------------
    * outdir: str: path to the directory where the subdirectory will be created.
    * zedges: list: list of redshift bin edges, e.g., [0.1, 0.3, 0.5] for a
                    2-bin analysis.
    * zmin_overall: float: minimum redshift over which to evaluate the N(z).
    * zmax_overall: float: maximum redshift over which to evaluate the N(z).
    * ilim: float: limiting i-band magnitude.

    Optional Inputs
    ---------------
    * redo_calc: bool: set to True to not read the results from disk, even if
                       if exists. Default: False
    * show_plot: bool: set to True to show the plot. Default: False.

    """
    sigma_z0 = 0.05
    z_bias = 0

    # set up
    print('\n## running get_dndz ... ')
    nzbins = len(zedges)
    subdir = '%s/dndz_%sbins/' % (outdir, nzbins)
    os.makedirs(subdir, exist_ok=True)
    print('## will be saving stuff in %s ' % subdir)
    # ------------------------------------------------------
    fname = 'dndz_%sbins.pickle' % (nzbins)
    file_path = '%s/%s' % (subdir, fname)
    if os.path.exists(file_path) and not redo_calc:
        print('## reading in the already-saved dndz from %s\n' % fname)
        with open(file_path, 'rb') as handle:
            dict = pickle.load(file=handle)
        z_arr = dict['z_arr']
        dndz = dict['dndz']
    else:
        print('## calculating dndz .. ')
        # ------------------------------------------------------
        # true n(z) distribution
        def nz_unnormed(z, alpha=2, beta=1.5, z0=0.5):
            # based on eq 13.10 in lsst science book
            return  z**alpha * np.exp( -(z/z0)**beta )

        # normalize the n(z) distribution
        def nz_normed(z_arr, i=25.3, minz=0, maxz=10):
            # normalization by eq 3.7: 46 \times 3600 \times 10^{0.31(i-25)} \mathrm{galaxies / deg^2} \\
            ngal_deg2 = 46 * 3600 * 10**(0.31*(i-25))
            val, _ = scipy.integrate.quad(func=nz_unnormed, a=minz, b=maxz)
            return  (ngal_deg2/val) * nz_unnormed(z=z_arr)

        # ------------------------------------------------------
        # photo-z distribution
        def photoz_dist_unnormed(z_phot, z, sigma_z0, z_bias # redshift bias
                                ):
            # based on eq 13.9 from the lsst science book
            sigma_z = sigma_z0 * (1 + z)
            if z_phot < 0:
                return 0
            else:
                return np.exp( -(z_phot - z - z_bias)**2 / (2 * sigma_z**2)) / (np.sqrt(2 * np.pi) * sigma_z)

        # ------------------------------------------------------
        z_arr = np.linspace(zmin_overall, zmax_overall, 1000)
        # initiate array for dndz
        nz = nz_normed(z_arr=z_arr, i=ilim)

        dndz = []
        for zi in range(nzbins):
            # double integral to take care of the normalization
            denominator, _ = scipy.integrate.dblquad(func=photoz_dist_unnormed,
                                                     a=0, b=10,
                                                     gfun=lambda x: zedges[zi][0],
                                                     hfun=lambda x: zedges[zi][1],
                                                     args=(sigma_z0, z_bias)
                                                    )

            pz_normed = []
            for z in z_arr:
                numerator, _ = scipy.integrate.quad(func=photoz_dist_unnormed,
                                                    a=zedges[zi][0], b=zedges[zi][1],
                                                    args=(z, sigma_z0, z_bias))
                pz_normed.append(numerator / denominator)
            pz_normed = np.array(pz_normed)
            dndz.append( nz * pz_normed )

        # --------------------------------
        # plot the N(z) for each bin
        summed = 0
        plt.clf()
        for i in range(nzbins):
            plt.plot(z_arr, dndz[i], '.-', label='bin %s' % (i+1))
            summed += dndz[i]
        plt.plot(z_arr, summed, 'k--')
        # plot the zbin edges
        ylims = plt.gca().get_ylim()
        for i, zj in enumerate( np.unique(zedges.flatten()) ):
            plt.plot( [ zj, zj ] , ylims, '-' , color='#7f7f7f', lw=2)
        #
        plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        plt.legend()
        plt.xlabel('z')
        plt.ylabel(r'$n(z) = d/dz (dN/d\Omega)$')
        fname = 'plot_dndz_%sbinsz.png' % (nzbins)
        plt.savefig('%s/%s' % (subdir, fname), format='png', bbox_inches='tight')
        print('## saved %s' % fname)
        if show_plot:
            plt.show()
        else:
            plt.close()
        # --------------------------------
        # save the dndz
        dict = {'z_arr': z_arr, 'dndz': dndz}
        with open(file_path, 'wb') as f_out:
            pickle.dump(dict, f_out)
        print('## saved dndz in %s\n' % fname)

    return z_arr, dndz

# ------------------------------------------------------------------------------
# from https://github.com/ACTCollaboration/xgal/blob/master/xgal/deltag.py#L11
def generate_correlated_alm(input_alm_f1,Clf1f1,Clf2f2,Clf1f2,seed=None):
    correlated = hp.almxfl(input_alm_f1, Clf1f2/Clf1f1)
    ps_noise = Clf2f2 - np.nan_to_num(Clf1f2**2/Clf1f1)

    assert np.all(ps_noise >= 0)

    if seed is not None: np.random.seed(seed)
    noise = hp.synalm(ps_noise, lmax=hp.Alm.getlmax(input_alm_f1.size))

    total = correlated + noise
    total[~np.isfinite(total)] = 0

    return total