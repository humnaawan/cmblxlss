import numpy as np
import healpy as hp
import numpy.ma as ma

from utils_plot import plot_power_spectrum, plot_skymap

def add_fakelss(theory_density_map, fakelss_map_path,
                show_plot=False, save_plot=False, outdir=None):
    """

    This function the input theory galaxy density map and add fake LSS due to
    MW dust and LSST observing strategy.

    Required Inputs
    ---------------
    * theory_density_map: array: healpix map of the theory galaxy density map.
    * fakelss_map_path: str: path to the file with fakelss (npz file).

    Optional Inputs
    ---------------
    * show_plot: bool: set to True to display plots. Default: False
    * save_plot: bool: set to True to save plots. Default: False
    * outdir: str: path to the output directory; must be specified if save_plot
                   is True. Default: None.

    Returns
    -------
    * masked array containing the density map with fake-lss with mask
      corresponding to the fakelss map.

    """
    if save_plot and outdir is None:
        raise ValueError('outdir cannot be None if plot(s) are to be saved.')
    # load the fakelss map
    fakelss_data = np.load(fakelss_map_path)
    fakelss_map = fakelss_data['data'].view(ma.MaskedArray)
    fakelss_map.mask = fakelss_data['mask']
    fakelss_map.fill_value = fakelss_data['fill_value']

    print('fakelss_map ', fakelss_map)
    # get the nsides for the two maps
    nside_theory = hp.npix2nside(npix=len(theory_density_map))
    nside_fakelss = hp.npix2nside(npix=len(fakelss_map.data))
    # see if need to plot things
    if show_plot or save_plot:
        # -----
        # input theory
        title = 'input theory'
        file_tag = 'input-theory-density'
        map_in = theory_density_map
        nside = nside_theory
        plot_skymap(map_in=map_in, title=title, data_label='density',
                    show_plot=show_plot, save_plot=save_plot,
                    outdir=outdir, file_tag=file_tag
                    )
        plot_power_spectrum(map_in=map_in, lmax=int(nside*2),
                            title=title, return_cls=False,
                            show_plot=show_plot, save_plot=save_plot,
                            outdir=outdir, file_tag=file_tag,
                            )

        # read-in fakelss map
        title = 'read in fakelss'
        file_tag = 'readin-fakelss'
        map_in = fakelss_map
        nside = nside_fakelss
        plot_skymap(map_in=map_in, title=title, data_label='density',
                    show_plot=show_plot, save_plot=save_plot,
                    outdir=outdir, file_tag=file_tag
                    )
        plot_power_spectrum(map_in=map_in, lmax=int(nside*2),
                            title=title, return_cls=False,
                            show_plot=show_plot, save_plot=save_plot,
                            outdir=outdir, file_tag=file_tag,
                            )

    # see if the fakelss map resolution would work with the input
    if nside_fakelss != nside_theory:
        print('Must change fakelss map resolution from nside %s to nside %s' % (nside_fakelss, nside_theory))
        fakelss_map = hp.ud_grade(map_in=fakelss_map, nside_out=nside_theory)
        # see if need to plot things
        if show_plot or save_plot:
            # read-in fakelss map
            title = 'udgraded fakelss'
            file_tag = 'udgraded-fakelss',
            map_in = fakelss_map
            nside = nside_theory
            plot_skymap(map_in=map_in, title=title, data_label='density',
                        show_plot=show_plot, save_plot=save_plot,
                        outdir=outdir, file_tag=file_tag
                        )
            plot_power_spectrum(map_in=map_in, lmax=int(nside*2),
                                title=title, return_cls=False,
                                show_plot=show_plot, save_plot=save_plot,
                                outdir=outdir, file_tag=file_tag,
                                )

    # add fakelss
    theory_wfakelss = np.zeros(len(theory_density_map))
    in_survey = np.where( fakelss_map.mask == False )[0]
    theory_wfakelss[in_survey] = ( ( 1 + theory_density_map[in_survey] ) * (1 + fakelss_map.data[in_survey]) ) - 1

    # set up the masked array for the output
    theory_wfakelss = theory_wfakelss.view(ma.MaskedArray)
    theory_wfakelss.mask = fakelss_map.mask.copy()
    theory_wfakelss.fill_value = fakelss_map.fill_value

    # see if need to plot things
    if show_plot or save_plot:
        # -----
        # input theory
        title = 'theory+fakelss'
        file_tag = title
        map_in = theory_wfakelss.copy()
        nside = nside_theory
        plot_skymap(map_in=map_in, title=title, data_label='density',
                    show_plot=show_plot, save_plot=save_plot,
                    outdir=outdir, file_tag=file_tag
                    )
        plot_power_spectrum(map_in=map_in, lmax=int(nside*2),
                            title=title, return_cls=False,
                            show_plot=show_plot, save_plot=save_plot,
                            outdir=outdir, file_tag=file_tag,
                            )

    return theory_wfakelss
