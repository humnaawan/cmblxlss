#!/bin/bash
################################################################################
# This script runs two python scripts: the first calculates and saves alms for
# kappa and g with different variations, while the second reads in the saved
# alms and calculates spectra.
#
# The script allows re-using saved alms for different spectra without having to
# run time-consuming calculations repeatedly.
################################################################################
# source the stack
source /global/common/software/lsst/cori-haswell-gcc/stack/setup_current_sims.sh
# ------------------------------------------------------------------------------
save_alms=1  # set to 0 to just use the saved alms and calculate spectra.
# ------------------------------------------------------------------------------
# set up inputs
lmax=3000
repo_path=/global/homes/a/awan/so/cmblxlss/
outdir=${repo_path}/outputs

# cmb related inputs
lensed_cmb_map_path=/global/cscratch1/sd/msyriac/so_lsst/lensed_cmb_nside_2048_uK.fits
cosmo_path=/global/cscratch1/sd/msyriac/so_lsst/cosmo2017_10K_acc3
kappa_norm_path=/global/cscratch1/sd/msyriac/so_lsst/norm.txt
kappa_alm_theory_path=/global/cscratch1/sd/msyriac/so_lsst/phi_alm.fits
gcls_path=/global/cscratch1/sd/msyriac/so_lsst/gcls.txt
fg_map_path=/global/cscratch1/sd/msyriac/so_lsst/galaxy_145.fits

# lsst related inputs
lsst_path='/global/cscratch1/sd/awan/soxlsst/gal_maps/'
lsst_path=${lsst_path}'fakelss_baseline2018a_Y10_nside256_i<25.3_allz_wDust_noPoisson_no0pt_egfootprint_normedNgal_directory/'
lsst_path=${lsst_path}'dNbyN_data/dNbyN_data_i_RandomDitherPerNight.npz'
lsst_data_tag=baseline2018a-allz-eg

# ------------------------------------------------------------------------------
if [ $save_alms == 1 ];
then
# run the code to save alms
python ${repo_path}/code/save_needed_alms.py \
                    --outdir=${outdir} \
                    --lensed-cmb-path=${lensed_cmb_map_path} \
                    --cosmology-path=${cosmo_path} \
                    --kappa-norm-path=${kappa_norm_path} \
                    --kappa-alm-theory_path=${kappa_alm_theory_path} \
                    --gcals-path=${gcls_path} \
                    --fg-map-path=${fg_map_path} \
                    --lsst-path=${lsst_path} \
                    --lsstdata-tag=${lsst_data_tag} \
                    --lmax=${lmax}
fi
# ------------------------------------------------------------------------------
echo Calculating spectra
outdir=${outdir}/output_lmax${lmax}_lsst-${lsst_data_tag}/
# run the code to calculate spectra
python ${repo_path}/code/get_spectra.py --outdir=${outdir}
