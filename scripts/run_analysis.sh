#!/bin/bash
################################################################################
# source the stack
source /global/common/software/lsst/cori-haswell-gcc/stack/setup_current_sims.sh
#--------------------------------------------------------------------------------------
# set up inputs
repo_path=/global/homes/a/awan/SO/CMBL-LSS-CrossCorr/
nside=1024
lmax=3000

# cmb related inputs
cmb_tag=sim0
conv_path=/global/cscratch1/sd/awan/soxlsst/takahashi_data/allskymap_nres12r000.zs66.mag.dat
halocat_path=/global/homes/a/awan/SO/CMBL-LSS-CrossCorr/data/haloCat_readcol_39773347gals.csv
#halocat_path=/global/cscratch1/sd/awan/soxlsst/takahashi_data//interm_debug//finalOutput_505gals.csv

# lsst related inputs
lsst_mask_path=${repo_path}/data/deltaNByNData_masked_i_RandomDitherFieldPerVisit.npz
lsst_tag=minion1016_i-band_RandomFPV_deltaNByN

# output directory
outdir=${repo_path}/outputs

# run the code
python ${repo_path}/code/analysis.py \
                    --outdir=${outdir} \
                    --nside=${nside} --lmax=${lmax} --halocat-path=${halocat_path} \
                    --convergence-path=${conv_path} --cmb-sim-tag=${cmb_tag} \
                    --lsst-mask-path=${lsst_mask_path} \
                    --lsstdata-tag=${lsst_tag}
