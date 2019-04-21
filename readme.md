## Goal: end-to-end simulation of LSST-LSS x SO-CMB Lensing cross-correlations.

We calculate the cross-correlations between CMB lensing reconstructions and a (correlated) galaxy density field. The reconstructions account for (all) foregrounds and the LSST survey mask based on a depth-limited LSST footprint, while the correlated galaxy density field accounts for the LSST survey mask and modulation seen due to the systematics induced by the depth non-uniformities.

## Analysis Details
[Need to be written out]

## Repo structure
- The main code is in `code/analysis.py`, which uses the rest of the files on that folder.
- `scripts/run_analysis.sh` runs the analysis python script with specific inputs.
- The `legacy` folder contains code/docs from initial analyses.

## Dependencies
The code developed here needs [orphics](https://github.com/msyriac/orphics ), [pixell](https://github.com/simonsobs/pixell ), and [falafel](https://github.com/ACTCollaboration/falafel ).
