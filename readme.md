## Goal: end-to-end simulation of LSST-LSS x SO-CMB Lensing cross-correlations.

We calculate the cross-correlations between CMB lensing reconstructions and a (correlated) galaxy density field. The lensing reconstructions currently account for (stacked) foregrounds and the LSST survey mask based on a depth-limited LSST footprint, while the correlated galaxy density field accounts for the LSST survey mask and modulation seen due to the systematics induced by the depth non-uniformities (with, without MW dust).

## Analysis Details
[Need to be written out, once finalized]

## Repo structure
- The main code is in two scripts: `code/save_needed_alms.py` calculates the various alms and saves them, while `code/get_spectra.py` reads in the saved alms to calculate and plot the spectra. These scripts which use the rest of the files in that folder.
- `scripts/run_analysis.sh` runs the python scripts with specific inputs.
- The `legacy` folder contains code/docs from initial analyses/exploration.

## Dependencies
The code developed here needs [orphics](https://github.com/msyriac/orphics ), [pixell](https://github.com/simonsobs/pixell ), and [falafel](https://github.com/ACTCollaboration/falafel ).
