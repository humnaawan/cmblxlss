## Goal: end-to-end simulation of LSST-LSSxSO-CMBLensing cross-correlations.

We use the Takahashi sims; see [here](http://cosmo.phys.hirosaki-u.ac.jp/takahasi/allsky_raytracing/nres12.html) for details for sims at Nside= 4096.

For each realization, we get three things from the sims linked above:
- CMB convergence map: `zs66` weak lensing map: a `.mag.dat` file that has 3 other maps (*strange to read, so needs a check*).
- Lensed CMB map: the `beta=0` link: a `.fits` file.
- The halo catalog: a 40G `ascii` file when unzipped. Impossible to read with `pandas` so use a python-analog of the IDL `readcol` routine; from [here](http://www.adamgginsburg.com/pyreadcol.htm).

#### LSST Galaxy Sample Cuts
These are the cuts applied to the halo catalog for now:
- parent_ID= -1 so that choose only the central haloes (i.e. not subhaloes).
- 5e12<Mvir<8e13

#### NERSC
The halo catalog file (for each realization) is big enough that needed to resort to NERSC to read it. Also, needed to compile `lenspix` (see below) which requires Fortran compiler; easy to handle on NERSC.
The files from the sims (for one realization only, rn) can be found on NERSC: `/global/cscratch1/sd/awan/takahashi_data/`

#### Lenspix
The [lenspix](https://github.com/cmbant/lenspix/tree/master) package contains the routine [`LensReconExample.f90`](https://github.com/cmbant/lenspix/blob/master/LensReconExample.f90) that does full sky lensing reconstruction when given a lensed CMB map. The code is compiled on NERSC using instructions written out [here](https://gist.github.com/humnaawan/60d742960060613b08e09a175073799f).

## Analysis Details
[Need to be written out]

## Repo structure
[Needs to be written out]
