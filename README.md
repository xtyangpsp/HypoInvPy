# HypoInvPy
*A Python interface for HypoInverse earthquake relocation software

![plot1](/figs/hypoinvpylogo.png)

This package is not intended to show documentation of how to use `hypoinverse`. Instead, we focus on buiding an user-friendly Python interface to run `hypoinverse`. For detailed documentation of `hypoinverse`, please read the documentation for hypoinverse in `hyp1.40/doc/hyp1.40.pdf`. The documentation for this interface package is in `doc/HypoInvPy_manual.pdf`.

## Credit note
The codes and workflow in this package are modified and simplified from Hypo-Interface-Py (https://github.com/YijianZhou/Hypo-Interface-Py). The hypoinverse version is 1.4 and can be downloaded from USGS (https://www.usgs.gov/software/hypoinverse-earthquake-location). 

## Installation
1. Create and activate a Python virtual environment `hypoinv`. This is optional but recommended. This will help isolate the computational needs from other packages. This will also help avoid version incompatibility in case new packages for some codes are updated in the `base` conda environment. However, since `HypoInvPy` doesn't currently need complicated packages, running under the `base` environment may just work fine. 
```
$ conda create -n hypoinv -c conda-forge jupyter numpy scipy pandas python cartopy obspy mpi4py
```
Then, activate the environment with: `$ conda activate hypoinv`

2. Install the package with the parameters in `setup.py`.
```
$ pip install .
```

3. Create jupyter notebook kernel with the environment (after activating the environment).
```
$ pip install --user ipykernel
$ python -m ipykernel install --user --name=hypoinv
```

## Run the hypoinverse example:
1. Make sure hypo1.40 has been installed on your computer. Download the version from the USGS website (https://www.usgs.gov/software/hypoinverse-earthquake-location). For your convinience, a copy of the hyp1.40 codes is available under folder `hyp1.40`.
2. In terminal under the package directory, type: python run_hypo.py
3. The output directory includes the resultant earthquake catalogs.

The `run_hypo.py` script already includes steps to reformat the phase and station information. 

## Upcoming changes:
* Incorporate the functionality into a standalone python package.
* Change mk_pha.py and mk_sta.py into functions.
* Embed the template parameters into functions.
* Add utility functions to convert data formats.
* Add visualization functionalities, e.g., earthquake locations and depth distribution and temporal sequence, etc.




