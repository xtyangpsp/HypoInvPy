# HypoInvPy
*A Python interface for HypoInverse earthquake relocation software*

![plot1](/figs/hypoinvpylogo.png)

This package is not intended to show documentation of how to use `hypoinverse`. Instead, we focus on buiding an user-friendly Python interface to run `hypoinverse`. For detailed documentation of `hypoinverse`, please read the documentation for hypoinverse in `hyp1.40/doc/hyp1.40.pdf`. The documentation for this interface package is in `doc/HypoInvPy_manual.pdf`.

## Credit note
The codes and workflow in this package are modified and simplified from Hypo-Interface-Py (https://github.com/YijianZhou/Hypo-Interface-Py). The hypoinverse version is 1.4 and can be downloaded from USGS (https://www.usgs.gov/software/hypoinverse-earthquake-location). 

## Installation
### Install as a pip package
HypoInvPy is available on pypl as a standalone package. It could be installed as regular pip package: `pip install hypoinvpy`. The latest version is always available on GitHub.

### Install with local copy from GitHub
This installation method will get the latest version from GitHub.

1. Create and activate a Python virtual environment `hypoinv`. This is optional but recommended. This will help isolate the computational needs from other packages. This will also help avoid version incompatibility in case new packages for some codes are updated in the `base` conda environment. However, since `HypoInvPy` doesn't currently need complicated packages, running under the `base` environment may just work fine. 
```
$ conda create -n hypoinv -c conda-forge jupyter numpy scipy pandas python cartopy obspy mpi4py
```
Then, activate the environment with: `$ conda activate hypoinv`

2. Clone the github repository, in terminal under the desired directory:
```
$ git clone https://github.com/xtyangpsp/HypoInvPy.git
```

3. `cd` to the repository folder in terminal and install the package with the parameters in `setup.py`.
```
$ pip install .
```

4. Create jupyter notebook kernel with the environment (after activating the environment).
```
$ pip install --user ipykernel
$ python -m ipykernel install --user --name=hypoinv
```

## Run the hypoinverse example:
1. Make sure hypo1.40 has been installed on your computer. Download the version from the USGS website (https://www.usgs.gov/software/hypoinverse-earthquake-location). For your convinience, a copy of the hyp1.40 codes is available under folder `hyp1.40`.
2. In terminal under the `example` directory, run the jupyter notebook `HypoInvPy_run_example.ipynb`
3. The output directory includes the resultant earthquake catalogs.

## Upcoming changes:
* Add visualization functionalities, e.g., earthquake locations and depth distribution and temporal sequence, etc.




