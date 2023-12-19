# HypoInvPy
Python interface for hypoinverse for relocating earthquakes. The codes and workflow are modified from Hypo-Interface-Py (https://github.com/YijianZhou/Hypo-Interface-Py). The hypoinverse version is 1.4 (https://www.usgs.gov/software/hypoinverse-earthquake-location). 

## How to run the hypoinverse example:
1. Make sure hypo1.40 has been installed on your computer. Download the version from the USGS website (https://www.usgs.gov/software/hypoinverse-earthquake-location).
2. In terminal under the package directory, type: python run_hypo.py
3. The output directory includes the resultant earthquake catalogs.

The `run_hypo.py` script already includes steps to reformat the phase and station information. 

## Upcoming changes:
* Incorporate the functionality into a standalone python package.
* Change mk_pha.py and mk_sta.py into functions.
* Embed the template parameters into functions.
* Add utility functions to convert data formats.
* Add visualization functionalities, e.g., earthquake locations and depth distribution and temporal sequence, etc.




