# N7

[![GPL](https://img.shields.io/badge/license-GNU%20GPLv3-brightgreen.svg)](http://choosealicense.com/licenses/gpl-3.0/)

N7 v0.1 is released!

N7 is a custom built code for modelling the NI lines in HST/COS data. The code models the NI lines found around 1160 and 1200 Angstrom. All the parameters used in the fitting and plotting are available in params.json. The code is under active development and will be featured in upcoming work by Wilson et al. (2017).

**fit.py** - Fits the nitrogen tripplet using voigt profiles and estimates column densities.

**MCMC_XA.py** - Code specifically designed to run 24 MCMC chains simultaneously on the exoatmos super computer.

**params.json** - Parameter file which contains all parameters used in the modelling of the NI lines.