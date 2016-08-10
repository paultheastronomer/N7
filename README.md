# sky

[![GPL](https://img.shields.io/badge/license-GNU%20GPLv3-brightgreen.svg)](http://choosealicense.com/licenses/gpl-3.0/)


Sky is the Norwegian word for cloud.

This code is currently being built and is not fit for use. The work using this code will be features in Wilson et al. (2016).

**fit.py** - Fits the nitrogen tripplet using voigt profiles and estimates column densities.

**plots/exocomets.py** - Shows a simple plot of spectra centered around the region of the NI line.

**convert2owens.py** - Converts a general text file into a **owens.f** friendly format.

**RotBroad.py** - Applies rotational broadening to a given spectrum using the formulae given in Gray’s “The Observation and Analysis of Stellar Photospheres”. It allows for limb darkening parameterized by the linear limb-darkening law. The code has been adapted from [pyatronomy](http://pyastronomy.readthedocs.io/en/latest/pyaslDoc/aslDoc/rotBroad.html).
