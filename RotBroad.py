from __future__ import print_function, division
import numpy as np
import matplotlib.pylab as plt
from PyAstronomy import funcFit as fuf
from PyAstronomy import pyasl
import scipy.integrate as sci
import json

def Initialise():
    with open('params.json') as param_file:    
		param = json.load(param_file)
    return param

W_H, F_H  = np.genfromtxt('/home/paw/science/betapic/data/stellar/BetaPic_flux_sum_2003-2015.txt',unpack=True,skip_header=1)


Wx  = np.genfromtxt('/home/paw/science/betapic/data/stellar/UVBLUE_wavelengths.dat',unpack=True)
Fx, Fy  = np.genfromtxt('/home/paw/science/betapic/data/stellar/t08000g45p00k2.flx',unpack=True,skip_header=3)
Wx = Wx[12200:14000]   # NI wavelengths
Fx = Fx[12200:14000]

#Wx = Wx[60900:63500]   # CaII wavelengths
#Fx = Fx[60900:63500]

#plt.plot(W_H+1.16, F_H)
#plt.plot(Wx, Fx)
#'''
Wx = np.linspace(Wx[0],Wx[-1],len(Wx))
# Obtain the broadened spectrum using
# vsini = 80 km/s and no limb-darkening
rflux = pyasl.rotBroad(Wx, Fx, 0.0, 130)

# Obtain the broadened spectrum using
# vsini = 80 km/s and strong limb-darkening
lflux = pyasl.rotBroad(Wx, Fx, 0.9, 132)

# Check that the area of the line did not change
# in response to the broadening
print("Initial EW [A]: ", 4. - sci.trapz(Fx, Wx))
print("After broadening without LD: ", 4. - sci.trapz(rflux, Wx))
print("After broadening with LD: ", 4. - sci.trapz(lflux, Wx))

# Read all parameters from params.json file.
param           = Initialise()

# Define the data directory
dat_directory   = param["directories"]["workdir"]
np.savetxt(dat_directory+"N_RotB132L09.dat",np.column_stack((Wx, lflux)))

# Plot the results
plt.title("Rotational broadening")
plt.xlabel("Wavelength [A]")
plt.ylabel("Normalized flux")
#plt.plot(Wx, Fx, 'b-')
plt.plot(Wx, rflux*0.75, 'r-')
plt.plot(Wx, lflux*0.75, 'g-')
#plt.xlim(3900,4200)
#'''
plt.show()
