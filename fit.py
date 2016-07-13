#
#   Sky v0.1 Alpha
#   Please the required parameters in params.json.
#

import numpy as np
import matplotlib.pyplot as plt
import json, sys
from scipy.optimize import leastsq

# Load all the def functions used
from src.calculations import Calc
from src.model import Model
from src.statistics import Stats

c   = Calc()
m   = Model()
s   = Stats()

def Initialise():
    with open('params.json') as param_file:    
		param = json.load(param_file)
    return param

def FindBestParams(params,F,E,Const,ModelType, param):
    best_P, success = leastsq(s.chi2_lm, params, args=(F,E,Const,ModelType, param), maxfev=1000)
    return best_P

def BasicPlot(param,W,F,E):
    if param["display"]["bin"] > 1:
        Wb, Fb, Eb  = c.BinData(W,F,E,param["display"]["bin"])
        plt.step(Wb,Fb)
    else:
        plt.step(W,F)
    plt.xlim(param["display"]["window"]["x1"],param["display"]["window"]["x2"])
    plt.ylim(param["display"]["window"]["y1"],param["display"]["window"]["y2"])
    plt.show()

def main():
    # Read all parameters from params.json file.
    param           = Initialise()
    
    # Define the data directory
    dat_directory   = param["directories"]["workdir"]

    # Select the model type
    ModelType       = param["fit"]["ModelType"]

    # Load the data file
    W, F, E         = np.genfromtxt(dat_directory+param["files"]["datafile"]
                      ,unpack=True,skip_header=7400,skip_footer= 8500)
                      
    # Create an array of RV measurements with a resolution of 1 km/s
    v               = np.arange(-len(W)-300,len(W)+300,1) # RV values
    
    # Calculate the corresponding wavelengths
    l               = (W[0]+W[-1])/2.*(1.0 + v/3e5)
    #l               = param["lines"]["line"]["N3"]["Wavelength"]*(1.0 + v/3e5)
    
    # Select the number of lines to model
    if param["lines"]["total"] == 3:
        L1              = param["lines"]["line"]["N1"]["Wavelength"]*(1.0 + v/3e5)
        L2              = param["lines"]["line"]["N2"]["Wavelength"]*(1.0 + v/3e5)
        L3              = param["lines"]["line"]["N3"]["Wavelength"]*(1.0 + v/3e5)
    
    # Select the model type
    if ModelType == 1:
                
        # Fixed paramteres
        Const   =   [W,F,E,l,L1,L2,L3,param["BetaPictoris"]["RV"],
                    
                    # Fixed ISM parameters
                    param["fit"]["ISM"]["log(H)"],
                    param["fit"]["ISM"]["RV"],
                    param["fit"]["ISM"]["b"],
                    param["fit"]["ISM"]["T"],
                    
                    # Fixed disk parameters
                    param["fit"]["disk"]["b"],
                    param["fit"]["disk"]["T"]]

                    # Free disk parameters
        Par     =   [param["fit"]["disk"]["log(H)"],
                    param["fit"]["disk"]["RV"]]
        
        #print "Calculating the best parameters..."
        
        #f_before_fit, f_abs_ism, f_abs_bp   = m.LyModel(Par,Const,ModelType,
        #param["fit"]["continuum"]["start"],param["fit"]["continuum"]["stop"])
        #'''
        #X = F, E, m.LyModel(Par, Const, ModelType)[0]
        print Par
        print "\nBest fit paramters:"
        P =  FindBestParams(Par, F, E, Const, ModelType, param)
        
        print P
        #sys.exit()
        
        f_before_fit, f_abs_ism, f_abs_bp   = m.LyModel(P,Const,ModelType,param)
        
        plt.errorbar(W,np.ones(len(W))*3e-14,yerr=E)
        plt.step(W,F,color="black")
        plt.plot(W,f_before_fit,lw=3,color='#FF281C',label=r'Best fit')
        plt.plot(l,f_abs_ism,color="green",lw=3)
        plt.plot(l,f_abs_bp,color="blue",lw=3)
        plt.xlim(param["display"]["window"]["x1"],param["display"]["window"]["x2"])
        plt.ylim(param["display"]["window"]["y1"],param["display"]["window"]["y2"])
        plt.show()
        #'''
    #BasicPlot(param, W, F, E)



if __name__ == '__main__':
    main()

