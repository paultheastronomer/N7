#
#   Sky v0.1 Alpha
#   Please the required parameters in params.json.
#

import numpy as np
import matplotlib.pyplot as plt
import json, sys
from scipy.optimize import leastsq

# Load all the custom def functions used
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

def BasicPlot(param,W,F,E,l,f_fit,f_abs_ism,f_abs_bp):

    c1  = param["fit"]["continuum"]["cut1"]
    c2  = param["fit"]["continuum"]["cut2"]

    if param["display"]["bin"] > 1:
        bin_size = param["display"]["bin"]
        Wb, Fb, Eb  = c.BinData(W,F,E,bin_size)
        plt.errorbar(Wb,np.ones(len(Wb))*2e-14,yerr=Eb)
        plt.step(Wb,Fb)
        plt.step(Wb[:c1/bin_size],Fb[:c1/bin_size],color="red")
        plt.step(Wb[-c2/bin_size:],Fb[-c2/bin_size:],color="red")
    else:
        plt.errorbar(W,np.ones(len(W))*2e-14,yerr=E)
        plt.step(W,F)
        plt.step(Wb[:c1],F[:c1],color="red")
        plt.step(W[-c2:],F[-c2:],color="red")
    
    plt.plot(W,f_fit,lw=3,color='#FF281C',label=r'Best fit')
    plt.plot(l,f_abs_ism,color="green",lw=3)
    plt.plot(l,f_abs_bp,color="blue",lw=3)    
    plt.xlim(param["display"]["window"]["x1"],param["display"]["window"]["x2"])
    plt.ylim(param["display"]["window"]["y1"],param["display"]["window"]["y2"])
    plt.show()

def PrintParams(P):
    print "\tlog(N/1cm^2)\t=\t",P[0]
    #print "\t   v\t\t=\t",P[1],"km/s"
    #print "\t   b\t\t=\t",P[2],"km/s"
    print "\n"

def main():
    # Read all parameters from params.json file.
    param           = Initialise()
    
    # Define the data directory
    dat_directory   = param["directories"]["workdir"]

    # Select the model type
    ModelType       = param["fit"]["ModelType"]

    # Load the data file
    #W, F, E         = np.genfromtxt(dat_directory+param["files"]["datafile"]
    #                  ,unpack=True,skip_header=7400,skip_footer= 8500)
    W, F, E         = np.genfromtxt(dat_directory+param["files"]["datafile"]
                      ,unpack=True)

    fit_start   = param["fit"]["continuum"]["start"]
    fit_end     = param["fit"]["continuum"]["stop"]
    
    s_i = []
    for i in range(len(W)):
        if fit_start <= W[i] <= fit_end:
            s_i.append(i)
    
    W    = W[s_i[0]:s_i[-1]]
    F    = F[s_i[0]:s_i[-1]]
    E    = E[s_i[0]:s_i[-1]]
                      
    # Create an array of RV measurements with a resolution of 1 km/s
    v               = np.arange(-len(W)-250,len(W)+250,1) # RV values
    
    # Calculate the corresponding wavelengths
    l               = (W[0]+W[-1])/2.*(1.0 + v/3e5)
    
    # Select the number of lines to model
    if param["lines"]["total"] == 3:
        L1  = param["lines"]["line"]["N1"]["Wavelength"]*(1.0 + v/3e5)
        L2  = param["lines"]["line"]["N2"]["Wavelength"]*(1.0 + v/3e5)
        L3  = param["lines"]["line"]["N3"]["Wavelength"]*(1.0 + v/3e5)
    
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
                    param["fit"]["disk"]["RV"],
                    param["fit"]["disk"]["b"],
                    param["fit"]["disk"]["T"]]

                    # Free disk parameters
        Par     =   [param["fit"]["disk"]["log(H)"]]
        
        print "\nStarting paramters:"        
        PrintParams(Par)
        P =  FindBestParams(Par, F, E, Const, ModelType, param)
        print "Best fit paramters:"
        PrintParams(P)
        
        f_fit, f_abs_ism, f_abs_bp   = m.LyModel(P,Const,ModelType,param)

        BasicPlot(param, W, F, E, l, f_fit, f_abs_ism, f_abs_bp)


if __name__ == '__main__':
    main()

