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

def BasicPlot(param,W,F,E,l,f_fit,f_abs_ism,f_abs_bp,f_abs_X):

    c1  = param["fit"]["windows"]["window1"]["cut1"]
    c2  = param["fit"]["windows"]["window1"]["cut2"]

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
    plt.plot(l,f_abs_X,color="purple",lw=3)  
    plt.xlim(param["display"]["window1"]["x1"],param["display"]["window1"]["x2"])
    plt.ylim(param["display"]["window1"]["y1"],param["display"]["window1"]["y2"])
    plt.show()

def PrintParams(P):
    print "\nISM:"
    print "-"*50
    print "\t     b\t\t=\t",P[0],"km/s"

    print "\nCS:"
    print "-"*50
    print "\tlog(N/1cm^2)\t=\t",P[1]
    print "\t     b\t\t=\t",P[2],"km/s"

    print "\nExocomet:"
    print "-"*50
    print "\tlog(N/1cm^2)\t=\t",P[3],"km/s"
    print "\t    RV\t\t=\t",P[4],"km/s"
    print "\t     b\t\t=\t",P[5],"km/s"
    print "-"*50,"\n"

def main():
    # Read all parameters from params.json file.
    param           = Initialise()
    
    # Define the data directory
    dat_directory   = param["directories"]["workdir"]

    # Select the model type
    ModelType       = 1#param["fit"]["ModelType"]
    
    Nwindows   = param["fit"]["windows"]["number"]

    # Load the data file
    W, F, E         = np.genfromtxt(dat_directory+param["files"]["datafile"]
                      ,unpack=True)

    if Nwindows == 2:
        fit_start_w1   = param["fit"]["windows"]["window1"]["start"]
        fit_end_w1     = param["fit"]["windows"]["window1"]["stop"]
    
        s1_i = []
        s2_i = []
        for i in range(len(W)):
            if fit_start_w1 <= W[i] <= fit_end_w1:
                s1_i.append(i)
            if fit_start_w2 <= W[i] <= fit_end_w2:
                s2_i.append(i)
        
        W1    = W1[s1_i[0]:s1_i[-1]]
        F1    = F1[s1_i[0]:s1_i[-1]]
        E1    = E1[s1_i[0]:s1_i[-1]]

        W2    = W2[s2_i[0]:s2_i[-1]]
        F2    = F2[s2_i[0]:s2_i[-1]]
        E2    = E2[s2_i[0]:s2_i[-1]]
                      
        # Create an array of RV measurements with a resolution of 1 km/s
        v1    = np.arange(-len(W1)-250,len(W1)+250,1) # RV values
        v2    = np.arange(-len(W2)-250,len(W2)+250,1) # RV values
    
        # Calculate the corresponding wavelengths
        l1    = (W1[0]+W1[-1])/2.*(1.0 + v1/3e5)
        l2    = (W2[0]+W2[-1])/2.*(1.0 + v2/3e5)
    
    # Select the number of lines to model
    if param["lines"]["total"] == 5:
        L1  = param["lines"]["line"]["N1"]["Wavelength"]*(1.0 + v1/3e5)
        L2  = param["lines"]["line"]["N2"]["Wavelength"]*(1.0 + v1/3e5)
        L3  = param["lines"]["line"]["N3"]["Wavelength"]*(1.0 + v1/3e5)
        
        L4  = param["lines"]["line"]["Nw1"]["Wavelength"]*(1.0 + v2/3e5)
        L5  = param["lines"]["line"]["Nw2"]["Wavelength"]*(1.0 + v2/3e5)
    
    # Select the model type
    if ModelType == 1:
                
        # Fixed paramteres
        Const   =   [W1,W2,F1,F2,E1,E2,l1,l2,L1,L2,L3,L4,L5,param["BetaPictoris"]["RV"],
                    
                    # Fixed ISM parameters
                    param["fit"]["ISM"]["log(H)"],
                    param["fit"]["ISM"]["RV"],
                    param["fit"]["ISM"]["T"],
                    
                    # Fixed CS parameters
                    param["fit"]["disk"]["RV"],
                    param["fit"]["disk"]["T"],
                    
                    param["fit"]["exocomet"]["T"]]

                    # Free ISM parameters
        Par     =   [param["fit"]["ISM"]["b"],
                    
                    # Free CS parameters
                    param["fit"]["disk"]["log(H)"],
                    param["fit"]["disk"]["b"],
                    
                    # Free exocomet parameters
                    param["fit"]["exocomet"]["log(H)"],
                    param["fit"]["exocomet"]["RV"],
                    param["fit"]["exocomet"]["b"]]
        
        print "\nStarting paramters:"        
        PrintParams(Par)
        P =  FindBestParams(Par, F, E, Const, ModelType, param)
        print "Best fit paramters:"
        PrintParams(P)
        
        f_fit, f_abs_ism, f_abs_bp, f_abs_X   = m.LyModel(P,Const,ModelType,param)

        BasicPlot(param, W, F, E, l, f_fit, f_abs_ism, f_abs_bp, f_abs_X)


if __name__ == '__main__':
    main()
