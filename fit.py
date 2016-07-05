import numpy as np
import matplotlib.pyplot as plt
import json
from scipy.optimize import leastsq
import sys

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
    param           = Initialise()
    dat_directory   = param["directories"]["workdir"]

    ModelType = 1

    W, F, E         = np.genfromtxt(dat_directory+param["files"]["datafile"]
                      ,unpack=True,skip_header=7300,skip_footer= 8500)
                      
    #print W[0],W[-1]
    #sys.exit()

    # RV Values
    v               = np.arange(-len(W),len(W),1) # RV values
    
    # Corresponding wavengths
    l               = param["lines"]["line"]["N1"]["Wavelength"]*(1.0 + v/3e5)
    
    if param["lines"]["total"] == 3:
        L1              = param["lines"]["line"]["N1"]["Wavelength"]*(1.0 + v/3e5)
        L2              = param["lines"]["line"]["N2"]["Wavelength"]*(1.0 + v/3e5)
        L3              = param["lines"]["line"]["N3"]["Wavelength"]*(1.0 + v/3e5)
    
    
    if ModelType == 1:
        
        # Free parameters
        Par     =   [param["fit"]["disk"]["log(H)"],
                    param["fit"]["disk"]["RV"]]
        
        # Fixed paramteres
        Const   =   [W,l,L1,L2,L3,param["BetaPictoris"]["RV"],
                    
                    # Fixed ISM parameters
                    param["fit"]["ISM"]["log(H)"],
                    param["fit"]["ISM"]["RV"],
                    param["fit"]["ISM"]["b"],
                    param["fit"]["ISM"]["T"],
                    
                    # Fixed disk parameters
                    param["fit"]["disk"]["b"],
                    param["fit"]["disk"]["T"]]
        
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

