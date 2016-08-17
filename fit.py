#
#   Sky v0.1 Alpha
#   Please change the required parameters in params.json.
#

import numpy as np
import matplotlib.pyplot as plt
import json, sys
from scipy.optimize import leastsq

from colorama import init
init()
from colorama import Fore, Back, Style

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
    best_P, success = leastsq(s.chi2_lm, params, args=(F,E,Const,ModelType, param), maxfev=10000)
    return best_P

def BasicPlot(param,window,W,F,E,l,f_fit,f_abs_ism,f_abs_bp,f_abs_X,unconvolved):

    x1  = param["display"][window]["x1"]
    x2  = param["display"][window]["x2"]
    y1  = param["display"][window]["y1"]
    y2  = param["display"][window]["y2"]
    
    c1  = param["fit"]["windows"][window]["cut1"]
    c2  = param["fit"]["windows"][window]["cut2"]

    fig = plt.figure(figsize=(12,6))

    if param["display"]["bin"] > 1:
        bin_size = param["display"]["bin"]
        Wb, Fb, Eb  = c.BinData(W,F,E,bin_size)
        plt.errorbar(Wb,np.ones(len(Wb))*2e-14,yerr=Eb)
        plt.step(Wb,Fb)
        plt.step(Wb[:c1/bin_size],Fb[:c1/bin_size],color="red")
        plt.step(Wb[-c2/bin_size:],Fb[-c2/bin_size:],color="red")
    else:
        line1 = param["lines"]["line"]["Nw1"]["Wavelength"]
        line2 = param["lines"]["line"]["Nw2"]["Wavelength"]
        plt.plot([line1,line1],[0.2e-14,0.3e-14],color='black')
        plt.plot([line2,line2],[0.2e-14,0.3e-14],color='black')
        plt.errorbar(W,np.ones(len(W))*2e-14,yerr=E)
        plt.step(W,F)
        plt.step(W[:c1],F[:c1],color="red")
        plt.step(W[-c2:],F[-c2:],color="red")
    
    plt.plot(l,unconvolved,color="cyan")
    plt.plot(l,f_abs_ism,color="green",lw=3)
    plt.plot(l,f_abs_bp,color="blue",lw=3)   
    plt.plot(l,f_abs_X,color="purple",lw=3)
    plt.plot(W,f_fit,lw=3,color='#FF281C',label=r'Best fit')
    plt.xlim(x1,x2)
    plt.ylim(y1,y2)
    plt.show()

def PrintParams(P, ConstB):
    print "\nISM:"
    print "-"*50,Fore.GREEN
    print "Free parameters:\n"
    print "\t     b\t\t=\t",P[0],"km/s"
    print(Fore.RED)
    print "\nConstant parameters:\n"
    print "\tlog(N/1cm^2)\t=\t",ConstB[1],"km/s"
    print "\t    RV\t\t=\t",ConstB[2],"km/s"
    print "\t     T\t\t=\t",ConstB[3],"K",Style.RESET_ALL
    print "-"*50,"\n"

    print "\nCS:"
    print "-"*50,Fore.GREEN
    print "Free parameters:\n"
    print "\tlog(N/1cm^2)\t=\t",P[1]
    print "\t     b\t\t=\t",P[2],"km/s"
    print(Fore.RED)
    print "\nConstant parameters:\n"
    print "\t    RV\t\t=\t",ConstB[4],"km/s"
    print "\t     T\t\t=\t",ConstB[5],"K",Style.RESET_ALL
    print "-"*50,"\n"

    print "\nExocomet:"
    print "-"*50,Fore.GREEN
    print "Free parameters:\n"
    print "\tlog(N/1cm^2)\t=\t",P[3],"km/s"
    print "\t    RV\t\t=\t",P[4],"km/s"
    print "\t     b\t\t=\t",P[5],"km/s"
    print(Fore.RED)
    print "\nConstant parameters:\n"
    print "\t     T\t\t=\t",ConstB[6],"K",Style.RESET_ALL
    print "-"*50,"\n\n\n"

def Window(param,W,F,E,WindowName):
    fit_start   = param["fit"]["windows"][WindowName]["start"]
    fit_end     = param["fit"]["windows"][WindowName]["stop"]

    s_i = []
    for i in range(len(W)):
        if fit_start <= W[i] <= fit_end:
            s_i.append(i)
    
    W    = W[s_i[0]:s_i[-1]]
    F    = F[s_i[0]:s_i[-1]]
    E    = E[s_i[0]:s_i[-1]]

    # Create an array of RV measurements with a resolution of 1 km/s
    v    = np.arange(-len(W)-500,len(W)+500,0.1) # RV values

    # Calculate the corresponding wavelengths
    l    = (W[0]+W[-1])/2.*(1.0 + v/3e5)

    return W, F, E, v, l

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

    if Nwindows == 1:
        W1, F1, E1, v1, l1  = Window(param,W,F,E,"window1")
        ConstA              = [W1,F1,E1,l1]
    
    if Nwindows == 2:
        W1, F1, E1, v1, l1  = Window(param,W,F,E,"window1")
        W2, F2, E2, v2, l2  = Window(param,W,F,E,"window2")
        ConstA              = [W1,W2,F1,F2,E1,E2,l1,l2]

    # The parameters listed in ConstB are not dependant on the number of windows used.
    # ConstA is.       
    ConstB   =  [param["BetaPictoris"]["RV"],
                
                # Fixed ISM parameters
                param["fit"]["ISM"]["log(H)"],
                param["fit"]["ISM"]["RV"],
                param["fit"]["ISM"]["T"],
                
                # Fixed CS parameters
                param["fit"]["disk"]["RV"],
                param["fit"]["disk"]["T"],
                
                # Fixed exocomet parameters
                param["fit"]["exocomet"]["T"]]

    Const   =   np.concatenate((ConstA,ConstB))
    
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
   
    PrintParams(Par, ConstB)
    P =  FindBestParams(Par, F1, E1, Const, ModelType, param)
    print "Best fit paramters:"
    PrintParams(P, ConstB)
    
    if Nwindows == 1:
        f_fit1, f_abs_ism1, f_abs_bp1, f_abs_X1, unconvolved1 = m.Model(P,Const,ModelType,param)
        BasicPlot(param, param["display"]["window1"]["name"], W1, F1, E1, l1, f_fit1, f_abs_ism1, f_abs_bp1, f_abs_X1, unconvolved1)    
    if Nwindows == 2:
        f_fit1, f_abs_ism1, f_abs_bp1, f_abs_X1, unconvolved1, f_fit2, f_abs_ism2, f_abs_bp2, f_abs_X2, unconvolved2 = m.Model(P,Const,ModelType,param)
        BasicPlot(param, param["display"]["window1"]["name"], W1, F1, E1, l1, f_fit1, f_abs_ism1, f_abs_bp1, f_abs_X1, unconvolved1) 
        BasicPlot(param, param["display"]["window2"]["name"], W2, F2, E2, l2, f_fit2, f_abs_ism2, f_abs_bp2, f_abs_X2, unconvolved2)

if __name__ == '__main__':
    main()
