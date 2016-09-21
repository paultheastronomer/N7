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

    fig = plt.figure(figsize=(8,5))

    fontlabel_size  = 18
    tick_size       = 18
    params = {'backend': 'wxAgg', 'lines.markersize' : 2, 'axes.labelsize': fontlabel_size, 'font.size': fontlabel_size, 'legend.fontsize': 15, 'xtick.labelsize': tick_size, 'ytick.labelsize': tick_size, 'text.usetex': True}
    plt.rcParams.update(params)
    plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['text.latex.unicode'] = True 

    plt.plot(l,unconvolved,color="cyan")
    plt.plot(l,f_abs_ism,color="#0386FF",lw=2)
    plt.plot(l,f_abs_bp,color="#00B233",lw=2)   
    plt.plot(l,f_abs_X,color="#FF9303",lw=2)
    plt.plot(W,f_fit,lw=2,color='#FF281C',label=r'Best fit')

    if param["display"]["bin"] > 1:
        bin_size = param["display"]["bin"]
        Wb, Fb, Eb  = c.BinData(W,F,E,bin_size)
        plt.errorbar(Wb,np.ones(len(Wb))*2e-14,yerr=Eb)
        plt.step(Wb,Fb,color="#333333")
        plt.step(Wb[:c1/bin_size],Fb[:c1/bin_size],color="black",lw=2)
        plt.step(Wb[-c2/bin_size:],Fb[-c2/bin_size:],color="black",lw=2)
    else:
        line1 = param["lines"]["line"]["Nw1"]["Wavelength"]
        line2 = param["lines"]["line"]["Nw2"]["Wavelength"]
        plt.plot([line1,line1],[0.2e-14,0.3e-14],color='black')
        plt.plot([line2,line2],[0.2e-14,0.3e-14],color='black')
        plt.errorbar(W,np.ones(len(W))*2e-14,yerr=E)
        plt.step(W,F,color="#333333")
        plt.step(W[:c1],F[:c1],color="black",lw=2)
        plt.step(W[-c2:],F[-c2:],color="black",lw=2)
    

    plt.xlim(x1,x2)
    plt.ylim(y1,y2)
    '''
    if window == 'window1':
    	x = [1199,1200,1201]
    	labels = ['1199','1200','1201']
    	plt.xticks(x, labels)
    if window == 'window2':
    	x = [1160,1161]
    	labels = ['1160','1161']
    	plt.xticks(x, labels)
    '''
    plt.xlabel(r'Wavelength (\AA)')
    plt.ylabel(r'Flux (erg/s/cm$^2$/\AA)')

    plt.minorticks_on()
    fig.tight_layout()
    fig.savefig("plots/"+window+".pdf")
    plt.show()

def PrintParams(P, ConstB):
    print "\n",Fore.GREEN,"Free parameters",Style.RESET_ALL
    print Fore.RED,"Constant paramters",Style.RESET_ALL,"\n"
    
    print "ISM parameters:"
    print "-"*50,Fore.GREEN
    print "\tlog(N/1cm^2)\t=\t",P[0],"km/s"
    print "\t     b\t\t=\t",    P[1],Fore.RED
    print "\t     RV\t\t=\t",   ConstB[1],"km/s"
    print "\t     T\t\t=\t",    ConstB[2],"K",Style.RESET_ALL
    print "-"*50
    
    print "\nCS:"
    print "-"*50,Fore.GREEN
    print "\tlog(N/1cm^2)\t=\t",P[1]
    print "\t     b\t\t=\t",    P[2],"km/s",Fore.RED
    print "\t    RV\t\t=\t",    ConstB[3],"km/s"
    print "\t     T\t\t=\t",    ConstB[4],"K",Style.RESET_ALL
    print "-"*50,"\n"

    print "\nExocomet:"
    print "-"*50,Fore.GREEN
    print "\tlog(N/1cm^2)\t=\t",P[3],"km/s"
    print "\t    RV\t\t=\t",    P[4],"km/s"
    print "\t     b\t\t=\t",    P[5],"km/s",Fore.RED
    print "\t     T\t\t=\t",    ConstB[5],"K",Style.RESET_ALL
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
    v    = np.arange(-len(W)-100,len(W)+100,1) # RV values

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
    
    Nwindows        = param["fit"]["windows"]["number"]

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
                param["fit"]["ISM"]["RV"],
                param["fit"]["ISM"]["T"],
                
                # Fixed CS parameters
                param["fit"]["disk"]["RV"],
                param["fit"]["disk"]["T"],
                
                # Fixed exocomet parameters
                param["fit"]["exocomet"]["T"]]

    Const   =   np.concatenate((ConstA,ConstB))
    
                # Free ISM parameters
    Par     =   [param["fit"]["ISM"]["log(H)"],
                param["fit"]["ISM"]["b"],                
                
                # Free CS parameters
                param["fit"]["disk"]["log(H)"],
                param["fit"]["disk"]["b"],
                
                # Free exocomet parameters
                param["fit"]["exocomet"]["log(H)"],
                param["fit"]["exocomet"]["RV"],
                param["fit"]["exocomet"]["b"]]
 
    '''
    PrintParams(Par, ConstB)
    P =  FindBestParams(Par, F1, E1, Const, ModelType, param)
    print "Best fit paramters:"
    PrintParams(P, ConstB)
    '''
    #X = [F1, E1, f_fit1]
    #print "DOF:\t\t",len(RV)-len(Par)
    #print "Chi2 reduced:\t",s.chi2(X)/(len(F1)-len(Par)),"\n"

  
    if Nwindows == 1:
        f_fit1, f_abs_ism1, f_abs_bp1, f_abs_X1, unconvolved1 = m.Model(Par,Const,ModelType,param)
        BasicPlot(param, param["display"]["window1"]["name"], W1, F1, E1, l1, f_fit1, f_abs_ism1, f_abs_bp1, f_abs_X1, unconvolved1)    
    if Nwindows == 2:
        f_fit1, f_abs_ism1, f_abs_bp1, f_abs_X1, unconvolved1, f_fit2, f_abs_ism2, f_abs_bp2, f_abs_X2, unconvolved2 = m.Model(Par,Const,ModelType,param)
        BasicPlot(param, param["display"]["window1"]["name"], W1, F1, E1, l1, f_fit1, f_abs_ism1, f_abs_bp1, f_abs_X1, unconvolved1) 
        BasicPlot(param, param["display"]["window2"]["name"], W2, F2, E2, l2, f_fit2, f_abs_ism2, f_abs_bp2, f_abs_X2, unconvolved2)


if __name__ == '__main__':
    main()
