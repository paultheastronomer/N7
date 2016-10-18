'''
    fit.py calculates an optimal fit to the NI lines at
    ~1160 Ang. and ~1200 Ang. using a least-squares minimisation
    routine

    To change the parameters of the fit please make changes to
    params.json (the parameter file) only.
'''

import numpy as np
import json, sys


# Load the package necessary for priting in colour in the terminal.
from colorama import Fore, Back, Style

# Load all the custom def functions used.
# These functions can be found in the src directory.
from src.calculations import Calc
from src.model import Model
from src.statistics import Stats
from src.plotting import Plotting

c   = Calc()
m   = Model()
s   = Stats()
p   = Plotting()

def Initialise():
    with open('params.json') as param_file:    
		param = json.load(param_file)
    return param


def main():
    # Read all parameters from params.json file.
    param           = Initialise()
    
    # Define the data directory.
    dat_directory   = param["directories"]["workdir"]
    
    # Select the model type
    # This option is yet to be implemented.
    ModelType       = 1#param["fit"]["ModelType"]
    
    # The number of spectral windows used for fitting.
    Nwindows        = param["fit"]["windows"]["number"]

    ''' Load the data file
    W = Wavelength
    F = Flux
    E = Uncertainty on F '''
    W, F, E         = np.genfromtxt(dat_directory+param["files"]["datafile"]
                      ,unpack=True)

    # Need to change later
    F = F*8.e-8

    if Nwindows == 1:
        W1, F1, E1, v1, l1  = c.Window(param,W,F,E,"window1")
        ConstA              = [W1,F1,E1,l1]
    
    if Nwindows == 2:
        W1, F1, E1, v1, l1  = c.Window(param,W,F,E,"window1")
        W2, F2, E2, v2, l2  = c.Window(param,W,F,E,"window2")
        ConstA              = [W1,W2,F1,F2,E1,E2,l1,l2]

    if Nwindows == 3:
        W1, F1, E1, v1, l1  = c.Window(param,W,F,E,"window1")
        W2, F2, E2, v2, l2  = c.Window(param,W,F,E,"window2")
        W3, F3, E3, v3, l3  = c.Window(param,W,F,E,"window3")
        ConstA              = [W1,W2,W3,F1,F2,F3,E1,E2,E3,l1,l2,l3]

    # The parameters listed in ConstB are not dependant on the number of windows used.
    # ConstA is.

    #           Radial velocity of the beta Pictoris system. Default 20.5 km/s.       
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
 
    if param["fit"]["dofit"] == "yes":
        P =  c.FindBestParams(Par, F1, E1, Const, ModelType, param)
    else:
        P = Par
        #'''
        #c.PrintParams(Par, ConstB)
        #P =  c.FindBestParams(Par, F1, E1, Const, ModelType, param)
        #print "Best fit paramters:"
        #c.PrintParams(P, ConstB)
        #'''

    # Using the parameters above create a model then plot the result.
    if Nwindows == 1:
        f_fit1, f_abs_con1, f_abs_ism1, f_abs_bp1, f_abs_X1, unconvolved1 = m.Model(P,Const,ModelType,param)
        p.BasicPlot(param, param["display"]["window1"]["name"], W, F, E, l1, f_abs_con1, f_abs_ism1, f_abs_bp1, f_abs_X1, unconvolved1)    
    if Nwindows == 2:
        f_fit1, f_abs_con1, f_abs_ism1, f_abs_bp1, f_abs_X1, unconvolved1, f_fit2, f_abs_con2, f_abs_ism2, f_abs_bp2, f_abs_X2, unconvolved2 = m.Model(P,Const,ModelType,param)
        p.BasicPlot(param, param["display"]["window1"]["name"], W, F, E, l1, f_abs_con1, f_abs_ism1, f_abs_bp1, f_abs_X1, unconvolved1) 
        p.BasicPlot(param, param["display"]["window2"]["name"], W, F, E, l2, f_abs_con2, f_abs_ism2, f_abs_bp2, f_abs_X2, unconvolved2)
    if Nwindows == 3:
        f_fit1, f_abs_con1, f_abs_ism1, f_abs_bp1, f_abs_X1, unconvolved1, f_fit2, f_abs_con2, f_abs_ism2, f_abs_bp2, f_abs_X2, unconvolved2, f_fit3, f_abs_con3, f_abs_ism3, f_abs_bp3, f_abs_X3, unconvolved3 = m.Model(P,Const,ModelType,param)
        p.BasicPlot(param, param["display"]["window1"]["name"], W, F, E, l1, f_abs_con1, f_abs_ism1, f_abs_bp1, f_abs_X1, unconvolved1) 
        p.BasicPlot(param, param["display"]["window2"]["name"], W, F, E, l2, f_abs_con2, f_abs_ism2, f_abs_bp2, f_abs_X2, unconvolved2)
        p.BasicPlot(param, param["display"]["window3"]["name"], W, F, E, l3, f_abs_con3, f_abs_ism3, f_abs_bp3, f_abs_X3, unconvolved3)

    # Compute the goodness of fit.
    #X = [F1, E1, f_fit1]
    #print "DOF:\t\t",len(W1)-len(P)
    #print "Chi2 reduced:\t",s.chi2(X)/(len(F1)-len(P)),"\n"

if __name__ == '__main__':
    main()