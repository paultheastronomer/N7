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
                
                # Fixed CS parameters
                param["fit"]["CS"]["RV"]]

    Const   =   np.concatenate((ConstA,ConstB))
    
                # ISM parameters which can be set free
    Par     =   [param["fit"]["ISM"]["log(N)"],
                param["fit"]["ISM"]["log(S)"],
                param["fit"]["ISM"]["b"],
                param["fit"]["ISM"]["T"],
                param["fit"]["ISM"]["xi"],                
                
                # CS parameters which can be set free
                param["fit"]["CS"]["log(N)"],
                param["fit"]["CS"]["log(S)"],
                param["fit"]["CS"]["b"],
                param["fit"]["CS"]["T"],
                param["fit"]["CS"]["xi"],
                
                # Exocomet parameters which can be set free
                param["fit"]["exocomet1"]["log(N)"],
                param["fit"]["exocomet1"]["log(S)"],
                param["fit"]["exocomet1"]["b"],
                param["fit"]["exocomet1"]["T"],
                param["fit"]["exocomet1"]["xi"],
                param["fit"]["exocomet1"]["RV"],

                param["fit"]["exocomet2"]["log(N)"],
                param["fit"]["exocomet2"]["log(S)"],
                param["fit"]["exocomet2"]["b"],
                param["fit"]["exocomet2"]["T"],
                param["fit"]["exocomet2"]["xi"],
                param["fit"]["exocomet2"]["RV"]]

    if param["fit"]["dofit"] == "yes":
        P =  c.FindBestParams(Par, F1, E1, Const, ModelType, param)
        print P
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
        f_abs_int_w1, f_abs_con_w1, f_abs_ism_w1, f_abs_bp_w1, f_abs_X1_w1, f_abs_X2_w1, unconvolved_w1 = m.Model(P,Const,ModelType,param)
        p.BasicPlot(param, param["display"]["window1"]["name"], W, F, E, l1, f_abs_con_w1, f_abs_ism_w1, f_abs_bp_w1, f_abs_X1_w1, f_abs_X2_w1, unconvolved_w1)    
    if Nwindows == 2:
        f_abs_int_w1, f_abs_con_w1, f_abs_ism_w1, f_abs_bp_w1, f_abs_X1_w1, f_abs_X2_w1, unconvolved_w1, f_abs_int_w2, f_abs_con_w2, f_abs_ism_w2, f_abs_bp_w2, f_abs_X1_w2, f_abs_X2_w2, unconvolved_w2 = m.Model(P,Const,ModelType,param)
        p.BasicPlot(param, param["display"]["window1"]["name"], W, F, E, l1, f_abs_con_w1, f_abs_ism_w1, f_abs_bp_w1, f_abs_X1_w1, f_abs_X2_w1, unconvolved_w1) 
        p.BasicPlot(param, param["display"]["window2"]["name"], W, F, E, l2, f_abs_con_w2, f_abs_ism_w2, f_abs_bp_w2, f_abs_X1_w2, f_abs_X2_w2, unconvolved_w2)
        X = [F2, E2, f_abs_int_w2]
        print "Chi2:\t",s.chi2(X)
        print "DOF:\t\t",len(W1)-len(P)
        print "Chi2 reduced:\t",s.chi2(X)/(len(F2)-len(P)),"\n"

    if Nwindows == 3:
        f_abs_int_w1, f_abs_con_w1, f_abs_ism_w1, f_abs_bp_w1, f_abs_X1_w1, f_abs_X2_w1, f_abs_X3_w1, unconvolved_w1, f_abs_int_w2, f_abs_con_w2, f_abs_ism_w2, f_abs_bp_w2, f_abs_X1_w2, f_abs_X2_w2, f_abs_X3_w2, unconvolved_w2, f_abs_int_w3, f_abs_con_w3, f_abs_ism_w3, f_abs_bp_w3,  f_abs_X1_w3, f_abs_X2_w3, f_abs_X3_w3, unconvolved_w3 = m.Model(P,Const,ModelType,param)
        p.BasicPlot(param, param["display"]["window1"]["name"], W, F, E, l1, f_abs_con_w1, f_abs_ism_w1, f_abs_bp_w1, f_abs_X1_w1, f_abs_X2_w1, f_abs_X3_w1, unconvolved_w1) 
        p.BasicPlot(param, param["display"]["window2"]["name"], W, F, E, l2, f_abs_con_w2, f_abs_ism_w2, f_abs_bp_w2, f_abs_X1_w2, f_abs_X2_w2, f_abs_X3_w2, unconvolved_w2)
        p.BasicPlot(param, param["display"]["window3"]["name"], W, F, E, l3, f_abs_con_w3, f_abs_ism_w3, f_abs_bp_w3, f_abs_X1_w3, f_abs_X2_w3, f_abs_X3_w3, unconvolved_w3)

    # Compute the goodness of fit.
    X = [F1, E1, f_abs_int_w1]
    print "Chi2:\t",s.chi2(X)
    print "DOF:\t\t",len(W1)-len(P)
    print "Chi2 reduced:\t",s.chi2(X)/(len(F1)-len(P)),"\n"



if __name__ == '__main__':
    main()