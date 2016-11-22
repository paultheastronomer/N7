#!/usr/bin/env python
import numpy as np
import json, sys

from src.calculations import Calc
from src.model import Model
from src.mcmc import MCMC

c   = Calc()
m   = Model()
mc  = MCMC() 

def Initialise():
    with open('params.json') as param_file:    
        param = json.load(param_file)
    return param

def main():    

    # Read all parameters from params.json file.
    param           = Initialise()
    
    # Define the data directory
    dat_directory   = param["directories"]["exoatmdir"]
    
    # Select the model type
    ModelType       = 1#param["fit"]["ModelType"]
    
    Nwindows        = param["fit"]["windows"]["number"]

    # Load the data file
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
                param["fit"]["exocomet"]["log(N)"],
                param["fit"]["exocomet"]["log(S)"],
                param["fit"]["exocomet"]["b"],
                param["fit"]["exocomet"]["T"],
                param["fit"]["exocomet"]["xi"],
                param["fit"]["exocomet"]["RV"]]

    X = F1, E1, m.Model(Par,Const,ModelType,param)[0]

    step = np.array([0.1, 0.0, 0.0, 0.0, 0.5,    0.2, 0.0, 0.0, 0.0, 0.2,    0.05, 0.0, 0.0, 0.0, 1.0, 2.0])
    #step = np.array([0.2, 0.0, 0.0, 0.0, 0.1,    0.05, 0.0, 0.0, 0.0, 0.1,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    chain, moves = mc.McMC(W,X,m.Model, ModelType, param, Par, Const, step,1.0e5)
    
    outfile = 'chains/chain_b7fast_'+sys.argv[1]
    np.savez(outfile, nN_ISM = chain[:,0], nS_ISM = chain[:,1], b_ISM = chain[:,2], T_ISM = chain[:,3], xi_ISM = chain[:,4],\
        nN_CS = chain[:,5], nS_CS = chain[:,6],  b_CS = chain[:,7], T_CS = chain[:,8], xi_CS = chain[:,9],\
        nN_X = chain[:,10], nS_X = chain[:,11],  b_X = chain[:,12], T_X = chain[:,13],xi_X = chain[:,14],RV_X = chain[:,15])

if __name__ == '__main__':
    main()