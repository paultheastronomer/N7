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
    
    Par     =   [
                # GEO parameters which can be set free
                param["fit"]["GEO"]["log(N)"],
                param["fit"]["GEO"]["log(S)"],
                param["fit"]["GEO"]["b"],
                param["fit"]["GEO"]["T"],
                param["fit"]["GEO"]["xi"],
                param["fit"]["GEO"]["RV"],

                # ISM parameters which can be set free
                param["fit"]["ISM"]["log(N)"],
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


    X = F1, E1, m.Model(Par,Const,ModelType,param)[0]

    #step = np.array([0.1, 0.0, 0.0, 0.0, 0.5, 0.0,    0.25, 0.0, 0.0, 0.0, 0.0,    0.25, 0.0, 0.0, 0.0, 0.1,    0.25, 0.0, 0.0, 0.0, 0.1, 1.0,     0.25, 0.0, 0.0, 0.0, 1.0, 1.0])
    #step = np.array([0.1, 0.0, 0.0, 0.0, 0.1, 0.0,    0.15, 0.0, 0.0, 0.0, 0.0,    0.09, 0.0, 0.0, 0.0, 0.1,    0.10, 0.0, 0.0, 0.0, 0.1, 1.0,     0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    step = np.array([0.3, 0.0, 0.0, 0.0, 0.3, 0.0,    0.3, 0.0, 0.0, 0.0, 0.0,    0.05, 0.0, 0.0, 0.0, 0.3,    0.3, 0.0, 0.0, 0.0, 0.3, 1.0,     0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    chain, moves, stats = mc.McMC(W,X,m.Model, ModelType, param, Par, Const, step,param["fit"]["MCMC"]["number_of_steps"])

    name = param["fit"]["MCMC"]["chain_name"]

    text_file = open('chains/'+name+'_'+sys.argv[1]+'moves.txt', 'w+')
    print >> text_file, moves, "Accepted steps: ",round(100.*(moves/param["fit"]["MCMC"]["number_of_steps"]),2),"%"
    text_file.close()
    
    outfile = 'chains/'+name+'_'+sys.argv[1]
    np.savez(outfile,
        nN_GEO = chain[:,0], nS_GEO = chain[:,1], b_GEO = chain[:,2], T_GEO = chain[:,3], xi_GEO = chain[:,4], RV_GEO = chain[:,5],\
        
        nN_ISM = chain[:,6], nS_ISM = chain[:,7], b_ISM = chain[:,8], T_ISM = chain[:,9], xi_ISM = chain[:,10],\
        nN_CS = chain[:,11], nS_CS = chain[:,12],  b_CS = chain[:,13], T_CS = chain[:,14], xi_CS = chain[:,15],\
        
        nN_X1 = chain[:,16], nS_X1 = chain[:,17],  b_X1 = chain[:,18], T_X1 = chain[:,19],xi_X1 = chain[:,20],RV_X1 = chain[:,21],
        nN_X2 = chain[:,22], nS_X2 = chain[:,23],  b_X2 = chain[:,24], T_X2 = chain[:,25],xi_X2 = chain[:,26],RV_X2 = chain[:,27])
    statsfile = 'chains/Stat_'+name+'_'+sys.argv[1]
    np.savez(statsfile, Chi2 = stats[:,0])

if __name__ == '__main__':
    main()