#!/usr/bin/env python
import numpy as np
import json, sys

from src.calculations import Calc
from src.statistics import Stats
from src.model import Model
from src.mcmc import MCMC

c   = Calc()
s   = Stats()
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

    # Need to change later
    #F = F*1e-7

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
                param["fit"]["ISM"]["T"],
                
                # Fixed CS parameters
                param["fit"]["CS"]["RV"],
                param["fit"]["CS"]["T"],
                
                # Fixed exocomet parameters
                param["fit"]["exocomet"]["T"]]

    Const   =   np.concatenate((ConstA,ConstB))
    
                # Free ISM parameters
    Par     =   [param["fit"]["ISM"]["log(N)"],
                param["fit"]["ISM"]["log(S)"],
                param["fit"]["ISM"]["T"],
                param["fit"]["ISM"]["xi"],                
                
                # Free CS parameters
                param["fit"]["CS"]["log(N)"],
                param["fit"]["CS"]["log(S)"],
                param["fit"]["CS"]["T"],
                param["fit"]["CS"]["xi"],
                
                # Free exocomet parameters
                param["fit"]["exocomet"]["log(N)"],
                param["fit"]["exocomet"]["log(S)"],
                param["fit"]["exocomet"]["T"],
                param["fit"]["exocomet"]["xi"],
                param["fit"]["exocomet"]["RV"]]

    X = F1, E1, m.Model(Par,Const,ModelType,param)[0]

    step = np.array([0.05,0.0,100.0,1.0,0.05,0.05,100,1.0,1.0,0.05,0.05,100,1])

    chain, moves = mc.McMC(W,X,m.Model, ModelType, param, Par, Const, step,1e4)
    
    outfile = 'chains/chain_B_'+sys.argv[1]
    np.savez(outfile, nN_ISM = chain[:,0], T_ISM = chain[:,1], xi_ISM = chain[:,2],\
        nN_CS = chain[:,3], nS_CS = chain[:,4], T_CS = chain[:,5], xi_CS = chain[:,6],\
        RV_X = chain[:,7], nN_X = chain[:,8], nS_X = chain[:,9], T_X = chain[:,10],\
        xi_X = chain[:,11])

    Pout = chain[moves,:]
    P_plot1 = [0,1]
    P_plot2 = [2,3]
    P_plot3 = [4,5]
    P_plot4 = [6,7]
    P_plot5 = [8,9]

    PU1 = mc.Median_and_Uncertainties(P_plot1,step,chain)
    PU2 = mc.Median_and_Uncertainties(P_plot2,step,chain)
    PU3 = mc.Median_and_Uncertainties(P_plot3,step,chain)
    PU4 = mc.Median_and_Uncertainties(P_plot4,step,chain)
    PU5 = mc.Median_and_Uncertainties(P_plot5,step,chain)
    
    print "log(N(N))_ISM\t=\t"  ,PU1[0][0],"\t+",PU1[1][0],"\t-",PU1[2][0]
    print "log(N(S))_ISM\t=\t"  ,PU1[0][1],"\t+",PU1[1][1],"\t-",PU1[2][1]
    print "xi_ISM\t\t=\t"       ,PU2[0][0],"\t+",PU2[1][0],"\t-",PU2[2][0]
    print "log(N(N))_CS\t=\t"   ,PU2[0][1],"\t+",PU2[1][1],"\t-",PU2[2][1]
    print "log(N(S))_CS\t=\t"   ,PU3[0][0],"\t+",PU3[1][0],"\t-",PU3[2][0]
    print "xi_CS\t\t=\t"        ,PU3[0][1],"\t+",PU3[1][1],"\t-",PU3[2][1]
    print "log(N(N))_X\t=\t"    ,PU4[0][0],"\t+",PU4[1][0],"\t-",PU4[2][0]
    print "log(N(S))_X\t=\t"    ,PU4[0][1],"\t+",PU4[1][1],"\t-",PU4[2][1]
    print "RV_X\t\t=\t"         ,PU5[0][0],"\t+",PU5[1][0],"\t-",PU5[2][0]
    print "xi_X\t\t=\t"         ,PU5[0][1],"\t+",PU5[1][1],"\t-",PU5[2][1]

if __name__ == '__main__':
    main()