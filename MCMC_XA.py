#!/usr/bin/env python
import numpy as np
import json, sys

from src.statistics import Stats
from src.model import Model
from src.mcmc import MCMC

s   = Stats()
m   = Model()
mc  = MCMC() 

def Initialise():
    with open('params.json') as param_file:    
        param = json.load(param_file)
    return param

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
    dat_directory   = param["directories"]["exoatmdir"]
    
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

    X = F1, E1, m.Model(Par,Const,ModelType,param)[0]

    step = np.array([0.1,0.0,0.1,0.0,0.06,1.0,0.0])
    chain, moves = mc.McMC(W,X,m.Model, ModelType, param, Par, Const, step,1e5)
    
    outfile = 'chains/chain_R_'+sys.argv[1]
    np.savez(outfile, nh_ISM = chain[:,0], b_ISM = chain[:,1], nh_CS = chain[:,2], b_CS = chain[:,3], nh_X = chain[:,4], RV_X = chain[:,5], b_X = chain[:,6])

    Pout = chain[moves,:]
    P_plot1 = [0,1]
    P_plot2 = [2,3]
    P_plot3 = [4,5]

    PU1 = mc.Median_and_Uncertainties(P_plot1,step,chain)
    PU2 = mc.Median_and_Uncertainties(P_plot2,step,chain)
    PU3 = mc.Median_and_Uncertainties(P_plot3,step,chain)
    
    print "log(N(H))_CS\t=\t" ,PU1[0][0],"\t+",PU1[1][0],"\t-",PU1[2][0]
    print "b_ISM\t\t=\t"      ,PU1[0][1],"\t+",PU1[1][1],"\t-",PU1[2][1]
    print "b_CS\t\t=\t"       ,PU2[0][0],"\t+",PU2[1][0],"\t-",PU2[2][0]
    print "log(N(H))_X\t=\t"  ,PU2[0][1],"\t+",PU2[1][1],"\t-",PU2[2][1]
    print "RV_X\t\t=\t"       ,PU3[0][0],"\t+",PU3[1][0],"\t-",PU3[2][0]
    print "b_X\t\t=\t"        ,PU3[0][1],"\t+",PU3[1][1],"\t-",PU3[2][1]

if __name__ == '__main__':
    main()