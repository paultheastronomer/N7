import pyfits, os
import matplotlib.pyplot as plt
import numpy as np
import json, sys

from src.calculations import Calc
from src.plotting import Plotting
c   = Calc()
p   = Plotting()

def Initialise():
    with open('params.json') as param_file:    
        param = json.load(param_file)
    return param

def main():

    # Read all parameters from params.json file.
    param           = Initialise()


    dat_directory   = param["directories"]["workdir"]
    stellar_dir     = param["directories"]["stellardir"]
    part            = param["BetaPictoris"]["RV"]   # A = red, B = blue
    bin_size        = param["display"]["bin"]
    x_lim1          = param["display"]["overview"]["Ax1"]
    x_lim2          = param["display"]["overview"]["Ax2"]   
    Line            = param["lines"]["line"]["N1"]["Wavelength"]
    RV_BP           = param["BetaPictoris"]["RV"]

    
    #c.BroadenSpec(1180, 1210, 1, 80)
    #sys.exit()

    # Load data visit 1 2014
    fits_location   = param["directories"]["2014"]
    w0_0,f0_0,e0_0,w_AG_0,f_AG_0,e_AG_0,NumFits_0                                               = c.GetData(fits_location,part)

    # Load data visit 1 2015
    fits_location   = param["directories"]["2015a"]
    w0_1,w1_1,w2_1,f0_1,f1_1,f2_1,e0_1,e1_1,e2_1,w_AG_1,f_AG_1,e_AG_1,NumFits_1                 = c.GetData(fits_location,part)
    
    # Load data visit 2 2015
    fits_location   = param["directories"]["2015b"]
    w0_2,w1_2,w2_2,w3_2,f0_2,f1_2,f2_2,f3_2,e0_2,e1_2,e2_2,e3_2,w_AG_2,f_AG_2,e_AG_2,NumFits_2  = c.GetData(fits_location,part)

    # Load data visit 3 2016
    fits_location   = param["directories"]["2016"]
    w0_3,w1_3,w2_3,w3_3,f0_3,f1_3,f2_3,f3_3,e0_3,e1_3,e2_3,e3_3,w_AG_3,f_AG_3,e_AG_3,NumFits_3  = c.GetData(fits_location,part)    
    
    # The next part of the code is aimed at finding quiet regions in the spectra
    # not affected by FEB activity. The binning of the data is done so that
    # the quiet regions can be better seen.
    
    # Bin the data by bin_size (defined under "display" in the params.json file).
    w0_bin, f0_bin, e0_bin   =  c.BinData(w0_0,f0_0,e0_0,bin_size)
    w1_bin, f1_bin, e1_bin   =  c.BinData(w0_1,f0_1,e0_1,bin_size)
    w2_bin, f2_bin, e2_bin   =  c.BinData(w0_2,f0_2,e0_2,bin_size)
    w3_bin, f3_bin, e3_bin   =  c.BinData(w0_3,f0_3,e0_3,bin_size)
    
    # To avoid dividing my 0 when calculating the ratios below, NaN and 0's are
    # replaced by the median. This is only for visual purposes.
    f0_bin = c.ReplaceWithMedian(f0_bin)
    f1_bin = c.ReplaceWithMedian(f1_bin)
    f2_bin = c.ReplaceWithMedian(f2_bin)
    f3_bin = c.ReplaceWithMedian(f3_bin)
    
    #AG0_bin = c.ReplaceWithMedian(AG0_bin)
    #AG1_bin = c.ReplaceWithMedian(AG1_bin)
    #AG2_bin = c.ReplaceWithMedian(AG2_bin)
    #AG3_bin = c.ReplaceWithMedian(AG3_bin)

    ratio1  = f0_bin/f1_bin
    ratio2  = f0_bin/f2_bin
    ratio3  = f0_bin/f3_bin

    # Define new start and stop values
    # which determine the region to be
    # cross correlated.
    if part == 'A':    
        s1 = 6100
        s2  = 11400  
    else:
        s1 = 10000
        s2  = 13473
        #s1 = 12550  # Region around a SII line
        #s2 = 12750
        #s1 = 7200  # Region just left of the NI line
        #s2 = 7400



    # Creates a figure showing the quiet regions of the spectra.
    #p.OverviewPlot(part, x_lim1, x_lim2, w0_0, w0_bin, ratio1, ratio2, ratio3, f0_bin, f1_bin, f2_bin, f3_bin, s1, s2)
    #sys.exit()
    # --------------------------------------------------------------------------
    # The heart of the code. This is where the the quiet region of the spectra
    # are used to cross correlate the spectra and calculate the shift.
    # --------------------------------------------------------------------------

    # These are dummy arrays used for the 10 Dec 2015 observations
    f_empty = []
    e_empty = []
    
    print "\n\nShifting 10 Dec 2015 observations:"
    W, F0_1, E0_1, F1_1, E1_1, F2_1, E2_1, AG1, AG1err, F_ave_w_1, E_ave_w_1                = c.ExportShitedSpectra(w0_0,f0_0,f0_1,f1_1,f2_1,f_empty,f_AG_1,e0_0,e0_1,e1_1,e2_1,e_empty,e_AG_1,NumFits_1,s1,s2,Line)
    
    print "\n\nShifting 24 Dec 2015 observations:"
    W, F0_2, E0_2, F1_2, E1_2, F2_2, E2_2, F3_2, E3_2, AG2, AG2err, F_ave_w_2, E_ave_w_2    = c.ExportShitedSpectra(w0_0,f0_0,f0_2,f1_2,f2_2,f3_2,f_AG_2,e0_0,e0_2,e1_2,e2_2,e3_2,e_AG_2,NumFits_2,s1,s2,Line)
    
    print "\n\nShifting 30 Jan 2016 observations:"
    W, F0_3, E0_3, F1_3, E1_3, F2_3, E2_3, F3_3, E3_3, AG3, AG3err, F_ave_w_3, E_ave_w_3    = c.ExportShitedSpectra(w0_0,f0_0,f0_3,f1_3,f2_3,f3_3,f_AG_3,e0_0,e0_3,e1_3,e2_3,e3_3,e_AG_3,NumFits_3,s1,s2,Line)

    # Save data into .dat file
    #np.savetxt(dat_directory+"/shifted_spectra/"+"V1.dat",np.column_stack((W, f0_0, e0_0, f_AG_0, e_AG_0)))
    #np.savetxt(dat_directory+"/shifted_spectra/"+"V2.dat",np.column_stack((W, F0_1, E0_1, F1_1, E1_1, F2_1, E2_1, AG1, AG1err, F_ave_w_1, E_ave_w_1)))
    #np.savetxt(dat_directory+"/shifted_spectra/"+"V3.dat",np.column_stack((W, F0_2, E0_2, F1_2, E1_2, F2_2, E2_2, F3_2, E3_2, AG2, AG2err, F_ave_w_2, E_ave_w_2)))
    #np.savetxt(dat_directory+"/shifted_spectra/"+"V4.dat",np.column_stack((W, F0_3, E0_3, F1_3, E1_3, F2_3, E2_3, F3_3, E3_3, AG3, AG3err, F_ave_w_3, E_ave_w_3)))

    F0_0    = f0_0
    E0_0    = e0_0

    n1 = 7700
    n2 = 8050

    # Uncomment to see region
    '''
    plt.plot(W,F0_0)
    plt.plot(W,F0_1)
    plt.plot(W[n1:n2],F0_1[n1:n2])
    plt.show()
    sys.exit()
    '''

    # Calculate the Correction Factor
    C   = [c.CF(F0_1,E0_1,F0_0,E0_0,n1,n2),c.CF(F1_1,E1_1,F0_0,E0_0,n1,n2),c.CF(F2_1,E2_1,F0_0,E0_0,n1,n2),\
    c.CF(F0_2,E0_2,F0_0,E0_0,n1,n2),c.CF(F1_2,E1_2,F0_0,E0_0,n1,n2),c.CF(F2_2,E2_2,F0_0,E0_0,n1,n2),c.CF(F3_2,E3_2,F0_0,E0_0,n1,n2),\
    c.CF(F0_3,E0_3,F0_0,E0_0,n1,n2),c.CF(F1_3,E1_3,F0_0,E0_0,n1,n2),c.CF(F2_3,E2_3,F0_0,E0_0,n1,n2),c.CF(F3_3,E3_3,F0_0,E0_0,n1,n2)]

    F  = [F0_1,F1_1,F2_1,F0_2,F1_2,F2_2,F3_2,F0_3,F1_3,F2_3,F3_3]
    E  = [E0_1,E1_1,E2_1,E0_2,E1_2,E2_2,E3_2,E0_3,E1_3,E2_3,E3_3]
    

    Fc = [[] for _ in range(len(C))]
    Ec = [[] for _ in range(len(C))]

    for i in range(len(C)):
        Fc[i] = F[i]*C[i]   # Correct for lower efficiency
        Ec[i] = E[i]*C[i]   # accordingly correct the tabulated error bars

    # -0.8" Ly-alpha wing
    #############################################################################################    
    # For all data uncomment the two lines below
    Flux = np.array([Fc[1],Fc[4],Fc[8]])
    Err  = np.array([Ec[1],Ec[4],Ec[8]])
    
    # For 10 Dec data uncomment the two lines below
    #Flux = np.array([Fc[1]])
    #Err  = np.array([Ec[1]])

    # For 24 Dec data uncomment the two lines below
    #Flux = np.array([Fc[4]])
    #Err  = np.array([Ec[4]])

    # For 30 Jan data uncomment the two lines below
    #Flux = np.array([Fc[8]])
    #Err  = np.array([Ec[8]])
    
    F1, F1_err    =  c.WeightedAvg(Flux,Err)         

    #############################################################################################    

    # 0.8" Ly-alpha wing
    #############################################################################################
    # For all data uncomment the two lines below
    Flux = np.array([Fc[2],Fc[9]])
    Err  = np.array([Ec[2],Ec[9]])
    
    # Fc[5] shows strange "emission feature" not consistent with the 1.1" offset data
    # and has thus been removed. To include it uncomment the two lines below.
    #Flux = np.array([Fc[2],Fc[5],Fc[9]])
    #Err  = np.array([Ec[2],Ec[5],Ec[9]])

    # For 10 Dec data uncomment the two lines below
    #Flux = np.array([Fc[2]])
    #Err  = np.array([Ec[2]])

    # For 24 Dec data uncomment the two lines below
    #Flux = np.array([Fc[5]])
    #Err  = np.array([Ec[5]])

    # For 30 Jan data uncomment the two lines below
    #Flux = np.array([Fc[9]])
    #Err  = np.array([Ec[9]])

    F2, F2_err    =  c.WeightedAvg(Flux,Err)         
    #############################################################################################

    # 1.1" Ly-alpha wing
    #############################################################################################
    # For all data uncomment the two lines below
    Flux = np.array([Fc[6],Fc[10]])
    Err  = np.array([Ec[6],Ec[10]])

    # For 24 Dec data uncomment the two lines below    
    #Flux = np.array([Fc[6]])
    #Err  = np.array([Ec[6]])

    # For 30 Jan data uncomment the two lines below    
    #Flux = np.array([Fc[10]])
    #Err  = np.array([Ec[10]])
    
    F3, F3_err    =  c.WeightedAvg(Flux,Err)         
    #############################################################################################

    Flux = np.array([F0_0,Fc[0],Fc[1],Fc[2],Fc[3],Fc[4],Fc[5],Fc[6],Fc[7],Fc[8],Fc[9],Fc[10]])
    Err  = np.array([E0_0,Ec[0],Ec[1],Ec[2],Ec[3],Ec[4],Ec[5],Ec[6],Ec[7],Ec[8],Ec[9],Ec[10]])

    Flux_no_2014 = np.array([Fc[0],Fc[1],Fc[2],Fc[3],Fc[4],Fc[5],Fc[6],Fc[7],Fc[8],Fc[9],Fc[10]])
    Err_no_2014  = np.array([Ec[0],Ec[1],Ec[2],Ec[3],Ec[4],Ec[5],Ec[6],Ec[7],Ec[8],Ec[9],Ec[10]])


    # For 2014 data see line straight after "F_tot, F_tot_err"
    
    # For 10 Dec data uncomment the two lines below
    Flux2 = np.array([Fc[0],Fc[1],Fc[2]])
    Err2  = np.array([Ec[0],Ec[1],Ec[2]])

    # For 24 Dec data uncomment the two lines below
    Flux3 = np.array([Fc[3],Fc[4],Fc[5],Fc[6]])
    Err3  = np.array([Ec[3],Ec[4],Ec[5],Ec[6]])

    # For 30 Jan data uncomment the two lines below
    Flux4 = np.array([Fc[7],Fc[8],Fc[9],Fc[10]])
    Err4  = np.array([Ec[7],Ec[8],Ec[9],Ec[10]])
    
    F_tot, F_tot_err    =  c.WeightedAvg(Flux_no_2014,Err_no_2014)
    #F_tot_no_V3, F_tot_err_no_V3    =  c.WeightedAvg(Flux_no_V3,Err_no_V3)

    F1, E1    =  F0_0, E0_0
    F2, E2    =  c.WeightedAvg(Flux2,Err2)
    F3, E3    =  c.WeightedAvg(Flux3,Err3)
    F4, E4    =  c.WeightedAvg(Flux4,Err4)
    
    #F_tot, F_tot_err    =  F0_0, E0_0
    #############################################################################################

    W_bin, F_tot_bin, E_tot_bin = c.BinData(w0_0, F_tot, F_tot_err,bin_size)
    W_bin, F1_bin, E1_bin       =   c.BinData(w0_0,F1,E2,bin_size)
    W_bin, F2_bin, E2_bin       =   c.BinData(w0_0,F2,E2,bin_size)
    W_bin, F3_bin, E3_bin       =   c.BinData(w0_0,F3,E3,bin_size)
    W_bin, F4_bin, E4_bin       =   c.BinData(w0_0,F4,E4,bin_size)

    #np.savetxt(dat_directory+"/shifted_spectra/"+"NI_2016_10_05.txt",np.column_stack((w0_0, F_tot, F_tot_err)))
    
    W80, F80                      = np.genfromtxt(stellar_dir+param["files"]["stellar_spectrum1"],unpack=True)
    W130, F130                      = np.genfromtxt(stellar_dir+param["files"]["stellar_spectrum2"],unpack=True)

    F80i                         = np.interp(w0_0,W80,F80)
    F130i                         = np.interp(w0_0,W130,F130)

    W_dummy, F80i_bin, E_dummy    = c.BinData(w0_0, F80i, F_tot_err,bin_size)
    W_dummy, F130i_bin, E_dummy    = c.BinData(w0_0, F130i, F_tot_err,bin_size)


    p.CombinedPlot(param, param["display"]["window1"]["name"], w0_0, W_bin, F_tot, F_tot_err, F_tot_bin, F80i, F80i_bin, F130i, F130i_bin, F1_bin, F2_bin, F3_bin, F4_bin)

if __name__ == '__main__':
    main()
