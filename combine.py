import pyfits, os
import matplotlib.pyplot as plt
import numpy as np
import json, sys

from src.calculations import Calc
from src.plotting import Plotting
from src.model import Model
c   = Calc()
p   = Plotting()
m   = Model()

def Initialise():
    with open('params.json') as param_file:    
        param = json.load(param_file)
    return param

def main():

    # Read all parameters from params.json file.
    param           = Initialise()

    part            = param["BetaPictoris"]["part"]   # A = red, B = blue
    bin_size        = param["display"]["bin"]
    Line            = param["lines"]["line"]["N1"]["Wavelength"]
    RV_BP           = param["BetaPictoris"]["RV"]


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

    ratio1  = f0_bin/f1_bin
    ratio2  = f0_bin/f2_bin
    ratio3  = f0_bin/f3_bin

    # Define new start and stop values
    # which determine the region to be
    # cross correlated.
    if part == 'A':    
        x_lim1          = param["display"]["overview"]["Ax1"]
        x_lim2          = param["display"]["overview"]["Ax2"]   
        y_lim1          = param["display"]["overview"]["Ay1"]
        y_lim2          = param["display"]["overview"]["Ay2"]  
        n1              = param["regions"]["norm"]["A1"]
        n2              = param["regions"]["norm"]["A2"]
        s1              = param["regions"]["align"]["A1"]
        s2              = param["regions"]["align"]["A2"]
    else:
        x_lim1          = param["display"]["overview"]["Bx1"]
        x_lim2          = param["display"]["overview"]["Bx2"]   
        y_lim1          = param["display"]["overview"]["By1"]
        y_lim2          = param["display"]["overview"]["By2"]  
        n1              = param["regions"]["norm"]["B1"]    # "B1"        : 7700,
        n2              = param["regions"]["norm"]["B2"]    # "B2"        : 8050
        s1              = param["regions"]["align"]["B1"]
        s2              = param["regions"]["align"]["B2"]
        #s1 = 10000     # Original region
        #s2  = 13473
        #s1 = 12550  # Region around a SII line
        #s2 = 12750
        #s1 = 7200  # Region just left of the NI line
        #s2 = 7400
        #s1 = 4800
        #s2 = 5020
        #s1 = 7205#7450#7205  # Region including NI triplet
        #s2 = 7647

    # Creates a figure showing the quiet regions of the spectra.
    #p.OverviewPlot(part, x_lim1, x_lim2, y_lim1, y_lim2, w0_0, w0_bin, ratio1, ratio2, ratio3, f0_bin, f1_bin, f2_bin, f3_bin, s1, s2)
    #sys.exit()

    # --------------------------------------------------------------------------
    # This is where the the quiet region of the spectra
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


    #'''
    # Save the shifted from each visit into .dat file
    #np.savetxt(param["directories"]["workdir"]+"V1.dat",np.column_stack((W, f0_0, e0_0, f_AG_0, e_AG_0)))
    #np.savetxt(param["directories"]["workdir"]+"V2.dat",np.column_stack((W, F0_1, E0_1, F1_1, E1_1, F2_1, E2_1, AG1, AG1err, F_ave_w_1, E_ave_w_1)))
    #np.savetxt(param["directories"]["workdir"]+"V3.dat",np.column_stack((W, F0_2, E0_2, F1_2, E1_2, F2_2, E2_2, F3_2, E3_2, AG2, AG2err, F_ave_w_2, E_ave_w_2)))
    #np.savetxt(param["directories"]["workdir"]+"V4.dat",np.column_stack((W, F0_3, E0_3, F1_3, E1_3, F2_3, E2_3, F3_3, E3_3, AG3, AG3err, F_ave_w_3, E_ave_w_3)))

    # For a better naming convention we will rename the following variables
    F0_0    = f0_0
    E0_0    = e0_0

    # Calculate the Correction Factor. Find which factor to multiply the spectra by to match the F0_0 flux level.
    # First define the region [n1:n2] where the flux from the individual spectra should be compared


    # To see the chosen region uncomment the lines below.
    '''
    plt.plot(W,F0_0)
    plt.plot(W,F0_1)
    plt.plot(W[n1:n2],F0_1[n1:n2])
    plt.show()
    sys.exit()
    '''

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

    # All data of the debris disk combined
    #Flux = np.array([F0_0,Fc[0],Fc[1],Fc[2],Fc[3],Fc[4],Fc[5],Fc[6],Fc[7],Fc[8],Fc[9],Fc[10]])
    #Err  = np.array([E0_0,Ec[0],Ec[1],Ec[2],Ec[3],Ec[4],Ec[5],Ec[6],Ec[7],Ec[8],Ec[9],Ec[10]])

    # All data combined except the 2014 data due to large amounts of airglow
    Flux_no_2014 = np.array([Fc[0],Fc[1],Fc[2],Fc[3],Fc[4],Fc[5],Fc[6],Fc[7],Fc[8],Fc[9],Fc[10]])
    Err_no_2014  = np.array([Ec[0],Ec[1],Ec[2],Ec[3],Ec[4],Ec[5],Ec[6],Ec[7],Ec[8],Ec[9],Ec[10]])
    
    # For 10 Dec data uncomment the two lines below
    Flux2 = np.array([Fc[0],Fc[1],Fc[2]])
    Err2  = np.array([Ec[0],Ec[1],Ec[2]])

    # For 24 Dec data uncomment the two lines below
    Flux3 = np.array([Fc[3],Fc[4],Fc[5],Fc[6]])
    Err3  = np.array([Ec[3],Ec[4],Ec[5],Ec[6]])

    # For 30 Jan data uncomment the two lines below
    Flux4 = np.array([Fc[7],Fc[8],Fc[9],Fc[10]])
    Err4  = np.array([Ec[7],Ec[8],Ec[9],Ec[10]])
    
    F_tot, F_tot_err = c.WeightedAvg(Flux_no_2014,Err_no_2014)

    F1, E1 =  F0_0, E0_0
    F2, E2 =  c.WeightedAvg(Flux2,Err2)
    F3, E3 =  c.WeightedAvg(Flux3,Err3)
    F4, E4 =  c.WeightedAvg(Flux4,Err4)
    
    #F_tot, F_tot_err    =  Fc[3], Ec[3]
    #############################################################################################

    W_bin, F_tot_bin, E_tot_bin =   c.BinData(w0_0, F_tot, F_tot_err,bin_size)
    W_bin, F1_bin, E1_bin       =   c.BinData(w0_0,F1,E2,bin_size)
    W_bin, F2_bin, E2_bin       =   c.BinData(w0_0,F2,E2,bin_size)
    W_bin, F3_bin, E3_bin       =   c.BinData(w0_0,F3,E3,bin_size)
    W_bin, F4_bin, E4_bin       =   c.BinData(w0_0,F4,E4,bin_size)

    # Save the combined spectrum
    #np.savetxt(param["directories"]["workdir"]+"NI_2016_11_07.txt",np.column_stack((w0_0, F_tot, F_tot_err)))
    #np.savetxt(param["directories"]["workdir"]+"2016_11_17_2014_B.txt",np.column_stack((w0_0, F1, E1)))
    #sys.exit()
    #'''
    
    # Temp code
    #w0_0, F_tot, F_tot_err = np.genfromtxt(param["directories"]["workdir"]+'NI_2016_10_12.txt',unpack=True)

    #p.CombinedPlot(param, "window1", W, F_tot, W_bin, F2_bin, F3_bin, F4_bin)

    # Stellar code part
    Wx      = np.genfromtxt(param["directories"]["stellardir"]+'UVBlue/wl.t08052g415p004.dat',unpack=True,usecols=0)
    H_nu    = np.genfromtxt(param["directories"]["stellardir"]+'UVBlue/sp.t08052g415p004.dat',unpack=True,skip_header=3,usecols=0)
    c_light = 299792458
    Fx      = H_nu*4*np.pi*c_light/Wx**2

    # To speed up the code, only broaden region around NI
    #Wx = Wx[13650:14000]   
    #Fx = Fx[13650:14000]

    #Shift to BetaPic RV
    Wx=Wx*(1.+param["BetaPictoris"]["RV"]/3e5)

    #Apply shift
    #for i in range(len(Wx)):
    #    Wx[i] = Wx[i]+0.11

    # Load broadening parameters
    eps             = param["BetaPictoris"]["eps"]
    vsini           = param["BetaPictoris"]["vsini"]

    # Interpolate synthetic spectrum onto the wavelength scale of the COS data
    Fint            = np.interp(w0_0,Wx,Fx)

    #s1 = 10000
    #s2  = 13473

    # Uncomment below to check shift in synthetic spectrum
    # relative to COS spectrum.
    #c.ShiftSpec(F_tot,Fint,F_tot_err,w0_0,s1,s2,Line)

    # Rotationally broaden spectrum
    RotBroadSpec    = c.RotBroad(w0_0, Fint, eps, vsini)

    # Perform LSF calculation
    kernel1         = m.LSF(param["lines"]["line"]["N1"]["Wavelength"], w0_0)
    Fx_con          = np.convolve(RotBroadSpec, kernel1, mode='valid')[:-1]

    synth_norm_region = []
    norm_region       = []
    for i in range(len(w0_0)):
        if 1198.5 < w0_0[i] < 1199.4:
            norm_region.append(F_tot[i])
            synth_norm_region.append(RotBroadSpec[i])
    norm_region         = np.array(norm_region)
    synth_norm_region   = np.array(synth_norm_region)
    median_synth_region = np.median(synth_norm_region)
    flux_ratio          = np.median(norm_region)/median_synth_region
 
    #np.savetxt(param["directories"]["workdir"]+"div_by_UVBLUE_2016_11_02.dat",np.column_stack((w0_0, F_tot/(RotBroadSpec/median_synth_region), F_tot_err)))

    fig = plt.figure(figsize=(8,5))
    #fig = plt.figure(figsize=(10,14))
    #fig = plt.figure(figsize=(14,10))

    fontlabel_size  = 18
    tick_size       = 18
    params = {'backend': 'wxAgg', 'lines.markersize' : 2, 'axes.labelsize': fontlabel_size, 'font.size': fontlabel_size, 'legend.fontsize': 15, 'xtick.labelsize': tick_size, 'ytick.labelsize': tick_size, 'text.usetex': True}
    plt.rcParams.update(params)
    plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['text.latex.unicode'] = True

    plt.step(w0_0, F_tot,color='black',lw=1.3)
    #plt.step(w0_0,AG1,lw=1.2,color="#FF9303",label='2015v1')
    #plt.step(w0_0,AG2,lw=1.2,color="#0386FF",label='2015v2')
    #plt.step(w0_0,AG3,lw=1.2,color="#00B233",label='2016v3')
    
    plt.step(w0_0,F_tot/(RotBroadSpec/median_synth_region),color='blue',lw=1)
    plt.plot(w0_0,RotBroadSpec*flux_ratio,color='red',lw=2)
    #plt.step(w0_0, 0.85e3*(F_tot/RotBroadSpec),color='black',lw=1.3)
    plt.xlim(1198.5,1201.5)
    plt.ylim(0,2e-14)
    fig.tight_layout()
    x = [1199,1199.5,1200,1200.5,1201,1201.5,1202]
    labels = ['1199.0','1199.5','1200.0','1200.5','1201','1201.5','1202.0']
    plt.xticks(x, labels)
    plt.xlabel(r'Wavelength [\AA]')
    #plt.savefig('spectra_comparison.pdf', bbox_inches='tight', pad_inches=0.1,dpi=300)
    plt.show()

if __name__ == '__main__':
    main()