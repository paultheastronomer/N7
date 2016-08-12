import numpy as np
import matplotlib.pyplot as plt
import json

from scipy.interpolate import interp1d

def Initialise():
    with open('../params.json') as param_file:    
		param = json.load(param_file)
    return param

def Bin_data(x,y1,e1,bin_pnts):
    bin_size    = int(len(x)/bin_pnts)
    bins        = np.linspace(x[0], x[-1], bin_size)
    digitized   = np.digitize(x, bins)
    bin_y       = np.array([y1[digitized == i].mean() for i in range(0, len(bins))])
    bin_e       = np.array([e1[digitized == i].mean() for i in range(0, len(bins))])
    return bins, bin_y, bin_e/np.sqrt(bin_pnts)
    
def main():    
    # Read all parameters from params.json file.
    param           = Initialise()

    # Define the data directory
    dat_directory   = param["directories"]["workdir"]

    # Load 2014 data
    W, RV, f0_0, e0_0, f_AG_0, e_AG_0   = np.genfromtxt(dat_directory+'B_N_2014.dat',unpack=True)
    
    # Load the data with little airglow
    W1, F1, E1  = np.genfromtxt(dat_directory+'N_AG_not_subtracted_2016_08_08_10Dec.txt',unpack=True)
    W2, F2, E2  = np.genfromtxt(dat_directory+'N_AG_not_subtracted_2016_08_08_26Dec.txt',unpack=True)
    W3, F3, E3  = np.genfromtxt(dat_directory+'N_AG_not_subtracted_2016_08_08_30Jan.txt',unpack=True)
    
    # Load all combined data (except 2014 data)
    Wc, Fc, Ec  = np.genfromtxt(dat_directory+'N_AG_not_subtracted_2016_08_08.txt',unpack=True)
    #print Wc
    #Wc  = Wc[12200:14000]
    #Fc  = Fc[12200:14000]
    #Ec  = Ec[12200:14000]
    
    #Fx, Fy  = np.genfromtxt('/home/paw/science/betapic/data/stellar/t08000g45p00k2.flx',unpack=True,skip_header=3)
    #Wx  = np.genfromtxt('/home/paw/science/betapic/data/stellar/UVBLUE_wavelengths.dat',unpack=True)
    Wx, Fx      = np.genfromtxt('/home/paw/science/betapic/data/stellar/N_RotB132L09.dat',unpack=True)
    
    # Load stellar model convolved with beta Pic broadening and shifted by +20.5 km/s (RV of Beta Pic).
    #Ws130, Fs130      = np.genfromtxt('/home/paw/science/betapic/data/stellar/BP_Phoenix_g4.5_v20.5_vsin130.dat',unpack=True)
    #Ws80, Fs80      = np.genfromtxt('/home/paw/science/betapic/data/stellar/BP_Phoenix_g4.5_v20.5_vsin80_shift0.045.dat',unpack=True)
    

    
    # Find the normalisation number used to show the divided spectra next to measured data.
    norm = np.median(Fc[7720:7850])
    
    # Configuring the plotting window
    #fig = plt.figure(figsize=(6.5,4.5))
    fig = plt.figure(figsize=(11,8))
    fontlabel_size  = 18
    tick_size       = 18
    params = {'backend': 'wxAgg', 'lines.markersize' : 2, 'axes.labelsize': fontlabel_size, 'font.size': fontlabel_size, 'legend.fontsize': 15, 'xtick.labelsize': tick_size, 'ytick.labelsize': tick_size, 'text.usetex': True}
    plt.rcParams.update(params)
    plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['text.latex.unicode'] = True
    
    #plt.step(Wx,Fx*1e-17,lw=2., color="red")

    
    #bin_pnts    = 3
    
    fS         = np.interp(Wc, Wx, Fx)
    #factor80    = np.median(Fc[7720:7850])/np.median(f80[7720:7850])
    norm      = np.median(Fc[7720:7850])
    DIV       = Fc/fS
    norm      = np.median(Fc[7720:7850])/np.median(DIV)
    DIV       = DIV*norm
    #DIV130      = norm*(Fc/f130)
    
    # Calculate the factor between the stellar model with data
    #factor80 = np.median(Fc[7720:7850]/Fs80[7720:7850])
    #factor130 = np.median(Fc[7720:7850]/Fs130[7720:7850])
    
    # Normalise the stellar model
    #Fs80    = Fs80*factor130
    #Fs130   = Fs130*factor130
    
        
    # Binning the combined data
    #Wcb, Fcb, Ecb           = Bin_data(Wc,Fc,Ec,bin_pnts)
    #Wcb, DIV80b, Ecb    = Bin_data(Wc,DIV80,Ec,bin_pnts)
    
    # Uncomment text below to see each epoch of data.
    '''
    w1b, f1b, e1b   = Bin_data(W,F1,E1,bin_pnts)
    w2b, f2b, e2b   = Bin_data(W,F2,E2,bin_pnts)
    w3b, f3b, e3b   = Bin_data(W,F3,E3,bin_pnts)
    plt.step(w1b, f1b, color="#FF9303", label=r'Method 1')
    plt.step(w2b, f2b, color="#0386FF", label=r'Method 1')
    plt.step(w3b, f3b, color="#00B233", label=r'Method 1')
    '''
    
    plt.step(Wc,DIV,lw=2., color="red")
    #plt.plot(Wc,f80,lw=2., color="red")
    plt.step(Wc,Fc,lw=2., color="black")

    
    #plt.step(Ws80,Fs80,lw=2., color="orange")
    #plt.step(Ws130,Fs130,lw=2., color="red")
    #plt.step(Wcb,Fcb,lw=2., color="black")
    
    

    #plt.step(Wcb,DIV80b,lw=2., color="orange")
    #plt.step(Wcb,DIV130b,lw=2., color="red")
    
    np.savetxt(dat_directory+"N_DIVB132L09.dat",np.column_stack((Wc, DIV, Ec)))
    
    #'''
    #plt.step(Wc[7720:7850],Fc[7720:7850],lw=2., color="blue")
    
    plt.xlim(1198.3,1202.3)
    plt.ylim(-1.0e-16,1.6e-14)
    plt.show()  
    
if __name__ == '__main__':
    main()
