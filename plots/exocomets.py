import numpy as np
import matplotlib.pyplot as plt
import json


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

    wc, fc, ec      = np.genfromtxt(dat_directory+'N_no2014_2016_07_27.txt',unpack=True)    
    
    w1, f1, e1      = np.genfromtxt(dat_directory+'N_2015v1_2016_07_27.txt',unpack=True)
    w2, f2, e2      = np.genfromtxt(dat_directory+'N_2015v2_2016_07_27.txt',unpack=True)
    w3, f3, e3      = np.genfromtxt(dat_directory+'N_2016v1_2016_07_27.txt',unpack=True)
    
    W, FitNoPS, F, E, C, Fit, N1, N2, N3      = np.genfromtxt(dat_directory+'NI_fit_2016_08_01.dat',unpack=True,skip_header=3)
    
    
    #fig = plt.figure(figsize=(6.5,4.5))
    fig = plt.figure(figsize=(11,8))
    fontlabel_size  = 18
    tick_size       = 18
    params = {'backend': 'wxAgg', 'lines.markersize' : 2, 'axes.labelsize': fontlabel_size, 'font.size': fontlabel_size, 'legend.fontsize': 15, 'xtick.labelsize': tick_size, 'ytick.labelsize': tick_size, 'text.usetex': True}
    plt.rcParams.update(params)
    plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['text.latex.unicode'] = True
    
    bin_pnts = 3
    
    wcb, fcb, ecb   = Bin_data(wc,fc,ec,bin_pnts)
    
    w1b, f1b, e1b   = Bin_data(w1,f1,e1,bin_pnts)
    w2b, f2b, e2b   = Bin_data(w2,f2,e2,bin_pnts)
    w3b, f3b, e3b   = Bin_data(w3,f3,e3,bin_pnts)
    
    Wb, Fb, Diffb    = Bin_data(W,F,F-Fit,bin_pnts)

    
    #plt.plot(W,N1,lw=1,color="blue")
    #plt.plot(W,N2,lw=1,color="green")
    #plt.plot(W,Fit,lw=2,color="red")
    
    #plt.step(wc, fc, lw=2, color="black", label=r'Method 1')
    #plt.step(W, F, lw=2, color="green", label=r'Method 1')
    
    plt.plot(W,Fit,lw=2., color="red")
    plt.step(Wb,Fb,lw=2., color="black")
    
    plt.plot(W,F-F,lw=2,color="gray")
    plt.step(Wb,Diffb,lw=2., color="blue")
    
    #plt.step(w1b, f1b, label=r'Method 1')
    #plt.step(w2b, f2b, label=r'Method 1')
    #plt.step(w3b, f3b, label=r'Method 1')
    plt.xlim(1198.3,1202.3)
    plt.ylim(-2.5e-14,2.2e-14)
    plt.show()  
    
if __name__ == '__main__':
    main()
