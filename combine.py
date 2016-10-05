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
    part            = 'B'   # A = red, B = blue
    bin_size        = param["display"]["bin"]
    x_lim1          = param["display"]["overview"]["Ax1"]
    x_lim2          = param["display"]["overview"]["Ax2"]   
    LyA             = param["lines"]["line"]["N1"]["Wavelength"]
    RV_BP           = param["BetaPictoris"]["RV"]

    # Load data visit 1 2014
    fits_location   = param["directories"]["2014"]
    w0_0,f0_0,e0_0,w_AG_0,f_AG_0,e_AG_0,NumFits_0   = c.GetData(fits_location,part)

    # Load data visit 1 2015
    fits_location   = param["directories"]["2015a"]
    w0_1,w1_1,w2_1,f0_1,f1_1,f2_1,e0_1,e1_1,e2_1,w_AG_1,f_AG_1,e_AG_1,NumFits_1 = c.GetData(fits_location,part)
    
    # Load data visit 2 2015
    fits_location   = param["directories"]["2015b"]
    w0_2,w1_2,w2_2,w3_2,f0_2,f1_2,f2_2,f3_2,e0_2,e1_2,e2_2,e3_2,w_AG_2,f_AG_2,e_AG_2,NumFits_2 = c.GetData(fits_location,part)

    # Load data visit 3 2016
    fits_location   = param["directories"]["2016"]
    w0_3,w1_3,w2_3,w3_3,f0_3,F0_3,F1_3,F2_3,e0_3,e1_3,e2_3,e3_3,w_AG_3,f_AG_3,e_AG_3,NumFits_3 = c.GetData(fits_location,part)    
    
    # The next part of the code is aimed at finding quiet regions in the spectra
    # not affected by FEB activity. The binning of the data is done so that
    # the quiet regions can be better seen.
    
    # Bin the data by bin_size (defined under "display" in the params.json file).
    w_bin, y0_bin, y1_bin       = c.BinData(w0_0,f0_0,f0_1,bin_size)
    w_bin, y2_bin, y3_bin       = c.BinData(w0_0,f0_2,f0_3,bin_size)
    w_AG_bin, AG0_bin, AG1_bin  = c.BinData(w_AG_0,f_AG_0,f_AG_1,bin_size)
    w_AG_bin, AG2_bin, AG3_bin  = c.BinData(w_AG_0,f_AG_2,f_AG_3,bin_size)
    
    # To avoid dividing my 0 when calculating the ratios below, NaN and 0's are
    # replaced by the median. This is only for visual purposes.
    y0_bin = c.ReplaceWithMedian(y0_bin)
    y1_bin = c.ReplaceWithMedian(y1_bin)
    y2_bin = c.ReplaceWithMedian(y2_bin)
    y3_bin = c.ReplaceWithMedian(y3_bin)
    
    AG0_bin = c.ReplaceWithMedian(AG0_bin)
    AG1_bin = c.ReplaceWithMedian(AG1_bin)
    AG2_bin = c.ReplaceWithMedian(AG2_bin)
    AG3_bin = c.ReplaceWithMedian(AG3_bin)

    ratio1  = y0_bin/y1_bin
    ratio2  = y0_bin/y2_bin
    ratio3  = y0_bin/y3_bin

    # Creates a figure showing the quiet regions of the spectra.
    p.OverviewPlot(part, x_lim1, x_lim2, w0_0, w_bin, ratio1, ratio2, ratio3, y0_bin, y1_bin, y2_bin, y3_bin)

    # --------------------------------------------------------------------------
    # The heart of the code. This is where the the quiet region of the spectra
    # are used to cross correlate the spectra and calculate the shift.
    # --------------------------------------------------------------------------
    
    # These are dummy arrays used for the 10 Dec 2015 observations
    f_empty = []
    e_empty = []
    
    '''
    print "\n\nShifting 10 Dec 2015 observations:"
    W, F0_1, E0_1, F1_1, E1_1, F2_1, E2_1, AG1, AG1err, F_ave_w_1, E_ave_w_1 = c.ExportShitedSpectra(w0_0,f0_0,f0_1,f1_1,f2_1,f_empty,f_AG_1,e0_0,e0_1,e1_1,e2_1,e_empty,e_AG_1,NumFits_1,s1,s2,LyA)
    
    print "\n\nShifting 24 Dec 2015 observations:"
    W, F0_2, E0_2, F1_2, E1_2, F2_2, E2_2, F3_2, E3_2, AG2, AG2err, F_ave_w_2, E_ave_w_2 = c.ExportShitedSpectra(w0_0,f0_0,f0_2,f1_2,f2_2,f3_2,f_AG_2,e0_0,e0_2,e1_2,e2_2,e3_2,e_AG_2,NumFits_2,s1,s2,LyA)
    
    print "\n\nShifting 30 Jan 2016 observations:"
    W, F0_3, E0_3, F1_3, E1_3, F2_3, E2_3, F3_3, E3_3, AG3, AG3err, F_ave_w_3, E_ave_w_3 = c.ExportShitedSpectra(w0_0,f0_0,f0_3,F0_3,F1_3,F2_3,f_AG_3,e0_0,e0_3,e1_3,e2_3,e3_3,e_AG_3,NumFits_3,s1,s2,LyA)
    '''

    # Convert from wavlength to radial velocity
    #RV = c.wave2RV(W,LyA,RV_BP)

    # Save data into .dat file
    #np.savetxt(dat_directory+part+"2_2014.dat",np.column_stack((W, RV, f0_0, e0_0, f_AG_0, e_AG_0)))
    #np.savetxt(dat_directory+part+"2_10Dec.dat",np.column_stack((W, RV, F0_1, E0_1, F1_1, E1_1, F2_1, E2_1, AG1, AG1err, F_ave_w_1, E_ave_w_1)))
    #np.savetxt(dat_directory+part+"2_24Dec.dat",np.column_stack((W, RV, F0_2, E0_2, F1_2, E1_2, F2_2, E2_2, F3_2, E3_2, AG2, AG2err, F_ave_w_2, E_ave_w_2)))
    #np.savetxt(dat_directory+part+"2_30Jan.dat",np.column_stack((W, RV, F0_3, E0_3, F1_3, E1_3, F2_3, E2_3, F3_3, E3_3, AG3, AG3err, F_ave_w_3, E_ave_w_3)))

if __name__ == '__main__':
    main()
