import numpy as np
from scipy.optimize import leastsq

from src.statistics import Stats
s   = Stats()

class Calc:
    '''
    A collection of calculations.
    '''

    def BinData(self,x,y1,e1,bin_pnts):
        bin_size    = int(len(x)/bin_pnts)
        bins        = np.linspace(x[0], x[-1], bin_size)
        digitized   = np.digitize(x, bins)
        bin_y       = np.array([y1[digitized == i].mean() for i in range(0, len(bins))])
        bin_e       = np.array([e1[digitized == i].mean() for i in range(0, len(bins))])
        return bins, bin_y, bin_e/np.sqrt(bin_pnts)

    def WeightedAvg(self,Flux, Err):
        """
        Return the weighted average and Error bars.
        """        
        Flux        = self.ReplaceWithMedian(Flux)
        Err         = self.ReplaceWithMedian(Err)
        weights     = 1./(Err**2)
        average     = np.average(Flux, axis=0, weights=weights)
        errorbars_2 = np.sum((weights*Err)**2, axis=0)
        return average, np.sqrt(errorbars_2)/ np.sum(weights, axis=0)

    def CF(self,flux,flux_err,ref,ref_err,n1,n2):
        flux        = self.ReplaceWithMedian(flux)
        flux_err    = self.ReplaceWithMedian(flux_err)
        ref_err     = self.ReplaceWithMedian(ref_err)
        ratio       = np.average(ref[n1:n2],  axis=0, weights=1./(ref_err[n1:n2]**2 ))/ \
                        np.average(flux[n1:n2], axis=0, weights=1./(flux_err[n1:n2]**2))                
        return ratio
        
    def Wave2RV(self,Wave,rest_wavelength,RV_BP):
        c = 299792458
        rest_wavelength = rest_wavelength*(RV_BP*1.e3)/c + rest_wavelength # Convert to beta pic reference frame
        delta_wavelength = Wave-rest_wavelength
        RV = ((delta_wavelength/rest_wavelength)*c)/1.e3	# km/s
        return RV

    def FindCenter(self,w,l):
        for i in range(len(w)):
          if w[i] > l:
            if abs(l - w[i-1]) < abs(l - w[i]):
              ans = i-1
              break
            else:
              ans = i
              break
        return ans

    def FindFactor(self, RV,D,E,A,l1,l2):
        region  = []
        err     = []
        E         = self.ReplaceWithMedian(E)
        for i in range(len(RV)):
            if l1 < RV[i] < l2:
                region.append(D[i])
                err.append(np.sqrt(E[i]**2+A[i]**2))
        region  = np.array(region)
        err     = np.array(err)
        
        factor, factor_err = self.WeightedAvg(region,err)

        print "Factor:",factor
        return factor, factor_err

    def ReplaceWithMedian(self, X):
        X[np.isnan(X)] = 0
        m = np.median(X[X > 0])
        X[X == 0] = m
        return X

    def ShiftAG(self, AG,units):
        zeros   = np.zeros(abs(units))
        if units > 0.:
            AG      = np.concatenate((zeros,AG))[:-units]
        else:
            AG      = np.concatenate((AG,zeros))[abs(units):]
        return AG

    def PrintParams(self, P, ConstB):
        print "\n",Fore.GREEN,"Free parameters",Style.RESET_ALL
        print Fore.RED,"Constant paramters",Style.RESET_ALL,"\n"
        
        print "ISM parameters:"
        print "-"*50,Fore.GREEN
        print "\tlog(N/1cm^2)\t=\t",P[0],"km/s"
        print "\t     b\t\t=\t",    P[1],Fore.RED
        print "\t     RV\t\t=\t",   ConstB[1],"km/s"
        print "\t     T\t\t=\t",    ConstB[2],"K",Style.RESET_ALL
        print "-"*50
        
        print "\nCS:"
        print "-"*50,Fore.GREEN
        print "\tlog(N/1cm^2)\t=\t",P[1]
        print "\t     b\t\t=\t",    P[2],"km/s",Fore.RED
        print "\t    RV\t\t=\t",    ConstB[3],"km/s"
        print "\t     T\t\t=\t",    ConstB[4],"K",Style.RESET_ALL
        print "-"*50,"\n"

        print "\nExocomet:"
        print "-"*50,Fore.GREEN
        print "\tlog(N/1cm^2)\t=\t",P[3],"km/s"
        print "\t    RV\t\t=\t",    P[4],"km/s"
        print "\t     b\t\t=\t",    P[5],"km/s",Fore.RED
        print "\t     T\t\t=\t",    ConstB[5],"K",Style.RESET_ALL
        print "-"*50,"\n\n\n"

    def FindBestParams(self, params,F,E,Const,ModelType, param):
        best_P, success = leastsq(s.chi2_lm, params, args=(F,E,Const,ModelType, param), maxfev=10000)
        return best_P

    def Window(self, param,W,F,E,WindowName):
        fit_start   = param["fit"]["windows"][WindowName]["start"]
        fit_end     = param["fit"]["windows"][WindowName]["stop"]

        s_i = []
        for i in range(len(W)):
            if fit_start <= W[i] <= fit_end:
                s_i.append(i)
        
        W    = W[s_i[0]:s_i[-1]]
        F    = F[s_i[0]:s_i[-1]]
        E    = E[s_i[0]:s_i[-1]]

        # Create an array of RV measurements with a resolution of 0.2 km/s

        v    = np.arange(-len(W)-100,len(W)+100,0.2) # RV values

        # Calculate the corresponding wavelengths
        l    = (W[0]+W[-1])/2.*(1.0 + v/3e5)

        return W, F, E, v, l