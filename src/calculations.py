import numpy as np

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
