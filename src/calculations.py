import numpy as np
from scipy.optimize import leastsq
import pyfits, os

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

    def ExportShitedSpectra(self, w0,f0,f1,f2,f3,f4,AG,e0,e1,e2,e3,e4,eAG,NumFits,start,stop,rest_wavelength):
       
        # Creating empty arrays to be filled with
        # shifted spectra
        F = [[] for _ in range(NumFits-1)]    # -1 because we want to avoid the airglow observation
        W = [[] for _ in range(NumFits-1)]
        E = [[] for _ in range(NumFits-1)]    

        W[0],F[0],E[0]  =   self.ShiftSpec(f0,f1,e1,w0,start,stop,rest_wavelength)
        W[1],F[1],E[1]  =   self.ShiftSpec(f0,f2,e2,w0,start,stop,rest_wavelength)
        W[2],F[2],E[2]  =   self.ShiftSpec(f0,f3,e3,w0,start,stop,rest_wavelength)
        
        #W[3],F[3],E[3] =   shift_spec(f0,f3,e3,w0,start,stop)

        if NumFits > 4:
            W[3],F[3],E[3]  =   self.ShiftSpec(f0,f4,e4,w0,start,stop,rest_wavelength)

        F = np.array(F)
        E = np.array(E)
        
        #F_ave_w =  np.average(F, axis=0,weights=1./E**2)
        
        F_ave_w, E_ave_w    = self.WeightedAvg(F, E)
        

        if NumFits > 4:
            return W[0], F[0], E[0], F[1], E[1], F[2], E[2], F[3], E[3], AG, eAG, F_ave_w, E_ave_w
        elif NumFits < 4:
            return W[0], F[0], E[0], AG, eAG, F_ave_w, E_ave_w
        else:
            return W[0], F[0], E[0], F[1], E[1], F[2], E[2], AG, eAG, F_ave_w, E_ave_w

    def ExtractData(self, fits_file, part):
        f           = pyfits.open(fits_file)
        tbdata      = f[1].data
        net         = tbdata['NET']
        gcounts     = tbdata['GCOUNTS']
        exptime     = np.array([tbdata['EXPTIME']])
        wavelength  = tbdata['WAVELENGTH']
        flux        = tbdata['FLUX'] 
        a           = net*exptime.T
        for i in range(len(a)):
            a[i]        = [1e-15 if x==0 else x for x in a[i]]
        err         = np.sqrt(gcounts+1)*(flux / (a))
        if part == 'A':
            return wavelength[0], flux[0], err[0]
        else:
            return wavelength[1], flux[1], err[1]   

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

    def GetData(self, fits_location,part):
        dir_contents    = os.listdir(fits_location)
        fits            = sorted([fn for fn in dir_contents if fn.startswith('l') and fn.endswith('sum.fits')])
        NumFits         = len(fits)

        # Extracting data from fits files.
        #
        
        if fits_location[-13:] == '2014/visit_1/':
            # The 2014 data had no shifts.
            wavelength0, flux0, err0    = self.ExtractData(fits_location+fits[0],part)
            wavelength_AG, flux_AG, err_AG  = self.ExtractData(fits_location+fits[1],part)
            return wavelength0, flux0 ,err0, wavelength_AG, flux_AG, err_AG, NumFits
        
        elif fits_location[-13:] == '2015/visit_1/':
            # 2015 visit 1 data had a shift +0.8.
            wavelength0, flux0, err0    = self.ExtractData(fits_location+fits[0],part)
            wavelength1, flux1, err1        = self.ExtractData(fits_location+fits[1],part)
            wavelength2, flux2, err2        = self.ExtractData(fits_location+fits[2],part)
            wavelength_AG, flux_AG, err_AG  = self.ExtractData(fits_location+fits[3],part)
            return wavelength0, wavelength1, wavelength2, flux0, flux1, flux2, err0, err1, err2, wavelength_AG, flux_AG, err_AG, NumFits
        
        else:
            # The 2015 v2 and 2016 data had multiple shifts.
            wavelength0, flux0, err0    = self.ExtractData(fits_location+fits[0],part)
            wavelength1, flux1, err1        = self.ExtractData(fits_location+fits[1],part)
            wavelength2, flux2, err2        = self.ExtractData(fits_location+fits[2],part)
            wavelength3, flux3, err3        = self.ExtractData(fits_location+fits[3],part)
            wavelength_AG, flux_AG, err_AG  = self.ExtractData(fits_location+fits[4],part)
            return wavelength0, wavelength1, wavelength2, wavelength3, flux0, flux1, flux2, flux3, err0, err1, err2, err3, wavelength_AG, flux_AG, err_AG, NumFits

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

    def ShiftSpec(self, ref,spec,error,wave,start,stop,rest_wavelength):
        # This routine correlates the spectrum: spec
        # with the reference spectrum: ref
        ref_c   = ref[start:stop]   # Reference spectrum
        spec_c  = spec[start:stop]  # Spectrum to be shifted
        error_c = error[start:stop] # Error on spec to be shifted
        wave_c  = wave[start:stop]  # Wavelength of spec to be shifted
        
        ref_c       = ref_c - np.mean(ref_c)
        spec_c      = spec_c - np.mean(spec_c)
        error_c       = error_c - np.mean(error_c)

        c           = np.correlate(spec_c,ref_c,mode='full')

        x           = np.arange(c.size)
        c_max       = np.argmax(c)      # Maximum correlation
        
        c_light = 299792458

        print "=================================="
        if ref_c.size-1 > c_max:        # If spectrum is redshifted
          #ref_c.size-1 because python starts index at 0
          shift = wave_c[ref_c.size-1]-wave_c[np.argmax(c)]
          units = (ref_c.size-1)-np.argmax(c)
          RV_shift = ((shift/rest_wavelength)*c_light)/1.e3
          print "Pixel Shift:\t",units
          print "Angstrom Shift:\t",shift
          print "RV Shift:\t",RV_shift
          zeros = np.zeros(units)
          if units != 0:
            spec    = np.concatenate((zeros, spec), axis=0)[:-units]
            error   = np.concatenate((zeros, error), axis=0)[:-units]
        else:                           # If spectrum is blueshifted
          c = np.correlate(ref_c,spec_c,mode='full')
          shift = wave_c[np.argmax(c)]-wave_c[ref_c.size-1]
          units = abs(np.argmax(c)-(ref_c.size-1))
          RV_shift = ((shift/rest_wavelength)*c_light)/1.e3
          print "Pixel Shift:\t",-units
          print "Angstrom Shift:\t",shift
          print "RV Shift:\t",RV_shift
          zeros     = np.zeros(units)
          if units != 0:
            spec    = np.concatenate((spec, zeros), axis=0)[units:]
            error   = np.concatenate((error, zeros), axis=0)[units:]
        print "=================================="

        return wave,spec,error

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