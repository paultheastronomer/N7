import numpy as np
from scipy.special import wofz

class Model:
    '''
    A collection of functions for modeling the line absorption.
    '''
        
    def voigt_wofz(self, a, u):
        ''' Compute the Voigt function using Scipy's wofz().

        # Code from https://github.com/nhmc/Barak/blob/\
        087602fb372812c0603441ca1ce6820e96963b88/barak/absorb/voigt.py
        
        Explanation: https://nhmc.github.io/Barak/generated/barak.voigt.voigt.html#barak.voigt.voigt

        Parameters
        ----------
        a: float
          Ratio of Lorentzian to Gaussian linewidths.
        u: array of floats
          The frequency or velocity offsets from the line centre, in units
          of the Gaussian broadening linewidth.

        See the notes for `voigt` for more details.
        '''
        try:
             from scipy.special import wofz
        except ImportError:
             s = ("Can't find scipy.special.wofz(), can only calculate Voigt "
                  " function for 0 < a < 0.1 (a=%g)" % a)  
             print(s)
        else:
             return wofz(u + 1j * a).real


    def K(self, W, l, sigma_kernel):
        ''' LSF
        Dispersion of the theoretical wavelength range
        np.roll is equivalent to the IDL shift function
        '''
        dispersion          =   9.97e-3                 # Ang /pix
        dl                  = np.mean((l-np.roll(l,1))[1:])
        # dl is the step size of the wavelength on which "kernel" is to be calculated.
        dwave               = np.median((W-np.roll(W,1))[1:])
        fwhm_cos_G130M_Ang  = sigma_kernel * dispersion # FWHM in Ang.
        fwhm_cos_G130M_dl   = fwhm_cos_G130M_Ang / dl
        sigma_kernel_dl     = fwhm_cos_G130M_dl / (2.*np.sqrt(2.*np.log(2.)))
        kernel              = np.arange(-len(W)/2.,len(W)/2.,1)
        kernel              = np.exp(-kernel**2/2./((sigma_kernel*dwave/dl)**2))
        kernel              = kernel/np.sum(kernel)     
        
        return kernel


    def Continuum(self, param, l, W, F, E):  
        
        c1  = param["fit"]["continuum"]["cut1"]
        c2  = param["fit"]["continuum"]["cut2"]
        
        W = np.concatenate((W[:c1],W[-c2:]))
        F = np.concatenate((F[:c1],F[-c2:]))
        E = np.concatenate((E[:c1],E[-c2:]))
        
        weights     = 1./E#**2 # Weights to apply to the y-coordinates of the sample points. For gaussian uncertainties, use 1/sigma (not 1/sigma**2). https://docs.scipy.org/doc/numpy/reference/generated/numpy.polyfit.html
        z           = np.polyfit(W, F, param["fit"]["continuum"]["order"], rcond=None, full=False)#, w=weights)
        pn          = np.poly1d(z)
        f           = pn(l)
        return f        

    def absorption(self, l,v_bp,nh,vturb,T,param):
        
        # [Hydrogen, Deuterium]   
        w       = [param["lines"]["line"]["N1"]["Wavelength"],
                  param["lines"]["line"]["N2"]["Wavelength"],
                  param["lines"]["line"]["N3"]["Wavelength"]]
        
        mass    = [param["lines"]["line"]["N1"]["Mass"],
                  param["lines"]["line"]["N2"]["Mass"],
                  param["lines"]["line"]["N3"]["Mass"]]
        
        fosc    = [param["lines"]["line"]["N1"]["Strength"],
                  param["lines"]["line"]["N2"]["Strength"],
                  param["lines"]["line"]["N3"]["Strength"]]
        
        delta   = np.array([param["lines"]["line"]["N1"]["Gamma"],
                  param["lines"]["line"]["N2"]["Gamma"],
                  param["lines"]["line"]["N3"]["Gamma"]]) /(4.*np.pi)
        
        N_col   = np.array([1.,1.,1.])*10**nh
        
        c       = 2.99793e14
        k       = 1.38064852e-23    # Boltzmann constant in J/K = m^2*kg/(s^2*K) in SI base units
        u       = 1.660539040e-27   # Atomic mass unit (Dalton) in kg

        abs_ism = np.ones(len(l))

        for i in range(len(w)):
            b_wid   = np.sqrt((T/mass[i]) + ((vturb/np.sqrt(2*k/u)/1e3)**2)) # non-thermal + thermal broadening
            b       = 4.30136955e-3*b_wid
            dnud    = b*c/w[i]
            xc      = l/(1.+v_bp*1.e9/c)
            v       = 1.e4*abs(((c/xc)-(c/w[i]))/dnud)
            tv      = 1.16117705e-14*N_col[i]*w[i]*fosc[i]/b_wid
            a       = delta[i]/dnud
            hav     = tv*self.voigt_wofz(a,v)
            
            # To avoid underflow which occurs when you have exp(small negative number)
            for j in range(len(hav)):
                if hav[j] < 20.:      
                    abs_ism[j]  =   abs_ism[j]*np.exp(-hav[j])       
                else:
                    abs_ism[j]  =   0.
                    
        return abs_ism


    def LyModel(self, params, Const, ModelType, param):
        
        '''
        ModelType refers to the kind of model you are interested in.
  
        ModelType = 1
        ========================================================================
        No extra components but the H absorption is not fixed to the beta pic
        reference frame, but is free to vary.
        ========================================================================                 
        '''
        sigma_kernel    = param["instrument"]["sigma_kernel"]
        
        if ModelType == 1:
            # Free parameters
            nh_bp     = params
 
            # Fixed parameters
            W,F,E,l,L1,L2,L3,BetaPicRV,nh_ism,v_ism,b_ism,T_ism,v_bp,b_bp,T_bp   = Const
        '''
        l = []
        #W = []
        for i in range(len(lw)):
            if s1 <= lw[i] <= s2:
                l.append(lw[i])
        #for i in range(len(w)):
        #    if s1 <= w[i] <= s2:
        #        W.append(w[i])
        
        l = np.array(l)
        W = w
        #W = np.array(W)
        '''
        
        kernel      =   self.K(W, l, sigma_kernel)

        # Calculates the ISM absorption
        abs_ism     =   self.absorption(l,v_ism,nh_ism,b_ism,T_ism,param)
        abs_bp      =   self.absorption(l,v_bp,nh_bp,b_bp,T_bp,param)


        # Continuum line
        f           =   self.Continuum(param, l, W, F, E)
       
        # Stellar spectral profile, as seen from Earth
        # after absorption by the ISM and BP CS disk.
        # Profile has been convolved with HST LSF
        #    -  in (erg cm-2 s-1 A-1)
        
        f_abs_con   =   np.convolve(f*abs_ism*abs_bp, kernel, mode='same')
        
        f_abs_ism   =   np.convolve(f*abs_ism, kernel, mode='same')
        
        # Absorption by beta Pictoris  
        f_abs_bp    =   np.convolve(f*abs_bp, kernel, mode='same')
        
        # Interpolation on COS wavelengths, relative to the star
        f_abs_int   =   np.interp(W,l,f_abs_con)
                    
        return f_abs_int, f_abs_ism, f_abs_bp
