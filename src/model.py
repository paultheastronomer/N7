import numpy as np
from scipy.special import wofz
import sys, time
import pandas as pd

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
        '''
        This code has been removed for speed purposes
        try:
             from scipy.special import wofz
        except ImportError:
             s = ("Can't find scipy.special.wofz(), can only calculate Voigt "
                  " function for 0 < a < 0.1 (a=%g)" % a)  
             print(s)
        else:
             return wofz(u + 1j * a).real
        '''
        return wofz(u + 1j * a).real

    def LSF(self,lsf_cen, W):
        ''' Tabulated Theoretical Line Spread Functions at Lifetime position 3
        Taken from http://www.stsci.edu/hst/cos/performance/spectral_resolution/
        '''
        # Load the wavelengths which have computed LSF values
        X =  np.genfromtxt('data/fuv_G130M_1291_lsf.dat', unpack=True).T[0]
        
        # Find the LSF closest to lsf_cen in Angstrom. See params.json.
        closest_LSF = min(X, key=lambda x:abs(x-lsf_cen))

        df      = pd.read_csv('data/fuv_G130M_1291_lsf.dat', delim_whitespace=True)
        Y       = df[0:][str(int(closest_LSF))]

        pix     = np.arange(-(len(Y))/2,len(Y)/2)+1

        kern    = np.arange(-len(W),len(W),1)

        LSF_kernel      = np.interp(kern,pix,Y)
        LSF_kernel      = LSF_kernel/np.sum(LSF_kernel)

        return LSF_kernel

    def K(self, W, l, sigma_kernel):
        ''' LSF
        Dispersion of the theoretical wavelength range.
        This LSF is based on a Gaussian with a FWHM given in the instrument parameter
        'sigma_kernel' found in params.json
        np.roll is equivalent to the IDL shift function
        '''
        # dl is the step size of the wavelength (l) in units of Angstrom
        # on which "kernel" is to be calculated.
        dl                  = np.mean((l-np.roll(l,1))[1:])
        dwave               = np.median((W-np.roll(W,1))[1:])   # Dispersion [Ang /pix]
        fwhm_cos_G130M_Ang  = sigma_kernel * dwave              # FWHM in Ang. 0.0648 Ang eq. to 6.5 pix.
        fwhm_cos_G130M_dl   = fwhm_cos_G130M_Ang / dl           # FWHM in Pix?
        c                   = fwhm_cos_G130M_dl/(2*np.sqrt(2*np.log(2.)))
        kern                = np.arange(-len(W),len(W),1)   # Same as l in fit.py
        kernel              = np.exp(-kern**2/(2*c**2))
        kernel              = kernel/np.sum(kernel)     
        return kernel

    def VoigtModel(self,params,Const):
        ''' Basic Voigt model to fit given data '''
        max_f, av, a, b         = params
        S,BetaPicRV,RV,kernel   = Const

        f       = max_f*self.voigt_wofz(av,RV)
        f       = f + (a*RV + b)
        f_star  = np.convolve(f,kernel,mode='same')
        
        return f_star

    def GaussianModel(self,params,Const):
        ''' Basic Gaussian model to fit given data '''
        A, sigma,  a, b     = params
        S,BetaPicRV,RV,kernel,mu        = Const

        u = RV + mu

        G = A*np.exp(-(u)**2/(2.*sigma**2))
        G = G + (a*u + b)

        f_star  =   np.convolve(G,kernel,mode='same')
        
        return f_star

    def GaussianModelDouble(self,params,Const):
        ''' Basic double Gaussian model to fit given data '''
        A, sigma1, sigma2, a, b = params
        S,BetaPicRV,RV,kernel,mu        = Const

        u = RV + mu

        G1 = np.exp(-(u)**2/(2.*sigma1**2))
        G2 = np.exp(-(u)**2/(2.*sigma2**2))

        Gtot = A*(G1+G2)+(a*u + b)

        f_star  =   np.convolve(Gtot,kernel,mode='same')

        return f_star

    def Continuum(self, param, window, l, W, F, E):  
        ''' Function used to model the continuum using only data outside
        of the lines being modeled. '''

        c1  = param["fit"]["windows"][window]["cut1"]
        c2  = param["fit"]["windows"][window]["cut2"]
        c3  = param["fit"]["windows"][window]["cut3"]
        c4  = param["fit"]["windows"][window]["cut4"]

        closest_W1  = min(W, key=lambda x:abs(x-c1))
        i1          = list(W).index(closest_W1)

        closest_W2  = min(W, key=lambda x:abs(x-c2))
        i2          = list(W).index(closest_W2)

        closest_W3  = min(W, key=lambda x:abs(x-c3))
        i3          = list(W).index(closest_W3)

        closest_W4  = min(W, key=lambda x:abs(x-c4))
        i4          = list(W).index(closest_W4)

        
        W = np.concatenate((W[:i1],W[i2:i3],W[i4:]))
        F = np.concatenate((F[:i1],F[i2:i3],F[i4:]))
        E = np.concatenate((E[:i1],E[i2:i3],E[i4:]))
 
        # Weights to apply to the y-coordinates of the sample points. For gaussian uncertainties, use 1/sigma (not 1/sigma**2).
        # https://docs.scipy.org/doc/numpy/reference/generated/numpy.polyfit.html
        weights     = 1./E
        
        z           = np.polyfit(W, F, param["fit"]["windows"][window]["order"], rcond=None, full=False, w=weights)
        pn          = np.poly1d(z)
        f           = pn(l)
        
        return f        

    def absorption(self, l,v_bp,nh,vturb,T,param,Nwindows):
        
        if Nwindows == 1:
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
        
        if Nwindows == 2:
            # [Hydrogen, Deuterium]   
            w       = [param["lines"]["line"]["N1"]["Wavelength"],
                      param["lines"]["line"]["N2"]["Wavelength"],
                      param["lines"]["line"]["N3"]["Wavelength"],
                      param["lines"]["line"]["Nw1"]["Wavelength"],
                      param["lines"]["line"]["Nw2"]["Wavelength"]]
            
            mass    = [param["lines"]["line"]["N1"]["Mass"],
                      param["lines"]["line"]["N2"]["Mass"],
                      param["lines"]["line"]["N3"]["Mass"],
                      param["lines"]["line"]["Nw1"]["Mass"],
                      param["lines"]["line"]["Nw2"]["Mass"]]
            
            fosc    = [param["lines"]["line"]["N1"]["Strength"],
                      param["lines"]["line"]["N2"]["Strength"],
                      param["lines"]["line"]["N3"]["Strength"],
                      param["lines"]["line"]["Nw1"]["Strength"],
                      param["lines"]["line"]["Nw2"]["Strength"]]
            
            delta   = np.array([param["lines"]["line"]["N1"]["Gamma"],
                      param["lines"]["line"]["N2"]["Gamma"],
                      param["lines"]["line"]["N3"]["Gamma"],
                      param["lines"]["line"]["Nw1"]["Gamma"],
                      param["lines"]["line"]["Nw2"]["Gamma"]]) /(4.*np.pi)                
            
            N_col   = np.array([1.,1.,1.,1.,1.])*10**nh
        
        c           = 2.99793e14        # Speed of light
        k           = 1.38064852e-23    # Boltzmann constant in J/K = m^2*kg/(s^2*K) in SI base units
        u           = 1.660539040e-27   # Atomic mass unit (Dalton) in kg
        absorption  = np.ones(len(l))

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
                    absorption[j]  =   absorption[j]*np.exp(-hav[j])       
                else:
                    absorption[j]  =   0.      
        return absorption

    def Absorptions(self,Const, params, param, sigma_kernel, Nwindows):
        
        if Nwindows == 1:
            W1,F1,E1,l1,BetaPicRV,v_ism,T_ism,v_bp,T_bp,T_X = Const
        if Nwindows == 2:
            W1,W2,F1,F2,E1,E2,l1,l2,BetaPicRV,v_ism,T_ism,v_bp,T_bp,T_X = Const

        nh_ism, b_ism, nh_bp, b_bp, nh_X, v_X, b_X     = params

        if param["fit"]["lsf"] == 'tabulated':
            kernel1      =   self.LSF(param["lines"]["line"]["N1"]["Wavelength"], W1)
        else:
            kernel1      =   self.K(W1, l1, sigma_kernel)

        
        # Calculates the ISM absorption
        abs_ism1     =   self.absorption(l1,v_ism,nh_ism,b_ism,T_ism,param,Nwindows)
        abs_bp1      =   self.absorption(l1,v_bp,nh_bp,b_bp,T_bp,param,Nwindows)
        abs_X1       =   self.absorption(l1,v_X,nh_X,b_X,T_X,param,Nwindows)
        
        
        # Continuum line
        f1           =   self.Continuum(param, param["display"]["window1"]["name"], l1, W1, F1, E1)
        
        # Stellar spectral profile, as seen from Earth
        # after absorption by the ISM and BP CS disk.
        # Profile has been convolved with HST LSF
        #    -  in (erg cm-2 s-1 A-1)
        
        f_abs_con1   =   np.convolve(f1*abs_ism1*abs_bp1*abs_X1, kernel1, mode='same')

        # Absorption by ISM
        f_abs_ism1   =   np.convolve(f1*abs_ism1, kernel1, mode='same')        

        # Absorption by beta Pictoris  
        f_abs_bp1    =   np.convolve(f1*abs_bp1, kernel1, mode='same')

        # Absorption by exocomets  
        f_abs_X1     =   np.convolve(f1*abs_X1, kernel1, mode='same')

        # Interpolation on COS wavelengths, relative to the star
        f_abs_int1   =   np.interp(W1,l1,f_abs_con1)

        unconvolved1             =  f1*abs_ism1*abs_bp1*abs_X1
                

        if Nwindows == 2:

            if param["fit"]["lsf"] == 'tabulated':
                kernel2      =   self.LSF(param["lines"]["line"]["Nw1"]["Wavelength"], W2)
            else:
                kernel2      =   self.K(W2, l2, sigma_kernel)

            kernel2      =   self.K(W2, l2, sigma_kernel)

            # Calculates the ISM absorption
            
            abs_ism2     =   self.absorption(l2,v_ism,nh_ism,b_ism,T_ism,param,Nwindows)
            abs_bp2      =   self.absorption(l2,v_bp,nh_bp,b_bp,T_bp,param,Nwindows)
            abs_X2       =   self.absorption(l2,v_X,nh_X,b_X,T_X,param,Nwindows)
            
            # Continuum line
            f2           =   self.Continuum(param, param["display"]["window2"]["name"], l2, W2, F2, E2)
            
            # Stellar spectral profile, as seen from Earth
            # after absorption by the ISM and BP CS disk.
            # Profile has been convolved with HST LSF
            #    -  in (erg cm-2 s-1 A-1)
            
            f_abs_con2   =   np.convolve(f2*abs_ism2*abs_bp2*abs_X2, kernel2, mode='same')

            # Absorption by ISM
            f_abs_ism2   =   np.convolve(f2*abs_ism2, kernel2, mode='same')        

            # Absorption by beta Pictoris  
            f_abs_bp2    =   np.convolve(f2*abs_bp2, kernel2, mode='same')

            # Absorption by exocomets  
            f_abs_X2     =    np.convolve(f2*abs_X2, kernel2, mode='same')

            # Interpolation on COS wavelengths, relative to the star
            f_abs_int2   =   np.interp(W2,l2,f_abs_con2)

            unconvolved2 =  f2*abs_ism2*abs_bp2*abs_X2
            
        if Nwindows == 1:
            return f_abs_int1, f_abs_ism1, f_abs_bp1, f_abs_X1, unconvolved1
        
        if Nwindows == 2:
            return f_abs_int1, f_abs_ism1, f_abs_bp1, f_abs_X1, unconvolved1, f_abs_int2, f_abs_ism2, f_abs_bp2, f_abs_X2, unconvolved2

    def Model(self, params, Const, ModelType, param):
        
        sigma_kernel    = param["instrument"]["sigma_kernel"]

        # Free parameters
        nh_ism, b_ism, nh_bp, b_bp, nh_X, v_X, b_X     = params

        Nwindows        = param["fit"]["windows"]["number"]

        if Nwindows == 1:
            # Fixed parameters
            W1,F1,E1,l1,BetaPicRV,v_ism,T_ism,v_bp,T_bp,T_X  = Const
            f_abs_int1, f_abs_ism1, f_abs_bp1, f_abs_X1, unconvolved1 = self.Absorptions(Const, params, param, sigma_kernel, Nwindows)      
            return f_abs_int1, f_abs_ism1, f_abs_bp1, f_abs_X1, unconvolved1
        if Nwindows == 2:
            # Fixed parameters
            W1,W2,F1,F2,E1,E2,l1,l2,BetaPicRV,v_ism,T_ism,v_bp,T_bp,T_X   = Const
            f_abs_int1, f_abs_ism1, f_abs_bp1, f_abs_X1, unconvolved1, f_abs_int2, f_abs_ism2, f_abs_bp2, f_abs_X2, unconvolved2  = self.Absorptions(Const, params, param,  sigma_kernel, Nwindows) 
            return f_abs_int1, f_abs_ism1, f_abs_bp1, f_abs_X1, unconvolved1, f_abs_int2, f_abs_ism2, f_abs_bp2, f_abs_X2, unconvolved2