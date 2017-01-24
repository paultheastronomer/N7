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

    def LSF(self, param, W):
        ''' Tabulated Theoretical Line Spread Functions at Lifetime position 3
        Taken from http://www.stsci.edu/hst/cos/performance/spectral_resolution/
        '''
        lsf_cen         = param["lines"]["line"]["N1"]["Wavelength"]
        # Load the wavelengths which have computed LSF values
        X               =  np.genfromtxt('data/fuv_G130M_1291_lsf.dat', unpack=True).T[0]
        
        # Find the LSF closest to lsf_cen in Angstrom. See params.json.
        closest_LSF     = min(X, key=lambda x:abs(x-lsf_cen))

        df              = pd.read_csv('data/fuv_G130M_1291_lsf.dat', delim_whitespace=True)
        Y               = df[0:][str(int(closest_LSF))]   

        LSF_pix_tabu    = np.arange(-(len(Y))/2,len(Y)/2)+0.5
        LSF_pix_kernel  = np.arange(-len(W)-100,len(W)+100,param["fit"]["pix_resolution"])       # In pixels. [-345:344]. Len tabulated LSF is 321

        dw_cos          = (W[-1]-W[0])/len(W)    # The size of a COS pixel in Angstrom: 1 pix = 0.00993320637681 Angstrom

        x_LSF           = LSF_pix_tabu*dw_cos   # The x-coordinates of the tabulated LSF in units of Angstrom. Length 321.
        y_LSF           = Y                     # The y-coordinates of the tabulated LSF ranging from 0 to 1. Length 321.
        x_val           = LSF_pix_kernel*dw_cos     # The x-coordinates of the interpolated values in RV
        
        LSF_kernel      = np.interp(x_val, x_LSF, y_LSF)
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
        kern                = np.arange(-len(W),len(W),1)       # Same as l in fit.py
        kernel              = np.exp(-kern**2/(2*c**2))
        kernel              = kernel/np.sum(kernel)     
        return kernel

    def Continuum(self, param, window, l, W, F, E):  
        ''' Function used to model the continuum using only data outside
        of the lines being modeled. '''

        # I am redefining W, F, E as I do wish to model a region different from the region
        # used to determine the continuum
        dat_directory   = param["directories"]["exoatmdir"]
        W, F, E         = np.genfromtxt(dat_directory+param["files"]["datafile"],unpack=True)

        c0  = param["fit"]["windows"][window]["cut0"]
        c1  = param["fit"]["windows"][window]["cut1"]
        c2  = param["fit"]["windows"][window]["cut2"]
        c3  = param["fit"]["windows"][window]["cut3"]
        c4  = param["fit"]["windows"][window]["cut4"]
        c5  = param["fit"]["windows"][window]["cut5"]

        closest_W0  = min(W, key=lambda x:abs(x-c0))
        i0          = list(W).index(closest_W0)

        closest_W1  = min(W, key=lambda x:abs(x-c1))
        i1          = list(W).index(closest_W1)

        closest_W2  = min(W, key=lambda x:abs(x-c2))
        i2          = list(W).index(closest_W2)

        closest_W3  = min(W, key=lambda x:abs(x-c3))
        i3          = list(W).index(closest_W3)

        closest_W4  = min(W, key=lambda x:abs(x-c4))
        i4          = list(W).index(closest_W4)

        closest_W5  = min(W, key=lambda x:abs(x-c5))
        i5          = list(W).index(closest_W5)

        W = np.concatenate((W[i0:i1],W[i2:i3],W[i4:i5]))
        F = np.concatenate((F[i0:i1],F[i2:i3],F[i4:i5]))
        E = np.concatenate((E[i0:i1],E[i2:i3],E[i4:i5]))
 
        # Weights to apply to the y-coordinates of the sample points. For gaussian uncertainties, use 1/sigma (not 1/sigma**2).
        # https://docs.scipy.org/doc/numpy/reference/generated/numpy.polyfit.html
        weights     = 1./E
        
        z           = np.polyfit(W, F, param["fit"]["windows"][window]["order"], rcond=None, full=False, w=weights)
        pn          = np.poly1d(z)
        f           = pn(l)
        
        return f        

    def absorption(self, l, v_comp, nN, nS, br, vturb, T, param, Nwindows):
        
        if Nwindows == 1:
            # [Hydrogen, Deuterium]   
            w       = [param["lines"]["line"]["N1"]["Wavelength"],
                      param["lines"]["line"]["N2"]["Wavelength"],
                      param["lines"]["line"]["N3"]["Wavelength"],
                      param["lines"]["line"]["SIII**"]["Wavelength"]]
            
            mass    = [param["lines"]["line"]["N1"]["Mass"],
                      param["lines"]["line"]["N2"]["Mass"],
                      param["lines"]["line"]["N3"]["Mass"],
                      param["lines"]["line"]["SIII**"]["Mass"]]
            
            fosc    = [param["lines"]["line"]["N1"]["Strength"],
                      param["lines"]["line"]["N2"]["Strength"],
                      param["lines"]["line"]["N3"]["Strength"],
                      param["lines"]["line"]["SIII**"]["Strength"]]
            
            delta   = np.array([param["lines"]["line"]["N1"]["Gamma"],
                      param["lines"]["line"]["N2"]["Gamma"],
                      param["lines"]["line"]["N3"]["Gamma"],
                      param["lines"]["line"]["SIII**"]["Gamma"]]) /(4.*np.pi)
            
            N_col   = np.array([1.,1.,1.])*10**nN
            N_S     = np.array([1.])*10**nS
        
        if Nwindows == 2:
            # [Hydrogen, Deuterium]   
            w       = [param["lines"]["line"]["N1"]["Wavelength"],
                      param["lines"]["line"]["N2"]["Wavelength"],
                      param["lines"]["line"]["N3"]["Wavelength"],
                      param["lines"]["line"]["Nw1"]["Wavelength"],
                      param["lines"]["line"]["Nw2"]["Wavelength"],
                      param["lines"]["line"]["SIII**"]["Wavelength"]]
            
            mass    = [param["lines"]["line"]["N1"]["Mass"],
                      param["lines"]["line"]["N2"]["Mass"],
                      param["lines"]["line"]["N3"]["Mass"],
                      param["lines"]["line"]["Nw1"]["Mass"],
                      param["lines"]["line"]["Nw2"]["Mass"],
                      param["lines"]["line"]["SIII**"]["Mass"]]
            
            fosc    = [param["lines"]["line"]["N1"]["Strength"],
                      param["lines"]["line"]["N2"]["Strength"],
                      param["lines"]["line"]["N3"]["Strength"],
                      param["lines"]["line"]["Nw1"]["Strength"],
                      param["lines"]["line"]["Nw2"]["Strength"],
                      param["lines"]["line"]["SIII**"]["Strength"]]
            
            delta   = np.array([param["lines"]["line"]["N1"]["Gamma"],
                      param["lines"]["line"]["N2"]["Gamma"],
                      param["lines"]["line"]["N3"]["Gamma"],
                      param["lines"]["line"]["Nw1"]["Gamma"],
                      param["lines"]["line"]["Nw2"]["Gamma"],
                      param["lines"]["line"]["SIII**"]["Gamma"]]) /(4.*np.pi)                
            
            N_col   = np.array([1.,1.,1.,1.,1.])*10**nN
            N_S     = np.array([1.])*10**nS

        if Nwindows == 3:
            # [Hydrogen, Deuterium]   
            w       = [param["lines"]["line"]["N1"]["Wavelength"],
                      param["lines"]["line"]["N2"]["Wavelength"],
                      param["lines"]["line"]["N3"]["Wavelength"],
                      param["lines"]["line"]["Nw1"]["Wavelength"],
                      param["lines"]["line"]["Nw2"]["Wavelength"],
                      param["lines"]["line"]["Nm1"]["Wavelength"],
                      param["lines"]["line"]["Nm2"]["Wavelength"],
                      param["lines"]["line"]["Nm3"]["Wavelength"],
                      param["lines"]["line"]["SIII**"]["Wavelength"]]
            
            mass    = [param["lines"]["line"]["N1"]["Mass"],
                      param["lines"]["line"]["N2"]["Mass"],
                      param["lines"]["line"]["N3"]["Mass"],
                      param["lines"]["line"]["Nw1"]["Mass"],
                      param["lines"]["line"]["Nw2"]["Mass"],
                      param["lines"]["line"]["Nm1"]["Mass"],
                      param["lines"]["line"]["Nm2"]["Mass"],
                      param["lines"]["line"]["Nm3"]["Mass"],
                      param["lines"]["line"]["SIII**"]["Mass"]]
            
            fosc    = [param["lines"]["line"]["N1"]["Strength"],
                      param["lines"]["line"]["N2"]["Strength"],
                      param["lines"]["line"]["N3"]["Strength"],
                      param["lines"]["line"]["Nw1"]["Strength"],
                      param["lines"]["line"]["Nw2"]["Strength"],
                      param["lines"]["line"]["Nm1"]["Strength"],
                      param["lines"]["line"]["Nm2"]["Strength"],
                      param["lines"]["line"]["Nm3"]["Strength"],
                      param["lines"]["line"]["SIII**"]["Strength"]]
            
            delta   = np.array([param["lines"]["line"]["N1"]["Gamma"],
                      param["lines"]["line"]["N2"]["Gamma"],
                      param["lines"]["line"]["N3"]["Gamma"],
                      param["lines"]["line"]["Nw1"]["Gamma"],
                      param["lines"]["line"]["Nw2"]["Gamma"],
                      param["lines"]["line"]["Nm1"]["Gamma"],
                      param["lines"]["line"]["Nm2"]["Gamma"],
                      param["lines"]["line"]["Nm3"]["Gamma"],
                      param["lines"]["line"]["SIII**"]["Gamma"]]) /(4.*np.pi)                
            
            N_col   = np.array([1.,1.,1.,1.,1.,1.,1.,1.])*10**nN
            N_S     = np.array([1.])*10**nS
        
        c_light     = 2.99793e14        # Speed of light
        k           = 1.38064852e-23    # Boltzmann constant in J/K = m^2*kg/(s^2*K) in SI base units
        u           = 1.660539040e-27   # Atomic mass unit (Dalton) in kg
        absorption  = np.ones(len(l))

        for i in range(len(w)):
            b_wid   = np.sqrt((T/mass[i]) + ((vturb/0.12895223)**2))
            b       = 4.30136955e-3*b_wid
            dnud    = b*c_light/w[i]
            xc      = l/(1.+v_comp*1.e9/c_light)    # In Angstrom
            v       = 1.e4*abs(((c_light/xc)-(c_light/w[i]))/dnud)
            if w[i] == param["lines"]["line"]["SIII**"]["Wavelength"]:
              tv      = 1.16117705e-14*N_S[0]*w[i]*fosc[i]/b_wid
            else:
              tv      = 1.16117705e-14*N_col[i]*w[i]*fosc[i]/b_wid
            a       = delta[i]/dnud
            hav     = tv*self.voigt_wofz(a,v)
            
            # To avoid calculating super tiny numbers
            for j in range(len(hav)):
                
                if hav[j] < 50:      
                    absorption[j]  =   absorption[j]*np.exp(-hav[j])       
                else:
                    absorption[j]  =   0.

        return absorption

    def Absorptions(self,Const, params, param, sigma_kernel, Nwindows):
        
        if Nwindows == 1:
            W1 ,F1, E1, l1, BetaPicRV, v_ISM, v_CS = Const
        if Nwindows == 2:
            W1, W2, F1, F2, E1, E2, l1, l2, BetaPicRV, v_ISM, v_CS = Const
        if Nwindows == 3:
            W1, W2, W3, F1, F2, F3, E1, E2, E3,l1,l2,l3, BetaPicRV, v_ISM, v_CS = Const

        nN_ISM, nS_ISM, b_ISM, T_ISM, xi_ISM, nN_CS, nS_CS, b_CS, T_CS, xi_CS, nN_X1, nS_X1, b_X1, T_X1, xi_X1, v_X1, nN_X2, nS_X2, b_X2, T_X2, xi_X2, v_X2     = params

        if param["fit"]["lsf"] == 'tabulated':
            kernel_w1  =   self.LSF(param, W1)
        else:
            kernel_w1  =   self.K(W1, l1, sigma_kernel)
        
        # Calculates the ISM absorption
        abs_ism_w1     =   self.absorption(l1, v_ISM, nN_ISM, nS_ISM, b_ISM, xi_ISM, T_ISM, param, Nwindows)
        abs_bp_w1      =   self.absorption(l1, v_CS, nN_CS, nS_CS, b_CS, xi_CS, T_CS, param, Nwindows)
        abs_X1_w1      =   self.absorption(l1, v_X1, nN_X1, nS_X1, b_X1, xi_X1, T_X1, param, Nwindows)
        abs_X2_w1      =   self.absorption(l1, v_X2, nN_X2, nS_X2, b_X2, xi_X2, T_X2, param, Nwindows)
        
        # Continuum line
        f_w1           =   self.Continuum(param, param["display"]["window1"]["name"], l1, W1, F1, E1)
        
        # Profile has been convolved with HST LSF
        #    -  in (erg cm-2 s-1 A-1)
        
        f_abs_con_w1    =   np.convolve(f_w1*abs_ism_w1*abs_bp_w1*abs_X1_w1*abs_X2_w1, kernel_w1, mode='same')

        # Absorption by ISM
        f_abs_ism_w1    =   np.convolve(f_w1*abs_ism_w1, kernel_w1, mode='same')        

        # Absorption by beta Pictoris  
        f_abs_bp_w1     =   np.convolve(f_w1*abs_bp_w1, kernel_w1, mode='same')

        # Absorption by exocomets  
        f_abs_X1_w1     =   np.convolve(f_w1*abs_X1_w1, kernel_w1, mode='same')
        f_abs_X2_w1     =   np.convolve(f_w1*abs_X2_w1, kernel_w1, mode='same')

        # Interpolation on COS wavelengths, relative to the star
        f_abs_int_w1    =   np.interp(W1,l1,f_abs_con_w1)

        unconvolved_w1  =  f_w1*abs_ism_w1*abs_bp_w1*abs_X1_w1*abs_X2_w1
                

        if Nwindows == 2:

            if param["fit"]["lsf"] == 'tabulated':
                kernel_w2      =   self.LSF(param, W2)
            else:
                kernel_w2      =   self.K(W2, l2, sigma_kernel)

            # Calculates the absorptions            
            abs_ism_w2     =   self.absorption(l2, v_ISM, nN_ISM, nS_ISM, b_ISM, xi_ISM, T_ISM, param, Nwindows)
            abs_bp_w2      =   self.absorption(l2, v_CS, nN_CS, nS_CS, b_CS, xi_CS, T_CS, param, Nwindows)
            abs_X1_w2      =   self.absorption(l2, v_X1, nN_X1, nS_X1, b_X1, xi_X1, T_X1, param, Nwindows)
            abs_X2_w2      =   self.absorption(l2, v_X2, nN_X2, nS_X2, b_X2, xi_X2, T_X2, param, Nwindows)
            
            # Continuum line
            f_w2           =   self.Continuum(param, param["display"]["window2"]["name"], l2, W2, F2, E2)
            
            f_abs_con_w2   =   np.convolve(f_w2*abs_ism_w2*abs_bp_w2*abs_X1_w2*abs_X2_w2, kernel_w2, mode='same')

            # Absorption by ISM
            f_abs_ism_w2   =   np.convolve(f_w2*abs_ism_w2, kernel_w2, mode='same')        

            # Absorption by beta Pictoris  
            f_abs_bp_w2    =   np.convolve(f_w2*abs_bp_w2, kernel_w2, mode='same')

            # Absorption by exocomets  
            f_abs_X1_w2    =   np.convolve(f_w2*abs_X1_w2, kernel_w2, mode='same')
            f_abs_X2_w2    =   np.convolve(f_w2*abs_X2_w2, kernel_w2, mode='same')

            # Interpolation on COS wavelengths, relative to the star
            f_abs_int_w2   =   np.interp(W2,l2,f_abs_con_w2)

            unconvolved_w2 =   f_w2*abs_ism_w2*abs_bp_w2*abs_X1_w2*abs_X2_w2

        if Nwindows == 3:

            if param["fit"]["lsf"] == 'tabulated':
                kernel_w2      =   self.LSF(param, W2)
            else:
                kernel_w2      =   self.K(W2, l2, sigma_kernel)

            # Calculates the absorptions            
            abs_ism_w2     =   self.absorption(l2, v_ISM, nN_ISM, nS_ISM, b_ISM, xi_ISM, T_ISM, param, Nwindows)
            abs_bp_w2      =   self.absorption(l2, v_CS, nN_CS, nS_CS, b_CS, xi_CS, T_CS, param, Nwindows)
            abs_X1_w2      =   self.absorption(l2, v_X1, nN_X1, nS_X1, b_X1, xi_X1, T_X1, param, Nwindows)
            abs_X2_w2      =   self.absorption(l2, v_X2, nN_X2, nS_X2, b_X2, xi_X2, T_X2, param, Nwindows)
            
            # Continuum line
            f_w2           =   self.Continuum(param, param["display"]["window2"]["name"], l2, W2, F2, E2)
            
            f_abs_con_w2   =   np.convolve(f_w2*abs_ism_w2*abs_bp_w2*abs_X1_w2*abs_X2_w2, kernel_w2, mode='same')

            # Absorption by ISM
            f_abs_ism_w2   =   np.convolve(f_w2*abs_ism_w2, kernel_w2, mode='same')        

            # Absorption by beta Pictoris  
            f_abs_bp_w2    =   np.convolve(f_w2*abs_bp_w2, kernel_w2, mode='same')

            # Absorption by exocomets  
            f_abs_X1_w2    =   np.convolve(f_w2*abs_X1_w2, kernel_w2, mode='same')
            f_abs_X2_w2    =   np.convolve(f_w2*abs_X2_w2, kernel_w2, mode='same')

            # Interpolation on COS wavelengths, relative to the star
            f_abs_int_w2   =   np.interp(W2,l2,f_abs_con_w2)

            unconvolved_w2 =   f_w2*abs_ism_w2*abs_bp_w2*abs_X1_w2*abs_X2_w2

            if param["fit"]["lsf"] == 'tabulated':
                kernel_w3  =   self.LSF(param, W2)
            else:
                kernel_w3  =   self.K(W3, l3, sigma_kernel)

            # Calculates the absorptions
            abs_ism_w3     =   self.absorption(l3, v_ISM, nN_ISM, nS_ISM, b_ISM, xi_ISM, T_ISM, param, Nwindows)
            abs_bp_w3      =   self.absorption(l3, v_CS, nN_CS, nS_CS, b_CS, xi_CS, T_CS, param, Nwindows)
            abs_X1_w3      =   self.absorption(l3, v_X1, nN_X1, nS_X1, b_X1, xi_X1, T_X1, param, Nwindows)
            abs_X2_w3      =   self.absorption(l3, v_X2, nN_X2, nS_X2, b_X2, xi_X2, T_X2, param, Nwindows)
            
            # Continuum line
            f_w3           =   self.Continuum(param, param["display"]["window3"]["name"], l3, W3, F3, E3)
            
            f_abs_con_w3   =   np.convolve(f_w3*abs_ism_w3*abs_bp_w3*abs_X1_w3*abs_X2_w3, kernel_w3, mode='same')

            # Absorption by ISM
            f_abs_ism_w3   =   np.convolve(f_w3*abs_ism_w3, kernel_w3, mode='same')        

            # Absorption by beta Pictoris  
            f_abs_bp_w3    =   np.convolve(f_w3*abs_bp_w3, kernel_w3, mode='same')

            # Absorption by exocomets
            f_abs_X1_w3     =    np.convolve(f_w3*abs_X1_w3, kernel_w3, mode='same')
            f_abs_X2_w3     =    np.convolve(f_w3*abs_X2_w3, kernel_w3, mode='same')  

            # Interpolation on COS wavelengths, relative to the star
            f_abs_int_w3   =   np.interp(W3,l3,f_abs_con_w3)

            unconvolved_w3 =   f_w3*abs_ism_w3*abs_bp_w3*abs_X1_w3*abs_X2_w3

            
        if Nwindows == 1:
            return f_abs_int_w1, f_abs_con_w1, f_abs_ism_w1, f_abs_bp_w1, f_abs_X1_w1, f_abs_X2_w1, unconvolved_w1
        
        if Nwindows == 2:
            return f_abs_int_w1, f_abs_con_w1, f_abs_ism_w1, f_abs_bp_w1, f_abs_X1_w1, f_abs_X2_w1, unconvolved_w1, f_abs_int_w2, f_abs_con_w2, f_abs_ism_w2, f_abs_bp_w2, f_abs_X1_w2, f_abs_X2_w2, unconvolved_w2

        if Nwindows == 3:
            return f_abs_int_w1, f_abs_con_w1, f_abs_ism_w1, f_abs_bp_w1, f_abs_X1_w1, f_abs_X2_w1, unconvolved_w1, f_abs_int_w2, f_abs_con_w2, f_abs_ism_w2, f_abs_bp_w2, f_abs_X1_w2, f_abs_X2_w2, unconvolved_w2, f_abs_int_w3, f_abs_con_w3, f_abs_ism_w3, f_abs_bp_w3,  f_abs_X1_w3, f_abs_X2_w3, unconvolved_w3

    def Model(self, params, Const, ModelType, param):
        
        sigma_kernel    = param["instrument"]["sigma_kernel"]

        # Free parameters
        nN_ISM, nS_ISM, b_ISM, T_ISM, xi_ISM, nN_CS, nS_CS, b_CS, T_CS, xi_CS, nN_X1, nS_X1, b_X1, T_X1, xi_X1, v_X1, nN_X2, nS_X2, b_X2, T_X2, xi_X2, v_X2    = params

        Nwindows        = param["fit"]["windows"]["number"]

        if Nwindows == 1:
            # Fixed parameters
            W1, F1, E1, l1, BetaPicRV, v_ISM, v_CS = Const
            return self.Absorptions(Const, params, param, sigma_kernel, Nwindows)

        if Nwindows == 2:
            # Fixed parameters
            W1, W2, F1, F2, E1, E2, l1, l2, BetaPicRV, v_ISM, v_CS   = Const
            return self.Absorptions(Const, params, param,  sigma_kernel, Nwindows)

        if Nwindows == 3:
            # Fixed parameters
            W1, W2, W3, F1, F2, F3, E1, E2, E3, l1, l2, l3, BetaPicRV, v_ISM, v_CS   = Const
            return self.Absorptions(Const, params, param,  sigma_kernel, Nwindows)