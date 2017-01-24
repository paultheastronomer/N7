import numpy as np
import matplotlib.pyplot as plt

from src.calculations import Calc
c   = Calc()

class Plotting:
    '''
    A collection of plotting routines.
    '''
    def FigParams(self):
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

        return fig

    def AGPlot(self, w_AG_0,f_AG_0,w_AG_1,f_AG_1,w_AG_2,f_AG_2,w_AG_3,f_AG_3):
        fig = self.FigParams()
        #fig = plt.figure(figsize=(10,14))


        # Top plot shows the ratio between the spectra. Flat regions are
        # indicative of low FEB activity.
        #ax1 = plt.subplot(211)
        # Plot the ratios of the specra
        plt.step(w_AG_0,f_AG_0)
        plt.step(w_AG_1,f_AG_1)
        plt.step(w_AG_2,f_AG_2)
        plt.step(w_AG_3,f_AG_3)

        plt.show()


    def OverviewPlot(self, part, x_lim1, x_lim2, y_lim1, y_lim2, w0_0, w_bin, ratio1, ratio2, ratio3, y0_bin, y1_bin, y2_bin, y3_bin, s1, s2):
        fig = self.FigParams()
        #fig = plt.figure(figsize=(10,14))


        # Top plot shows the ratio between the spectra. Flat regions are
        # indicative of low FEB activity.
        ax1 = plt.subplot(211)
        # Plot the ratios of the specra
        plt.step(w_bin,ratio1,label='2014/2015v1',color="#FF32B1") 
        plt.step(w_bin,ratio2+1,label='2014/2015v2',color="#13C1CC")
        plt.step(w_bin,ratio3+2,label='2014/2016v3',color="#BDAC42")
        plt.plot([w0_0[s1],w0_0[s1]],[0,10],'-k')
        plt.plot([w0_0[s2],w0_0[s2]],[0,10],'-k')
        plt.xlim(x_lim1,x_lim2)
        plt.ylim(0,4)
        plt.legend(loc='upper left', numpoints=1)


        #'''            
        ax2 = plt.subplot(212, sharex=ax1)
        # Plot of the individual spectra
        plt.step(w_bin,y2_bin,lw=1.2,color="#0386FF",label='2015v2')
        plt.step(w_bin,y3_bin,lw=1.2,color="#00B233",label='2016v3')
        plt.step(w_bin,y1_bin,lw=1.2,color="#FF9303",label='2015v1')
        plt.step(w_bin,y0_bin,lw=1.2,color="#FF281C",label='2014')
        plt.plot([w0_0[s1],w0_0[s1]],[0,3],'-k')
        plt.plot([w0_0[s2],w0_0[s2]],[0,3],'-k')
        plt.xlim(x_lim1,x_lim2)
        plt.ylim(y_lim1,y_lim2)          

        plt.legend(loc='upper left', numpoints=1)
        
        fig.tight_layout()
        plt.xlabel(r'Wavelength \AA')
        #plt.savefig('FEB_quiet_regions.pdf', bbox_inches='tight', pad_inches=0.1,dpi=300)
        plt.show()

    def CombinedPlot(self, param, window, W, F_tot, W_bin, F2_bin, F3_bin, F4_bin):
        
        fig = self.FigParams()

        x1  = param["display"][window]["x1"]
        x2  = param["display"][window]["x2"]
        y1  = param["display"][window]["y1"]
        y2  = param["display"][window]["y2"]

        plt.step(W_bin, F2_bin, color='#FF9303', lw=1.3)
        plt.step(W_bin, F3_bin, color='#0386FF', lw=1.3)
        plt.step(W_bin, F4_bin, color='#00B233', lw=1.3)
        plt.step(W, F_tot, color='black', lw=1.3)
        plt.xlim(x1,x2)
        plt.ylim(y1,y2)
        plt.xlabel(r'Wavelength (\AA)')
        plt.ylabel(r'Flux (erg/s/cm$^2$/\AA)')
        plt.minorticks_on()
        fig.tight_layout()
        plt.show()

    def BasicPlot(self,param,window,W, F, E,l,f_fit,f_abs_ism,f_abs_bp,f_abs_X1,f_abs_X2,unconvolved):

        fig = self.FigParams()
        
        x1  = param["display"][window]["x1"]
        x2  = param["display"][window]["x2"]
        y1  = param["display"][window]["y1"]
        y2  = param["display"][window]["y2"]

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

        plt.plot(l,unconvolved,color="cyan")
        plt.plot(l,f_abs_X1,color="#FF9303",lw=2)
        plt.plot(l,f_abs_X2,color="#FF9303",lw=2)
        plt.plot(l,f_abs_ism,color="#0386FF",lw=2)
        plt.plot(l,f_abs_bp,color="#00B233",lw=2)   
        plt.plot(l,f_fit,lw=2,color='#FF281C',label=r'Best fit')

        if param["display"]["bin"] > 1:
            bin_size = param["display"]["bin"]
            Wb, Fb, Eb      = c.BinData(W,F,E,bin_size)
            #plt.errorbar(Wb,np.ones(len(Wb))*2e-14,yerr=Eb)
            plt.step(Wb,Fb,color="black",lw=1)
            plt.step(Wb[i0/bin_size:i1/bin_size],Fb[i0/bin_size:i1/bin_size],color="black",lw=2)
            plt.step(Wb[i2/bin_size:i3/bin_size],Fb[i2/bin_size:i3/bin_size],color="black",lw=2)
            plt.step(Wb[i4/bin_size:i5/bin_size],Fb[i4/bin_size:i5/bin_size],color="black",lw=2)
        else:
            line1 = param["lines"]["line"]["Nw1"]["Wavelength"]
            line2 = param["lines"]["line"]["Nw2"]["Wavelength"]
            plt.plot([line1,line1],[0.2e-14,0.3e-14],color='black')
            plt.plot([line2,line2],[0.2e-14,0.3e-14],color='black')
            #plt.errorbar(W,np.ones(len(W))*0.2e-14,yerr=E)
            plt.step(W,F,color="black",lw=1)
            plt.step(W[i0:i1],F[i0:i1],color="black",lw=2)
            plt.step(W[i2:i3],F[i2:i3],color="black",lw=2)
            plt.step(W[i4:i5],F[i4:i5],color="black",lw=2)


        plt.xlim(x1,x2)
        plt.ylim(y1,y2)
        #'''
        if window == 'window1':
            x = [1199,1199.5,1200,1200.5,1201,1201.5,1202]
            labels = ['1199.0','1199.5','1200.0','1200.5','1201','1201.5','1202.0']
            plt.xticks(x, labels)
        if window == 'window2':
            x = [1160.0,1160.5,1161.0,1161.5]
            labels = ['1160.0','1160.5','1161.0','1161.5']
            plt.xticks(x, labels)
        if window == 'window3':
            x = [1133.5,1134.0,1134.5,1135.0,1135.5]
            labels = ['1133.5','1134.0','1134.5','1135.0','1135.5']
            plt.xticks(x, labels)
        #'''
        plt.xlabel(r'Wavelength (\AA)')
        plt.ylabel(r'Flux (erg/s/cm$^2$/\AA)')

        plt.minorticks_on()
        fig.tight_layout()
        #fig.savefig("plots/"+window+".png")
        plt.show()

    def OwensPlot(self, param, window, W, NoPSF, F, E, Continuum, Fit, NI_1, SIII_1, NI_2, SIII_2, NI_3, SIII_3):

        fig = self.FigParams()
        
        x1  = param["display"][window]["x1"]
        x2  = param["display"][window]["x2"]
        y1  = param["display"][window]["y1"]
        y2  = param["display"][window]["y2"]


        if param["display"]["bin"] > 1:
            bin_size = param["display"]["bin"]
            Wb, Fb, Eb      = c.BinData(W,F,E,bin_size)
            #plt.errorbar(Wb,np.ones(len(Wb))*2e-14,yerr=Eb)
            plt.step(Wb,Fb,color="black",lw=1.2)
        else:
            line1 = param["lines"]["line"]["Nw1"]["Wavelength"]
            line2 = param["lines"]["line"]["Nw2"]["Wavelength"]
            plt.plot([line1,line1],[0.2e-14,0.3e-14],color='black')
            plt.plot([line2,line2],[0.2e-14,0.3e-14],color='black')
            #plt.errorbar(W,np.ones(len(W))*0.2e-14,yerr=E)
            plt.step(W,F,color="black",lw=1.2)

        plt.plot(W,NoPSF,color="cyan")
        plt.plot(W,SIII_3,color="#FF9303",lw=2)
        plt.plot(W,SIII_1,color="#0386FF",lw=2)
        plt.plot(W,SIII_2,color="#00B233",lw=2)   
        plt.plot(W,NI_3,color="#FF9303",lw=2)
        plt.plot(W,NI_1,color="#0386FF",lw=2)
        plt.plot(W,NI_2,color="#00B233",lw=2)   
        plt.plot(W,Fit,lw=2,color='#FF281C',label=r'Best fit')

        plt.xlim(x1,x2)
        plt.ylim(y1,y2)

        if window == 'window1':
            x = [1199,1199.5,1200,1200.5,1201,1201.5,1202]
            labels = ['1199.0','1199.5','1200.0','1200.5','1201','1201.5','1202.0']
            plt.xticks(x, labels)
        if window == 'window2':
            x = [1160.0,1160.5,1161.0,1161.5]
            labels = ['1160.0','1160.5','1161.0','1161.5']
            plt.xticks(x, labels)
        if window == 'window3':
            x = [1133.5,1134.0,1134.5,1135.0,1135.5]
            labels = ['1133.5','1134.0','1134.5','1135.0','1135.5']
            plt.xticks(x, labels)

        plt.xlabel(r'Wavelength (\AA)')
        plt.ylabel(r'Flux (erg/s/cm$^2$/\AA)')

        plt.minorticks_on()
        fig.tight_layout()
        fig.savefig("plots/"+window+"_NI_diff_align.pdf")
        plt.show()

    def CompareSpec(self, param, window, W, F1, F2):

        fig = self.FigParams()
        
        x1  = param["display"][window]["x1"]
        x2  = param["display"][window]["x2"]
        y1  = param["display"][window]["y1"]
        y2  = param["display"][window]["y2"]


        if param["display"]["bin"] > 1:
            bin_size = param["display"]["bin"]
            Wb, F1b, Db      = c.BinData(W,F1,F1,bin_size)
            Wb, F2b, Db      = c.BinData(W,F2,F2,bin_size)
            #plt.errorbar(Wb,np.ones(len(Wb))*2e-14,yerr=Eb)
            plt.step(Wb,F1b,color="black",lw=1.2)
            plt.step(Wb,F2b,color="red",lw=1.2)
        else:
            line1 = param["lines"]["line"]["Nw1"]["Wavelength"]
            line2 = param["lines"]["line"]["Nw2"]["Wavelength"]
            plt.plot([line1,line1],[0.2e-14,0.3e-14],color='black')
            plt.plot([line2,line2],[0.2e-14,0.3e-14],color='black')
            #plt.errorbar(W,np.ones(len(W))*0.2e-14,yerr=E)
            plt.step(W,F1,color="black",lw=1.2)
            plt.step(W,F2,color="red",lw=1.2)

        plt.xlim(x1,x2)
        plt.ylim(y1,y2)

        if window == 'window1':
            x = [1199,1199.5,1200,1200.5,1201,1201.5,1202]
            labels = ['1199.0','1199.5','1200.0','1200.5','1201','1201.5','1202.0']
            plt.xticks(x, labels)
        if window == 'window2':
            x = [1160.0,1160.5,1161.0,1161.5]
            labels = ['1160.0','1160.5','1161.0','1161.5']
            plt.xticks(x, labels)
        if window == 'window3':
            x = [1133.5,1134.0,1134.5,1135.0,1135.5]
            labels = ['1133.5','1134.0','1134.5','1135.0','1135.5']
            plt.xticks(x, labels)

        plt.xlabel(r'Wavelength (\AA)')
        plt.ylabel(r'Flux (erg/s/cm$^2$/\AA)')

        plt.minorticks_on()
        fig.tight_layout()
        fig.savefig("plots/"+window+"_comparison.pdf")
        plt.show()