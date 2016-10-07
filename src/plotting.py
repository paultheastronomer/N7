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

    def OverviewPlot(self, part, x_lim1, x_lim2, w0_0, w_bin, ratio1, ratio2, ratio3, y0_bin, y1_bin, y2_bin, y3_bin, s1, s2):
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
        if part == 'A':
          plt.ylim(0.25,2.3)
        else:
          plt.ylim(0.7,4.0)
        plt.legend(loc='upper left', numpoints=1)


        #'''            
        ax2 = plt.subplot(212, sharex=ax1)
        # Plot of the individual spectra
        plt.step(w_bin,y0_bin,lw=1.2,color="#FF281C",label='2014')
        plt.step(w_bin,y1_bin,lw=1.2,color="#FF9303",label='2015v1')
        plt.step(w_bin,y2_bin,lw=1.2,color="#0386FF",label='2015v2')
        plt.step(w_bin,y3_bin,lw=1.2,color="#00B233",label='2016v3') 
        plt.plot([w0_0[s1],w0_0[s1]],[0,3],'-k')
        plt.plot([w0_0[s2],w0_0[s2]],[0,3],'-k')
        plt.xlim(x_lim1,x_lim2)
          
        if part == 'A':
          plt.ylim(0.,1.2e-12) 
        else:
          plt.ylim(0.,1.0e-13) 

        plt.legend(loc='upper left', numpoints=1)
        
        fig.tight_layout()
        plt.xlabel(r'Wavelength \AA')
        #plt.savefig('FEB_quiet_regions.pdf', bbox_inches='tight', pad_inches=0.1,dpi=300)
        plt.show()

    def CombinedPlot(self, param, window, w, W, F_tot, F_tot_err, F_tot_bin, FS1, FS1_bin, FS2, FS2_bin, F1, F2, F3, F4):
        
        fig = self.FigParams()

        x1  = param["display"][window]["x1"]
        x2  = param["display"][window]["x2"]
        y1  = param["display"][window]["y1"]
        y2  = param["display"][window]["y2"]

        #plt.step(W,F1/F_tot,lw=1.2,color="#FF281C",label='2014')
        #plt.step(W,F2,lw=1.2,color="#FF9303",label='2015v1')
        #plt.step(W,F3,lw=1.2,color="#0386FF",label='2015v2')
        #plt.step(W,F4,lw=1.2,color="#00B233",label='2016v3')
        plt.step(w+0.12,FS1*9.8e-18,lw=1.2,color="blue",label='Total')
        plt.step(w+0.12,FS2*9.8e-18,lw=1.2,color="red",label='Total')
        plt.step(W,F_tot_bin,lw=1.2,color="black",label='Total')
        #plt.xlim(x1,x2)
        plt.ylim(y1,y2*15)
        #plt.ylim(0.,2.e-15)
        plt.show()

    def BasicPlot(self,param,window,W,F,E,l,f_fit,f_abs_ism,f_abs_bp,f_abs_X,unconvolved):

        fig = self.FigParams()
        
        x1  = param["display"][window]["x1"]
        x2  = param["display"][window]["x2"]
        y1  = param["display"][window]["y1"]
        y2  = param["display"][window]["y2"]

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

        '''
        plt.figure(figsize=(8,5))

        fontlabel_size  = 18
        tick_size       = 18
        params = {'backend': 'wxAgg', 'lines.markersize' : 2, 'axes.labelsize': fontlabel_size, 'font.size': fontlabel_size, 'legend.fontsize': 15, 'xtick.labelsize': tick_size, 'ytick.labelsize': tick_size, 'text.usetex': True}
        plt.rcParams.update(params)
        plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
        plt.rcParams['text.usetex'] = True
        plt.rcParams['text.latex.unicode'] = True
        '''

        plt.plot(l,unconvolved,color="cyan")
        plt.plot(l,f_abs_ism,color="#0386FF",lw=2)
        plt.plot(l,f_abs_bp,color="#00B233",lw=2)   
        plt.plot(l,f_abs_X,color="#FF9303",lw=2)
        plt.plot(W,f_fit,lw=2,color='#FF281C',label=r'Best fit')

        #plt.step(W,F-f_fit,lw=2,color='#FF281C',label=r'Best fit')
        #plt.show()
        #sys.exit()

        if param["display"]["bin"] > 1:
            bin_size = param["display"]["bin"]
            Wb, Fb, Eb  = c.BinData(W,F,E,bin_size)
            plt.errorbar(Wb,np.ones(len(Wb))*2e-14,yerr=Eb)
            plt.step(Wb,Fb,color="#333333")
            plt.step(Wb[:i1/bin_size],Fb[:i1/bin_size],color="black",lw=2)
            plt.step(Wb[i2/bin_size:i3/bin_size],Fb[i2/bin_size:i3/bin_size],color="black",lw=2)
            plt.step(Wb[i4/bin_size:],Fb[i4/bin_size:],color="black",lw=2)
        else:
            line1 = param["lines"]["line"]["Nw1"]["Wavelength"]
            line2 = param["lines"]["line"]["Nw2"]["Wavelength"]
            plt.plot([line1,line1],[0.2e-14,0.3e-14],color='black')
            plt.plot([line2,line2],[0.2e-14,0.3e-14],color='black')
            plt.errorbar(W,np.ones(len(W))*2e-14,yerr=E)
            plt.step(W,F,color="#333333")
            plt.step(W[:i1],F[:i1],color="black",lw=2)
            plt.step(W[i2:i3],F[i2:i3],color="black",lw=2)
            plt.step(W[i4:],F[i4:],color="black",lw=2)


        plt.xlim(x1,x2)
        plt.ylim(y1,y2)
        '''
        if window == 'window1':
            x = [1199,1200,1201]
            labels = ['1199','1200','1201']
            plt.xticks(x, labels)
        if window == 'window2':
            x = [1160,1161]
            labels = ['1160','1161']
            plt.xticks(x, labels)
        '''
        plt.xlabel(r'Wavelength (\AA)')
        plt.ylabel(r'Flux (erg/s/cm$^2$/\AA)')

        plt.minorticks_on()
        fig.tight_layout()
        fig.savefig("plots/"+window+".pdf")
        plt.show()