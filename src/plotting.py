import numpy as np
import matplotlib.pyplot as plt

from src.calculations import Calc
c   = Calc()

class Plotting:
    '''
    A collection of plotting routines.
    '''

    def BasicPlot(self,param,window,W,F,E,l,f_fit,f_abs_ism,f_abs_bp,f_abs_X,unconvolved):

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

        fig = plt.figure(figsize=(8,5))

        fontlabel_size  = 18
        tick_size       = 18
        params = {'backend': 'wxAgg', 'lines.markersize' : 2, 'axes.labelsize': fontlabel_size, 'font.size': fontlabel_size, 'legend.fontsize': 15, 'xtick.labelsize': tick_size, 'ytick.labelsize': tick_size, 'text.usetex': True}
        plt.rcParams.update(params)
        plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
        plt.rcParams['text.usetex'] = True
        plt.rcParams['text.latex.unicode'] = True 

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