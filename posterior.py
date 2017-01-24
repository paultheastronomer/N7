import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import corner
import pandas as pd
import os, sys, json

def Distrib(x):
   '''Finds median and 68% interval of array x.'''
   y    = sorted(x)
   up   = y[int(0.8413*len(y))]
   down = y[int(0.1587*len(y))]
   med  = y[int(0.5*len(y))]
   
   return med,up,down   

def Uncertainties(x):
   '''Finds median and 68% interval of array x.'''
   y    = sorted(x)
   up   = y[int(0.8413*len(y))]
   down = y[int(0.1587*len(y))]
   med  = y[int(0.5*len(y))]
   
   return med,up-med,med-down   

def Initialise():
    with open('params.json') as param_file:    
        param = json.load(param_file)
    return param

param           = Initialise()
# Select which posterior distributions to use.
letter = param["fit"]["MCMC"]["chain_name"]#'3EXO_Tab_1'
chains = [f for f in os.listdir('chains') if f.startswith(letter)]
stats  = [f for f in os.listdir('chains') if f.startswith("Stat_"+letter)]
number_of_chains = len(chains)

chain = []
stat  = []
for i in range(1,number_of_chains):
  chain.append(np.load('chains/'+chains[i]))
  stat.append(np.load('chains/'+stats[i]))

Chi2    = stat[0]['Chi2']

print "\nNumber of accepted steps: ",round((len(Chi2)*number_of_chains)/1.0e6,2)," million"
print "Best Chi2 value obtained: ",round(Chi2.min(),2),"\n"



nN_ISM  = chain[0]['nN_ISM']
b_ISM   = chain[0]['b_ISM']
T_ISM   = chain[0]['T_ISM']
xi_ISM  = chain[0]['xi_ISM']

nN_CS   = chain[0]['nN_CS']
nS_CS   = chain[0]['nS_CS']
b_CS    = chain[0]['b_CS']
T_CS    = chain[0]['T_CS']
xi_CS   = chain[0]['xi_CS']

nN_X1   = chain[0]['nN_X1']
nS_X1   = chain[0]['nS_X1']
b_X1    = chain[0]['b_X1']
T_X1    = chain[0]['T_X1']
xi_X1   = chain[0]['xi_X1']
RV_X1   = chain[0]['RV_X1']

nN_X2   = chain[0]['nN_X2']
nS_X2   = chain[0]['nS_X2']
b_X2    = chain[0]['b_X2']
T_X2    = chain[0]['T_X2']
xi_X2   = chain[0]['xi_X2']
RV_X2   = chain[0]['RV_X2']

for i in range(number_of_chains-1):
  Chi2    = np.concatenate((Chi2,stat[i]['Chi2']))

  nN_ISM  = np.concatenate((nN_ISM,chain[i]['nN_ISM']))
  b_ISM   = np.concatenate((b_ISM,chain[i]['b_ISM']))
  T_ISM   = np.concatenate((T_ISM,chain[i]['T_ISM']))
  xi_ISM  = np.concatenate((xi_ISM,chain[i]['xi_ISM']))
  
  nN_CS   = np.concatenate((nN_CS,chain[i]['nN_CS']))
  nS_CS   = np.concatenate((nS_CS,chain[i]['nS_CS']))
  b_CS    = np.concatenate((b_CS,chain[i]['b_CS']))
  T_CS    = np.concatenate((T_CS,chain[i]['T_CS']))
  xi_CS   = np.concatenate((xi_CS,chain[i]['xi_CS']))
  
  nN_X1   = np.concatenate((nN_X1,chain[i]['nN_X1']))
  nS_X1   = np.concatenate((nS_X1,chain[i]['nS_X1']))
  b_X1    = np.concatenate((b_X1,chain[i]['b_X1']))
  T_X1    = np.concatenate((T_X1,chain[i]['T_X1']))
  xi_X1   = np.concatenate((xi_X1,chain[i]['xi_X1']))
  RV_X1   = np.concatenate((RV_X1,chain[i]['RV_X1']))
  
  nN_X2   = np.concatenate((nN_X2,chain[i]['nN_X2']))
  nS_X2   = np.concatenate((nS_X2,chain[i]['nS_X2']))
  b_X2    = np.concatenate((b_X2,chain[i]['b_X2']))
  T_X2    = np.concatenate((T_X2,chain[i]['T_X2']))
  xi_X2   = np.concatenate((xi_X2,chain[i]['xi_X2']))
  RV_X2   = np.concatenate((RV_X2,chain[i]['RV_X2']))


k         = 1.38064852e-23    # Boltzmann constant in J/K = m^2*kg/(s^2*K) in SI base units
u         = 1.660539040e-27   # Atomic mass unit (Dalton) in kg

b_ISM     = 0.12895223*np.sqrt(T_ISM / 14.007  + (xi_ISM/0.12895223)**2)
b_CS      = 0.12895223*np.sqrt(T_CS / 14.007   + (xi_CS/0.12895223)**2)
b_X1      = 0.12895223*np.sqrt(T_X1 / 14.007    + (xi_X1/0.12895223)**2)
b_X2      = 0.12895223*np.sqrt(T_X2 / 14.007    + (xi_X2/0.12895223)**2)



# Arange the data into pandas format to be compatible with corner.py
'''
print Uncertainties(nN_ISM)
print Uncertainties(nN_CS)
print Uncertainties(nN_X1)
print Uncertainties(nN_X2)
print Uncertainties(nN_X3)

print "\nb:"
print Uncertainties(b_ISM)
print Uncertainties(b_CS)
print Uncertainties(b_X1)
print Uncertainties(b_X2)
print Uncertainties(b_X3)

#print "\nT:"
#print Uncertainties(T_ISM)
#print Uncertainties(T_CS)
#print Uncertainties(T_X)
'''
print "\nxi:"
print Uncertainties(xi_ISM)
print Uncertainties(xi_CS)
print Uncertainties(xi_X1)
print Uncertainties(xi_X2)


'''
print "\nSiIII**:"
print Uncertainties(nS_CS)
print Uncertainties(nS_X1)
print Uncertainties(nS_X2)
'''

print "\nRV:"
print Uncertainties(RV_X1)
print Uncertainties(RV_X2)
#sys.exit()

#data = np.array([nN_ISM,nN_CS,nN_X1,nN_X2,b_ISM,b_CS,b_X1,b_X2]).T
#columns = ['nN_ISM','nN_CS','nN_X1','nN_X2','b_ISM','b_CS','b_X1','b_X2']

data = np.array([nN_ISM,nN_CS,nN_X1,nN_X2,b_CS,b_X1,b_X2]).T
columns = ['nN_ISM','nN_CS','nN_X1','nN_X2','b_CS','b_X1','b_X2']



df = pd.DataFrame(data,columns=columns)

#data = np.array([nN_ISM,nN_CS,nN_X1,nN_X2,nN_X3,b_CS,b_X1,b_X2,b_X3,Chi2]).T
#columns = ['nN_ISM','nN_CS','nN_X1','nN_X2','nN_X3','b_CS','b_X1','b_X2','b_X3','Chi2']
#dfS = pd.DataFrame(data,columns=columns)


#df = df.drop(df[df.nN_ISM < 12].index)

#before = float(len(df['nN_ISM']))

#df = df[df['nN_ISM'] < 16.0]
#df = df[df['nN_CS'] < 16.0]
#df = df[df['b_ISM'] < 10.0]
#df = df[df['b_X2'] < 80.0]
#df = df[df['nN_X2'] > 11.5]
#df = df[df['nN_X3'] < 14.0]
#df = df[df['nN_CS'] < 16.0]
#df = df[df['b_X2'] < 20.0]
#df = df[df['b_X3'] < 20.0]

# Plot the posterior distributions.

fontlabel_size  = 14
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
params = {'backend': 'wxAgg', 'lines.markersize' : 2, 'axes.labelsize': fontlabel_size, 'font.size': fontlabel_size, 'legend.fontsize': fontlabel_size, 'font.family': 'Computer Modern'}
plt.rcParams.update(params)

# I'd like to have TeX font on the axis. Sadly the above line does not work.
#plt.rc('text', usetex=True)
#'''
figure = corner.corner(df,labels=[
r"$\log(N_{\mathrm{NI}})_{\mathrm{ISM}}$",\
r"$\log(N_{\mathrm{NI}})_{\mathrm{CS}}$",\
r"$\log(N_{\mathrm{NI}})_{\mathrm{X1}}$",\
r"$\log(N_{\mathrm{NI}})_{\mathrm{X2}}$",\
#r"$\log(N_{\mathrm{NI}})_{\mathrm{X3}}$",\
#r"$b_\mathrm{ISM}$",\
r"$b_\mathrm{CS}$",\
r"$b_\mathrm{X1}$",\
r"$b_\mathrm{X2}$"],
                                     quantiles=[0.16, 0.5,0.8413],
                                     levels=(1-np.exp(-0.5),),
                                     #truths=[13.49,17.15,14.03,0.77,1.33,2.32,43.70],
                                     #range=[(12,17),(13.40,17.8),(13.8,14.30),(0,2),(0.0,4.0),(0,36),(0,14)],
                                     #fill_contours=True,
                                     ret=True,
                                     bins=50,
                                     smooth=0.8,
                                     show_titles=True, title_kwargs={"fontsize": 13},
                                     label_kwargs = {"fontsize": 22},
                                     plot_contours=True,
                                     verbose=False,
                                     use_math_text=True)

figure.savefig("plots/mcmc_"+letter+".pdf")