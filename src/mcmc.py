import numpy as np

from src.statistics import Stats
from src.model import Model

s   = Stats()
m   = Model()

class MCMC:

    def Distrib(self, x):
       '''Finds median and 68% interval of array x.'''
       y    = sorted(x)
       up   = y[int(0.8413*len(y))]
       down = y[int(0.1587*len(y))]
       med  = y[int(0.5*len(y))]
       
       return med,up,down   
        
    def Median_and_Uncertainties(self, P, S, chain):
      ''' Prints the median parameter values
      with their associated uncertainties '''
      param_ans = []
      param_u   = []
      param_l   = []

      for i in range(len(P)):
        up_err      = self.Distrib(chain[:,P[i]])[1]-self.Distrib(chain[:,P[i]])[0]
        median_val  = self.Distrib(chain[:,P[i]])[0]
        low_err     = self.Distrib(chain[:,P[i]])[0]-self.Distrib(chain[:,P[i]])[2]

        param_ans.append(median_val)
        param_u.append(up_err)
        param_l.append(low_err)

      return param_ans,param_u,param_l

    def McMC(self, x, X, F, ModelType, param, P, Const, S, C):
      '''
      x => x-axis values (In this case wavelength)
      X => Data (y,yerr,model)
      F => Function used
      P => Parameters
      S => Scale
      C => Chain length
      '''
      L         = s.Merit(X)
      moves     = 0
      chain     = np.zeros(shape=(int(C),len(P)))
      L_chain   = np.zeros(shape=(int(C),1))
      for i in range(int(C)):
        if i%100 == 0.:
          print (i/C)*100.," % done"
        jump        = np.random.normal(0.,1.,len(S)) * S
        P           = P + jump
        new_fit     = m.Model(P, Const, ModelType, param)[0]
        X           = X[0],X[1],new_fit
        L_new       = s.Merit(X)
        L_chain[i]  = L_new
        ratio       = L_new/L

        if (np.random.random() > ratio):
          P     = P - jump
          moved = 0
        else:
          L     = L_new
          moved = 1
        moves  += moved
        chain[i,:] = np.array(P)
      print "\nAccepted steps: ",round(100.*(moves/C),2),"%"
      
      return chain, moves