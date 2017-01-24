import numpy as np

from src.statistics import Stats
from src.model import Model

s   = Stats()
m   = Model()

class MCMC:

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
      stats     = np.zeros(shape=(int(C),1))
      
      for i in range(int(C)):
        if i%100 == 0.:
          print (i/C)*100.," % done"
        jump        = np.random.normal(0.,1.,len(S)) * S
        P           = P + jump
        
        while P[0] < 12.0 or P[5] < 12.0:
          P         = P - jump
          jump      = np.random.normal(0.,1.,len(S)) * S
          P         = P + jump


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

        # Calculate Chi2
        stats[i] = s.chi2(X)
      
      print "\nAccepted steps: ",round(100.*(moves/C),2),"%"
      
      return chain, moves, stats