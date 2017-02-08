'''
This code is used to submit tasks to exoatmos.
Written by: Dr. Paul A. Wilson (paul.wilson@iap.fr)
'''

from condorpy import Job, Templates
import numpy as np

cores = np.arange(20)						# Run on 20 cores as there is 20 cores per node
nodes = ['"exa01.iap.fr"','"exa02.iap.fr"','"exa03.iap.fr"','"exa04.iap.fr"',\
'"exa11.iap.fr"','"exa12.iap.fr"','"exa13.iap.fr"','"exa14.iap.fr"',\
'"exa21.iap.fr"','"exa22.iap.fr"','"exa23.iap.fr"','"exa24.iap.fr"']

job = Job("MCMC")
for i in range(len(nodes)):
	for j in range(len(cores)):
		print "Running chain #",j+1,"\ton ",nodes[i]
		job.executable = "MCMC_XA.py"
		job.arguments = "c"+str(i)+"_n"+str(j+1)
		job.requirements = "Machine=="+str(nodes[i])
		job.submit()
