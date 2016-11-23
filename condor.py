'''
This code is used to submit tasks to exoatmos.
Written by: Dr. Paul A. Wilson (paul.wilson@iap.fr)
'''

from condorpy import Job, Templates
import numpy as np

# A list of cores used
cores = np.arange(20)+1

# A list of nodes used:
#nodes = ['exa01','exa02','exa03','exa04','exa11','exa12','exa13','exa14','exa21','exa22','exa23','exa24']
nodes = ['exa14','exa21']

for i in range(len(nodes)):
	for j in range(len(cores)):
		print "running core ",j+1," on node ",nodes[i]+".iap.fr"
		job = Job('MCMC_NI_'+str(nodes[i])+"_"+str(j+1),attributes=Templates.base)
		job.executable = "MCMC_XA.py"
		job.arguments = str(j)
		job.requirements = "Machine=="+nodes[i]+".iap.fr"
		job.submit()