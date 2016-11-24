'''
This code is used to submit tasks to exoatmos.
Written by: Dr. Paul A. Wilson (paul.wilson@iap.fr)
'''

from condorpy import Job, Templates
import numpy as np

number_of_chains = 100

job = Job("JobName")
for i in range(number_of_chains):
	print "Running chain #",i+1
	job.executable = "MCMC_XA.py"
	job.arguments = str(i+1)
	job.submit()