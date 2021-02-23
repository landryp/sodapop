#!/usr/bin/env python

import numpy as np
from scipy.stats import uniform, norm

### BASIC PRIOR DISTRIBUTIONS

def flat(size=1,lb=0.,ub=1.):

	return uniform.rvs(loc=lb,scale=ub-lb,size=int(size))

def normal(size=1,med=0.,std=1.):

	return norm.rvs(loc=med,scale=std,size=int(size))
	
### PRIOR LOOKUP AND SAMPLING FUNCTIONS
	
param_priors = {'flat': flat, 'norm': normal}

def get_param_prior(key):

	try: prior_func = param_priors[key]
	except KeyError:
		print('Invalid prior specification, accepted keys are as follows:\n{0}\n'.format(param_priors.keys()))
		raise KeyError

	return prior_func

def get_param_samples(size=1,distr='flat',params=None):

	hyperparams = [float(param) for param in params]
	
	param_prior = get_param_prior(distr)	
	samps = param_prior(size,*hyperparams)

	return samps

