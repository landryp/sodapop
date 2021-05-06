#!/usr/bin/env python

import numpy as np
from scipy.stats import uniform, norm

### BASIC PRIOR DISTRIBUTIONS

def flat_prior(x,lb=0.,ub=1.):

	if x > ub or x < lb: val = 0.
	else: val = 1./(ub-lb)

	return val
	
def flat12_prior(x,y,lbx=0.,ubx=1.,lby=0.,uby=1.):

	if x > ubx or x < lbx or y > uby or y < lby or x > y: val = 0.
	else: val = 1./((ubx-lbx)*(uby-lby))

	return val

def normal_prior(x,mu=0.,sigma=1.):

	val = np.exp(-((m-mu)/(np.sqrt(2)*sigma))**2)/(sigma*np.sqrt(2*np.pi))
    
	return val
	
### BASIC PRIOR DISTRIBUTIONS SAMPLERS

def flat(size=1,lb=0.,ub=1.):

	return uniform.rvs(loc=lb,scale=ub-lb,size=int(size))

def normal(size=1,med=0.,std=1.):

	return norm.rvs(loc=med,scale=std,size=int(size))
	
### PRIOR LOOKUP AND SAMPLING FUNCTIONS
	
param_priors = {'flat': flat, 'norm': normal}
param_prior_funcs = {'flat': flat_prior, 'flat12': flat12_prior, 'norm': normal_prior}

def get_param_prior_func(key):

	try: prior_func = param_prior_funcs[key]
	except KeyError:
		print('Invalid prior specification, accepted keys are as follows:\n{0}\n'.format(param_prior_funcs.keys()))
		raise KeyError

	return prior_func

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

