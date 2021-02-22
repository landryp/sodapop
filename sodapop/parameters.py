#!/usr/bin/env python

import numpy as np
from scipy.stats import uniform, norm

def flat(size=1,lb=0.,ub=1.):

	return uniform.rvs(loc=lb,scale=ub-lb,size=int(size))

def normal(size=1,med=0.,std=1.):

	return norm.rvs(loc=med,scale=std,size=int(size))
	
param_priors = {'flat': flat, 'norm': normal}

def get_param_samples(size=1,distr='flat',params=None):

	hyperparams = [float(param) for param in params]
	
	samps = param_priors[distr](size,*hyperparams)

	return samps

