#!/usr/bin/env python

import numpy as np
from scipy.stats import uniform, norm, powerlaw, gaussian_kde
import os

### LOAD EOS-INFORMED MMAX DATA

MMAX_PATH = './sodapop/etc/LCEHL_NS_observables_samples_5k.csv'
MMAX_KEY = 'Mmax'

if not os.path.isfile(MMAX_PATH):

	import requests as rq
	import zipfile

	MTOV_PATH = 'https://zenodo.org/record/5397808/files/NS_samples.zip?download=1'
	MTOV_NAME = 'LCEHL_NS_observables_samples'
	
	mtov_file = rq.get(MTOV_PATH,allow_redirects=True)
	with open(MTOV_NAME+'.zip','wb') as outfile:
		outfile.write(mtov_file.content)

	with zipfile.ZipFile(MTOV_NAME+'.zip') as infile:
		infile.extractall('./')
        
	os.remove(MTOV_NAME+'.zip')
	MMAX_PATH = 'NS_samples.csv'
	MMAX_KEY = 'M_max'

### BASIC PRIOR DISTRIBUTIONS

def flat_prior(x,lb=0.,ub=1.):

	if x > ub or x < lb: val = 0.
	else: val = 1./(ub-lb)

	return val
	
def flat12_prior(x,y,lbx=0.,ubx=1.,lby=0.,uby=1.):

	if x > ubx or x < lbx or y > uby or y < lby or x > y: val = 0.
	else: val = 1./((ubx-lbx)*(uby-lby))

	return val

def flat123_prior(x,y,z,lbx=0.,ubx=1.,lby=0.,uby=1.,lbz=0.,ubz=1.):

	if x > ubx or x < lbx or y > uby or y < lby or x > y or z > min(ubz,y) or z < max(lbz,x): val = 0.
	else: val = 1./((ubx-lbx)*(uby-lby)*(ubz-lbz))

	return val
	
def flat1234_prior(x,y,z,w,lbx=0.,ubx=1.,lby=0.,uby=1.,lbz=0.,ubz=1.,lbw=0.,ubw=1.):

	if x > ubx or x < lbx or y > uby or y < lby or x > y or z > min(ubz,y) or z < max(lbz,x) or w > min(ubw,y) or w < max(lbw,x) or z > w: val = 0.
	else: val = 1./((ubx-lbx)*(uby-lby)*(ubz-lbz)*(ubw-lbw))

	return val

def normal_prior(x,mu=0.,sigma=1.):

	val = np.exp(-((m-mu)/(np.sqrt(2)*sigma))**2)/(sigma*np.sqrt(2*np.pi))
    
	return val
	
def quad_prior(x,lb=0.,ub=1.):

	if np.isscalar(x): x = np.array([x])
	else: x = np.array(x)
	z = np.zeros(len(x))
	
	p = 3.*((x-lb)/(ub-lb))**2/(ub-lb)
	
	return np.where((x < lb) | (x > ub), z, p)
	
### BASIC PRIOR DISTRIBUTIONS SAMPLERS

def flat(size=1,lb=0.,ub=1.):

	return uniform.rvs(loc=lb,scale=ub-lb,size=int(size))

def normal(size=1,med=0.,std=1.):

	return norm.rvs(loc=med,scale=std,size=int(size))
	
def quad(size=1,lb=0.,ub=1.):

	return powerlaw.rvs(3.,loc=lb,scale=ub-lb,size=int(size))

### EOS-INFORMED MMAX DISTRIBUTION

mmax_dat = np.genfromtxt(MMAX_PATH,names=True,dtype=None,encoding=None,delimiter=',')
	
def mmax_prior(x,mtovs=mmax_dat[MMAX_KEY],min_mmax=1.5,max_mtov=5.):

	if x < min_mmax or x > max_mtov: val = 0.

	else:
		mmax_grid = np.arange(min_mmax,max_mtov,0.01)
		val = len(mtovs[mtovs >= x])/len(mtovs) # require mmax <= mtov

	return val
	
def flatmmin_mmax_prior(x,y,lbx=0.,ubx=1.,mtovs=mmax_dat[MMAX_KEY],min_mmax=1.5,max_mtov=5.):

	if x > ubx or x < lbx or x > y or y < min_mmax or y > max_mtov: val = 0.
	else: val = mmax_prior(y,mtovs,min_mmax,max_mtov)/(ubx-lbx)

	return val

def flatmminmu_mmax_prior(x,y,z,lbx=0.,ubx=1.,lbz=0.,ubz=1.,mtovs=mmax_dat[MMAX_KEY],min_mmax=1.5,max_mtov=5.):

	if x > ubx or x < lbx or x > y or z > min(ubz,y) or z < max(lbz,x) or y < min_mmax or y > max_mtov: val = 0.
	else: val = mmax_prior(y,mtovs,min_mmax,max_mtov)/((ubx-lbx)*(ubz-lbz))

	return val
	
def flatmminmu1mu2_mmax_prior(x,y,z,w,lbx=0.,ubx=1.,lbz=0.,ubz=1.,lbw=0.,ubw=1.,mtovs=mmax_dat[MMAX_KEY],min_mmax=1.5,max_mtov=5.):

	if x > ubx or x < lbx or x > y or z > min(ubz,y) or z < max(lbz,x) or w > min(ubw,y) or w < max(lbw,x) or z > w or y < min_mmax or y > max_mtov: val = 0.
	else: val = mmax_prior(y,mtovs,min_mmax,max_mtov)/((ubx-lbx)*(ubz-lbz)*(ubw-lbw))

	return val
	
### PRIOR LOOKUP AND SAMPLING FUNCTIONS
	
param_priors = {'flat': flat, 'norm': normal, 'quad': quad}
param_prior_funcs = {'flat': flat_prior, 'flat12': flat12_prior, 'flat123': flat123_prior, 'flat1234': flat1234_prior, 'norm': normal_prior, 'quad': quad_prior, 'mmax': mmax_prior, 'flatmmin_mmax': flatmmin_mmax_prior, 'flatmminmu_mmax': flatmminmu_mmax_prior, 'flatmminmu1mu2_mmax': flatmminmu1mu2_mmax_prior}

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
	
def load_dist_prior(dist_prior_str):

	dist_shape = (dist_prior_str).split(',')[0]
	dist_func = get_param_prior_func(dist_shape)
	dist_params = [float(val) for val in (dist_prior_str).split(',')[1:]]

	return dist_func, dist_params
	
def create_prior_dict(prior_str_list):

	prior_dict = {}

	for prior_str in prior_str_list:

		prior_key = prior_str.split('+')[0]
		prior_func = get_param_prior_func((prior_str.split('+')[1]).split(',')[0])
		hyp_params = [float(val) for val in (prior_str.split('+')[1]).split(',')[1:]]
		prior_dict[prior_key] = [prior_func,hyp_params]

	return prior_dict
	
def create_fixed_dict(fixed_str_list):	

	fixed_dict = {}
	fixed_params = [fixed_str.split(',')[0] for fixed_str in fixed_str_list]
	fixed_vals = [fixed_str.split(',')[1] for fixed_str in fixed_str_list]
	for i,fixed_param in enumerate(fixed_params):
		fixed_dict[fixed_param] = float(fixed_vals[i])

	return fixed_params, fixed_dict

def create_lambda_dict(pop_param_names,fixed_params=[],fixed_dict={}):

	lambda_dict = {}
	j = 0
	for param in pop_param_names:
		if param in fixed_params:
			lambda_dict[param] = lambda x, param=param : fixed_dict[param]
		else:
			lambda_dict[param] = lambda x, j=j : x[j]
			j += 1

	return lambda_dict
	
