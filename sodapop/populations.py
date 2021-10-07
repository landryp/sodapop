#!/usr/bin/env python

import numpy as np
import scipy.special
#import scipy.integrate

### BASIC ANALYTIC DISTRIBUTIONS

def uniform(m,lambdaa):

	mmin, mmax = lambdaa[:2]
	
	if np.isscalar(m): m = np.array([m])
	else: m = np.array(m)
	z = np.zeros(len(m))
	
	p = np.full(len(m),1./(mmax-mmin))
	
	return np.where((m > mmax) | (m < mmin), z, p)

def powerlaw(m,lambdaa):

	alpha, mmin, mmax = lambdaa[:3]
	
	if np.isscalar(m): m = np.array([m])
	else: m = np.array(m)
	z = np.zeros(len(m))
	
	if alpha != -1.: p = (1.+alpha)*m**alpha/(mmax**(1.+alpha)-mmin**(1.+alpha))
	else: p = m**alpha/(np.log(mmax)-np.log(mmin))
	
	return np.where((m > mmax) | (m < mmin), z, p)

def gaussian(m,lambdaa):

	mu, sigma = lambdaa[:2]
	
	if np.isscalar(m): m = np.array([m])
	else: m = np.array(m)
	
	p = np.exp(-((m-mu)/(np.sqrt(2)*sigma))**2)/(sigma*np.sqrt(2*np.pi))
	
	return p

def smooth(m,lambdaa):

	mmin, delta = lambdaa[:2]
	
	if np.isscalar(m): m = np.array([m])
	else: m = np.array(m)
	z = np.zeros(len(m))
	o = np.ones(len(m))
	
	p = 1./(np.exp(delta/(m-mmin)+delta/(m-mmin-delta))+1.)
	
	return np.select([m < mmin, (m >= mmin) & (m < mmin + delta), m >= mmin + delta], [z, p, o])

### INDIVIDUAL MASS DISTRIBUTIONS

def unif_mass(m,lambdaa): # uniform mass distribution

	p = uniform(m,lambdaa)

	return p
	
def peak_mass(m,lambdaa): # gaussian mass distribution

	p = gaussian(m,lambdaa)

	return p
	
def bimod_mass(m,lambdaa): # double gaussian mass distribution

	mu1, sigma1, mu2, sigma2, w = lambdaa[:5]
	
	norm1 = 1.
	norm2 = 1.
	p = w*gaussian(m,(mu1,sigma1))/norm1 + (1.-w)*gaussian(m,(mu2,sigma2))/norm2

	return p
	
def power_mass(m,lambdaa): # power law mass distribution

	p = powerlaw(m,lambdaa)

	return p
	
def peakcut_mass(m,lambdaa): # gaussian mass distribution with high- and low-mass cutoffs

	mu, sigma, mmin, mmax = lambdaa[:4]
	
	if np.isscalar(m): m = np.array([m])
	else: m = np.array(m)
	z = np.zeros(len(m))
	
	norm = 0.5*(scipy.special.erf((mmax-mu)/(np.sqrt(2)*sigma))-scipy.special.erf((mmin-mu)/(np.sqrt(2)*sigma)))
	p = gaussian(m,(mu,sigma))/norm
	
	return np.where((m > mmax) | (m < mmin), z, p)
	
def bimodcut_mass(m,lambdaa): # double gaussian mass distribution with high- and low-mass cutoffs

	mu1, sigma1, mu2, sigma2, w, mmin, mmax = lambdaa[:7]
	
	if np.isscalar(m): m = np.array([m])
	else: m = np.array(m)
	z = np.zeros(len(m))
	
	norm1 = 0.5*(scipy.special.erf((mmax-mu1)/(np.sqrt(2)*sigma1))-scipy.special.erf((mmin-mu1)/(np.sqrt(2)*sigma1)))
	norm2 = 0.5*(scipy.special.erf((mmax-mu2)/(np.sqrt(2)*sigma2))-scipy.special.erf((mmin-mu2)/(np.sqrt(2)*sigma2)))
	p = w*gaussian(m,(mu1,sigma1))/norm1 + (1.-w)*gaussian(m,(mu2,sigma2))/norm2
	
	return np.where((m > mmax) | (m < mmin), z, p)

### BINARY MASS DISTRIBUTIONS

def unif_m1m2(m1,m2,lambdaa): # uniform distribution in masses, subject to m1 >= m2 convention

	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	z = np.zeros(len(m1))

	p = unif_mass(m1,lambdaa)*unif_mass(m2,lambdaa)
	
	return np.where(m1 < m2, z, p)
	
def peak_m1m2(m1,m2,lambdaa): # gaussian distribution in masses, subject to m1 >= m2 convention

	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	z = np.zeros(len(m1))

	p = peak_mass(m1,lambdaa)*peak_mass(m2,lambdaa)

	return np.where(m1 < m2, z, p)
	
def bimod_m1m2(m1,m2,lambdaa): # double gaussian distribution in masses, subject to m1 >= m2 convention

	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	z = np.zeros(len(m1))
	
	p = bimod_mass(m1,lambdaa)*bimod_mass(m2,lambdaa)
	
	return np.where(m1 < m2, z, p)

def power_m1m2(m1,m2,lambdaa): # power law in masses, subject to m1 >= m2 convention

	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	z = np.zeros(len(m1))
	
	p = power_mass(m1,lambdaa)*power_mass(m2,lambdaa)
	
	return np.where(m1 < m2, z, p)

def peakcut_m1m2(m1,m2,lambdaa): # gaussian distribution in masses with high- and low-mass cutoffs, subject to m1 >= m2 convention

	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	z = np.zeros(len(m1))

	p = peakcut_mass(m1,lambdaa)*peakcut_mass(m2,lambdaa)

	return np.where(m1 < m2, z, p)
	
def bimodcut_m1m2(m1,m2,lambdaa): # double gaussian distribution in masses with high- and low-mass cutoffs, subject to m1 >= m2 convention

	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	z = np.zeros(len(m1))
	
	p = bimodcut_mass(m1,lambdaa)*bimodcut_mass(m2,lambdaa)
	
	return np.where(m1 < m2, z, p)
	
def unif_m1m2_qpair(m1,m2,lambdaa): # uniform distribution in masses, subject to m1 >= m2 convention and q-dependent pairing

	beta = lambdaa[-1]
	
	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	
	p = unif_m1m2(m1,m2,lambdaa)*(m2/m1)**beta
    
	return p
	
def power_m1m2_qpair(m1,m2,lambdaa): # uniform distribution in masses, subject to m1 >= m2 convention and q-dependent pairing

	beta = lambdaa[-1]
	
	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	
	p = power_m1m2(m1,m2,lambdaa)*(m2/m1)**beta
    
	return p
	
def peakcut_m1m2_qpair(m1,m2,lambdaa): # gaussian distribution in masses, subject to m1 >= m2 convention and q-dependent pairing

	beta = lambdaa[-1]
	
	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	
	p = peakcut_m1m2(m1,m2,lambdaa)*(m2/m1)**beta
    
	return p
	
def bimodcut_m1m2_qpair(m1,m2,lambdaa): # double gaussian distribution in masses, subject to m1 >= m2 convention and q-dependent pairing

	beta = lambdaa[-1]
	
	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	
	p = bimodcut_m1m2(m1,m2,lambdaa)*(m2/m1)**beta
    
	return p
	
def power_m1_unif_m2_qpair(m1,m2,lambdaa): # power-law in m1, uniform in m2, subject to m1 >= m2 convention and q-dependent pairing

	beta = lambdaa[-1]
	
	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	z = np.zeros(len(m1))
	
	p = power_mass(m1,lambdaa)*unif_mass(m2,lambdaa[1:])*(m2/m1)**beta
    
	return np.where(m1 < m2, z, p)
	
def bimodcut_m1_unif_m2_qpair(m1,m2,lambdaa): # bimodal in m1, uniform in m2, subject to m1 >= m2 convention and q-dependent pairing

	beta = lambdaa[-1]
	
	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	z = np.zeros(len(m1))
	
	p = bimodcut_mass(m1,lambdaa)*unif_mass(m2,lambdaa[5:])*(m2/m1)**beta
    
	return np.where(m1 < m2, z, p)
	
# NSBH MASS DISTRIBUTIONS

def unif_m1_unif_m2(m1,m2,lambdaa): # uniform distributions in bh and ns masses, subject to m1 >= m2 convention

	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	z = np.zeros(len(m1))

	p = unif_mass(m1,lambdaa[2:])*unif_mass(m2,lambdaa)
	
	return np.where(m1 < m2, z, p)
	
def unif_m1_power_m2(m1,m2,lambdaa): # uniform distributions in bh mass and power law in ns mass, subject to m1 >= m2 convention

	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	z = np.zeros(len(m1))

	p = unif_mass(m1,lambdaa[3:])*power_mass(m2,lambdaa)
	
	return np.where(m1 < m2, z, p)
	
def unif_m1_peakcut_m2(m1,m2,lambdaa): # uniform distribution in bh masses and gaussian distribution in ns masses, subject to m1 >= m2 convention

	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	z = np.zeros(len(m1))

	p = unif_mass(m1,lambdaa[4:])*peakcut_mass(m2,lambdaa)
	
	return np.where(m1 < m2, z, p)
	
def unif_m1_bimodcut_m2(m1,m2,lambdaa): # uniform distribution in bh masses and double gaussian distribution in ns masses, subject to m1 >= m2 convention

	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	z = np.zeros(len(m1))

	p = unif_mass(m1,lambdaa[7:])*bimodcut_mass(m2,lambdaa)
	
	return np.where(m1 < m2, z, p)
	
def unif_m1_unif_m2_qpair(m1,m2,lambdaa): # uniform distributions in bh and ns masses, subject to m1 >= m2 convention and q-dependent pairing

	beta = lambdaa[-1]
	
	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	
	p = unif_m1_unif_m2(m1,m2,lambdaa)*(m2/m1)**beta
    
	return p
	
def unif_m1_power_m2_qpair(m1,m2,lambdaa): # uniform distributions in bh and ns masses, subject to m1 >= m2 convention and q-dependent pairing

	beta = lambdaa[-1]
	
	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	
	p = unif_m1_power_m2(m1,m2,lambdaa)*(m2/m1)**beta
    
	return p
	
def unif_m1_peakcut_m2_qpair(m1,m2,lambdaa): # uniform distribution in bh masses and gaussian distribution in ns masses, subject to m1 >= m2 convention and q-dependent pairing

	beta = lambdaa[-1]
	
	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	
	p = unif_m1_peakcut_m2(m1,m2,lambdaa)*(m2/m1)**beta
	
	return p
	
def unif_m1_bimodcut_m2_qpair(m1,m2,lambdaa): # uniform distribution in bh masses and double gaussian distribution in ns masses, subject to m1 >= m2 convention and q-dependent pairing

	beta = lambdaa[-1]
	
	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	
	p = unif_m1_bimodcut_m2(m1,m2,lambdaa)*(m2/m1)**beta
	
	return p

def power2_m1_unif_m2_qpair(m1,m2,lambdaa): # power-law in m1 with exponent -2, uniform in m2, subject to m1 >= m2 convention and q-dependent pairing

        beta = lambdaa[-1]

        if np.isscalar(m1): m1 = np.array([m1])
        else: m1 = np.array(m1)
        if np.isscalar(m2): m2 = np.array([m2])
        else: m2 = np.array(m2)
        z = np.zeros(len(m1))

        p = power_mass(m1,[-2.]+lambdaa[2:])*unif_mass(m2,lambdaa)*(m2/m1)**beta

        return np.where(m1 < m2, z, p)
	
# LOOKUP FUNCTIONS

pop_priors = {'unif_mass': unif_mass, 'peak_mass': peak_mass, 'bimod_mass': bimod_mass, 'peakcut_mass': peakcut_mass, 'bimodcut_mass': bimodcut_mass, 'unif_m1m2': unif_m1m2, 'peak_m1m2': peak_m1m2, 'bimod_m1m2': bimod_m1m2, 'power_m1m2': power_m1m2, 'peakcut_m1m2': peakcut_m1m2, 'bimodcut_m1m2': bimodcut_m1m2, 'unif_m1m2_qpair': unif_m1m2_qpair, 'power_m1m2_qpair': power_m1m2_qpair, 'peakcut_m1m2_qpair': peakcut_m1m2_qpair, 'bimodcut_m1m2_qpair': bimodcut_m1m2_qpair, 'unif_m1_unif_m2': unif_m1_unif_m2, 'unif_m1_unif_m2_qpair': unif_m1_unif_m2_qpair, 'unif_m1_power_m2': unif_m1_power_m2, 'unif_m1_power_m2_qpair': unif_m1_power_m2_qpair, 'unif_m1_peakcut_m2': unif_m1_peakcut_m2, 'unif_m1_peakcut_m2_qpair': unif_m1_peakcut_m2_qpair, 'unif_m1_bimodcut_m2': unif_m1_bimodcut_m2, 'unif_m1_bimodcut_m2_qpair': unif_m1_bimodcut_m2_qpair, 'power_m1_unif_m2_qpair': power_m1_unif_m2_qpair, 'bimodcut_m1_unif_m2_qpair': bimodcut_m1_unif_m2_qpair,'power2_m1_unif_m2_qpair': power2_m1_unif_m2_qpair}

pop_params = {'unif_mass': 'mmin,mmax', 'peak_mass': 'mu,sigma', 'bimod_mass': 'mu1,sigma1,mu2,sigma2,w', 'peakcut_mass': 'mu,sigma,mmin,mmax', 'bimodcut_mass': 'mu1,sigma1,mu2,sigma2,w,mmin,mmax', 'unif_m1m2': 'mmin,mmax', 'peak_m1m2': 'mu,sigma', 'bimod_m1m2': 'mu1,sigma1,mu2,sigma2,w', 'power_m1m2': 'alpha,mmin,mmax', 'peakcut_m1m2': 'mu,sigma,mmin,mmax', 'bimodcut_m1m2': 'mu1,sigma1,mu2,sigma2,w,mmin,mmax', 'unif_m1m2_qpair': 'mmin,mmax,beta', 'power_m1m2_qpair': 'alpha,mmin,mmax,beta', 'peakcut_m1m2_qpair': 'mu,sigma,mmin,mmax,beta', 'bimodcut_m1m2_qpair': 'mu1,sigma1,mu2,sigma2,w,mmin,mmax,beta', 'unif_m1_unif_m2': 'mmin,mmax', 'unif_m1_unif_m2_qpair': 'mmin,mmax,beta', 'unif_m1_power_m2': 'alpha,mmin,mmax', 'unif_m1_power_m2_qpair': 'alpha,mmin,mmax,beta', 'unif_m1_peakcut_m2': 'mu,sigma,mmin,mmax', 'unif_m1_peakcut_m2_qpair': 'mu,sigma,mmin,mmax,beta', 'unif_m1_bimodcut_m2': 'mu1,sigma1,mu2,sigma2,w,mmin,mmax', 'unif_m1_bimodcut_m2_qpair': 'mu1,sigma1,mu2,sigma2,w,mmin,mmax,beta', 'power_m1_unif_m2_qpair': 'alpha,mmin,mmax,beta', 'bimodcut_m1_unif_m2_qpair': 'mu1,sigma1,mu2,sigma2,w,mmin,mmax,beta', 'power2_m1_unif_m2_qpair': 'mmin,mmax,beta'}

def get_pop_prior(key):

	try: prior_func = pop_priors[key]
	except KeyError:
		print('Invalid prior specification, accepted keys are as follows:\n{0}\n'.format(pop_priors.keys()))
		raise KeyError

	return prior_func

def get_pop_params(key):

	try: prior_func = pop_params[key]
	except KeyError:
		print('Invalid prior specification, accepted keys are as follows:\n{0}\n'.format(pop_params.keys()))
		raise KeyError

	return prior_func
