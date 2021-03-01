#!/usr/bin/env python

import numpy as np
import scipy.special

### BASIC ANALYTIC DISTRIBUTIONS

def powerlaw(m,alpha,mmax):

	if m > mmax: val = 0.
	else: val = m**alpha
        
	return val

def gaussian(m,mu,sigma):

	val = np.exp(-((m-mu)/(np.sqrt(2)*sigma))**2)/(sigma*np.sqrt(2*np.pi))
    
	return val 

def smooth(m,mmin,delta):

	if m < mmin: val =  0.
	elif m >= mmin + delta: val = 1.
	else: val = 1./(np.exp(delta/(m-mmin)+delta/(m-mmin-delta))+1.)
    
	return val

### INDIVIDUAL MASS DISTRIBUTIONS

def unif_mass(m,mmin=1.,mmax=3.): # uniform mass distribution between mmin and mmax

	if m > mmax or m < mmin: val = 0.
	else: val = 1./(mmax-mmin)

	return val
	
def peak_mass(m,mu=1.34,sigma=0.02): # gaussian mass distribution

	val = gaussian(m,mu,sigma)

	return val
	
def bimod_mass(m,mu1=1.34,sigma1=0.02,mu2=1.47,sigma2=0.15,alpha=0.68): # double Gaussian fit to Galactic double neutron star population from Farrow+ arXiv:1902.03300

	val = alpha*gaussian(m,mu1,sigma1) + (1.-alpha)*gaussian(m,mu2,sigma2)

	return val
	
def peakcut_mass(m,mu=1.34,sigma=0.02,mmin=1.,mmax=3.): # gaussian mass distribution with high- and low-mass cutoffs

	if m > mmax or m < mmin: val = 0.
	else: val = gaussian(m,mu,sigma)
	
	return val
	
def bimodcut_mass(m,mu1=1.34,sigma1=0.02,mu2=1.47,sigma2=0.15,alpha=0.68,mmin=1.,mmax=3.): # double gaussian mass distribution with high- and low-mass cutoffs

	if m > mmax or m < mmin: val = 0.
	else: val = alpha*gaussian(m,mu1,sigma1) + (1.-alpha)*gaussian(m,mu2,sigma2)
	
	return val

### BINARY MASS DISTRIBUTIONS

def unif_m1m2(m1,m2,mmin=5.,mmax=2e2): # uniform distribution in source frame masses, subject to m1 >= m2 convention

	if m1 < m2 or m1 > mmax or m2 < mmin: val = 0.
	else: val = 1./(mmax-mmin)**2
    
	return val
	
def peak_m1m2(m1,m2,mu=1.34,sigma=0.02): # gaussian distribution in source frame masses, subject to m1 >= m2 convention

	if m1 < m2: val = 0.
	else: val = gaussian(m1,mu,sigma)*gaussian(m2,mu,sigma)
	
	return val
	
def bimod_m1m2(m1,m2,mu1=1.34,sigma1=0.02,mu2=1.47,sigma2=0.15,alpha=0.68): # double gaussian distribution in source frame masses, subject to m1 >= m2 convention

	if m1 < m2: val = 0.
	else: val = (alpha*gaussian(m1,mu1,sigma1) + (1.-alpha)*gaussian(m1,mu2,sigma2))*(alpha*gaussian(m2,mu1,sigma1) + (1.-alpha)*gaussian(m2,mu2,sigma2))
	
	return val
	
def peakcut_m1m2(m1,m2,mu=1.34,sigma=0.02,mmin=1.,mmax=3.): # gaussian distribution in source frame masses with high- and low-mass cutoffs, subject to m1 >= m2 convention

	if m1 < m2 or m1 > mmax or m2 < mmin: val = 0.
	else: val = gaussian(m1,mu,sigma)*gaussian(m2,mu,sigma)
	
	return val
	
def bimodcut_m1m2(m1,m2,mu1=1.34,sigma1=0.02,mu2=1.47,sigma2=0.15,alpha=0.68,mmin=1.,mmax=3.): # double gaussian mass distribution with high- and low-mass cutoffs

	if m1 < m2 or m1 > mmax or m2 < mmin: val = 0.
	else: val = (alpha*gaussian(m1,mu1,sigma1) + (1.-alpha)*gaussian(m1,mu2,sigma2))*(alpha*gaussian(m2,mu1,sigma1) + (1.-alpha)*gaussian(m2,mu2,sigma2))
	
	return val
	
def peak_m1_unif_m2(m1,m2,mu=1.34,sigma=0.02,mmin=1.,mmax=3.): # gaussian distribution in source frame masses with high- and low-mass cutoffs, subject to m1 >= m2 convention

	if m1 < m2: val = 0.
	else: val = gaussian(m1,mu,sigma)*unif_mass(m2,mmin,mmax)
	
	return val
	
def bimod_m1_unif_m2(m1,m2,mu1=1.34,sigma1=0.02,mu2=1.47,sigma2=0.15,alpha=0.68,mmin=1.,mmax=3.): # double gaussian mass distribution with high- and low-mass cutoffs

	if m1 < m2: val = 0.
	else: val = (alpha*gaussian(m1,mu1,sigma1) + (1.-alpha)*gaussian(m1,mu2,sigma2))*unif_mass(m2,mmin,mmax)
	
	return val
	
def unif_m1m2_qpair(m1,m2,mmin=5.,mmax=2e2,beta=0.): # uniform distribution in source frame masses, subject to m1 >= m2 convention and q-dependent pairing

	if m1 < m2 or m1 > mmax or m2 < mmin: val = 0.
	else: val = (m2/m1)**beta
    
	return val
	
def peakcut_m1m2_qpair(m1,m2,mu=1.34,sigma=0.02,mmin=1.,mmax=3.,beta=0.): # gaussian distribution in source frame masses, subject to m1 >= m2 convention and q-dependent pairing

	if m1 < m2 or m1 > mmax or m2 < mmin: val = 0.
	else: val = gaussian(m1,mu,sigma)*gaussian(m2,mu,sigma)*(m2/m1)**beta
	
	return val
	
def bimodcut_m1m2_qpair(m1,m2,mu1=1.34,sigma1=0.02,mu2=1.47,sigma2=0.15,alpha=0.68,mmin=1.,mmax=3.,beta=0.): # double gaussian distribution in source frame masses, subject to m1 >= m2 convention and q-dependent pairing

	if m1 < m2 or m1 > mmax or m2 < mmin: val = 0.
	else: val = (alpha*gaussian(m1,mu1,sigma1) + (1.-alpha)*gaussian(m1,mu2,sigma2))*(alpha*gaussian(m2,mu1,sigma1) + (1.-alpha)*gaussian(m2,mu2,sigma2))*q**beta
	
	return val
	
def o3a_powerpeak_m1m2_qpair(m1,m2,lpeak=0.10,alpha=2.62,beta=1.26,mmin=4.53,mmax=86.73,delta=4.88,mu=33.49,sigma=5.09): # power law + peak mass model from O3a populations paper
    
	val = ((1.-lpeak)*powerlaw(m1,-alpha,mmax)+lpeak*gaussian(m1,mu,sigma))*smooth(m1,mmin,delta)*(m2/m1)**beta*smooth(m2,mmin,delta)/m1
    
	return val

def o3a_powerbreak_m1m2_qpair(m1,m2,alpha1=1.58,alpha2=5.59,beta=1.40,mmin=4.83,mmax=3.96,delta=87.14,b=0.43): # broken power law mass model from O3a populations paper
    
    mbreak = mmin+b*(mmax-mmin)
    
    if m1 >= mmin and m1 < mbreak: val = m1**(-alpha1)*smooth(m1,mmin,delta)*(m2/m1)**beta*s(m2,mmin,delta)/m1
    elif m1 >= mbreak and m1 <= mmax: val = m1**(-alpha2)*smooth(m1,mmin,delta)*(m2/m1)**beta*s(m2,mmin,delta)/m1
    else: val = 0.
    
    return val
    
### PRIOR LOOKUP FUNCTIONS
	
pop_priors = {'unif_mass': unif_mass, 'peak_mass': peak_mass, 'bimod_mass': bimod_mass, 'peakcut_mass': peakcut_mass, 'bimodcut_mass': bimodcut_mass, 'unif_m1m2': unif_m1m2, 'peak_m1m2': peak_m1m2, 'bimod_m1m2': bimod_m1m2, 'peakcut_m1m2': peakcut_m1m2, 'bimodcut_m1m2': bimodcut_m1m2, 'peak_m1_unif_m2': peak_m1_unif_m2, 'bimod_m1_unif_m2': bimod_m1_unif_m2, 'unif_m1m2_qpair': unif_m1m2_qpair, 'peakcut_m1m2_qpair': peakcut_m1m2_qpair, 'bimodcut_m1m2_qpair': bimodcut_m1m2_qpair, 'o3a_powerpeak_m1m2_qpair': o3a_powerpeak_m1m2_qpair, 'o3a_powerbreak_m1m2_qpair': o3a_powerbreak_m1m2_qpair}

pop_params = {'unif_mass': 'mmin,mmax', 'peak_mass': 'mu,sigma', 'bimod_mass': 'mu1,sigma1,mu2,sigma2,alpha', 'peakcut_mass': 'mu,sigma,mmin,mmax', 'bimodcut_mass': 'mu1,sigma1,mu2,sigma2,alpha,mmin,mmax', 'unif_m1m2': 'mmin,mmax', 'peak_m1m2': 'mu,sigma', 'bimod_m1m2': 'mu1,sigma1,mu2,sigma2,alpha', 'peakcut_m1m2': 'mu,sigma,mmin,mmax', 'bimodcut_m1m2': 'mu1,sigma1,mu2,sigma2,alpha,mmin,mmax', 'peak_m1_unif_m2': 'mu,sigma,mmin,mmax', 'bimod_m1_unif_m2': 'mu1,sigma1,mu2,sigma2,alpha,mmin,mmax', 'unif_m1m2_qpair': 'mmin,mmax,beta', 'peakcut_m1m2_qpair': 'mu,sigma,mmin,mmax,beta', 'bimodcut_m1m2_qpair': 'mu1,sigma1,mu2,sigma2,alpha,mmin,mmax,beta', 'o3a_powerpeak_m1m2_qpair': 'lpeak,alpha,beta,mmin,mmax,delta,mu,sigma', 'o3a_powerbreak_m1m2_qpair': 'alpha1,alpha2,beta,mmin,mmax,delta,b'}

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
