#!/usr/bin/env python

import numpy as np
import scipy.special
import scipy.integrate

### BASIC ANALYTIC DISTRIBUTIONS

def powerlaw(m,alpha,mmin,mmax):

	if m > mmax or m < mmin: val = 0.
	else: val = (1.+alpha)*m**alpha/(mmax**(1.+alpha)-mmin**(1.+alpha))
        
	return val

def gaussian(m,mu,sigma):

	val = np.exp(-((m-mu)/(np.sqrt(2)*sigma))**2)/(sigma*np.sqrt(2*np.pi))
    
	return val 

def smooth(m,mmin,delta):

	if m < mmin: val =  0.
	elif m >= mmin + delta: val = 1.
	else: val = 1./(np.exp(delta/(m-mmin)+delta/(m-mmin-delta))+1.)
    
	return val

### SELECTION FUNCTIONS

def snrcut(m1,m2,dl,snr0=8.,r0=176.): # snr selection function from https://arxiv.org/abs/1108.5161 w/ Advanced LIGO 3-IFO estimate of r0 [Mpc]

	Mc = (m1*m2)**(3./5.)/(m1+m2)**(1./5.)
	z = 67.4*dl/(2.998e5) # assume local Hubble's law w/ Planck 2018 cosmology
	Mcdet = Mc*(1.+z)
	x = (snr0/8.)*(dl/r0)*(1.2/Mcdet)**(5./6.)

	if x <= 0: C = 1.
	elif x > 0. and x <= 4.: C = (1.+x)*(4.-x)**4./256.
	else: C = 0.

	return C

select_funcs = {'snrcut': snrcut}

def get_select_func(key):

	try: select_func = select_funcs[key]
	except KeyError:
		print('Invalid selection function specification, accepted keys are as follows:\n{0}\n'.format(select_funcs.keys()))
		raise KeyError

	return select_func

### INDIVIDUAL MASS DISTRIBUTIONS

def unif_mass(m,mmin=1.,mmax=3.): # uniform mass distribution between mmin and mmax

	if m > mmax or m < mmin: val = 0.
	else: val = 1./(mmax-mmin)

	return val
	
def peak_mass(m,mu=1.34,sigma=0.02): # gaussian mass distribution

	val = gaussian(m,mu,sigma)

	return val
	
def bimod_mass(m,mu1=1.34,sigma1=0.02,mu2=1.47,sigma2=0.15,alpha=0.68,norms=False): # double Gaussian fit to Galactic double neutron star population from Farrow+ arXiv:1902.03300
	
	if norms:
		norm1 = 0.5*(scipy.special.erf(mu1/(np.sqrt(2)*sigma1))+1.)
		norm2 = 0.5*(scipy.special.erf(mu2/(np.sqrt(2)*sigma2))+1.)
	else:
		norm1 = 1.
		norm2 = 1.
	val = alpha*gaussian(m,mu1,sigma1)/norm1 + (1.-alpha)*gaussian(m,mu2,sigma2)/norm2

	return val
	
def peakcut_mass(m,mu=1.34,sigma=0.02,mmin=1.,mmax=3.): # gaussian mass distribution with high- and low-mass cutoffs

	if m > mmax or m < mmin: val = 0.
	else:
		norm = 0.5*(scipy.special.erf((mmax-mu)/(np.sqrt(2)*sigma))-scipy.special.erf((mmin-mu)/(np.sqrt(2)*sigma)))
		val = gaussian(m,mu,sigma)/norm
	
	return val
	
def bimodcut_mass(m,mu1=1.34,sigma1=0.07,mu2=1.80,sigma2=0.21,alpha=0.65,mmin=0.9,mmax=2.12,norms=True): # double gaussian mass distribution with high- and low-mass cutoffs from Alsing+

	if m > mmax or m < mmin: val = 0.
	else:
		if norms:
			norm1 = 0.5*(scipy.special.erf((mmax-mu1)/(np.sqrt(2)*sigma1))-scipy.special.erf((mmin-mu1)/(np.sqrt(2)*sigma1)))
			norm2 = 0.5*(scipy.special.erf((mmax-mu2)/(np.sqrt(2)*sigma2))-scipy.special.erf((mmin-mu2)/(np.sqrt(2)*sigma2)))
		else:
			norm1 = 1.
			norm2 = 1.
		val = alpha*gaussian(m,mu1,sigma1)/norm1 + (1.-alpha)*gaussian(m,mu2,sigma2)/norm2
	
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
	
def bimod_m1m2(m1,m2,mu1=1.34,sigma1=0.02,mu2=1.47,sigma2=0.15,alpha=0.68,norms=False): # double gaussian distribution in source frame masses, subject to m1 >= m2 convention

	if m1 < m2: val = 0.
	else:
		if norms:
			norm1 = 0.5*(scipy.special.erf(mu1/(np.sqrt(2)*sigma1))+1.)
			norm2 = 0.5*(scipy.special.erf(mu2/(np.sqrt(2)*sigma2))+1.)
		else:
			norm1 = 1.
			norm2 = 1.
		val = (alpha*gaussian(m1,mu1,sigma1)/norm1 + (1.-alpha)*gaussian(m1,mu2,sigma2)/norm2)*(alpha*gaussian(m2,mu1,sigma1)/norm1 + (1.-alpha)*gaussian(m2,mu2,sigma2)/norm2)
	
	return val
	
def peakcut_m1m2(m1,m2,mu=1.34,sigma=0.02,mmin=1.,mmax=3.): # gaussian distribution in source frame masses with high- and low-mass cutoffs, subject to m1 >= m2 convention

	if m1 < m2 or m1 > mmax or m2 < mmin: val = 0.
	else:
		norm = 0.5*(scipy.special.erf((mmax-mu)/(np.sqrt(2.)*sigma))-scipy.special.erf((mmin-mu)/(np.sqrt(2.)*sigma)))
		val = gaussian(m1,mu,sigma)*gaussian(m2,mu,sigma)/norm**2
	
	return val
	
def bimodcut_m1m2(m1,m2,mu1=1.34,sigma1=0.07,mu2=1.80,sigma2=0.21,alpha=0.65,mmin=0.9,mmax=2.12,norms=True): # double gaussian mass distribution with high- and low-mass cutoffs from Alsing+

	if m1 < m2 or m1 > mmax or m2 < mmin: val = 0.
	else:
		if norms:
			norm1 = 0.5*(scipy.special.erf((mmax-mu1)/(np.sqrt(2)*sigma1))-scipy.special.erf((mmin-mu1)/(np.sqrt(2)*sigma1)))
			norm2 = 0.5*(scipy.special.erf((mmax-mu2)/(np.sqrt(2)*sigma2))-scipy.special.erf((mmin-mu2)/(np.sqrt(2)*sigma2)))
		else:
			norm1 = 1.
			norm2 = 1.
		val = (alpha*gaussian(m1,mu1,sigma1)/norm1 + (1.-alpha)*gaussian(m1,mu2,sigma2)/norm2)*(alpha*gaussian(m2,mu1,sigma1)/norm1 + (1.-alpha)*gaussian(m2,mu2,sigma2)/norm2)
	
	return val
	
def peakcut_m1_unif_m2(m1,m2,mu=1.34,sigma=0.02,mmin=1.,mmax=3.): # gaussian distribution in source frame masses with high- and low-mass cutoffs, subject to m1 >= m2 convention

	if m1 < m2 or m1 > mmax or m2 < mmin: val = 0.
	else:
		norm = 0.5*(scipy.special.erf((mmax-mu)/(np.sqrt(2)*sigma))-scipy.special.erf((mmin-mu)/(np.sqrt(2)*sigma)))
		val = gaussian(m1,mu,sigma)*unif_mass(m2,mmin,mmax)/norm
	
	return val
	
def bimodcut_m1_unif_m2(m1,m2,mu1=1.34,sigma1=0.02,mu2=1.47,sigma2=0.15,alpha=0.68,mmin=1.,mmax=3.,norms=False): # double gaussian mass distribution with high- and low-mass cutoffs

	if m1 < m2 or m1 > mmax or m2 < mmin: val = 0.
	else:
		if norms:
			norm1 = 0.5*(scipy.special.erf(mu1/(np.sqrt(2)*sigma1))+1.)
			norm2 = 0.5*(scipy.special.erf(mu2/(np.sqrt(2)*sigma2))+1.)
		else:
			norm1 = 1.
			norm2 = 1.
		val = (alpha*gaussian(m1,mu1,sigma1)/norm1 + (1.-alpha)*gaussian(m1,mu2,sigma2)/norm2)*unif_mass(m2,mmin,mmax)
	
	return val
	
def unif_m1m2_qpair(m1,m2,mmin=5.,mmax=2e2,beta=0.): # uniform distribution in source frame masses, subject to m1 >= m2 convention and q-dependent pairing

	if m1 < m2 or m1 > mmax or m2 < mmin: val = 0.
	else:
		norm = (mmax*mmin)**beta*(mmin**(beta+1)-mmax**(beta+1))*(mmax*mmin**beta-mmin*mmax**beta)/(beta**2-1.)
		val = (m2/m1)**beta/norm
    
	return val
	
def peakcut_m1m2_qpair(m1,m2,mu=1.34,sigma=0.02,mmin=1.,mmax=3.,beta=0.): # gaussian distribution in source frame masses, subject to m1 >= m2 convention and q-dependent pairing

	if m1 < m2 or m1 > mmax or m2 < mmin: val = 0.
	else:
		norm = scipy.integrate.dblquad(lambda m1,m2 : gaussian(m1,mu,sigma)*gaussian(m2,mu,sigma)*(m2/m1)**beta, mmin, mmax, mmin, mmax)[0]
		val = gaussian(m1,mu,sigma)*gaussian(m2,mu,sigma)*(m2/m1)**beta/norm
	
	return val
	
def bimodcut_m1m2_qpair(m1,m2,mu1=1.34,sigma1=0.07,mu2=1.80,sigma2=0.21,alpha=0.65,mmin=0.9,mmax=2.12,beta=0.,norms=True): # double gaussian Alsing+ distribution in source frame masses, subject to m1 >= m2 convention and q-dependent pairing

	if m1 < m2 or m1 > mmax or m2 < mmin: val = 0.
	else:
		if norms:
			norm1 = 0.5*(scipy.special.erf((mmax-mu1)/(np.sqrt(2)*sigma1))-scipy.special.erf((mmin-mu1)/(np.sqrt(2)*sigma1)))
			norm2 = 0.5*(scipy.special.erf((mmax-mu2)/(np.sqrt(2)*sigma2))-scipy.special.erf((mmin-mu2)/(np.sqrt(2)*sigma2)))
		else:
			norm1 = 1.
			norm2 = 1.
			
		norm = scipy.integrate.dblquad(lambda m1,m2 : (alpha*gaussian(m1,mu1,sigma1)/norm1 + (1.-alpha)*gaussian(m1,mu2,sigma2)/norm2)*(alpha*gaussian(m2,mu1,sigma1)/norm1 + (1.-alpha)*gaussian(m2,mu2,sigma2)/norm2)*(m2/m1)**beta, mmin, mmax, mmin, mmax)[0]
			
		val = (alpha*gaussian(m1,mu1,sigma1)/norm1 + (1.-alpha)*gaussian(m1,mu2,sigma2)/norm2)*(alpha*gaussian(m2,mu1,sigma1)/norm1 + (1.-alpha)*gaussian(m2,mu2,sigma2)/norm2)*(m2/m1)**beta/norm
	
	return val
	
def o3a_powerpeak_m1m2_qpair(m1,m2,lpeak=0.10,alpha=2.62,beta=1.26,mmin=4.53,mmax=86.73,delta=4.88,mu=33.49,sigma=5.09): # power law + peak mass model from O3a populations paper
    
	norm = 0.5*(scipy.special.erf((mmax-mu)/(np.sqrt(2)*sigma))-scipy.special.erf((mmin-mu)/(np.sqrt(2)*sigma)))
	val = ((1.-lpeak)*powerlaw(m1,-alpha,mmin,mmax)+lpeak*gaussian(m1,mu,sigma)/norm)*smooth(m1,mmin,delta)*(m2/m1)**beta*smooth(m2,mmin,delta)/m1
    
	return val # FIXME: normalize me!

def o3a_powerbreak_m1m2_qpair(m1,m2,alpha1=1.58,alpha2=5.59,beta=1.40,mmin=4.83,mmax=3.96,delta=87.14,b=0.43): # broken power law mass model from O3a populations paper
    
    mbreak = mmin+b*(mmax-mmin)
    
    if m1 >= mmin and m1 < mbreak: val = m1**(-alpha1)*smooth(m1,mmin,delta)*(m2/m1)**beta*s(m2,mmin,delta)/m1
    elif m1 >= mbreak and m1 <= mmax: val = m1**(-alpha2)*smooth(m1,mmin,delta)*(m2/m1)**beta*s(m2,mmin,delta)/m1
    else: val = 0.
    
    return val # FIXME: normalize me!

def o3a_powerpeak_m1_unif_m2(m1,m2,mu=1.34,sigma=0.02,mmin=1.,mmax=3.):

	if m1 < m2 or m2 > mmax or m2 < mmin: val = 0.
	else:
		norm = 0.5*(scipy.special.erf((mmax-mu)/(np.sqrt(2)*sigma))-scipy.special.erf((mmin-mu)/(np.sqrt(2)*sigma)))
		val = ((1.-lpeak)*powerlaw(m1,-alpha,mmin,mmax)+lpeak*gaussian(m1,mu,sigma)/norm)*smooth(m1,mmin,delta)*unif_mass(m2,mmin,mmax)
	
	return val
	
# NSBH MASS DISTRIBUTIONS

def unif_m1_unif_m2(m1,m2,mmin=1.,mmax=3.,mmin_bh=3.,mmax_bh=23.): # uniform distribution in source frame masses, subject to m1 >= m2 convention

	if m1 < m2 or m2 > mmax or m2 < mmin or m1 < mmin_bh or m1 > mmax_bh: val = 0.
	else: val = 1./((mmax-mmin)*(mmax_bh-mmin_bh))
    
	return val
	
def unif_m1_unif_m2_qpair(m1,m2,mmin=1.,mmax=3.,beta=2.,mmin_bh=3.,mmax_bh=23.): # uniform distribution in source frame masses, subject to m1 >= m2 convention

	if m1 < m2 or m2 > mmax or m2 < mmin or m1 < mmin_bh or m1 > mmax_bh: val = 0.
	else:
		norm = -(mmax_bh*mmin_bh)**(-beta)*(mmax**(beta+1)-mmin**(beta+1))*(mmax_bh*mmin_bh**beta-mmin_bh*mmax_bh**beta)/(beta**2-1.)
		val = (m2/m1)**beta/norm
    
	return val
	
# FIXME: NEED TO FILL IN MODELS BELOW
	
def unif_m1_peakcut_m2(m1,m2,mmin=1.,mmax=3.,mmin_bh=3.,mmax_bh=23.): # uniform distribution in source frame masses, subject to m1 >= m2 convention

	if m1 < m2 or m2 > mmax or m2 < mmin or m1 < mmin_bh or m1 > mmax_bh: val = 0.
	else: val = 1./((mmax-mmin)*(mmax_bh-mmin_bh))
    
	return val
	
def unif_m1_peakcut_m2_qpair(m1,m2,mmin=1.,mmax=3.,beta=2.,mmin_bh=3.,mmax_bh=23.): # uniform distribution in source frame masses, subject to m1 >= m2 convention

	if m1 < m2 or m2 > mmax or m2 < mmin or m1 < mmin_bh or m1 > mmax_bh: val = 0.
	else:
		norm = -(mmax_bh*mmin_bh)**(-beta)*(mmax**(beta+1)-mmin**(beta+1))*(mmax_bh*mmin_bh**beta-mmin_bh*mmax_bh**beta)/(beta**2-1.)
		val = (m2/m1)**beta/norm
    
	return val
	
def unif_m1_bimodcut_m2(m1,m2,mmin=1.,mmax=3.,mmin_bh=3.,mmax_bh=23.): # uniform distribution in source frame masses, subject to m1 >= m2 convention

	if m1 < m2 or m2 > mmax or m2 < mmin or m1 < mmin_bh or m1 > mmax_bh: val = 0.
	else: val = 1./((mmax-mmin)*(mmax_bh-mmin_bh))
    
	return val
	
def unif_m1_bimodcut_m2_qpair(m1,m2,mmin=1.,mmax=3.,beta=2.,mmin_bh=3.,mmax_bh=23.): # uniform distribution in source frame masses, subject to m1 >= m2 convention

	if m1 < m2 or m2 > mmax or m2 < mmin or m1 < mmin_bh or m1 > mmax_bh: val = 0.
	else:
		norm = -(mmax_bh*mmin_bh)**(-beta)*(mmax**(beta+1)-mmin**(beta+1))*(mmax_bh*mmin_bh**beta-mmin_bh*mmax_bh**beta)/(beta**2-1.)
		val = (m2/m1)**beta/norm
    
	return val
	
# FIXME: NEED TO FILL IN MODELS ABOVE
	
# LOOKUP FUNCTIONS

pop_priors = {'unif_mass': unif_mass, 'peak_mass': peak_mass, 'bimod_mass': bimod_mass, 'peakcut_mass': peakcut_mass, 'bimodcut_mass': bimodcut_mass, 'unif_m1m2': unif_m1m2, 'peak_m1m2': peak_m1m2, 'bimod_m1m2': bimod_m1m2, 'peakcut_m1m2': peakcut_m1m2, 'bimodcut_m1m2': bimodcut_m1m2, 'peakcut_m1_unif_m2': peakcut_m1_unif_m2, 'bimodcut_m1_unif_m2': bimodcut_m1_unif_m2, 'unif_m1m2_qpair': unif_m1m2_qpair, 'peakcut_m1m2_qpair': peakcut_m1m2_qpair, 'bimodcut_m1m2_qpair': bimodcut_m1m2_qpair, 'o3a_powerpeak_m1m2_qpair': o3a_powerpeak_m1m2_qpair, 'o3a_powerbreak_m1m2_qpair': o3a_powerbreak_m1m2_qpair, 'o3a_powerpeak_m1_unif_m2': o3a_powerpeak_m1_unif_m2, 'unif_m1_unif_m2': unif_m1_unif_m2}

pop_params = {'unif_mass': 'mmin,mmax', 'peak_mass': 'mu,sigma', 'bimod_mass': 'mu1,sigma1,mu2,sigma2,alpha', 'peakcut_mass': 'mu,sigma,mmin,mmax', 'bimodcut_mass': 'mu1,sigma1,mu2,sigma2,alpha,mmin,mmax', 'unif_m1m2': 'mmin,mmax', 'peak_m1m2': 'mu,sigma', 'bimod_m1m2': 'mu1,sigma1,mu2,sigma2,alpha', 'peakcut_m1m2': 'mu,sigma,mmin,mmax', 'bimodcut_m1m2': 'mu1,sigma1,mu2,sigma2,alpha,mmin,mmax', 'peakcut_m1_unif_m2': 'mu,sigma,mmin,mmax', 'bimodcut_m1_unif_m2': 'mu1,sigma1,mu2,sigma2,alpha,mmin,mmax', 'unif_m1m2_qpair': 'mmin,mmax,beta', 'peakcut_m1m2_qpair': 'mu,sigma,mmin,mmax,beta', 'bimodcut_m1m2_qpair': 'mu1,sigma1,mu2,sigma2,alpha,mmin,mmax,beta', 'o3a_powerpeak_m1m2_qpair': 'lpeak,alpha,beta,mmin,mmax,delta,mu,sigma', 'o3a_powerbreak_m1m2_qpair': 'alpha1,alpha2,beta,mmin,mmax,delta,b', 'o3a_powerpeak_m1_unif_m2': 'm1,m2,mu,sigma,mmin,mmax', 'unif_m1_unif_m2': 'm1,m2,mmin,mmax'}

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
