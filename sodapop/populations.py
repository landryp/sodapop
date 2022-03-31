#!/usr/bin/env python

import numpy as np
import scipy.special

### DEFINE CONSTANTS AND USEFUL RELATIONS

c = 2.998e10
G = 6.674e-8
Msun = 1.988e33

def ILove(m,Lambda): # empirical I-Love relation from arXiv:1608.02582

	Ibar = np.exp(1.496+0.05951*np.log(Lambda)+0.02238*np.log(Lambda)**2-6.953e-4*np.log(Lambda)**3+8.345e-6*np.log(Lambda)**4) 

	return Ibar*G**2*m**3/c**4
	
def fmaxLove(m,Lambda): # empirical fit to Kepler frequency in Hz from arXiv:0901.1268

	fmaxkHz = 1.08e3*np.sqrt(m)*(G*m*Msun*c**(-2)*(Lambda/0.0093)**(1./6.)/1e4)**(-1.5)
	
	return fmaxkHz*1e3
	
def chiI(m,I,f): # definition of chi in terms of rotational frequency and moment of inertia

	return 2.*np.pi*f*I*c/(G*m*Msun**2)

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

### INDIVIDUAL SPIN DISTRIBUTIONS

def unif_chi(chi,lambdaa): # uniform spin magnitude distribution

	p = uniform(chi,lambdaa)

	return p
	
### BINARY SPIN DISTRIBUTIONS

def unif_chi1_unif_chi2(chi1,chi2,lambdaa): # uniform distribution in component spin magnitudes

	if np.isscalar(chi1): chi1 = np.array([chi1])
	else: chi1 = np.array(chi1)
	if np.isscalar(chi2): chi2 = np.array([chi2])
	else: chi2 = np.array(chi2)
	
	p = unif_chi(m1,lambdaa)*unif_chi(m2,lambdaa)
	
	return p
	
### BINARY MASS & SPIN DISTRIBUTIONS

def unif_m1m2_common_chi1chi2(m1,m2,chi1,chi2,lambdaa): # uniform distribution in component masses and spins, subject to common-EOS dependent maximum spin

	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	if np.isscalar(chi1): chi1 = np.array([chi1])
	else: chi1 = np.array(chi1)
	if np.isscalar(chi2): chi2 = np.array([chi2])
	else: chi2 = np.array(chi2)
	z = np.zeros(len(m1))
	
	I1 = ILove(m1,Lambda1)
	I2 = ILove(m2,Lambda2)
	f1max = fmaxLove(m1,Lambda1)
	f2max = fmaxLove(m2,Lambda2)
	chi1max = chiI(m1,I1,f1max)
	chi2max = chiI(m2,I2,f2max)
	
	common1 = (lambdaa[2],min(lambdaa[3],chi1max))
	common2 = (lambdaa[2],min(lambdaa[3],chi2max))
	
	p = unif_m1m2(m1,m2,lambdaa)*unif_chi(chi1,common1)*unif_chi(chi2,common2)
	
	return np.where(m1 < m2, z, p)
	
def peakcut_m1m2_common_chi1chi2(m1,m2,chi1,chi2,lambdaa): # uniform distribution in component masses and spins, subject to common-EOS dependent maximum spin

	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	if np.isscalar(chi1): chi1 = np.array([chi1])
	else: chi1 = np.array(chi1)
	if np.isscalar(chi2): chi2 = np.array([chi2])
	else: chi2 = np.array(chi2)
	z = np.zeros(len(m1))
	
	I1 = ILove(m1,Lambda1)
	I2 = ILove(m2,Lambda2)
	f1max = fmaxLove(m1,Lambda1)
	f2max = fmaxLove(m2,Lambda2)
	chi1max = chiI(m1,I1,f1max)
	chi2max = chiI(m2,I2,f2max)
	
	common1 = (lambdaa[4],min(lambdaa[5],chi1max))
	common2 = (lambdaa[4],min(lambdaa[5],chi2max))
	
	p = peakcut_m1m2(m1,m2,lambdaa)*unif_chi(chi1,common1)*unif_chi(chi2,common2)
	
	return np.where(m1 < m2, z, p)
	
def bimodcut_m1m2_common_chi1chi2(m1,m2,chi1,chi2,lambdaa): # uniform distribution in component masses and spins, subject to common-EOS dependent maximum spin

	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	if np.isscalar(chi1): chi1 = np.array([chi1])
	else: chi1 = np.array(chi1)
	if np.isscalar(chi2): chi2 = np.array([chi2])
	else: chi2 = np.array(chi2)
	z = np.zeros(len(m1))
	
	I1 = ILove(m1,Lambda1)
	I2 = ILove(m2,Lambda2)
	f1max = fmaxLove(m1,Lambda1)
	f2max = fmaxLove(m2,Lambda2)
	chi1max = chiI(m1,I1,f1max)
	chi2max = chiI(m2,I2,f2max)
	
	common1 = (lambdaa[7],min(lambdaa[8],chi1max))
	common2 = (lambdaa[7],min(lambdaa[8],chi2max))
	
	p = bimodcut_m1m2(m1,m2,lambdaa)*unif_chi(chi1,common1)*unif_chi(chi2,common2)
	
	return np.where(m1 < m2, z, p)

### INDIVIDUAL TIDAL DEFORMABILITY DISTRIBUTIONS

def unif_Lambda(Lambda,lambdaa): # uniform tidal deformability distribution

	p = uniform(Lambda,lambdaa)

	return p
	
def peak_Lambda(Lambda,lambdaa): # gaussian tidal deformability distribution

	p = gaussian(Lambda,lambdaa)

	return p
	
### BINARY TIDAL DEFORMABILITY DISTRIBUTIONS

def unif_Lambda1_unif_Lambda2(Lambda1,Lambda2,lambdaa): # uniform distribution in component tidal deformabilities

	if np.isscalar(Lambda1): Lambda1 = np.array([Lambda1])
	else: Lambda1 = np.array(Lambda1)
	if np.isscalar(Lambda2): Lambda2 = np.array([Lambda2])
	else: Lambda2 = np.array(Lambda2)
	
	p = unif_Lambda(Lambda1,lambdaa)*unif_Lambda(Lambda2,lambdaa)
	
	return p
	
def unif_Lambda1Lambda2(Lambda1,Lambda2,lambdaa): # uniform distribution in component tidal deformabilities

	if np.isscalar(Lambda1): Lambda1 = np.array([Lambda1])
	else: Lambda1 = np.array(Lambda1)
	if np.isscalar(Lambda2): Lambda2 = np.array([Lambda2])
	else: Lambda2 = np.array(Lambda2)
	z = np.zeros(len(Lambda1))
	
	p = unif_Lambda(Lambda1,lambdaa)*unif_Lambda(Lambda2,lambdaa)
	
	return np.where(Lambda2 < Lambda1, z, p)

### BINARY MASS & TIDAL DEFORMABILITY DISTRIBUTIONS

def unif_m1m2_common_Lambda1Lambda2(m1,m2,Lambda1,Lambda2,lambdaa): # uniform distribution in component masses and tidal deformabilities, subject to common-EOS approximation

	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	if np.isscalar(Lambda1): Lambda1 = np.array([Lambda1])
	else: Lambda1 = np.array(Lambda1)
	if np.isscalar(Lambda2): Lambda2 = np.array([Lambda2])
	else: Lambda2 = np.array(Lambda2)
	z = np.zeros(len(m1))
	
	common_Lambda1 = Lambda2*(m2/m1)**6
	common_std = 10.
	common = (common_Lambda1,common_std)
	
	p = unif_m1m2(m1,m2,lambdaa)*peak_Lambda(Lambda1,common)*unif_Lambda(Lambda2,lambdaa[2:])
	
	return np.where((m1 < m2) | (Lambda2 < Lambda1), z, p)

def peakcut_m1m2_common_Lambda1Lambda2(m1,m2,Lambda1,Lambda2,lambdaa): # uniform distribution in component masses and tidal deformabilities, subject to common-EOS approximation

	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	if np.isscalar(Lambda1): Lambda1 = np.array([Lambda1])
	else: Lambda1 = np.array(Lambda1)
	if np.isscalar(Lambda2): Lambda2 = np.array([Lambda2])
	else: Lambda2 = np.array(Lambda2)
	z = np.zeros(len(m1))
	
	common_Lambda1 = Lambda2*(m2/m1)**6
	common_std = 10.
	common = (common_Lambda1,common_std)
	
	p = peakcut_m1m2(m1,m2,lambdaa)*peak_Lambda(Lambda1,common)*unif_Lambda(Lambda2,lambdaa[4:])
	
	return np.where((m1 < m2) | (Lambda2 < Lambda1), z, p)
	
def bimodcut_m1m2_common_Lambda1Lambda2(m1,m2,Lambda1,Lambda2,lambdaa): # uniform distribution in component masses and tidal deformabilities, subject to common-EOS approximation

	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	if np.isscalar(Lambda1): Lambda1 = np.array([Lambda1])
	else: Lambda1 = np.array(Lambda1)
	if np.isscalar(Lambda2): Lambda2 = np.array([Lambda2])
	else: Lambda2 = np.array(Lambda2)
	z = np.zeros(len(m1))
	
	common_Lambda1 = Lambda2*(m2/m1)**6
	common_std = 10.
	common = (common_Lambda1,common_std)
	
	p = bimodcut_m1m2(m1,m2,lambdaa)*peak_Lambda(Lambda1,common)*unif_Lambda(Lambda2,lambdaa[7:])
	
	return np.where((m1 < m2) | (Lambda2 < Lambda1), z, p)
	
### BINARY MASS, SPIN & TIDAL DEFORMABILITY DISTRIBUTIONS

def unif_m1m2_unif_chi1_unif_chi2_unif_Lambda1Lambda2(m1,m2,chi1,chi2,Lambda1,Lambda2,lambdaa): # bimodal mass model with random pairing, uniform uncorrelated spins and common-EOS approximated tidal deformabilities

	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	if np.isscalar(chi1): chi1 = np.array([chi1])
	else: chi1 = np.array(chi1)
	if np.isscalar(chi2): chi2 = np.array([chi2])
	else: chi2 = np.array(chi2)
	if np.isscalar(Lambda1): Lambda1 = np.array([Lambda1])
	else: Lambda1 = np.array(Lambda1)
	if np.isscalar(Lambda2): Lambda2 = np.array([Lambda2])
	else: Lambda2 = np.array(Lambda2)
	z = np.zeros(len(m1))
	
	p = unif_m1m2(m1,m2,lambdaa)*unif_chi1_unif_chi2(chi1,chi2,lambdaa[2:4])*unif_Lambda1Lambda2(Lambda1,Lambda2,lambdaa[4:])

	return np.where((m1 < m2) | (Lambda2 < Lambda1), z, p)
	
def peakcut_m1m2_unif_chi1_unif_chi2_unif_Lambda1Lambda2(m1,m2,chi1,chi2,Lambda1,Lambda2,lambdaa): # bimodal mass model with random pairing, uniform uncorrelated spins and common-EOS approximated tidal deformabilities

	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	if np.isscalar(chi1): chi1 = np.array([chi1])
	else: chi1 = np.array(chi1)
	if np.isscalar(chi2): chi2 = np.array([chi2])
	else: chi2 = np.array(chi2)
	if np.isscalar(Lambda1): Lambda1 = np.array([Lambda1])
	else: Lambda1 = np.array(Lambda1)
	if np.isscalar(Lambda2): Lambda2 = np.array([Lambda2])
	else: Lambda2 = np.array(Lambda2)
	z = np.zeros(len(m1))
	
	p = peakcut_m1m2(m1,m2,lambdaa)*unif_chi1_unif_chi2(chi1,chi2,lambdaa[4:6])*unif_Lambda1Lambda2(Lambda1,Lambda2,lambdaa[6:])

	return np.where((m1 < m2) | (Lambda2 < Lambda1), z, p)

def bimodcut_m1m2_unif_chi1_unif_chi2_unif_Lambda1Lambda2(m1,m2,chi1,chi2,Lambda1,Lambda2,lambdaa): # bimodal mass model with random pairing, uniform uncorrelated spins and common-EOS approximated tidal deformabilities

	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	if np.isscalar(chi1): chi1 = np.array([chi1])
	else: chi1 = np.array(chi1)
	if np.isscalar(chi2): chi2 = np.array([chi2])
	else: chi2 = np.array(chi2)
	if np.isscalar(Lambda1): Lambda1 = np.array([Lambda1])
	else: Lambda1 = np.array(Lambda1)
	if np.isscalar(Lambda2): Lambda2 = np.array([Lambda2])
	else: Lambda2 = np.array(Lambda2)
	z = np.zeros(len(m1))
	
	p = bimodcut_m1m2(m1,m2,lambdaa)*unif_chi1_unif_chi2(chi1,chi2,lambdaa[7:9])*unif_Lambda1Lambda2(Lambda1,Lambda2,lambdaa[9:])

	return np.where((m1 < m2) | (Lambda2 < Lambda1), z, p)
	
def unif_m1m2_unif_chi1_unif_chi2_common_Lambda1Lambda2(m1,m2,chi1,chi2,Lambda1,Lambda2,lambdaa): # bimodal mass model with random pairing, uniform uncorrelated spins and common-EOS approximated tidal deformabilities

	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	if np.isscalar(chi1): chi1 = np.array([chi1])
	else: chi1 = np.array(chi1)
	if np.isscalar(chi2): chi2 = np.array([chi2])
	else: chi2 = np.array(chi2)
	if np.isscalar(Lambda1): Lambda1 = np.array([Lambda1])
	else: Lambda1 = np.array(Lambda1)
	if np.isscalar(Lambda2): Lambda2 = np.array([Lambda2])
	else: Lambda2 = np.array(Lambda2)
	z = np.zeros(len(m1))
	
	p = unif_m1m2_common_Lambda1Lambda2(m1,m2,Lambda1,Lambda2,lambdaa[:2]+lambdaa[4:])*unif_chi1_unif_chi2(chi1,chi2,lambdaa[2:4])

	return np.where((m1 < m2) | (Lambda2 < Lambda1), z, p)
	
def peakcut_m1m2_unif_chi1_unif_chi2_common_Lambda1Lambda2(m1,m2,chi1,chi2,Lambda1,Lambda2,lambdaa): # bimodal mass model with random pairing, uniform uncorrelated spins and common-EOS approximated tidal deformabilities

	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	if np.isscalar(chi1): chi1 = np.array([chi1])
	else: chi1 = np.array(chi1)
	if np.isscalar(chi2): chi2 = np.array([chi2])
	else: chi2 = np.array(chi2)
	if np.isscalar(Lambda1): Lambda1 = np.array([Lambda1])
	else: Lambda1 = np.array(Lambda1)
	if np.isscalar(Lambda2): Lambda2 = np.array([Lambda2])
	else: Lambda2 = np.array(Lambda2)
	z = np.zeros(len(m1))
	
	p = peakcut_m1m2_common_Lambda1Lambda2(m1,m2,Lambda1,Lambda2,lambdaa[:4]+lambdaa[6:])*unif_chi1_unif_chi2(chi1,chi2,lambdaa[4:6])

	return np.where((m1 < m2) | (Lambda2 < Lambda1), z, p)

def bimodcut_m1m2_unif_chi1_unif_chi2_common_Lambda1Lambda2(m1,m2,chi1,chi2,Lambda1,Lambda2,lambdaa): # bimodal mass model with random pairing, uniform uncorrelated spins and common-EOS approximated tidal deformabilities

	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	if np.isscalar(chi1): chi1 = np.array([chi1])
	else: chi1 = np.array(chi1)
	if np.isscalar(chi2): chi2 = np.array([chi2])
	else: chi2 = np.array(chi2)
	if np.isscalar(Lambda1): Lambda1 = np.array([Lambda1])
	else: Lambda1 = np.array(Lambda1)
	if np.isscalar(Lambda2): Lambda2 = np.array([Lambda2])
	else: Lambda2 = np.array(Lambda2)
	z = np.zeros(len(m1))
	
	p = bimodcut_m1m2_common_Lambda1Lambda2(m1,m2,Lambda1,Lambda2,lambdaa[:7]+lambdaa[9:])*unif_chi1_unif_chi2(chi1,chi2,lambdaa[7:9])

	return np.where((m1 < m2) | (Lambda2 < Lambda1), z, p)
	
def unif_m1m2_common_chi1chi2_common_Lambda1Lambda2(m1,m2,chi1,chi2,Lambda1,Lambda2,lambdaa): # uniform mass model with random pairing, uniform uncorrelated spins and common-EOS approximated tidal deformabilities

	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	if np.isscalar(chi1): chi1 = np.array([chi1])
	else: chi1 = np.array(chi1)
	if np.isscalar(chi2): chi2 = np.array([chi2])
	else: chi2 = np.array(chi2)
	if np.isscalar(Lambda1): Lambda1 = np.array([Lambda1])
	else: Lambda1 = np.array(Lambda1)
	if np.isscalar(Lambda2): Lambda2 = np.array([Lambda2])
	else: Lambda2 = np.array(Lambda2)
	z = np.zeros(len(m1))
	
	I1 = ILove(m1,Lambda1)
	I2 = ILove(m2,Lambda2)
	f1max = fmaxLove(m1,Lambda1)
	f2max = fmaxLove(m2,Lambda2)
	chi1max = chiI(m1,I1,f1max)
	chi2max = chiI(m2,I2,f2max)
	
	common1 = (lambdaa[2],min(lambdaa[3],chi1max))
	common2 = (lambdaa[2],min(lambdaa[3],chi2max))
	
	p = unif_m1m2_common_Lambda1Lambda2(m1,m2,Lambda1,Lambda2,lambdaa[:2]+lambdaa[4:])*unif_chi(chi1,common1)*unif_chi(chi2,common2)

	return np.where((m1 < m2) | (Lambda2 < Lambda1), z, p)
	
def peakcut_m1m2_common_chi1chi2_common_Lambda1Lambda2(m1,m2,chi1,chi2,Lambda1,Lambda2,lambdaa): # Gaussian mass model with random pairing, uniform uncorrelated spins and common-EOS approximated tidal deformabilities

	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	if np.isscalar(chi1): chi1 = np.array([chi1])
	else: chi1 = np.array(chi1)
	if np.isscalar(chi2): chi2 = np.array([chi2])
	else: chi2 = np.array(chi2)
	if np.isscalar(Lambda1): Lambda1 = np.array([Lambda1])
	else: Lambda1 = np.array(Lambda1)
	if np.isscalar(Lambda2): Lambda2 = np.array([Lambda2])
	else: Lambda2 = np.array(Lambda2)
	z = np.zeros(len(m1))
	
	I1 = ILove(m1,Lambda1)
	I2 = ILove(m2,Lambda2)
	f1max = fmaxLove(m1,Lambda1)
	f2max = fmaxLove(m2,Lambda2)
	chi1max = chiI(m1,I1,f1max)
	chi2max = chiI(m2,I2,f2max)
	
	common1 = (lambdaa[4],min(lambdaa[5],chi1max))
	common2 = (lambdaa[4],min(lambdaa[5],chi2max))
	
	p = peakcut_m1m2_common_Lambda1Lambda2(m1,m2,Lambda1,Lambda2,lambdaa[:4]+lambdaa[6:])*unif_chi(chi1,common1)*unif_chi(chi2,common2)

	return np.where((m1 < m2) | (Lambda2 < Lambda1), z, p)

def bimodcut_m1m2_common_chi1chi2_common_Lambda1Lambda2(m1,m2,chi1,chi2,Lambda1,Lambda2,lambdaa): # bimodal mass model with random pairing, uniform uncorrelated spins and common-EOS approximated tidal deformabilities

	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	if np.isscalar(chi1): chi1 = np.array([chi1])
	else: chi1 = np.array(chi1)
	if np.isscalar(chi2): chi2 = np.array([chi2])
	else: chi2 = np.array(chi2)
	if np.isscalar(Lambda1): Lambda1 = np.array([Lambda1])
	else: Lambda1 = np.array(Lambda1)
	if np.isscalar(Lambda2): Lambda2 = np.array([Lambda2])
	else: Lambda2 = np.array(Lambda2)
	z = np.zeros(len(m1))
	
	I1 = ILove(m1,Lambda1)
	I2 = ILove(m2,Lambda2)
	f1max = fmaxLove(m1,Lambda1)
	f2max = fmaxLove(m2,Lambda2)
	chi1max = chiI(m1,I1,f1max)
	chi2max = chiI(m2,I2,f2max)
	
	common1 = (lambdaa[7],min(lambdaa[8],chi1max))
	common2 = (lambdaa[7],min(lambdaa[8],chi2max))
	
	p = bimodcut_m1m2_common_Lambda1Lambda2(m1,m2,Lambda1,Lambda2,lambdaa[:7]+lambdaa[9:])*unif_chi(chi1,common1)*unif_chi(chi2,common2)

	return np.where((m1 < m2) | (Lambda2 < Lambda1), z, p)

# LOOKUP FUNCTIONS

pop_priors = {'unif_mass': unif_mass, 'peak_mass': peak_mass, 'bimod_mass': bimod_mass, 'peakcut_mass': peakcut_mass, 'bimodcut_mass': bimodcut_mass, 'unif_m1m2': unif_m1m2, 'peak_m1m2': peak_m1m2, 'bimod_m1m2': bimod_m1m2, 'power_m1m2': power_m1m2, 'peakcut_m1m2': peakcut_m1m2, 'bimodcut_m1m2': bimodcut_m1m2, 'unif_m1m2_qpair': unif_m1m2_qpair, 'power_m1m2_qpair': power_m1m2_qpair, 'peakcut_m1m2_qpair': peakcut_m1m2_qpair, 'bimodcut_m1m2_qpair': bimodcut_m1m2_qpair, 'unif_m1_unif_m2': unif_m1_unif_m2, 'unif_m1_unif_m2_qpair': unif_m1_unif_m2_qpair, 'unif_m1_power_m2': unif_m1_power_m2, 'unif_m1_power_m2_qpair': unif_m1_power_m2_qpair, 'unif_m1_peakcut_m2': unif_m1_peakcut_m2, 'unif_m1_peakcut_m2_qpair': unif_m1_peakcut_m2_qpair, 'unif_m1_bimodcut_m2': unif_m1_bimodcut_m2, 'unif_m1_bimodcut_m2_qpair': unif_m1_bimodcut_m2_qpair, 'power_m1_unif_m2_qpair': power_m1_unif_m2_qpair, 'bimodcut_m1_unif_m2_qpair': bimodcut_m1_unif_m2_qpair, 'unif_chi': unif_chi, 'unif_chi1_unif_chi2': unif_chi1_unif_chi2, 'unif_m1m2_common_chi1chi2': unif_m1m2_common_chi1chi2, 'peakcut_m1m2_common_chi1chi2': peakcut_m1m2_common_chi1chi2, 'bimodcut_m1m2_common_chi1chi2': bimodcut_m1m2_common_chi1chi2, 'unif_Lambda': unif_Lambda, 'peak_Lambda': peak_Lambda, 'unif_Lambda1_unif_Lambda2': unif_Lambda1_unif_Lambda2, 'unif_Lambda1Lambda2': unif_Lambda1Lambda2, 'unif_m1m2_common_Lambda1Lambda2': unif_m1m2_common_Lambda1Lambda2, 'peakcut_m1m2_common_Lambda1Lambda2': peakcut_m1m2_common_Lambda1Lambda2, 'bimodcut_m1m2_common_Lambda1Lambda2': bimodcut_m1m2_common_Lambda1Lambda2, 'unif_m1m2_unif_chi1_unif_chi2_unif_Lambda1Lambda2': unif_m1m2_unif_chi1_unif_chi2_unif_Lambda1Lambda2, 'peakcut_m1m2_unif_chi1_unif_chi2_unif_Lambda1Lambda2': peakcut_m1m2_unif_chi1_unif_chi2_unif_Lambda1Lambda2, 'bimodcut_m1m2_unif_chi1_unif_chi2_unif_Lambda1Lambda2': bimodcut_m1m2_unif_chi1_unif_chi2_unif_Lambda1Lambda2, 'unif_m1m2_unif_chi1_unif_chi2_common_Lambda1Lambda2': unif_m1m2_unif_chi1_unif_chi2_common_Lambda1Lambda2, 'peakcut_m1m2_unif_chi1_unif_chi2_common_Lambda1Lambda2': peakcut_m1m2_unif_chi1_unif_chi2_common_Lambda1Lambda2, 'bimodcut_m1m2_unif_chi1_unif_chi2_common_Lambda1Lambda2': bimodcut_m1m2_unif_chi1_unif_chi2_common_Lambda1Lambda2, 'unif_m1m2_common_chi1chi2_common_Lambda1Lambda2': unif_m1m2_common_chi1chi2_common_Lambda1Lambda2, 'peakcut_m1m2_common_chi1chi2_common_Lambda1Lambda2': peakcut_m1m2_common_chi1chi2_common_Lambda1Lambda2, 'bimodcut_m1m2_common_chi1chi2_common_Lambda1Lambda2': bimodcut_m1m2_common_chi1chi2_common_Lambda1Lambda2}

pop_params = {'unif_mass': 'mmin,mmax', 'peak_mass': 'mu,sigma', 'bimod_mass': 'mu1,sigma1,mu2,sigma2,w', 'peakcut_mass': 'mu,sigma,mmin,mmax', 'bimodcut_mass': 'mu1,sigma1,mu2,sigma2,w,mmin,mmax', 'unif_m1m2': 'mmin,mmax', 'peak_m1m2': 'mu,sigma', 'bimod_m1m2': 'mu1,sigma1,mu2,sigma2,w', 'power_m1m2': 'alpha,mmin,mmax', 'peakcut_m1m2': 'mu,sigma,mmin,mmax', 'bimodcut_m1m2': 'mu1,sigma1,mu2,sigma2,w,mmin,mmax', 'unif_m1m2_qpair': 'mmin,mmax,beta', 'power_m1m2_qpair': 'alpha,mmin,mmax,beta', 'peakcut_m1m2_qpair': 'mu,sigma,mmin,mmax,beta', 'bimodcut_m1m2_qpair': 'mu1,sigma1,mu2,sigma2,w,mmin,mmax,beta', 'unif_m1_unif_m2': 'mmin,mmax', 'unif_m1_unif_m2_qpair': 'mmin,mmax,beta', 'unif_m1_power_m2': 'alpha,mmin,mmax', 'unif_m1_power_m2_qpair': 'alpha,mmin,mmax,beta', 'unif_m1_peakcut_m2': 'mu,sigma,mmin,mmax', 'unif_m1_peakcut_m2_qpair': 'mu,sigma,mmin,mmax,beta', 'unif_m1_bimodcut_m2': 'mu1,sigma1,mu2,sigma2,w,mmin,mmax', 'unif_m1_bimodcut_m2_qpair': 'mu1,sigma1,mu2,sigma2,w,mmin,mmax,beta', 'power_m1_unif_m2_qpair': 'alpha,mmin,mmax,beta', 'bimodcut_m1_unif_m2_qpair': 'mu1,sigma1,mu2,sigma2,w,mmin,mmax,beta', 'unif_chi': 'chimin,chimax', 'unif_chi1_unif_chi2': 'chimin,chimax', 'unif_m1m2_common_chi1chi2': 'mmin,mmax,chimin,chimax', 'peakcut_m1m2_common_chi1chi2': 'mu,sigma,mmin,mmax,chimin,chimax', 'bimodcut_m1m2_common_chi1chi2': 'mu1,sigma1,mu2,sigma2,w,mmin,mmax,chimin,chimax', 'unif_Lambda': 'Lmin,Lmax', 'peak_Lambda': 'Lmu,Lsigma', 'unif_Lambda1_unif_Lambda2': 'Lmin,Lmax', 'unif_Lambda1Lambda2': 'Lmin,Lmax', 'unif_m1m2_common_Lambda1Lambda2': 'mmin,mmax,Lmin,Lmax', 'peakcut_m1m2_common_Lambda1Lambda2': 'mu,sigma,mmin,mmax,Lmin,Lmax', 'bimodcut_m1m2_common_Lambda1Lambda2': 'mu1,sigma1,mu2,sigma2,w,mmin,mmax,Lmin,Lmax', 'unif_m1m2_unif_chi1_unif_chi2_unif_Lambda1Lambda2': 'mmin,mmax,chimin,chimax,Lmin,Lmax', 'peakcut_m1m2_unif_chi1_unif_chi2_unif_Lambda1Lambda2': 'mu,sigma,mmin,mmax,chimin,chimax,Lmin,Lmax', 'bimodcut_m1m2_unif_chi1_unif_chi2_unif_Lambda1Lambda2': 'mu1,sigma1,mu2,sigma2,w,mmin,mmax,chimin,chimax,Lmin,Lmax', 'unif_m1m2_unif_chi1_unif_chi2_common_Lambda1Lambda2': 'mmin,mmax,chimin,chimax,Lmin,Lmax', 'peakcut_m1m2_unif_chi1_unif_chi2_common_Lambda1Lambda2': 'mu,sigma,mmin,mmax,chimin,chimax,Lmin,Lmax', 'bimodcut_m1m2_unif_chi1_unif_chi2_common_Lambda1Lambda2': 'mu1,sigma1,mu2,sigma2,w,mmin,mmax,chimin,chimax,Lmin,Lmax', 'unif_m1m2_common_chi1chi2_common_Lambda1Lambda2': 'mmin,mmax,chimin,chimax,Lmin,Lmax', 'peakcut_m1m2_common_chi1chi2_common_Lambda1Lambda2': 'mu,sigma,mmin,mmax,chimin,chimax,Lmin,Lmax', 'bimodcut_m1m2_common_chi1chi2_common_Lambda1Lambda2': 'mu1,sigma1,mu2,sigma2,w,mmin,mmax,chimin,chimax,Lmin,Lmax'}

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
