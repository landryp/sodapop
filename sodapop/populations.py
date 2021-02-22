#!/usr/bin/env python

import numpy as np

# basic functions and distributions

def powerlaw(m,alpha,mmax):
    
    if m <= mmax: val = m**alpha
    else: val = 0.
        
    return val

def gaussian(m,mu,sigma):
    
    return np.exp(-((m-mu)/(np.sqrt(2)*sigma))**2)/(sigma*np.sqrt(2*np.pi))

def smooth(m,mmin,delta):
    
    if m < mmin: val =  0.
    elif m >= mmin + delta: val = 1.
    else: val = 1./(np.exp(delta/(m-mmin)+delta/(m-mmin-delta))+1.)
    
    return val

def double_gaussian(m,mu1=1.34,sigma1=0.02,mu2=1.47,sigma2=0.15,alpha=0.68): # double Gaussian fit to Galactic double neutron star population from Farrow+ arXiv:1902.03300
    
    return alpha*gaussian(m,mu1,sigma1) + (1.-alpha)*gaussian(m,mu2,sigma2)
    
def src_to_det(m,z): # convert to detector frame mass from source frame
    
    return m*(1.+z)
    
def dL_to_z(dL,H0=67.4): # convert distance to redshift, assuming a cosmology
    
    return H0*dL/2.998e5

def flat(x): # truly flat prior
    
    return 1.
    
# define available prior distributions in (m1,m2)

def flat_m1m2(m1,m2,dL,mmin=0.,mmax=1e10): # flat in source frame masses, subject to m1 >= m2 convention
    
    if m1 >= m2 and m2 >= mmin and m1 <= mmax : prior = 1.
    else: prior = 0.
    
    return prior

def flat_m1m2det(m1,m2,dL,dL_prior=flat): # flat in detector frame masses, subject to m1 >= m2 convention
    
    if m1 >= m2: prior = (1.+dL_to_z(dL))**2*dL_prior(dL)
    else: prior = 0.
    
    return prior

def flat_mceta(m1,m2,dL): # flat in chirp mass and symmetric mass ratio
    
    if m1 >= m2: prior = (m1-m2)*(m1*m2)**0.6/(m1+m2)**3.2
    else: prior = 0.
    
    return prior

def flat_mcetadet(m1,m2,dL,dL_prior=flat): # flat in chirp mass and symmetric mass ratio
    
    if m1 >= m2: prior = (1.+dL_to_z(dL))**2*dL_prior(dL)*(m1-m2)*(m1*m2)**0.6/(m1+m2)**3.2
    else: prior = 0.
    
    return prior
    
def doublegaussian_m1m2(m1,m2,dL,mu1=1.34,sigma1=0.02,mu2=1.47,sigma2=0.15,w=0.68): # double Gaussian for m1, m2, subject to m1 >= m2 convention
    
    if m1 >= m2: prior = double_gaussian(m1,mu1,sigma1,mu2,sigma2,w)*double_gaussian(m2,mu1,sigma1,mu2,sigma2,w)
    else: prior = 0.
    
    return prior
    
def doublegaussian_m1_flat_m2(m1,m2,dL): # double Gaussian for m1, flat for m2, subject to m1 >= m2 convention

    if m1 >= m2: prior = double_gaussian(m1)
    else: prior = 0.
    
    return prior

def doublegaussian_m1m2_qpair(m1,m2,dL,beta=4.): # double Gaussian for m1, m2 with mass-ratio dependent pairing, subject to m1 >= m2 convention
    
    if m1 >= m2: prior = double_gaussian(m1)*double_gaussian(m2)*(m2/m1)**beta
    else: prior = 0.
    
    return prior

def doublegaussian_m1_flat_m2_qpair(m1,m2,dL,beta=4.): # double Gaussian for m1, flat for m2 with mass-ratio dependent pairing, subject to m1 >= m2 convention

    if m1 >= m2: prior = double_gaussian(m1)*(m2/m1)**beta
    else: prior = 0.
    
    return prior

def o3a_powerpeak_m1m2_qpair(m1,m2,dL,lpeak=0.10,alpha=2.62,beta=1.26,mmin=4.53,mmax=86.73,delta=4.88,mu=33.49,sigma=5.09): # power law + peak mass model from O3a populations paper
    
    prior = ((1.-lpeak)*powerlaw(m1,-alpha,mmax)+lpeak*gaussian(m1,mu,sigma))*smooth(m1,mmin,delta)*(m2/m1)**beta*smooth(m2,mmin,delta)/m1
    
    return prior

def o3a_powerbreak_m1m2_qpair(m1,m2,dL,alpha1=1.58,alpha2=5.59,beta=1.40,mmin=4.83,mmax=3.96,delta=87.14,b=0.43): # broken power law mass model from O3a populations paper
    
    mbreak = mmin+b*(mmax-mmin)
    
    if m1 >= mmin and m1 < mbreak: prior = m1**(-alpha1)*smooth(m1,mmin,delta)*(m2/m1)**beta*s(m2,mmin,delta)/m1
    elif m1 >= mbreak and m1 <= mmax: prior = m1**(-alpha2)*smooth(m1,mmin,delta)*(m2/m1)**beta*s(m2,mmin,delta)/m1
    else: prior = 0.
    
    return prior
