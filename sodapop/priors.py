#!/usr/bin/env python

import numpy as np
import astropy.units as u
from astropy.coordinates import Distance
from astropy.cosmology import Planck15 as cosmo

### BASIC FUNCTIONS

def src_to_det(m,z): # convert to detector frame mass from source frame
    
	return m*(1.+z)
    
def dL_to_z(dL): # convert distance to redshift, assuming a cosmology
    
	return Distance(dL,unit=u.Mpc).compute_z(cosmology=cosmo)
	
### BINARY MASS PRIORS

def flat_m1m2(m1,m2,dL): # uniform prior in source frame masses, subject to m1 >= m2 convention

	if m1 < m2: val = 0.
	else: val = 1.
    
	return val

def flat_m1m2_quad_dL(m1,m2,dL): # flat in source frame masses, subject to m1 >= m2 convention and quadratic prior in dL
    
	if m1 < m2: val = 0.
	else: val = dL**2
    
	return val

def flat_m1m2det(m1,m2,dL): # flat in detector frame masses, subject to m1 >= m2 convention and flat prior in dL
    
	if m1 < m2: val = 0.
	else: val = (1.+dL_to_z(dL))**2
    
	return val
	
def flat_m1m2det_quad_dL(m1,m2,dL): # flat in detector frame masses, subject to m1 >= m2 convention and quadratic prior in dL
    
	if m1 < m2: val = 0.
	else: val = dL**2*(1.+dL_to_z(dL))**2
    
	return val

def flat_mceta(m1,m2,dL): # flat in chirp mass and symmetric mass ratio
    
	if m1 < m2: val = 0.
	else: val = (m1-m2)*(m1*m2)**0.6/(m1+m2)**3.2
    
	return val

def flat_mcetadet(m1,m2,dL): # flat in chirp mass and symmetric mass ratio, assuming flat prior in dL
    
	if m1 < m2: val = 0.
	else: val = (1.+dL_to_z(dL))**2*(m1-m2)*(m1*m2)**0.6/(m1+m2)**3.2
    
	return val
    
def flat_mcq(m1,m2,dL): # flat in chirp mass and mass ratio
    
	if m1 < m2: val = 0.
	else: val = (m1*m2)**0.6/(m1**2*(m1+m2)**0.2)
    
	return val

def flat_mcqdet(m1,m2,dL): # flat in chirp mass and mass ratio, assuming flat prior in dL
    
	if m1 < m2: val = 0.
	else: val = (1.+dL_to_z(dL))*(m1*m2)**0.6/(m1**2*(m1+m2)**0.2)
    
	return val
	
def flat_mcqdet_quad_dL(m1,m2,dL): # flat in chirp mass and mass ratio, assuming quadratic prior in dL
    
	if m1 < m2: val = 0.
	else: val = dL**2*(1.+dL_to_z(dL))*(m1*m2)**0.6/(m1**2*(m1+m2)**0.2)
    
	return val

# BINARY MASS & SPIN PRIORS

def flat_mcqdet_quad_dL_flat_chieff_aligned(m1,m2,dL,chi1,chi2):

	chimax = 0.89
	chieff = (chi1+(m2/m1)*chi2)/(1.+(m2/m1))

	if (m1 < m2) or (np.abs(chi1) > chimax) or (np.abs(chi2) > chimax): val = 0.
	elif np.abs(chieff) <= chimax*(1.-(m2/m1))/(1.+(m2/m1)): val = 1.
	else: val = (1.-(1.-(m2/m1)/(1.+(m2/m1))))/(1.-np.abs(chieff)/chimax)*dL**2*(1.+dL_to_z(dL))*(m1*m2)**0.6/(m1**2*(m1+m2)**0.2)

	return val
	
# BINARY MASS & TIDAL DEFORMABILITY PRIORS

def flat_mcqdet_quad_dL_flat_LambdaT(m1,m2,dL,Lambda1,Lambda2):

	if (m1 < m2) or (Lambda2 < Lambda1): val = 0.
	else: val = np.exp(-((Lambda1-Lambda2*(m2/m1)**6)/(np.sqrt(2)*10.))**2)/(10.*np.sqrt(2*np.pi))*dL**2*(1.+dL_to_z(dL))*(m1*m2)**0.6/(m1**2*(m1+m2)**0.2)

	return val
	
# BINARY MASS, SPIN & TIDAL DEFORMABILITY PRIORS

def flat_mcqdet_quad_dL_flat_chieff_aligned_flat_LambdaT(m1,m2,dL,chi1,chi2,Lambda1,Lambda2):

	chimax = 0.89
	chieff = (chi1+(m2/m1)*chi2)/(1.+(m2/m1))

	if (m1 < m2) or (np.abs(chi1) > chimax) or (np.abs(chi2) > chimax) or (Lambda2 < Lambda1): val = 0.
	elif np.abs(chieff) <= chimax*(1.-(m2/m1))/(1.+(m2/m1)): val = np.exp(-((Lambda1-Lambda2*(m2/m1)**6)/(np.sqrt(2)*10.))**2)/(10.*np.sqrt(2*np.pi))*dL**2*(1.+dL_to_z(dL))*(m1*m2)**0.6/(m1**2*(m1+m2)**0.2)
	else: val = (1.-(1.-(m2/m1)/(1.+(m2/m1))))/(1.-np.abs(chieff)/chimax)*np.exp(-((Lambda1-Lambda2*(m2/m1)**6)/(np.sqrt(2)*10.))**2)/(10.*np.sqrt(2*np.pi))*dL**2*(1.+dL_to_z(dL))*(m1*m2)**0.6/(m1**2*(m1+m2)**0.2)

	return val
	
def flat_mcqdet_quad_dL_flat_chieff_flat_LambdaT(m1,m2,dL,chi1,chi2,cost1,cost2,Lambda1,Lambda2):

	chimax = 0.89
	chieff = (chi1*cost1+(m2/m1)*chi2*cost2)/(1.+(m2/m1))

	if (m1 < m2) or (np.abs(chi1) > chimax) or (np.abs(chi2) > chimax) or (np.abs(cost1) > 1.) or (np.abs(cost2) > 1.) or (Lambda2 < Lambda1): val = 0.
	elif np.abs(chieff) <= chimax*(1.-(m2/m1))/(1.+(m2/m1)): val = np.exp(-((Lambda1-Lambda2*(m2/m1)**6)/(np.sqrt(2)*10.))**2)/(10.*np.sqrt(2*np.pi))*dL**2*(1.+dL_to_z(dL))*(m1*m2)**0.6/(m1**2*(m1+m2)**0.2)
	else: val = (1.-(1.-(m2/m1)/(1.+(m2/m1))))/(1.-np.abs(chieff)/chimax)*np.exp(-((Lambda1-Lambda2*(m2/m1)**6)/(np.sqrt(2)*10.))**2)/(10.*np.sqrt(2*np.pi))*dL**2*(1.+dL_to_z(dL))*(m1*m2)**0.6/(m1**2*(m1+m2)**0.2)

	return val*4.*chimax*np.log(np.abs(chi1*cost1)/chimax)**(-1)*np.log(np.abs(chi2*cost2)/chimax)**(-1)

### PRIOR LOOKUP FUNCTIONS
    
binary_mass_priors = {'flat_m1m2': flat_m1m2, 'flat_m1m2det': flat_m1m2det, 'flat_mceta': flat_mceta, 'flat_mcetadet': flat_mcetadet, 'flat_m1m2det_quad_dL': flat_m1m2det_quad_dL, 'flat_m1m2_quad_dL': flat_m1m2_quad_dL, 'flat_mcq': flat_mcq , 'flat_mcqdet': flat_mcqdet, 'flat_mcqdet_quad_dL': flat_mcqdet_quad_dL, 'flat_mcqdet_quad_dL_flat_chieff_aligned': flat_mcqdet_quad_dL_flat_chieff_aligned, 'flat_mcqdet_quad_dL_flat_LambdaT': flat_mcqdet_quad_dL_flat_LambdaT, 'flat_mcqdet_quad_dL_flat_chieff_aligned_flat_LambdaT': flat_mcqdet_quad_dL_flat_chieff_aligned_flat_LambdaT, 'flat_mcqdet_quad_dL_flat_chieff_flat_LambdaT': flat_mcqdet_quad_dL_flat_chieff_flat_LambdaT}

def get_binary_mass_prior(key):

	try: prior_func = binary_mass_priors[key]
	except KeyError:
		print('Invalid prior specification, accepted keys are as follows:\n{0}\n'.format(binary_mass_priors.keys()))
		raise KeyError

	return prior_func
	
