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

def flat_mcqdet_quad_dL_flat_chieff(m1,m2,dL,chi1,chi2):

	if (m1 < m2) or (chi1 > 1.) or (chi1 < 0.) or (chi2 > 1.) or (chi2 < 0.): val = 0.
	else: val = (1.)*dL**2*(1.+dL_to_z(dL))*(m1*m2)**0.6/(m1**2*(m1+m2)**0.2) #FIXME: add p(chi1,chi2) given flat chieff and aligned/isotropic spins

	return val
	
# BINARY MASS & TIDAL DEFORMABILITY PRIORS

def flat_mcqdet_quad_dL_flat_LambdaT(m1,m2,dL,Lambda1,Lambda2):

	if (m1 < m2) or (chi1 > 1.) or (chi1 < 0.) or (chi2 > 1.) or (chi2 < 0.): val = 0.
	else: val = (1.)*dL**2*(1.+dL_to_z(dL))*(m1*m2)**0.6/(m1**2*(m1+m2)**0.2) #FIXME: add p(Lambda1,Lambda2) given flat LambdaT

	return val

### PRIOR LOOKUP FUNCTIONS
    
binary_mass_priors = {'flat_m1m2': flat_m1m2, 'flat_m1m2det': flat_m1m2det, 'flat_mceta': flat_mceta, 'flat_mcetadet': flat_mcetadet, 'flat_m1m2det_quad_dL': flat_m1m2det_quad_dL, 'flat_m1m2_quad_dL': flat_m1m2_quad_dL, 'flat_mcq': flat_mcq , 'flat_mcqdet': flat_mcqdet, 'flat_mcqdet_quad_dL': flat_mcqdet_quad_dL}

def get_binary_mass_prior(key):

	try: prior_func = binary_mass_priors[key]
	except KeyError:
		print('Invalid prior specification, accepted keys are as follows:\n{0}\n'.format(binary_mass_priors.keys()))
		raise KeyError

	return prior_func
	
