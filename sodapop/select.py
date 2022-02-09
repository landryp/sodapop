#!/usr/bin/env python

import numpy as np
import numpy.random
import h5py
import astropy.units as u
from astropy.coordinates import Distance
from astropy.cosmology import Planck15 as cosmo
import os

### SELECTION FUNCTIONS
	
def chirpmass52(m1,m2,dl,z): # snr selection proportional to (chirp mass)^5/2

	if np.isscalar(m1): m1 = np.array([m1])
	else: m1 = np.array(m1)
	if np.isscalar(m2): m2 = np.array([m2])
	else: m2 = np.array(m2)
	if np.isscalar(z): z = np.array([z])
	else: z = np.array(z)
	if np.isscalar(dl): dl = np.array([dl])
	else: dl = np.array(dl)

	Mc = (m1*m2)**(3./5.)/(m1+m2)**(1./5.)
	Mcdet = Mc*(1.+z)
	x = Mcdet**(5./2.)

	return x

select_funcs = {'chirpmass52': chirpmass52}
	
### TOOLS FOR HANDLING SELECTION FUNCTIONS

def get_select_func(key):

	try: select_func = select_funcs[key]
	except KeyError:
		print('Invalid selection function specification, accepted keys are as follows:\n{0}\n'.format(select_funcs.keys()))
		raise KeyError

	return select_func

def load_select_effect(sel_func_name,sel_priors,num_sel):

	if sel_func_name == False or sel_func_name == 'False': sel_func = lambda *x : 1.
	else: sel_func = get_select_func(sel_func_name)
	
	m1lb,m1ub,m2lb,m2ub,dllb,dlub = [float(val) for val in (sel_priors).split(',')]
	
	m1_u = np.random.uniform(m1lb,m1ub,num_sel)
	m2_u = np.random.uniform(m2lb,m2ub,num_sel)
	dl_select = np.random.uniform(dllb,dlub,num_sel)
	z_select = [Distance(dL,unit=u.Mpc).compute_z(cosmology=cosmo) for dL in dl_select]
	
	m1_select = []
	m2_select = []
	for m1u,m2u in zip(m1_u,m2_u):
		while m1u < m2u:
			m1u = np.random.uniform(m1lb,m1ub,1)[0]
			m2u = np.random.uniform(m2lb,m2ub,1)[0]
			
		m1_select += [m1u]
		m2_select += [m2u]
	
	sel_samps = [np.array(m1_select), np.array(m2_select), np.array(dl_select), np.array(z_select)] # uniform m1,m2,dl samples subject to m1 > m2

	return sel_func, sel_samps
	
