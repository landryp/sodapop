#!/usr/bin/env python

import numpy as np

### PRIOR, POSTERIOR AND LIKELIHOOD FOR EMCEE

def logprior(lambdaa, lambda_dict, prior_dict): # log of prior probability of population parameters lambda
	
	log_prior = 0.
	for prior_key in prior_dict.keys():
	
		prior_func, hyp_params = prior_dict[prior_key]
		lambda_subset = [lambda_dict[param](lambdaa) for param in prior_key.split(',')]
		
		log_prior += np.log(prior_func(*lambda_subset,*hyp_params))

	return log_prior
	
def loglikelihood(lambdaa, lambda_dict, lambdabh, obs_data, obs_classes, pop_models, sel_samps, sel_func, dist_func, dist_params, detfrac_dict, likelihood_dict): # log of likelihood of population parameters lambda

	num_obs = len(obs_data)

	bns_lambda = [lambda_dict[param](lambdaa) for param in lambda_dict.keys()]
	if 'beta' in lambda_dict.keys(): nsbh_lambda = bns_lambda[:-1]+lambdabh
	else: nsbh_lambda = bns_lambda+lambdabh
	
	lambdas = {'bns': bns_lambda, 'nsbh': nsbh_lambda}
	
	log_like = 0.
	for i,likedata in enumerate(obs_data):
	
		obs_class = obs_classes[i]
		pop_model = pop_models[obs_class]
		lambda_subset = lambdas[obs_class]
		
		m1s,m2s,dls,zs = likedata
		like_num = np.sum(pop_model(m1s,m2s,lambda_subset) * dist_func(dls,*dist_params))
			
		log_like += np.log(like_num)
		
		if not np.isfinite(log_like): return -np.inf, -np.inf

	m1s,m2s,dls,zs = sel_samps
	nsbh_pop_probs = pop_models['nsbh'](m1s,m2s,lambdas['nsbh'])
	bns_pop_probs = pop_models['bns'](m1s,m2s,lambdas['bns'])
	
	pop_probs = np.where(m1s > lambdabh[0], nsbh_pop_probs, bns_pop_probs) # use NSBH model if sample in NSBH region
			
	like_denom = np.sum(sel_func(m1s,m2s,dls,zs) * pop_probs * dist_func(dls,*dist_params))

	log_like -= num_obs*np.log(like_denom)
	
	if not np.isfinite(log_like): return -np.inf, np.log(like_denom)

	return log_like, np.log(like_denom)
	
def logposterior(lambdaa, lambda_dict, lambdabh, obs_data, obs_classes, pop_models, sel_samps, sel_func, prior_dict, dist_func, dist_params, detfrac_dict, likelihood_dict): # log of posterior probability of population parameters lambda

	log_prior = logprior(lambdaa, lambda_dict, prior_dict)
	
	if not np.isfinite(log_prior): return -np.inf, list(lambdaa), -np.inf, -np.inf
	
	log_like, log_detfrac = loglikelihood(lambdaa, lambda_dict, lambdabh, obs_data, obs_classes, pop_models, sel_samps, sel_func, dist_func, dist_params, detfrac_dict, likelihood_dict)

	return log_prior + log_like, list(lambdaa), log_like, log_detfrac
	
