#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns
import pandas as pd

### DIAGNOSTIC PLOTS

def corner_plot(samples,pop_param_names,fixed_params,outpath):
	
	params = [param for param in pop_param_names if param not in fixed_params]
	
	samps_frame = pd.DataFrame(samples,columns=params)
	
	fig = sns.pairplot(samps_frame, kind='scatter', markers='.', diag_kind='kde', corner=True, diag_kws=dict(lw=1,common_norm=False))
	fig.map_lower(sns.kdeplot, levels=[0.1,0.5])
	
	plt.savefig(outpath)

	return 0
	
def trace_plot(raw_samples,pop_param_names,fixed_params,outpath,num_burn=False):
	
	params = [param for param in pop_param_names if param not in fixed_params]
	num_params = len(params)
	
	fig = plt.figure(figsize=(8.,6.))
	gs = gridspec.GridSpec(num_params, 1)
	axs = [plt.subplot(gs[i]) for i in range(num_params)]
	plt.subplots_adjust(hspace=0.05)
	
	for i,ax in enumerate(axs):
		for chain_samples in raw_samples: ax.plot(chain_samples[:,i], c=sns.color_palette()[0], alpha=0.2)
		ax.set_ylabel(params[i])
		if num_burn: ax.axvline(num_burn,lw=1,ls='--',c='k')
		if i < len(axs)-1: ax.tick_params(labelbottom=False)
		
	axs[-1].set_xlabel('steps')
	
	plt.savefig(outpath)

	return 0
	
def ppd_plot(ppds,outpath,med=False,logscale=True):
	
	fig = plt.figure(figsize=(8.,6.))
	gs = gridspec.GridSpec(3, 1)
	axs = [plt.subplot(gs[i]) for i in range(3)]
	plt.subplots_adjust(hspace=0.05)
	
	if med: cent = [2,6,10,14,18,22]
	else: cent = [1,5,9,13,17,21]
	
	axs[0].plot(ppds[:,0],ppds[:,cent[0]],c=sns.color_palette()[0],label='pop')
	axs[0].fill_between(ppds[:,0],ppds[:,3],ppds[:,4],color=sns.color_palette()[0],alpha=0.1,lw=0.)
	
	axs[0].plot(ppds[:,0],ppds[:,cent[1]],c=sns.color_palette()[0],ls='--',label='obs')
	axs[0].plot(ppds[:,0],ppds[:,7],c=sns.color_palette()[0],lw=0.5,ls='--')
	axs[0].plot(ppds[:,0],ppds[:,8],c=sns.color_palette()[0],lw=0.5,ls='--')
	
	axs[1].plot(ppds[:,0],ppds[:,cent[2]],c=sns.color_palette()[0])
	axs[1].fill_between(ppds[:,0],ppds[:,11],ppds[:,12],color=sns.color_palette()[0],alpha=0.1,lw=0.)
	
	axs[1].plot(ppds[:,0],ppds[:,cent[3]],c=sns.color_palette()[0],ls='--')
	axs[1].plot(ppds[:,0],ppds[:,15],c=sns.color_palette()[0],lw=0.5,ls='--')
	axs[1].plot(ppds[:,0],ppds[:,16],c=sns.color_palette()[0],lw=0.5,ls='--')
	
	axs[2].plot(ppds[:,0],ppds[:,cent[4]],c=sns.color_palette()[0])
	axs[2].fill_between(ppds[:,0],ppds[:,19],ppds[:,20],color=sns.color_palette()[0],alpha=0.1,lw=0.)
	
	axs[2].plot(ppds[:,0],ppds[:,cent[5]],c=sns.color_palette()[0],ls='--')
	axs[2].plot(ppds[:,0],ppds[:,23],c=sns.color_palette()[0],lw=0.5,ls='--')
	axs[2].plot(ppds[:,0],ppds[:,24],c=sns.color_palette()[0],lw=0.5,ls='--')
	
	if logscale:
		axs[0].set_yscale('log')
		axs[1].set_yscale('log')
		axs[2].set_yscale('log')

		axs[0].set_ylim(0.1,5.)
		axs[1].set_ylim(0.1,5.)
		axs[2].set_ylim(0.1,5.)
	
	axs[0].set_ylabel('ppd(m)')
	axs[1].set_ylabel('ppd(m1)')
	axs[2].set_ylabel('ppd(m2)')
	
	axs[0].tick_params(labelbottom=False)
	axs[1].tick_params(labelbottom=False)	
	axs[2].set_xlabel('m')
	
	axs[0].legend()
	
	plt.savefig(outpath)

	return 0
	
