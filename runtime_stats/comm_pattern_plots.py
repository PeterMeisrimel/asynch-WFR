#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 13:10:54 2019

@author: Peter Mesirimel, Lund University
"""

from __future__ import division
import numpy as np
import pylab as pl
import json

def plot_comm_pattern(path):
    parameters = {}
    with open(path + 'parameters.txt', 'r') as myfile:
        parameters = json.load(myfile)
    
    data = {}
    with open(path + 'comm.txt', 'r') as myfile:
        data = json.load(myfile)
        
    ## create necessary directories
    import os
    plot_path = 'plots_' + path[path.find('/') + 1:-1]
    dirs = [plot_path, plot_path + '/comm']
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)
    
    data = data['NEW']
    for data_set in data:
        tol = data_set['tol']
        
        msteps = parameters['macrosteps']
        iters_per_macro = data_set['mlog']
        log0 = data_set['log0']
        log1 = data_set['log1']
        
        # get averages
        steps_per_macro_0 = int(parameters['timesteps1']/msteps)
        steps_per_macro_1 = int(parameters['timesteps2']/msteps)
        
        log0_avg = [np.zeros(steps_per_macro_0) for i in range(msteps)]
        log1_avg = [np.zeros(steps_per_macro_1) for i in range(msteps)]
        for t in range(parameters['times']):
            for i, iters in enumerate(iters_per_macro[t]):
                for j in range(iters):
                    log0_avg[i] += np.array(log0[t][i][j])/iters
                    log1_avg[i] += np.array(log1[t][i][j])/iters
        
        for i in range(msteps):
            log0_avg[i] /= parameters['times']
            log1_avg[i] /= parameters['times']
        
        fig, ax = pl.subplots(2, figsize = (12, 10))
        for i in range(msteps):
            ax[0].bar(range(steps_per_macro_0*i, steps_per_macro_0*(i+1)), log0_avg[i], color = 'blue')
            ax[0].plot([steps_per_macro_0*(i+1)]*2, [0, 1], '--k')
            ax[1].bar(range(steps_per_macro_1*i, steps_per_macro_1*(i+1)), log1_avg[i], color = 'red')
            ax[1].plot([steps_per_macro_1*(i+1)]*2, [0, 1], '--k')
            
        for i in range(2):
            ax[i].axhline(0, linestyle = '-', color = 'black')
            ax[i].axhline(1, linestyle = '--', color = 'black')
            ax[i].set_ylim([-0.1, 1.1])
            ax[i].set_xlim([0, max(parameters['timesteps1'], parameters['timesteps2'])])
            ax[i].set_ylabel('% use new data')
            ax[i].set_xlabel('timestep')
        ax[0].set_title('Processor 0: Dirichlet, WF tol = {}'.format(tol))
        ax[1].set_title('Processor 1: Neumann, WF tol = {}'.format(tol))
        fig.savefig(plot_path + '/comm/tol_{}.png'.format(tol), dpi = 100)
        