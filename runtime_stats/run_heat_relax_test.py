#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 15:56:25 2019

@author: Peter Mesirimel, Lund University
"""

from __future__ import division
from run_tolerances import run_tolerances
from postprocessing_cpp_output import process_output, delete_raw_data
from produce_plotting_data import produce_plotting_data
import subprocess

times = 3
tolerances = [10**(-i) for i in range(7)]
parameters = {'timesteps' : 40, 'macrosteps': 2, 'maxiter': 1000,
              'gridsize': 32, 'alpha': 1, 'lambda': 0.1, 'tend': 0.2, 'u0': 2, 'nconv': 3,
              'w_relax_jac': 0.2, 'w_relax_gs': 0.6, 'times_only_new': 1, 'match_which_conv_relax': 1}

print('Starting run...')
print('times: ', times)
print('tolerances: ', tolerances)
print('parameters: ', parameters)
path = run_tolerances('heat_DN', 'heat_DN', 'heat_CN', times = times,
                      tolerances = tolerances, parameters = parameters,
                      runmodes = ['GS', 'GS', 'JAC', 'NEW', 'NEW'],
                      run_names = ['GS_DN', 'GS_ND', 'JAC', 'NEW_single', 'NEW_opt'],
                      first = [1]+[0]*4,
                      ref_run_name = 'GS_DN',
                      labels = ['GS_DN', 'GS_ND', 'JAC', 'NEW_base', 'NEW_opt'],
                      w_relax = [0.6, 0.6, 0.2, 0.4, 1],
                      new_relax_opt = [0]*4 + [1])

print('...processing output')
process_output(path)
print('...deleting raw data')
#delete_raw_data(path)
print('...producing plotting data')
produce_plotting_data(path)
print('...plotting')
subprocess.call('python3 plotting.py {}'.format(path), shell = True)
print('...done')