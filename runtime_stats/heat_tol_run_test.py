#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 15:56:25 2019

@author: Peter Mesirimel, Lund University
"""

from __future__ import division
from run_tolerances import run_tolerances
from postprocessing_cpp_output import process_output_errlog, delete_raw_data
from produce_plotting_data import produce_plotting_data
import subprocess

from heat_relax_aux import get_parameters, get_relax_parameters

times = 3

tolerances = [10**(-i) for i in range(9)]
which = 'test'
heat_para = get_parameters(which)
## adjust tend depending on problem!!!

parameters = {'timesteps1' : 20, 'timesteps2': 20, 'macrosteps': 1, 'maxiter': 100,
              'gridsize': 32, 'tend': 1, 'u0': 2, 'errlog': 1, 
              'times_only_new': True, ## only let NEW method run <times> times
              'nsteps_conv_check': 1, ## steps below tolerance for convergence check
              'var_relax': 1,
              **heat_para}

relax = get_relax_parameters(**parameters)
parameters.update(relax)

print('Starting run...')
print('times: ', times)
print('tolerances: ', tolerances)
print('parameters: ', parameters)
path = run_tolerances('heat_DN', 'heat_DN', 'heat_CN_{}'.format(which), times = times,
                      tolerances = tolerances, parameters = parameters,
#                      run_prefix = 'nohup',
                      runmodes = ['GS1', 'JAC', 'NEW', 'NEW'],
                      run_names = ['GS_DN', 'JAC', 'NEW0', 'NEW2'],
                      ref_run_name = 'GS_DN',
                      labels = ['GS_DN', 'JAC', 'NEW0', 'NEW2'],
                      var_relax = [0, 0, 0, 2],
                      errorlogging = True
                      )

print('...processing output')
process_output_errlog(path, which = 2)
#print('...deleting raw data')
##delete_raw_data(path)
print('...producing plotting data')
produce_plotting_data(path)
print('...plotting')
subprocess.call('python3 plotting.py {}'.format(path), shell = True)
print('...done')
