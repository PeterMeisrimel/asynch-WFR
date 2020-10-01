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

from heat_relax_aux import get_parameters, heat_prob

times = 5
## need to use larger tol for air_steel, 11, 7 otherwise
#tolerances = [10**(-i) for i in range(11)]
#heat_para = get_parameters('air_steel')

#tolerances = [10**(-i) for i in range(11)]
#which = 'air_steel'
#heat_para = get_parameters(which)

tolerances = [10**(-i) for i in range(6)]
which = 'test'
heat_para = get_parameters(which)

parameters = {'timesteps' : 40, 'macrosteps': 1, 'maxiter': 50,
              'gridsize': 32, 'tend': 10000, 'u0': 2,
              'times_only_new': True, ## only let NEW method run <times> times
              'nsteps_conv_check': 1, ## steps below tolerance for convergence check
              **heat_para}

heat_relax_class = heat_prob(parameters['gridsize'], **heat_para)
dt = parameters['tend']/parameters['timesteps']

theta_gs = heat_relax_class.GS_theta_opt(dt, dt)
theta_jac = heat_relax_class.JACOBI_theta_opt(dt, dt)

relax = {'theta_relax1': theta_jac, 'theta_relax2': theta_jac, # jacobi relax
         'var_relax': 1,
         'theta_relax_gs_a_1': 1, 'theta_relax_gs_a_2': theta_gs, # 1->2 gs
         'theta_relax_gs_b_1': theta_gs, 'theta_relax_gs_b_2': 1  # 2->1 gs
         }

parameters.update(relax)

print('Starting run...')
print('times: ', times)
print('tolerances: ', tolerances)
print('parameters: ', parameters)
#path = run_tolerances('heat_DN', 'heat_DN', 'heat_CN_{}'.format(which), times = times,
#                      tolerances = tolerances, parameters = parameters,
##                      run_prefix = 'nohup',
#                      runmodes = ['GS1', 'GS2', 'JAC', 'NEW', 'NEW', 'NEW'],
#                      run_names = ['GS_DN', 'GS_ND', 'JAC', 'NEW0', 'NEW1', 'NEW2'],
#                      ref_run_name = 'GS_DN',
#                      labels = ['GS_DN', 'GS_ND', 'JAC', 'NEW0', 'NEW1', 'NEW2'],
#                      var_relax = [0, 0, 0, 0, 1, 2]
#                      )
path = run_tolerances('heat_DN', 'heat_DN', 'heat_CN_{}'.format(which), times = times,
                      tolerances = tolerances, parameters = parameters,
#                      run_prefix = 'nohup',
                      runmodes = ['GS1', 'NEW'],
                      run_names = ['GS_DN', 'NEW2'],
                      ref_run_name = 'GS_DN',
                      labels = ['GS_DN', 'NEW2'],
                      var_relax = [0, 2]
                      )
print('...processing output')
process_output(path)
print('...deleting raw data')
#delete_raw_data(path)
print('...producing plotting data')
produce_plotting_data(path)
print('...plotting')
subprocess.call('python3 plotting.py {}'.format(path), shell = True)
print('...done')
