#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 15:56:25 2019

@author: Peter Mesirimel, Lund University
"""

from __future__ import division
import subprocess
from run_tolerances import run_tolerances
from postprocessing_cpp_output import process_output, delete_raw_data
from produce_plotting_data import produce_plotting_data

times = 100
tolerances = [10**(-i) for i in range(10)]
parameters = {'timesteps' : 100, 'macrosteps': 1, 'maxiter': 200, 'nsteps_conv_check': 3}

print('Starting run...')
print('times: ', times)
print('tolerances: ', tolerances)
print('parameters: ', parameters)

path = run_tolerances('toy', 'TOY', 'toy_test', times = times,
                      tolerances = tolerances, parameters = parameters,
                      runmodes = ['GS1', 'GS2', 'JAC', 'NEW'],
                      run_names = ['GS_1', 'GS_2', 'JAC', 'NEW'],
                      ref_run_name = 'GS_1',
                      labels = ['GS 1', 'GS 2', 'JAC', 'NEW'])

print('...processing output')
process_output(path)
print('...deleting raw data')
#delete_raw_data(path)
print('...producing plotting data')
produce_plotting_data(path)
print('...plotting')
subprocess.call('python3 plotting.py {} False'.format(path), shell = True)
print('...done')
