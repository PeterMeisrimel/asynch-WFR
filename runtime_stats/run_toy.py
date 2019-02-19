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

times = 5
tolerances = [10**(-i) for i in range(8)]
parameters = {'timesteps1' : 50, 'timesteps2' : 50, 'macrosteps': 5, 'maxiter': 20}

print('Starting run...')
print('times: ', times)
print('tolerances: ', tolerances)
print('parameters: ', parameters)

path = run_tolerances('toy', 'TOY', 'toy_test', times = times, tolerances = tolerances, parameters = parameters)

print('...processing output')
process_output(path)
print('...deleting raw data')
delete_raw_data(path)
print('...producing plotting data')
produce_plotting_data(path)
print('...plotting')
subprocess.call('python3 plotting.py {}'.format(path), shell = True)
print('...done')