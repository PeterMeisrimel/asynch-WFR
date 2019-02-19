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
from comm_pattern_plots import plot_comm_pattern
import subprocess

times = 10
tolerances = [10**(-i) for i in range(3, 8)]

gs = [32, 64, 128]

for steps in [50, 100]:
    for g in gs:
        parameters = {'timesteps1' : steps, 'timesteps2' : steps, 'macrosteps': 5, 'maxiter': 20,
                      'gridsize': g, 'alpha': 1, 'gamma': 0.01, 'logging': 1}
    
        print('Starting run...')
        print('times: ', times)
        print('tolerances: ', tolerances)
        print('parameters: ', parameters)
        
        path = run_tolerances('heat_DN', 'heat_DN', 'heat', times = times, tolerances = tolerances, parameters = parameters)
        print('...processing output')
        process_output(path)
        print('...deleting raw data')
        delete_raw_data(path)
        if parameters['logging']:
            print('...creating communication pattern log')
            plot_comm_pattern(path)
        print('...producing plotting data')
        produce_plotting_data(path)
        print('...plotting')
        subprocess.call('python3 plotting.py {}'.format(path), shell = True)
        print('...done')