#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 13:31:49 2019

@author: Peter Mesirimel, Lund University
"""

from __future__ import division
from postprocessing_cpp_output import split_and_float
import subprocess
import os
import json
import numpy as np
import pylab as pl
import datetime
pl.close('all')

tol = 1e-6
gridsize = 16
name = 'heat_{}'.format(gridsize)
folder = 'heat_DN'
exe = 'heat_DN'

time_string = str(datetime.datetime.now()).replace(' ', '_').replace(':', '_').replace('.', '_')

if not os.path.exists('output'):
    os.makedirs('output')
output_dir = 'output/order_verification_{}_{}'.format(name, time_string)
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
#########
output_parameters = '{}/parameters.txt'.format(output_dir)
out_files = {'files': []}

print('Starting order verification run...')
steps_list = [2**i for i in range(10)]
for steps in steps_list:
    parameters = {'timesteps' : steps, 'macrosteps': 1, 'maxiter': 1000,
                  'gridsize': gridsize, 'alpha': 1, 'gamma': 0.01,
                  'logging': 0, 'runmode': 'GS', 'wftol': tol}
    
    parameter_string = ' '.join(['-' + str(key) + ' ' + str(parameters[key]) for key in parameters.keys()])
    out_f = output_dir + '/' + str(steps) + '.txt'
    out_files['files'].append(out_f)
            
    run_string = ('mpirun -np 2 ../src/' + folder + '/' + exe + ' ' + parameter_string + ' >> ' + out_f + '\n')
    print('running with {} steps'.format(steps))
    subprocess.call(run_string, shell = True)
                
parameters['steps_list'] = steps_list
with open(output_parameters, 'w') as myfile:
    myfile.write(json.dumps(parameters, indent = 4, sort_keys = True))
""" process data """

dic = {}
path = output_dir + '/'
with open(path + 'parameters.txt', 'r') as myfile:
    dic = json.load(myfile)
files = [output_dir + '/' + str(i) for i in dic['steps_list']]

sol1 = []
sol2 = []
    
for steps in steps_list:
    data = 0
    with open(output_dir + '/' + str(steps) + '.txt', 'r') as myfile:
        data = myfile.read().split('\n')[:-1]
    _, _, vec1, vec2 = split_and_float(data[0], data[1])
    sol1.append(np.array(vec1))
    sol2.append(np.array(vec2))

""" calculate errors """
errors1, errors2 = [], []
ref1, ref2 = sol1[-1], sol2[-1]

for i in range(len(steps_list) - 1):
    errors1.append(np.linalg.norm(sol1[i] - ref1, 2))
    errors2.append(np.linalg.norm(sol2[i] - ref2, 2))
    
""" plotting """
plot_path = 'plots_' + path[path.find('/') + 1:-1]
if not os.path.exists(plot_path):
    os.makedirs(plot_path)
fig, ax = pl.subplots(figsize = (12, 10))

fs = 32
style = {'linestyle': '-', 'lw': 3, 'markersize': 17}

ax.loglog(steps_list[:-1], errors1, color = 'g', marker = 'o', label = 'Flux error', **style)
ax.loglog(steps_list[:-1], errors2, color = 'b', marker = 'd', label = 'Interface temp error', **style)
ax.loglog(steps_list, [1/i for i in steps_list], '--k', label = 'Order 1')
ax.loglog(steps_list, [1/i**2 for i in steps_list], ':k', label = 'Order 2')
ax.axhline(tol, color = 'r', label = 'WFR Tol', **style)
ax.grid(b = True, which = 'major')
        
ax.set_title('Order verification', fontsize = fs + 2)
ax.set_xlabel('Steps', fontsize = fs)
ax.set_ylabel('Error', fontsize = fs)

ax.tick_params(labelsize = 20)
ax.legend(fontsize = fs - 12, loc = 0)
fig.savefig(plot_path + '/tolerance_vs_iters.eps', dpi = 100)
fig.savefig(plot_path + '/tolerance_vs_iters.png', dpi = 100)