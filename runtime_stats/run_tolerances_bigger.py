#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 12:53:31 2019

@author: Peter Mesirimel, Lund University
"""

from __future__ import division
import numpy as np
import os
import subprocess
import datetime
import json

#########
problem_name = 'bigger' # for output file only
problem_folder = 'bigger'
problem_executable = 'BIGGER'

times = 50
runmodes = ['GS', 'JAC', 'NEW']
tolerances = [10**(-i) for i in range(9)]
parameters = {'timesteps1' : 50, 'timesteps2' : 50, 'macrosteps': 5, 'maxiter': 50}
#########

time_string = str(datetime.datetime.now()).replace(' ', '_')
output_file_tmp = 'output_{}'.format(time_string)
out_files = []
output_file = 'results_{}_{}.txt'.format(problem_name, time_string)

results_list = {}
for r in runmodes:
    results_list[r] = {'tols': tolerances, 'iter': [], 't': [], 'sol': []}

for v, tol in enumerate(tolerances):
    parameters['wftol'] = tol
    parameter_string = ' '.join(['-' + str(key) + ' ' + str(parameters[key]) for key in parameters.keys()])
    for r in runmodes:
        if v == len(tolerances) - 1 and r != 'GS':
            continue
        out_f = output_file_tmp + str(r) + str(tol) + '.txt'
        out_files.append(out_f)
        
        run_string = ('mpirun -np 2 ../src/' + problem_folder + '/' + problem_executable +
                     ' -runmode ' + r + ' ' + parameter_string + ' >> ' + out_f + '\n')
        print('running', r, 'for tolerance of', tol)
        for i in range(times):
            subprocess.call(run_string, shell = True)
        # process output file
        content = 0
        with open(out_f, 'r') as f:
            content = f.read().split('\n')[:-1]
        p0, p1 = [], []
        for line in content:
            if line[0] == '0':
                p0.append(line.split())
            else:
                p1.append(line.split())
        results_list[r]['iter'].append(np.mean([float(p0[i][1]) for i in range(times)]))
        results_list[r]['t']   .append(np.mean([max(float(p0[i][2]), float(p1[i][2])) for i in range(times)]))
        sol = []
        for i, j in enumerate(p0[0][3:]):
            sol.append(np.mean([float(p0[k][i+3]) for k in range(times)]))
        for i, j in enumerate(p1[0][3:]):
            sol.append(np.mean([float(p1[k][i+3]) for k in range(times)]))
        results_list[r]['sol'].append(sol)
        
with open(output_file, 'w') as myfile:
    myfile.write(json.dumps(results_list, sort_keys = True, indent = 4))

for i in out_files:
    os.remove(i)