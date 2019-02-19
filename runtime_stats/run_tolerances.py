#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 12:53:31 2019

@author: Peter Mesirimel, Lund University
"""

from __future__ import division
import os
import subprocess
import datetime
import json

def run_tolerances(folder, exe, name, times = 10, tolerances = None, parameters = None, runmodes = None):
    
    if runmodes is None:
        runmodes = ['GS', 'JAC', 'NEW']
    if tolerances is None:
        tolerances = [10**(-i) for i in range(6)]
    if parameters is None:
        parameters = {'timesteps1' : 100, 'timesteps2' : 100, 'macrosteps': 5, 'maxiter': 50}
        
    tt = times
    #########
    time_string = str(datetime.datetime.now()).replace(' ', '_')
    time_string = time_string.replace(':', '_')
    #########
    # create output directory
    if not os.path.exists('output'):
        os.makedirs('output')
    if "gridsize" in parameters.keys():
        output_dir = 'output/{}_{}_{}'.format(name, parameters['gridsize'], time_string)
    else:
        output_dir = 'output/{}_{}'.format(name, time_string)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    #########
    output_content = '{}/content.txt'.format(output_dir)
    output_parameters = '{}/parameters.txt'.format(output_dir)
    out_files = {'files': []}
    
    for v, tol in enumerate(tolerances):
        parameters['wftol'] = tol
        parameter_string = ' '.join(['-' + str(key) + ' ' + str(parameters[key]) for key in parameters.keys()])
        for r in runmodes:
            if v == len(tolerances) - 1 :
                if r != 'GS':
                    continue
                else:
                    times = 1
            out_f = output_dir + '/' + str(r) + '_' + str(tol) + '.txt'
            out_files['files'].append(out_f)
            
            run_string = ('mpirun -np 2 ../src/' + folder + '/' + exe +
                         ' -runmode ' + r + ' ' + parameter_string + ' >> ' + out_f + '\n')
            print('running', r, 'for tolerance of', tol)
            for i in range(times):
                subprocess.call(run_string, shell = True)
                
    with open(output_content, 'w') as myfile:
        myfile.write(json.dumps(out_files, indent = 4))

    parameters['runmodes'] = runmodes
    parameters['tolerances'] = tolerances
    parameters['times'] = tt        
    del parameters['wftol']
    
    with open(output_parameters, 'w') as myfile:
        myfile.write(json.dumps(parameters, indent = 4, sort_keys = True))
    return output_dir + '/'
