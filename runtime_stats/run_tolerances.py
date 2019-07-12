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

def run_tolerances(folder, exe, name, times = 10, tolerances = None,
                   parameters = None, runmodes = None,
                   first = None, run_names = None, ref_run_name = 'GS',
                   labels = None, times_only_new = False):
    # set times_only_new = true to only run new method <times> many times
    # other methodss are deterministic, several runs only give mean runtime
    
    if runmodes is None:
        runmodes = ['GS', 'JAC', 'NEW']
    if tolerances is None:
        tolerances = [10**(-i) for i in range(6)]
    if parameters is None:
        parameters = {'timesteps1' : 100, 'timesteps2' : 100, 'macrosteps': 5, 'maxiter': 50}
    if first is None:
        first = [True]*len(runmodes)
    if run_names is None:
        run_names = runmodes
    if labels is None:
        labels = run_names
        
        
    for rr in run_names:
        if ' ' in run_names:
            raise ValueError('spaces in run_names not permitted')
        
    parameters_run = parameters.copy()
    parameters['runmodes'] = runmodes
    parameters['folder'] = folder
    parameters['executable'] = exe
    parameters['name'] = name
    parameters['run_names'] = run_names
    parameters['reference_sol'] = ref_run_name
    parameters['labels'] = labels
        
    tt = times
    #########
    time_string = str(datetime.datetime.now()).replace(' ', '_').replace(':', '_').replace('.', '_')
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
        parameters_run['wftol'] = tol
        parameter_string = ' '.join(['-' + str(key) + ' ' + str(parameters_run[key]) for key in parameters_run.keys()])
        for w, r in enumerate(runmodes):
            parameter_string_full = parameter_string + ' -first ' + str(int(first[w]))
            if v == len(tolerances) - 1 :
                if run_names[w] != ref_run_name:
                    continue
                else:
                    times = 1
            out_f = output_dir + '/' + str(run_names[w]) + '_' + str(tol) + '.txt'
            out_files['files'].append(out_f)
            
            run_string = ('mpirun -np 2 ../src/' + folder + '/' + exe +
                         ' -runmode ' + r + ' ' + parameter_string_full + ' >> ' + out_f + '\n')
            #print(run_string)
            print('running', r, '({})'.format(run_names[w]), 'for tolerance of', tol)
            for i in range(times):
                subprocess.call(run_string, shell = True)
                if times_only_new and 'NEW' not in r:
                    break
                
    with open(output_content, 'w') as myfile:
        myfile.write(json.dumps(out_files, indent = 4))

    parameters['runmodes'] = runmodes
    parameters['tolerances'] = tolerances
    parameters['times'] = tt        
    
    with open(output_parameters, 'w') as myfile:
        myfile.write(json.dumps(parameters, indent = 4, sort_keys = True))
    return output_dir + '/'
