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

def run_tolerances(folder, exe, name, parameters, times = 10, run_prefix = '', # base inputs
                   tolerances = None, run_names = None, runmodes = None, ref_run_name = None, # mandatory list based inputs
                   errorlogging = False, ## run only last tolerances and extract other results from logs
                   **kwargs): # optional list based inputs
    # name as in base name for output and such
    # non-list based inputs go into "parameters"
    # list based (optional) parameters go into **kwargs
    if runmodes is None:   raise ValueError('Runmodes not specified')
    if tolerances is None: raise ValueError('Tolerances not specified')
    if run_names is None:  raise ValueError('run_names (name identifiers) not specified')
    for rr in run_names:
        if ' ' in run_names: raise ValueError('spaces in run_names not permitted')
    if len(run_names) != len(set(run_names)):
        raise ValueError('run_names must be destinct')
    
    if 'labels' not in kwargs.keys():
        labels = run_names.copy()
    else:
        labels = kwargs['labels'].copy()
        if len(labels) != len(run_names): raise ValueError('number of label names does not match number of run_names')
        del kwargs['labels']
    
    # verify all input lists are of same length
    s = set([len(kwargs[key]) for key in kwargs.keys()] + [len(i) for i in [run_names, runmodes]])
    if len(s) > 1: raise ValueError('Lengths of input lists do not match up')
    
    ## cut down tolerances if errorlogging
    tolerances_org = tolerances.copy()
    if errorlogging:
        tolerances = tolerances[-2:]
    
    time_string = str(datetime.datetime.now()).replace(' ', '_').replace(':', '_').replace('.', '_')[:-7]
    
    # create output directory
    if not os.path.exists('output'): os.makedirs('output')
    if "gridsize" in parameters.keys(): 
        output_dir = 'output/{}_{}_{}'.format(name, parameters['gridsize'], time_string)
    else:
        output_dir = 'output/{}_{}'.format(name, time_string)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    output_content = '{}/content.txt'.format(output_dir)
    output_parameters = '{}/parameters.txt'.format(output_dir)
    out_files = {'files': []}
    
    import threading
    shared_par = [folder, exe, output_dir, out_files, parameters, times, tolerances, ref_run_name, run_prefix]
    # SERIAL, 1 processor, 'GS' in runmode as indicator
    serial_kwargs = {'runmodes': [r for r in runmodes if 'GS' in r],
                     'run_names': [run_names[i] for i, r in enumerate(runmodes) if 'GS' in r],
                     'errorlogging': errorlogging}
    for key in kwargs.keys():
        serial_kwargs[key] = [kwargs[key][i] for i, r in enumerate(runmodes) if 'GS' in r]
        
    thread_serial = threading.Thread(target = run_tolerances_thread,
                                     args = shared_par + [1], kwargs = serial_kwargs)
    
    # PARALLEL, 2 processors
    par_kwargs = {'runmodes': [r for r in runmodes if 'GS' not in r],
                  'run_names': [run_names[i] for i, r in enumerate(runmodes) if 'GS' not in r],
                  'errorlogging': errorlogging}
    for key in kwargs.keys():
        par_kwargs[key] = [kwargs[key][i] for i, r in enumerate(runmodes) if 'GS' not in r]
        
    thread_parallel = threading.Thread(target = run_tolerances_thread,
                                       args = shared_par + [2], kwargs = par_kwargs)
    
    thread_serial.start()
    thread_parallel.start()
    thread_serial.join()
    thread_parallel.join() #to get most relevant 
                
    with open(output_content, 'w') as myfile:
        myfile.write(json.dumps(out_files, indent = 4))
    
    parameters['runmodes'] = runmodes
    parameters['folder'] = folder
    parameters['executable'] = exe
    parameters['name'] = name
    parameters['run_names'] = run_names
    parameters['reference_sol'] = ref_run_name
    parameters['labels'] = labels
    parameters['tolerances'] = tolerances_org
    parameters['times'] = times
    
    with open(output_parameters, 'w') as myfile:
        myfile.write(json.dumps(parameters, indent = 4, sort_keys = True))
    return output_dir + '/'

def run_tolerances_thread(folder, exe, output_dir, out_files, parameters, times, tolerances,
                          ref_run_name, run_prefix, num_proc, runmodes, run_names, errorlogging = False, **kwargs):
    for i, tol in enumerate(tolerances):
        last = False
        parameters['wftol'] = tol
        parameter_string = ' '.join([f'-{key} {parameters[key]}' for key in parameters.keys()])
        for j, r in enumerate(runmodes):
            string_full = parameter_string + ' ' + ' '.join([f'-{key} {kwargs[key][j]}' for key in kwargs.keys()])
            
            ## if errorlogging, skip first tolerance for reference solution
            if errorlogging and (i == 0) and (run_names[j] == ref_run_name):
                continue
            
            ## last tolerance, skip for non reference solution
            if i == len(tolerances) - 1: 
                if run_names[j] != ref_run_name:
                    continue
                else:
                    last = True
            
            out_f = output_dir + '/' + str(run_names[j]) + '_' + str(tol) + '.txt'
            out_files['files'].append(out_f)
            
            run_string = (f'{run_prefix} mpirun -np {num_proc} -quiet ../src/{folder}/{exe} -runmode {r} {string_full} >> {out_f}\n')
#            print(run_string)
            print(f'running {r} ({run_names[j]}) for tolerance of {tol}')
            for _ in range(times):
                subprocess.call(run_string, shell = True)
                if last: 
                    break
                try:
                    if parameters['times_only_new'] and 'NEW' not in r:
                        break
                except KeyError: # in case times_only_new not set, noting bad should happen
                    pass