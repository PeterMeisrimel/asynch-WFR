#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 13:21:10 2019

@author: Peter Mesirimel, Lund University
"""

from __future__ import division
import os
import json
import numpy as np
                
def cut_path(f):
    return f[-f[::-1].find('/'):]

def split_and_float(f1, f2):
    a = [float(i) for i in f1.split(' ')[:-1]]
    b = [float(i) for i in f2.split(' ')[:-1]]
    assert a[0] != b[0], 'file output incorrect or corrupted, {}, {}'.format(a[0], b[0])
    assert a[1] == b[1], 'file output incorrect or corrupted, {}, {}'.format(a[1], b[1])
    return a[1], (a[2] + b[2])/2, a[3:], b[3:]

def process_log(log):
    log = [int(i) for i in log[4:-1].split(' ')]
    p = log[0]
    timesteps = log[1] # total
    macro_steps = log[2]
    macro_list = [i for i in log[3:3+macro_steps]]
    log = log[3+macro_steps:]
    timesteps = int(timesteps/macro_steps) # per macro step
    
    log_list = []
    s = 0
    for m in macro_list:
        log_list.append([log[(s+i)*timesteps:(s+i+1)*timesteps] for i in range(m)])
        s += m
    return p, macro_steps, macro_list, log_list

def dune_preprocess(data, exceptions = []):
    data_new = []
    for line in data:
        ## skip empty lines
        if line == '':
            continue
        ## skip lines that contain non-numbers
        for l in line.split():
            try:
                float(l)
            except:
                if l in exceptions:
                    continue
                else:
                    break
        else:
            data_new.append(line)
    return data_new
        
def process_output(path):
    dic = {}
    with open(path + 'content.txt', 'r') as myfile:
        dic = json.load(myfile)
    files = dic['files']
    
    parameters = {}
    with open(path + 'parameters.txt', 'r') as myfile:
        parameters = json.load(myfile)
        
    results, results_comm = {}, {}
    run_names = parameters['run_names'] + ['REF']
    for r in run_names:
        results[r] = []
        if 'NEW' not in r:
            results_comm[r] = []
    
    for file in files:
        res = {}
        f = cut_path(file[:-4])
        splitted = f.split('_')
        rm = '_'.join(splitted[:-1]) # in case name contains _
        tol = float(splitted[-1])
        
        if tol == parameters['tolerances'][-1]:
            rm = 'REF'
        
        data = 0
        with open(file, 'r') as myfile:
            data = myfile.read().split('\n')[:-1]
            
        data = dune_preprocess(data) ## dune preprocessing, i.e., filtering out string output from alugrid
            
        res = {'tol': tol, 'iter': [], 'times': [], 'sol1': [], 'sol2': []}
        res_comm_log = {'tol': tol, 'msteps': [], 'mlog': [], 'log0': [], 'log1': []}
        for i in range(min(int(parameters['times'])*4, int(len(data)/2))):
            if data[2*i][:3] == 'LOG':
                if rm in ['JAC', 'GS', 'REF']: # no logs for dterministic processes
                    continue
                p, macro_steps, macro_list, log_list = process_log(data[2*i])
                assert p == 0
                
                res_comm_log['msteps'].append(macro_steps)
                res_comm_log['mlog'].append(macro_list)
                res_comm_log['log0'].append(log_list)
                
                p, macro_steps, macro_list, log_list = process_log(data[2*i + 1])
                assert p == 1
                res_comm_log['log1'].append(log_list)
            else:
                it, time, vec0, vec1 = split_and_float(data[2*i], data[2*i+1])
                res['iter'].append(it)
                res['times'].append(time)
                if rm not in ['JAC', 'GS', 'REF'] or i == 0: # skip storing of vector, if result is deterministic anyways
                    res['sol1'].append(vec0)
                    res['sol2'].append(vec1)
                else:
                    continue
        results[rm].append(res)
        if 'NEW' not in rm:
            results_comm[rm].append(res_comm_log)
    with open(path + 'results.txt', 'w') as myfile:
        myfile.write(json.dumps(results, sort_keys = True, indent = 4))
    if 'commlog' in parameters.keys() and parameters['commlog']:
        with open(path + 'comm.txt', 'w') as myfile:
            myfile.write(json.dumps(results_comm))
            
def process_output_errlog(path, which = 2):
    dic = {}
    with open(path + 'content.txt', 'r') as myfile:
        dic = json.load(myfile)
    files = dic['files']
    
    parameters = {}
    with open(path + 'parameters.txt', 'r') as myfile:
        parameters = json.load(myfile)
        
    results, results_comm = {}, {}
    run_names = parameters['run_names'] + ['REF']
    for r in run_names:
        results[r] = []
        if 'NEW' not in r:
            results_comm[r] = []
            
    for file in files:
        res = {}
        f = cut_path(file[:-4])
        splitted = f.split('_')
        rm = '_'.join(splitted[:-1]) # in case name contains _
        print(rm)
        tol = float(splitted[-1])
        
        data = 0
        with open(file, 'r') as myfile:
            data = myfile.read().split('\n')[:-1]
            
        data = dune_preprocess(data) ## dune preprocessing, i.e., filtering out string output from alugrid
        
        sol1_list, sol2_list, times_list, updates_list = [], [], [], []
        ## store results in lists of lists
        i = 0
        while True: ## this loop is to handle multiple outputs in one file
            sols1, sols2, times, updates = [], [], [], []
            ## process lines
            while int(data[i][:2]) != -1: ## until line with updates is found
                it, time, vec1, vec2 = split_and_float(data[i], data[i+1])
                sols1.append(vec1)
                sols2.append(vec2)
                times.append(time)
                i += 2
            up = [float(a) for a in data[i].split()] ## process updates
            i+= 1
            rel_tol = up[1]
            
            ## store list of result in list
            sol1_list.append(sols1)
            sol2_list.append(sols2)
            times_list.append(times)
            updates_list.append(np.array(up[2:])/rel_tol)
            if i >= len(data): ## file done
                break

        ## iterate over all tolerances        
        for tol in parameters['tolerances'][:-1]:
            res = {'tol': tol, 'iter': [], 'times': [], 'sol1': [], 'sol2': []}
            ## iterate over each result in the results list
            for sols1, sols2, times, updates in zip(sol1_list, sol2_list, times_list, updates_list):
                idx = len(sols1) - 3 ## default, if not converged
                try:
                    ## find idx corresponding to first update meeting the tolerance requirement
                    idx = np.where(updates < tol)[0][0]
                except:
                    pass
                ## filter out actual solution
                res['iter'].append(float(idx + 2))
                res['times'].append(times[idx + 2])
                res['sol1'].append(sols1[idx + 2])
                res['sol2'].append(sols2[idx + 2])
            ## append to results data for producing plotting data
            results[rm].append(res)
        if rm == parameters['reference_sol']:
            res = {'tol': parameters['tolerances'][-1], 'iter': [it], 'times': [times[-1]], 'sol1': [sols1[-1]], 'sol2': [sols2[-1]]}
            results['REF'].append(res)

    with open(path + 'results.txt', 'w') as myfile:
        myfile.write(json.dumps(results, sort_keys = True, indent = 4))

def delete_raw_data(path):
    dic = {}
    with open(path + 'content.txt', 'r') as myfile:
        dic = json.load(myfile)
    files = dic['files']
    for i in files:
        os.remove(i)
    os.remove(path + 'content.txt')