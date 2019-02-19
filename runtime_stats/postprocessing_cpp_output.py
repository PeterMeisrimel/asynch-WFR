#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 13:21:10 2019

@author: Peter Mesirimel, Lund University
"""

from __future__ import division
import os
import json
                
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
        
def process_output(path):
    dic = {}
    with open(path + 'content.txt', 'r') as myfile:
        dic = json.load(myfile)
    files = dic['files']
    
    parameters = {}
    with open(path + 'parameters.txt', 'r') as myfile:
        parameters = json.load(myfile)
        
    results, results_comm = {}, {}
    runmodes = parameters['runmodes'] + ['REF']
    for r in runmodes:
        results[r] = []
        if r not in ['JAC', 'GS', 'REF']:
            results_comm[r] = []
    
    for file in files:
        res = {}
        f = cut_path(file[:-4])
        rm, tol = f.split('_')
        tol = float(tol)
        
        if tol == parameters['tolerances'][-1]:
            rm = 'REF'
        
        data = 0
        with open(file, 'r') as myfile:
            data = myfile.read().split('\n')[:-1]
        
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
        if rm not in ['JAC', 'GS', 'REF']:
            results_comm[rm].append(res_comm_log)
    with open(path + 'results.txt', 'w') as myfile:
        myfile.write(json.dumps(results))
    if 'logging' in parameters.keys() and parameters['logging']:
        with open(path + 'comm.txt', 'w') as myfile:
            myfile.write(json.dumps(results_comm))

def delete_raw_data(path):
    dic = {}
    with open(path + 'content.txt', 'r') as myfile:
        dic = json.load(myfile)
    files = dic['files']
    for i in files:
        os.remove(i)
    os.remove(path + 'content.txt')