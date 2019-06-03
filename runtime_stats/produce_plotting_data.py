#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 13:20:46 2019

@author: Peter Mesirimel, Lund University
"""

from __future__ import division
import numpy as np
import json

def mean_var_med(f):
    return np.mean(f), np.var(f), np.median(f)

def produce_plotting_data(path):
    parameters = {}
    with open(path + 'parameters.txt', 'r') as myfile:
        parameters = json.load(myfile)
    
    data = {}
    with open(path + 'results.txt', 'r') as myfile:
        data = json.load(myfile)
        
    results = {}
    for r in parameters['run_names']:
        results[r] = {}
        
        ## tolerances
        tols = [i['tol'] for i in data[r]]
        results[r]['tols'] = tols
        
        ## times
        times, times_var, times_min, times_max, times_med = [], [], [], [], []
        for i in data[r]:
            f = i['times']
            m, v, med = mean_var_med(f)
            times.append(m)
            times_var.append(v)
            times_med.append(med)
            times_min.append(np.min(f))
            times_max.append(np.max(f))
        results[r]['times'] = times
        results[r]['times_var'] = times_var
        results[r]['times_med'] = times_med
        results[r]['times_min'] = times_min
        results[r]['times_max'] = times_max
            
        ## iterations
        iters, iters_var, iters_min, iters_max, iters_med = [], [], [], [], []
        for i in data[r]:
            f = i['iter']
            m, v, med = mean_var_med(f)
            iters.append(m)
            iters_var.append(v)
            iters_med.append(med)
            iters_min.append(np.min(f))
            iters_max.append(np.max(f))
        results[r]['iters'] = iters
        results[r]['iters_var'] = iters_var
        results[r]['iters_med'] = iters_med
        results[r]['iters_min'] = iters_min
        results[r]['iters_max'] = iters_max
            
        ## errors, 1
        ref = np.array(data['REF'][0]['sol1'][0])
        
        errors1, errors1_var, errors1_min, errors1_max, errors1_med = [], [], [], [], []
        for i in data[r]:
            errs = []
            f = i['sol1']
            for sol in f:
                errs.append(np.linalg.norm(np.array(sol) - ref, 2))
            m, v, med = mean_var_med(errs)
            errors1.append(m)
            errors1_var.append(v)
            errors1_med.append(med)
            errors1_min.append(np.min(errs))
            errors1_max.append(np.max(errs))
        results[r]['errors1'] = errors1
        results[r]['errors1_var'] = errors1_var
        results[r]['errors1_med'] = errors1_med
        results[r]['errors1_min'] = errors1_min
        results[r]['errors1_max'] = errors1_max
            
        ## errors, 2
        ref = np.array(data['REF'][0]['sol2'][0])
        
        errors2, errors2_var, errors2_min, errors2_max, errors2_med = [], [], [], [], []
        for i in data[r]:
            errs = []
            f = i['sol2']
            for sol in f:
                errs.append(np.linalg.norm(np.array(sol) - ref, 2))
            m, v, med = mean_var_med(errs)
            errors2.append(m)
            errors2_var.append(v)
            errors2_med.append(med)
            errors2_min.append(np.min(errs))
            errors2_max.append(np.max(errs))
        results[r]['errors2'] = errors2
        results[r]['errors2_var'] = errors2_var
        results[r]['errors2_med'] = errors2_med
        results[r]['errors2_min'] = errors2_min
        results[r]['errors2_max'] = errors2_max
            
        ## errors combined
        ref = np.array(data['REF'][0]['sol1'][0] + data['REF'][0]['sol2'][0])
        
        errors_comb, errors_comb_var, errors_comb_min, errors_comb_max, errors_comb_med = [], [], [], [], []
        for i in data[r]:
            errs = []
            f = [i['sol1'][j] + i['sol2'][j] for j in range(len(i['sol1']))]
            for sol in f:
                errs.append(np.linalg.norm(np.array(sol) - ref, 2))
            m, v, med = mean_var_med(errs)
            errors_comb.append(m)
            errors_comb_var.append(v)
            errors_comb_med.append(med)
            errors_comb_min.append(np.min(errs))
            errors_comb_max.append(np.max(errs))
        results[r]['errors'] = errors_comb
        results[r]['errors_var'] = errors_comb_var
        results[r]['errors_med'] = errors_comb_med
        results[r]['errors_min'] = errors_comb_min
        results[r]['errors_max'] = errors_comb_max
        
    with open(path + 'plotting_data.txt', 'w') as myfile:
        myfile.write(json.dumps(results, sort_keys = True))