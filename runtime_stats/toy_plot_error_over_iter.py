#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 13:31:49 2019

@author: Peter Mesirimel, Lund University
"""

from __future__ import division
from postprocessing_cpp_output import cut_path
from run_tolerances import run_tolerances
import json
import numpy as np
import pylab as pl
pl.close('all')

def split_and_float(f):
    a = [float(i) for i in f.split(' ')[:-1]]
    return a[0], a[1], a[2:]

tol = 1e-4
gridsize = 32
name = 'toy_test'
folder = 'toy'
exe = 'TOY'

parameters = {'timesteps' : 100, 'macrosteps': 1, 'maxiter': 200, 'errlog': 1, 'nsteps_conv_check': 3}
path = run_tolerances(folder, exe, name, times = 1, tolerances = [tol, tol/100], parameters = parameters,
                      runmodes = ['GS1', 'GS2', 'JAC', 'NEW'], run_names = ['GS_1', 'GS_2', 'JAC', 'NEW'],
                      ref_run_name = 'GS_1',
                      labels = ['GS 1', 'GS 2', 'JAC', 'NEW'])
""" process data """
results, dic = {}, {}
with open(path + 'content.txt', 'r') as myfile:
    dic = json.load(myfile)
files = dic['files']
    
with open(path + 'parameters.txt', 'r') as myfile:
    parameters = json.load(myfile)
        
run_names = parameters['run_names'] + ['REF']
for r in run_names:
    results[r] = []
    
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
        
    res = {'tol': tol, 'iters': 0, 'sol1': [], 'sol2': []}
    for line in data:
        p, it, vec = split_and_float(line)
        if p == 0:
            res['sol1'].append(vec)
        elif p == 1:
            res['sol2'].append(vec)
    res['iters'] = it
    if rm == 'REF':
        res['sol1'] = res['sol1'][-1]
        res['sol2'] = res['sol2'][-1]
    results[rm] = res
with open(path + 'results.txt', 'w') as myfile:
    myfile.write(json.dumps(results, sort_keys = True, indent = 4))
    
""" produce plotting data """
data = {}
with open(path + 'results.txt', 'r') as myfile:
    data = json.load(myfile)
    
ref_sol = np.array(data['REF']['sol2'])
results = {}
for r in parameters['run_names']:
    results[r] = {}
    results[r]['errs2'] = [np.linalg.norm(ref_sol - np.array(ii), 2) for ii in data[r]['sol2']]
    results[r]['iters'] = data[r]['iters']
            
with open(path + 'plotting_data.txt', 'w') as myfile:
    myfile.write(json.dumps(results, sort_keys = True))
    
""" plotting """
pl.close('all')
pl.rcParams['lines.linewidth'] = 2
pl.rcParams['font.size'] = 18
pl.rcParams['lines.markersize'] = 9

import itertools
def reset_markers():
    return itertools.cycle('o^D*s'), itertools.cycle('gbrckmy')

def reset_fc():
    return itertools.cycle(['green', 'blue', 'red', 'cyan', 'black', 'magenta', 'yellow'])

data = {}
with open(path + 'plotting_data.txt', 'r') as myfile:
    data = json.load(myfile)
label = parameters['labels']

## create necessary directories
import os
plot_path = 'plots_iterations_' + path[path.find('/') + 1:-1]
dirs = [plot_path, plot_path + '/eps', plot_path + '/png']
for d in dirs:
    if not os.path.exists(d):
        os.makedirs(d)
        
from shutil import copy2
copy2(path + 'parameters.txt', plot_path + '/eps')
copy2(path + 'parameters.txt', plot_path + '/png')
        
fig, ax = pl.subplots()
markers, colors = reset_markers()

for i, r in enumerate(parameters['run_names']):
    c = next(colors)
    m = next(markers)
    ax.semilogy(range(int(data[r]['iters'])+1), data[r]['errs2'], label = label[i], marker = m, color = c)
    
ax.set_xlabel('k', labelpad = -20, position = (1.08, -1), fontsize = 20)
ax.set_ylabel('Err', rotation = 0, labelpad = -20, position = (2., 1.05), fontsize = 20)
ax.grid(b = True, which = 'major')
ax.legend(loc = 0)
fig.savefig(plot_path + '/eps/error_vs_iters.eps', dpi = 100)
fig.savefig(plot_path + '/png/error_vs_iters.png', dpi = 100)