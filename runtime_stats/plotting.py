#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 16:09:23 2019

@author: Peter Mesirimel, Lund University
"""

from __future__ import division
import pylab as pl
import numpy as np
import scipy as sp
pl.close('all')
pl.rcParams['lines.linewidth'] = 3
pl.rcParams['lines.markersize'] = 17
fs = 32

import itertools
def reset_markers():
    return itertools.cycle('o^D*s'), itertools.cycle('gbrckmy')

import sys
input_file = sys.argv[1] # 0-th is filename

import json
data = {}
with open(input_file, 'r') as myfile:
    data = json.load(myfile)
## create necessary directories
import os
plot_path = 'plots_' + input_file
dirs = [plot_path, plot_path + '/eps', plot_path + '/png']
for d in dirs:
    if not os.path.exists(d):
        os.makedirs(d)

# get reference solution from smallest GS tolerance        
ref_sol = np.array(data['GS']['sol'][-1])
# cut away 
data['GS']['sol']  = data['GS']['sol'][:-1]
data['GS']['tols'] = data['GS']['tols'][:-1]
data['GS']['iter'] = data['GS']['iter'][:-1]
data['GS']['t']    = data['GS']['t'][:-1]
data['NEW']['tols'] = data['NEW']['tols'][:-1]
if 'JAC' in data.keys():
    data['JAC']['tols'] = data['JAC']['tols'][:-1]

# calculate errors
for k in data.keys():
    data[k]['errs'] = []
    for sol in data[k]['sol']:
        data[k]['errs'].append(np.linalg.norm(np.array(sol) - ref_sol, 2))
        
fig, ax = pl.subplots(figsize = (12, 10))
markers, colors = reset_markers()

for k in data.keys():
    pl.semilogx([float(x) for x in data[k]['tols']],
                [float(x) for x in data[k]['iter']],
                label = str(k), marker = next(markers), color = next(colors))
ax.set_title('Tolerance vs. #Waveform iterations', fontsize = fs + 2)
ax.set_xlabel('Waveform update tolerance', fontsize = fs)
ax.set_ylabel('Iterations', fontsize = fs)
ax.grid(b = True, which = 'major')
ax.tick_params(labelsize = 20)
ax.legend(fontsize = fs - 12, loc = 0)
fig.savefig(plot_path + '/eps/tolerance_vs_iters.eps', dpi = 100)
fig.savefig(plot_path + '/png/tolerance_vs_iters.png', dpi = 100)
###########################################################
fig, ax = pl.subplots(figsize = (12, 10))
markers, colors = reset_markers()

pl.loglog([float(x) for x in data['GS']['tols']],
          [float(x) for x in data['GS']['tols']],
          label = 'tol', marker = None, color = 'k', linestyle = '--')
for k in data.keys():
    pl.loglog([float(x) for x in data[k]['tols']],
               data[k]['errs'],
               label = str(k), marker = next(markers), color = next(colors))
ax.set_title('Tolerance vs. Error', fontsize = fs + 2)
ax.set_xlabel('Waveform update tolerance', fontsize = fs)
ax.set_ylabel('Err', fontsize = fs)
ax.grid(b = True, which = 'major')
ax.tick_params(labelsize = 20)
ax.legend(fontsize = fs-12, loc = 0)
fig.savefig(plot_path + '/eps/tolerance_vs_err.eps', dpi = 100)
fig.savefig(plot_path + '/png/tolerance_vs_err.png', dpi = 100)
###########################################################

fig, ax = pl.subplots(figsize = (12, 10))
markers, colors = reset_markers()

for k in data.keys():
    pl.loglog(data[k]['errs'],
              [float(x) for x in data[k]['t']],
              label = str(k), marker = next(markers), color = next(colors))
ax.set_title('Error vs. Computational time', fontsize = fs + 2)
ax.set_xlabel('Error', fontsize = fs)
ax.set_ylabel('Time', fontsize = fs)
#ax.set_xlim(left = 1e-16) # to-do temporary fix only
ax.grid(b = True, which = 'major')
ax.tick_params(labelsize = 20)
ax.legend(fontsize = fs-12, loc = 0)
fig.savefig(plot_path + '/eps/err_vs_time.eps', dpi = 100)
fig.savefig(plot_path + '/png/err_vs_time.png', dpi = 100)
###########################################################

fig, ax = pl.subplots(figsize = (12, 10))
markers, colors = reset_markers()

for k in data.keys():
    pl.loglog([float(x) for x in data[k]['tols']],
              [float(x) for x in data[k]['t']],
              label = str(k), marker = next(markers), color = next(colors))
ax.set_title('Tolerance vs. Computational time', fontsize = fs + 2)
ax.set_xlabel('Waveform update tolerance', fontsize = fs)
ax.set_ylabel('Time', fontsize = fs)
ax.grid(b = True, which = 'major')
ax.tick_params(labelsize = 20)
ax.legend(fontsize = fs-12, loc = 0)
fig.savefig(plot_path + '/eps/tolerance_vs_time.eps', dpi = 100)
fig.savefig(plot_path + '/png/tolerance_vs_time.png', dpi = 100)
###########################################################
"""
fig, ax = pl.subplots(figsize = (12, 10))
markers, colors = reset_markers()

for i in range(num):
    pl.loglog(data[i]['tols'], np.array(data[i]['times'])/np.array(data[i]['iters']), label = labels[i], marker = next(markers), color = next(colors))
for i in range(sd_idx):
    pl.semilogx(data_sd[i]['tols'], np.array(data_sd[i]['times'])/np.array(data_sd[i]['iters']), marker = None, color = 'k', linestyle = '--', linewidth = 1)
ax.set_title('Tolerance vs. Time/Iteration', fontsize = fs + 2)
ax.set_xlabel('Waveform update tolerance', fontsize = fs)
ax.set_ylabel('Time/Iter', fontsize = fs)
ax.grid(b = True, which = 'major')
ax.tick_params(labelsize = 20)
ax.set_xlim(left = max_tol, right = min_tol)
ax.legend(fontsize = fs-12, loc = 0)
fig.savefig(plot_path + '/eps/tolerance_vs_time_per_iteration.eps', dpi = 100)
fig.savefig(plot_path + '/png/tolerance_vs_time_per_iteration.png', dpi = 100)
"""