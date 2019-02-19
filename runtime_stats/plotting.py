#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 16:09:23 2019

@author: Peter Mesirimel, Lund University
"""

from __future__ import division
import matplotlib.pylab as pl
import json

pl.close('all')
pl.rcParams['lines.linewidth'] = 3
pl.rcParams['lines.markersize'] = 17
fs = 32

import itertools
def reset_markers():
    return itertools.cycle('o^D*s'), itertools.cycle('gbrckmy')

def reset_fc():
    return itertools.cycle(['green', 'blue', 'red', 'cyan', 'black', 'magenta', 'yellow'])

import sys
path = sys.argv[1] # 0-th is filename

data = {}
with open(path + 'plotting_data.txt', 'r') as myfile:
    data = json.load(myfile)
    
parameters = {}
with open(path + 'parameters.txt', 'r') as myfile:
    parameters = json.load(myfile)
    
## create necessary directories
import os
plot_path = 'plots_' + path[path.find('/') + 1:-1]
dirs = [plot_path, plot_path + '/eps', plot_path + '/png']
for d in dirs:
    if not os.path.exists(d):
        os.makedirs(d)
        
from shutil import copy2
copy2(path + 'parameters.txt', plot_path + '/eps')
copy2(path + 'parameters.txt', plot_path + '/png')
        
fig, ax = pl.subplots(figsize = (12, 10))
markers, colors = reset_markers()

for k in data.keys():
    c = next(colors)
    m = next(markers)
    pl.semilogx(data[k]['tols'], data[k]['iters_med'], label = str(k), marker = m, color = c)
    if k == 'NEW':
        pl.semilogx(data[k]['tols'], data[k]['iters_min'], color = c, lw = 2)
        pl.semilogx(data[k]['tols'], data[k]['iters_max'], color = c, lw = 2)
        
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
    c = next(colors)
    m = next(markers)
    pl.semilogx(data[k]['tols'], data[k]['errors2_med'], label = str(k), marker = m, color = c)
    if k == 'NEW':
        pl.semilogx(data[k]['tols'], data[k]['errors2_min'], color = c, lw = 2)
        pl.semilogx(data[k]['tols'], data[k]['errors2_max'], color = c, lw = 2)
        
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
    c = next(colors)
    m = next(markers)
    pl.loglog(data[k]['tols'], data[k]['times_med'], label = str(k), marker = m, color = c)
    if k == 'NEW':
        pl.loglog(data[k]['tols'], data[k]['times_min'], color = c, lw = 2)
        pl.loglog(data[k]['tols'], data[k]['times_max'], color = c, lw = 2)
    
ax.set_title('Tolerance vs. Computational time', fontsize = fs + 2)
ax.set_xlabel('Waveform update tolerance', fontsize = fs)
ax.set_ylabel('Time (Wall Clock)', fontsize = fs)
ax.grid(b = True, which = 'major')
ax.tick_params(labelsize = 20)
ax.legend(fontsize = fs - 12, loc = 0)
fig.savefig(plot_path + '/eps/tolerance_vs_time.eps', dpi = 100)
fig.savefig(plot_path + '/png/tolerance_vs_time.png', dpi = 100)

###########################################################
fig, ax = pl.subplots(figsize = (12, 10))
markers, colors = reset_markers()
facecolor = reset_fc()

for k in data.keys():
    c = next(colors)
    m = next(markers)
    f = next(facecolor)
    ax.loglog(data[k]['times_med'], data[k]['errors2_med'], label = str(k), marker = m, color = c)
    for n, t in enumerate(data[k]['times']):
        val_t = data[k]['times_med'][n]
        min_t = data[k]['times_min'][n]
        max_t = data[k]['times_max'][n]
        
        val_e = data[k]['errors2_med'][n]
        min_e = data[k]['errors2_min'][n]
        max_e = data[k]['errors2_max'][n]
        ax.fill_between([min_t, val_t, max_t], [val_e, max_e, val_e], [val_e, min_e, val_e], facecolor=f, interpolate=True, alpha = 0.3)
    
ax.set_title('Computational time vs. Error', fontsize = fs + 2)
ax.set_xlabel('Time (Wall Clock)', fontsize = fs)
ax.set_ylabel('Error', fontsize = fs)
ax.grid(b = True, which = 'major')
ax.tick_params(labelsize = 20)
ax.legend(fontsize = fs-12, loc = 0)
fig.savefig(plot_path + '/eps/err_vs_time.eps', dpi = 100)
fig.savefig(plot_path + '/png/err_vs_time.png', dpi = 100)

###########################################################
fig, ax = pl.subplots(figsize = (12, 10))
markers, colors = reset_markers()
facecolor = reset_fc()

for k in data.keys():
    c = next(colors)
    m = next(markers)
    f = next(facecolor)
    pl.semilogx(data[k]['errors2_med'], data[k]['iters_med'], label = str(k), marker = m, color = c)
    for n, t in enumerate(data[k]['iters_med']):
        val_e = data[k]['errors2_med'][n]
        min_e = data[k]['errors2_min'][n]
        max_e = data[k]['errors2_max'][n]
        
        val_i = data[k]['iters_med'][n]
        min_i = data[k]['iters_min'][n]
        max_i = data[k]['iters_max'][n]
        ax.fill_between([min_e, val_e, max_e], [val_i, max_i, val_i], [val_i, min_i, val_i], facecolor = f, interpolate=True, alpha = 0.3)
        
ax.set_title('Error vs. Iterations', fontsize = fs + 2)
ax.set_xlabel('Error', fontsize = fs)
ax.set_ylabel('Iterations', fontsize = fs)
ax.grid(b = True, which = 'major')
ax.tick_params(labelsize = 20)
ax.legend(fontsize = fs-12, loc = 0)
fig.savefig(plot_path + '/eps/err_vs_iter.eps', dpi = 100)
fig.savefig(plot_path + '/png/err_vs_iter.png', dpi = 100)