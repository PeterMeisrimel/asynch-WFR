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
pl.rcParams['lines.linewidth'] = 2
pl.rcParams['font.size'] = 18
pl.rcParams['lines.markersize'] = 9

import itertools
def reset_markers():
    return itertools.cycle('o^D*s'), itertools.cycle('gbrckmy')

def reset_fc():
    return itertools.cycle(['green', 'blue', 'red', 'cyan', 'black', 'magenta', 'yellow'])

import sys
path = sys.argv[1] # 0-th is filename

# option to deactivate 
plot_variances = True
if len(sys.argv) == 3:
    if sys.argv[2] == 'False':
        plot_variances = False

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
        
fig, ax = pl.subplots()
markers, colors = reset_markers()

for num, k in enumerate(data.keys()):
    c = next(colors)
    m = next(markers)
    pl.semilogx(data[k]['tols'], data[k]['iters_med'], label = data[k]['label'], marker = m, color = c)
    if plot_variances and ('NEW' in k):
        pl.semilogx(data[k]['tols'], data[k]['iters_min'], color = c)
        pl.semilogx(data[k]['tols'], data[k]['iters_max'], color = c)
        
#ax.set_title('Tolerance vs. #Waveform iterations', fontsize = fs + 2)
ax.set_xlabel('TOL', labelpad = -20, position = (1.05, -1), fontsize = 20)
ax.set_ylabel('k', rotation = 0, labelpad = -40, position = (2., 1.05), fontsize = 20)
ax.grid(b = True, which = 'major')
ax.tick_params(labelsize = 20)
ax.legend(loc = 0)
fig.savefig(plot_path + '/eps/tolerance_vs_iters.eps', dpi = 100)
fig.savefig(plot_path + '/png/tolerance_vs_iters.png', dpi = 100)
###########################################################

fig, ax = pl.subplots()
markers, colors = reset_markers()

pl.loglog([float(x) for x in parameters['tolerances']],
          [float(x) for x in parameters['tolerances']],
          label = 'TOL', marker = None, color = 'k', linestyle = '--')
for num, k in enumerate(data.keys()):
    c = next(colors)
    m = next(markers)
    pl.semilogx(data[k]['tols'], data[k]['errors2_med'], label = data[k]['label'], marker = m, color = c)
    if plot_variances and ('NEW' in k):
        pl.semilogx(data[k]['tols'], data[k]['errors2_min'], color = c)
        pl.semilogx(data[k]['tols'], data[k]['errors2_max'], color = c)
        
#ax.set_title('Tolerance vs. Error', fontsize = fs + 2)
ax.set_xlabel('TOL', labelpad = -20, position = (1.08, -1), fontsize = 20)
ax.set_ylabel('Err', rotation = 0, labelpad = -50, position = (2., 1.05), fontsize = 20)
ax.grid(b = True, which = 'major')
ax.tick_params(labelsize = 20)
ax.legend(loc = 0)
fig.savefig(plot_path + '/eps/tolerance_vs_err.eps', dpi = 100)
fig.savefig(plot_path + '/png/tolerance_vs_err.png', dpi = 100)
###########################################################

fig, ax = pl.subplots()
markers, colors = reset_markers()

for num, k in enumerate(data.keys()):
    c = next(colors)
    m = next(markers)
    pl.loglog(data[k]['tols'], data[k]['times_med'], label = data[k]['label'], marker = m, color = c)
    if plot_variances and ('NEW' in k):
        pl.loglog(data[k]['tols'], data[k]['times_min'], color = c)
        pl.loglog(data[k]['tols'], data[k]['times_max'], color = c)
    
#ax.set_title('Tolerance vs. Computational time', fontsize = fs + 2)
ax.set_xlabel('TOL', labelpad = -20, position = (1.05, -1), fontsize = 20)
ax.set_ylabel('Time', rotation = 0, labelpad = -70, position = (2., 1.05), fontsize = 20)
ax.grid(b = True, which = 'major')
ax.tick_params(labelsize = 20)
ax.legend(loc = 0)
fig.savefig(plot_path + '/eps/tolerance_vs_time.eps', dpi = 100)
fig.savefig(plot_path + '/png/tolerance_vs_time.png', dpi = 100)

###########################################################
fig, ax = pl.subplots()
markers, colors = reset_markers()
facecolor = reset_fc()

for num, k in enumerate(data.keys()):
    c = next(colors)
    m = next(markers)
    f = next(facecolor)
    ax.loglog(data[k]['times_med'], data[k]['errors2_med'], label = data[k]['label'], marker = m, color = c)
    if plot_variances and ('NEW' in k):
        for n, t in enumerate(data[k]['times']):
            val_t = data[k]['times_med'][n]
            min_t = data[k]['times_min'][n]
            max_t = data[k]['times_max'][n]
            
            val_e = data[k]['errors2_med'][n]
            min_e = data[k]['errors2_min'][n]
            max_e = data[k]['errors2_max'][n]
            ax.fill_between([min_t, val_t, max_t], [val_e, max_e, val_e], [val_e, min_e, val_e], facecolor=f, interpolate=True, alpha = 0.3)
    
ax.set_xlabel('Time', labelpad = -20, position = (1.05, -1), fontsize = 20)
ax.set_ylabel('Err', rotation = 0, labelpad = -50, position = (2., 1.05), fontsize = 20)
ax.grid(b = True, which = 'major')
ax.tick_params(labelsize = 20)
ax.legend(loc = 0)
fig.savefig(plot_path + '/eps/err_vs_time.eps', dpi = 100)
fig.savefig(plot_path + '/png/err_vs_time.png', dpi = 100)

###########################################################
fig, ax = pl.subplots()
markers, colors = reset_markers()
facecolor = reset_fc()

for num, k in enumerate(data.keys()):
    c = next(colors)
    m = next(markers)
    f = next(facecolor)
    pl.semilogy(data[k]['iters_med'], data[k]['errors2_med'], label = data[k]['label'], marker = m, color = c)
    if plot_variances and ('NEW' in k):
        for n, t in enumerate(data[k]['iters_med']):
            val_e = data[k]['errors2_med'][n]
            min_e = data[k]['errors2_min'][n]
            max_e = data[k]['errors2_max'][n]
            
            val_i = data[k]['iters_med'][n]
            min_i = data[k]['iters_min'][n]
            max_i = data[k]['iters_max'][n]
            ax.fill_between([min_i, val_i, max_i], [val_e, max_e, val_e], [val_e, min_e, val_e], facecolor = f, interpolate=True, alpha = 0.3)
        
#ax.set_title('Error vs. Iterations')
ax.set_xlabel('k', labelpad = -20, position = (1.08, -1), fontsize = 20)
ax.set_ylabel('Err', rotation = 0, labelpad = -50, position = (2., 1.05), fontsize = 20)
ax.grid(b = True, which = 'major')
ax.tick_params(labelsize = 20)
ax.legend(loc = 0)
fig.savefig(plot_path + '/eps/err_vs_iter.eps', dpi = 100)
fig.savefig(plot_path + '/png/err_vs_iter.png', dpi = 100)