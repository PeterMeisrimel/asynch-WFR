#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 16:43:51 2020

@author: Peter Meisrimel, Lund University
"""

import numpy as np
import pylab as pl
from heat_dune import Problem_heat
from scipy.interpolate import interp1d

pl.rcParams['lines.linewidth'] = 2
pl.rcParams['font.size'] = 14

def solve_exact(tf, Nx = None, Ny = None, dx = None, **kwargs):
    xa, xb = kwargs['xa'], kwargs['xb']
    if Nx is None: Nx = kwargs['gridsize'] + 1
    if Ny is None: Ny = kwargs['gridsize'] + 1
    if dx is None: dx = 1./(kwargs['gridsize'] + 1)
    
    prob = Problem_heat(**kwargs)
    prob.get_ex_sol(tf)
    sol = prob.get_sol((xb - xa)*Nx + 1, Ny + 1, dx, xx = xa)
    prob.get_ex_flux(tf)
    flux = prob.get_u_gamma(prob.unew)
    return sol, flux

def verify_time(prob = None, D = True, tf = 1, k = 10, kmin = 0, order = 1, savefig = None, **kwargs):
    xa, xb = kwargs['xa'], kwargs['xb']
    Nx = kwargs['gridsize'] + 1 ## nodes per unit length, - 1
    dx = 1./Nx
    
    ref, ref_f = solve_exact(tf, **kwargs)
    if D:
        ref = ref[:(Nx + 1)*(abs(xa)*Nx + 1)]
    else:
        ref = ref[-(Nx + 1)*(abs(xb)*Nx + 1):]
    
    L2_fac, L2_fac_intf = dx, np.sqrt(dx)
    L2_fac_intf = np.sqrt(L2_fac)
    p = prob(**kwargs, order = order)
    
    dts, errs, errs_f = [], [], []
    for n_steps in [2**i for i in range(kmin, k)]:
        p.reset()
        dts.append(tf/n_steps)
        p.solve(tf, n_steps)
        if D:
            u = p.get_sol(abs(xa)*Nx + 1, Nx + 1, dx, xa)
            errs.append(np.linalg.norm(u - ref, 2)*L2_fac)
            f = -p.get_flux()
            errs_f.append(np.linalg.norm(f - ref_f, 2)*L2_fac_intf)
        else:
            u = p.get_sol(Nx + 1, Nx + 1, dx, 0)
            errs.append(np.linalg.norm(u - ref, 2)*L2_fac)
    for i in range(k-1-kmin):
        print(np.log2(errs[i]/errs[i+1]))
    if D:
        print('flux')
        for i in range(k-1-kmin):
            print(np.log2(errs_f[i]/errs_f[i+1]))
     
    dts = np.array(dts)
    pl.figure()
    pl.loglog(dts, errs, label = 'err', marker = 'o')
    if D:
        pl.loglog(dts, errs_f, label = 'flux', marker = 'o')    
    pl.loglog(dts, dts, label = '1 st order', linestyle = '--')
    pl.loglog(dts, dts**2, label = '2 nd order', linestyle = '--')
    pl.legend()
    pl.grid(True, which = 'major')
    if savefig is not None:
        s = f'verify_time_ord_{order}_{D}.png'
        pl.savefig(savefig  + s, dpi = 100)
        
def verify_space(prob = None, D = True, tf = 1, N_steps = 100, k = 10, kmin = 0, order = 1, savefig = None, **kwargs):
    xa, xb = kwargs['xa'], kwargs['xb']
    
    dxs, errs, errs_f = [], [], []
    for gridsize in [2**i for i in range(kmin, k)]:
        Nx = gridsize + 1
        ref, ref_f = solve_exact(tf, gridsize = gridsize, **kwargs)
        
        if D:
            ref = ref[:(Nx + 1)*(abs(xa)*Nx + 1)]
        else:
            ref = ref[-(Nx + 1)*(abs(xb)*Nx + 1):]
            
        dxs.append(1./(gridsize + 1))
        L2_fac, L2_fac_intf = dxs[-1], np.sqrt(dxs[-1])
        
        p = prob(**kwargs, gridsize = gridsize, order = order)
        p.solve(tf, N_steps)
        if D:
            u = p.get_sol(abs(xa)*Nx + 1, Nx + 1, dxs[-1], xa)
            errs.append(np.linalg.norm(u - ref, 2)*L2_fac)
            f = -p.get_flux()
            errs_f.append(np.linalg.norm(f - ref_f, 2)*L2_fac_intf)
        else:
            u = p.get_sol(Nx + 1, Nx + 1, dxs[-1], 0)
            errs.append(np.linalg.norm(u - ref, 2)*L2_fac)
    for i in range(k-1-kmin):
        print(np.log2(errs[i]/errs[i+1]))
    if D:
        print('flux')
        for i in range(k-1-kmin):
            print(np.log2(errs_f[i]/errs_f[i+1]))
     
    dxs = np.array(dxs)
    pl.figure()
    pl.loglog(dxs, errs, label = 'err', marker = 'o')
    if D:
        pl.loglog(dxs, errs_f, label = 'flux', marker = 'o')    
    pl.loglog(dxs, dxs, label = '1 st order', linestyle = '--')
    pl.loglog(dxs, dxs**2, label = '2 nd order', linestyle = '--')
    pl.legend()
    pl.grid(True, which = 'major')
    if savefig is not None:
        s = f'verify_space_ord_{order}_{D}.png'
        pl.savefig(savefig  + s, dpi = 100)