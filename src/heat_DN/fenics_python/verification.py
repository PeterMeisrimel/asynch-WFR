#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 16:43:51 2020

@author: Peter Meisrimel, Lund University
"""

import numpy as np
import pylab as pl
from heat_fenics import Problem_heat
from scipy.interpolate import interp1d

pl.rcParams['lines.linewidth'] = 2
pl.rcParams['font.size'] = 18
pl.rcParams['lines.markersize'] = 12

def get_parameters(which = 'test'):
    if which == 'test':
        return {'alpha': 1., 'lambda_diff': 0.1}
    else:
        return ValueError('invalid parameters')
    
def solve_monolithic(tf, N, Nx = None, Ny = None, dx = None, **kwargs):
    xa, xb = kwargs['xa'], kwargs['xb']
    if Nx is None: Nx = kwargs['gridsize'] + 1
    if Ny is None: Ny = kwargs['gridsize'] + 1
    if dx is None: dx = 1./(kwargs['gridsize'] + 1)
    
    prob = Problem_heat(**kwargs)
    prob.solve(tf, N)
    return prob.get_sol((xb - xa)*Nx + 1, Ny + 1, dx, xx = xa)

def solve_exact(tf, Nx = None, Ny = None, dx = None, **kwargs):
    xa, xb = kwargs['xa'], kwargs['xb']
    if Nx is None: Nx = kwargs['gridsize'] + 1
    if Ny is None: Ny = kwargs['gridsize'] + 1
    if dx is None: dx = 1./(kwargs['gridsize'] + 1)
    
    prob = Problem_heat(**kwargs)
    prob.get_ex_sol(tf)
    return prob.get_sol((xb - xa)*Nx + 1, Ny + 1, dx, xx = xa)

def verify_space_error(tf = 1, n_min = 2, N_steps = 100, k = 8, savefig = None, **kwargs):
    ## verifiy order of space discretization by picking a sufficiently fine grid in time and an exact solution
    errs = []
    
    n_list = np.array([n_min*(2**i) for i in range(k)])
    for n in n_list:
        # L2 factor for 2D case, inner points
        L2_fac = 1/(n + 1)
        ref = solve_exact(tf, gridsize = n, **kwargs)
        sol = solve_monolithic(tf, N_steps, gridsize = n, **kwargs)
        errs.append(np.linalg.norm(ref - sol, 2)*L2_fac)
    for i in range(k-1):
        print(np.log2(errs[i]/errs[i+1]))
    
    n_list = np.array(n_list)
    pl.figure()
    pl.loglog(n_list, errs, label = 'error', marker = 'o')
    pl.loglog(n_list, 1/n_list, label = '1st', linestyle = '--')
    pl.loglog(n_list, 1/(n_list**2), label = '2nd', linestyle = '--')
    pl.legend()
    pl.xlabel('gridsize'); pl.ylabel('err')
    pl.grid(True, which = 'major')
    pl.title('error in space')
    if savefig is not None:
        s = f'err_space_steps_{N_steps}.png'
        pl.savefig(savefig + s, dpi = 100)
        
def verify_mono_time(tf = 1, k = 5, savefig = None, **kwargs):
    ## verify time-integration order of monolithic solution with itself
    ref = solve_monolithic(tf, 2**(k+1), **kwargs)

    L2_fac = 1./(kwargs['gridsize'] + 1)
    errs = []
    steps = [2**i for i in range(k)]
    for s in steps:
        sol = solve_monolithic(tf, s, **kwargs)
        errs.append(np.linalg.norm(sol - ref, 2)*L2_fac)
    for i in range(k-1):
        print(np.log2(errs[i]/errs[i+1]))
        
    pl.figure()
    dts = np.array([tf/s for s in steps])
    pl.loglog(dts, errs, label = 'err', marker = 'o')
    pl.loglog(dts, dts, label = '1 st order', linestyle = '--')
    pl.loglog(dts, dts**2, label = '2 nd order', linestyle = '--')
    pl.legend()
    pl.xlabel('dt'); pl.ylabel('Err')
    pl.grid(True, which = 'major')
    pl.title('time-integration error, monolithic')
    if savefig is not None:
        order = kwargs['order']
        s = f'mono_time_ord_{order}.png'
        pl.savefig(savefig + s, dpi = 100)
        
def verify_with_monolithic(solve_WR = None, tf = 1, N_steps = 20, k = 8, theta = 1, maxiter = 100, savefig = None, **kwargs):
    ## verify convergence of coupling scheme for decreasing tolerances with monolithic solution, for fixed delta t
    ref = solve_monolithic(tf, N_steps, **kwargs)
    
    L2_fac = 1./(kwargs['gridsize'] + 1)
    errs = []
    tols = np.array([10**(-i) for i in range(k)])
    for tol in tols:
        u1, u2, ug, _, _ = solve_WR(tf, N_steps, N_steps, maxiter, tol, theta, **kwargs)
        errs.append(np.linalg.norm(ref - np.hstack((u1, ug, u2)), 2)*L2_fac)
    for i in range(k-1):
        print(np.log10(errs[i]/errs[i+1]))
        
    pl.figure()
    pl.loglog(tols, errs, label = 'err', marker = 'o')
    pl.loglog(tols, tols, label = '1 st order', linestyle = '--')
    pl.loglog(tols, tols**2, label = '2 nd order', linestyle = '--')
    pl.title('Verification with monolithic solution')
    pl.legend()
    pl.xlabel('TOL'); pl.ylabel('Err')
    pl.grid(True, which = 'major')
    if savefig is not None:
        order = kwargs['order']
        s = f'verify_mono_time_steps_{N_steps}_ord_{order}.png'
        pl.savefig(savefig + s, dpi = 100)
        
def verify_comb_error(solve_WR = None, tf = 1, k = 10, kmin = 0, theta = 1, order = 1, savefig = None, maxiter = 100, **kwargs):
    ## sum of splitting and time-integration order
    ## calculating time-integration error( + splitting error) for a small tolerance
    ref = solve_monolithic(tf, 2**(k+1), **kwargs)
    
    L2_fac = 1/(kwargs['gridsize'] + 1)
    
    dts, errs = [], []
    for n_steps in [2**i for i in range(kmin, k)]:
        dts.append(tf/n_steps)
        u1, u2, ug, _, _ = solve_WR(tf, n_steps, n_steps, maxiter, 1e-10, theta, order = order, **kwargs)
        errs.append(np.linalg.norm(ref - np.hstack((u1, ug, u2)), 2)*L2_fac)
    for i in range(k-1-kmin):
        print(np.log2(errs[i]/errs[i+1]))
     
    dts = np.array(dts)
    pl.figure()
    pl.loglog(dts, errs, label = 'err', marker = 'o')
    pl.loglog(dts, dts, label = '1 st order', linestyle = '--')
    pl.loglog(dts, dts**2, label = '2 nd order', linestyle = '--')
    pl.legend()
    pl.title('Splitting + time int error test')
    pl.xlabel('dt'); pl.ylabel('Err')
    pl.grid(True, which = 'major')
    if savefig is not None:
        s = f'verify_comb_error_ord_{order}.png'
        pl.savefig(savefig  + s, dpi = 100)
        
def verify_comb_error_space(solve_WR = None, tf = 1, N_steps = 100, k = 10, kmin = 0, theta = 1, order = 1, savefig = None, maxiter = 100, TOL = 1e-10, **kwargs):
    # for each dx calculate sum of splitting and time-integration error.
    dxs, errs = [], []
    for gridsize in [2**i for i in range(kmin, k)]:
        ref = solve_monolithic(tf, N_steps, gridsize = gridsize, order = order, **kwargs)
        dxs.append(1./(gridsize + 1))
        L2_fac = dxs[-1]
        u1, u2, ug, _, _ = solve_WR(tf, N_steps, N_steps, maxiter, TOL, theta, gridsize = gridsize, order = order, **kwargs)
        errs.append(np.linalg.norm(ref - np.hstack((u1, ug, u2)), 2)*L2_fac)
    for i in range(k-1-kmin):
        print(np.log2(errs[i]/errs[i+1]))
     
    dxs = np.array(dxs)
    pl.figure()
    pl.loglog(dxs, errs, label = 'err', marker = 'o')
    pl.loglog(dxs, dxs, label = '1 st order', linestyle = '--')
    pl.loglog(dxs, dxs**2, label = '2 nd order', linestyle = '--')
    pl.legend()
    pl.xlabel('dx'); pl.ylabel('Err')
    pl.grid(True, which = 'major')
    if savefig is not None:
        s = f'verify_comb_error_space_ord_{order}.png'
        pl.savefig(savefig  + s, dpi = 100)
        
def verify_self_time(solve_WR = None, tf = 1, k = 10, kmin = 0, theta = 1, order = 1, savefig = None, maxiter = 100, TOL = 1e-10, **kwargs):
    ## sum of splitting and time-integration order
    ## calculating time-integration error( + splitting error) for a small tolerance
    refu1, refu2, refug, _, _ = solve_WR(tf, 2**(k+1), 2**(k+1), maxiter, TOL, theta, order = order, **kwargs)
    ref = np.hstack((refu1, refug, refu2))
    
    L2_fac = 1/(kwargs['gridsize'] + 1)
    
    dts, errs, errs_ug = [], [], []
    for n_steps in [2**i for i in range(kmin, k)]:
        dts.append(tf/n_steps)
        u1, u2, ug, _, _ = solve_WR(tf, n_steps, n_steps, maxiter, TOL, theta, order = order, **kwargs)
        errs.append(np.linalg.norm(ref - np.hstack((u1, ug, u2)), 2)*L2_fac)
        errs_ug.append(np.linalg.norm(refug - ug, 2)*np.sqrt(L2_fac))
    print('err full')
    for i in range(k-1-kmin):
        print(np.log2(errs[i]/errs[i+1]))
        
    print('err ug')
    for i in range(k-1-kmin):
        print(np.log2(errs_ug[i]/errs_ug[i+1]))
     
    dts = np.array(dts)
    pl.figure()
    pl.loglog(dts, errs, label = 'err', marker = 'o')
    pl.loglog(dts, errs_ug, label = 'err ug', marker = 'o')
    pl.loglog(dts, dts, label = '1 st order', linestyle = '--')
    pl.loglog(dts, dts**2, label = '2 nd order', linestyle = '--')
    pl.legend()
    pl.xlabel('dt'); pl.ylabel('Err')
    pl.grid(True, which = 'major')
    if savefig is not None:
        s = f'verify_self_time_ord_{order}.png'
        pl.savefig(savefig  + s, dpi = 100)

def plot_theta(solve_WR = None, tf = 1, N_steps = 50, order = 1, savefig = None, kmax = 6, res = 20, TOL = 1e-6, **kwargs):
    thetas = np.linspace(0, 1, res + 1)[1:]
    rates = []
    for th in thetas:
        print('theta', th)
        _, _, _, up, _ = solve_WR(tf, N_steps, N_steps, kmax, TOL, th, Nx = 1, Ny = 1, dx = 1, order = order, **kwargs)
        x = np.array(up[:-1])
        rates.append(np.mean(x[1:]/x[:-1]))
        
    pl.figure()
    pl.plot(thetas, rates, label = 'rate')
    pl.legend()
    if savefig is not None:
        s = f'rates_ord_{order}.png'
        pl.savefig(savefig  + s, dpi = 100)