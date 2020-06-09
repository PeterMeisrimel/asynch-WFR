#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 11:12:33 2020

@author: Peter Meisrimel, Lund University
"""

import numpy as np
import ufl
from numpy import sin, cos, exp, pi
import pylab as pl

from dune.grid import cartesianDomain
from dune.alugrid import aluConformGrid
from dune.fem.space import lagrange as solutionSpace
from dune.fem.utility import lineSample
from dune.ufl import expression2GF

pl.rcParams['lines.linewidth'] = 2
pl.rcParams['font.size'] = 14

def flux_test(alpha, lambda_diff, gs, t = 0, xa = -2, xb = 1):
    dx = 1./(gs + 1)
    length = xb - xa
    
    x = ufl.SpatialCoordinate(ufl.triangle)
    u0_expr = 500*ufl.sin(ufl.pi*(x[0] - xa)/length)*ufl.sin(ufl.pi*x[1])
    f0_expr = lambda_diff*500*ufl.cos(ufl.pi*(x[0] - xa)/length)*ufl.sin(ufl.pi*x[1])*ufl.pi/length
    
    domain = cartesianDomain([xa, 0.], [xb, 1.], [length*(gs + 1), gs + 1])
    mesh = aluConformGrid(domain)
    space = solutionSpace(mesh, order = 1)
    
    uref = space.interpolate(u0_expr, name = 'uref')
    uref_flux = space.interpolate(f0_expr, name = 'uref_flux')
    
    u_flux = expression2GF(uref.space.grid, ufl.grad(uref), uref.space.order)
    
    return (lineSample(uref, [0., 0.], [0., 1.], gs + 2)[1],
            lineSample(uref_flux, [0., 0.], [0., 1.], gs + 2)[1],
            lambda_diff*lineSample(u_flux, [0., 0.], [0., 1.], gs + 2)[1])
    
def py_test(alpha, lambda_diff, gs, t = 0, xa = -2, xb = 1):
    yy = np.linspace(0, 1, gs + 2) ## meshpoints for interface values
    length = xb - xa
    
    u0_expr = lambda x, y: 500*sin((x - xa)*pi/length)*sin(y*pi)
    f0_expr = lambda x, y: lambda_diff*500*cos((x - xa)*pi/length)*sin(y*pi)*pi/length
    
    return (np.array([u0_expr(0., y) for y in yy]),
            np.array([f0_expr(0., y) for y in yy]))
    
if __name__ == '__main__':
    pl.close('all')
    pp = {'alpha': 1, 'lambda_diff': 0.1, 't': 0, 'xa': -2, 'xb': 1}
    k = 8
    
    dxs, errs_u, errs_f1, errs_f2 = [], [], [], []
    for gs in [2**i for i in range(k)]:
        dxs.append(1./(gs + 1))
        L2_fac, L2_fac_intf = dxs[-1], np.sqrt(dxs[-1])
        
        u, u_flux, flux = flux_test(**pp, gs = gs)
#        print('u', u)
#        print('u_flux', u_flux)
#        print('flux', flux)
        uref, fref = py_test(**pp, gs = gs)
        
        errs_u.append(np.linalg.norm(u - uref, 2)*L2_fac_intf)
        errs_f1.append(np.linalg.norm(u_flux - fref, 2)*L2_fac_intf)
        errs_f2.append(np.linalg.norm(flux - fref, 2)*L2_fac_intf)
        
    dxs = np.array(dxs)
    pl.figure()
    pl.loglog(dxs, errs_u, label = 'intf val')
    pl.loglog(dxs, errs_f1, label = 'f1 val')
    pl.loglog(dxs, errs_f2, label = 'f2 val')
    
    pl.loglog(dxs, dxs, label = '1st', linestyle = '--')
    pl.loglog(dxs, dxs**2, label = '2nd', linestyle = ':')
    pl.legend()
    pl.grid(True, 'major')
    
    yys, sols = [], []
    for gs in [2**i for i in range(k)]:
        u, u_flux, flux = flux_test(**pp, gs = gs)
        sols.append(flux)
        yys.append(np.linspace(0, 1, gs + 2))
    
    pl.figure()
    pl.plot(yys[-1], u_flux, label = 'ref', color = 'k')
    for y, f in zip(yys, sols):
        pl.plot(y, f, label = str(len(y) - 2))
    pl.legend()