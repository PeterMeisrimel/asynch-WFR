#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 11:12:33 2020

@author: Peter Meisrimel, Lund University
"""

import numpy as np
import dolfin as dol
from numpy import sin, cos, exp, pi
import pylab as pl

pl.rcParams['lines.linewidth'] = 2
pl.rcParams['font.size'] = 14

## testing spatial accuracy of flux computation

def flux_test(alpha, lambda_diff, gs, t = 0, xa = -2, xb = 1):
    yy = np.linspace(0, 1, gs + 2) ## meshpoints for interface values
    dx = 1./(gs + 1)
    length = xb - xa
    
    u0_expr = "500*sin((x[0] - {})*M_PI/{})*sin(M_PI*x[1])".format(xa, length)
    f0_expr = "{}*500*cos((x[0] - {})*M_PI/{})*sin(M_PI*x[1])*M_PI/{}".format(lambda_diff, xa, length, length)
    
    mesh = dol.RectangleMesh(dol.MPI.comm_self, dol.Point((xa, 0.)), dol.Point((0, 1.)), length*(gs + 1), gs + 1)
    
    space = dol.FunctionSpace(mesh, 'CG', 1)
    
    v = dol.TestFunction(space)
    
    uref = dol.Function(space)
    uref_flux = dol.Function(space)
    uref.interpolate(dol.Expression(u0_expr, degree = 2))
    uref_flux.interpolate(dol.Expression(f0_expr, degree = 2))
    
    u_flux = dol.Function(space)
    
    n = dol.Constant((1.0,0.0))
    
    u_flux_res = dol.assemble(lambda_diff*dol.dot(dol.grad(uref), n)*v*dol.ds)
    u_flux.vector().set_local(u_flux_res)
    
    
    return (np.array([uref(0., y) for y in yy]),
            np.array([uref_flux(0., y) for y in yy]),
            np.array([u_flux(0., y) for y in yy])/dx,
            uref_flux(0., 0.5),
            u_flux(0., 0.5)/dx)
        
def py_test(alpha, lambda_diff, gs, t = 0, xa = -2, xb = 1):
    yy = np.linspace(0, 1, gs + 2) ## meshpoints for interface values
    dx = 1./(gs + 1)
    length = xb - xa
    
    u0_expr = lambda x, y: 500*sin((x - xa)*pi/length)*sin(y*pi)
    f0_expr = lambda x, y: lambda_diff*500*cos((x - xa)*pi/length)*sin(y*pi)*pi/length
    
    return (np.array([u0_expr(0., y) for y in yy]),
            np.array([f0_expr(0., y) for y in yy]),
            lambda_diff*np.array([(u0_expr(-2*dx, y) - 4*u0_expr(-dx, y) + 3*u0_expr(0, y))/(2*dx) for y in yy]))
    
if __name__ == '__main__':
    pl.close('all')
    pp = {'alpha': 1, 'lambda_diff': 0.1, 't': 0, 'xa': -2, 'xb': 1}
    k = 8
    
    dxs, errs_u, errs_f1, errs_f2, errs_f_py, errs_p = [], [], [], [], [], []
    for gs in [2**i for i in range(k)]:
        dxs.append(1./(gs + 1))
        L2_fac, L2_fac_intf = dxs[-1], np.sqrt(dxs[-1])
        
        u, u_flux, flux, ref_p, flux_p = flux_test(**pp, gs = gs)
        uref, fref, f_approx_py = py_test(**pp, gs = gs)
        
        errs_u.append(np.linalg.norm(u - uref, 2)*L2_fac_intf)
        errs_f1.append(np.linalg.norm(u_flux - fref, 2)*L2_fac_intf)
        errs_f2.append(np.linalg.norm(flux - fref, 2)*L2_fac_intf)
        errs_f_py.append(np.linalg.norm(f_approx_py - fref, 2)*L2_fac_intf)
        errs_p.append(abs(ref_p - flux_p))
        
    print('flux_rate')
    for i in range(k-1):
        print(np.log2(errs_f2[i]/errs_f2[i+1]))
        
    print('point_rate')
    for i in range(k-1):
        print(np.log2(errs_p[i]/errs_p[i+1]))
        
    dxs = np.array(dxs)
    pl.figure()
#    pl.loglog(dxs, errs_u, label = 'intf val')
#    pl.loglog(dxs, errs_f1, label = 'f1 val')
    pl.loglog(dxs, errs_f2, label = 'flux grad')
    pl.loglog(dxs, errs_f_py, label = 'flux python fin diff')
    pl.loglog(dxs, errs_p, label = 'midpoint')
    
    pl.loglog(dxs, dxs, label = '1st', linestyle = '--')
    pl.loglog(dxs, dxs**2, label = '2nd', linestyle = ':')
    pl.legend()
    pl.title('spatial errors')
    pl.grid(True, 'major')
    
    yys, sols = [], []
    for gs in [2**i for i in range(k)]:
        u, u_flux, flux, _, _ = flux_test(**pp, gs = gs)
        sols.append(flux)
        yys.append(np.linspace(0, 1, gs + 2))
    
    pl.figure()
    pl.plot(yys[-1], u_flux, label = 'ref', color = 'k')
    for y, f in zip(yys, sols):
        pl.plot(y, f, label = str(len(y) - 2))
    pl.legend()
    pl.xlabel('y')
    pl.title('flux value visiualized for different dx')