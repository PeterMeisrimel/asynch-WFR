#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 11:10:56 2020

@author: Peter Meisrimel, with lots of help from Robert Klöfkorn
"""

import numpy as np
import pylab as pl
from scipy.interpolate import interp1d

# doesn't find dune.fem on second run?
import ufl
from dune.fem.space import lagrange as solutionSpace
from dune.ufl import DirichletBC, Constant
from dune.fem.scheme import galerkin as solutionScheme
from dune.fem.utility import lineSample
from dune.grid import cartesianDomain
from dune.alugrid import aluConformGrid

## standard heat equation: alpha u_t + lamda_diff \Delta u = 0
## assume interface to be at zero, spatial domain marked by xa and xb
class Problem_heat:
    eps = 1e-8
    ## gridsize = number of internal variables
    def __init__(self, gridsize, alpha, lambda_diff, order = 2, xa = -1, xb = 1, mono = True):
        ## n = gridsize = number of internal unknowns
        ## => n + 2 nodes per unit length
        ## => n + 1 cells per unit length
        self.NN = gridsize + 2
        self.yy = np.linspace(0, 1, self.NN) ## meshpoints for interface values
        self.dx = 1./(gridsize + 1)
        
        self.a, self.lam = Constant(alpha), Constant(lambda_diff)
        self.dt = Constant(0.)
        
        if mono:
            self.xa, self.xb = xa, xb
            self.len = self.xb - self.xa ## length of x domain
            self.domain = cartesianDomain([self.xa, 0.], [self.xb, 1.], [self.len*(gridsize + 1), gridsize + 1])
        else:
            self.len = self.xb - self.xa ## length of x domain
            self.domain = cartesianDomain([xa, 0.], [xb, 1.], [self.len*(gridsize + 1), gridsize + 1])
        
        self.mesh = aluConformGrid(self.domain)
        self.space = solutionSpace(self.mesh, order = 1)
        
        self.x = ufl.SpatialCoordinate(ufl.triangle)
        self.u0 = 500*ufl.sin(ufl.pi*(self.x[0] - self.xa)/self.len)*ufl.sin(ufl.pi*self.x[1])
        self.f0_expr = lambda_diff*500*ufl.cos(ufl.pi*(self.x[0] - self.xa)/self.len)*ufl.sin(ufl.pi*self.x[1])*ufl.pi/self.len
#        exact_sol = ufl.exp(-(5*lambda_diff*ufl.pi**2*tend)/(4*alpha))*u0

        self.u = ufl.TrialFunction(self.space)
        self.v = ufl.TestFunction(self.space)
        self.uold = self.space.interpolate(self.u0, name = 'uold')
        self.unew = self.space.interpolate(self.u0, name = 'unew')
        self.u_checkpoint = self.uold.copy(name = 'u_checkpoint')
        
        if order == 1:
            self.A = self.a*self.u*self.v*ufl.dx + self.dt*self.lam*ufl.dot(ufl.grad(self.u), ufl.grad(self.v))*ufl.dx
            self.b = self.a*self.uold*self.v*ufl.dx
        elif order == 2:
            ## weak form, Crank-Nicolson method
            self.A = self.a*self.u*self.v*ufl.dx + 0.5*self.dt*self.lam*ufl.dot(ufl.grad(self.u), ufl.grad(self.v))*ufl.dx
            self.b = self.a*self.uold*self.v*ufl.dx - 0.5*self.dt*self.lam*ufl.dot(ufl.grad(self.uold), ufl.grad(self.v))*ufl.dx
        
        self.bc_bottom = DirichletBC(self.space, Constant(0.), self.x[1] < self.eps)
        self.bc_top = DirichletBC(self.space, Constant(0.), self.x[1] > 1 - self.eps)
        self.bcs = [self.bc_bottom, self.bc_top]
        #### add manually for D, resp N
        if mono:
            self.left = DirichletBC(self.space, Constant(0.), self.x[0] < self.xa + self.eps)
            self.right = DirichletBC(self.space, Constant(0.), self.x[0] > self.xb - self.eps)
            self.bcs.append(self.left); self.bcs.append(self.right)
            
            self.scheme = solutionScheme([self.A == self.b, *self.bcs], solver = 'cg')

        self.create_checkpoint()
        self.reset()
        
        
        
        self.u_gamma_f = lambda x: lineSample(x, [0., 0.], [0., 1.], self.NN)[1]
        
    def reset(self):
        self.uold.interpolate(self.u_checkpoint)
        self.unew.interpolate(self.u_checkpoint)
        self.dt.value = 0.
    def create_checkpoint(self):
        self.u_checkpoint.interpolate(self.uold)
        
    def get_u_gamma(self, ugamma):
#        return lineSample(ugamma, [0., 0.], [0., 1.], self.NN)[1]
        return self.u_gamma_f(ugamma)
    
    def do_step(self, dt, ug):
        self.dt.value = dt
        
        intp = interp1d(self.yy, ug)
        self.u_gamma.interpolate(lambda x: intp(x[1]))
        
        self.scheme.solve(target = self.unew)
        flux = self.get_flux(dt, self.unew, self.uold)
        self.uold.assign(self.unew)
        return flux  
        
    def solve(self, tf, n_steps):
        self.dt.value = tf/n_steps
        for _ in range(n_steps):
            self.scheme.solve(target = self.unew)
            self.uold.assign(self.unew)
    
    def get_sol(self, Nx, Ny, dx, xx = 0, offset = 0):
        res = np.zeros(Nx*Ny)
        for i in range(offset,  Nx):
#            res[(i-offset)*Ny:(i-offset+1)*Ny] = lineSample(self.unew, [xx + i*dx, 0.], [xx + i*dx, 1.], self.NN)[1]
            res[(i-offset)*Ny:(i-offset+1)*Ny] = lineSample(self.unew, [xx + i*dx, 0.], [xx + i*dx, 1.], Ny)[1]
        return res
    
    def get_ex_sol(self, t = 1):
        t_fac = np.exp(-(self.len**2 + 1)/(self.len**2)*np.pi**2*t*float(self.lam)/float(self.a))
        self.u_ex = self.u0*t_fac
        self.unew.interpolate(self.u_ex)
        
    def get_ex_flux(self, t = 1):
        t_fac = np.exp(-(self.len**2 + 1)/(self.len**2)*np.pi**2*t*float(self.lam)/float(self.a))
        self.u_flux_ex = self.f0_expr*t_fac
        self.unew.interpolate(self.u_flux_ex)

def get_solve_WR(Problem_heat_D, Problem_heat_N):
    def solve_WR(tf, N1, N2, maxiter, TOL, th, Nx = None, Ny = None, dx = None, **kwargs):
        xa = kwargs['xa']
        xb = kwargs['xb']
        if Nx is None: Nx = kwargs['gridsize'] + 1 ## nodes per unit length, - 1
        if Ny is None: Ny = kwargs['gridsize'] + 1 ## nodes per unit length, - 1
        if dx is None: dx = 1./(kwargs['gridsize'] + 1)
        
        pD, pN = Problem_heat_D(**kwargs), Problem_heat_N(**kwargs)
        tt1, tt2 = np.linspace(0, tf, N1 + 1), np.linspace(0, tf, N2 + 1)
        dt1, dt2 = tf/N1, tf/N2
        norm = lambda x: np.linalg.norm(x, 2)*np.sqrt(pN.dx)
        
        ## run Waveform iteration 
        ug0, f0 = pN.get_u0(), pD.get_u0()
        ugold = [np.copy(ug0) for i in range(N2 + 1)]
        ugnew = [np.copy(ug0) for i in range(N2 + 1)]
        flux_WF = [np.copy(f0) for i in range(N1 + 1)]
        
        pD.create_checkpoint(); pN.create_checkpoint()
        updates, rel_tol_fac = [], norm(ug0)
        for j in range(maxiter):
            ugold = [np.copy(uu) for uu in ugnew]
            pD.reset(); pN.reset()
            
            ## Dirichlet
            ug_WF_f = interp1d(tt2, ugold, kind = 'linear', axis = 0, fill_value = 'extrapolate')
            for i, t in enumerate(tt1[:-1]):
                flux_WF[i+1] = pD.do_step(dt1, ug_WF_f(t + dt1))
                
            ## Neumann
            flux_WF_f = interp1d(tt1, flux_WF, kind = 'linear', axis = 0, fill_value = 'extrapolate')
            for i, t in enumerate(tt2[:-1]):
                ugnew[i+1] = pN.do_step(dt2, flux_WF_f(t), flux_WF_f(t + dt2))
                
            ## relaxation
            tmp = np.copy(ugold[-1])
            for i in range(N2 + 1):
                ugnew[i] = (1-th)*ugold[i] + th*ugnew[i]
            
            updates.append(norm(ugnew[-1] - tmp))
            print(j, ' update = ', updates[-1])
            if updates[-1]/rel_tol_fac < TOL: # STOPPING CRITERIA FOR FIXED POINT ITERATION
                break
        return (pD.get_sol(Nx * abs(xa), Ny + 1, dx, xx = xa),
                pN.get_sol(Nx * abs(xb), Ny + 1, dx, xx = 0, offset = 1),
                ugnew[-1], updates, j+1)
    return solve_WR

if __name__ == '__main__':
    from verification import get_parameters, verify_space_error, verify_mono_time
    savefig = 'mono_'
    ## verify space order of monolithic solution
    verify_space_error(1., k = 6, order = 2, **get_parameters(), xa = -2, xb = 1, savefig = savefig)
    
    ## verify time order with itself
    pp = {'tf': 1., **get_parameters(), 'gridsize': 64, 'xa': -2, 'xb': 1}
    verify_mono_time(k = 7, order = 1, **pp, savefig = savefig)
    verify_mono_time(k = 5, order = 2, **pp, savefig = savefig)