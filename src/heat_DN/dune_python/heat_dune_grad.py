#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 22:31:31 2020

@author: Peter Meisrimel, with lots of help from Robert Klöfkorn
"""

import pylab as pl
import numpy as np
from scipy.interpolate import interp1d
from heat_dune import Problem_heat

import ufl
from dune.ufl import DirichletBC, Constant
from dune.fem.scheme import galerkin as solutionScheme
from dune.fem.utility import lineSample
from dune.ufl import expression2GF

class Problem_heat_D(Problem_heat):
    ## Dirichlet Problem, non-zero boundary to the right
    def __init__(self, gridsize, alpha, lambda_diff, order = 2, xa = -1, xb = 1):
        self.xa, self.xb, self.x_not_gamma = xa, xb, xa
        super(Problem_heat_D, self).__init__(gridsize, alpha, lambda_diff, order, xa = xa, xb = 0, mono = False)
        
        self.bc_side = DirichletBC(self.space, Constant(0.), self.x[0] < xa + self.eps)
        self.bcs.append(self.bc_side)
        
        self.u_gamma = self.space.interpolate(self.u0, name = 'u_gamma')
        self.bc_gamma = DirichletBC(self.space, self.u_gamma, self.x[0] > - self.eps)
        self.bcs.append(self.bc_gamma)
        
        self.scheme = solutionScheme([self.A == self.b, *self.bcs], solver = 'cg')
        
        self.normal = Constant((1., 0.))
        self.flux_grad = expression2GF(self.unew.space.grid, ufl.dot(ufl.grad(self.unew), self.normal), self.unew.space.order)
        self.flux_f = lambda x: float(self.lam)*lineSample(x, [0., 0.], [0., 1.], self.NN)[1]
        
        self.ug_interp = interp1d(self.yy, np.zeros(self.NN))
        self.ug_f = lambda x: self.ug_interp(x[1])

    def get_u0(self): return self.get_flux()
		
    def get_flux(self):
        return self.flux_f(self.flux_grad)
    
    def do_step(self, dt, ug):
        self.dt.value = dt
        
        self.ug_interp.y = ug
        self.u_gamma.interpolate(self.ug_f)
        
        self.scheme.solve(target = self.unew)
        self.uold.assign(self.unew)
        return self.get_flux()
    
    def solve(self, tf, n_steps):
        dt = tf/n_steps
        for i, t in enumerate(np.linspace(0, tf, n_steps + 1)[:-1]):
            self.get_ex_sol(t + dt) ## ex sol now rests in unew
            self.do_step(dt, self.get_u_gamma(self.unew))

class Problem_heat_N(Problem_heat):
    ## Neumann Problem, boundary to the left
    def __init__(self, gridsize, alpha, lambda_diff, order = 2, xa = -1, xb = 1):
        self.xa, self.xb, self.x_not_gamma = xa, xb, xb
        super(Problem_heat_N, self).__init__(gridsize, alpha, lambda_diff, order, xa = 0, xb = xb, mono = False)
        
        self.fold = self.space.interpolate(0., name = 'fold')
        self.fnew = self.space.interpolate(0., name = 'fnew')
        if order == 1:
            self.b -= self.dt*self.fnew*ufl.conditional(self.x[0] < self.eps, 1, 0)*self.v*ufl.ds
        elif order == 2:
            self.b -= 0.5*self.dt*(self.fnew + self.fold)*ufl.conditional(self.x[0] < self.eps, 1, 0)*self.v*ufl.ds
        else:
            raise ValueError('no implementation for this order')
        
        self.bc_side = DirichletBC(self.space, Constant(0.), self.x[0] > xb - self.eps)
        self.bcs.append(self.bc_side)
        
        self.scheme = solutionScheme([self.A == self.b, *self.bcs], solver = 'cg')
        
        
        self.f_old_interpol = interp1d(self.yy, np.zeros(self.NN))
        self.f_old_f = lambda x: self.f_old_interpol(x[1])
        
        self.f_new_interpol = interp1d(self.yy, np.zeros(self.NN))
        self.f_new_f = lambda x: self.f_new_interpol(x[1])
        
    def get_u0(self): return self.get_u_gamma(self.unew)
    
    def do_step(self, dt, flux_1, flux_2):
        self.dt.value = dt
        
        self.f_old_interpol.y = flux_1
        self.fold.interpolate(self.f_old_f)
        
        self.f_new_interpol.y = flux_2
        self.fnew.interpolate(self.f_new_f)

        self.scheme.solve(target = self.unew)
        self.uold.assign(self.unew)
        return self.get_u_gamma(self.unew)
    
    def solve(self, tf, n_steps):
        dt = tf/n_steps
        for i, t in enumerate(np.linspace(0, tf, n_steps + 1)[:-1]):
            self.get_ex_flux(t) # ex flux now sits in unew
            flux_old = np.copy(self.get_u_gamma(self.unew))
            
            self.get_ex_flux(t + dt) # ex flux now sits in unew
            flux_new = np.copy(self.get_u_gamma(self.unew))
            
            self.do_step(dt, flux_old, flux_new)

if __name__ == '__main__':
    pl.close('all')
    
    from heat_dune import get_solve_WR
    from verification import get_parameters
    solver = get_solve_WR(Problem_heat_D, Problem_heat_N)
    savefig = 'grad_'
    
    ## verify dirichlet and neumann solvers on their own
    from verification_D_N import verify_time
    pp = {'tf': 1., **get_parameters(), 'gridsize': 32, 'xa': -2, 'xb': 1}
#    verify_time(Problem_heat_D, True, k = 8, order = 1, **pp, savefig = savefig)
#    verify_time(Problem_heat_D, True, k = 8, order = 2, **pp, savefig = savefig)
    
#    verify_time(Problem_heat_N, False, k = 5, order = 1, **pp, savefig = savefig)
#    verify_time(Problem_heat_N, False, k = 5, order = 2, **pp, savefig = savefig)
    
    from verification_D_N import verify_space
    pp = {'tf': 1., **get_parameters(), 'xa': -2, 'xb': 1}
    
#    verify_space(Problem_heat_D, True, k = 6, order = 1, N_steps = 100, **pp, savefig = savefig)
#    verify_space(Problem_heat_D, True, k = 8, order = 2, N_steps = 50, **pp, savefig = savefig)
    
#    verify_space(Problem_heat_N, False, k = 6, order = 1, N_steps = 100, **pp, savefig = savefig)
#    verify_space(Problem_heat_N, False, k = 6, order = 2, N_steps = 50, **pp, savefig = savefig)

    ## verify WR converges against monolithic solution    
    from verification import verify_with_monolithic
    pp = {'tf': 1., **get_parameters(), 'gridsize': 32, 'xa': -2, 'xb': 1, 'theta': 0.5}

#    verify_with_monolithic(solve_WR = solver, k = 8, order = 1, **pp, savefig = savefig)
#    verify_with_monolithic(solve_WR = solver, k = 8, order = 2, **pp, savefig = savefig)
#    
    from verification import verify_comb_error
    ## verify combined error, splitting + time int for decreasing dt
#    verify_comb_error(solve_WR = solver, k = 8, order = 1, **pp, savefig = savefig)
#    verify_comb_error(solve_WR = solver, k = 8, order = 2, **pp, savefig = savefig)
    
    ## verify combined error, splitting + time int for decreasing dx
    from verification import verify_comb_error_space
    pp = {'tf': 1., **get_parameters(), 'xa': -2, 'xb': 1, 'theta': 0.5}
#    verify_comb_error_space(solve_WR = solver, k = 8, order = 1, **pp, savefig = savefig, TOL = 1e-3, N_steps = 100)
#    verify_comb_error_space(solve_WR = solver, k = 6, order = 2, **pp, savefig = savefig, TOL = 1e-6, N_steps = 50)
    
    ## verify convergence rate
    from verification import plot_theta
    pp = {'tf': 1., **get_parameters(), 'gridsize': 32, 'xa': -2, 'xb': 1}
#    plot_theta(solve_WR = solver, savefig = savefig, order = 1, **pp)
#    plot_theta(solve_WR = solver, savefig = savefig, order = 2, **pp)
    
    ## verify time-integration order with itself
    from verification import verify_self_time
    pp = {'tf': 1., **get_parameters(), 'gridsize': 32, 'xa': -2, 'xb': 1, 'theta': 0.5}
#    verify_self_time(solve_WR = solver, savefig = savefig, k = 6, order = 1, **pp)
#    verify_self_time(solve_WR = solver, savefig = savefig, k = 6, order = 2, **pp)