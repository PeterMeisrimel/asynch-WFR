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
from dune.fem.utility import Sampler
from dune.ufl import expression2GF
from dune.fem.operator import galerkin

class Problem_heat_D(Problem_heat):
    ## Dirichlet Problem, non-zero boundary to the right
    def __init__(self, gridsize, alpha, lambda_diff, order = 2, xa = -1, xb = 1, init = True):
        self.xa, self.xb, self.x_not_gamma = xa, xb, xa
        super(Problem_heat_D, self).__init__(gridsize, alpha, lambda_diff, order, xa = xa, xb = 0, mono = False)

        self.bc_side = DirichletBC(self.space, Constant(0.), self.x[0] < xa + self.eps)
        self.bcs.append(self.bc_side)

        self.u_gamma = self.space.interpolate(self.u0, name = 'u_gamma')
        self.bc_gamma = DirichletBC(self.space, self.u_gamma, self.x[0] > - self.eps)
        self.bcs.append(self.bc_gamma)

        self.scheme = solutionScheme([self.A == self.b, *self.bcs], solver = 'cg')
        
        self.normal = ufl.as_vector([1.0, 0.0])
        self.flux_grad = expression2GF(self.unew.space.grid, ufl.dot(ufl.grad(self.unew), self.normal), self.unew.space.order)
        self.sampler = Sampler(self.flux_grad)
        self.flux_f = lambda : float(self.lam)*self.sampler.lineSample([0., 0.], [0., 1.], self.NN)[1]

        self.ug_interp = interp1d(self.yy, np.zeros(self.NN))
        self.ug_f = lambda x: self.ug_interp(x[1])
        
        if init:
            ## first timestep may invoke some loading/compiling
            ## do this up front in case of timed-runs
            self.create_checkpoint()
            ## not sure what a suitable stepsize would actually be here
            self.do_step(1, np.zeros(self.NN))
            self.reset()
       
    def get_u0(self): 
        return self.flux_f()

    def get_flux(self):
        return self.flux_f()

    def do_step(self, dt, ug):
        self.scheme.model.dt = dt

        self.ug_interp.y = ug
        self.u_gamma.interpolate(self.ug_f)

        self.scheme.solve(target = self.unew)
        flux = self.get_flux()
        self.uold.assign(self.unew)
        return flux

    def solve(self, tf, n_steps):
        dt = tf/n_steps
        for i, t in enumerate(np.linspace(0, tf, n_steps + 1)[:-1]):
            self.get_ex_sol(t + dt) ## ex sol now rests in unew
            self.do_step(dt, self.get_u_gamma())
            
class Problem_heat_D_weak(Problem_heat_D):
    def __init__(self, gridsize, alpha, lambda_diff, order = 2, xa = -1, xb = 1):            
        ## required for initial reset, values do not matter
        self.flux0 = np.array([0])
        self.flux_old = np.array([0])
        super(Problem_heat_D_weak, self).__init__(gridsize, alpha, lambda_diff, order, xa, xb, init = False)
        self.setup_flux_expression()
        
        ## first timestep may invoke some loading/compiling
        ## do this up front in case of timed-runs
        self.create_checkpoint()
        ## not sure what a suitable stepsize would actually be here
        self.do_step(1, np.zeros(self.NN))
        self.reset()
        
    def setup_flux_expression(self):
        self.flux_lhs = self.u*self.v*ufl.dx
        if self.order == 1:
            self.flux_expr = (self.a*(self.u - self.uold)*self.v/self.dt*ufl.dx + 
                              self.lam*ufl.dot(ufl.grad(self.u), ufl.grad(self.v))*ufl.dx)
            self.flux_expr_operator = galerkin(self.flux_expr)
            self.get_flux = self.get_flux_weak1
        elif self.order == 2:
            self.flux_expr = (self.a*(self.u - self.uold)*self.v/self.dt*ufl.dx + 
                              0.5*self.lam*ufl.dot(ufl.grad(self.u + self.uold), ufl.grad(self.v))*ufl.dx)
            self.flux_expr_operator = galerkin(self.flux_expr)
            self.get_flux = self.get_flux_weak2
        else:
            raise KeyError('Order not available')
        self.flux_sol = self.space.interpolate(self.u0, name = 'flux_sol')
        
        self.flux_sol_expr = expression2GF(self.flux_sol.space.grid, self.flux_sol, self.flux_sol.space.order)
        self.sampler_weak_flux = Sampler(self.flux_sol_expr)
        self.flux_f_weak = lambda : self.sampler_weak_flux.lineSample([0., 0.], [0., 1.], self.NN)[1]
            
    def get_u0(self):
        flux = super(Problem_heat_D_weak, self).get_u0()
        self.flux0 = np.copy(flux)
        self.flux_old = np.copy(flux)
        return flux
    
    def reset(self):
        super(Problem_heat_D, self).reset()
        self.flux_old = np.copy(self.flux0)
    
    def create_checkpoint(self):
        super(Problem_heat_D, self).create_checkpoint()
        self.flux0 = np.copy(self.flux_old)
            
    def get_flux_weak1(self):
        ## compute flux, solution goes into self.flux_sol
        self.flux_expr_operator(self.unew, self.flux_sol)
        flux = self.flux_f_weak()/self.dx
        flux[0], flux[-1] = 0, 0
        self.flux_old = flux
        return flux
    
    def get_flux_weak2(self):
        self.flux_expr_operator(self.unew, self.flux_sol)
        flux = 2*self.flux_f_weak()/self.dx - self.flux_old
        flux[0], flux[-1]  = 0, 0
        self.flux_old = flux
        return flux
    
    def do_step(self, dt, ug):
        self.flux_expr_operator.model.dt = dt
        self.scheme.model.dt = dt

        self.ug_interp.y = ug
        self.u_gamma.interpolate(self.ug_f)

        self.scheme.solve(target = self.unew)
        flux = self.get_flux()
        self.uold.assign(self.unew)
        return flux
    
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
        
        self.do_initialization_step()
        
    def do_initialization_step(self):
        ## first timestep may invoke some loading/compiling
        ## do this up front in case of timed-runs
        self.create_checkpoint()
        ## not sure what a suitable stepsize would actually be here
        f = np.zeros(self.NN)
        self.do_step(1, f, f)
        self.reset()

    def get_u0(self): return self.get_u_gamma()

    def do_step(self, dt, flux_1, flux_2):
        self.scheme.model.dt = dt

        self.f_old_interpol.y = flux_1
        self.fold.interpolate(self.f_old_f)

        self.f_new_interpol.y = flux_2
        self.fnew.interpolate(self.f_new_f)

        self.scheme.solve(target = self.unew)
        self.uold.assign(self.unew)
        return self.get_u_gamma()

    def solve(self, tf, n_steps):
        dt = tf/n_steps
        for i, t in enumerate(np.linspace(0, tf, n_steps + 1)[:-1]):
            self.get_ex_flux(t) # ex flux now sits in unew
            flux_old = np.copy(self.get_u_gamma())

            self.get_ex_flux(t + dt) # ex flux now sits in unew
            flux_new = np.copy(self.get_u_gamma())

            self.do_step(dt, flux_old, flux_new)

if __name__ == '__main__':
    pl.close('all')

    from heat_dune import get_solve_WR
    from verification import get_parameters, verify_with_monolithic, verify_comb_error, verify_comb_error_space, plot_theta, verify_self_time
    from verification_D_N import verify_time, verify_space
    
    for i, (D_prob, savefig) in enumerate(zip([Problem_heat_D, Problem_heat_D_weak], ['grad_', 'weak_'])):
        if i == 0:
            continue
        solver = get_solve_WR(D_prob, Problem_heat_N)

        ## verify dirichlet and neumann solvers on their own
        pp = {'tf': 1., **get_parameters(), 'gridsize': 64, 'xa': -1, 'xb': 1}
#        verify_time(D_prob, True, k = 8, order = 1, **pp, savefig = savefig)
        verify_time(D_prob, True, k = 8, order = 2, **pp, savefig = savefig)
    
#        if i == 0:
#            verify_time(Problem_heat_N, False, k = 5, order = 1, **pp, savefig = savefig)
#            verify_time(Problem_heat_N, False, k = 5, order = 2, **pp, savefig = savefig)
    
        pp = {'tf': 1., **get_parameters(), 'xa': -1, 'xb': 1}
#        verify_space(D_prob, True, k = 6, order = 1, N_steps = 100, **pp, savefig = savefig)
        verify_space(D_prob, True, k = 8, order = 2, N_steps = 50, **pp, savefig = savefig)
    
#        if i == 0:
#            verify_space(Problem_heat_N, False, k = 6, order = 1, N_steps = 100, **pp, savefig = savefig)
#            verify_space(Problem_heat_N, False, k = 6, order = 2, N_steps = 50, **pp, savefig = savefig)
    
        ## verify WR converges against monolithic solution
        pp = {'tf': 1., **get_parameters(), 'gridsize': 64, 'xa': -1, 'xb': 1, 'theta': 0.5}
#        verify_with_monolithic(solve_WR = solver, k = 8, order = 1, **pp, savefig = savefig)
        verify_with_monolithic(solve_WR = solver, k = 8, order = 2, **pp, savefig = savefig)
    
        # verify combined error, splitting + time int for decreasing dt
#        verify_comb_error(solve_WR = solver, k = 8, order = 1, **pp, savefig = savefig)
        verify_comb_error(solve_WR = solver, k = 8, order = 2, **pp, savefig = savefig)
    
        ## verify combined error, splitting + time int for decreasing dx
        pp = {'tf': 1., **get_parameters(), 'xa': -1, 'xb': 1, 'theta': 0.5}
#        verify_comb_error_space(solve_WR = solver, k = 8, order = 1, **pp, savefig = savefig, TOL = 1e-3, N_steps = 100)
        verify_comb_error_space(solve_WR = solver, k = 6, order = 2, **pp, savefig = savefig, TOL = 1e-6, N_steps = 50)
    
        ## verify convergence rate
        pp = {'tf': 1., **get_parameters(), 'gridsize': 64, 'xa': -1, 'xb': 1}
#        plot_theta(solve_WR = solver, savefig = savefig, order = 1, **pp)
        plot_theta(solve_WR = solver, savefig = savefig, order = 2, **pp)
    
        ## verify time-integration order with itself
        pp = {'tf': 1., **get_parameters(), 'gridsize': 64, 'xa': -1, 'xb': 1, 'theta': 0.5}
#        verify_self_time(solve_WR = solver, savefig = savefig, k = 6, order = 1, **pp)
        verify_self_time(solve_WR = solver, savefig = savefig, k = 6, order = 2, **pp)
        
#    pp = {'tf': 1., **get_parameters(), 'gridsize': 512, 'xa': -1, 'xb': 1, 'theta': 0.5}
#    verify_comb_error(solve_WR = get_solve_WR(Problem_heat_D_weak, Problem_heat_N),
#                      k = 10, order = 1, **pp, savefig = '512_weak_')
#    verify_comb_error(solve_WR = get_solve_WR(Problem_heat_D, Problem_heat_N),
#                      k = 10, order = 1, **pp, savefig = '512_grad_')