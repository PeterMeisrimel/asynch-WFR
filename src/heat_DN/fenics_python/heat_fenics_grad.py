#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 22:31:31 2020

@author: Peter Meisrimel
"""

import dolfin as dol
import pylab as pl
import numpy as np
from heat_fenics import Problem_heat, Custom_Expr

class Problem_heat_D(Problem_heat):
    ## Dirichlet Problem, boundary to the right, standard computation of flux via gradient of solution
    def __init__(self, gridsize, alpha, lambda_diff, order = 2, xa = -1, xb = 1):
        self.xa, self.xb, self.x_not_gamma = xa, xb, xa
        self.f_old, self.f_new = dol.Constant(0.), dol.Constant(0.)
        
        super(Problem_heat_D, self).__init__(gridsize, alpha, lambda_diff, order, xa = xa, xb = 0, mono = False)
        
        self.n = dol.Constant((1.0, 0.0))
        self.F_flux = (self.lam*dol.dot(dol.grad(self.usol), self.n)*self.vtest*dol.ds)
        
        self.flux_f = dol.Function(self.V)
        self.ugamma = Custom_Expr() ## for discrete Dirichlet boundary data
        
        self.bc_gamma = dol.DirichletBC(self.V, self.ugamma, self.bc_gamma)
        self.bcs = [self.bc_D, self.bc_gamma]
        
    def get_flux(self):
        flux_sol = dol.assemble(self.F_flux)
        self.flux_f.vector().set_local(flux_sol)
#        flux = self.get_u_gamma(self.flux_f)
        flux = self.get_u_gamma(self.flux_f)/self.dx
        flux[0], flux[-1] = 0, 0 # manually enforce zero flux at endpoints
        print('fluxxxx1', flux)
        return flux
    
    def get_u0(self):
        return self.get_flux()
    
    def do_step(self, dt, ug):
        print('do D step ', dt)
        self.dt.assign(dt)
        self.ugamma.update_vals(self.yy, ug) ## set boundary data
        dol.solve(self.lhs == self.rhs, self.usol, self.bcs)
        flux = self.get_flux() # potentially based on uold & usol
        self.uold.assign(self.usol)
        return flux

    def solve(self, tf, n_steps):
        self.get_u0() ## for Dirichlet only verification
        dt = tf/n_steps
        for i, t in enumerate(np.linspace(0, tf, n_steps + 1)[:-1]):
            self.get_ex_sol(t + dt) # get solution for u_gamma into self.usol
            self.do_step(dt, self.get_u_gamma(self.usol))
            
class Problem_heat_D_weak(Problem_heat_D):
    ## Dirichlet Problem, boundary to the right, compute flux via weak form
    def __init__(self, gridsize, alpha, lambda_diff, order = 2, xa = -1, xb = 1):            
        ## required for initial reset, values do not matter
        self.flux0, self.flux_old = np.array([0]), np.array([0])
        super(Problem_heat_D_weak, self).__init__(gridsize, alpha, lambda_diff, order, xa, xb)
        
        if order == 1:
            self.F_flux = (self.a*(self.usol - self.uold)*self.vtest/self.dt*dol.dx + 
                           self.lam*dol.dot(dol.grad(self.usol), dol.grad(self.vtest))*dol.dx)
        elif order == 2:
            self.F_flux = (self.a*(self.usol - self.uold)*self.vtest/self.dt*dol.dx + 
                           0.5*self.lam*dol.dot(dol.grad(self.usol + self.uold), dol.grad(self.vtest))*dol.dx)
            self.get_flux = self.get_flux2
            
    def get_u0(self):
        ## essentially the one from Problem_heat_D
        F_flux0 = dol.assemble(self.lam*dol.dot(dol.grad(self.usol), dol.Constant((1.0, 0.0)))*self.vtest*dol.ds)
        self.flux_f.vector().set_local(F_flux0)
        self.flux0 = self.get_u_gamma(self.flux_f)/self.dx
        self.flux_old = np.copy(self.flux0)
        return self.flux_old
            
    def get_flux2(self):
        flux_sol = dol.assemble(self.F_flux)
        self.flux_f.vector().set_local(flux_sol)
        flux_new = 2*self.get_u_gamma(self.flux_f)/self.dx - self.flux_old
        flux_new[0], flux_new[-1] = 0, 0
        self.flux_old = flux_new
        print('fluxxx2', flux_new)
        return flux_new
    
    def reset(self):
        super(Problem_heat_D, self).reset()
        self.flux_old = np.copy(self.flux0)
        
    def create_checkpoint(self):
        super(Problem_heat_D, self).create_checkpoint()
        self.flux0 = np.copy(self.flux_old)
    
class Problem_heat_N(Problem_heat):
    ## Neumann Problem, boundary to the left
    def __init__(self, gridsize, alpha, lambda_diff, order = 2, xa = -1, xb = 1):
        self.xa, self.xb, self.x_not_gamma = xa, xb, xb

        self.f_old, self.f_new = Custom_Expr(), Custom_Expr()
        super(Problem_heat_N, self).__init__(gridsize, alpha, lambda_diff, order, xa = 0, xb = xb, mono = False)

    def get_u0(self): return self.get_u_gamma(self.usol)
    
    def do_step(self, dt, flux_1, flux_2):
        self.dt.assign(dt)
        ## update fluxes
        self.f_old.update_vals(self.yy, flux_1)
        self.f_new.update_vals(self.yy, flux_2)

        dol.solve(self.lhs == self.rhs, self.usol, self.bc_D)
        self.uold.assign(self.usol)
        return self.get_u_gamma(self.usol)
    
    def solve(self, tf, n_steps):
        dt = tf/n_steps
        for i, t in enumerate(np.linspace(0, tf, n_steps + 1)[:-1]):
            self.get_ex_flux(t) # get solution for flux into self.usol
            flux_old = np.copy(self.get_u_gamma(self.usol))
            
            self.get_ex_flux(t + dt) # get solution for flux into self.usol
            flux_new = np.copy(self.get_u_gamma(self.usol))
            
            self.do_step(dt, flux_old, flux_new)
    
if __name__ == '__main__':
    pl.close('all')
    
    from heat_fenics import get_solve_WR
    from verification import get_parameters, verify_with_monolithic, verify_comb_error, verify_comb_error_space, plot_theta, verify_self_time
    from verification_D_N import verify_time, verify_space
    
    #solv = get_solve_WR(Problem_heat_D_weak, Problem_heat_N)
    #print(solv(1, 20, 20, 100, 1e-10, 0.5, xa = -1, xb = 1, gridsize = 32, alpha = 1, lambda_diff = 0.01)[2:])
    
    for i, (D_prob, savefig) in enumerate(zip([Problem_heat_D, Problem_heat_D_weak], ['grad_', 'weak_'])):
        solver = get_solve_WR(D_prob, Problem_heat_N)
        if i == 0: 
            continue
    
        ## verify dirichlet and neumann solvers on their own
        ## for refinement in time
        pp = {'tf': 1., **get_parameters(), 'gridsize': 64, 'xa': -1, 'xb': 1}
#        verify_time(D_prob, True, k = 8, order = 1, **pp, savefig = savefig)
#        verify_time(D_prob, True, k = 8, order = 2, **pp, savefig = savefig)
        
#        if i == 0:
#            verify_time(Problem_heat_N, False, k = 6, order = 1, **pp, savefig = savefig)
#            verify_time(Problem_heat_N, False, k = 6, order = 2, **pp, savefig = savefig)
        
        ## for refinement in space
        pp = {'tf': 1., **get_parameters(), 'xa': -1, 'xb': 1}
        
#        verify_space(D_prob, True, k = 8, order = 1, N_steps = 100, **pp, savefig = savefig)
#        verify_space(D_prob, True, k = 8, order = 2, N_steps = 100, **pp, savefig = savefig)
        
#        if i == 0:
#            verify_space(Problem_heat_N, False, k = 8, order = 1, N_steps = 200, **pp, savefig = savefig)
#            verify_space(Problem_heat_N, False, k = 8, order = 2, N_steps = 50, **pp, savefig = savefig)
    
        pp = {'tf': 1., **get_parameters(), 'gridsize': 64, 'xa': -1, 'xb': 1, 'theta': 0.5}
        ## verify WR converges against monolithic solution
        verify_with_monolithic(solve_WR = solver, k = 8, order = 1, **pp, savefig = savefig)
#        verify_with_monolithic(solve_WR = solver, k = 8, order = 2, **pp, savefig = savefig)
        
        ## verify combined error, splitting + time int for decreasing dt
#        verify_comb_error(solve_WR = solver, k = 8, order = 1, **pp, savefig = savefig)
#        verify_comb_error(solve_WR = solver, k = 8, order = 2, **pp, savefig = savefig)
        
        ## verify combined error, splitting + time int for decreasing dx
        pp = {'tf': 1., **get_parameters(), 'xa': -1, 'xb': 1, 'theta': 0.5}
#        verify_comb_error_space(solve_WR = solver, k = 8, order = 1, **pp, savefig = savefig, TOL = 1e-10, N_steps = 100)
#        verify_comb_error_space(solve_WR = solver, k = 8, order = 2, **pp, savefig = savefig, TOL = 1e-10, N_steps = 100)
        
        ## verify convergence rate
        pp = {'tf': 1., **get_parameters(), 'gridsize': 64, 'xa': -1, 'xb': 1}
#        plot_theta(solve_WR = solver, savefig = savefig, order = 1, **pp)
#        plot_theta(solve_WR = solver, savefig = savefig, order = 2, **pp)
        
        ## verify time-integration order with itself
        pp = {'tf': 1., **get_parameters(), 'gridsize': 64, 'xa': -1, 'xb': 1, 'theta': 0.5}
#        verify_self_time(solve_WR = solver, savefig = savefig, k = 6, order = 1, **pp)
#        verify_self_time(solve_WR = solver, savefig = savefig, k = 6, order = 2, **pp)

    pp = {'tf': 1., **get_parameters(), 'gridsize': 512, 'xa': -1, 'xb': 1, 'theta': 0.5}
#    verify_comb_error(solve_WR = get_solve_WR(Problem_heat_D_weak, Problem_heat_N),
#                      k = 10, order = 1, **pp, savefig = '512_weak_')
#    verify_comb_error(solve_WR = get_solve_WR(Problem_heat_D, Problem_heat_N),
#                      k = 10, order = 1, **pp, savefig = '512_grad_')