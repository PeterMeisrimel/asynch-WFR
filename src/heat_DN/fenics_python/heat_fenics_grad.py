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
    ## Dirichlet Problem, boundary to the right
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

    def get_u0(self): 
        return self.get_flux()
		
    def get_flux(self):
        flux_sol = dol.assemble(self.F_flux)
        self.flux_f.vector().set_local(flux_sol/self.dx)
        return self.get_u_gamma(self.flux_f)
    
    def do_step(self, dt, ug):
        self.dt.assign(dt)
        self.ugamma.update_vals(self.yy, ug) ## set boundary data
        dol.solve(self.lhs == self.rhs, self.usol, self.bcs)
        self.uold.assign(self.usol)
        return self.get_flux() # based on usol

    def solve(self, tf, n_steps):
        dt = tf/n_steps
        for i, t in enumerate(np.linspace(0, tf, n_steps + 1)[:-1]):
            self.get_ex_sol(t + dt) # get solution for u gamma into self.usol
            self.do_step(dt, self.get_u_gamma(self.usol))
    
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
    from verification import get_parameters, verify_with_monolithic, verify_comb_error, verify_comb_error_space, plot_theta
    solver = get_solve_WR(Problem_heat_D, Problem_heat_N)
    savefig = 'grad_'
    
    ## verify dirichlet and neumann solvers on their own
    from verification_D_N import verify_time
    pp = {'tf': 1., **get_parameters(), 'gridsize': 128, 'xa': -2, 'xb': 1}
#    verify_time(Problem_heat_D, True, k = 8, order = 1, **pp, savefig = savefig)
#    verify_time(Problem_heat_D, True, k = 8, order = 2, **pp, savefig = savefig)
    
#    verify_time(Problem_heat_N, False, k = 6, order = 1, **pp, savefig = savefig)
#    verify_time(Problem_heat_N, False, k = 6, order = 2, **pp, savefig = savefig)
    
    from verification_D_N import verify_space
    pp = {'tf': 1., **get_parameters(), 'xa': -2, 'xb': 1}
    
#    verify_space(Problem_heat_D, True, k = 6, order = 1, N_steps = 100, **pp, savefig = savefig)
#    verify_space(Problem_heat_D, True, k = 8, order = 2, N_steps = 50, **pp, savefig = savefig)
    
#    verify_space(Problem_heat_N, False, k = 8, order = 1, N_steps = 200, **pp, savefig = savefig)
#    verify_space(Problem_heat_N, False, k = 8, order = 2, N_steps = 50, **pp, savefig = savefig)
    
    pp = {'tf': 1., **get_parameters(), 'gridsize': 32, 'xa': -2, 'xb': 1, 'theta': 0.5}
    ## verify WR converges against monolithic solution
#    verify_with_monolithic(solve_WR = solver, k = 8, order = 1, **pp, savefig = savefig)
#    verify_with_monolithic(solve_WR = solver, k = 8, order = 2, **pp, savefig = savefig)
    
    ## verify combined error, splitting + time int for decreasing dt
#    verify_comb_error(solve_WR = solver, k = 8, order = 1, **pp, savefig = savefig)
#    verify_comb_error(solve_WR = solver, k = 8, order = 2, **pp, savefig = savefig)
    
    ## verify combined error, splitting + time int for decreasing dx
    pp = {'tf': 1., **get_parameters(), 'xa': -2, 'xb': 1, 'theta': 0.5}
#    verify_comb_error_space(solve_WR = solver, k = 8, order = 1, **pp, savefig = savefig, TOL = 1e-3, N_steps = 100)
#    verify_comb_error_space(solve_WR = solver, k = 6, order = 2, **pp, savefig = savefig, TOL = 1e-6, N_steps = 50)
    
    ## verify convergence rate
    pp = {'tf': 1., **get_parameters(), 'gridsize': 32, 'xa': -2, 'xb': 1}
#    plot_theta(solve_WR = solver, savefig = savefig, order = 1, **pp)
#    plot_theta(solve_WR = solver, savefig = savefig, order = 2, **pp)
    
    ## verify time-integration order with itself
    from verification import verify_self_time
    pp = {'tf': 1., **get_parameters(), 'gridsize': 16, 'xa': -2, 'xb': 1, 'theta': 0.5}
#    verify_self_time(solve_WR = solver, savefig = savefig, k = 6, order = 1, **pp)
#    verify_self_time(solve_WR = solver, savefig = savefig, k = 6, order = 2, **pp)