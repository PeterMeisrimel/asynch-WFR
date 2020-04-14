#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 22:31:31 2020

@author: Peter Meisrimel
"""
## visually verified for correctness

import numpy as np
import dolfin as dol
from scipy.interpolate import interp1d

"""
Careful with indentations here, using the standard ones from Spyder might break something for the C++ code
"""

class Custom_Expr(dol.UserExpression):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
    def update_vals(self, grid, values):
        self.interpolated_data = interp1d(grid, values)
    def eval(self, value, x):
        value[0] = self.interpolated_data(x[1])
    def value_shape(self):
        return ()
    
class Problem_heat:
    def __init__(self, gridsize, alpha, lam_diff):
        dol.set_log_level(30) # supress non critical ouput
        self.a = dol.Constant(alpha)
        self.lam = dol.Constant(lam_diff)
        self.dt = dol.Constant(1.)
        self.yy = np.linspace(0, 1, gridsize + 2) ## meshpoints for interface values
        
        self.mesh = dol.UnitSquareMesh(dol.MPI.comm_self, gridsize + 1, gridsize + 1)
        
        self.V = dol.FunctionSpace(self.mesh, "CG", 1)

        self.u0 = dol.interpolate(dol.Expression("500*sin((x[0] + {})*M_PI/2)*sin(M_PI*x[1])".format(self.x_init), degree = 2), self.V)

        self.unew = dol.TrialFunction(self.V)
        self.uold = dol.Function(self.V)
        self.uold.interpolate(self.u0)
        self.vtest = dol.TestFunction(self.V)
        self.usol = dol.Function(self.V)
        self.usol.interpolate(self.u0)

        self.F = (self.a*(self.unew - self.uold)*self.vtest*dol.dx
                  + 0.5*self.dt*self.lam*dol.dot(dol.grad(self.uold + self.unew), dol.grad(self.vtest))*dol.dx
                  + 0.5*self.dt*(self.f_old + self.f_new)*self.vtest*dol.ds)
        self.lhs = dol.lhs(self.F)
        self.rhs = dol.rhs(self.F)

        def boundary_zero(x, on_boundary):
            return on_boundary and (dol.near(x[0], self.x_not_gamma) or dol.near(x[1], 0) or dol.near(x[1], 1))
        self.bc_D = dol.DirichletBC(self.V, dol.Constant(0.), boundary_zero)

        def boundary_gamma(x, on_boundary):
            return on_boundary and dol.near(x[0], self.x_gamma)
        self.bc_gamma = boundary_gamma
        
        self.create_checkpoint()
        self.reset()
        
    def reset(self):
        self.uold.vector()[:] = self.u_checkpoint.vector()
        self.usol.vector()[:] = self.u_checkpoint.vector()
        self.dt.assign(1.)
    def create_checkpoint(self):
        self.u_checkpoint = self.uold.copy()

class Problem_heat_D(Problem_heat):
    ## Dirichlet Problem, boundary to the right
    def __init__(self, gridsize, alpha, lam_diff):
        self.x_init = 0
        self.x_gamma, self.x_not_gamma = 1., 0.
        self.f_old = dol.Constant(0.)
        self.f_new = dol.Constant(0.)
        super(Problem_heat_D, self).__init__(gridsize, alpha, lam_diff)
        
        self.F_flux = (self.a*(self.usol - self.uold)/self.dt*self.vtest*dol.dx
                       + self.lam*dol.dot(dol.grad(self.usol), dol.grad(self.vtest))*dol.dx)
        self.ugamma = Custom_Expr()

    def get_u0(self):
        return self.get_flux(1)
		
    def get_flux(self, dt):
#        dx = self.yy[1] - self.yy[0]
#        return float(self.lam)*np.array([0.] + [(self.usol(self.x_gamma, y) - self.usol(self.x_gamma  - dx, y))/dx for y in self.yy[1:-1]] + [0.])
        self.dt.assign(dt)
        flux_sol = dol.assemble(self.F_flux)
        flux_f = dol.Function(self.V)
        flux_f.vector().set_local(flux_sol)
        return -np.array([0.] + [flux_f(self.x_gamma, y) for y in self.yy[1:-1]] + [0.]) # manually enfore zero fluxes at upper and lower bound

    def do_step(self, dt, ug):
        self.dt.assign(dt)
        self.ugamma.update_vals(self.yy, ug)

        bc_gamma = dol.DirichletBC(self.V, self.ugamma, self.bc_gamma)
        dol.solve(self.lhs == self.rhs, self.usol, [self.bc_D, bc_gamma])
        flux = self.get_flux(dt)
        self.uold.assign(self.usol)
        return flux      
    
class Problem_heat_N(Problem_heat):
    ## Neumann Problem, boundary to the left
    def __init__(self, gridsize, alpha, lam_diff):
        self.x_init = 1
        self.x_gamma, self.x_not_gamma = 0, 1

        self.f_old = Custom_Expr()
        self.f_new = Custom_Expr()
        super(Problem_heat_N, self).__init__(gridsize, alpha, lam_diff)

    def get_u0(self):
        return np.array([0.] + [self.uold(self.x_gamma, y) for y in self.yy[1:-1]] + [0.])
		
    def do_step(self, dt, flux_old, flux_new):
        self.dt.assign(dt)
        self.f_old.update_vals(self.yy, flux_old)
        self.f_new.update_vals(self.yy, flux_new)

        dol.solve(self.lhs == self.rhs, self.usol, self.bc_D)
        self.uold.assign(self.usol)
        return np.array([0.] + [self.usol(self.x_gamma, y) for y in self.yy[1:-1]] + [0.])

if __name__ == '__main__':
    print('running python mode, normally')
    p = {'alpha': 1., 'lam_diff': 0.01, 'gridsize': 32}
    tf = 1
    N = 100
    dt = tf/N
    maxiter = 40

    pD = Problem_heat_D(**p)
    pN = Problem_heat_N(**p)
    
    ug0 = pN.get_u0()
    print('ug0', pN.get_u0())
    print('flux0', pD.get_u0())

    ## run Waveform iteration 
    ug = [np.copy(ug0) for i in range(N + 1)]
    pD.create_checkpoint()
    pN.create_checkpoint()
    for k in range(maxiter):
        pD.reset()
        pN.reset()
        f = [np.copy(pD.get_flux(1.))]
        for i in range(N):
            f.append(pD.do_step(dt, ug[i+1]))
        tmp = np.copy(ug[-1])
        for i in range(N):
            ug[i+1] = pN.do_step(dt, f[i], f[i+1])
        update = np.linalg.norm(ug[-1] - tmp, 2)
        print(k, ' update = ', update)
        if update < 1e-8:
            break
    print(f[-1])
    print(ug[-1])

    ## plotting
    plotting_points = 64
    plotting_dx = 1/(plotting_points + 1)

    plot_mat_D = np.zeros((plotting_points+1, plotting_points+1))
    plot_mat_N = np.zeros((plotting_points+1, plotting_points+1))
    for i in range(plotting_points + 1):
        for j in range(plotting_points + 1):
            plot_mat_D[i, j] = pD.uold(i*plotting_dx, j*plotting_dx)
            plot_mat_N[i, j] = pN.uold(i*plotting_dx, j*plotting_dx)
   
    plot_total = np.hstack((plot_mat_D.T, plot_mat_N[1:].T))
    import pylab as pl
    pl.close('all')
    fig, ax = pl.subplots(figsize = (20, 10))
    p = ax.pcolor(plot_total)
    pl.colorbar(p)
    fig.savefig('test_total.png', dpi = 100)