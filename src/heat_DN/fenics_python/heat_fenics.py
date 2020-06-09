#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 22:31:31 2020

@author: Peter Meisrimel
"""
import numpy as np
import dolfin as dol
import pylab as pl
from scipy.interpolate import interp1d

## Class to make discrete data accessible to the weak forms, via linear interpolation
class Custom_Expr(dol.UserExpression):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
    def update_vals(self, grid, values):
        self.interpolated_data = interp1d(grid, values)
    def eval(self, value, x):
        value[0] = self.interpolated_data(x[1])
    def value_shape(self):
        return ()

## standard heat equation: alpha u_t + lamda_diff \Delta u = 0
## assume interface to be at zero, spatial domain marked by xa and xb
class Problem_heat:
    ## gridsize = number of internal variables
    def __init__(self, gridsize, alpha, lambda_diff, order = 2, xa = -1, xb = 1, mono = True):
        dol.set_log_level(30) # supress non critical ouput
        self.a, self.lam = dol.Constant(alpha), dol.Constant(lambda_diff)
        self.dt = dol.Constant(1.)
        ## n = gridsize = number of internal unknowns
        ## => n + 2 nodes per unit length
        ## => n + 1 cells per unit length
        self.yy = np.linspace(0, 1, gridsize + 2) ## meshpoints for interface values
        self.dx = 1./(gridsize + 1)
        
        ## logic for spatial domain: self.x marks the whole domain, while inputs to __init__ mark local domain
        if mono:
            self.xa, self.xb = xa, xb
            self.len = self.xb - self.xa ## length of x domain
            self.mesh = dol.RectangleMesh(dol.MPI.comm_self, dol.Point((self.xa, 0.)), dol.Point((self.xb, 1.)), self.len*(gridsize + 1), gridsize + 1)
            self.f_old, self.f_new = dol.Constant(0.), dol.Constant(0.)
        else:
            self.len = self.xb - self.xa ## length of x domain
            self.mesh = dol.RectangleMesh(dol.MPI.comm_self, dol.Point((xa, 0.)), dol.Point((xb, 1.)), (xb - xa)*(gridsize + 1), gridsize + 1)
        self.V = dol.FunctionSpace(self.mesh, "CG", 1)

        # initial conditions, fulfills Dirichlet boundaries
        self.u0_expr = "500*sin((x[0] - {})*M_PI/{})*sin(M_PI*x[1])".format(self.xa, self.len)
        self.u0 = dol.interpolate(dol.Expression(self.u0_expr, degree = 2), self.V)
        
        self.f0_expr = "{}*500*cos((x[0] - {})*M_PI/{})*sin(M_PI*x[1])*M_PI/{}".format(lambda_diff, self.xa, self.len, self.len)

        self.unew, self.vtest = dol.TrialFunction(self.V), dol.TestFunction(self.V)
        self.uold, self.usol = dol.Function(self.V), dol.Function(self.V)
        self.uold.interpolate(self.u0); self.usol.interpolate(self.u0)

        if order == 1: ## Implicit euler method
            self.F = (self.a*(self.unew - self.uold)*self.vtest*dol.dx
                      + self.dt*self.lam*dol.dot(dol.grad(self.unew), dol.grad(self.vtest))*dol.dx
                      - self.dt*self.f_new*self.vtest*dol.ds)
        elif order == 2: ## Crank-Nicolson method
            self.F = (self.a*(self.unew - self.uold)*self.vtest*dol.dx
                      + 0.5*self.dt*self.lam*dol.dot(dol.grad(self.uold + self.unew), dol.grad(self.vtest))*dol.dx
                      - 0.5*self.dt*(self.f_old + self.f_new)*self.vtest*dol.ds)
        self.lhs, self.rhs = dol.lhs(self.F), dol.rhs(self.F)

        if mono:
            ## marking zero-boundaries
            def boundary_zero(x, on_boundary):
                return on_boundary and (dol.near(x[0], self.xa) or dol.near(x[0], self.xb) or dol.near(x[1], 0) or dol.near(x[1], 1))
        else:            
            ## marking zero-boundaries
            def boundary_zero(x, on_boundary):
                return on_boundary and (dol.near(x[0], self.x_not_gamma) or dol.near(x[1], 0) or dol.near(x[1], 1))
        self.bc_D = dol.DirichletBC(self.V, dol.Constant(0.), boundary_zero)

        ## marking interface
        def boundary_gamma(x, on_boundary):
            return on_boundary and dol.near(x[0], 0.)
        self.bc_gamma = boundary_gamma
        
        self.create_checkpoint()
        self.reset()
        
    def reset(self):
        self.uold.vector()[:] = self.u_checkpoint.vector()
        self.usol.vector()[:] = self.u_checkpoint.vector()
        self.dt.assign(1.)
    
    def create_checkpoint(self):
        self.u_checkpoint = self.uold.copy()
        
    def get_u_gamma(self, ugamma):
        return np.array([ugamma(0., y) for y in self.yy])
        
    def solve(self, tf, n_steps):
        self.dt.assign(tf/n_steps)
        for _ in range(n_steps):
            dol.solve(self.lhs == self.rhs, self.usol, self.bc_D)
            self.uold.assign(self.usol)
    
    def get_sol(self, Nx, Ny, dx, xx = 0, offset = 0):
        res = np.zeros(Nx*Ny)
        for i in range(offset,  Nx):
            for j in range(Ny):
                res[(i - offset)*Ny + j] = self.usol(xx + i*dx, j*dx)
        return res
    
    def get_ex_sol(self, t = 1):
        t_fac = np.exp(-(self.len**2 + 1)/(self.len**2)*np.pi**2*t*float(self.lam)/float(self.a))
        self.usol.interpolate(dol.Expression(self.u0_expr + "*" + str(t_fac), degree = 2))
        
    def get_ex_flux(self, t = 1):
        t_fac = np.exp(-(self.len**2 + 1)/(self.len**2)*np.pi**2*t*float(self.lam)/float(self.a))
        self.usol.interpolate(dol.Expression(self.f0_expr + "*" + str(t_fac), degree = 2))
     
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
        print('rel tolerance factor', rel_tol_fac)
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
    pl.close('all')
    from verification import get_parameters, verify_space_error, verify_mono_time
    savefig = 'mono_'
    ## verify space order of monolithic solution
#    verify_space_error(1., k = 6, order = 2, **get_parameters(), xa = -2, xb = 1, savefig = savefig)
    
    ## verify time order with itself
    pp = {'tf': 1., **get_parameters(), 'gridsize': 64, 'xa': -2, 'xb': 1}
#    verify_mono_time(k = 7, order = 1, **pp, savefig = savefig)
#    verify_mono_time(k = 5, order = 2, **pp, savefig = savefig)