"""
Created on Thu Jun 11 14:52:41 2020

modified from dune-fempy demo: 
https://www.dune-project.org/sphinx/content/sphinx/dune-fem/euler_nb.html (May 2020)
@author: Peter Meisrimel, Lund University
"""

## Euler System of Gas Dynamics

import numpy as np
import matplotlib
matplotlib.rc( 'image', cmap='jet' )
from matplotlib import pyplot
import ufl
from dune.grid import reader
from dune.fem.space import finiteVolume
from dune.fem.function import gridFunction
from dune.femdg import femDGOperator
from dune.femdg.rk import femdgStepper
from dune.alugrid import aluSimplexGrid as simplexGrid
from dune.ufl import expression2GF
from dune.fem.utility import lineSample
from scipy.interpolate import interp1d

# Basic model for hyperbolic conservation law
## has not been touched, just merged the velocity, physical and jump functions into the class
#class Model:
#    ## U = 4 variables
#    ## primitive version = 3 variables with velocity being a vector!
#    gamma = 1.4
#    R = 287.058
#    # helper function
#    def toPrim(U):
#        v = ufl.as_vector([U[i]/U[0] for i in range(1,3)])
#        kin = ufl.dot(v,v) * U[0]/2
#        pressure = (Model.gamma - 1)*(U[3] - kin)
#        return U[0], v, pressure
#
#    # interface methods for model
#    def F_c(t, x, U):
#        rho, v, p = Model.toPrim(U)
#        return ufl.as_matrix([
#                  [rho*v[0], rho*v[1]],
#                  [rho*v[0]*v[0] + p, rho*v[0]*v[1]],
#                  [rho*v[0]*v[1], rho*v[1]*v[1] + p],
#                  [(U[3] + p)*v[0], (U[3] + p)*v[1]]])
#    # simple 'outflow' boundary conditions on all boundaries
#    boundary = {1: lambda t, x, U: [U[0], 0, 0, U[3]],
#                2: lambda t, x, U: [U[0], 0, 0, U[3]],
#                3: lambda t, x, U: [U[0], 0, 0, U[3]],
#                4: lambda t, x, U: [U[0], 0, 0, U[3]]}
#    ## which is which? possibly fix these via ufl.conditionals? 
##    boundary = {range(1, 5): lambda t, x, U: [10, 0, 0, 10]}
#
#    # interface method needed for LLF and time step control
#    def maxLambda(t, x, U, n):
#        rho, v, p = Model.toPrim(U)
#        return abs(ufl.dot(v,n)) + ufl.sqrt(Model.gamma*p/rho)
#
#    ## methods for limiting
#    def velocity(t, x, U):
#        _, v, _ = Model.toPrim(U)
#        return v
#    def physical(t, x, U):
#        rho, _, p = Model.toPrim(U)
#        return ufl.conditional(rho > 1e-8, ufl.conditional(p > 1e-8, 1, 0), 0)
#    def jump(t, x, U, V):
#        _, _, pL = Model.toPrim(U)
#        _, _, pR = Model.toPrim(V)
#        return (pL - pR)/(0.5*(pL + pR))
#    
#    def temp(U):
#        rho, _, p = Model.toPrim(U)
#        return p/(rho*Model.R)

class Model:
    ## U = 4 variables
    ## primitive version = 3 variables with velocity being a vector!
    gamma = 1.4
    R = 287.058
    # helper function
    def toPrim(self, U):
        v = ufl.as_vector([U[i]/U[0] for i in range(1,3)])
        kin = ufl.dot(v,v) * U[0]/2
        pressure = (self.gamma - 1)*(U[3] - kin)
        return U[0], v, pressure

    # interface methods for model
    def F_c(self, t, x, U):
        rho, v, p = self.toPrim(U)
        return ufl.as_matrix([
                  [rho*v[0], rho*v[1]],
                  [rho*v[0]*v[0] + p, rho*v[0]*v[1]],
                  [rho*v[0]*v[1], rho*v[1]*v[1] + p],
                  [(U[3] + p)*v[0], (U[3] + p)*v[1]]])
    # simple 'outflow' boundary conditions on all boundaries
    boundary = {1: lambda t, x, U: [U[0], -U[1], -U[2], U[3]],
                2: lambda t, x, U: [U[0], -U[1], -U[2], U[3]],
                3: lambda t, x, U: [U[0], U[1], U[2], U[3]],
                4: lambda t, x, U: [U[0], U[1], U[2], U[3]]}
    ## which is which? possibly fix these via ufl.conditionals? 
#    boundary = {range(1, 5): lambda t, x, U: [10, 0, 0, 10]}

    # interface method needed for LLF and time step control
    def maxLambda(self, t, x, U, n):
        rho, v, p = self.toPrim(U)
        return abs(ufl.dot(v,n)) + ufl.sqrt(self.gamma*p/rho)

    ## methods for limiting
    def velocity(self, t, x, U):
        _, v, _ = self.toPrim(U)
        return v
    def physical(self, t, x, U):
        rho, _, p = self.toPrim(U)
        return ufl.conditional(rho > 1e-8, ufl.conditional(p > 1e-8, 1, 0), 0)
    def jump(self, t, x, U, V):
        _, _, pL = self.toPrim(U)
        _, _, pR = self.toPrim(V)
        return (pL - pR)/(0.5*(pL + pR))
    
    def temp(self, U):
        rho, _, p = self.toPrim(U)
        return p/(rho*self.R)
    
## TODO
"""
how to disable alugrid messages
Implicit time integration?
Change domain/grid
Change boundaries
input of boundary data
Fluid solver takes heat flux, i.e., energy flux as input and outputs boundary temperature
"""
    
class EulerSolver:
    ## get some more inputs, e.g. various orders
    def __init__(self, order_space = 2, order_time = 2):
        self.x = ufl.SpatialCoordinate(ufl.triangle)
        self.u0 = ufl.conditional(ufl.sqrt(ufl.dot(self.x, self.x)) < 0.5, ## circle with radius 0.05
                                  ufl.as_vector([1, 0, 0, 2.5]), ## inside, high density + pressure
                                  ufl.as_vector([0.125, 0, 0, 0.25])) ## outside, low density  + pressure
        
        self.domain = (reader.dgf, "triangle.dgf") ## read domain specifications from file
        self.gridView = simplexGrid(self.domain, dimgrid = 2) ## grid with triangles, using alugrid module
        self.space = finiteVolume(self.gridView, dimRange = 4) ## solution space
        self.uh = self.space.interpolate(self.u0, name = "solution")
        
        self.local_uh = self.uh.localFunction()
        
        ## this one causes the crash upon deconstruction of sorts?
        ## current error, something is trying to be freed, that was not allocated by malloc.
        self.rho_gf = gridFunction(self.space.grid, name = "rho", order = order_space)(self.rho)
        
        self.Model = Model()
#        self.operator = femDGOperator(Model, self.space, limiter = "MinMod") # note: that u_h.space fails since not ufl_space
        self.operator = femDGOperator(self.Model, self.space, limiter = "MinMod") # note: that u_h.space fails since not ufl_space
        ## does not seem to work, resp. take way too long for implicit?
        self.stepper = femdgStepper(order = order_time, operator = self.operator, rkType = 'EX')
        self.operator.applyLimiter(self.uh)
        
        ## plotting
        self.c, self.saveStep = 0, 0.01
        self.fig = pyplot.figure(figsize = (30, 10))
        self.rho_gf.plot(gridLines = "white", figure = (self.fig, 131 + self.c), colorbar = False, clim = [0.125, 1])
        
        ## flux computation
        ## how to compute flux only for temperature? resp. only for rho + pressure?
#        self.flux_grad = expression2GF(self.uh.space.grid, ufl.grad(self.uh[3]), self.uh.space.order) ## this would be the energy gradient?
#        self.flux_grad = expression2GF(self.uh.space.grid, ufl.grad(Model.temp(self.uh)), self.uh.space.order)
        self.flux_grad = expression2GF(self.uh.space.grid, ufl.grad(self.Model.temp(self.uh)), self.uh.space.order)
        self.flux_f = lambda x: lineSample(x, [-1., -1.], [-1., 1.], self.NN)[1]
        
    def rho(self, e, x):
        self.local_uh.bind(e)
        return self.local_uh(x)[0]
    
    ## this one might just work?
    def get_flux(self):
#        f = self.flux_f(self.uh) ## currently assume this would be fluxes of [rho, rho v1, rho v2, rho E]
#        return Model.temp(f)
        return self.flux_f(self.flux_grad)
        
    def do_step(self, t, dt):
        self.operator.setTime(t)
        ## this doesn't appear to actually do anything?
        self.Model.boundary[1] = lambda t, x, U: [-(t + 100)*U[0], -200*U[1], -200*U[2], -(t + 100)*U[3]]
        self.t += self.stepper(self.uh, dt = dt)
        ## do plotting every self.saveStep time units
        if self.t > self.saveStep:
            self.saveStep += 0.1
            self.c += 1
            self.rho_gf.plot(gridLines = "white", figure = (self.fig, 131 + self.c), colorbar = False, clim = [0.125, 1])
            
    def solve(self, tf, n):
        dt = tf/n
        self.t = 0
        for _ in range(n):
            print('doing step', self.t, dt, flush = True)
            self.do_step(self.t, dt)
        pyplot.show()
        
#    def __del__(self):
#        print('deconstructor')
#        pass
            
if __name__ == '__main__':
#    tf = 0.3
#    n = 100
    tf = 0.6
    n = 200
    solver = EulerSolver()
    solver.solve(tf, n)
    print('solving over')

# time integration
#def evolve(space, u_h):
#    lu_h = u_h.localFunction()
#    @gridFunction(space.grid, name = "rho", order = space.order)
#    def rho(e, x):
#        lu_h.bind(e)
#        return lu_h(x)[0]
#    operator = femDGOperator(Model, space, limiter = "MinMod") # note: that u_h.space fails since not ufl_space
#    stepper  = femdgStepper(order = space.order if space.order else 2, operator = operator)
#    operator.applyLimiter(u_h)
#    t = 0
#    saveStep = 0.01
#    c = 0
#    fig = pyplot.figure(figsize = (30, 10))
#    rho.plot(gridLines = "white", figure = (fig, 131 + c), colorbar = False, clim = [0.125, 1])
#    while t < 0.3:
#        operator.setTime(t)
#        t += stepper(u_h)
#        print(t)
#        if t > saveStep:
#            saveStep += 0.1
#            c += 1
#            rho.plot(gridLines = "white", figure=(fig, 131 + c), colorbar = False, clim = [0.125, 1])
#    # return lineSample(u_h,[0,0],[1,1],200)
#    res = numpy.zeros((2, space.grid.size(0)))
#    for i,e in enumerate(space.grid.elements):
#        x = e.geometry.center
#        res[0][i] = x.two_norm
#        res[1][i] = rho(e, e.geometry.toLocal(x))
#    pyplot.show()
#    return res


## initial condition
#x = ufl.SpatialCoordinate(ufl.triangle)
##initial = ufl.conditional(ufl.dot(x, x) < 0.1, ## circle with radius 0.1
##                          ufl.as_vector([1, 0, 0, 2.5]), ## inside, high density + pressure
##                          ufl.as_vector([0.125, 0, 0, 0.25])) ## outside, low density  + pressure
#initial = ufl.conditional(ufl.sqrt(ufl.dot(x, x)) < 0.5, ## circle with radius 0.05
#                          ufl.as_vector([1, 0, 0, 2.5]), ## inside, high density + pressure
#                          ufl.as_vector([0.125, 0, 0, 0.25])) ## outside, low density  + pressure
#
#domain = (reader.dgf, "triangle.dgf") ## read domain specifications from file
#gridView = simplexGrid(domain, dimgrid = 2) ## grid with triangles, using alugrid module
#
#space = finiteVolume(gridView, dimRange = 4) ## solution space
#
#u_h   = space.interpolate(initial, name = "solution")
#res   = [evolve(space, u_h)]

## not sure what I am actually seeing here?
## Solutions along the diagonal
#pyplot.scatter(res[0][0], res[0][1])
#pyplot.scatter(res[1][0], res[1][1])
# pyplot.scatter(res[2][0],res[2][1])
#pyplot.show()
