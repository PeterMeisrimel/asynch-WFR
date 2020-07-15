"""
Created on Thu Jun 11 14:52:41 2020

modified from dune-fempy demo: 
https://www.dune-project.org/sphinx/content/sphinx/dune-fem/euler_nb.html (May 2020)
@author: Peter Meisrimel, Lund University, with lots of help from Robert Klöfkorn
"""

## Euler System of Gas Dynamics

import numpy as np
import matplotlib
matplotlib.rc('image', cmap = 'jet')
import ufl
from dune.grid import reader
from dune.fem.space import finiteVolume
from dune.fem.function import gridFunction
from dune.femdg import femDGOperator
from dune.femdg.rk import femdgStepper
from dune.alugrid import aluSimplexGrid as simplexGrid
from dune.ufl import expression2GF, Constant
from dune.fem.utility import lineSample
from scipy.interpolate import interp1d

xx = np.linspace(-1, 1, 100)
aa = np.random.rand(100)*20
interp = interp1d(xx, aa)
def bound(x):
#    print('aaaa', x, type(x), type(x[0]))
#    return float(interp([x[0]])[0])
    return float(interp([x])[0])
#bound = lambda x: float(interp(x[0]))
fff = lambda t: (t + 1)**20

def fxx(t, x):
    print(t, x, type(t), type(x))
    return (t+1)

## boundary takes whatever the input function is and converts it to a ufl function, no evaluations of this particular one except for initialization

class Model:
    ## U = 4 variables
    ## primitive version = 3 variables with velocity being a vector!
    
    ## source https://www.ohio.edu/mechanical/thermo/property_tables/air/air_Cp_Cv.html
    cp = 1005
    cv = 718
    gamma = 1.4 # adiabatic exponent cp/cv
#    R = 287.058 # specific gas constant
    R = cp - cv # specific gas constant
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
    
    ## by ID, see triangle.dgf file, given in conservative form
    boundary = {1: lambda t, x, U: U, ## right
                2: lambda t, x, U: [1, 1, 0, 1], ## upper 
                3: lambda t, x, U: [1, 1, 0, 1], ## left
                4: lambda t, x, U: [1, 1, 0, 1], ## lower left
                5: lambda t, x, U: [3, 1, 0, 1], ## lower centre
                6: lambda t, x, U: [1, 1, 0, 1]  ## lower right
                } 

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
Implicit time integration?
Change domain/grid
input of boundary data
Fluid solver takes heat flux, i.e., energy flux as input and outputs boundary temperature
"""
    
class EulerSolver:
    ## get some more inputs, e.g. various orders
    def __init__(self, order_space = 2, order_time = 2):
        self.x = ufl.SpatialCoordinate(ufl.triangle)
        self.u0 = ufl.as_vector([1, 1, 0, 1])
        
        self.domain = (reader.dgf, "triangle.dgf") ## read domain specifications from file
        self.gridView = simplexGrid(self.domain, dimgrid = 2) ## grid with triangles, using alugrid module
        self.space = finiteVolume(self.gridView, dimRange = 4) ## solution space
        self.uh = self.space.interpolate(self.u0, name = "solution")
        
        self.local_uh = self.uh.localFunction()
        
        ## this one causes the crash upon deconstruction of sorts?
        ## current error, something is trying to be freed, that was not allocated by malloc.
        self.rho_gf = gridFunction(self.space.grid, name = "rho", order = order_space)(self.rho)
        
        self.Model = Model()
        self.operator = femDGOperator(self.Model, self.space, limiter = "MinMod") # note: that u_h.space fails since not ufl_space
        ## does not seem to work, resp. take way too long for implicit?
        self.stepper = femdgStepper(order = order_time, operator = self.operator, rkType = 'EX')
        self.operator.applyLimiter(self.uh)
        
        ## flux computation
        ## can one write it like this, i.e., the self.Model.temp(self.uh) ? 
#        self.normal = Constant((0., -1.))
#        self.flux_grad = expression2GF(self.uh.space.grid, ufl.dot(ufl.grad(self.Model.temp(self.uh)), self.normal), self.uh.space.order)
        self.normal = Constant((1., 0., 0., 0.))
        self.flux_grad = expression2GF(self.uh.space.grid, ufl.dot(self.uh, self.normal), self.uh.space.order)
        self.NN = 32 # number of sampling points?
        ## always all zeros?
        self.flux_f = lambda x: lineSample(x, [-1., 0.], [1., 0.], self.NN)[1] ## flux on the horizontal middle line
#        self.flux_f = lambda x: lineSample(x, [-1., -1.], [1., -1.], self.NN)[1] ## flux at bottom boundary
        
    def rho(self, e, x):
        self.local_uh.bind(e)
        return self.local_uh(x)[0]
    
    def get_flux(self):
        return self.flux_f(self.flux_grad)
        
    ## variable "i" is for plotting purposes only
    def do_step(self, t, dt, i):
        self.operator.setTime(t)
        ## this doesn't do anything, unphysical values and no complain
        self.Model.boundary[1] = lambda t, x, U: [-10, 0, 0, -10]
        self.t += self.stepper(self.uh, dt = dt)
        self.gridView.writeVTK('solution', pointdata={'solution': self.uh}, number = i)
        
#        print(dir(self.operator.models[0]))
#        print(dir(self.operator.models))
#        assert False
            
    def solve(self, tf, n):
        dt = tf/n
        self.t = 0
        for i in range(n):
            print('doing step', self.t, dt, flush = True)
            self.do_step(self.t, dt, i)
#            print(self.get_flux())
        
#    def __del__(self):
#        print('deconstructor')
#        pass
            
if __name__ == '__main__':
    tf = 0.6
    n = 200
    solver = EulerSolver()
    solver.solve(tf, n)
    
    print('all done')