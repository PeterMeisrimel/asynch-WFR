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

#xx = np.linspace(-1, 1, 100)
#aa = np.random.rand(100)*20
#interp = interp1d(xx, aa)
#def bound(x):
##    print('aaaa', x, type(x), type(x[0]))
##    return float(interp([x[0]])[0])
#    return float(interp([x])[0])
##bound = lambda x: float(interp(x[0]))
#fff = lambda t: (t + 1)**20

#def fxx(t, x):
#    ## this needs to be an ufl expression or similar
#    return (t+1) * ufl.exp(x[0])

#class interpol_exp(ufl.classes.Expr):
#    def __init__(self):
#        pass
#    
#    def evaluate(self, x, a, b, c):
#        return 2
    

#class interpol_exp(ufl.core.expr.Expr):
#@ufl.core.ufl_type.Expr()
#class interpol_exp(ufl.classes.Operator):
#    def __init__(self, *args):
#        super().__init__(*args)
##        self.ufl_operands = None
#        pass
#    def f(self, *args):
#        return 2
#    
#    def evaluate(self, *args):
#        return 100
#    
#    def _call(self, *args):
#        return 100

#@ufl.core.ufl_type.Expr()
#@ufl.core.ufl_type.ufl_type(is_abstract = False)
##class interpol_exp(ufl.core.expr.Expr):
#class interpol_exp(ufl.core.operator.Operator):
#    _ufl_noslots_ = True
#    ufl_shape = (1,)
##    def __init__(self, *args):
##        super().__init__(*args)
###        self.ufl_operands = None
##        pass
#    def evaluate(self, x, mapping, component, index_values):
#        a = self.ufl_operands[0].evaluate(x, mapping, component, index_values)
#        return 2
import math

@ufl.core.ufl_type.ufl_type()
class interpol_exp(ufl.mathfunctions.MathFunction):
    __slots__ = ()
    
    def set_vals(self, xx, yy):
        self.interp = interp1d(xx.copy(), yy.copy())

## cls functions have access to class parameters, but not parameters specific to the instance
#    def __new__(cls, argument):
    def __new__(self, argument):
#        print(argument)
#        print(type(argument))
        return ufl.constantvalue.FloatValue(math.exp(float(argument)))

    def __init__(self, argument):
        ufl.mathfunctions.MathFunction.__init__(self, "exp", argument)

#@ufl.core.ufl_type.ufl_type()
#class interpol_exp(ufl.core.expr.Expr):
#    _ufl_noslots_ = True
#    is_terminal = True
#    def evaluate(self, x, mapping, component, index_values):
#        print('evaluation!!!')
#        return 1
    
#@ufl.core.ufl_type
#class interpol_exp(ufl.classes.Operator):
#    def __init__(self, *args):
#        print(*args)
#        super().__init__(*args)
#        self.ufl_operands = None
##        pass
#    def f(self, *args):
#        return 2
#    
#    def evaluate(self, t, x, *args):
#        return 2

class Model:
    def __init__(self):
        ## material parameters
        ## source https://www.ohio.edu/mechanical/thermo/property_tables/air/air_Cp_Cv.html
        self.cp, self.cv = 1005, 718
        self.gamma = 1.4
        self.R = 287.058 ## R = cp - cv ## specific gas constant
    
        ## initial conditions
        self.rho0 = 1.225 # kg/m**3
        self.T0 = 293.15 # 20 C
        self.T1 = self.T0 + 400
        self.vx0, self.vy0 = 1, 0 # m/s
        
        self.E0 = self.rho0 * (self.R*self.T0*(self.gamma - 1) + 0.5 * (self.vx0**2 + self.vy0**2))
        self.E1 = self.rho0 * (self.R*self.T1*(self.gamma - 1) + 0.5 * (self.vx0**2 + self.vy0**2))
        ## something wrong about these boundary conditions/initial conditions
        ## sends quite a wave (not shockwave though) upwards at beginning
        
        ## by ID, see triangle.dgf file
        # given in conservative form
        self.boundary = {1: lambda t, x, U: U, ## right
                         2: lambda t, x, U: [self.rho0, self.vx0, self.vy0, self.E0], ## upper 
                         3: lambda t, x, U: [self.rho0, self.vx0, self.vy0, self.E0], ## left
                         4: lambda t, x, U: [U[0], U[1], 0,
                                             U[0]*(self.R*(self.gamma-1)*self.T1 + 0.5*(U[1]/U[0])**2)], ## lower left
                         5: lambda t, x, U: [U[0], U[1], 0, 
                                             U[0]*(self.R*(self.gamma-1)*self.T1 + 0.5*(U[1]/U[0])**2)], ## lower centre
                         6: lambda t, x, U: [U[0], U[1], 0, 
                                             U[0]*(self.R*(self.gamma-1)*self.T1 + 0.5*(U[1]/U[0])**2)]  ## lower right
                         } 
        
        self.eps = 1e-13
        self.x = ufl.SpatialCoordinate(ufl.triangle)
        # ufl.conditional(cond, True, False)
        ## it doesn't like the ufl.And?
#        self.u0 = ufl.conditional(ufl.And(ufl.And(ufl.le(self.x[1], -1 + self.eps),
#                                  ufl.ge(self.x[0], -0.5 - self.eps)),
#                                  ufl.le(self.x[0], 0.5 + self.eps)), ## condition for lower boundary
#                                  ufl.as_vector([self.rho0, self.vx0, self.vy0, self.E0]), 
#                                  ufl.as_vector([self.rho0, self.vx0, self.vy0, self.E1]))
        self.u0 = ufl.conditional(ufl.le(self.x[1], -1 + self.eps), ## condition for lower boundary
                          ufl.as_vector([self.rho0, self.vx0, self.vy0, self.E1]), 
                          ufl.as_vector([self.rho0, self.vx0, self.vy0, self.E0]))
                
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
    
    def set_temp(self, T, U):
        rho, v, p = self.toPrim(U)
        return ufl.as_vector([rho, rho*v[0], 0, rho*(self.R*(self.gamma - 1)*T + 0.5*v[0]*v[0])])
    
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
        self.Model = Model()
        self.u0 = self.Model.u0
        
        self.domain = (reader.dgf, "triangle.dgf") ## read domain specifications from file
        self.gridView = simplexGrid(self.domain, dimgrid = 2) ## grid with triangles, using alugrid module
        self.space = finiteVolume(self.gridView, dimRange = 4) ## solution space
        self.uh = self.space.interpolate(self.u0, name = "solution")
        
        self.uh_bd = self.space.interpolate(Constant((1.225, 1., 0., self.Model.E0)), name = "bd")
        
        for i in range(1, 7):
            ## this seems to work as a start? 
            ## next steps:
            # get values at boundary, I guess on a line, do temperature modification, interpolate from pyhton function to function space
            self.Model.boundary[i] = lambda t, x, U: self.u0
        
        self.local_uh = self.uh.localFunction()
        
        ## this one causes the crash upon deconstruction of sorts?
        ## current error, something is trying to be freed, that was not allocated by malloc.
        self.rho_gf = gridFunction(self.space.grid, name = "rho", order = order_space)(self.rho)
        
#        yy = np.linspace(-0.5, 0.5, 20)
#        zz = np.array(range(1, 21))
#        ii = interp1d(yy, zz)
        
#        self.Model = Model()
        
#        x = ufl.SpatialCoordinate(ufl.triangle)
#        aa = interpol_exp(1)
#        aa = interpol_exp()
#        self.Model.boundary[5] = lambda t, x, U, n: [ii(x[0]), 1, 0, 1]
#        self.Model.boundary[5] = lambda t, x, U, *args: ufl.as_vector([aa.evaluate, 1, 0, 1])
#        self.Model.boundary[5] = lambda t, x, U, *args: self.Model.set_temp(10000, U)
#        self.Model.boundary[5] = lambda t, x, U, *args: [3, 1, 0, 1]
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
#        self.Model.boundary[1] = lambda t, x, U: [-10, 0, 0, -10]
        self.t += self.stepper(self.uh, dt = dt)
#        self.gridView.writeVTK('solution', pointdata={'solution': self.uh, 'T': self.Model.temp(self.uh)}, number = i)
        
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
#    tf = 1.2
#    n = 400
    tf = 0.05
    n = 2000
    solver = EulerSolver()
    solver.solve(tf, n)
    
    print('all done')