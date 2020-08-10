#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 17:37:25 2020

@author: Peter Meisrimel, Lund University
"""

import numpy as np
import pylab as pl
import scipy as sp
pl.close('all')

def get_parameters(which):
    ## parameter selection
    lambda_air = 0.0243
    lambda_water = 0.58
    lambda_steel = 48.9
    
    alpha_air = 1.293*1005
    alpha_water = 999.7*4192.1
    alpha_steel = 7836*443
    if which == 'air_water':
        return {'alpha1': alpha_air, 'alpha2': alpha_water, 'lambda1': lambda_air, 'lambda2': lambda_water}
    elif which == 'air_steel':
        return {'alpha1': alpha_air, 'alpha2': alpha_steel, 'lambda1': lambda_air, 'lambda2': lambda_steel}
    elif which == 'water_steel':
        return {'alpha1': alpha_water, 'alpha2': alpha_steel, 'lambda1': lambda_water, 'lambda2': lambda_steel}
    elif which == 'test':
        return {'alpha1': 1., 'alpha2': 1., 'lambda1': 0.1, 'lambda2': 0.1}
    
from numpy import sin, cos, pi
#############################################################
#  DNWR Relaxation
#############################################################
## see NNWR paper formulas (9.2) (9.3)
## Monge, Azahar, and Philipp Birken.
## "A Multirate Neumann--Neumann Waveform Relaxation Method for Heterogeneous Coupled Heat Equations."
## SIAM Journal on Scientific Computing 41.5 (2019): S86-S105.
class heat_prob:
    def __init__(self, n = 31, alpha1 = 1, alpha2 = 1, lambda1 = 0.1, lambda2 = 0.1):
        self.n = n
        self.alpha1, self.alpha2 = alpha1, alpha2
        self.lambda1, self.lambda2  = lambda1, lambda2
        self.dx = 1/(n+1)
        
    def w_i(self, a, lam, dt):
        dx = self.dx
        return sum([3*dt*dx**2 * sin(i*pi*dx)**2/(2*a*dx**2 + 6*lam*dt + (a*dx**2 - 6*lam*dt) * cos(i*pi*dx)) for i in range(1, self.n + 1)])
    def S_i(self, a, lam, dt):
        dx = self.dx
        return (6*dt*dx*(a*dx**2 + 3*lam*dt) - (a*dx**2 - 6*lam*dt)**2 * self.w_i(a, lam, dt))/(18*dt**2*dx**3)
    
    def GS_theta_opt(self, dt1, dt2):
        dt = max(dt1, dt2)
        S1 = self.S_i(self.alpha1, self.lambda1, dt)
        S2 = self.S_i(self.alpha2, self.lambda2, dt)
        return abs(1/(1 + S1/S2))
    
    def GS_theta_opt2(self, dt1, dt2):
        dt = max(dt1, dt2)
        S1 = self.S_i(self.alpha2, self.lambda2, dt)
        S2 = self.S_i(self.alpha1, self.lambda1, dt)
        return abs(1/(1 + S1/S2))
    
    def JACOBI_theta_opt(self, dt1, dt2):
        dt = max(dt1, dt2)
        S1 = self.S_i(self.alpha1, self.lambda1, dt)
        S2 = self.S_i(self.alpha2, self.lambda2, dt)
        
        if -S1/S2 < 0:
            return 1/((S1/S2)**2 + 1)
        else:
            return 1