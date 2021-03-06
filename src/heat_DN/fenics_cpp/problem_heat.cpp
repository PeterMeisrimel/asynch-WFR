/*
Authors: Peter Meisrimel, Benjamin Rueth
September 2018
*/

#ifndef PROBLEM_HEAT_CPP_
#define PROBLEM_HEAT_CPP_

#include "problem_heat.h"
#include "dolfin.h"
#include "mpi.h"
#include "heat.h" // generated from heat.cfl

Problem_heat::Problem_heat(int gridsize, double a, double g, double bx, int which_u0){
    _length = gridsize + 2; // number of "unknowns", including boundaries
    _dx    = 1/float(gridsize + 1);
    _dxinv = float(gridsize + 1);
    _bx    = bx;
    
    _dt    = std::make_shared<dolfin::Constant>(0);
    _alpha = std::make_shared<dolfin::Constant>(a);
    _lambda = std::make_shared<dolfin::Constant>(g);
    
    _mesh = std::make_shared<dolfin::UnitSquareMesh>(MPI_COMM_SELF, gridsize + 1, gridsize + 1);
    _V    = std::make_shared<heat::FunctionSpace>(_mesh);
    
    _uold = std::make_shared<dolfin::Function>(_V);
    _unew = std::make_shared<dolfin::Function>(_V);
    _ucheckpoint = std::make_shared<dolfin::Function>(_V);
    
    _a = new heat::BilinearForm(_V, _V);
    _a -> dt          = _dt;
    _a -> alpha       = _alpha;
    _a -> lambda_diff = _lambda;
    
    _L = new heat::LinearForm(_V);
    _L -> dt          = _dt;
    _L -> alpha       = _alpha;
    _L -> lambda_diff = _lambda;
    _L -> u0 = _uold;
    
    InitialConditions u_init(_bx, which_u0);
    _uold -> interpolate(u_init);
    _unew -> interpolate(u_init);
}

#endif //PROBLEM_HEAT_CPP_