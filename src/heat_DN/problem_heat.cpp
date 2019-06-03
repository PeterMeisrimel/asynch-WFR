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

using namespace dolfin;

Problem_heat::Problem_heat(int gridsize, double a, double g, double const_c, double bx, double by, int which_u0){
    _N     = gridsize;
    _dx    = 1/float(gridsize);
    _c     = const_c;
    _bx    = bx;
    _by    = by;

    _alpha = std::make_shared<Constant>(a);
    _lambda = std::make_shared<Constant>(g);
    _dt    = std::make_shared<Constant>(0);

    _mesh = std::make_shared<UnitSquareMesh>(MPI_COMM_SELF, _N, _N);
    _V    = std::make_shared<heat::FunctionSpace>(_mesh);

    _a = new heat::BilinearForm(_V, _V);
    _L = new heat::LinearForm(_V);

    _L -> alpha       = _alpha;
    _L -> lambda_diff = _lambda;
    _L -> dt          = _dt;
    _a -> alpha       = _alpha;
    _a -> lambda_diff = _lambda;
    _a -> dt          = _dt;
	
    InitialConditions u_init(const_c, _bx, _by, which_u0);
    _uold        = std::make_shared<Function>(_V);
    _ucheckpoint = std::make_shared<Function>(_V);
    _uold -> interpolate(u_init);
    _unew = new Function(_V);   
}

#endif //PROBLEM_HEAT_CPP_
