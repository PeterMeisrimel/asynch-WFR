/*
Authors: Peter Meisrimel, Benjamin Rueth
September 2018
*/

#ifndef PROBLEM_HEAT_D_H_
#define PROBLEM_HEAT_D_H_

#include "waveform.h"
#include "problem_heat.h"
#include "boundaries.h"
#include "int_expression.h"
#include "dolfin.h"
#include "mpi.h"
#include "heat.h" // generated from heat.cfl
#include "heat_flux.h" // generated from heat_flux.cfl

// Dirichlet part of a coupled heat equation
// interface on right-hand side 

class Problem_heat_D : public Problem_heat{
private:
    std::shared_ptr<Boundaries_D> _dirichlet_boundary;
    std::shared_ptr<Interface_D> _dirichlet_boundary_inft;
    const dolfin::DirichletBC * _BC_inft;
  
    std::shared_ptr<InterpolatedExpression> _interface_vals;

    void get_flux(double, double*);

    heat_flux::LinearForm * _heat_F;
    std::shared_ptr<dolfin::Function> _flux_function;
    dolfin::Vector * _flux_vec;
public:
    // gridsize, alpha, lambda
    Problem_heat_D(int, double, double, int = 0);

    void init_other(int);
    // t, dt, unew (out), WF
    void do_step(double, double, double *, Waveform *);
};

#endif //PROBLEM_HEAT_D_H_
