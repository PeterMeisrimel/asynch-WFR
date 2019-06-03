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
//#include "heat_flux.h" // generated from heat_flux.cfl

using namespace dolfin;

// Dirichlet part of a coupled heat equation
// interface on right-hand side 

class Problem_heat_D : public Problem_heat{
private:
    std::shared_ptr<Constant> _dirichlet_boundary_val;
    std::shared_ptr<Boundaries_D> _dirichlet_boundary;
    std::shared_ptr<Interface_D> _dirichlet_boundary_inft;
    const DirichletBC * _BC;
    const DirichletBC * _BC_inft;
    std::vector<const DirichletBC *> _bcs;
  
    std::shared_ptr<Function> _unew_flux;
    Function * _flux_function;
    std::shared_ptr<InterpolatedExpression> _interface_vals;

    //heat_flux::LinearForm * _heat_flux_form;

    Vector * _flux_vec;

    void get_flux(double, double*);
public:
    // gridsize, alpha, lambda, const for u0
    Problem_heat_D(int, double, double, double = 500, int = 0);

    void init_other(int);
    // t, dt, unew (out), WF
    void do_step(double, double, double *, Waveform *);
};

#endif //PROBLEM_HEAT_D_H_
