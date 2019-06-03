/*
Authors: Peter Meisrimel, Benjamin Rueth
September 2018
*/

#ifndef PROBLEM_HEAT_N_H_
#define PROBLEM_HEAT_N_H_

#include "waveform.h"
#include "problem_heat.h"
#include "boundaries.h"
#include "int_expression.h"
#include <dolfin.h>
#include "mpi.h"
#include "heat.h" // generated from heat.cfl

using namespace dolfin;

// Dirichlet part of a coupled heat equation
// interface on right-hand side 

class Problem_heat_N : public Problem_heat{
private:
    std::shared_ptr<Boundaries_N> _dirichlet_boundary;
  
    std::shared_ptr<InterpolatedExpression> _fluxx_old;
    std::shared_ptr<InterpolatedExpression> _fluxx_new;
public:
    // gridsize, alpha, lambda, const for u0
    Problem_heat_N(int, double, double, double = 500, int = 0);

    void init_other(int);

    // t, dt, unew (out), WF
    void do_step(double, double, double *, Waveform *);
};

#endif // PROBLEM_HEAT_N_H_
