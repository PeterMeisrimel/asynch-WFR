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
  std::shared_ptr<Constant> _dirichlet_boundary_val;
  std::shared_ptr<Boundaries_N> _dirichlet_boundary;
  const DirichletBC * _BC;
  std::vector<const DirichletBC *> _bcs;
  
  //std::shared_ptr<Function> _rhs_f;
  std::shared_ptr<InterpolatedExpression> _fluxx;
public:
  // gridsize, alpha, gamma, const for u0
  Problem_heat_N(int, double, double, double = 500);

  void init_other(int);

  // t, dt, unew (out), WF
  void do_step(double, double, double *, Waveform *);
};

#endif // PROBLEM_HEAT_N_H_
