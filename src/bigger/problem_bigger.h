/*
Authors: Peter Meisrimel
December 2018
*/

#ifndef PROBLEM_BIGGER_H_
#define PROBLEM_BIGGER_H_

#include "problem.h"

/*
Underlying system of ODEs:

u' = A u, u(0) = u0, where

A = 
[-1  1  0  0  0,
  0 -1  1  0  0,
  0  0 -1  1  0,
  0  0  0 -1  1,
  1  0  0  0  -1]

u0 = [1, 2, 3, 4, 5]^T

Part 1 takes the first three variables and Part 2 the latter 2
Time-integration is done using the explicit Euler scheme
*/

class Problem_bigger : public Problem{
protected:
  double *_uold, *_unew, *_ucheck;
  Waveform * _WF; // pointer to Waveform, for internal handling in rhs
  virtual void rhs(double, double*, double*) = 0;
public:
  Problem_bigger(int);

  void do_step(double, double, double*, Waveform *);
  void create_checkpoint();
  void reset_to_checkpoint();
};

class Problem_bigger_part_1: public Problem_bigger{
private:
  void rhs(double, double*, double*);
public:
  Problem_bigger_part_1();
};

class Problem_bigger_part_2: public Problem_bigger{
private:
  void rhs(double, double*, double*);
public:
  Problem_bigger_part_2();
};

#endif //PROBLEM_BIGGER_H_
