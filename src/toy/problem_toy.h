/*
Authors: Peter Meisrimel
December 2018
*/

#ifndef PROBLEM_TOY_H_
#define PROBLEM_TOY_H_

#include "problem.h"

/*
Underlying system of ODEs:

u' = A u, u(0) = u0, where

A = 
[-1  0.5,
 -0.5 -1]

u0 = [1, 3]^T

Part 1 takes the first variable, Part 2 the second one
Time-integration is done using the explicit Euler scheme
*/


class Problem_toy : public Problem{
protected:
    double *_uold, *_unew, *_ucheck;
    Waveform * _WF; // pointer to Waveform, for internal handling in rhs
    virtual void rhs(double, double*, double*) = 0;
public:
    Problem_toy();

    void do_step(double, double, double*, Waveform *);
    void create_checkpoint();
    void reset_to_checkpoint();
};

class Problem_toy_part_1: public Problem_toy{
private:
    void rhs(double, double*, double*);
public:
    Problem_toy_part_1();
};

class Problem_toy_part_2: public Problem_toy{
private:
    void rhs(double, double*, double*);
public:
    Problem_toy_part_2();
};

#endif //PROBLEM_TOY_H_
