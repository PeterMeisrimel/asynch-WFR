/*
Authors: Peter Meisrimel
December 2018
*/

#ifndef REF_SOL_H_
#define REF_SOL_H_

#include "dolfin.h"
#include "math.h"

// Reference solution for the given initial conditions, geometry and constant parameters alpha and lambda
class ReferenceSolution : public dolfin::Expression{
private:
    double alpha, lambda, cons, t_end;
public:
    ReferenceSolution(double a, double g, double c, double t){
        alpha = a;
        lambda = g;
        cons  = c;
        t_end = t;
    }

    void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const{
        values[0]= cons*sin(M_PI*x[0]/2)*sin(M_PI*x[1])*exp(-5.0/4.0*lambda/alpha*M_PI*M_PI*t_end);
    }
};

#endif //REF_SOL_H_
