/*
Authors: Peter Meisrimel
March 2020
*/

#ifndef PROBLEM_HEAT_PY_H_
#define PROBLEM_HEAT_PY_H_

#include "Python.h"
#include "waveform.h"
#include "problem.h"
#include "iostream"
//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION // disables from warning from outdated API
//#include "numpy/arrayobject.h"

class Problem_heat_python : public Problem{
protected:
    PyObject *pClass_obj, *Pdt;
    double _dx;
public:
    Problem_heat_python(int, double, double, const char*);
    void create_checkpoint();
    void reset_to_checkpoint();
    void get_u0(double*);
    double get_norm_factor(){
        return std::sqrt(_dx);
    }
    void init_other(int len_other){
        _length_other = len_other;
    };
};

class Problem_heat_python_D : public Problem_heat_python{
public:
    Problem_heat_python_D(int gridsize, double a, double g):Problem_heat_python(gridsize, a, g, "Problem_heat_D_weak"){};
    void do_step(double, double, double *, Waveform *);
};

class Problem_heat_python_N : public Problem_heat_python{
public:
    Problem_heat_python_N(int gridsize, double a, double g) :Problem_heat_python(gridsize, a, g, "Problem_heat_N"){};
    void do_step(double, double, double *, Waveform *);
};

#endif //PROBLEM_HEAT_PY_H_