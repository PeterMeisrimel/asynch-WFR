/*
Authors: Peter Meisrimel
March 2020
*/

#ifndef PROBLEM_HEAT_PY_H_
#define PROBLEM_HEAT_PY_H_

#include "waveform.h"
#include "problem.h"
#include "Python.h"
#include "iostream"
//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION // disables from warning from outdated API
//#include "numpy/arrayobject.h"

class Problem_heat_python : public Problem{
protected:
    int _N; // meshsize, number of grid cells in one direction, thus points per side = _N + 1
	PyObject *pClass_obj;
public:
    Problem_heat_python(int, double, double, const char*);
    void create_checkpoint();
    void reset_to_checkpoint();
	void get_u0(double*);
    void init_other(int len_other){
        _length_other = len_other;
    };
};

class Problem_heat_python_D : public Problem_heat_python{
public:
	Problem_heat_python_D(int gridsize, double a, double g):Problem_heat_python(gridsize, a, g, "Problem_heat_D"){};
	void do_step(double, double, double *, Waveform *);
};

class Problem_heat_python_N : public Problem_heat_python{
public:
	Problem_heat_python_N(int gridsize, double a, double g) :Problem_heat_python(gridsize, a, g, "Problem_heat_N"){};
	void do_step(double, double, double *, Waveform *);
};

#endif //PROBLEM_HEAT_PY_H_
