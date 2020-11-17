/*
Authors: Peter Meisrimel
June 2020
*/

#ifndef PROBLEM_HEAT_DUNE_H_
#define PROBLEM_HEAT_DUNE_H_

#include "Python.h"
#include "waveform.h"
#include "problem.h"
#include "iostream"
//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION // disables from warning from outdated API
//#include "numpy/arrayobject.h"

class Problem_heat_dune : public Problem{
protected:
    PyObject *pClass_obj, *Pdt;
    double _dx;
public:
    Problem_heat_dune(){};
    Problem_heat_dune(int, double, double, const char*, const char*);
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

class Problem_heat_dune_D : public Problem_heat_dune{
public:
    Problem_heat_dune_D(){};
    Problem_heat_dune_D(int gridsize, double a, double g,
                        const char* py_file = "heat_dune_grad", const char* py_class = "Problem_heat_D_weak")
                        :Problem_heat_dune(gridsize, a, g, py_file, py_class){};
    void do_step(double, double, double *, Waveform *);
};

class Problem_heat_dune_N : public Problem_heat_dune{
public:
    Problem_heat_dune_N(int gridsize, double a, double g,
                        const char* py_file = "heat_dune_grad", const char* py_class = "Problem_heat_N")
                        :Problem_heat_dune(gridsize, a, g, py_file, py_class){};
    void do_step(double, double, double *, Waveform *);
};

class Problem_heat_dune_N_euler : public Problem_heat_dune_N{
public:
    Problem_heat_dune_N_euler(int gridsize, double a, double g,
                              const char* py_file = "heat_dune_grad", const char* py_class = "Problem_heat_N_euler")
                              :Problem_heat_dune_N(gridsize, a, g, py_file, py_class){};
};

class Problem_heat_dune_euler_fluid : public Problem_heat_dune_D{
protected:
    PyObject *Pt;
public:
    Problem_heat_dune_euler_fluid(int);
    void do_step(double, double, double *, Waveform *);
};

#endif //PROBLEM_HEAT_DUNE_H_