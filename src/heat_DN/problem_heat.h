/*
Authors: Peter Meisrimel, Benjamin Rueth
September 2018
*/

#ifndef PROBLEM_HEAT_H_
#define PROBLEM_HEAT_H_

#include "waveform.h"
#include "problem.h"
#include "u0.h"
#include "dolfin.h"
#include "mpi.h"
#include "heat.h" // generated from heat.cfl

using namespace dolfin;

class Problem_heat : public Problem{
protected:
    int _N; // meshsize, number of grid cells in one direction, thus points per side = _N + 1
    double _dx;
    double _c, _bx, _by; // parameters for building initial condition
    std::shared_ptr<Constant> _dt;
    std::shared_ptr<Constant> _alpha, _gamma; // equation parameters

    std::shared_ptr<UnitSquareMesh> _mesh;
    std::shared_ptr<FunctionSpace> _V; // FunctionSpace
  
    std::shared_ptr<Function> _uold;
    std::shared_ptr<Function> _ucheckpoint;
    Function * _unew;
    std::shared_ptr<Function> _rhs_fold;
    std::shared_ptr<Function> _rhs_fnew;

    heat::BilinearForm * _a;
    heat::LinearForm   * _L;

    double * _uother_old, * _uother_new;

public:
    Problem_heat(int, double, double, double, double, double);

    void create_checkpoint(){
        *_ucheckpoint -> vector() = *_uold -> vector();
    }
    void reset_to_checkpoint(){
        *_uold -> vector() = *_ucheckpoint -> vector();
    }
    double * get_uother_old_p(){return _uother_old;};
    double * get_uother_new_p(){return _uother_new;};
};

#endif //PROBLEM_HEAT_H_
