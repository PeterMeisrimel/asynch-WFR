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

//using namespace dolfin;

class Problem_heat : public Problem{
protected:
    double _dx, _bx; // parameters for building initial condition
    std::shared_ptr<dolfin::Constant> _dt;
    std::shared_ptr<dolfin::Constant> _alpha, _lambda; // equation parameters
    
    std::shared_ptr<dolfin::UnitSquareMesh> _mesh;
    std::shared_ptr<dolfin::FunctionSpace> _V; // FunctionSpace
    
    // boundary stuff
    std::shared_ptr<dolfin::Constant> _dirichlet_boundary_val;
    const dolfin::DirichletBC * _BC;
    std::vector<const dolfin::DirichletBC *> _bcs;
      
    std::shared_ptr<dolfin::Function> _uold;
    std::shared_ptr<dolfin::Function> _ucheckpoint;
    //Function * _unew;
    std::shared_ptr<dolfin::Function> _unew;
    std::shared_ptr<dolfin::Function> _rhs_fnew;
    
    heat::BilinearForm * _a;
    heat::LinearForm   * _L;
    
    double * _uother_old, * _uother_new;
public:
    // gridsize, a, g, (bx, which_u0), last 3 are for initial conditions
    // gridsize = number of internal unknowns on interface
    Problem_heat(int, double, double, double, int);
    
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
