/*
Authors: Peter Meisrimel, Benjamin Rueth
September 2018
*/

#include "problem_heat_N.h"
#include "problem_heat.h"

Problem_heat_N::Problem_heat_N(int gridsize, double a, double g, int which) : Problem_heat(gridsize, a, g, 1.0, which){
    other_init_done = false;
    _u0 = new double[_length];

    for(int i = 1; i < _length - 1; i++)
        _u0[i] = (*_uold)(0., i*_dx);
    _u0[0] = 0.; // manually enforce boundaries
    _u0[_length - 1] = 0.;
    
    _dirichlet_boundary_val = std::make_shared<dolfin::Constant>(0.0);
    _dirichlet_boundary     = std::make_shared<Boundaries_N>();
    _BC = new dolfin::DirichletBC(_V, _dirichlet_boundary_val, _dirichlet_boundary);

    _bcs.push_back(_BC);

    _rhs_fold = std::make_shared<dolfin::Function>(_V);
    _rhs_fnew = std::make_shared<dolfin::Function>(_V);
    _L -> fold = _rhs_fold;
    _L -> fnew = _rhs_fnew;
}

void Problem_heat_N::init_other(int len_other){
    if(not other_init_done){
        _fluxx = std::make_shared<InterpolatedExpression>(len_other);
        _length_other = len_other;
        _uother_old = new double[len_other];
        _uother_new = new double[len_other];
        other_init_done = true;
    }
}

void Problem_heat_N::do_step(double t, double dt, double * u_output, Waveform * WF_in){
    *_dt = dt;
    WF_in -> eval(t, _uother_old);
//    std::cout << " N eval at t = " << t << std::endl;
//    for (int i = 0; i < 32; i++)
//        std::cout << _uother_old[i] << " ";
//    std::cout << std::endl;
    _fluxx -> set_vals(get_uother_old_p());
    _rhs_fold -> interpolate(*_fluxx);
    
    WF_in -> eval(t + dt, _uother_new);
//    std::cout << " N eval at t = " << t + dt << std::endl;
//    for (int i = 0; i < 32; i++)
//        std::cout << _uother_new[i] << " ";
//    std::cout << std::endl;
    _fluxx -> set_vals(get_uother_new_p());
    _rhs_fnew -> interpolate(*_fluxx);

	dolfin::solve(*_a == *_L, *_unew, _bcs);
    *_uold->vector() = *(_unew -> vector());

    // evalute solution and return boundary values
    for(int i = 1; i < _length - 1; i++)
        u_output[i] = (*_uold)(0, i*_dx);
    u_output[0] = 0.; // manually enforce boundaries
    u_output[_length - 1] = 0.;
}
