/*
Authors: Peter Meisrimel, Benjamin Rueth
September 2018
*/

#include "problem_heat_N.h"
#include "problem_heat.h"

Problem_heat_N::Problem_heat_N(int gridsize, double a, double g, double const_c) : Problem_heat(gridsize, a, g, const_c, 1.0, 0.0){
    other_init_done = false;
    _length = _N + 1;
    _u0 = new double[_length];

    for(int i = 0; i < _length; i++)
        _u0[i] = (*_uold)(0, i*_dx);

    _dirichlet_boundary_val = std::make_shared<Constant>(0.0);
    _dirichlet_boundary     = std::make_shared<Boundaries_N>();
    _BC = new DirichletBC(_V, _dirichlet_boundary_val, _dirichlet_boundary);

    _bcs.push_back(_BC);

    _rhs_f = std::make_shared<Function>(_V);
    _L -> f = _rhs_f;
}

void Problem_heat_N::init_other(int len_other){
    if(not other_init_done){
        _fluxx = std::make_shared<InterpolatedExpression>(len_other);
        _length_other = len_other;
        _uother = new double[len_other];
        other_init_done = true;
    }
}

// implicit Euler
void Problem_heat_N::do_step(double t, double dt, double * unew, Waveform * WF_in){
    // t + dt as we do implicit euler
    WF_in -> eval(t + dt, _uother);

    // change sign of flux
    for(int i = 0; i < _length_other; i++)
        _uother[i] *= -1;
	_uother[0] = 0;
	_uother[_length_other-1] = 0;

    _fluxx -> set_vals(get_uother_p());
    _rhs_f -> interpolate(*_fluxx);

	_L -> u0 = _uold;
	*_dt = dt;
	solve(*_a == *_L, *_unew, _bcs);

    *_uold->vector() = *(_unew -> vector());

    // evalute solution and return boundary values
    for(int i = 0; i < _length; i++)
        unew[i] = (*_uold)(0, i*_dx);
}
