/*
Authors: Peter Meisrimel, Benjamin Rueth
September 2018
*/

#include "problem_heat_N.h"
#include "problem_heat.h"

Problem_heat_N::Problem_heat_N(int gridsize, double a, double g, int which) : Problem_heat(gridsize, a, g, 1.0, 0.0, which){
    other_init_done = false;
    _length = _N + 1;
    _u0 = new double[_length];

    for(int i = 0; i < _length; i++)
        _u0[i] = (*_uold)(0, i*_dx);

    _dirichlet_boundary_val = std::make_shared<Constant>(0.0);
    _dirichlet_boundary     = std::make_shared<Boundaries_N>();
    _BC = new DirichletBC(_V, _dirichlet_boundary_val, _dirichlet_boundary);

    _bcs.push_back(_BC);

    _rhs_fold = std::make_shared<Function>(_V);
    _rhs_fnew = std::make_shared<Function>(_V);
    _L -> fold = _rhs_fold;
    _L -> fnew = _rhs_fnew;
}

void Problem_heat_N::init_other(int len_other){
    if(not other_init_done){
        _fluxx_old = std::make_shared<InterpolatedExpression>(len_other);
        _fluxx_new = std::make_shared<InterpolatedExpression>(len_other);
        _length_other = len_other;
        _uother_old = new double[len_other];
        _uother_new = new double[len_other];
        other_init_done = true;
    }
}

void Problem_heat_N::do_step(double t, double dt, double * unew, Waveform * WF_in){
    WF_in -> eval(t, _uother_old);
    // force boundary points of flux to zero, coefficients and sign have been set accordingly in heat.ufl
	_uother_old[0] = 0;
	_uother_old[_length_other-1] = 0;
    _fluxx_old -> set_vals(get_uother_old_p());
    _rhs_fold -> interpolate(*_fluxx_old);

    WF_in -> eval(t + dt, _uother_new);
	_uother_new[0] = 0;
	_uother_new[_length_other-1] = 0;
    _fluxx_new -> set_vals(get_uother_new_p());
    _rhs_fnew -> interpolate(*_fluxx_new);

	_L -> u0 = _uold;
	*_dt = dt;
	solve(*_a == *_L, *_unew, _bcs);

    *_uold->vector() = *(_unew -> vector());

    // evalute solution and return boundary values
    for(int i = 0; i < _length; i++)
        unew[i] = (*_uold)(0, i*_dx);
}
