/*
Authors: Peter Meisrimel, Benjamin Rueth
September 2018
*/

#include "problem_heat_D.h"
#include "problem_heat.h"

Problem_heat_D::Problem_heat_D(int gridsize, double a, double g, int which) : Problem_heat(gridsize, a, g, 0.0, 0.0, which){
    other_init_done = false;
    _length = _N + 1;
    _u0 = new double[_length];

    _dirichlet_boundary_val  = std::make_shared<Constant>(0.0);
    _dirichlet_boundary      = std::make_shared<Boundaries_D>();
    _dirichlet_boundary_inft = std::make_shared<Interface_D>();
    _BC = new DirichletBC(_V, _dirichlet_boundary_val, _dirichlet_boundary);

    _bcs.push_back(_BC);

    _rhs_fold = std::make_shared<Function>(_V);
    _rhs_fnew = std::make_shared<Function>(_V);
    Constant const_0(0.0);
    _rhs_fold -> interpolate(const_0);
    _rhs_fnew -> interpolate(const_0);
    _L -> fold = _rhs_fold;
    _L -> fnew = _rhs_fnew;

    _unew->vector() = _uold->vector();

    get_flux(1., _u0);
}

void Problem_heat_D::init_other(int len_other){
    if(not other_init_done){
        _interface_vals = std::make_shared<InterpolatedExpression>(len_other);
		_BC_inft = new DirichletBC(_V, _interface_vals, _dirichlet_boundary_inft);
		_bcs.push_back(_BC_inft);
        _length_other = len_other;
        _uother_old = new double[len_other];
        other_init_done = true;
    }
};

// compute fluxes
void Problem_heat_D::get_flux(double dt, double * flux_out){
    for(int i = 0; i < _length; i++)
        flux_out[i] = (*_lambda)*((*_uold)(1, i*_dx) - (*_uold)(1-_dx, i*_dx))/_dx;
}

void Problem_heat_D::do_step(double t, double dt, double * unew, Waveform * WF_in){
    WF_in -> eval(t + dt, _uother_old);    

    _interface_vals -> set_vals(get_uother_old_p());
	_L -> u0 = _uold;
	*_dt = dt;
	solve(*_a == *_L, *_unew, _bcs);

    *_uold->vector() = *(_unew -> vector());
    get_flux(dt, unew);
}
