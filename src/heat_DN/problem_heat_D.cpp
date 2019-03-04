/*
Authors: Peter Meisrimel, Benjamin Rueth
September 2018
*/

#include "problem_heat_D.h"
#include "problem_heat.h"
#include <algorithm>    // std::min, std::max
#include <iomanip> // set precision

Problem_heat_D::Problem_heat_D(int gridsize, double a, double g, double const_c) : Problem_heat(gridsize, a, g, const_c, 0.0, 0.0){
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

    // setting up everything for computing fluxes
    _unew_flux = std::make_shared<Function>(_V);

    _heat_flux_form = new heat_flux::LinearForm(_V);
    _heat_flux_form -> dt    = _dt;
    _heat_flux_form -> alpha = _alpha;
    _heat_flux_form -> gamma = _gamma;
    _heat_flux_form -> uold  = _uold;
    _heat_flux_form -> unew  = _unew_flux;

    _unew->vector() = _uold->vector();

    _flux_function = new Function(_V);
    _flux_vec = new Vector(MPI_COMM_SELF, _length*_length);
    _flux_vec -> zero();

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
        flux_out[i] = (*_gamma)*((*_uold)(1, i*_dx) - (*_uold)(1-_dx, i*_dx))/_dx;
}

// implicit Euler
void Problem_heat_D::do_step(double t, double dt, double * unew, Waveform * WF_in){
    // dirichlet boundary condition linear in time, thus use midpoint for CN
    WF_in -> eval(t + 0.5*dt, _uother_old);
    _interface_vals -> set_vals(get_uother_old_p());
	_L -> u0 = _uold;
	*_dt = dt;
	solve(*_a == *_L, *_unew, _bcs);

    get_flux(dt, unew);

    *_uold->vector() = *(_unew -> vector());
}
