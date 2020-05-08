/*
Authors: Peter Meisrimel, Benjamin Rueth
September 2018
*/

#include "problem_heat_D.h"
#include "problem_heat.h"

Problem_heat_D::Problem_heat_D(int gridsize, double a, double g, int which) : Problem_heat(gridsize, a, g, 0.0, which){
    other_init_done = false;
    _u0 = new double[_length];

    _dirichlet_boundary_val  = std::make_shared<dolfin::Constant>(0.0);
    _dirichlet_boundary      = std::make_shared<Boundaries_D>();
    _dirichlet_boundary_inft = std::make_shared<Interface_D>();
    _BC = new dolfin::DirichletBC(_V, _dirichlet_boundary_val, _dirichlet_boundary);
    
    _bcs.push_back(_BC);

    _rhs_fnew = std::make_shared<dolfin::Function>(_V);
    dolfin::Constant const_0(0.0);
    _rhs_fnew -> interpolate(const_0);
    _L -> fnew = _rhs_fnew;
    
    _unew->vector() = _uold->vector();

    _heat_F = new heat_flux::LinearForm(_V);
    _heat_F -> alpha = _alpha;
    _heat_F -> lambda_diff = _lambda;
    _heat_F -> dt = _dt;
    _heat_F -> uold = _uold;
    _heat_F -> unew = _unew;
    
    _flux_function = std::make_shared<dolfin::Function>(_V);
    _flux_vec = new dolfin::Vector(MPI_COMM_SELF, _length*_length);
    _flux_vec -> zero();
    
    get_flux(1., _u0);
}

void Problem_heat_D::init_other(int len_other){
    if(not other_init_done){
        _interface_vals = std::make_shared<InterpolatedExpression>(len_other);
    _BC_inft = new dolfin::DirichletBC(_V, _interface_vals, _dirichlet_boundary_inft);
    _bcs.push_back(_BC_inft);
    _length_other = len_other;
    _uother_old = new double[len_other];
    other_init_done = true;
    }
};

// compute fluxes
void Problem_heat_D::get_flux(double dt, double * flux_out){
    *_dt = dt;
    dolfin::assemble(*_flux_vec, *_heat_F);
    *(_flux_function -> vector()) = *_flux_vec;
    for(int i = 1; i < _length - 1; i++)
        flux_out[i] = -(*_flux_function)(1., i*_dx); // negative flux
    flux_out[0] = 0.; //enforcing boundary manually
    flux_out[_length - 1] = 0.;
}

void Problem_heat_D::do_step(double t, double dt, double * unew, Waveform * WF_in){
    WF_in -> eval(t + dt, _uother_old);    
    
    _interface_vals -> set_vals(get_uother_old_p());
    _L -> u0 = _uold;
    *_dt = dt;
    dolfin::solve(*_a == *_L, *_unew, _bcs);
    
    get_flux(dt, unew);
    *_uold->vector() = *(_unew -> vector());
}