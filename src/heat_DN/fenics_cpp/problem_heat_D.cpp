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

    _rhs_fold = std::make_shared<dolfin::Function>(_V);
    _rhs_fnew = std::make_shared<dolfin::Function>(_V);
    _rhs_fnew -> interpolate(dolfin::Constant(0.0));
    _rhs_fold -> interpolate(dolfin::Constant(0.0));
    _L -> fold = _rhs_fold;
    _L -> fnew = _rhs_fnew;
    
    _unew->vector() = _uold->vector();

    _heat_Flux_grad = new heat_flux_grad::LinearForm(_V);
    _heat_Flux_grad -> lambda_diff = _lambda;
    _heat_Flux_grad -> unew = _unew;
    _heat_Flux_grad -> n = std::make_shared<dolfin::Constant>(1.0, 0.0);
    
    _heat_Flux_weak = new heat_flux_weak::LinearForm(_V);
    _heat_Flux_weak -> alpha = _alpha;
    _heat_Flux_weak -> lambda_diff = _lambda;
    _heat_Flux_weak -> dt = _dt;
    _heat_Flux_weak -> uold = _uold;
    _heat_Flux_weak -> unew = _unew;
    
    _flux_function = std::make_shared<dolfin::Function>(_V);
    _flux_vec = new dolfin::Vector(MPI_COMM_SELF, _length*_length);
    _flux_vec -> zero();
    _flux_vec_old = new double[_length];
    
    // initial flux is computed via gradient
    // for a parallel method we need an intial value and the weak form relies on a difference with previous flux
    dolfin::assemble(*_flux_vec, *_heat_Flux_grad);
    *(_flux_function -> vector()) = *_flux_vec;
    for(int i = 1; i < _length - 1; i++)
        _u0[i] = (*_flux_function)(1., i*_dx)*_dxinv;
    // manually enforce boundary values
    _u0[0] = 0;
    _u0[_length - 1] = 0; // rounding might result in i*_dx to be outside of domain
    
    // store old flux
    for(int i = 0; i < _length; i++)
        _flux_vec_old[i] = _u0[i];
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

void Problem_heat_D::get_flux(double * flux_out){
    dolfin::assemble(*_flux_vec, *_heat_Flux_weak);
    *(_flux_function -> vector()) = *_flux_vec;
    for(int i = 1; i < _length - 1; i++)
        flux_out[i] = 2*(*_flux_function)(1., i*_dx)*_dxinv - _flux_vec_old[i]; // negative flux
    flux_out[0] = 0.; //enforcing boundary manually
    flux_out[_length - 1] = 0.;
    
    for(int i = 1; i < _length - 1; i++) // store old flux
        _flux_vec_old[i] = flux_out[i];
}

void Problem_heat_D::reset_to_checkpoint(){
    // also need to reset initial flux
    *_uold -> vector() = *_ucheckpoint -> vector();
    
    // restore initial flux
    for(int i = 0; i < _length; i++)
        _flux_vec_old[i] = _u0[i];
}

void Problem_heat_D::do_step(double t, double dt, double * flux_out, Waveform * WF_in){
    *_dt = dt;

    WF_in -> eval(t + dt, _uother_old);    
    _interface_vals -> set_vals(get_uother_old_p());
    
    dolfin::solve(*_a == *_L, *_unew, _bcs);
    
    get_flux(flux_out);
    *_uold->vector() = *(_unew -> vector());
}