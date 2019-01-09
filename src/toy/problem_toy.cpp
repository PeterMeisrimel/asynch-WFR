/*
Authors: Peter Meisrimel
December 2018
*/

#include "problem_toy.h"

/////////////////// Problem_toy

Problem_toy::Problem_toy(){
    other_init_done = false;
    _length  = 1;
    _u0     = new double[1];
    _uold   = new double[1];
    _unew   = new double[1];
    _ucheck = new double[1];
}

void Problem_toy::do_step(double t, double dt, double * uout, Waveform * WF_in){
    _WF = WF_in;
    // simple explicit Euler step
    rhs(t, _uold, _unew);
    _unew[0] = _uold[0] + dt*_unew[0];
    uout[0]  = _unew[0];
    _uold[0] = _unew[0];
}

void Problem_toy::create_checkpoint(){
    _ucheck[0] = _uold[0];
}

void Problem_toy::reset_to_checkpoint(){
    _uold[0] = _ucheck[0];
}

/////////////////// Part 1

Problem_toy_part_1::Problem_toy_part_1() : Problem_toy(){
    _u0[0]   = 1.;
    _uold[0] = _u0[0];
    _unew[0] = _u0[0];
}

void Problem_toy_part_1::rhs(double t, double *uold, double *uout){
    _WF->eval(t, _uother);
    uout[0] = -uold[0] + 0.5*_uother[0];
}

/////////////////// Part 2

Problem_toy_part_2::Problem_toy_part_2(){
    _u0[0]   = 3.;
    _uold[0] = _u0[0];
    _unew[0] = _u0[0];
}

void Problem_toy_part_2::rhs(double t, double *uold, double *uout){
    _WF->eval(t, _uother);
    uout[0] = -uold[0] - 0.5*_uother[0];
}
