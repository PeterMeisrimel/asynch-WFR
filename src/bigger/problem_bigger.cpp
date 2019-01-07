/*
Authors: Peter Meisrimel
December 2018
*/

#include "problem_bigger.h"

/////////////////// Problem_bigger

Problem_bigger::Problem_bigger(int len){
  other_init_done = false;
  _length  = len;
  _u0     = new double[len];
  _uold   = new double[len];
  _unew   = new double[len];
  _ucheck = new double[len];
}

void Problem_bigger::do_step(double t, double dt, double * uout, Waveform * WF_in){
  _WF = WF_in;
  // simple explicit Euler step
  rhs(t, _uold, _unew);
  for(int i = 0; i < _length; i++){
    _unew[i] = _uold[i] + dt*_unew[i];
    uout [i] = _unew[i];
    _uold[i] = _unew[i];
  }
}

void Problem_bigger::create_checkpoint(){
  for(int i = 0; i < _length; i++)
    _ucheck[i] = _uold[i];
}

void Problem_bigger::reset_to_checkpoint(){
  for(int i = 0; i < _length; i++)
    _uold[i] = _ucheck[i];
}

/////////////////// Part 1

Problem_bigger_part_1::Problem_bigger_part_1() : Problem_bigger(3){
  _u0[0]   = 1.;
  _u0[1]   = 2.;
  _u0[2]   = 3.;
  for(int i = 0; i < _length; i++){
    _uold[i] = _u0[i];
    _unew[i] = _u0[i];
  }
}

void Problem_bigger_part_1::rhs(double t, double *uold, double *uout){
  _WF->eval(t, _uother);
  uout[0] = -uold[0] + uold[1];
  uout[1] = -uold[1] + uold[2];
  uout[2] = -uold[2] + _uother[0];
}

/////////////////// Part 2

Problem_bigger_part_2::Problem_bigger_part_2(): Problem_bigger(2){
  _u0[0]   = 4.;
  _u0[1]   = 5.;
  for(int i = 0; i < _length; i++){
    _uold[i] = _u0[i];
    _unew[i] = _u0[i];
  }
}

void Problem_bigger_part_2::rhs(double t, double *uold, double *uout){
  _WF->eval(t, _uother);
  uout[0] = -uold[0] + uold[1];
  uout[1] = -uold[1] + _uother[0];
}
