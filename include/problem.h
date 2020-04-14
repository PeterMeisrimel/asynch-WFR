/*
Authors: Peter Meisrimel
December 2018
*/

#ifndef PROBLEM_H_
#define PROBLEM_H_

#include "waveform.h"

// general form of a problem to be coupled
class Problem{
protected:
    double * _u0, * _uother; // initial value and storage for evaluating waveforms
    int _length, _length_other; // corresponding lenghtes of the problems
    bool other_init_done; // marker that init_other has been called
public:
    int get_length(){ return _length;}; // 
    // u (inout)
    virtual void get_u0(double *uout){
        for(int i = 0; i < _length; i++)
            uout[i] = _u0[i];
    };

    // initialization of _uother which has the correct size for evaluating the waveform, may be overwritten for further initialization
    virtual void init_other(int len_other){
        _length_other = len_other;
        if(not other_init_done){
            _length_other = len_other;
            _uother = new double[len_other];
            other_init_done = true;
        }
    };

    // t, dt, unew (out), WF
    virtual void do_step(double, double, double *, Waveform *) = 0;

    virtual void create_checkpoint()   = 0; // to create backup of possible internal variables at start of macrostep
    virtual void reset_to_checkpoint() = 0; // if macro step needs to be repeated, resetting to previous start
};

#endif //PROBLEM_H_