/*
Authors: Peter Meisrimel
December 2018
*/

#ifndef MPI_TAGS_
#define MPI_TAGS_

const int TAG_IDX  = 1;
const int TAG_TIME = 2;
const int TAG_DATA = 3;
const int TAG_DONE = 4;

#endif // MPI_TAGS_

#ifndef WFR_H_
#define WFR_H_

#include "problem.h"
#include <iostream>

class WFR{
protected:
    bool first_iter;
    // which norm to use for checking for convergence
    // 0: normal 2 norm, 1 (default): weighted scaled norm
    // -1: 2-norm of first system, -2: 2-norm of second system
    int conv_which;
    double _t_end;
    Problem * prob;

    int       ID_SELF     ,  ID_OTHER;
    int       DIM_SELF    ,  DIM_OTHER;
    double   *u0_self     , *u0_other;
    Waveform *WF_self     , *WF_other;  
    double   *WF_self_data, *WF_other_data;
    double   *WF_self_last, *WF_other_last;
    double   *times_self  , *times_other;

    double up_self, up_other, update;
    int WF_iters;

    double runtime;
    bool log_errors;
    int err_log_counter;
    double * error_log;

    double rel_update_fac;
    int nconv, nconvmax; // convergence check requires nconv updates below tolerance to register convergence

    virtual void do_WF_iter        (double, int, int, int) = 0; // tol, maxiter, steps_self, steps_other
    virtual void integrate_window  (int)                   = 0; // steps_self
public:
    virtual void run               (double, int, int, int, int = 1) = 0; // tol, maxiter, macro_steps, steps_self, conv_which
    virtual bool check_convergence (double); // checking convergence, input: tolerance, ouput: true/false
    virtual void get_relative_tol();
    virtual void get_sol(double * out){ return WF_self->get_last(out);}; // obtain solution
    int get_WF_iters(){ return WF_iters;}; // obtain total number of waveform iteration
    double get_runtime(){ return runtime;}; // obtain runtime, not including initialization

    void write_results();
    void init_error_log(int, int);
    void update_error_log();
    void write_error_log();
    virtual void write_log(int a, int b){
        if (ID_SELF == 0)
            std::cout << "LOG NONE" << std::endl << "LOG NONE" << std::endl;
    };
};

#endif // WFR_GS_H_
