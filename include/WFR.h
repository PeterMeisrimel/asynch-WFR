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
const int TAG_MISC = 5;

#endif // MPI_TAGS_

#ifndef WFR_H_
#define WFR_H_

#include "problem.h"
#include <iostream>
#include "mpi.h"

/*
General viewing point is to consider a problem by its in -and outputs (waveforms)
as such SELF = output of own, OTHER = input

for classical, constant splititng methods, relaxation is done on the sending end
for NEW method, with variable splittings, relaxation is done on the receiving end

termination check is done on both processes either way, use an output based notation here
i.e., using only a single output for the termination check, we mark "2" as checking the output of second problem

do relaxation, even if relaxation parameter might be one, for consistency?
enable use of 2 different relaxation parameters, as to enable relaxation on only a single in/output
=> pass single relaxation parameter into each run function, each for their own ouput
=> new method, needs 3 input parameters
exchange them via mpi calls

=> need to eliminate "match relax to conv" parameter and replace with something more general and sensible
*/

class WFR{
protected:
    double _t_end;
    Problem * prob_self;

    int np;

    int       ID_SELF     ,  ID_OTHER;
    int       DIM_SELF    ,  DIM_OTHER;
    double   *u0_self     , *u0_other;
    Waveform *WF_self     , *WF_other;
    double   *WF_self_data, *WF_other_data;
    double   *WF_self_last, *WF_other_last;
    double   *times_self  , *times_other;

    // relaxation
//    double theta_self, theta_other; // relaxation parameters
    bool RELAX;
    double w_relax;
    double * relax_aux_vec;
    
    double theta_relax;

    // used in termination criterion
    bool first_iter;
    // which norm to use for checking for convergence
    // 0: normal 2 norm, 1 (default): weighted scaled norm
    // negative for single check, defined via outputs
    // -1: 2-norm of output of first system, self on p0, other on p1
    // -2: 2-norm of output of second system, other on p0, self on p1
    int conv_which, steps_converged, steps_converged_required;

    // for convergence checking
    double rel_update_fac, update, up_self, up_other, norm_factor;

    // gathering and printing of results
    int WF_iters, sol_size, sol_other_size, iters;
    double * sol, *sol_other;
    double runtime, runtime_self, runtime_other;

    // logging of error over iterations
    bool log_errors;
    int err_log_counter;
    double * error_log, * error_other_log;

    virtual void do_WF_iter(double WF_TOL, int WF_MAX_ITER, int steps_self, int steps_other) = 0;
    virtual void integrate_window(Waveform * WF_calc, Waveform * WF_src, int steps_self, Problem * p);
    virtual void integrate_window(int start, Waveform * WF_calc, Waveform * WF_src, int steps_self, Problem * p);
public:
    WFR(){
        MPI_Comm_size(MPI_COMM_WORLD, &np);
    };
    virtual void run(double WF_TOL, int WF_MAX_ITER, int steps_macro, int steps_self, int steps_other,
                     int conv_check = 1, int steps_converged_in = 1, bool errorlogging = false) = 0;
    virtual void set_relax_params(double t){theta_relax = t;};
    virtual void set_relax_params(double, double){};
    virtual void set_relax_params(double, double, double){};

    virtual bool check_convergence(double WF_TOL);
    virtual void get_relative_tol();

    int get_WF_iters(){ return WF_iters;};
    double get_runtime(){ return runtime;};

    virtual void write_results();

    virtual void init_error_log(int steps_macro, int WF_MAX_ITER) = 0;
    virtual void update_error_log() = 0;
    virtual void write_error_log();
};

class WFR_serial: public WFR{
protected:
    Problem * prob_other;
    bool FIRST; // to define the ordering for coupling two problems, true = 1 -> 2 ordering
    bool RELAX_0, RELAX_1; // need more relaxation flags for possibility of only relaxing a single interface
    
    double theta_relax_self, theta_relax_other;
public:
    WFR_serial() : WFR(){
        ID_SELF = 0;
        ID_OTHER = 1;
    };
    virtual void set_relax_params(double t1, double t2){
        theta_relax_self = t1;
        theta_relax_other = t2;
    };
    
    void write_results();

    void init_error_log(int, int);
    void update_error_log();
};

class WFR_parallel: public WFR{
public:
    WFR_parallel(int id_in_self, int id_in_other) : WFR(){
        ID_SELF = id_in_self;
        ID_OTHER = id_in_other;
    };
    void write_results();

    void init_error_log(int, int);
    void update_error_log();
    void write_error_log();
};

#endif // WFR_GS_H_
