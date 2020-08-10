/*
Authors: Peter Meisrimel
December 2018
*/

#ifndef WFR_NEW_H_
#define WFR_NEW_H_

#include "problem.h"
#include "mpi.h"
#include "WFR.h"

class WFR_NEW : public WFR_parallel{
// only using put operations here
protected:
    int msg_sent;

    MPI_Win WIN_data;

    virtual void do_WF_iter(double WF_TOL, int WF_MAX_ITER, int steps_per_window_self, int steps_per_window_other);
    virtual void integrate_window(Waveform * WF_calc, Waveform * WF_src, int steps_self, Problem * p, double theta_relax);

    // for creating logs of communication patterns
    bool log_pattern;
    int * iter_per_macro; // store number of iterations per macro step
    int log_p, log_m; // idx for writing into comm_pattern and iter_per_macro
public:
    WFR_NEW(int id_in_self, int id_in_other, double tend, Problem * p);

    virtual void run(double WF_TOL, int WF_MAX_ITER, int steps_macro, int steps_self, int steps_other,
                     int conv_check, int steps_converged_required_in, bool errlogging);
};

class WFR_NEW_var_relax : public WFR_NEW{
// check optimal relaxation stepwise, assumes same number of steps for both problems
private:
    // gs_self = own process is last, rationale: this step can be done mostly independent
    // gs_other = own process was first
    double theta_relax_self_ahead, theta_relax_other_ahead;

    bool * WF_other_data_recv_flag, * relax_self_done_flag, * relax_other_done_flag, * relax_self_done_flag_jac;
    MPI_Win WIN_recv_flag;
    double * WF_recv_data;
    Waveform * WF_recv;

    void do_WF_iter(double WF_TOL, int WF_MAX_ITER, int steps_per_window_self, int steps_per_window_other);
    void integrate_window(Waveform * WF_calc, Waveform * WF_src, int steps_self, Problem * p);
    const bool TRUE_SEND = true;
public:
//    virtual void set_conv_check_WF_ptr(int conv_which, bool match_which_conv_relax);

    WFR_NEW_var_relax(int id_in_self, int id_in_other, double tend, Problem * p);
    
    virtual void set_relax_params(double t1, double t2, double t3){
        theta_relax = t1; // jacobi (of other)
        theta_relax_self_ahead = t2; // 
        theta_relax_other_ahead = t3;
    };

    void run(double WF_TOL, int WF_MAX_ITER, int steps_macro, int steps_self, int steps_other,
             int conv_check, int steps_converged_required_in, bool errlogging);
};

#endif // WFR_NEW_H_
