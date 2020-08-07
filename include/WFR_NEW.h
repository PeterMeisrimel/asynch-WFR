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
    virtual void integrate_window(Waveform * WF_calc, Waveform * WF_src, int steps_self, Problem * p);

    // for creating logs of communication patterns
    bool log_pattern;
    int * iter_per_macro; // store number of iterations per macro step
    int log_p, log_m; // idx for writing into comm_pattern and iter_per_macro
public:
    WFR_NEW(int id_in_self, int id_in_other, double tend, Problem * p, bool errlogging = false);

    virtual void run(double WF_TOL, int WF_MAX_ITER, int steps_macro, int steps_self, int steps_other, int conv_check, int steps_converged_required_in, double relax_param);
};

class WFR_NEW_relax_opt : public WFR_NEW{
// check optimal relaxation stepwise, assumes same number of steps for both problems
private:
    // gs_self = own process is last, rationale: this step can be done mostly independent
    // gs_other = own process was first
    double w_relax_jac, w_relax_gs_self, w_relax_gs_other;

    bool * WF_other_data_recv_flag, * relax_self_done_flag, * relax_other_done_flag, * relax_self_done_flag_jac;
    MPI_Win WIN_recv_flag;
    double * WF_recv_data;
    Waveform * WF_recv;

    void do_WF_iter(double WF_TOL, int WF_MAX_ITER, int steps_per_window_self, int steps_per_window_other);
    void integrate_window(Waveform * WF_calc, Waveform * WF_src, int steps_self, Problem * p);
    const bool TRUE_SEND = true;
public:
//    virtual void set_conv_check_WF_ptr(int conv_which, bool match_which_conv_relax);

    WFR_NEW_relax_opt(int id_in_self, int id_in_other, double tend, Problem * p, bool errlogging = false, double w_relax_gs = 1);

    void run(double WF_TOL, int WF_MAX_ITER, int steps_macro, int steps_self, int steps_other, int conv_check, int steps_converged_required_in, double relax_param);
};

#endif // WFR_NEW_H_
