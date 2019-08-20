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
private:
    int msg_sent;

    MPI_Win WIN_data;

    void do_WF_iter(double WF_TOL, int WF_MAX_ITER, int steps_per_window_self, int steps_per_window_other);
    void integrate_window(Waveform * WF_calc, Waveform * WF_src, int steps_self, Problem * p);

    // for creating logs of communication patterns
    bool log_pattern;
    bool * comm_pattern; // creating a log of the used communication patterns used
    int * iter_per_macro; // store number of iterations per macro step
    int log_p, log_m; // idx for writing into comm_pattern and iter_per_macro
public:
    WFR_NEW(int id_in_self, int id_in_other, double tend, Problem * p, bool errlogging = false, bool commlogging = false);

    void run(double WF_TOL, int WF_MAX_ITER, int steps_macro, int steps_self, int steps_other, int conv_check, int steps_converged_required_in);

    void write_log(int macro, int steps);
};

#endif // WFR_NEW_H_
