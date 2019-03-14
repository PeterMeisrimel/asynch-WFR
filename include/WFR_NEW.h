/*
Authors: Peter Meisrimel
December 2018
*/

#ifndef WFR_NEW_H_
#define WFR_NEW_H_

#include "problem.h"
#include "mpi.h"
#include "WFR.h"

class WFR_NEW : public WFR{
// only using put operations here
private:
    int msg_sent, msg_to_recv;
    int IDX, IDX_aux, * IDX_win;

    double *WF_other_data_win;
    Waveform *WF_win;

    MPI_Win WIN_data, WIN_idx;

    void check_new_data(); // check for and recieve new data

    void do_WF_iter        (double, int, int, int); // tol, maxiter, steps_self , steps_other
    void integrate_window  (int); // steps_self
    bool check_convergence (double); // tol

    // for creating logs of communication patterns
    bool log_pattern;
    bool * comm_pattern; // creating a log of the used communication patterns used
    int * iter_per_macro; // store number of iterations per macro step
    int log_p, log_m; // idx for writing into comm_pattern and iter_per_macro
public:
    WFR_NEW(int, int, double, Problem *, bool = false); // bool = pattern log

    void run(double, int, int, int, int = 1); // tol, maxiter, macro_steps, steps_self, conv_check

    void write_log(int, int);
};

#endif // WFR_NEW_H_
