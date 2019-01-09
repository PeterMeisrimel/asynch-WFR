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
private:
    int msg_sent, msg_to_recv;
    int IDX, IDX_aux, * IDX_win;

    double *WF_other_data_new, *WF_other_data_win;
    Waveform *WF_other_new, *WF_win, *WF_curr;

    MPI_Win WIN_data, WIN_idx;

    void check_new_data(); // check for and recieve new data

    void do_WF_iter        (double, int, int, int); // tol, maxiter, steps_self , steps_other
    void integrate_window  (int); // steps_self
    bool check_convergence (double); // tol
public:
    WFR_NEW(int, int, double, Problem *);

    void run               (double, int, int, int, int = 1); // tol, maxiter, macro_steps, steps_self, conv_check
};

#endif // WFR_NEW_H_
