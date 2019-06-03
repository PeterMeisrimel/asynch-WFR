/*
Authors: Peter Meisrimel
December 2018
*/

#include "WFR_JAC.h"
#include "waveform.h"
#include "problem.h"
#include "math.h" // for sqrt
#include "mpi.h"

WFR_JAC::WFR_JAC(int id_self, int id_other, double t_end, Problem * p, bool errlogging){
    ID_SELF  = id_self;
    ID_OTHER = id_other;
    _t_end   = t_end;
    prob    = p;
    WF_iters = 0;
    log_errors = errlogging;
    err_log_counter = 0;
}

void WFR_JAC::do_WF_iter(double WF_TOL, int WF_MAX_ITER, int steps_per_window_self, int steps_per_window_other){
    first_iter = true;
    for(int i = 0; i < WF_MAX_ITER; i++){ // WF iter loop
        WF_iters++;
        
        integrate_window(steps_per_window_self);
	    MPI_Sendrecv(WF_self_data , (steps_per_window_self  + 1) * DIM_SELF , MPI_DOUBLE, ID_OTHER, TAG_DATA, 
                     WF_other_data, (steps_per_window_other + 1) * DIM_OTHER, MPI_DOUBLE, ID_OTHER, TAG_DATA,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (check_convergence(WF_TOL)){
            break;
        }else{
            WF_self ->get_last(WF_self_last);
            WF_other->get_last(WF_other_last);
            prob    ->reset_to_checkpoint();
        }
    } // END WF iter loop
}
