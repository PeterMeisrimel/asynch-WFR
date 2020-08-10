/*
Authors: Peter Meisrimel
December 2018
*/

#include "WFR_GS.h"
#include "waveform.h"
#include "problem.h"
#include "math.h" // for sqrt
#include "mpi.h"
#include <stdexcept>
#include <iostream>

WFR_GS::WFR_GS(double t_end, Problem * p1, Problem * p2, bool first) : WFR_serial(){
    _t_end   = t_end;
    prob_self  = p1;
    prob_other = p2;
    FIRST    = first;
    WF_iters = 0;
}

void WFR_GS::run(double WF_TOL, int WF_MAX_ITER, int steps_macro, int steps_self, int steps_other, int conv_check, int nsteps_conv_check, bool errlogging, double relax_param){
    steps_converged = 0;
    steps_converged_required = nsteps_conv_check;
    conv_which = conv_check;

    w_relax = relax_param;
    if (w_relax == 1)
        RELAX = false;
    else
        RELAX = true;

    DIM_SELF = prob_self->get_length();
    DIM_OTHER = prob_other->get_length();
    prob_self->init_other(DIM_OTHER);
    prob_other->init_other(DIM_SELF);

    u0_self  = new double[DIM_SELF];
    u0_other = new double[DIM_OTHER];
    prob_self->get_u0(u0_self);
    prob_other->get_u0(u0_other);

    WF_self_last  = new double[DIM_SELF];
    WF_other_last = new double[DIM_OTHER];

    int WF_LEN_SELF  = steps_self/steps_macro + 1;
    int WF_LEN_OTHER = steps_other/steps_macro + 1;

    // initiliaze own waveform
    times_self = new double[WF_LEN_SELF];
    double dt_self = _t_end/steps_self;
    for (int i = 0; i < WF_LEN_SELF; i++)
        times_self[i] = i*dt_self;

    WF_self_data = new double[WF_LEN_SELF * DIM_SELF];
    WF_self      = new Waveform(WF_LEN_SELF, DIM_SELF, times_self, WF_self_data);
    WF_self->set_last(u0_self);

    relax_aux_vec = new double[DIM_SELF];
    
    // initialize other waveform
    times_other = new double[WF_LEN_OTHER];
    double dt_other = _t_end/steps_other;
    for (int i = 0; i < WF_LEN_OTHER; i++)
        times_other[i] = i*dt_other;

    WF_other_data = new double[WF_LEN_OTHER * DIM_OTHER];
    WF_other      = new Waveform(WF_LEN_OTHER, DIM_OTHER, times_other, WF_other_data);
    WF_other->set_last(u0_other);

    log_errors = errlogging;
    err_log_counter = 0;
    init_error_log(steps_macro, WF_MAX_ITER);

    double window_length = _t_end/steps_macro;

    norm_factor = prob_self -> get_norm_factor(); // implicitly assumed to be identical for both subproblems

    MPI_Barrier(MPI_COMM_WORLD);
    runtime = MPI_Wtime(); // runtime measurement start
    for(int i = 0; i < steps_macro; i++){
        prob_self->create_checkpoint();
        prob_other->create_checkpoint();

        WF_self ->get_last(u0_self);
        WF_self ->set(0, u0_self);

        WF_other->get_last(u0_other);
        WF_other->set(0, u0_other);

        if (FIRST){ // first == true, self first, then other
            WF_other->init_by_last();
        }else{
            WF_self->init_by_last();
        }

        get_relative_tol(); // get tolerance for relative update termination check
      
        do_WF_iter(WF_TOL, WF_MAX_ITER, WF_LEN_SELF - 1, WF_LEN_OTHER - 1);
        
        if(i != steps_macro - 1){
            WF_self ->time_shift(window_length);
            WF_other->time_shift(window_length);
        }
    }
    runtime = MPI_Wtime() - runtime; // runtime measurement end
}

void WFR_GS::do_WF_iter(double WF_TOL, int WF_MAX_ITER, int steps_per_window_self, int steps_per_window_other){
    first_iter = true;
    if (RELAX){
        WF_other->init_by_last();
        WF_self->init_by_last();
    }
    for(int i = 0; i < WF_MAX_ITER; i++){ // WF iter loop
        WF_iters++;
        
        if (FIRST){
            RELAX = RELAX_0;
            integrate_window(WF_self , WF_other, steps_per_window_self , prob_self);
            RELAX = RELAX_1;
            integrate_window(WF_other, WF_self , steps_per_window_other, prob_other);
        }else{
            RELAX = RELAX_1;
            integrate_window(WF_other, WF_self , steps_per_window_other, prob_other);
            RELAX = RELAX_0;
            integrate_window(WF_self , WF_other, steps_per_window_self , prob_self);
        }
        
        if (check_convergence(WF_TOL)){
            steps_converged++;
            if (steps_converged >= steps_converged_required) // sufficiently many steps registered convergence, stop
                break;
        }else{
            steps_converged = 0;
        }
        WF_self    ->get_last(WF_self_last);
        WF_other   ->get_last(WF_other_last);
        prob_self  ->reset_to_checkpoint();
        prob_other ->reset_to_checkpoint();
    } // END WF iter loop
}
