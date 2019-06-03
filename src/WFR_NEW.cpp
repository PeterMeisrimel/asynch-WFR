/*
Authors: Peter Meisrimel
December 2018
*/

#include "WFR_NEW.h"
#include "waveform.h"
#include "waveform_locking.h"
#include "problem.h"
#include <stdexcept>
#include "math.h"
#include "mpi.h"
#include "unistd.h"
#include <cassert>
#include <iostream>

WFR_NEW::WFR_NEW(int id_self, int id_other, double tend, Problem * p, bool errlogging, bool commlogging){
    ID_SELF  = id_self;
    ID_OTHER = id_other;
    _t_end = tend;
    prob    = p;
    WF_iters = 0;
    //log_pattern = commlogging;
    log_pattern = false;
    log_errors = errlogging;
    err_log_counter = 0;
}

void WFR_NEW::run(double WF_TOL, int WF_MAX_ITER, int num_macro, int steps_self, int conv_check){
    conv_which = conv_check;
    /***************************
    initialize stuff
    ****************************/
    // initialize own waveform

    DIM_SELF = prob->get_length();
    // Get vectors length from other problem
    MPI_Sendrecv(&DIM_SELF , 1, MPI_INT, ID_OTHER, TAG_DATA,
                 &DIM_OTHER, 1, MPI_INT, ID_OTHER, TAG_DATA,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    prob->init_other(DIM_OTHER);

    u0_self  = new double[DIM_SELF];
    u0_other = new double[DIM_OTHER];
    prob->get_u0(u0_self);

    WF_self_last  = new double[DIM_SELF];
    WF_other_last = new double[DIM_OTHER];

    int steps_other;
    MPI_Sendrecv(&steps_self , 1, MPI_INT, ID_OTHER, TAG_DATA,
                 &steps_other, 1, MPI_INT, ID_OTHER, TAG_DATA,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    int WF_LEN_SELF  = steps_self/num_macro + 1;
    int WF_LEN_OTHER = steps_other/num_macro + 1;

    // initiliaze own waveform
    times_self = new double[WF_LEN_SELF];
    double dt_self = _t_end/steps_self;
    for (int i = 0; i < WF_LEN_SELF; i++)
        times_self[i] = i*dt_self;

    WF_self_data = new double[WF_LEN_SELF * DIM_SELF];
    WF_self      = new Waveform(WF_LEN_SELF, DIM_SELF, times_self, WF_self_data);
    WF_self->set_last(u0_self);

    // initialize other waveform
    MPI_Sendrecv(u0_self , DIM_SELF, MPI_DOUBLE, ID_OTHER, TAG_DATA,
                 u0_other, DIM_OTHER, MPI_DOUBLE, ID_OTHER, TAG_DATA,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    times_other = new double[WF_LEN_OTHER];
    double dt_other = _t_end/steps_other;
    for (int i = 0; i < WF_LEN_OTHER; i++)
        times_other[i] = i*dt_other;

	MPI_Win_allocate(WF_LEN_OTHER * DIM_OTHER * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &WF_other_data, &WIN_data);
    WF_other = new Waveform_locking(WF_LEN_OTHER, DIM_OTHER, times_other, WF_other_data, &WIN_data, ID_SELF);
    WF_other->set_last(u0_other);

    // if logging patterns, allocate necessary memory
    /*
    if (log_pattern){
        log_p = 0;
        log_m = 0;
        iter_per_macro = new int[num_macro];
        comm_pattern = new bool[WF_MAX_ITER*steps_self];
    }
    */

    init_error_log(num_macro, WF_MAX_ITER);

    double window_length = _t_end/num_macro;

    MPI_Barrier(MPI_COMM_WORLD);
	runtime = MPI_Wtime(); // runtime measurement start
	for(int i = 0; i < num_macro; i++){ // Macro step loop
        prob -> create_checkpoint();
        WF_self ->get_last(u0_self);
        WF_self ->set(0, u0_self);

        WF_other->init_by_last();

        get_relative_tol(); // get tolerance for relative update termination check

        MPI_Barrier(MPI_COMM_WORLD); // not quite sure if needed
        do_WF_iter(WF_TOL, WF_MAX_ITER, WF_LEN_SELF - 1, WF_LEN_OTHER - 1);
        if(i != num_macro - 1){
            WF_self  -> time_shift(window_length);
            WF_other -> time_shift(window_length);
        }
	} // endfor macrostep loop
	runtime = MPI_Wtime() - runtime; // runtime measurement end

	MPI_Win_free(&WIN_data);
}

void WFR_NEW::do_WF_iter(double WF_TOL, int WF_MAX_ITER, int steps_per_window, int dummy){
	first_iter = true;
	for(int i = 0; i < WF_MAX_ITER; i++){ // Waveform loop
        WF_iters++;

        integrate_window(steps_per_window);

        // all outstanding RMA operations done up to this point
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Win_lock(MPI_LOCK_SHARED, ID_SELF, 0, WIN_data); // Shared as it is read only
        MPI_Win_sync(WIN_data); // lock neccessary in intel MPI
        MPI_Win_unlock(ID_SELF, WIN_data);

        if (check_convergence(WF_TOL) or (i == (WF_MAX_ITER - 1))){ // convergence or maximum number of iterations reached
            /*
            if (log_pattern){
                iter_per_macro[log_m] = i+1;
                log_m++;
            }
            */
            break;
        }else{
            WF_self -> get_last(WF_self_last);
            WF_other-> get_last(WF_other_last);
            prob    -> reset_to_checkpoint();
        }
    } // endfor waveform loop
}

void WFR_NEW::integrate_window(int steps){
    double t, dt;
	for(int i = 0; i < steps; i++){ // timestepping loop
        MPI_Win_lock(MPI_LOCK_SHARED, ID_SELF, 0, WIN_data); // Shared as it is read only
        MPI_Win_sync(WIN_data); // lock neccessary in intel MPI
        MPI_Win_unlock(ID_SELF, WIN_data);

        // make the choice of which waveform to choose from and use linear interpolation			
        // i = number of steps done by own process
        t = WF_self->get_time(i);
        dt = WF_self->get_time(i+1) - t;
        
        // designed for implicit time_integration in current form
        // TODO: figure out better method to accurately show percentage of new data being involved?
        /*
        if (log_pattern){
            if (times_other[IDX] >= t + dt)
                comm_pattern[log_p] = 1;
            else
                comm_pattern[log_p] = 0;
            log_p++;
        }
        */
        
        // actual timestep
        prob->do_step(t, dt, (*WF_self)[i+1], WF_other);

		msg_sent = i+1;
      	// Send out new data to other process
      	MPI_Win_lock(MPI_LOCK_EXCLUSIVE, ID_OTHER, 0, WIN_data);
        MPI_Put((*WF_self)[msg_sent], DIM_SELF, MPI_DOUBLE, ID_OTHER, DIM_SELF * msg_sent, DIM_SELF, MPI_DOUBLE, WIN_data);
		MPI_Win_unlock(ID_OTHER, WIN_data);
    }
}

void WFR_NEW::write_log(int macro, int steps){
    assert(log_pattern == true);

    int log_size = log_p; 
    int log_size_other = 0;
    
    MPI_Sendrecv(&log_size, 1, MPI_INT, ID_OTHER, 0,
                 &log_size_other, 1, MPI_INT, ID_OTHER, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    int timesteps = steps;
    int timesteps_other = 0;

    MPI_Sendrecv(&timesteps, 1, MPI_INT, ID_OTHER, 1,
                 &timesteps_other, 1, MPI_INT, ID_OTHER, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if(ID_SELF == 0){
        std::cout << "LOG " << ID_SELF << " " << timesteps << " " << macro << " ";
        for(int i = 0; i < macro; i++)
            std::cout << iter_per_macro[i] << " ";
        for(int i = 0; i < log_size; i++)
            std::cout << comm_pattern[i] << " ";
        std::cout << std::endl;

        MPI_Recv(comm_pattern, log_size_other, MPI_C_BOOL, ID_OTHER, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        std::cout << "LOG " << ID_OTHER << " " << timesteps_other << " " << macro << " ";
        for(int i = 0; i < macro; i++)
            std::cout << iter_per_macro[i] << " ";
        for(int i = 0; i < log_size_other; i++)
            std::cout << comm_pattern[i] << " ";
        std::cout << std::endl;
    }else{ // ID_SELF == 1
        MPI_Send(comm_pattern, log_size, MPI_C_BOOL, ID_OTHER, 2, MPI_COMM_WORLD);
    }
}
