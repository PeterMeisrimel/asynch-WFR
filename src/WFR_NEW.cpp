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
#include <iostream>

WFR_NEW::WFR_NEW(int id_in_self, int id_in_other, double tend, Problem * p, bool errlogging, bool commlogging) : WFR_parallel(id_in_self, id_in_other){
    _t_end = tend;
    prob_self    = p;
    WF_iters = 0;
    //log_pattern = commlogging;
    log_pattern = false;
    log_errors = errlogging;
    err_log_counter = 0;
}

void WFR_NEW::run(double WF_TOL, int WF_MAX_ITER, int steps_macro, int steps_self, int steps_other, int conv_check, int steps_converged_required_in, double relax_param, bool match_which_conv_relax){
    conv_which = conv_check;
    steps_converged = 0;
    steps_converged_required = steps_converged_required_in;

    w_relax = relax_param;
    if (w_relax == 1)
        RELAX = false;
    else
        RELAX = true;

    DIM_SELF = prob_self->get_length();
    // Get vectors length from other problem
    MPI_Sendrecv(&DIM_SELF , 1, MPI_INT, ID_OTHER, TAG_DATA,
                 &DIM_OTHER, 1, MPI_INT, ID_OTHER, TAG_DATA,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    prob_self->init_other(DIM_OTHER);

    u0_self  = new double[DIM_SELF];
    u0_other = new double[DIM_OTHER];
    prob_self->get_u0(u0_self);

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

    if (RELAX)
        relax_aux_vec = new double[DIM_SELF];

    // initialize other waveform
    MPI_Sendrecv(u0_self , DIM_SELF, MPI_DOUBLE, ID_OTHER, TAG_DATA,
                 u0_other, DIM_OTHER, MPI_DOUBLE, ID_OTHER, TAG_DATA,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    times_other = new double[WF_LEN_OTHER];
    double dt_other = _t_end/steps_other;
    for (int i = 0; i < WF_LEN_OTHER; i++){
        times_other[i] = i*dt_other;
	}

	MPI_Win_allocate(WF_LEN_OTHER * DIM_OTHER * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &WF_other_data, &WIN_data);
    WF_other = new Waveform_locking(WF_LEN_OTHER, DIM_OTHER, times_other, WF_other_data, &WIN_data, ID_SELF);
    WF_other->set_last(u0_other);

    // if logging patterns, allocate necessary memory
    /*
    if (log_pattern){
        log_p = 0;
        log_m = 0;
        iter_per_macro = new int[steps_macro];
        comm_pattern = new bool[WF_MAX_ITER*steps_self];
    }
    */

    init_error_log(steps_macro, WF_MAX_ITER);

    double window_length = _t_end/steps_macro;
    norm_factor = prob_self -> get_norm_factor(); // implicitly assumed to be identical for both subproblems
    set_conv_check_WF_ptr(conv_which, match_which_conv_relax);

    MPI_Barrier(MPI_COMM_WORLD);
	runtime = MPI_Wtime(); // runtime measurement start
	for(int i = 0; i < steps_macro; i++){ // Macro step loop
        prob_self -> create_checkpoint();
        WF_self ->get_last(u0_self);
        WF_self ->set(0, u0_self);

        WF_other->init_by_last();

        get_relative_tol(); // get tolerance for relative update termination check

        MPI_Barrier(MPI_COMM_WORLD); // not quite sure if needed
        do_WF_iter(WF_TOL, WF_MAX_ITER, WF_LEN_SELF - 1, WF_LEN_OTHER - 1);
        if(i != steps_macro - 1){
            WF_self  -> time_shift(window_length);
            WF_other -> time_shift(window_length);
        }
	} // endfor macrostep loop
	runtime = MPI_Wtime() - runtime; // runtime measurement end

	MPI_Win_free(&WIN_data);
}

void WFR_NEW::do_WF_iter(double WF_TOL, int WF_MAX_ITER, int steps_per_window_self, int steps_per_window_other){
	first_iter = true;

    if (RELAX){
        WF_self->init_by_last(); // need a base for relaxation
    }

    steps_converged = 0;
	for(int i = 0; i < WF_MAX_ITER; i++){ // Waveform loop
        WF_iters++;

        integrate_window(WF_self, WF_other, steps_per_window_self, prob_self);

        // all outstanding RMA operations done up to this point
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, ID_SELF, 0, WIN_data); // Shared as it is read only // excl since one is writing in own window of sorts?
        MPI_Win_sync(WIN_data); // lock neccessary in intel MPI
        MPI_Win_unlock(ID_SELF, WIN_data);

        if (check_convergence(WF_TOL)){ // convergence or maximum number of iterations reached
            steps_converged++;
            /*
            if (log_pattern){
                iter_per_macro[log_m] = i+1;
                log_m++;
            }
            */
            if (steps_converged >= steps_converged_required) // sufficiently many steps registered convergence, stop
                break;
        }else{
            steps_converged = 0;
        }
        WF_self   -> get_last(WF_self_last);
        WF_other  -> get_last(WF_other_last);
        prob_self -> reset_to_checkpoint();
    } // endfor waveform loop
}

void WFR_NEW::integrate_window(Waveform * WF_calc, Waveform * WF_src, int steps, Problem * p){
    double t, dt;
    
    if (RELAX){
	    for(int i = 0; i < steps; i++){ // timestepping loop
            MPI_Win_lock(MPI_LOCK_EXCLUSIVE, ID_SELF, 0, WIN_data); // excl since it involves writing
            MPI_Win_sync(WIN_data); // lock neccessary in intel MPI
            MPI_Win_unlock(ID_SELF, WIN_data);

            // make the choice of which waveform to choose from and use linear interpolation			
            // i = number of steps done by own process
            t = WF_self->get_time(i);
            dt = WF_self->get_time(i+1) - t;

            WF_calc -> get(i+1, relax_aux_vec); // backup for self of new timestep

            // actual timestep
            p->do_step(t, dt, (*WF_calc)[i+1], WF_src); // do timestep, write into WF_self

            WF_calc -> relax_by_single(w_relax, i+1, relax_aux_vec); // do interpol in place, relax_aux_vec being previous iterate at new timepoint

		    msg_sent = i+1;
          	// Send out new data to other process
          	MPI_Win_lock(MPI_LOCK_EXCLUSIVE, ID_OTHER, 0, WIN_data);
            // directly send relaxed value
            MPI_Put((*WF_self)[msg_sent], DIM_SELF, MPI_DOUBLE, ID_OTHER, DIM_SELF * msg_sent, DIM_SELF, MPI_DOUBLE, WIN_data);
		    MPI_Win_unlock(ID_OTHER, WIN_data);
        }
    }else{
	    for(int i = 0; i < steps; i++){ // timestepping loop
            MPI_Win_lock(MPI_LOCK_SHARED, ID_SELF, 0, WIN_data); // Shared as it is read only
            MPI_Win_sync(WIN_data); // lock neccessary in intel MPI
            MPI_Win_unlock(ID_SELF, WIN_data);

            // make the choice of which waveform to choose from and use linear interpolation			
            // i = number of steps done by own process
            t = WF_self->get_time(i);
            dt = WF_self->get_time(i+1) - t;
            
            // designed for implicit time_integration in current form
            /// \todo figure out better method to accurately show percentage of new data being involved?
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
            p->do_step(t, dt, (*WF_calc)[i+1], WF_src);

		    msg_sent = i+1;
          	// Send out new data to other process
          	MPI_Win_lock(MPI_LOCK_EXCLUSIVE, ID_OTHER, 0, WIN_data);
            MPI_Put((*WF_self)[msg_sent], DIM_SELF, MPI_DOUBLE, ID_OTHER, DIM_SELF * msg_sent, DIM_SELF, MPI_DOUBLE, WIN_data);
		    MPI_Win_unlock(ID_OTHER, WIN_data);
        }
    }
}

void WFR_NEW::write_log(int macro, int steps){
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

/////////////////////////////////
// OPT RELAX TESTING
/////////////////////////////////
void WFR_NEW_relax_opt::set_conv_check_WF_ptr(int conv_which, bool match_which_conv_relax){
    // important note: relaxation is now done at receiving end, RELAX flags will be opposite of previous case
    switch(conv_which){
        case -1:{ // use ouput of first model (ID_SELF == 0) as measurement for convergence
            if (ID_SELF == 0){
                WF_conv_check = WF_self;
                WF_conv_check_last = WF_self_last;
                if (match_which_conv_relax) // relaxation for output of first model, done in second model
                    RELAX = false; // => deactivate relaxation in first model
            }else{
                WF_conv_check = WF_other;
                WF_conv_check_last = WF_other_last;
            }
            break;   
        }
        case -2:{ // use ouput of second model (ID_SELF == 1) as measurement for convergence
            if (ID_SELF == 0){
                WF_conv_check = WF_other;
                WF_conv_check_last = WF_other_last;
            }else{
                WF_conv_check = WF_self;
                WF_conv_check_last = WF_self_last;
                if (match_which_conv_relax) // relaxation of output of second model, done in first model
                    RELAX = false; // => deactive relaxation in second model
            }
            break;
        }
    }
}

WFR_NEW_relax_opt::WFR_NEW_relax_opt(int id_in_self, int id_in_other, double tend, Problem * p, bool errlogging, bool commlogging, double w_relax_gs) : WFR_NEW(id_in_self, id_in_other, tend, p, errlogging, commlogging){
    /*
    w_relax logic, a given problem would know what relax parameter it should take if it goes first,
    but relaxation is done at receiving end. Thus a given problem/process knows the relaxation parameter for the other problem
    */
    w_relax_gs_other = w_relax_gs;
    RELAX = true; 
}

void WFR_NEW_relax_opt::run(double WF_TOL, int WF_MAX_ITER, int steps_macro, int steps_self, int steps_other, int conv_check, int steps_converged_required_in, double relax_param, bool match_which_conv_relax){
    conv_which = conv_check;
    steps_converged = 0;
    steps_converged_required = steps_converged_required_in;

    w_relax_jac = relax_param;

    DIM_SELF = prob_self->get_length();
    // Get vectors length from other problem
    MPI_Sendrecv(&DIM_SELF , 1, MPI_INT, ID_OTHER, TAG_DATA,
                 &DIM_OTHER, 1, MPI_INT, ID_OTHER, TAG_DATA,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    prob_self->init_other(DIM_OTHER);

    u0_self  = new double[DIM_SELF];
    u0_other = new double[DIM_OTHER];
    prob_self->get_u0(u0_self);

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

    // initialize other waveform
    MPI_Sendrecv(u0_self , DIM_SELF, MPI_DOUBLE, ID_OTHER, TAG_DATA,
                 u0_other, DIM_OTHER, MPI_DOUBLE, ID_OTHER, TAG_DATA,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    times_other = new double[WF_LEN_OTHER];
    double dt_other = _t_end/steps_other;
    for (int i = 0; i < WF_LEN_OTHER; i++)
        times_other[i] = i*dt_other;

    // WF_other contains relaxed data, does not need to be in window now
    WF_other_data = new double[WF_LEN_OTHER*DIM_OTHER];
    WF_other = new Waveform(WF_LEN_OTHER, DIM_OTHER, times_other, WF_other_data);
    WF_other -> set_last(u0_other);

    // Window & Waveform for receiving data
	MPI_Win_allocate(WF_LEN_OTHER * DIM_OTHER * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &WF_recv_data, &WIN_data);
    WF_recv = new Waveform_locking(WF_LEN_OTHER, DIM_OTHER, times_other, WF_recv_data, &WIN_data, ID_SELF);
    // flags for which data has already been received
    MPI_Win_allocate(WF_LEN_OTHER * sizeof(bool), sizeof(bool), MPI_INFO_NULL, MPI_COMM_WORLD, &WF_other_data_recv_flag, &WIN_recv_flag);
    // init flags
    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, ID_SELF, 0, WIN_recv_flag); // Excl, writing
    for (int i = 1; i < WF_LEN_OTHER; i++)
        WF_other_data_recv_flag[i] = false;
    MPI_Win_unlock(ID_SELF, WIN_recv_flag);

    // flags for when other process did GS relax
    relax_other_done_flag = new bool[WF_LEN_OTHER];
    // flags for marking which entry already received relaxation (GS)
    relax_self_done_flag = new bool[WF_LEN_OTHER];
    // flags for marking which entry already received relaxation (jac)
    relax_self_done_flag_jac = new bool[WF_LEN_OTHER];
    // init flags
    for (int i = 1; i < WF_LEN_OTHER; i++){
        relax_self_done_flag[i] = false;
        relax_other_done_flag[i] = false;
        relax_self_done_flag_jac[i] = false;
    }
    
    init_error_log(steps_macro, WF_MAX_ITER);

    double window_length = _t_end/steps_macro;
    norm_factor = prob_self -> get_norm_factor(); // implicitly assumed to be identical for both subproblems
    set_conv_check_WF_ptr(conv_which, match_which_conv_relax);
    
    MPI_Sendrecv(&w_relax_gs_other, 1, MPI_DOUBLE, ID_OTHER, TAG_MISC,
                 &w_relax_gs_self , 1, MPI_DOUBLE, ID_OTHER, TAG_MISC,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // deactivate relaxation by setting all parameters to 1
    if (not RELAX){
        w_relax_gs_self = 1;
        w_relax_gs_other = 1;
        w_relax_jac = 1;
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
	runtime = MPI_Wtime(); // runtime measurement start
	for(int i = 0; i < steps_macro; i++){ // Macro step loop
        prob_self -> create_checkpoint();
        WF_self ->get_last(u0_self);
        WF_self ->set(0, u0_self);

        WF_other->init_by_last();

        get_relative_tol(); // get tolerance for relative update termination check

        MPI_Barrier(MPI_COMM_WORLD); // not quite sure if needed
        do_WF_iter(WF_TOL, WF_MAX_ITER, WF_LEN_SELF - 1, WF_LEN_OTHER - 1);
        if(i != steps_macro - 1){
            WF_self  -> time_shift(window_length);
            WF_other -> time_shift(window_length);
        }
	} // endfor macrostep loop
	runtime = MPI_Wtime() - runtime; // runtime measurement end

	MPI_Win_free(&WIN_data);
}

void WFR_NEW_relax_opt::do_WF_iter(double WF_TOL, int WF_MAX_ITER, int steps_per_window_self, int steps_per_window_other){
	first_iter = true;
    steps_converged = 0;
	for(int i = 0; i < WF_MAX_ITER; i++){ // Waveform loop
        MPI_Barrier(MPI_COMM_WORLD);
        WF_iters++;

        integrate_window(WF_self, WF_other, steps_per_window_self, prob_self);

        // relax_self_done_flag are the flags for relaxation on the data of the other process on the own process
        // thus the size of relax_self_done_flag is determined by the steps of the other process
        // this seems unnecessarily confusing, possibly find better names?
        MPI_Sendrecv(relax_self_done_flag, steps_per_window_other + 1, MPI_C_BOOL, ID_OTHER, TAG_MISC,
                     relax_other_done_flag, steps_per_window_self + 1, MPI_C_BOOL, ID_OTHER, TAG_MISC,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // all outstanding RMA (send) operations done up to this point
        //MPI_Barrier(MPI_COMM_WORLD); // previous sendrecv acts as barrier?
        // sync data
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, ID_SELF, 0, WIN_data); // Excl since it is updating self
        MPI_Win_sync(WIN_data); // lock neccessary in intel MPI
        MPI_Win_unlock(ID_SELF, WIN_data);
        // necessary to also synch flags?
        // WF_recv and WF_other_data_recv_flag should now be up to date

        // do all outstanding relaxation (jacobi)
        MPI_Win_lock(MPI_LOCK_SHARED, ID_SELF, 0, WIN_data);
        for (int i = 0; i < steps_per_window_other; i++){
            if (relax_self_done_flag[i+1] or relax_self_done_flag_jac[i+1]){
                continue;
            }else{
                if (relax_other_done_flag[i+1]){
                    for (int j = 0; j < DIM_OTHER; j++)
                        WF_other_data[(i+1)*DIM_OTHER + j] = (1 - w_relax_gs_other)*WF_other_data[(i+1)*DIM_OTHER + j] + w_relax_gs_other*WF_recv_data[(i+1)*DIM_OTHER + j];
                }else{
                    for (int j = 0; j < DIM_OTHER; j++)
                        WF_other_data[(i+1)*DIM_OTHER + j] = (1 - w_relax_jac)*WF_other_data[(i+1)*DIM_OTHER + j] + w_relax_jac*WF_recv_data[(i+1)*DIM_OTHER + j];
                }
            }
        }
        MPI_Win_unlock(ID_SELF, WIN_data);

        // reset flags
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, ID_SELF, 0, WIN_recv_flag); // Excl, writing
        for (int i = 0; i < steps_per_window_other; i++)
            WF_other_data_recv_flag[i+1] = false;
        MPI_Win_unlock(ID_SELF, WIN_recv_flag);

        for (int i = 1; i < steps_per_window_other + 1; i++){
            relax_self_done_flag[i] = false;
            relax_self_done_flag_jac[i] = false;
        }

        // since own data is discarded after each iteration over a given time-window, relaxation is only done at receiver on other
        // the exception is the convergence check, for which one need the relaxed data
        MPI_Sendrecv((*WF_other)[steps_per_window_other], DIM_OTHER , MPI_DOUBLE, ID_OTHER, TAG_DATA, 
                     (*WF_self)[steps_per_window_self], DIM_SELF, MPI_DOUBLE, ID_OTHER, TAG_DATA,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (check_convergence(WF_TOL)){ // convergence or maximum number of iterations reached
            steps_converged++;
            if (steps_converged >= steps_converged_required) // sufficiently many steps registered convergence, stop
                break;
        }else{
            steps_converged = 0;
        }
        WF_self   -> get_last(WF_self_last);
        WF_other  -> get_last(WF_other_last);
        prob_self -> reset_to_checkpoint();
    } // endfor waveform loop
}

void WFR_NEW_relax_opt::integrate_window(Waveform * WF_calc, Waveform * WF_src, int steps, Problem * p){
    double t, dt;

	for(int i = 0; i < steps; i++){ // timestepping loop
        // first sync flags, other way around a flag might be true, but data might not be there yet (flag updated during data sync.)
        // alternative: enclose locks?
        // synch flags marking which data is new
        MPI_Win_lock(MPI_LOCK_SHARED, ID_SELF, 0, WIN_recv_flag); // Shared as it is read only
        MPI_Win_sync(WIN_recv_flag); // lock neccessary in intel MPI
        MPI_Win_unlock(ID_SELF, WIN_recv_flag);
        // sync data exposed by window, i.e. catch new updates
        MPI_Win_lock(MPI_LOCK_SHARED, ID_SELF, 0, WIN_data); // Shared as it is read only
        MPI_Win_sync(WIN_data); // lock neccessary in intel MPI
        MPI_Win_unlock(ID_SELF, WIN_data);

        // \todo: determine the truth values of the if-cases beforehand to reduce total locking time on flags
        MPI_Win_lock(MPI_LOCK_SHARED, ID_SELF, 0, WIN_recv_flag); // read only
        if (WF_other_data_recv_flag[i+1]){
            // other process is already done with corresponding (the one we are calculating right now) timestep, do GS relaxation, 
            // do relaxation, GS
            MPI_Win_lock(MPI_LOCK_SHARED, ID_SELF, 0, WIN_data);
            for (int j = 0; j < DIM_OTHER; j++)
                WF_other_data[(i+1)*DIM_OTHER + j] = (1 - w_relax_gs_self)*WF_other_data[(i+1)*DIM_OTHER + j] + w_relax_gs_self*WF_recv_data[(i+1)*DIM_OTHER + j];
            MPI_Win_unlock(ID_SELF, WIN_data);
            relax_self_done_flag[i+1] = true; // mark timestep as relaxated via GS

            //if ((not relax_self_done_flag[i]) and WF_other_data_recv_flag[i]){ // data not relaxed, but new values arrived in the meantime (during previous timestep)
            if (not relax_self_done_flag[i]){ // data not relaxed, but new values arrived in the meantime (during previous timestep)
                MPI_Win_lock(MPI_LOCK_SHARED, ID_SELF, 0, WIN_data);
                for (int j = 0; j < DIM_OTHER; j++)
                    WF_other_data[i*DIM_OTHER + j] = (1 - w_relax_jac)*WF_other_data[i*DIM_OTHER + j] + w_relax_jac*WF_recv_data[i*DIM_OTHER + j];
                MPI_Win_unlock(ID_SELF, WIN_data);
                relax_self_done_flag_jac[i] = true;
            }
        }
        MPI_Win_unlock(ID_SELF, WIN_recv_flag);

        // make the choice of which waveform to choose from and use linear interpolation			
        // i = number of steps done by own process
        t = WF_self->get_time(i);
        dt = WF_self->get_time(i+1) - t;

        // actual timestep
        p->do_step(t, dt, (*WF_calc)[i+1], WF_other);

        msg_sent = i+1;
        // Send out new data to other process
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, ID_OTHER, 0, WIN_data);
        MPI_Put((*WF_self)[msg_sent], DIM_SELF, MPI_DOUBLE, ID_OTHER, DIM_SELF * msg_sent, DIM_SELF, MPI_DOUBLE, WIN_data);
        MPI_Win_unlock(ID_OTHER, WIN_data);

        // Send flag for new data
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, ID_OTHER, 0, WIN_recv_flag);
        MPI_Put(&TRUE_SEND, 1, MPI_C_BOOL, ID_OTHER, msg_sent, 1, MPI_C_BOOL, WIN_recv_flag);
        MPI_Win_unlock(ID_OTHER, WIN_recv_flag);        
    }
}
