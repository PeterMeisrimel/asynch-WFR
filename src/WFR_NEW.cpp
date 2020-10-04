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
#include "utils.h"

WFR_NEW::WFR_NEW(int id_in_self, int id_in_other, double tend, Problem * p) : WFR_parallel(id_in_self, id_in_other){
    _t_end = tend;
    prob_self = p;
    WF_iters = 0;
}

void WFR_NEW::run(double WF_TOL, int WF_MAX_ITER, int steps_macro, int steps_self, int steps_other,
                  int conv_check, int nsteps_conv_check, bool errlogging){
    conv_which = conv_check;
    steps_converged = 0;
    steps_converged_required = nsteps_conv_check;

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

//    if (RELAX)
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

    log_errors = errlogging;
    init_error_log(steps_macro, WF_MAX_ITER);

    double window_length = _t_end/steps_macro;
    norm_factor = prob_self -> get_norm_factor(); // implicitly assumed to be identical for both subproblems

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

    WF_self->init_by_last(); // need a base for relaxation

    steps_converged = 0;
	for(int i = 0; i < WF_MAX_ITER; i++){ // Waveform loop
        WF_iters++;

        integrate_window(WF_self, WF_other, steps_per_window_self, prob_self, theta_relax);

        // all outstanding RMA operations done up to this point
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, ID_SELF, 0, WIN_data); // Shared as it is read only // excl since one is writing in own window of sorts?
        MPI_Win_sync(WIN_data); // lock neccessary in intel MPI
        MPI_Win_unlock(ID_SELF, WIN_data);

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

void WFR_NEW::integrate_window(Waveform * WF_calc, Waveform * WF_src, int steps, Problem * p, double theta_relax){
    double t, dt;
    bool RELAX = (theta_relax != 1);
    
    for(int i = 0; i < steps; i++){ // timestepping loop
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, ID_SELF, 0, WIN_data); // excl since it involves writing
        MPI_Win_sync(WIN_data); // lock neccessary in intel MPI
        MPI_Win_unlock(ID_SELF, WIN_data);

        // make the choice of which waveform to choose from and use linear interpolation			
        // i = number of steps done by own process
        t = WF_self->get_time(i);
        dt = WF_self->get_time(i+1) - t;

        if (RELAX)
            WF_calc -> get(i+1, relax_aux_vec); // backup for self of new timestep

        // actual timestep
        p->do_step(t, dt, (*WF_calc)[i+1], WF_src); // do timestep, write into WF_self

        if (RELAX)
            WF_calc -> relax_by_single(theta_relax, i+1, relax_aux_vec); // do interpol in place, relax_aux_vec being previous iterate at new timepoint

        msg_sent = i+1;
        // Send out new data to other process
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, ID_OTHER, 0, WIN_data);
        // directly send relaxed value
        MPI_Put((*WF_self)[msg_sent], DIM_SELF, MPI_DOUBLE, ID_OTHER, DIM_SELF * msg_sent, DIM_SELF, MPI_DOUBLE, WIN_data);
        MPI_Win_unlock(ID_OTHER, WIN_data);
    }
}

/////////////////////////////////
// OPT RELAX TESTING
/////////////////////////////////

WFR_NEW_var_relax::WFR_NEW_var_relax(int id_in_self, int id_in_other, double tend, Problem * p) : WFR_NEW(id_in_self, id_in_other, tend, p){}

void WFR_NEW_var_relax::run(double WF_TOL, int WF_MAX_ITER, int steps_macro, int steps_self, int steps_other, 
                            int conv_check, int nsteps_conv_check, bool errlogging){
    conv_which = conv_check;
    steps_converged = 0;
    steps_converged_required = nsteps_conv_check;
    
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
    // \todo TODO why would this be a locking waveform?
//    WF_recv = new Waveform_locking(WF_LEN_OTHER, DIM_OTHER, times_other, WF_recv_data, &WIN_data, ID_SELF);
    WF_recv = new Waveform(WF_LEN_OTHER, DIM_OTHER, times_other, WF_recv_data);
    // flags for which data has already been received
    MPI_Win_allocate(WF_LEN_OTHER * sizeof(bool), sizeof(bool), MPI_INFO_NULL, MPI_COMM_WORLD, &WF_other_data_recv_flag, &WIN_recv_flag);
    // init flags
    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, ID_SELF, 0, WIN_recv_flag); // Excl, writing
    for (int i = 1; i < WF_LEN_OTHER; i++)
        WF_other_data_recv_flag[i] = false;
    MPI_Win_unlock(ID_SELF, WIN_recv_flag);

    // flags for when other process did GS relax
    relax_other_done_flag = new bool[WF_LEN_OTHER];
    for (int i = 1; i < WF_LEN_OTHER; i++)
        relax_other_done_flag[i] = false;
    // flags for marking which entry already received relaxation (GS)
    relax_self_done_flag = new bool[WF_LEN_OTHER];
    for (int i = 1; i < WF_LEN_OTHER; i++)
        relax_self_done_flag[i] = false;
    // flags for marking which entry already received relaxation (jac)
    relax_self_done_flag_jac = new bool[WF_LEN_OTHER];
    for (int i = 1; i < WF_LEN_OTHER; i++)
        relax_self_done_flag_jac[i] = false;
    
    log_errors = errlogging;
    init_error_log(steps_macro, WF_MAX_ITER);

    double window_length = _t_end/steps_macro;
    norm_factor = prob_self -> get_norm_factor(); // implicitly assumed to be identical for both subproblems
    
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

void WFR_NEW_var_relax::do_WF_iter(double WF_TOL, int WF_MAX_ITER, int steps_per_window_self, int steps_per_window_other){
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
                    for (int j = 0; j < DIM_OTHER; j++) // self ahead
                        WF_other_data[(i+1)*DIM_OTHER + j] = (1 - theta_relax_self_ahead)*WF_other_data[(i+1)*DIM_OTHER + j] + theta_relax_self_ahead*WF_recv_data[(i+1)*DIM_OTHER + j];
                }else{
                    for (int j = 0; j < DIM_OTHER; j++)
                        WF_other_data[(i+1)*DIM_OTHER + j] = (1 - theta_relax)*WF_other_data[(i+1)*DIM_OTHER + j] + theta_relax*WF_recv_data[(i+1)*DIM_OTHER + j];
                }
            }
        }
        MPI_Win_unlock(ID_SELF, WIN_data);

        // reset flags
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, ID_SELF, 0, WIN_recv_flag); // Excl, writing
        for (int i = 0; i < steps_per_window_other; i++)
            WF_other_data_recv_flag[i+1] = false;
        MPI_Win_unlock(ID_SELF, WIN_recv_flag);

        for (int i = 1; i < steps_per_window_other + 1; i++)
            relax_self_done_flag[i] = false;
        for (int i = 1; i < steps_per_window_other + 1; i++)
            relax_self_done_flag_jac[i] = false;

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

void WFR_NEW_var_relax::integrate_window(Waveform * WF_calc, Waveform * WF_src, int steps, Problem * p){
    double t, dt;

	for(int i = 0; i < steps; i++){ // timestepping loop
        // first sync flags, other way around a flag might be true, but data might not be there yet (flag updated during data sync.)
        // alternative: enclose locks?
        // synch flags marking which data is new
        MPI_Win_lock(MPI_LOCK_SHARED, ID_SELF, 0, WIN_recv_flag);  
        MPI_Win_sync(WIN_recv_flag); // lock neccessary in intel MPI
        MPI_Win_unlock(ID_SELF, WIN_recv_flag);
        // sync data exposed by window, i.e. catch new updates
        MPI_Win_lock(MPI_LOCK_SHARED, ID_SELF, 0, WIN_data);
        MPI_Win_sync(WIN_data); // lock neccessary in intel MPI
        MPI_Win_unlock(ID_SELF, WIN_data);

        // \todo: determine the truth values of the if-cases beforehand to reduce total locking time on flags
        MPI_Win_lock(MPI_LOCK_SHARED, ID_SELF, 0, WIN_recv_flag); // read only
        if (WF_other_data_recv_flag[i+1]){
            // other process is already done with corresponding (the one we are calculating right now) timestep, do GS relaxation, 
            // do relaxation, GS
            MPI_Win_lock(MPI_LOCK_SHARED, ID_SELF, 0, WIN_data);
            for (int j = 0; j < DIM_OTHER; j++) // other process is ahead
                WF_other_data[(i+1)*DIM_OTHER + j] = (1 - theta_relax_other_ahead)*WF_other_data[(i+1)*DIM_OTHER + j] + theta_relax_other_ahead*WF_recv_data[(i+1)*DIM_OTHER + j];
            MPI_Win_unlock(ID_SELF, WIN_data);
            relax_self_done_flag[i+1] = true; // mark timestep as relaxated via GS

            //if ((not relax_self_done_flag[i]) and WF_other_data_recv_flag[i]){ // data not relaxed, but new values arrived in the meantime (during previous timestep)
            if (not relax_self_done_flag[i]){ // data not relaxed, but new values arrived in the meantime (during previous timestep)
                MPI_Win_lock(MPI_LOCK_SHARED, ID_SELF, 0, WIN_data);
                for (int j = 0; j < DIM_OTHER; j++)
                    WF_other_data[i*DIM_OTHER + j] = (1 - theta_relax)*WF_other_data[i*DIM_OTHER + j] + theta_relax*WF_recv_data[i*DIM_OTHER + j];
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


/////////////////////////////////
// OPT RELAX TESTING + MULTIRATE
/////////////////////////////////

WFR_NEW_var_relax_MR::WFR_NEW_var_relax_MR(int id_in_self, int id_in_other, double tend, Problem * p) : WFR_NEW(id_in_self, id_in_other, tend, p){}

void WFR_NEW_var_relax_MR::run(double WF_TOL, int WF_MAX_ITER, int steps_macro, int steps_self, int steps_other, 
                            int conv_check, int nsteps_conv_check, bool errlogging){
    conv_which = conv_check;
    steps_converged = 0;
    steps_converged_required = nsteps_conv_check;
    
//    size_shared_grid = lcm(steps_self, steps_other) + 1;
//    int steps_gcd = gcd(steps_self, steps_other);
//    size_shared_grid = steps_self + steps_other - gcd(steps_self, steps_other) + 1;
    
//    std::cout << "shared grid size " << size_shared_grid << std::endl;
    
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
    
    WF_LEN_SELF  = steps_self/steps_macro + 1;
    WF_LEN_OTHER = steps_other/steps_macro + 1;
    
    size_shared_grid = (WF_LEN_SELF - 1) + (WF_LEN_OTHER - 1) - gcd(WF_LEN_SELF - 1, WF_LEN_OTHER - 1) + 1;
    
//    std::cout << ID_SELF << " shared grid size " << size_shared_grid << std::endl;
//    std::cout << ID_SELF << " steps " << steps_self << " " << steps_other << std::endl;
    
    // initiliaze own waveform
    times_self = new double[WF_LEN_SELF];
    double dt_self = _t_end/steps_self;
    for (int i = 0; i < WF_LEN_SELF; i++)
        times_self[i] = i*dt_self;
    
    WF_self_data = new double[WF_LEN_SELF * DIM_SELF];
    WF_self      = new Waveform(WF_LEN_SELF, DIM_SELF, times_self, WF_self_data);
    WF_self->set_last(u0_self);
    
    // initialize other waveform
    // contain relaxed data, thus needs to be of length size_shared_grid
    MPI_Sendrecv(u0_self , DIM_SELF, MPI_DOUBLE, ID_OTHER, TAG_DATA,
                 u0_other, DIM_OTHER, MPI_DOUBLE, ID_OTHER, TAG_DATA,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                 
    // times for receiving data window; required in construction of shared time grid
    times_other = new double[WF_LEN_OTHER];
    double dt_other = _t_end/steps_other;
    for (int i = 0; i < WF_LEN_OTHER; i++)
        times_other[i] = i*dt_other;
                 
    // time grid for shared grid
    times_shared = new double[size_shared_grid];
    times_shared[0] = 0;
    times_shared[size_shared_grid - 1] = _t_end;
    
    // flags for marking which grid a point belongs to
    // checked via & FLAG_GRID_SELF or FLAG_GRID_OTHER
    times_shared_flags = new int[size_shared_grid];
    times_shared_flags[0] = FLAG_GRID_SELF + FLAG_GRID_OTHER; // first
    times_shared_flags[size_shared_grid - 1] = FLAG_GRID_SELF + FLAG_GRID_OTHER; // last
    
    // mapping from receving grid to shared grid
    mapping_recv_to_shared = new int[WF_LEN_OTHER];
    mapping_recv_to_shared[0] = 0;
    mapping_recv_to_shared[WF_LEN_OTHER - 1] = size_shared_grid - 1;
    
    mapping_self_to_shared = new int[WF_LEN_SELF];
    mapping_self_to_shared[0] = 0;
    mapping_self_to_shared[WF_LEN_SELF - 1] = size_shared_grid - 1;
    
    int i_self = 1;
    int i_other = 1;
    // iterate through both grids, pick minimum time point for shared grid
    for (int i = 1; i < size_shared_grid - 1; i++){
//        if (ID_SELF == 0)
//            std::cout << i << " " << times_self[i_self] << " " << times_other[i_other] << std::endl;
        // shared time point
        if (std::abs(times_self[i_self] - times_other[i_other]) < 1e-10){ // time-point on both grids
            times_shared[i] = times_self[i_self];
            times_shared_flags[i] = FLAG_GRID_SELF + FLAG_GRID_OTHER;
            
            mapping_recv_to_shared[i_other] = i;
            i_other++;
            
            mapping_self_to_shared[i_self] = i;
            i_self++;
        }else{
            // smaller one on times_self
            if (times_self[i_self] < times_other[i_other]){
                times_shared[i] = times_self[i_self];
                times_shared_flags[i] = FLAG_GRID_SELF;
                mapping_self_to_shared[i_self] = i;
                i_self++;
            }else{ // times_self[i_self] > times_other[i_other]
                times_shared[i] = times_other[i_other];
                times_shared_flags[i] = FLAG_GRID_OTHER;
                mapping_recv_to_shared[i_other] = i;
                i_other++;
            }
        }
    }
    
    // auxiliary vector for interpolation on receiving waveform before relaxation
    relax_interpol_aux = new double[DIM_OTHER];
    
//    
//    
//    std::cout << "shared grid" << std::endl;
//    for (int i = 0; i < size_shared_grid; i++)
//        std::cout << times_shared[i] << " ";
//    std::cout << std::endl;
//    
//    std::cout << "shared grid SELF " << ID_SELF << " " << FLAG_GRID_SELF << " " << FLAG_GRID_OTHER<< std::endl;
//    for (int i = 0; i < size_shared_grid; i++)
//        std::cout << times_shared_flags[i] << " ";
////        std::cout << (times_shared_flags[i]&FLAG_GRID_SELF) << " ";
//    std::cout << std::endl;
//    
//    std::cout << ID_SELF << " mapping" << std::endl;
//    for (int i = 0; i < WF_LEN_OTHER; i++)
//        std::cout << mapping_recv_to_shared[i] << " ";
//    std::cout << std::endl;
    
    
    // WF_other contains relaxed data, does not need to be in window now
    WF_other_data = new double[size_shared_grid*DIM_OTHER];
    WF_other = new Waveform(size_shared_grid, DIM_OTHER, times_shared, WF_other_data);
    WF_other -> set_last(u0_other);

    // Window & Waveform for receiving data
    MPI_Win_allocate(WF_LEN_OTHER * DIM_OTHER * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &WF_recv_data, &WIN_data);
    WF_recv = new Waveform(WF_LEN_OTHER, DIM_OTHER, times_other, WF_recv_data);
//    WF_recv = new Waveform_locking(WF_LEN_OTHER, DIM_OTHER, times_other, WF_recv_data, &WIN_data, ID_SELF);
    // flags for which data has already been received
    MPI_Win_allocate(WF_LEN_OTHER * sizeof(bool), sizeof(bool), MPI_INFO_NULL, MPI_COMM_WORLD, &WF_other_data_recv_flag, &WIN_recv_flag);
    // init flags
    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, ID_SELF, 0, WIN_recv_flag); // Excl, writing
    WF_other_data_recv_flag[0] = true;
    for (int i = 1; i < WF_LEN_OTHER; i++)
        WF_other_data_recv_flag[i] = false;
    MPI_Win_unlock(ID_SELF, WIN_recv_flag);

    // relaxation flags now live on shared grid
    // flags for when other process did GS relax
    relax_other_done_flag = new bool[size_shared_grid];
    for (int i = 0; i < size_shared_grid; i++)
        relax_other_done_flag[i] = false;
    // flags for marking which entry already received relaxation (GS)
    relax_self_done_flag = new bool[size_shared_grid];
    for (int i = 0; i < size_shared_grid; i++)
        relax_self_done_flag[i] = false;
    // flags for marking which entry already received relaxation (jac)
    relax_self_done_flag_jac = new bool[size_shared_grid];
    relax_self_done_flag_jac[0] = true;
    for (int i = 1; i < size_shared_grid; i++)
        relax_self_done_flag_jac[i] = false;
        
    relax_flag_jac_mark = new bool[size_shared_grid];
    for (int i = 0; i < size_shared_grid; i++)
        relax_flag_jac_mark[i] = false;
    
    
    log_errors = errlogging;
    init_error_log(steps_macro, WF_MAX_ITER);

    double window_length = _t_end/steps_macro;
    norm_factor = prob_self -> get_norm_factor(); // implicitly assumed to be identical for both subproblems
    
    MPI_Barrier(MPI_COMM_WORLD);
    runtime = MPI_Wtime(); // runtime measurement start
    for(int i = 0; i < steps_macro; i++){ // Macro step loop
//        std::cout << ID_SELF << " starting macrostepping loop, " << i << std::endl;
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

void WFR_NEW_var_relax_MR::do_WF_iter(double WF_TOL, int WF_MAX_ITER, int steps_per_window_self, int steps_per_window_other){
//    std::cout << ID_SELF << "relaxation parameters " << theta_relax << " " << theta_relax_self_ahead << " " << theta_relax_other_ahead << std::endl;
    first_iter = true;
    steps_converged = 0;
    for(int i = 0; i < WF_MAX_ITER; i++){ // Waveform loop
        MPI_Barrier(MPI_COMM_WORLD);
        WF_iters++;

        /////////////////// STOP MARKER
//        std::cout << ID_SELF << " integrate window " << i << std::endl;
        integrate_window(WF_self, WF_other, steps_per_window_self, prob_self);
        
//        std::cout << ID_SELF << " integrate window done " << i << std::endl;

        // relax_self_done_flag == true marks GS relaxation on the other process
        // having this be of length size_shared_grid might not be necessary, but it is simple for now
        MPI_Sendrecv(relax_self_done_flag, size_shared_grid, MPI_C_BOOL, ID_OTHER, TAG_MISC,
                     relax_other_done_flag, size_shared_grid, MPI_C_BOOL, ID_OTHER, TAG_MISC,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // all outstanding RMA (send) operations done up to this point, since SendRecv is passed
        // sync data
//        std::cout << ID_SELF << " integrate window done 0001" << std::endl;
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, ID_SELF, 0, WIN_data); // Excl since it is updating self
        MPI_Win_sync(WIN_data); // lock neccessary in intel MPI
        MPI_Win_unlock(ID_SELF, WIN_data);
        // necessary to also synch flags? Should not, since flags not used any longer
        // WF_recv and WF_other_data_recv_flag should now be up to date

//        std::cout << ID_SELF << " integrate window done 0002" << std::endl;
        // do all outstanding relaxation
        MPI_Win_lock(MPI_LOCK_SHARED, ID_SELF, 0, WIN_data);
        // do relaxation on shared grid
        int idx_recv_grid = 0; // idx for keeping track of corresponding index on receiving grid
//        double theta_tmp = 1; 
        for (int i = 1; i < size_shared_grid; i++){
//            std::cout << ID_SELF << " error hunting 0001 " << i << std::endl;
            if (times_shared_flags[i] & FLAG_GRID_OTHER) // keeping track of index on receiving grid
                idx_recv_grid++;
            // skip steps that where relaxation was already done
//            std::cout << ID_SELF << " error hunting 0002 " << i << std::endl;
            if (relax_self_done_flag[i] or relax_self_done_flag_jac[i]){
//                std::cout << ID_SELF << " error hunting 0003a " << i << std::endl;
                continue;
            }else{ // relaxation not done yet
//                std::cout << ID_SELF << " error hunting 0003b " << i << std::endl;
                // TODO: a bit of redundancy here, instead check relax type first and set temporary relax parameter instead
//                std::cout << ID_SELF << " error hunting 0004b " << i << std::endl;
                if (relax_other_done_flag[i]) // self ahead -> GS relax
                    theta_tmp = theta_relax_self_ahead;
                else // jacobi relax
                    theta_tmp = theta_relax;
//                std::cout << ID_SELF << " CLEANUP RELAX " << i << " theta_tmp " << theta_tmp << std::endl;
//                std::cout << ID_SELF << " error hunting 0005b " << i << std::endl; 
                if (times_shared_flags[i] & FLAG_GRID_OTHER){ // point on other time-grid, no interpolation necessary
                    for (int j = 0; j < DIM_OTHER; j++)
                        WF_other_data[i*DIM_OTHER + j] = (1 - theta_tmp)*WF_other_data[i*DIM_OTHER + j] + theta_tmp*WF_recv_data[idx_recv_grid*DIM_OTHER + j];
                }else{ // interpolation required
                    // write interpolation result into temporary vector
    //                std::cout << ID_SELF << " error hunting 0006b " << i << std::endl;
    //                std::cout << ID_SELF << " max len " << WF_LEN_OTHER << " i " << i << " times_shared " << times_shared[i] << std::endl;
                    WF_recv->eval(times_shared[i], relax_interpol_aux);
    //                std::cout << ID_SELF << " error hunting 0007b " << i << std::endl; 
                    // do relaxation
                    for (int j = 0; j < DIM_OTHER; j++)
                        WF_other_data[i*DIM_OTHER + j] = (1 - theta_tmp)*WF_other_data[i*DIM_OTHER + j] + theta_tmp*relax_interpol_aux[j];
                }
            }
        }
        MPI_Win_unlock(ID_SELF, WIN_data);

        // reset flags
//        std::cout << ID_SELF << " integrate window done 0003" << std::endl;
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, ID_SELF, 0, WIN_recv_flag); // Excl, writing
        for (int i = 1; i < WF_LEN_OTHER; i++)
            WF_other_data_recv_flag[i] = false;
        MPI_Win_unlock(ID_SELF, WIN_recv_flag);

        for (int i = 0; i < size_shared_grid; i++)
            relax_self_done_flag[i] = false;
        for (int i = 1; i < size_shared_grid; i++)
            relax_self_done_flag_jac[i] = false;

        // since own data is discarded after each iteration over a given time-window, relaxation is only done at receiver on other
        // the exception is the convergence check, for which one need the relaxed data
//        std::cout << ID_SELF << " integrate window done 0004" << std::endl;
        MPI_Sendrecv((*WF_other)[size_shared_grid-1], DIM_OTHER , MPI_DOUBLE, ID_OTHER, TAG_DATA, 
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

void WFR_NEW_var_relax_MR::integrate_window(Waveform * WF_calc, Waveform * WF_src, int steps, Problem * p){
// WF_self, WR_other(shared grid), steps_self, prob_self
    double t, dt;

    int idx_latest_timepoint_recv = 0;
    
    // on receiving grid
    int idx_relevant_interval_left = 0;
    int idx_relevant_interval_right = 0;
//    std::cout << ID_SELF << " WF_other_data_recv_flags: ";
//    for( int j = 0; j < WF_LEN_OTHER; j++)
//        std::cout << WF_other_data_recv_flag[j] << " ";
//    std::cout << std::endl;
    
	for(int i = 0; i < steps; i++){ // timestepping loop
//        std::cout << ID_SELF << " timestep " << i << std::endl;
	
        t = WF_self->get_time(i);
        dt = WF_self->get_time(i+1) - t;
        // first sync flags, other way around a flag might be true, but data might not be there yet (flag updated during data sync.)
        // alternative: enclose locks?
        // synch flags marking which data is new
//        std::cout << ID_SELF << " timestep " << i << "  001"<< std::endl;
        MPI_Win_lock(MPI_LOCK_SHARED, ID_SELF, 0, WIN_recv_flag); // Shared as it is read only
        MPI_Win_sync(WIN_recv_flag); // lock neccessary in intel MPI
        MPI_Win_unlock(ID_SELF, WIN_recv_flag);
        // sync data exposed by window, i.e. catch new updates
//        std::cout << ID_SELF << " timestep " << i << "  002"<< std::endl;
        MPI_Win_lock(MPI_LOCK_SHARED, ID_SELF, 0, WIN_data); // Shared as it is read only
        MPI_Win_sync(WIN_data); // lock neccessary in intel MPI
        MPI_Win_unlock(ID_SELF, WIN_data);

        // \todo: determine the truth values of the if-cases beforehand to reduce total locking time on flags
//        MPI_Win_lock(MPI_LOCK_SHARED, ID_SELF, 0, WIN_recv_flag); // read only
        // determine if other process is ahead, assuming implicit time-integration
        
        // 1. find time-point corresponding to latest available data point from receiving WF
//        int idx_latest_timepoint_recv = 1;
        MPI_Win_lock(MPI_LOCK_SHARED, ID_SELF, 0, WIN_recv_flag); // Shared as it is read only
//        std::cout << ID_SELF << " idx latest checking " << idx_latest_timepoint_recv << " WF_LEN_OTHER " << WF_LEN_OTHER << std::endl;
//        std::cout << ID_SELF << " WF_other_data_recv_flags: ";
//        for( int j = 0; j < WF_LEN_OTHER; j++)
//            std::cout << WF_other_data_recv_flag[j] << " ";
//        std::cout << std::endl;
        while (WF_other_data_recv_flag[idx_latest_timepoint_recv] and (idx_latest_timepoint_recv < WF_LEN_OTHER))
            idx_latest_timepoint_recv++;
        idx_latest_timepoint_recv--;
        MPI_Win_unlock(ID_SELF, WIN_recv_flag);
            
        // 2. find endpoint of current timesteps
        int idx_new_timepoint = i + 1;
        
        // determine relevant interval
        while (times_shared[mapping_recv_to_shared[idx_relevant_interval_left]] <= t)
            idx_relevant_interval_left++;
        idx_relevant_interval_left--;
            
        idx_relevant_interval_right = idx_relevant_interval_left + 1;
        while (times_shared[mapping_recv_to_shared[idx_relevant_interval_right]] < t + dt)
            idx_relevant_interval_right++;
        int idx_shared_grid_startpoint = mapping_recv_to_shared[idx_relevant_interval_left];
        int idx_shared_grid_endpoint = mapping_recv_to_shared[idx_relevant_interval_right];
        
        // compare on shared grid to see if other process is ahead
//        std::cout << ID_SELF << " idx latest recv " << idx_latest_timepoint_recv << " at t = " << mapping_recv_to_shared[idx_latest_timepoint_recv] << " new t = " << mapping_self_to_shared[idx_new_timepoint] << std::endl;
        bool other_process_ahead = (mapping_recv_to_shared[idx_latest_timepoint_recv] >= mapping_self_to_shared[idx_new_timepoint]);
        if (other_process_ahead){
//            std::cout << ID_SELF << " other process ahead " << std::endl;
        
            // determine relevant interval, by indices
            // left side, maximal point of recv grid <= t
//            while (times_shared[mapping_recv_to_shared[idx_relevant_interval_left]] <= t)
//                idx_relevant_interval_left++;
//            idx_relevant_interval_left--;
//                
//            idx_relevant_interval_right = idx_relevant_interval_left + 1;
//            while (times_shared[mapping_recv_to_shared[idx_relevant_interval_right]] < t + dt)
//                idx_relevant_interval_right++;
                
//            std::cout << ID_SELF << " limits " << idx_relevant_interval_left << " " << idx_relevant_interval_right << std::endl;
            
            // a) check relaxation status of previous received timestep, a.k.a. left limit of relevant interval
//            int idx_shared_grid_startpoint = mapping_recv_to_shared[idx_relevant_interval_left];
            
//            bool previous_recv_timestep_relaxed = relax_self_done_flag[idx_shared_grid_startpoint] || relax_self_done_flag_jac[idx_shared_grid_startpoint];
            bool previous_recv_timestep_relaxed = relax_self_done_flag[idx_shared_grid_startpoint] || relax_self_done_flag_jac[idx_shared_grid_startpoint];
            if (not previous_recv_timestep_relaxed){ // repeat JAC relaxation for said timepoint
//                std::cout << ID_SELF << " JAC RELAX of previous step " << idx_shared_grid_startpoint << std::endl;
                MPI_Win_lock(MPI_LOCK_SHARED, ID_SELF, 0, WIN_data); // lock receving data window
//                int idx_jac_relax_recv = idx_latest_timepoint_recv - 1;
//                int idx_jac_relax_shared = mapping_recv_to_shared[idx_jac_relax_recv];
                for (int j = 0; j < DIM_OTHER; j++)
                    WF_other_data[idx_shared_grid_startpoint*DIM_OTHER + j] = (1 - theta_relax)*WF_other_data[idx_shared_grid_startpoint*DIM_OTHER + j] + theta_relax*WF_recv_data[idx_shared_grid_startpoint*DIM_OTHER + j];
                MPI_Win_unlock(ID_SELF, WIN_data);
                relax_self_done_flag_jac[idx_shared_grid_startpoint] = true;
            } // else: relaxation for said time-point already done
            
            // b) do interpolation -> relaxation for all remaining points
//            int idx_shared_grid_endpoint = mapping_recv_to_shared[idx_relevant_interval_right];
            MPI_Win_lock(MPI_LOCK_SHARED, ID_SELF, 0, WIN_data);
            
            // points from [idx_shared_grid_startpoint+1, idx_shared_grid_endpoint) need interpolation
            // loop might be empty
            for (int j = idx_shared_grid_startpoint + 1; j < idx_shared_grid_endpoint; j++){
//                std::cout << ID_SELF << " GS RELAX OF STEP " << j << std::endl;
                // write interpolation result into temporary vector
                WF_recv->eval(times_shared[j], relax_interpol_aux);
                if (relax_flag_jac_mark[j])
                    theta_tmp = theta_relax;
                else
                    theta_tmp = theta_relax_other_ahead;
                // do relaxation
                for (int k = 0; k < DIM_OTHER; k++)
                    WF_other_data[j*DIM_OTHER + k] = (1 - theta_tmp)*WF_other_data[j*DIM_OTHER + k] + theta_tmp*relax_interpol_aux[k];
//                    WF_other_data[j*DIM_OTHER + k] = (1 - theta_relax_other_ahead)*WF_other_data[j*DIM_OTHER + k] + theta_relax_other_ahead*relax_interpol_aux[k];
                if (relax_flag_jac_mark[j])
                    relax_self_done_flag_jac[j] = true;
                else
                    relax_self_done_flag[j] = true;
            }
            // right endpoint does not require interpolation, since the time-point is on the shared grid
//            std::cout << ID_SELF << " GS RELAX OF STEP (endpoint)" << idx_shared_grid_endpoint << std::endl;

            if (relax_flag_jac_mark[idx_shared_grid_endpoint])
                theta_tmp = theta_relax;
            else
                theta_tmp = theta_relax_other_ahead;
            for (int k = 0; k < DIM_OTHER; k++)
                WF_other_data[idx_shared_grid_endpoint*DIM_OTHER + k] = (1 - theta_tmp)*WF_other_data[idx_shared_grid_endpoint*DIM_OTHER + k] + theta_tmp*WF_recv_data[idx_relevant_interval_right*DIM_OTHER + k];
//                WF_other_data[idx_shared_grid_endpoint*DIM_OTHER + k] = (1 - theta_relax_other_ahead)*WF_other_data[idx_shared_grid_endpoint*DIM_OTHER + k] + theta_relax_other_ahead*WF_recv_data[idx_relevant_interval_right*DIM_OTHER + k];
                // mark point as relaxed, on shared grid
            if (relax_flag_jac_mark[idx_shared_grid_endpoint])
                relax_self_done_flag_jac[idx_shared_grid_endpoint] = true;
            else
                relax_self_done_flag[idx_shared_grid_endpoint] = true;
//            relax_self_done_flag[idx_shared_grid_endpoint] = true;
            MPI_Win_unlock(ID_SELF, WIN_data);
        }else{ // other process not ahead, mark relevant time-point for tentative Jacobi relaxation
            for (int j = idx_shared_grid_startpoint; j < idx_shared_grid_endpoint; j++)
                relax_flag_jac_mark[j] = true;
        }
        
        // make the choice of which waveform to choose from and use linear interpolation			
        // i = number of steps done by own process
//        t = WF_self->get_time(i);
//        dt = WF_self->get_time(i+1) - t;
        
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