/*
Authors: Peter Meisrimel
January 2019
*/

#include "WFR.h"
#include "mpi.h"
#include <iostream>
#include <iomanip> // set precision
#include "math.h" // sqrt

void WFR_serial::set_conv_check_WF_ptr(int conv_which, bool match_which_conv_relax){
    RELAX_0 = RELAX;
    RELAX_1 = RELAX;
    switch(conv_which){
        case -1:{ // use ouput of first model (ID_SELF == 0) as measurement for convergence
            WF_conv_check = WF_self;
            WF_conv_check_last = WF_self_last;
            if (match_which_conv_relax) // only relaxation on output of first model
                RELAX_1 = false; // => deactivate relaxation for 1st processor
            break;
        }
        case -2:{ // use ouput of second model (ID_SELF == 1) as measurement for convergence
            WF_conv_check = WF_other;
            WF_conv_check_last = WF_other_last;
            if (match_which_conv_relax) // only relaxation on output of second model
                RELAX_0 = false; // => deactive relaxation for 2nd processor
            break;
        }
    }
}

void WFR_parallel::set_conv_check_WF_ptr(int conv_which, bool match_which_conv_relax){
    switch(conv_which){
        case -1:{ // use ouput of first model (ID_SELF == 0) as measurement for convergence
            if (ID_SELF == 0){
                WF_conv_check = WF_self;
                WF_conv_check_last = WF_self_last;
            }else{
                WF_conv_check = WF_other;
                WF_conv_check_last = WF_other_last;
                if (match_which_conv_relax) // only relaxation for output of first model
                    RELAX = false; // deactivate relaxation for ID_SELF == 1
            }
            break;   
        }
        case -2:{ // use ouput of second model (ID_SELF == 1) as measurement for convergence
            if (ID_SELF == 0){
                WF_conv_check = WF_other;
                WF_conv_check_last = WF_other_last;
                if (match_which_conv_relax) // only relaxation on output of second model
                    RELAX = false; // => deactive relaxation for output of first model
                }else{
                    WF_conv_check = WF_self;
                    WF_conv_check_last = WF_self_last;
                }
            break;
        }
    }
}

void WFR::integrate_window(Waveform * WF_calc, Waveform * WF_src, int steps, Problem * p){
    double t;
    if (RELAX){
        for(int i = 0; i < steps; i++){
            t = WF_calc->get_time(i);
            WF_calc -> get(i+1, relax_aux_vec); // store previous iterate at new timepoint 
            p->do_step(t, WF_calc->get_time(i+1) - t, (*WF_calc)[i+1], WF_src); // overwrite previous iterate at new timepoint
            WF_calc -> relax_by_single(w_relax, i+1, relax_aux_vec); // do relaxation, relax_aux_vec being the old iterate
        }
    }else{
        for(int i = 0; i < steps; i++){
            t = WF_calc->get_time(i);
            p->do_step(t, WF_calc->get_time(i+1) - t, (*WF_calc)[i+1], WF_src);
        }
    }
}

void WFR::integrate_window(int start, Waveform * WF_calc, Waveform * WF_src, int steps, Problem * p){
    double t;
    if (RELAX){
        for(int i = start; i < start + steps; i++){
            t = WF_calc->get_time(i);
            WF_calc -> get(i+1, relax_aux_vec); // store previous iterate at new timepoint 
            p->do_step(t, WF_calc->get_time(i+1) - t, (*WF_calc)[i+1], WF_src); // overwrite previous iterate at new timepoint
            WF_calc -> relax_by_single(w_relax, i+1, relax_aux_vec); // do relaxation, relax_aux_vec being the old iterate
        }
    }else{
        for(int i = start; i < start + steps; i++){
            t = WF_calc->get_time(i);
            p->do_step(t, WF_calc->get_time(i+1) - t, (*WF_calc)[i+1], WF_src);
        }
    }
}

void WFR::get_relative_tol(){
    switch(conv_which){
        case 0:{
            double val0 = WF_self ->get_norm_sq_last();
            double val1 = WF_other->get_norm_sq_last();
            rel_update_fac = sqrt(val0 + val1);
            break;
        }
        case 1:{
            double w_self  = float(DIM_SELF)  / float(DIM_SELF + DIM_OTHER);
            double w_other = float(DIM_OTHER) / float(DIM_SELF + DIM_OTHER);
            double val0 = WF_self ->get_norm_sq_last();
            double val1 = WF_other->get_norm_sq_last();
            rel_update_fac = sqrt(w_self*val0 + w_other*val1);
            break;
        }
        default:{
            double val;
            if ((conv_which == -1) || (conv_which == -2)){
                val = WF_conv_check -> get_norm_sq_last();
                rel_update_fac = sqrt(val);
                break;
            }else
                throw std::invalid_argument("No method to check for convergence implemented for this input");
        }
    }
}

bool WFR::check_convergence(double WF_TOL){
    update_error_log();
    if (first_iter){
        first_iter = false;
    }else{
        switch(conv_which){
            case 0:{ // usual 2-norm
                up_self  = WF_self ->get_err_norm_sq_last(WF_self_last);
                up_other = WF_other->get_err_norm_sq_last(WF_other_last);
                update   = sqrt(up_self + up_other);
                break;
            }
            case 1:{ // weighted scale norm
                double w_self  = float(DIM_SELF)  / float(DIM_SELF + DIM_OTHER);
                double w_other = float(DIM_OTHER) / float(DIM_SELF + DIM_OTHER);
                up_self  = WF_self ->get_err_norm_sq_last(WF_self_last);
                up_other = WF_other->get_err_norm_sq_last(WF_other_last);
                update   = sqrt(w_self*up_self + w_other*up_other);
                break;
            }
            default:{
                if ((conv_which == -1) || (conv_which == -2)){
                    up_self  = WF_conv_check -> get_err_norm_sq_last(WF_conv_check_last);
                    update = sqrt(up_self);
                    break;
                }else
                    throw std::invalid_argument("No method to check for convergence implemented for this input");
            }
        }
//        std::cout << "update " << update << std::endl;
        if (update/rel_update_fac < WF_TOL)
            return true;
    }// END ELSE
    return false;
}

void WFR_serial::write_results(){
    sol_size = prob_self->get_length();
    sol = new double[sol_size];
    WF_self -> get_last(sol);

    sol_other_size = prob_other->get_length();
    sol_other = new double[sol_other_size];
    WF_other -> get_last(sol_other);

    runtime_self = get_runtime();
    runtime_other = runtime_self;
    
    WFR::write_results();
}

void WFR_parallel::write_results(){
    sol_size = prob_self->get_length();
    sol = new double[sol_size];
    WF_self -> get_last(sol);

    MPI_Sendrecv(&sol_size, 1, MPI_INT, ID_OTHER, 0,
                 &sol_other_size, 1, MPI_INT, ID_OTHER, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    sol_other = new double[sol_other_size];
    MPI_Sendrecv(sol, sol_size, MPI_DOUBLE, ID_OTHER, 1,
                 sol_other, sol_other_size, MPI_DOUBLE, ID_OTHER, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    runtime_self = get_runtime();
    runtime_other = 0;
    MPI_Sendrecv(&runtime_self,  1, MPI_DOUBLE, ID_OTHER, 2,
                 &runtime_other, 1, MPI_DOUBLE, ID_OTHER, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    WFR::write_results();
}


void WFR::write_results(){
    iters = get_WF_iters();
	if(ID_SELF == 0){
        std::cout << std::setprecision(14) << ID_SELF << " " << iters << " " << runtime_self << " ";
	    for(int i = 0; i < sol_size; i++)
            std::cout << std::setprecision(14) << sol[i] << " ";
        std::cout << std::endl;

        std::cout << std::setprecision(14) << ID_OTHER << " " << iters << " " << runtime_other << " ";
	    for(int i = 0; i < sol_other_size; i++)
            std::cout << std::setprecision(14) << sol_other[i] << " ";
        std::cout << std::endl;
    }
}

void WFR_serial::init_error_log(int steps_macro, int WF_MAX_ITER){
    if (log_errors){
        if (steps_macro != 1)
            MPI_Abort(MPI_COMM_WORLD, 2);
        error_log = new double[(WF_MAX_ITER+1)*DIM_SELF];
        WF_self -> get_last(error_log);

        error_other_log = new double[(WF_MAX_ITER+1)*DIM_OTHER];
        WF_other -> get_last(error_other_log);

        err_log_counter++; // only one needed since error logging is based on macro steps rather than timesteps
    }
}

void WFR_parallel::init_error_log(int steps_macro, int WF_MAX_ITER){
    if (log_errors){
        if (steps_macro != 1)
            MPI_Abort(MPI_COMM_WORLD, 2);

        error_log = new double[(WF_MAX_ITER+1)*DIM_SELF];
        WF_self -> get_last(error_log);
        err_log_counter++;
    }
}

void WFR_serial::update_error_log(){
    if (log_errors){
        WF_self  -> get_last(error_log       + DIM_SELF *err_log_counter);
        WF_other -> get_last(error_other_log + DIM_OTHER*err_log_counter);
        err_log_counter++;
    }
}

void WFR_parallel::update_error_log(){
    if (log_errors){
        WF_self -> get_last(error_log + DIM_SELF*err_log_counter);
        err_log_counter++;
    }
}

// void WFR_serial::write_error_log() // not required

void WFR_parallel::write_error_log(){
    if (log_errors){
        error_other_log = new double[err_log_counter*DIM_OTHER];

        MPI_Sendrecv(error_log, err_log_counter*DIM_SELF, MPI_DOUBLE, ID_OTHER, 1,
                     error_other_log, err_log_counter*DIM_OTHER, MPI_DOUBLE, ID_OTHER, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        WFR::write_error_log();
    }
}

void WFR::write_error_log(){
    if (log_errors){
	    if(ID_SELF == 0){
            for (int j = 0; j < err_log_counter; j++){
                std::cout << std::setprecision(14) << ID_SELF << " " << j << " ";
	            for(int i = 0; i < DIM_SELF; i++)
                    std::cout << std::setprecision(14) << error_log[i + j*DIM_SELF] << " ";
                std::cout << std::endl;

                std::cout << std::setprecision(14) << ID_OTHER << " " << j << " ";
	            for(int i = 0; i < DIM_OTHER; i++)
                    std::cout << std::setprecision(14) << error_other_log[i + j*DIM_OTHER] << " ";
                std::cout << std::endl;
            }
        }
    }
}
