/*
Authors: Peter Meisrimel
January 2019
*/

#include "WFR.h"
#include "mpi.h"
#include <iostream>
#include <iomanip> // set precision
#include <cassert>
#include "math.h" // sqrt

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
        case -1:{
            double val;
            if (ID_SELF == 0)
                val = WF_self  ->get_norm_sq_last();
            else
                val = WF_other ->get_norm_sq_last();
            rel_update_fac = sqrt(val);
            break;
        }
        case -2:{
            double val;
            if (ID_SELF == 1)
                val = WF_self  ->get_norm_sq_last();
            else
                val = WF_other ->get_norm_sq_last();
            rel_update_fac = sqrt(val);
            break;
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
            case -1:{
                if (ID_SELF == 0)
                  up_self  = WF_self  ->get_err_norm_sq_last(WF_self_last);
                else
                  up_self  = WF_other ->get_err_norm_sq_last(WF_other_last);
                update = sqrt(up_self);
                break;
            }
            case -2:{
                if (ID_SELF == 1)
                  up_self  = WF_self  ->get_err_norm_sq_last(WF_self_last);
                else
                  up_self  = WF_other ->get_err_norm_sq_last(WF_other_last);
                update = sqrt(up_self);
                break;
            }
            default:{
                throw std::invalid_argument("No method to check for convergence implemented for this input");
            }
        }
        if (update/rel_update_fac < WF_TOL)
            return true;
    }// END ELSE
    return false;
}

void WFR::write_results(){
    int sol_size, sol_other_size, iters;
    double * sol, *sol_other;
    double runtime_self, runtime_other;


    sol_size = prob->get_length();
    sol = new double[sol_size];
    get_sol(sol);

    MPI_Sendrecv(&sol_size, 1, MPI_INT, ID_OTHER, 0,
                 &sol_other_size, 1, MPI_INT, ID_OTHER, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    sol_other = new double[sol_other_size];
    MPI_Sendrecv(sol, sol_size, MPI_DOUBLE, ID_OTHER, 1,
                 sol_other, sol_other_size, MPI_DOUBLE, ID_OTHER, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    runtime_self = get_runtime();
    runtime_other = 0;
    MPI_Sendrecv(&runtime_self, 1, MPI_DOUBLE, ID_OTHER, 2,
                 &runtime_other, 1, MPI_DOUBLE, ID_OTHER, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

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

void WFR::init_error_log(int num_macro, int WF_MAX_ITER){
    if (log_errors){
        assert(num_macro == 1);
        error_log = new double[(WF_MAX_ITER+1)*DIM_SELF];
		
        WF_self -> get_last(error_log);
        err_log_counter++;
    }
}

void WFR::update_error_log(){
    if (log_errors){
        WF_self -> get_last(error_log + DIM_SELF*err_log_counter);
        err_log_counter++;
    }
}

void WFR::write_error_log(){
    if (log_errors){
        int sol_other_size;
        double * sol_other = new double[err_log_counter*DIM_OTHER];

        MPI_Sendrecv(error_log, err_log_counter*DIM_SELF, MPI_DOUBLE, ID_OTHER, 1,
                     sol_other, err_log_counter*DIM_OTHER, MPI_DOUBLE, ID_OTHER, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    if(ID_SELF == 0){
            for (int j = 0; j < err_log_counter; j++){
                std::cout << std::setprecision(14) << ID_SELF << " " << j << " ";
	            for(int i = 0; i < DIM_SELF; i++)
                    std::cout << std::setprecision(14) << error_log[i + j*DIM_SELF] << " ";
                std::cout << std::endl;

                std::cout << std::setprecision(14) << ID_OTHER << " " << j << " ";
	            for(int i = 0; i < DIM_OTHER; i++)
                    std::cout << std::setprecision(14) << sol_other[i + j*DIM_OTHER] << " ";
                std::cout << std::endl;
            }
        }
    }
}
