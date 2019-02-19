/*
Authors: Peter Meisrimel
December 2018
*/

#include "WFR_synch.h"
#include "waveform.h"
#include "problem.h"
#include <stdexcept>
#include "math.h" // for sqrt
#include "mpi.h"

void WFR_synch::run(double WF_TOL, int WF_MAX_ITER, int num_macro, int steps_self, int conv_check){
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

    WF_other_data = new double[WF_LEN_OTHER * DIM_OTHER];
    WF_other      = new Waveform(WF_LEN_OTHER, DIM_OTHER, times_other, WF_other_data);
    WF_other->set_last(u0_other);

    double window_length = _t_end/num_macro;

	MPI_Barrier(MPI_COMM_WORLD);
	runtime = MPI_Wtime(); // runtime measurement start
    for(int i = 0; i < num_macro; i++){
        prob->create_checkpoint();
        WF_self ->get_last(u0_self);
        WF_self ->set(0, u0_self);
        WF_other->init_by_last();
      
        MPI_Barrier(MPI_COMM_WORLD);
        do_WF_iter(WF_TOL, WF_MAX_ITER, WF_LEN_SELF - 1, WF_LEN_OTHER - 1);
        if(i != num_macro - 1){
            WF_self ->time_shift(window_length);
            WF_other->time_shift(window_length);
        }
    }
    runtime = MPI_Wtime() - runtime; // runtime measurement end
}

void WFR_synch::integrate_window(int steps){
    double t;
    for(int i = 0; i < steps; i++){
        t = WF_self->get_time(i);
        prob->do_step(t, WF_self->get_time(i+1) - t, (*WF_self)[i+1], WF_other);
    }
}

bool WFR_synch::check_convergence(double WF_TOL){
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
        //std::cout << "id: " << ID_SELF << " " << update << std::endl;
        if (update < WF_TOL)
            return true;
    }// END ELSE
    return false;
}
