#include "WFR_GS.h"
#include "waveform.h"
#include "problem.h"
#include "math.h" // for sqrt
#include "mpi.h"

WFR_GS::WFR_GS(int id_self, int id_other, double t_end, Problem * p, bool first){
  ID_SELF  = id_self;
  ID_OTHER = id_other;
  _t_end   = t_end;
  prob    = p;
  FIRST    = first;
  WF_iters = 0;
}

void WFR_GS::do_WF_iter(double WF_TOL, int WF_MAX_ITER, int steps_per_window_self, int steps_per_window_other){
  first_iter = true;
  for(int i = 0; i < WF_MAX_ITER; i++){ // WF iter loop
    WF_iters++;
    
    if (FIRST){
      integrate_window(steps_per_window_self);
		  MPI_Sendrecv(WF_self_data , (steps_per_window_self  + 1) * DIM_SELF , MPI_DOUBLE, ID_OTHER, TAG_DATA, 
								   WF_other_data, (steps_per_window_other + 1) * DIM_OTHER, MPI_DOUBLE, ID_OTHER, TAG_DATA,
								   MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }else{ // not FIRST
      MPI_Recv(WF_other_data, (steps_per_window_other + 1) * DIM_OTHER, MPI_DOUBLE, ID_OTHER, TAG_DATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      integrate_window(steps_per_window_self);
      MPI_Send(WF_self_data, (steps_per_window_self + 1) * DIM_SELF, MPI_DOUBLE, ID_OTHER, TAG_DATA, MPI_COMM_WORLD);
    }
    
    if (check_convergence(WF_TOL)){
      break;
    }else{
      WF_self ->get_last(WF_self_last);
      WF_other->get_last(WF_other_last);
      prob    ->reset_to_checkpoint();
    }
  } // END WF iter loop
}
