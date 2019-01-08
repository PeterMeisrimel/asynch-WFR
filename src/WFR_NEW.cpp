/*
Authors: Peter Meisrimel
December 2018
*/

#include "WFR_NEW.h"
#include "waveform.h"
#include "problem.h"
#include <stdexcept>
#include "math.h"
#include "mpi.h"

WFR_NEW::WFR_NEW(int id_self, int id_other, double tend, Problem * p){
  ID_SELF  = id_self;
  ID_OTHER = id_other;
  _t_end = tend;
  prob    = p;
  WF_iters = 0;
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

  WF_other_data = new double[WF_LEN_OTHER * DIM_OTHER];
  WF_other      = new Waveform(WF_LEN_OTHER, DIM_OTHER, times_other, WF_other_data);
  WF_other->set_last(u0_other);

  double window_length = _t_end/num_macro;

	// Waveform iteration specific stuff
  WF_other_data_new = new double[WF_LEN_OTHER * DIM_OTHER];
  WF_other_new      = new Waveform(WF_LEN_OTHER, DIM_OTHER, times_other, WF_other_data_new);
  WF_other_new->set_last(u0_other);

	MPI_Win_allocate(sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &IDX_win, &WIN_idx);
	MPI_Win_allocate(WF_LEN_SELF * DIM_SELF * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &WF_other_data_win, &WIN_data);
  WF_win = new Waveform(WF_LEN_OTHER, DIM_OTHER, times_other, WF_other_data_win); // init with times not necessary in theory

  MPI_Barrier(MPI_COMM_WORLD);
	runtime = MPI_Wtime(); // runtime measurement start
	for(int i = 0; i < num_macro; i++){ // Macro step loop
    prob -> create_checkpoint();
    WF_self ->get_last(u0_self);
    WF_self ->set(0, u0_self);

    WF_other->init_by_last();

    WF_other_new -> get_last(u0_other);
    WF_other_new -> set(0, u0_other);

    MPI_Barrier(MPI_COMM_WORLD);
    do_WF_iter(WF_TOL, WF_MAX_ITER, WF_LEN_SELF - 1, WF_LEN_OTHER - 1);
    if(i != num_macro - 1){
      WF_self  -> time_shift(window_length);
      WF_other -> time_shift(window_length); // NOTE FOR FUTURE: WF_other and WF_other_new share same time vector, hence only 1 shift here
    }
	} // endfor macrostep loop
	runtime = MPI_Wtime() - runtime; // runtime measurement end

	MPI_Win_free(&WIN_data);
	MPI_Win_free(&WIN_idx);
}

void WFR_NEW::do_WF_iter(double WF_TOL, int WF_MAX_ITER, int steps_per_window, int dummy){
	first_iter = true;
	for(int i = 0; i < WF_MAX_ITER; i++){ // Waveform loop
    WF_iters++;
		IDX = 0; // marker how many timesteps other process has done
    // reset own marker, exclusive as writing is done
		MPI_Win_lock(MPI_LOCK_EXCLUSIVE, ID_SELF, 0, WIN_idx);
		IDX_win[0] = 0;
		MPI_Win_unlock(ID_SELF, WIN_idx);

    MPI_Barrier(MPI_COMM_WORLD); // make sure both markers have been properly reset, probably can be put somewhere else?

		msg_to_recv = 0;
    integrate_window(steps_per_window);

    // Yes, this part is actually necessary, for some reason, to force synchronization
		MPI_Sendrecv(&msg_sent, 1, MPI_INT, ID_OTHER, TAG_DONE,
		             &msg_to_recv, 1, MPI_INT, ID_OTHER, TAG_DONE,
		             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		while(true){
  		if(msg_to_recv == IDX){
  			//MPI_Win_lock(MPI_LOCK_SHARED, ID_SELF, 0, WIN_idx); // reset IDX in window
  			//IDX_win[0] = 0;
				//MPI_Win_unlock(ID_SELF, WIN_idx); 
				break;
			}
      check_new_data();
		}

    if (check_convergence(WF_TOL)){
      break;
    }else{
      WF_self      -> get_last(WF_self_last);
      WF_other_new -> get_last(WF_other_last);
      WF_swap_data_pointers(WF_other, WF_other_new);
      prob -> reset_to_checkpoint();
    }
	} // endfor waveform loop
}

void WFR_NEW::integrate_window(int steps){
  double t;
	for(int i = 0; i < steps; i++){ // timestepping loop, explicit euler

    check_new_data();

		// make the choice of which waveform to choose from and use linear interpolation			
    // IDX = number of timesteps done by other process, i = number of steps done by own process
    if(IDX >= i) // other one is further ahead
      WF_curr = WF_other_new;
    else
      WF_curr = WF_other;
		// actual timestep
    t = WF_self->get_time(i);
    prob->do_step(t, WF_self->get_time(i+1) - t, (*WF_self)[i+1], WF_curr);

		msg_sent = i+1;
  	// Send out new data to other process
  	// actual data first, IDX last, Exclusive locks, as data is written
  	MPI_Win_lock(MPI_LOCK_EXCLUSIVE, ID_OTHER, 0, WIN_data);
    MPI_Put((*WF_self)[msg_sent], DIM_SELF, MPI_DOUBLE, ID_OTHER, DIM_SELF * msg_sent, DIM_SELF, MPI_DOUBLE, WIN_data);
		MPI_Win_unlock(ID_OTHER, WIN_data);
		// IDX
		MPI_Win_lock(MPI_LOCK_EXCLUSIVE, ID_OTHER, 0, WIN_idx);
		MPI_Put(&msg_sent, 1, MPI_INT, ID_OTHER, 0, 1, MPI_INT, WIN_idx);
		MPI_Win_unlock(ID_OTHER, WIN_idx);
	}
}

void WFR_NEW::check_new_data(){
 	// check for new messages by index, lock own window and check for updates
	IDX_aux = IDX;
	MPI_Win_lock(MPI_LOCK_SHARED, ID_SELF, 0, WIN_idx); // Shared as it is read only
	IDX = IDX_win[0];
	MPI_Win_unlock(ID_SELF, WIN_idx);
	// if new data available, lock data window and copy data to local
	if(IDX > IDX_aux){
		MPI_Win_lock(MPI_LOCK_SHARED, ID_SELF, 0, WIN_data); // Shared as it is read only
    for(int h = (IDX_aux + 1); h <= IDX; h++) // start at old max IDX + 1 until (inclusive) new index
      WF_other_new -> set(h, (*WF_win)[h]);
		MPI_Win_unlock(ID_SELF, WIN_data);
	}
}

bool WFR_NEW::check_convergence(double WF_TOL){
  if(first_iter){ // update waveforms
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
		if(update < WF_TOL)
			return true;
	} // end elif first_iter
  return false;
}
