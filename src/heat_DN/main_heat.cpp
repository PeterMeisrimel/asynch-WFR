/*
Authors: Peter Meisrimel
December 2018
*/

#include "WFR_GS.h"
#include "WFR_JAC.h"
#include "WFR_NEW.h"
#include "input_reader.h"
#include "problem_heat_D.h"
#include "problem_heat_N.h"
#include "mpi.h"
#include "dolfin.h"
#include <iostream>
#include <iomanip> // set precision
#include "unistd.h"

int main(int argc, char *argv[]){
  init(argc, argv);
  set_log_level(30); // fenics output level: "things that may go boom later"

	int np, ID_SELF, ID_OTHER;
  MPI_Comm_rank(MPI_COMM_WORLD, &ID_SELF);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
	ID_OTHER = (ID_SELF + 1)%np;

  // require exactly 2 processes
  if (np != 2)
    return 1;

  // default init parameters
  double t_end = 1;
	int runmode = 1;
  int timesteps = 100;
	int timesteps1 = timesteps;
	int timesteps2 = timesteps;
  Problem * prob;
  WFR * wfr_method;
  bool FIRST = false; // for GS

  // default running parameters
  double WF_TOL = 1e-6;
  int WF_MAXITER = 30;
  int num_macro = 5;

  int gridsize = 32;
	int gridsize1 = gridsize, gridsize2 = gridsize;
  double alpha = 1;
	double alpha1 = alpha;
	double alpha2 = alpha;
  double gamma = 0.01;
	double gamma1 = gamma;
	double gamma2 = gamma;

	process_inputs(argc, argv, runmode, WF_TOL, timesteps1, timesteps2, num_macro, WF_MAXITER);
	process_inputs_heat(argc, argv, alpha1, alpha2, gamma1, gamma2, gridsize1, gridsize2);

  if (ID_SELF == 0){
    timesteps = timesteps1;
		alpha = alpha1;
		gamma = gamma1;
		gridsize = gridsize1;
    prob = new Problem_heat_D(gridsize, alpha, gamma);
    FIRST = true;
  }else{
    timesteps = timesteps2;
		alpha = alpha2;
		gamma = gamma2;
		gridsize = gridsize2;
    prob = new Problem_heat_N(gridsize, alpha, gamma);
  }
  double * sol = new double[prob->get_length()];

	switch(runmode){
		case 1:
		  wfr_method = new WFR_GS(ID_SELF, ID_OTHER, t_end, prob, FIRST);
			break;
		case 2:
			wfr_method = new WFR_JAC(ID_SELF, ID_OTHER, t_end, prob);
			break;
		case 3:
			wfr_method = new WFR_NEW(ID_SELF, ID_OTHER, t_end, prob);
			break;
	}

  wfr_method -> run(WF_TOL, WF_MAXITER, num_macro, timesteps, -2);
  std::cout << std::setprecision(14) << "id " << ID_SELF << ", Iterations: " << wfr_method -> get_WF_iters() << " runtime: " << wfr_method->get_runtime() << std::endl;

	return 0;
}
