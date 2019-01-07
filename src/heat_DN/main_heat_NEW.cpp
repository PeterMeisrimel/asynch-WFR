#include "WFR_NEW.h"
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

  // Init parameters
  double t_end = 1;
  int timesteps;
  Problem * prob;
  WFR * wfr_method;

  int gridsize = 32;
  double alpha = 1;
  double gamma = 0.01;

  if (ID_SELF == 0){
    timesteps = 100;
    prob = new Problem_heat_D(gridsize, alpha, gamma);
  }else{
    timesteps = 100;
    prob = new Problem_heat_N(gridsize, alpha, gamma);
  }
  double * sol = new double[prob->get_length()];

  // running parameters
  double WF_TOL = 1e-6;
  int WF_MAXITER = 30;
  int num_macro = 5;

  wfr_method = new WFR_NEW(ID_SELF, ID_OTHER, t_end, prob);
  wfr_method -> run(WF_TOL, WF_MAXITER, num_macro, timesteps, -2);
  //wfr_method -> get_sol(sol);
  //std::cout << std::setprecision(14) << "id " << ID_SELF << ", solution: "<< *sol << ", Iterations: " << wfr_method -> get_WF_iters() << std::endl;
  std::cout << std::setprecision(14) << "id " << ID_SELF << ", Iterations: " << wfr_method -> get_WF_iters() << std::endl;

	return 0;
}
