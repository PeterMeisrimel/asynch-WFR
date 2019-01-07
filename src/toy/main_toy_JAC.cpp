/*
Authors: Peter Meisrimel
December 2018
*/

#include "WFR_JAC.h"
#include "problem_toy.h"
#include "mpi.h"
#include <iostream>
#include <iomanip> // set precision
#include "unistd.h"

int main(int argc, char *argv[]){
	MPI_Init(&argc, &argv);

	int np, ID_SELF, ID_OTHER;
  MPI_Comm_rank(MPI_COMM_WORLD, &ID_SELF);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
	ID_OTHER = (ID_SELF + 1)%np;

  // require exactly 2 processes
  if (np != 2)
    return 1;

  // Init parameters
  double t_end = 2;
  int timesteps;
  Problem * prob;
  WFR * wfr_method;

  if (ID_SELF == 0){
    timesteps = 100;
    prob = new Problem_toy_part_1();
  }else{
    timesteps = 100;
    prob = new Problem_toy_part_2();
  }
  double * sol = new double[prob->get_length()];

  // running parameters
  double WF_TOL = 1e-10;
  int WF_MAXITER = 50;
  int num_macro = 5;

  wfr_method = new WFR_JAC(ID_SELF, ID_OTHER, t_end, prob);
  wfr_method -> run(WF_TOL, WF_MAXITER, num_macro, timesteps);
  wfr_method -> get_sol(sol);
  std::cout << std::setprecision(14) << "id " << ID_SELF << ", solution: "<< *sol << ", Iterations: " << wfr_method -> get_WF_iters() << std::endl;

  MPI_Finalize();
	return 0;
}
