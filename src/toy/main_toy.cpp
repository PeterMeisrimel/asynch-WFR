/*
Authors: Peter Meisrimel
December 2018
*/

#include "WFR_GS.h"
#include "WFR_JAC.h"
#include "WFR_NEW.h"
#include "problem_toy.h"
#include "input_reader.h"
#include "mpi.h"
#include <iostream>
#include <iomanip> // set precision
#include "unistd.h"

int main(int argc, char *argv[]){
	// the usual MPI initialization
	MPI_Init(&argc, &argv);
	int np, ID_SELF, ID_OTHER;
    MPI_Comm_rank(MPI_COMM_WORLD, &ID_SELF);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
	ID_OTHER = (ID_SELF + 1)%np;

    // require exactly 2 processes
    if (np != 2)
        return 1;

    // Default init parameters for problems
    double t_end = 2;
    int timesteps = 100;
	int timesteps1 = timesteps;
	int timesteps2 = timesteps;
	int runmode = 1;
    Problem * prob;
    WFR * wfr_method;
    bool FIRST = false; // for GS

    // default running parameters
    double WF_TOL = 1e-10;
    int WF_MAXITER = 50;
    int num_macro = 5;

	process_inputs(argc, argv, runmode, WF_TOL, timesteps1, timesteps2, num_macro, WF_MAXITER);

    if (ID_SELF == 0){
        timesteps = timesteps1;
        prob = new Problem_toy_part_1();
        FIRST = true;
    }else{
        timesteps = timesteps2;
        prob = new Problem_toy_part_2();
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
    wfr_method -> run(WF_TOL, WF_MAXITER, num_macro, timesteps, 1);
    wfr_method -> get_sol(sol);
    std::cout << std::setprecision(14) << ID_SELF << " " << wfr_method -> get_WF_iters() << " " << wfr_method -> get_runtime() << " " << *sol << std::endl;

    MPI_Finalize();
	return 0;
}
