/*
Authors: Peter Meisrimel
December 2018
*/

#include "WFR.h"
#include "WFR_GS.h"
#include "WFR_JAC.h"
#include "WFR_NEW.h"
#include "input_reader.h"
#include "input_reader_heat.h"
#include "problem_heat_D.h"
#include "problem_heat_N.h"
#include "mpi.h"
#include "dolfin.h"
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
    bool FIRST = true; // for GS
    bool commlogging = false;
    bool errlogging = false;
	int which_u0 = 0;
    int nconv = 1;

    // default running parameters
    double WF_TOL = 1e-6;
    int WF_MAXITER = 200;
    int num_macro = 5;

    int gridsize = 32;
	int gridsize1 = gridsize, gridsize2 = gridsize;
    double alpha = 1;
	double alpha1 = alpha;
	double alpha2 = alpha;
    double lambda = 0.01;
	double lambda1 = lambda;
	double lambda2 = lambda;

	process_inputs(argc, argv, runmode, WF_TOL, t_end, timesteps1, timesteps2, num_macro, WF_MAXITER, FIRST, errlogging, commlogging, nconv);
	process_inputs_heat(argc, argv, alpha1, alpha2, lambda1, lambda2, gridsize1, gridsize2, which_u0);

    if (ID_SELF == 0){
        timesteps = timesteps1;
		alpha = alpha1;
		lambda = lambda1;
		gridsize = gridsize1;
        prob = new Problem_heat_D(gridsize, alpha, lambda, 500, which_u0);
    }else{
        timesteps = timesteps2;
		alpha = alpha2;
		lambda = lambda2;
		gridsize = gridsize2;
        prob = new Problem_heat_N(gridsize, alpha, lambda, 500, which_u0);
        FIRST = not FIRST;
    }

	switch(runmode){
		case 1:
		    wfr_method = new WFR_GS(ID_SELF, ID_OTHER, t_end, prob, FIRST, errlogging, nconv);
			break;
		case 2:
			wfr_method = new WFR_JAC(ID_SELF, ID_OTHER, t_end, prob, errlogging, nconv);
			break;
		case 3:
			wfr_method = new WFR_NEW(ID_SELF, ID_OTHER, t_end, prob, errlogging, commlogging, nconv);
			break;
	}

    wfr_method -> run(WF_TOL, WF_MAXITER, num_macro, timesteps, -2);

    if (not errlogging){
        wfr_method -> write_results();
    }else{
        MPI_Barrier(MPI_COMM_WORLD);
        wfr_method -> write_error_log();
    }

	return 0;
}
