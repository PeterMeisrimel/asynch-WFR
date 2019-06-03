/*
Authors: Peter Meisrimel
May 2019
*/

#include "WFR.h"
#include "WFR_GS.h"
#include "WFR_JAC.h"
#include "WFR_NEW.h"
#include "problem_toy.h"
#include "input_reader.h"
#include "mpi.h"
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
    bool FIRST = true; // for GS
    // Note, logging formally does not make a lot of sense here, as it is build around implicit time-integration
    bool logging = false;

    // default running parameters
    double WF_TOL = 1e-10;
    int WF_MAXITER = 50;
    int num_macro = 5;
    int u0 = 0;

	process_inputs(argc, argv, runmode, WF_TOL, t_end, timesteps1, timesteps2, num_macro, WF_MAXITER, logging, FIRST, u0);

    t_end /= num_macro;
    num_macro = 1;

    if (ID_SELF == 0){
        timesteps = timesteps1;
        prob = new Problem_toy_part_1();
    }else{
        timesteps = timesteps2;
        prob = new Problem_toy_part_2();
        FIRST = not FIRST;
    }

	switch(runmode){
		case 1:
		    wfr_method = new WFR_GS(ID_SELF, ID_OTHER, t_end, prob, FIRST, true);
			break;
		case 2:
			wfr_method = new WFR_JAC(ID_SELF, ID_OTHER, t_end, prob, true);
			break;
		case 3:
			wfr_method = new WFR_NEW(ID_SELF, ID_OTHER, t_end, prob, true);
			break;
	}
    wfr_method -> run(WF_TOL, WF_MAXITER, num_macro, timesteps, 0);
    
    wfr_method -> write_error_log();

    MPI_Finalize();
	return 0;
}