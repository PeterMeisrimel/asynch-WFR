/*
Authors: Peter Meisrimel
January 2019
*/

#include "utils.h"
#include <string>  // string handling
#include <stdlib.h> // string to double conversion
#include <stdexcept>

#include "WFR.h"
#include "WFR_GS.h"
#include "WFR_JAC.h"
#include "WFR_NEW.h"
#include "mpi.h"

void process_inputs(int argc, char **argv, int& runmode, double& WF_TOL, double& t_end, int& timesteps1, int& timesteps2, int& macrosteps, int& maxiter, bool& FIRST, bool& error_logging, bool& comm_logging, int& nconv){
    if (argc % 2 != 1)
        throw std::invalid_argument("invalid number of input arguments, needs to be even, check for accidental spaces");

	for(int i = 1; i < argc; i+=2){
		std::string arg = argv[i];
		if (arg == "-runmode"){
			std::string arg2 = argv[i+1];
			if (arg2 == "GS")
				runmode = 1;
			else if (arg2 == "JAC")
				runmode = 2;
			else if (arg2 == "NEW")
				runmode = 3;
			else
				throw std::invalid_argument("invalid runmode, check your spelling");
		}
		else if(arg == "-tend")
			t_end = atof(argv[i+1]);
		else if(arg == "-wftol")
			WF_TOL = atof(argv[i+1]);
		else if (arg == "-timesteps1")
			timesteps1 = atoi(argv[i+1]);
		else if (arg == "-timesteps2")
			timesteps2 = atoi(argv[i+1]);
		else if (arg == "-timesteps"){
			timesteps1 = atoi(argv[i+1]);
			timesteps2 = atoi(argv[i+1]);
        }
		else if (arg == "-macrosteps")
			macrosteps = atoi(argv[i+1]);
		else if (arg == "-maxiter")
			maxiter = atoi(argv[i+1]);
        else if (arg == "-commlog")
            comm_logging = bool(atoi(argv[i+1]));
        else if (arg == "-first")
            FIRST = bool(atoi(argv[i+1]));
        else if (arg == "-errlog")
            error_logging = bool(atoi(argv[i+1]));
        else if (arg == "-nconv")
            nconv = atoi(argv[i+1]);
	}
}

void setup_and_run_WFR(Problem * prob1, Problem * prob2, int which_conv, double t_end, int timesteps, int argc, char *argv[]){
	int np, ID_SELF, ID_OTHER;
    MPI_Comm_rank(MPI_COMM_WORLD, &ID_SELF);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
	ID_OTHER = (ID_SELF + 1)%np;

    // default init parameters
	int runmode = 1;
    bool FIRST = true; // for GS
    bool commlogging = false;
    bool errlogging = false;
    int steps_converged_required = 3;

	int timesteps1 = timesteps, timesteps2 = timesteps;

    // default running parameters
    double WF_TOL = 1e-6;
    int WF_MAXITER = 200;
    int num_macro = 5;

	process_inputs(argc, argv, runmode, WF_TOL, t_end, timesteps1, timesteps2, num_macro, WF_MAXITER, FIRST, errlogging, commlogging, steps_converged_required);

    // parallel methods only get one problem as input, pick correct one if there are multiple processes
    Problem * prob;
    if (ID_SELF == 0)
        prob = prob1;
    else
        prob = prob2;

    // build WFR method based on input parameters
    WFR * wfr_method;
	switch(runmode){
		case 1:{ // Gauss-Seidel (GS)
            if (np != 1)
                MPI_Abort(MPI_COMM_WORLD, 1);
		    wfr_method = new WFR_GS(t_end, prob1, prob2, FIRST, errlogging);
			break;
        }
		case 2:{ // Jacobi (JAC)
            if (np != 2)
                MPI_Abort(MPI_COMM_WORLD, 1);
			wfr_method = new WFR_JAC(ID_SELF, ID_OTHER, t_end, prob, errlogging);
			break;
        }
		case 3:{ // NEW method, still needs better name
            if (np != 2)
                MPI_Abort(MPI_COMM_WORLD, 1);
			wfr_method = new WFR_NEW(ID_SELF, ID_OTHER, t_end, prob, errlogging, commlogging);
			break;
        }
	}
    
    wfr_method -> run(WF_TOL, WF_MAXITER, num_macro, timesteps1, timesteps2, which_conv, steps_converged_required);

    if (not errlogging){
        wfr_method -> write_results();
    }else{
        MPI_Barrier(MPI_COMM_WORLD);
        wfr_method -> write_error_log();
    }
}