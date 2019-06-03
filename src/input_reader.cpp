/*
Authors: Peter Meisrimel
January 2019
*/

#include "input_reader.h"
#include <cassert> // assert
#include <string>  // string handling
#include <stdlib.h> // string to double conversion
#include <stdexcept>

void process_inputs(int argc, char **argv, int& runmode, double& WF_TOL, double& t_end, int& timesteps1, int& timesteps2, int& macrosteps, int& maxiter, bool& FIRST, bool& error_logging, bool& comm_logging){
	assert( argc % 2 == 1);

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
	}
}
