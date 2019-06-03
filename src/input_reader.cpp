/*
Authors: Peter Meisrimel
January 2019
*/

#include "input_reader.h"
#include <cassert> // assert
#include <string>  // string handling
#include <stdlib.h> // string to double conversion
#include <stdexcept>

void process_inputs(int argc, char **argv, int& runmode, double& WF_TOL, double& t_end, int& timesteps1, int& timesteps2, int& macrosteps, int& maxiter, bool& logging, bool& FIRST, int& which_u0){
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
        else if (arg == "-logging")
            logging = atoi(argv[i+1]);
        else if (arg == "-first")
            FIRST = bool(atoi(argv[i+1]));
        else if (arg == "-u0")
            which_u0 = atoi(argv[i+1]);
	}
}

void process_inputs_heat(int argc, char **argv , double& alpha1, double& alpha2, double& lambda1, double& lambda2, int& gridsize1, int& gridsize2){
	assert( argc % 2 == 1);

	for(int i = 1; i < argc; i+=2){
		std::string arg = argv[i];
		if (arg == "-alpha"){
			alpha1 = atof(argv[i+1]);
			alpha2 = atof(argv[i+1]);
		}
        else if (arg == "-alpha1")
			alpha1 = atof(argv[i+1]);
		else if (arg == "-alpha2")
			alpha2 = atof(argv[i+1]);
		else if (arg == "-lambda"){
			lambda1 = atof(argv[i+1]);
			lambda2 = atof(argv[i+1]);
		}
        else if (arg == "-lambda1")
			lambda1 = atof(argv[i+1]);
		else if (arg == "-lambda2")
			lambda2 = atof(argv[i+1]);
		else if (arg == "-gridsize"){
			gridsize1 = atoi(argv[i+1]);
			gridsize2 = atoi(argv[i+1]);
        }
        else if (arg == "-gridsize1")
			gridsize1 = atoi(argv[i+1]);
		else if (arg == "-gridsize2")
			gridsize2 = atoi(argv[i+1]);
	}
}
