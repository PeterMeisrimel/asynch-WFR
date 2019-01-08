/*
Authors: Peter Meisrimel
January 2019
*/

#include "input_reader.h"
#include <cassert> // assert
#include <string>  // string handling
#include <stdlib.h> // string to double conversion
#include <stdexcept>

void process_inputs(int argc, char **argv, int& runmode, double& WF_TOL, int& timesteps1, int& timesteps2, int& macrosteps, int& maxiter){
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
		else if(arg == "-wftol")
			WF_TOL = atof(argv[i+1]);
		else if (arg == "-timesteps1")
			timesteps1 = atoi(argv[i+1]);
		else if (arg == "-timesteps2")
			timesteps2 = atoi(argv[i+1]);
		else if (arg == "-macrosteps")
			macrosteps = atoi(argv[i+1]);
		else if (arg == "-maxiter")
			maxiter = atoi(argv[i+1]);
		else
			throw std::invalid_argument("invalid input parameter, check your spelling");
	}
}

void process_inputs_heat(int argc, char **argv , double& alpha1, double& alpha2, double& gamma1, double& gamma2, int& gridsize1, int& gridsize2){
	assert( argc % 2 == 1);

	for(int i = 1; i < argc; i+=2){
		std::string arg = argv[i];
		if(arg == "-alpha1")
			alpha1 = atof(argv[i+1]);
		else if (arg == "-alpha2")
			alpha2 = atof(argv[i+1]);
		else if (arg == "-gamma1")
			gamma1 = atof(argv[i+1]);
		else if (arg == "-gamma2")
			gamma2 = atof(argv[i+1]);
		else if (arg == "-gridsize1")
			gridsize1 = atoi(argv[i+1]);
		else if (arg == "-gridsize2")
			gridsize2 = atoi(argv[i+1]);
		else
			std::invalid_argument("invalid input parameter, check your spelling");
	}
}
