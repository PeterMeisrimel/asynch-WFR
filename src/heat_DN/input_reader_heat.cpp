/*
Authors: Peter Meisrimel
June 2019
*/

#include "input_reader_heat.h"
#include <cassert> // assert
#include <string>  // string handling
#include <stdlib.h> // string to double conversion
#include <stdexcept>

void process_inputs_heat(int argc, char **argv , double& alpha1, double& alpha2, double& lambda1, double& lambda2, int& gridsize1, int& gridsize2, int& which_u0){
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
        else if (arg == "-u0")
            which_u0 = atoi(argv[i+1]);
	}
}
