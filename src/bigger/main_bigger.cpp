/*
Authors: Peter Meisrimel
December 2018
*/

#include "problem_bigger.h"
#include "utils.h"
#include "mpi.h"

int main(int argc, char *argv[]){
	// the usual MPI initialization
	MPI_Init(&argc, &argv);
	int np, ID_SELF, ID_OTHER;
    MPI_Comm_rank(MPI_COMM_WORLD, &ID_SELF);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
	ID_OTHER = (ID_SELF + 1)%np;

    // Default init parameters for problems
    double t_end = 2;
    int timesteps = 100;

    Problem * prob1, * prob2;

    int which_conv = 0; 

    // IMPORTANT TO INITIALIZE PROBLEMS LIKE THIS
    // ensures that for parallel methods each processor initializes only ones problem and for serial ones, both are initialized
    if (ID_SELF == 0)
        prob1 = new Problem_bigger_part_1();
    if (ID_OTHER == 0)
        prob2 = new Problem_bigger_part_2();

    setup_and_run_WFR(prob1, prob2, which_conv, t_end, timesteps, argc, argv);

    MPI_Finalize();
	return 0;
}
