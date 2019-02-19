/*
Authors: Peter Meisrimel
January 2019
*/

#include "WFR.h"
#include "mpi.h"
#include <iostream>
#include <iomanip> // set precision

void WFR::write_results(){
    int sol_size, sol_other_size, iters;
    double * sol, *sol_other;
    double runtime_self, runtime_other;


    sol_size = prob->get_length();
    sol = new double[sol_size];
    get_sol(sol);

    MPI_Sendrecv(&sol_size, 1, MPI_INT, ID_OTHER, 0,
                 &sol_other_size, 1, MPI_INT, ID_OTHER, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    sol_other = new double[sol_other_size];
    MPI_Sendrecv(sol, sol_size, MPI_DOUBLE, ID_OTHER, 1,
                 sol_other, sol_other_size, MPI_DOUBLE, ID_OTHER, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    runtime_self = get_runtime();
    runtime_other = 0;
    MPI_Sendrecv(&runtime_self, 1, MPI_DOUBLE, ID_OTHER, 2,
                 &runtime_other, 1, MPI_DOUBLE, ID_OTHER, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    iters = get_WF_iters();
	if(ID_SELF == 0){
        std::cout << std::setprecision(14) << ID_SELF << " " << iters << " " << runtime_self << " ";
	    for(int i = 0; i < sol_size; i++)
            std::cout << std::setprecision(14) << sol[i] << " ";
        std::cout << std::endl;

        std::cout << std::setprecision(14) << ID_OTHER << " " << iters << " " << runtime_other << " ";
	    for(int i = 0; i < sol_other_size; i++)
            std::cout << std::setprecision(14) << sol_other[i] << " ";
        std::cout << std::endl;
    }
}
