/*
Authors: Peter Meisrimel
December 2018
*/

#include "Python.h"
/*
https://docs.python.org/3/extending/extending.html: 
"Since Python may define some pre-processor definitions which affect the standard headers on some systems, 
you must include Python.h before any standard headers are included"
*/

#include "utils.h"
#include "mpi.h"

#include "dune_python/problem_heat_dune.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION // disables from warning from outdated API
#define PY_ARRAY_UNIQUE_SYMBOL cool_ARRAY_API
#include "numpy/arrayobject.h"
#include "unistd.h"

/*
Problem with current setup: One requires installations of fenics & (dune*) to actually be able to run this code 
*: fenics, resp. dolfin c++ is needed while dolfin and dune in python are only needed if called
*/

int main(int argc, char *argv[]){
//    dolfin::init(argc, argv);
//    dolfin::set_log_level(30); // fenics output level: "things that may go boom later"
    MPI_Init(&argc, &argv);
    
    int np, ID_SELF, ID_OTHER;
    MPI_Comm_rank(MPI_COMM_WORLD, &ID_SELF);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    ID_OTHER = (ID_SELF + 1)%np;
    
    Problem * prob1, * prob2;
    // problem specific default parameters
    double t_end = 0.001;
    int timesteps = 100;
    
    // heat eq specific init parameters
    int gridsize = 32;
    double alpha = 7836*443;
    double lambda = 48.9;
    
    int which_conv = -2; // temperature at interface as cancelation criterion
    
    Py_Initialize();
    import_array(); // numpy stuff, ignore warning here
    
    // IMPORTANT TO INITIALIZE PROBLEMS LIKE THIS
    // ensures that for parallel methods each processor initializes only one problem and for serial ones, both are initialized
    if (ID_SELF == 0)
        prob1 = new Problem_heat_dune_euler_fluid(gridsize);
    if (ID_OTHER == 0)
        prob2 = new Problem_heat_dune_N_euler(gridsize, alpha, lambda);
        
    setup_and_run_WFR(prob1, prob2, which_conv, t_end, timesteps, argc, argv); // ignore warning here
    
    Py_Finalize();
    
    return 0;
}