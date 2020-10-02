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
#include "input_reader_heat.h"
#include "mpi.h"
#include "dolfin.h"

#include "fenics_cpp/problem_heat_D.h"
#include "fenics_cpp/problem_heat_N.h"

#include "fenics_python/problem_heat_python.h"
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
    dolfin::init(argc, argv);
    dolfin::set_log_level(30); // fenics output level: "things that may go boom later"
    
    int np, ID_SELF, ID_OTHER;
    MPI_Comm_rank(MPI_COMM_WORLD, &ID_SELF);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    ID_OTHER = (ID_SELF + 1)%np;
    
    Problem * prob1, * prob2;
    // problem specific default parameters
    double t_end = 1;
    int timesteps = 100;
    
    // heat eq specific init parameters
    int gridsize = 32;
    double alpha = 1;
    double lambda = 0.01;
    int which_u0 = 0;
    // which solver to use for given problem, 1 = left, Dirichlet, 2 = right, Neumann
    // 0 = fenics cpp, 1 = fenics python, 2 = dune python
    int solver1 = 0, solver2 = 0;
    bool using_py = false; // flag if python is used, will cause call of various init/finalize functions
    
    int gridsize1 = gridsize, gridsize2 = gridsize;
    double alpha1 = alpha, alpha2 = alpha;
    double lambda1 = lambda, lambda2 = lambda;
    process_inputs_heat(argc, argv, alpha1, alpha2, lambda1, lambda2, gridsize1, gridsize2, which_u0, solver1, solver2);
    
    if (solver1 || solver2) // true if either of them not zero
        using_py = true;
    
    int which_conv = -2; // temperature at interface as cancelation criterion
    
    if (using_py){
        Py_Initialize();
        import_array(); // numpy stuff, ignore warning here
    }
    
    // IMPORTANT TO INITIALIZE PROBLEMS LIKE THIS
    // ensures that for parallel methods each processor initializes only one problem and for serial ones, both are initialized
    if (ID_SELF == 0){
        if (solver1 == 0)
            prob1 = new Problem_heat_D(gridsize1, alpha1, lambda1, which_u0);
        else if (solver1 == 1)
            prob1 = new Problem_heat_python_D(gridsize1, alpha1, lambda1);
        else if (solver1 == 2)
            prob1 = new Problem_heat_dune_D(gridsize1, alpha1, lambda1);
    }
    if (ID_OTHER == 0){
        if (solver2 == 0)
            prob2 = new Problem_heat_N(gridsize2, alpha2, lambda2, which_u0);
        else if (solver2 == 1)
            prob2 = new Problem_heat_python_N(gridsize2, alpha2, lambda2);
        else if (solver2 == 2)
            prob2 = new Problem_heat_dune_N(gridsize2, alpha2, lambda2);
    }
    setup_and_run_WFR(prob1, prob2, which_conv, t_end, timesteps, argc, argv); // ignore warning here
    
    if (using_py)
        Py_Finalize();
    
    return 0;
}