/*
Authors: Peter Meisrimel
March 2020
*/

#include "utils.h"
#include "input_reader_heat.h"
#include "problem_heat_python.h"
#include "mpi.h"
#include "Python.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION // disables from warning from outdated API
#define PY_ARRAY_UNIQUE_SYMBOL cool_ARRAY_API
#include "numpy/arrayobject.h"
#include "unistd.h"

int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);
    //wchar_t *program = Py_DecodeLocale(argv[0], NULL);
    //Py_SetProgramName(program);  /* optional but recommended */
    
    int np, ID_SELF, ID_OTHER;
    MPI_Comm_rank(MPI_COMM_WORLD, &ID_SELF);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    ID_OTHER = (ID_SELF + 1)%np;
    
    Py_Initialize();
    import_array(); // numpy stuff
    
    Problem * prob1, * prob2;
    // problem specific default parameters
    double t_end = 1;
    int timesteps = 100;
    
    // heat eq specific init parameters
    int gridsize = 32;
    double alpha = 1;
    double lambda = 0.01;
    int which_u0 = 0; // not used here
    
    int gridsize1 = gridsize, gridsize2 = gridsize;
    double alpha1 = alpha, alpha2 = alpha;
    double lambda1 = lambda, lambda2 = lambda;
    process_inputs_heat(argc, argv, alpha1, alpha2, lambda1, lambda2, gridsize1, gridsize2, which_u0);
    
    int which_conv = -2; // temperature at interface as cancelation criterion
    
    // IMPORTANT TO INITIALIZE PROBLEMS LIKE THIS
    // ensures that for parallel methods each processor initializes only one problem and for serial ones, both are initialized
    if (ID_SELF == 0){
        prob1 = new Problem_heat_python_D(gridsize1, alpha1, lambda1);
    }
    if (ID_OTHER == 0){
        prob2 = new Problem_heat_python_N(gridsize2, alpha2, lambda2);
    }
    setup_and_run_WFR(prob1, prob2, which_conv, t_end, timesteps, argc, argv);
    
    Py_Finalize();
    return 0; 
}
