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

/*
need total of 6 relaxation inputs for new method?
*/

void process_inputs(int argc, char **argv,
                    double& t_end,
                    int& timesteps1, int& timesteps2,
                    int& runmode, int& macrosteps, int& maxiter, double& WF_TOL, bool& error_logging, int& nsteps_conv_check,
                    double &theta_relax1, double &theta_relax2,
                    bool &var_relax,
                    double & theta_relax_gs_a_1, double & theta_relax_gs_a_2,
                    double & theta_relax_gs_b_1, double & theta_relax_gs_b_2){
    if (argc % 2 != 1){
        throw std::invalid_argument("invalid number of input arguments, needs to be even, check for accidental spaces");
    }

    for(int i = 1; i < argc; i+=2){
        std::string arg = argv[i];
//        problem parameters
        if(arg == "-tend")
            t_end = atof(argv[i+1]);
//        time-integration
        else if (arg == "-timesteps1")
            timesteps1 = atoi(argv[i+1]);
        else if (arg == "-timesteps2")
            timesteps2 = atoi(argv[i+1]);
        else if (arg == "-timesteps"){
            timesteps1 = atoi(argv[i+1]);
            timesteps2 = atoi(argv[i+1]);
        }
//        WR specific settings
        else if (arg == "-runmode"){
            std::string arg2 = argv[i+1];
            if ((arg2 == "GS") || (arg2 == "GS1")) // GS in 1 -> 2 ordering
                runmode = 1;
            else if (arg2 == "GS2") // GS in 2 -> 1 ordering
                runmode = 2;
            else if (arg2 == "JAC")
                runmode = 3;
            else if (arg2 == "NEW")
                runmode = 4;
            else
                throw std::invalid_argument("invalid runmode, check your spelling");
        }
        else if (arg == "-macrosteps")
            macrosteps = atoi(argv[i+1]);
        else if (arg == "-maxiter") // max number of iterations
            maxiter = atoi(argv[i+1]);
        else if (arg == "-wftol") // tolerance
            WF_TOL = atof(argv[i+1]);
        else if (arg == "-errlog") // logging of errors for error over iter plots
            error_logging = bool(atoi(argv[i+1]));            
        else if (arg == "-nsteps_conv_check") // termination check
            nsteps_conv_check = atoi(argv[i+1]);
//        relaxation, fixed splittings
        else if (arg == "-theta_relax1") // output of prob 1
            theta_relax1 = atof(argv[i+1]);
        else if (arg == "-theta_relax2") // output of prob 2
            theta_relax2 = atof(argv[i+1]);
        else if (arg == "-theta_relax"){ // both using the same
            theta_relax1 = atof(argv[i+1]);
            theta_relax2 = atof(argv[i+1]);
        }
//        relaxation, variable splittings
        else if (arg == "-var_relax")
            var_relax = bool(atoi(argv[i+1]));
//        take previous parameters as jacobi relaxation parameters
//        gs_a_n  a: 1 -> 2 ordering, n = 1 first problem, n = 2 second probem
        else if (arg == "-theta_relax_gs_a_1")
            theta_relax_gs_a_1 = atof(argv[i+1]);
        else if (arg == "-theta_relax_gs_a_2")
            theta_relax_gs_a_2 = atof(argv[i+1]);
        else if (arg == "-theta_relax_gs_a"){ // same for both
            theta_relax_gs_a_1 = atof(argv[i+1]);
            theta_relax_gs_a_2 = atof(argv[i+1]);
        }
//        gs_b_n  b: 2 -> 1 ordering, n = 1 first problem, n = 2 second probem
        else if (arg == "-theta_relax_gs_b_1")
            theta_relax_gs_b_1 = atof(argv[i+1]);
        else if (arg == "-theta_relax_gs_b_2")
            theta_relax_gs_b_2 = atof(argv[i+1]);
        else if (arg == "-theta_relax_gs_b"){ // same for both
            theta_relax_gs_b_1 = atof(argv[i+1]);
            theta_relax_gs_b_2 = atof(argv[i+1]);
        }
    }
}

void setup_and_run_WFR(Problem * prob1, Problem * prob2, int which_conv, double t_end, int timesteps, int argc, char *argv[]){
    int np, ID_SELF, ID_OTHER;
    MPI_Comm_rank(MPI_COMM_WORLD, &ID_SELF);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    ID_OTHER = (ID_SELF + 1)%np;
    
    // default init parameters
    int runmode = 1; // GS, JACOBI or NEW, see below for more details
    bool errlogging = false;
    int nsteps_conv_check = 3;
    
    int timesteps1 = timesteps, timesteps2 = timesteps;
    
    // default running parameters
    double WF_TOL = 1e-6;
    int WR_MAXITER = 50;
    int num_macro = 1;
    
    double theta_relax1 = 1, theta_relax2 = 1; // basic ones for fixed splittings, jacobi for variable
    
    bool var_relax = false; // new scheme only with relaxation parameters
    double theta_relax_a_1 = 1, theta_relax_a_2 = 1; // gs 1->2 ordering
    double theta_relax_b_1 = 1, theta_relax_b_2 = 1; // gs 2->1 ordering
    
    process_inputs(argc, argv,
                   t_end,
                   timesteps1, timesteps2,
                   runmode, num_macro, WR_MAXITER, WF_TOL, errlogging, nsteps_conv_check,
                   theta_relax1, theta_relax2, 
                   var_relax, 
                   theta_relax_a_1, theta_relax_a_2,
                   theta_relax_b_1, theta_relax_b_2);
    
    // parallel methods only get one problem as input, pick correct one if there are multiple processes
    Problem * prob;
    if (ID_SELF == 0)
        prob = prob1;
    else
        prob = prob2;

    // build WFR method based on input parameters
    WFR * wfr_method;
    switch(runmode){
        case 1:{ // Gauss-Seidel (GS), 1 -> 2 ordering
            if (np != 1)
                MPI_Abort(MPI_COMM_WORLD, 1);
            wfr_method = new WFR_GS(t_end, prob1, prob2, true);
            
            wfr_method -> set_relax_params(theta_relax_a_1, theta_relax_a_2);
            break;
        }
        case 2:{ // Gauss-Seidel (GS), 2 -> 1 ordering
            if (np != 1)
                MPI_Abort(MPI_COMM_WORLD, 1);
            wfr_method = new WFR_GS(t_end, prob1, prob2, false);
            
            wfr_method -> set_relax_params(theta_relax_b_1, theta_relax_b_2);
            break;
        }
        case 3:{ // Jacobi (JAC)
            if (np != 2)
                MPI_Abort(MPI_COMM_WORLD, 1);
            wfr_method = new WFR_JAC(ID_SELF, ID_OTHER, t_end, prob);
            
            if (ID_SELF == 0)
                wfr_method -> set_relax_params(theta_relax1);
            else
                wfr_method -> set_relax_params(theta_relax2);
            break;
        }
        case 4:{ // NEW method, still needs better name
            if (np != 2)
                MPI_Abort(MPI_COMM_WORLD, 1);
            // variable relaxation
            if (var_relax){
                wfr_method = new WFR_NEW_var_relax(ID_SELF, ID_OTHER, t_end, prob);
                // relaxation parameters reversed here, since relaxation is done on the receiving end
                if (ID_SELF == 0)
                    wfr_method -> set_relax_params(theta_relax2, theta_relax_a_2, theta_relax_b_2);
                else
                    wfr_method -> set_relax_params(theta_relax1, theta_relax_a_1, theta_relax_b_1);
            }else{ // fixed relaxation
                wfr_method = new WFR_NEW(ID_SELF, ID_OTHER, t_end, prob);
                if (ID_SELF == 0)
                    wfr_method -> set_relax_params(theta_relax1);
                else
                    wfr_method -> set_relax_params(theta_relax2);
            }
            break;
        }
        default:{
            throw std::invalid_argument("invalid runmode");
            break;
        }
    }
    
    wfr_method -> run(WF_TOL, WR_MAXITER, num_macro, timesteps1, timesteps2, which_conv, nsteps_conv_check, errlogging);

    if (not errlogging){
        wfr_method -> write_results();
    }else{
        MPI_Barrier(MPI_COMM_WORLD);
        wfr_method -> write_error_log();
    }
}
