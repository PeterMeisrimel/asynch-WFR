/*
Authors: Peter Meisrimel
January 2019
*/

#ifndef UTILS_H_
#define UTILS_H_

#include "problem.h"

void process_inputs(int argc, char **argv,
                    double& t_end,
                    int& timesteps1, int& timesteps2,
                    int& runmode, int& macrosteps, int& maxiter, double& WF_TOL, bool& error_logging, int& nsteps_conv_check,
                    double &theta_relax1, double &theta_relax2,
                    bool &var_relax,
                    double & theta_relax_gs_a_1, double & theta_relax_gs_a_2,
                    double & theta_relax_gs_b_1, double & theta_relax_gs_b_2);

void setup_and_run_WFR(Problem * prob1, Problem * prob2, int which_conv, double t_end, int timesteps, int argc, char *argv[]);

#endif // UTILS_H_
