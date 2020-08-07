/*
Authors: Peter Meisrimel
January 2019
*/

#ifndef UTILS_H_
#define UTILS_H_

#include "problem.h"

void process_inputs(int argc, char **argv, int& runmode, double& WF_TOL, double& t_end,
                    int& timesteps1, int& timesteps2, int& macrosteps, int& maxiter,
                    bool& FIRST, bool& error_logging, int& nsteps_conv_check,
                    double &w_relax, bool &new_relax_opt, double &w_relax_jac,
                    double &w_relax_gs1, double &w_relax_gs2);

void setup_and_run_WFR(Problem * prob1, Problem * prob2, int which_conv, double t_end, int timesteps, int argc, char *argv[]);

#endif // UTILS_H_
