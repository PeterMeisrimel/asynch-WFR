/*
Authors: Peter Meisrimel
January 2019
*/

#ifndef UTILS_H_
#define UTILS_H_

#include "problem.h"

void process_inputs(int, char **, int&, double&, double&, int&, int&, int&, int&, bool&, bool&, bool&, int&);

void setup_and_run_WFR(Problem * prob1, Problem * prob2, int which_conv, double t_end, int timesteps, int argc, char *argv[]);

#endif // UTILS_H_
