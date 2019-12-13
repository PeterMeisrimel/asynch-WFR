/*
Authors: Peter Meisrimel
December 2018
*/

#ifndef WFR_JAC_H_
#define WFR_JAC_H_

#include "problem.h"
#include "WFR.h"

class WFR_JAC: public WFR_parallel{
private:
    void do_WF_iter(double WF_TOL, int WF_MAX_TOL, int steps_self, int steps_other);
public:
    WFR_JAC(int id_in_self, int id_in_other, double t_end, Problem * p, bool errlogging = false);

    void run(double WF_TOL, int WF_MAX_ITER, int steps_macro, int steps_self, int steps_other, int conv_check, int steps_converged_required_in, double relax_param, bool match_which_conv_relax);
};

#endif // WFR_JAC_H_
