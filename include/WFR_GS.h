/*
Authors: Peter Meisrimel
December 2018
*/

#ifndef WFR_GS_H_
#define WFR_GS_H_

#include "problem.h"
#include "WFR.h"

class WFR_GS: public WFR_serial{
private:
    void do_WF_iter(double WF_TOL, int WF_MAX_ITER, int steps_self, int steps_other);
public:
    WFR_GS(double t_end, Problem * p1, Problem * p2, bool first, bool errlogging = false);

    void run(double WF_TOL, int WF_MAX_ITER, int steps_macro, int steps_self, int steps_other, int conv_check, int steps_converged_required_in, double relax_param, bool match_which_conv_relax);
};

#endif // WFR_GS_H_
