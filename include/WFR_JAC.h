#ifndef WFR_JAC_H_
#define WFR_JAC_H_

#include "problem.h"
#include "WFR_synch.h"

class WFR_JAC: public WFR_synch{
private:
  void do_WF_iter        (double, int, int, int); // tol, maxiter, steps_self, steps_other
public:
	WFR_JAC(int, int, double, Problem *);// ID_SELF, ID_OTHER, t_end, Problem
};

#endif // WFR_JAC_H_
