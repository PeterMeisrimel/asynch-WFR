/*
Authors: Peter Meisrimel
December 2018
*/

#ifndef WFR_GS_H_
#define WFR_GS_H_

#include "problem.h"
#include "WFR_synch.h"

class WFR_GS: public WFR_synch{
private:
  bool FIRST;

  void do_WF_iter        (double, int, int, int); // tol, maxiter, steps_self, steps_other
public:
	WFR_GS(int, int, double, Problem *, bool);// ID_SELF, ID_OTHER, t_end, Problem, FIRST
};

#endif // WFR_GS_H_
