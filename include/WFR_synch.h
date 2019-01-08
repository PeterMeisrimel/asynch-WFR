/*
Authors: Peter Meisrimel
December 2018
*/

#ifndef WFR_SYNCH_H_
#define WFR_SYNCH_H_

#include "problem.h"
#include "WFR.h"

class WFR_synch : public WFR{
protected:
  void integrate_window  (int); // steps_self
  bool check_convergence (double); // tol
public:
  void run               (double, int, int, int, int = 1); // tol, maxiter, macro_steps, steps_self, conv_check
};

#endif // WFR_SYNCH_H_
