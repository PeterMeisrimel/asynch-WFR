#ifndef U0_H_
#define U0_H_

#include "dolfin.h"
#include "math.h"

using namespace dolfin;

// Initial condition
class InitialConditions : public Expression{
private:
  double cons, bx, by;
public:
  InitialConditions(double c, double base_x, double base_y){
    cons = c;
    bx   = base_x;
    by   = base_y;
  }
  void eval(Array<double>& values, const Array<double>& x) const{
    values[0]= cons*sin(M_PI*(bx + x[0])/2)*sin(M_PI*(by + x[1]));
  }
};

#endif //U0_H_
