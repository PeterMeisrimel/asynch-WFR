/*
Authors: Peter Meisrimel
December 2018
*/

#ifndef U0_H_
#define U0_H_

#include "dolfin.h"
#include "math.h"
#include <iostream>

//using namespace dolfin;

// Initial condition
class InitialConditions : public dolfin::Expression{
private:
    double cons, bx, by;
	int which;
public:
    InitialConditions(double base_x, double base_y, int which_u0){
        cons = 500;
        bx   = base_x;
        by   = base_y;
		which = which_u0;
    }
    void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const{
		switch (which){
			case 0:{
				values[0]= cons*sin(M_PI*(bx + x[0])/2)*sin(M_PI*(by + x[1]));
				break;		
			}
			case 1:{
				values[0]= cons*sin(0.25*M_PI*(bx + x[0])*(bx + x[0]))*sin(M_PI*(by + x[1]));
				break;		
			}
			case 2:{
				values[0]= cons*sin(0.25*M_PI*(2 - bx - x[0])*(2 - bx - x[0]))*sin(M_PI*(by + x[1]));
				break;		
			}
			default:{
				std::cout << "invalid InitialCondition selected" << std::endl;
			}		
		}
    }
};

#endif //U0_H_
