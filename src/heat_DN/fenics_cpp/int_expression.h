/*
Authors: Peter Meisrimel, Benjamin Rueth
September 2018
*/

#ifndef INT_EXPRESSION_H_
#define INT_EXPRESSION_H_

// Expression class used to interpolate discrete data, used for exchange over interface
class InterpolatedExpression : public dolfin::Expression{
private:
    int _N;
	double _dx;
    double * _vals;
public:
	InterpolatedExpression(int N){
        _N = N;
        _dx = 1/float(_N - 1);
	}
    void set_vals(double * in){ 
        _vals = in;
    };
    void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const{
        int i;
        // find smallest i such that i*_dx > x[1] (yes, strict)
        for(i = 1; i < _N; i++)
            if(i*_dx > x[1])
                break;
        if(i == _N)
            i--;
        i--;
        double fac = (x[1] - i*_dx)/_dx;
        values[0] = (1 - fac) * _vals[i] + fac * _vals[i+1];
	}
};

#endif //INT_EXPRESSION_H_
