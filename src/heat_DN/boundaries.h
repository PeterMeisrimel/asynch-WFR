/*
Authors: Peter Meisrimel, Benjamin Rueth
September 2018
*/

#ifndef BOUNDARIES_H_
#define BOUNDARIES_H_

// left, top and bottom
class Boundaries_D : public SubDomain{
    bool inside(const Array<double>& x, bool on_boundary) const{
        return x[0] < DOLFIN_EPS or x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS;
    }
};

// right boundary
class Interface_D : public SubDomain{
    bool inside(const Array<double>& x, bool on_boundary) const{
        return x[0] > 1.0 - DOLFIN_EPS;
    }
};

// Dirichlet boundary on top, bottom and left side
class Boundaries_N : public SubDomain{
    bool inside(const Array<double>& x, bool on_boundary) const{
        return x[0] > 1 - DOLFIN_EPS or x[1] < DOLFIN_EPS     or x[1] > 1.0 - DOLFIN_EPS;
    }
};

#endif // BOUNDARIES_H_
