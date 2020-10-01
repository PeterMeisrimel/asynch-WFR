cd src/toy
mkdir obj
make

cd ../bigger
mkdir obj
make

cd ../heat_DN
cd fenics_cpp
ffc -l dolfin heat.ufl
ffc -l dolfin heat_flux_grad.ufl
ffc -l dolfin heat_flux_weak.ufl
cd ..
cmake .
make
