cd src/toy
mkdir obj
make

cd ../bigger
mkdir obj
make

cd ../heat_DN
cd fenics_cpp
ffc -l dolfin heat.ufl
ffc -l dolfin heat_flux.ufl
cd ..
cmake .
make
