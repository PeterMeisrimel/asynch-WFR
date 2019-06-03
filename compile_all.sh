cd src/toy
mkdir obj
make

cd ../bigger
mkdir obj
make

cd ../heat_DN
ffc -l dolfin heat.ufl
cmake .
make
