cd src/toy
make clean

cd ../bigger
make clean

cd ../heat_DN
rm fenics_cpp/heat.h
rm fenics_cpp/heat_flux.h
make clean
rm -r CMakeFiles
rm CMakeCache.txt
