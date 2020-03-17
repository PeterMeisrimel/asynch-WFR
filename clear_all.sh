cd src/toy
make clean

cd ../bigger
make clean

cd ../heat_DN
rm heat.h
make clean
rm -r CMakeFiles
rm CMakeCache.txt

cd ../heat_DN_python
make clean
