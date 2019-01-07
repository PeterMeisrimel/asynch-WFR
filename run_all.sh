# simple functionality test of all current examples
mpirun -np 2 src/toy/GS
mpirun -np 2 src/toy/JAC
mpirun -np 2 src/toy/NEW

mpirun -np 2 src/bigger/GS
mpirun -np 2 src/bigger/JAC
mpirun -np 2 src/bigger/NEW

mpirun -np 2 src/heat_DN/GS
mpirun -np 2 src/heat_DN/JAC
mpirun -np 2 src/heat_DN/NEW
