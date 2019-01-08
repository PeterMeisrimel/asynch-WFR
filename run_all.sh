# simple functionality test of all current examples
mpirun -np 2 src/toy/TOY -runmode GS
mpirun -np 2 src/toy/TOY -runmode JAC
mpirun -np 2 src/toy/TOY -runmode NEW

mpirun -np 2 src/bigger/BIGGER -runmode GS
mpirun -np 2 src/bigger/BIGGER -runmode JAC
mpirun -np 2 src/bigger/BIGGER -runmode NEW

mpirun -np 2 src/heat_DN/heat_DN -runmode GS
mpirun -np 2 src/heat_DN/heat_DN -runmode JAC
mpirun -np 2 src/heat_DN/heat_DN -runmode NEW
