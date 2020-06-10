# Waveform Relaxation using asynchronous time-integration

This is an C++ implementation of a Waveform Relaxation based coupled problems solver (for exactly two problems as of now). Aside from the classical Jacobi and Gauss-Seidel Waveform Relaxation (WR) algorithms it also includes a new method using asynchronous communication in the time-integration, which is being researched and developed.

## Authors

Peter Meisrimel, Lund University, Sweden, peter.meisrimel@na.lu.se

Benjamin Rüth, Munich University, Germany

Robert Klöfkorn, NORCE, Norway

Philipp Birken, Lund University, Sweden

## License

Published under the GNU General Public License v3.0 License

## Software requirements

At least MPI 2.x, respectively an implementation of the MPI 3 Standard. Has been developed using Open MPI 2.1.1
For the heat equation test problem: FEniCs 2019.1.0
DUNE

## Publications

Meisrimel, Peter, and Philipp Birken. "Waveform Iteration with asynchronous communication." PAMM 19.1 (2019)

## Related Literature

Parallel multigrid waveform relaxation for parabolic problems, Stefan Vandevalle, as a good introduction to WR. 
A proceedings article on the new method using asynchronous communication is be published soon.

## Useage

A given problem class needs to be an instance of problem.h and provide the according functions.

Run compile_all.sh to compile all problems.

Run run_test_all.sh in runtime_stats (parameters for run_heat.py may need to be tweaked to a smaller resolution in space) to make a test run of all implemented problems.

The python scripts in runtime_stats provide a framework to run the problems for various input paramters settings and plots results.

More documentation will come soon.

## Test Problems

Included test problems:
toy - A coupling of two scalar equations which was for development and verification.

bigger - A coupling of a scalar equation with 2 other scalar equations, also used for verification and testing.

heat_DN - Coupling of two heat problems via a Dirichlet-Neumann coupling over the interface. Uses a FE discretization using FEniCs.

## Acknowledgments

Joachim Hein, Lund University, Sweden, for support on many technical questions regarding MPI.

Azahar Monge, DeustoTech, Spain, usage of older code to verify correctness of heat equation example.
