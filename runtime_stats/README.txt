Author: Peter Meisrimel
January 2019

Usage of runtime measurement scripts:

Create measurements by executing run_tolerances_xyz.py, where you can adapt the parameters accordingly in the beginning. This does:
1. Runs the given runmodes for the given tolerances <times> number of times
2. Creates output files to record results (number of iterations, runtime, solution)
3. Processes said output files, taking averages of runtime, iterations and solution
4. Delete previous output files
5. Store processes results in a results_<name>_<timestamp>.txt file

Important notes:
The processing of data is based on the output format of the main_files, which need to be in the format <processor id> <num iterations> <runtime> <solution>, see the toy problem as example.
Solution means the final values of the data that the problem communicate to each other.
Any runmode that is not GS (Gauss-Seidel) will omit the smallest tolerance. The rationale is that one will use GS as a reference solution and thus needs it for a smaller tolerance.
As the solution is at least a vector of size 2, it here becomes a of all entries (sorted, first vector, second vector).

For plotting, simply run: python3 plotting.py <resultsfile.txt>, to get plots in both .png and .eps
