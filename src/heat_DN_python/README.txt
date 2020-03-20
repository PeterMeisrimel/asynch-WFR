Current status:

It works, barely. Always appears to finish on error messages due to something in MPI_Finalize, which I could not figure out thus far.

1. The root cause of the error messages appears to be the dolfin expression that is used for the initial value. Possibly some error because it goes back to C in some way to evaluate this?
Replacing this by a UserExpression didn't work out so far

2. Despite the fact that Python is called from a single processor, FEniCs somehow registers that there are multiple processors and tries to make use of these. One can restrict the grid to MPI_COMM_SELF, i.e. to not have a distributed grid, and this will enable compuations and give reasonable results. However, this does NOT fix the underlying issue that Python is called with multiple processors in the first place. This might also be part of the reason why the initialization appears to be very slow (or it is just importing dolfin).

Either way, the proof of concept for calling Python FEniCs from C++ works and we will not try to fix this mess at this points.

Good luck if anyone tries to clean up this mess in the future.
