SETUP: 

Open problem_heat_python.cpp
Change the line "PyRun_SimpleString("sys.path.append('/home/peter/asynch-WFR/src/heat_DN_python/')");" to correctly represent wherever you saved this




Current status: It works, somewhat.

Always appears to finish on error messages due to something in MPI_Finalize, which I could not figure out thus far. See also minimal working example (folder: min_ex) that can reproduce this problem. The root cause appears to be the dolfin expression. Replacing it by a UserExpression does not solve the problem. Possibly, the issue is that it goes back to C in some way to evaluate this? 

See also: https://fenicsproject.discourse.group/t/expression-class-causing-trouble-for-python-embedded-in-c/2851

This can be "fixed" skipping the MPI_Finalize and running with -quiet to hide the resulting error message at the end. However, this will also hide any legit MPI error messages.

Maybe this stuff will be fixed in future FEniCs versions?

This mostly serves as proof of concept.
