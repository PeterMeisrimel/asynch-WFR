Authors: Peter Meisrimel (Lund University), Benjamin Rüth (Munich University, Germany), Philipp Birken (Lund University)
Contact: peter.meisrimel@na.lu.se
Published under the GNU General Public License v3.0 License

Requirement: at least mpi 2.x


This is a work in progress solver for 2 coupled systems of (time dependent) ordinary differenential equations using Waveform relaxation (WFR).
Currently supported Waveform relaxation schemes are Gauss-Seidel, Jacobi and a new, asynchronous method we are developing. Here, iterations can be done using windows of equidistant sizes and a varying number of total timesteps for each problem.

General assumed problem setting:
1. u' = f(u, a(v)),
2. v' = g(b(u), v),
where u, v are time-dependent functions. a() and b() are potential mappings to the interfaces, which would usually be of lower dimension

Overview:
The Waveform class is used to store and organize discrete Waveforms and is used internally in the various WFR schemes, interpolation is linear. It can be evaluated using the eval(t) function.

A given problems needs to provide a do_step(t, dt, uout, WF) function that performs a single time-integration step. Here, t and dt are time and stepsize. Assuming the first problem (u' = f(u, a(v))), uout needs to be b(u) and WF would be the discrete Waveform representation of a(v). Storage and management of internal variables must be handled by the problem classes themselves. 
The create_checkpoint and reset_to_checkpoint functions need to provide backups of the current solution and the possibility to resote to these.

Current examples are a two coupled scalar equations (problem_toy) and a slightly larger problem (problem_bigger), coupling a 2x2 system with a 3x3 system. These serve as proof of concept and more sophisticated examples are in the making.

More documentation will be done as well.
