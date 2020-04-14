import dolfin as dol

mesh = dol.UnitSquareMesh(dol.MPI.comm_self, 10, 10)
V = dol.FunctionSpace(mesh, "CG", 1)
u0 = dol.interpolate(dol.Expression("1", degree = 1, mpi_comm = dol.MPI.comm_self), V)
print("interpolation successful", u0)
