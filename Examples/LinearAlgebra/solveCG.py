import sys, petsc4py
print "init"
PetscOptions =  sys.argv
PetscOptions.append("-ksp_monitor")
PetscOptions.append("-ksp_converged_reason")
PetscOptions.append("-ksp_rtol")
PetscOptions.append("1.e-20")
petsc4py.init(sys.argv)

print "PETSC"
from petsc4py import PETSc

print "numpy"
import numpy


n = 1000
print "create dense"
unsymmetric = numpy.random.rand(n,n)
symmetric = unsymmetric + unsymmetric.transpose()
A = PETSc.Mat().createDense([n, n], array=symmetric, comm=PETSc.COMM_WORLD)
A.assemblyBegin();A.assemblyEnd()
A.shift(n/39.0)
#A.view()
A.assemblyBegin();A.assemblyEnd()

# create linear solver
ksp = PETSc.KSP()
ksp.create(PETSc.COMM_WORLD)
# use conjugate gradients
ksp.setType('cg')
# and incomplete Cholesky
ksp.getPC().setType('none')
#ksp.getPC().setType('jacobi')
# obtain sol & rhs vectors
x, b = A.getVecs()
x.set(0)
b.set(1)
# and next solve
ksp.setOperators(A)
ksp.setFromOptions()
ksp.view()
ksp.solve(b, x)
