"""
   default LU Facto
"""
import sys, petsc4py
PetscOptions =  sys.argv
PetscOptions.append("-log_summary")
petsc4py.init(PetscOptions)

from petsc4py import PETSc

import numpy
import scipy.io as scipyio

# create stages for logging
LogSetup = PETSc.Log.Stage("Setup")
LogFactor= PETSc.Log.Stage("Factor")
LogSolve = PETSc.Log.Stage("Solve")

LogSetup.push()
# read in matrix
#sumcov = scipyio.loadmat("/home/fuentes/DDDAS/trunk/dddas/SumCov.mat" )
#n = len(sumcov['Mat_0'])
#AMat = PETSc.Mat().createDense([n, n], array=sumcov['Mat_0'],comm=PETSc.COMM_WORLD)
n = 7000
AMat = PETSc.Mat().createDense([n, n], array=numpy.random.rand(n,n), comm=PETSc.COMM_WORLD)
AMat.assemblyBegin();AMat.assemblyEnd()

# storage space for solution
AInv = PETSc.Mat().createDense([n, n], array=numpy.zeros((n,n)), comm=PETSc.COMM_WORLD)
AInv.assemblyBegin();AInv.assemblyEnd()

# setup identity matrix
eye  = PETSc.Mat().createDense([n, n], array=numpy.zeros((n,n)), comm=PETSc.COMM_WORLD)
#eye  = PETSc.Mat().createAIJ([n, n], comm=PETSc.COMM_WORLD)
eye.assemblyBegin();eye.assemblyEnd()
eye.shift(1.0)
eye.assemblyBegin();eye.assemblyEnd()
LogSetup.pop()

#Factor
LogFactor.push()
(isrow,iscol) = AMat.getOrdering("nd")
AMat.factorLU(isrow,iscol)
LogFactor.pop()

# Solve
LogSolve.push()
AMat.matSolve(eye,AInv)
LogSolve.pop()
    
print "n = %d"%n
print AInv
#AMat.view()
print AInv.norm(0)
print AInv.norm(2)
print AInv.norm(3)
