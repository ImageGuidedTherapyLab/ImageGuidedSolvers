import sys, petsc4py
PetscOptions =  sys.argv
PetscOptions.append("-log_summary")
#PetscOptions.append("-info")
petsc4py.init(PetscOptions)

from petsc4py import PETSc

import numpy
import scipy.io as scipyio

# break processors into separate communicators
petscRank = PETSc.COMM_WORLD.getRank()
petscSize = PETSc.COMM_WORLD.Get_size()
sys.stdout.write("petsc rank %d petsc nproc %d\n" % (petscRank, petscSize))

# create stages for logging
LogSetup      = PETSc.Log.Stage("Setup")
LogMUMPSFactor= PETSc.Log.Stage("MUMPSFactor")
LogMUMPSSolve = PETSc.Log.Stage("MUMPSSolve")
LogSUPERLUFactor= PETSc.Log.Stage("SUPERLUFactor")
LogSUPERLUSolve = PETSc.Log.Stage("SUPERLUSolve")
LogDFLTFactor = PETSc.Log.Stage("DFLTFactor")
LogDFLTSolve  = PETSc.Log.Stage("DFLTSolve")

LogSetup.push()
#read petsc matrix
matviewer = PETSc.Viewer().createBinary("/home/fuentes/DDDAS/trunk/dddas/SystemMatrix.bin",mode=PETSc.Viewer.Mode.READ,comm=PETSc.COMM_SELF) 
AMat = PETSc.Mat().create(comm = PETSc.COMM_SELF) 
AMat.load(matviewer,PETSc.Mat.Type.SEQAIJ) 
(M,N) = AMat.getSize()

# storage space for solution
nloc = N / petscSize
AInv = PETSc.Mat().createDense([M, nloc], array=numpy.zeros((M,nloc)), comm=PETSc.COMM_SELF)
AInv.assemblyBegin();AInv.assemblyEnd()

# setup identity matrix
RHS  = PETSc.Mat().createDense([M, nloc], array=numpy.random.rand(M,nloc), comm=PETSc.COMM_SELF)
LogSetup.pop()

print "nloc = %d"%nloc
import femLibrary
# Mumps Factor
LogMUMPSFactor.push()
AMatFactored = femLibrary.GetFactorMat(AMat,"mumps","nd")
LogMUMPSFactor.pop()

# Mumps Solve
LogMUMPSSolve.push()
AMatFactored.matSolve(RHS,AInv)
LogMUMPSSolve.pop()

# superlu Factor
LogSUPERLUFactor.push()
AMatSuperLUFactored = femLibrary.GetFactorMat(AMat,"superlu_dist","nd")
LogSUPERLUFactor.pop()

# Mumps Solve
LogSUPERLUSolve.push()
AMatSuperLUFactored.matSolve(RHS,AInv)
LogSUPERLUSolve.pop()

# Default Factor
LogDFLTFactor.push()
(isrow,iscol) = AMat.getOrdering("nd")
AMat.factorLU(isrow,iscol)
LogDFLTFactor.pop()

# Default Factor
LogDFLTSolve.push()
AMat.matSolve(RHS,AInv)
LogDFLTSolve.pop()
    
print AInv
#AMat.view()
print AInv.norm(0)
print AInv.norm(2)
print AInv.norm(3)
