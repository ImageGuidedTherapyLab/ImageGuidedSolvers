"""
iterative reconstruction
"""

# import needed modules
import petsc4py, numpy, sys
petsc4py.init(sys.argv)

from petsc4py import PETSc
from mpi4py import MPI

# break processors into separate communicators
petscRank = PETSc.COMM_WORLD.getRank()
petscSize = PETSc.COMM_WORLD.Get_size()
sys.stdout.write("petsc rank %d petsc nproc %d\n" % (petscRank, petscSize))

# break up processors into communicators
NumProcsPerSubComm = 4
color = petscRank/NumProcsPerSubComm 
NumSubCommunicators = petscSize / NumProcsPerSubComm + 1
subcomm = PETSc.Comm(MPI.COMM_WORLD.Split(color))
subRank = subcomm.Get_rank()
subSize = subcomm.Get_size()
sys.stdout.write("number of sub communictors %d \n" % ( NumSubCommunicators ))
sys.stdout.write("subcomm rank %d subcomm nproc %d\n" % (subRank, subSize))

# get image data and store as petsc vec
NumCoils = 10 
#import vtk 
##import vtk.util.numpy_support as vtkNumPy 
#vtkReader = vtk.vtkXMLImageDataReader() 
#vtkReader.SetFileName('tmapstd.0060.vti') 
#vtkReader.Update() 
#dimensions = vtkReader.GetOutput().GetDimensions()
dimensions = (10,10,10)
#print "dimensions ",dimensions 
#numberPointsImage =  vtkReader.GetOutput().GetNumberOfPoints()
numberPointsImage = 3000
#print "#points",numberPointsImage 
##image_cells = vtkReader.GetOutput().GetCellData() 
##data_array = vtkNumPy.vtk_to_numpy(image_cells.GetArray('phi')) 
data_array = numpy.zeros(numberPointsImage, dtype=PETSc.ScalarType)
v1 = PETSc.Vec().createWithArray(data_array, comm=PETSc.COMM_SELF)

# setup linear system matrix
A = PETSc.Mat()
A.create(subcomm)
#A.create(PETSc.COMM_WORLD)
A.setSizes([numberPointsImage, numberPointsImage])
#MATPYTHON should be a python version of MATSHELL
A.setType('python')
# set shell context
import signalmodel
shell = signalmodel.PySignalModel(4.0) # shell context
shell.SetuplibMesh( subcomm ) # initialize libMesh data structures

# initialize libMesh data structures
shell.SetupStructuredGrid( (10,10,4) ,[0.0,1.0],[0.0,1.0],[0.0,1.0]) 

# setup imaging to interpolate onto FEM mesh
shell.SetupImaging( dimensions ,[0.0,0.0,0.0],[1.0,1.0,1.0]) 

# add systems to hold data and initialize
shell.AddExplicitSystem( "Spin"               ,0 ) 
shell.AddExplicitSystem( "RealCoilSensitivity",0 ) 
shell.AddExplicitSystem( "ImagCoilSensitivity",0 ) 
shell.AddExplicitSystem( "T2Star"             ,0 ) 
shell.InitializeEquationSystems( ) 
shell.SetupLinearGradientModel( [3.0,4.0,9.0], 63.0) 

# print info
shell.printSelf() 

# Project imaging onto libMesh data structures
shell.ProjectImagingToFEMMesh("RealCoilSensitivity", v1)  
shell.ProjectImagingToFEMMesh("ImagCoilSensitivity", v1)  
shell.ProjectImagingToFEMMesh("T2Star"          , v1)  

A.setPythonContext(shell)
print A.getLocalSize() 

timeValues = [2.0,8.0,4.0]
print shell.AssembleSignal(timeValues)  

shell.WriteEquationSystems("fem_data.e")  

## setup linear system vectors
#x, b = A.getVecs()
#x.set(0.0)
#b.set(1.0)
#
## setup Krylov solver
#ksp = PETSc.KSP().create()
#pc = ksp.getPC()
#ksp.setType('cg')
#pc.setType('none')
#
## iteratively solve linear
## system of equations A*x=b
#ksp.setOperators(A)
#ksp.setFromOptions()
#ksp.solve(b, x)
#
## scale solution vector to
## account for grid spacing
#x.scale(h**2)
#

#try:
#    from matplotlib import pylab
#except ImportError:
#    raise SystemExit("matplotlib not available")
#from numpy import mgrid
#X, Y =  mgrid[0:1:1j*n,0:1:1j*n]
#Z = x[...].reshape(n,n,n)[:,:,n/2-2]
#pylab.contourf(X, Y, Z)
#pylab.axis('equal')
#pylab.colorbar()
#pylab.show()
#
