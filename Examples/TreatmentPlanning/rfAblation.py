"""
RF treatment planning model 
"""
# import needed modules
import petsc4py, numpy, sys
PetscOptions =  sys.argv
PetscOptions.append("-ksp_monitor")
PetscOptions.append("-ksp_converged_reason")
PetscOptions.append("-ksp_rtol")
PetscOptions.append("1.0e-12")
#PetscOptions.append("-help")
petsc4py.init(PetscOptions)

from petsc4py import PETSc
from mpi4py import MPI

# break processors into separate communicators
petscRank = PETSc.COMM_WORLD.getRank()
petscSize = PETSc.COMM_WORLD.Get_size()
sys.stdout.write("petsc rank %d petsc nproc %d\n" % (petscRank, petscSize))

# break up processors into communicators
NumProcsPerSubComm = 1000000
color = petscRank/NumProcsPerSubComm 
NumSubCommunicators = petscSize / NumProcsPerSubComm + 1
subcomm = PETSc.Comm(MPI.COMM_WORLD.Split(color))
subRank = subcomm.Get_rank()
subSize = subcomm.Get_size()
sys.stdout.write("number of sub communictors %d \n" % ( NumSubCommunicators ))
sys.stdout.write("subcomm rank %d subcomm nproc %d\n" % (subRank, subSize))

# set shell context
import femLibrary
fem = femLibrary.PyFEMInterface(4.0) 
fem.SetuplibMesh( subcomm ) # initialize libMesh data structures

# inifile = None ==> setup from command line
inifile=None
fem.SetupIni( inifile ) 
# set any additional parameters needed
fem.SetIniValue( "electric_conductivity/s_0_healthy","7.0" ) 
fem.SetIniValue( "perfusion/w_0_healthy","6.0" ) 
fem.SetIniValue( "perfusion/w_0_tumor"  ,"0.0" ) 
fem.SetIniValue( "bc/u_dirichletid","2" )    #apply dirichlet data on last domain
fem.SetIniValue( "bc/volt_dirichletid","2" ) #apply dirichlet data on last domain
fem.SetIniValue( "bc/v_newton_coeff","10000.0" ) #apply dirichlet data on last domain

# initialize FEM Mesh
# must setup Ini File first
# Identity transform
RotationMatrix = [[1.,0.,0.],
                  [0.,1.,0.],
                  [0.,0.,1.]]
Translation =     [0.,0.,0.]
fem.SetupUnStructuredGrid( "clusterVessel.e",0,RotationMatrix, Translation  ) 
#fem.SetupStructuredGrid( (10,10,4) ,[0.0,1.0],[0.0,1.0],[0.0,1.0]) 

# add the data structures for the Background System Solve
# set deltat, number of time steps, power profile, and add system
deltat = 0.5
ntime  = 40 
RFSystem=0
fem.AddPennesSystem(RFSystem,deltat,ntime,[[2,28,46,78,119],[0.0,15.0,0.0,0.0,0.0]]) 

# initialize libMesh data structures
fem.InitializeEquationSystems( ) 

# print info
fem.printSelf() 

# write IC
fem.WriteTimeStep("fem_data.e" , 1, 0.0 )  

# loop over time steps and solve
for timeID in range(1,ntime):
   print "time step = " ,timeID
   fem.UpdateTransientSystemTimeStep("StateSystem",timeID ) 
   fem.SystemSolve( "StateSystem" ) 
   #fem.StoreTransientSystemTimeStep("StateSystem",timeID ) 
   fem.WriteTimeStep("fem_data.e" , timeID+1, timeID*deltat )  
