"""
treatment planning model 
"""
# import needed modules
import petsc4py, numpy, sys
PetscOptions =  sys.argv
PetscOptions.append("-ksp_monitor")
PetscOptions.append("-ksp_rtol")
PetscOptions.append("1.0e-15")
#PetscOptions.append("-help")
petsc4py.init(PetscOptions)

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

# set shell context
# TODO import vtk should be called after femLibrary ???? 
# FIXME WHY IS THIS????
import femLibrary
fem = femLibrary.PyFEMInterface(4.0) 
fem.SetuplibMesh( subcomm ) # initialize libMesh data structures

# inifile = None ==> setup from command line
inifile=None
fem.SetupIni( inifile ) 
# from Duck table 2.15
fem.SetIniValue( "material/specific_heat","3840.0" ) 
# set ambient temperature 
fem.SetIniValue( "initial_condition/u_init","21.0" ) 
# from Duck
fem.SetIniValue( "thermal_conductivity/k_0_healthy","0.6" ) 
# from Welch
fem.SetIniValue( "optical/mu_a_healthy","3.0" ) 
# FIXME large mu_s (> 30) in agar causing negative fluence to satisfy BC 
fem.SetIniValue( "optical/mu_s_healthy",".3" ) 
# from AE paper
#http://scitation.aip.org/journals/doc/MPHYA6-ft/vol_36/iss_4/1351_1.html#F3
fem.SetIniValue( "optical/guass_radius","0.0025" ) 
fem.SetIniValue( "optical/mu_a_tumor","51.0" ) 
fem.SetIniValue( "optical/mu_s_tumor","192.0" ) 
#fem.SetIniValue( "optical/mu_a_tumor","71.0" ) 
#fem.SetIniValue( "optical/mu_s_tumor","89.0" ) 
fem.SetIniValue( "optical/anfact","0.9" ) 
#
#  given the original orientation as two points along the centerline z = x2 -x1
#     the transformed orienteation would be \hat{z} = A x2 + b - A x1 - b = A z
#  ie transformation w/o translation which is exactly w/ vtk has implemented w/ TransformVector
#  TransformVector = TransformPoint - the transation
#Setup Affine Transformation for registration
RotationMatrix = [[1.,0.,0.],
                  [0.,1.,0.],
                  [0.,0.,1.]]
Translation =     [0.,0.,0.]

import vtk
import vtk.util.numpy_support as vtkNumPy 
# FIXME  notice that order of operations is IMPORTANT
# FIXME   translation followed by rotation will give different results
# FIXME   than rotation followed by translation
# FIXME  Translate -> RotateZ -> RotateY -> RotateX -> Scale seems to be the order of paraview
AffineTransform = vtk.vtkTransform()
# should be in meters
AffineTransform.Translate([0.0,0.0,0.00000001])
AffineTransform.RotateZ(  0.0 )
AffineTransform.RotateY(  0.0 )
AffineTransform.RotateX(  0.0 )
AffineTransform.Scale([1.,1.,1.])
# get homogenius 4x4 matrix  of the form
#               A | b
#    matrix =   -----
#               0 | 1
#   
matrix = AffineTransform.GetConcatenatedTransform(0).GetMatrix()
#print matrix 
RotationMatrix = [[matrix.GetElement(0,0),matrix.GetElement(0,1),matrix.GetElement(0,2)],
                  [matrix.GetElement(1,0),matrix.GetElement(1,1),matrix.GetElement(1,2)],
                  [matrix.GetElement(2,0),matrix.GetElement(2,1),matrix.GetElement(2,2)]]
Translation =     [matrix.GetElement(0,3),matrix.GetElement(1,3),matrix.GetElement(2,3)] 
#print RotationMatrix ,Translation 

laserTip         =  AffineTransform.TransformPoint( [0.,0.,0.035] )
laserOrientation =  AffineTransform.TransformVector( [0.,0.,-1.0 ] )
fem.SetIniValue( "probe/x_0","%f" % laserTip[0]) 
fem.SetIniValue( "probe/y_0","%f" % laserTip[1]) 
fem.SetIniValue( "probe/z_0","%f" % laserTip[2]) 
fem.SetIniValue( "probe/x_orientation","%f" % laserOrientation[0] ) 
fem.SetIniValue( "probe/y_orientation","%f" % laserOrientation[1] ) 
fem.SetIniValue( "probe/z_orientation","%f" % laserOrientation[2] ) 

# initialize FEM Mesh
# must setup Ini File first
fem.SetupUnStructuredGrid( "phantomMesh.e",0,RotationMatrix, Translation  ) 
#fem.SetupStructuredGrid( (10,10,4) ,[0.0,1.0],[0.0,1.0],[0.0,1.0]) 

# add the data structures for the Background System Solve
# set deltat, number of time steps, power profile, and add system
nsubstep = 1
acquisitionTime = 5.0
deltat = acquisitionTime / nsubstep
ntime  = 60 
fem.AddPennesSystem(2,deltat,nsubstep,[[7,30,ntime],[0.0,4.0,0.0]]) 

# initialize libMesh data structures
fem.InitializeEquationSystems( ) 

# print info
fem.printSelf() 

# write IC
fem.WriteTimeStep("fem_data.e" , 1, 0.0 )  

# read imaging data geometry that will be used to project FEM data onto
vtkReader = vtk.vtkXMLImageDataReader() 
vtkReader.SetFileName("temperature.0000.vti" ) 
vtkReader.Update()
templateImage = vtkReader.GetOutput()
dimensions = templateImage.GetDimensions()
spacing    = templateImage.GetSpacing()
origin     = templateImage.GetOrigin()

# loop over time steps and solve
for timeID in range(1,ntime):
   print "time step = " ,timeID
   fem.UpdateTransientSystemTimeStep("StateSystem",timeID ) 
   fem.SystemSolve( "StateSystem" ) 
   #fem.StoreTransientSystemTimeStep("StateSystem",timeID ) 
   fem.WriteTimeStep("fem_data.e" , timeID+1, timeID*deltat )  

   # Interpolate FEM onto imaging data structures
   if ( petscRank == -1 ):
      vtkExodusIIReader = vtk.vtkExodusIIReader()
      vtkExodusIIReader.SetFileName("fem_data.e")
      vtkExodusIIReader.SetPointResultArrayStatus("u0",1)
      vtkExodusIIReader.SetTimeStep(timeID-1) 
      vtkExodusIIReader.Update()
 
      # reuse ShiftScale Geometry
      vtkResample = vtk.vtkCompositeDataProbeFilter()
      vtkResample.SetInput( templateImage )
      vtkResample.SetSource( vtkExodusIIReader.GetOutput() ) 
      vtkResample.Update()

      # FIXME this is prob the longest round about way possible...
      # FIXME convert to binary first then read back in a single component
      # For VTK to be able to use the data, it must be stored as a VTK-image. This can be done by the vtkImageImport-class which
      # imports raw data and stores it.
      dataImporter = vtk.vtkImageImport()
      imageData = vtkResample.GetOutput().GetPointData().GetArray('u0')
      data_array = vtkNumPy.vtk_to_numpy(imageData) 
      # write as binary also 
      #data_array.tofile("background.%04d.raw" % timeID )
      # The previously created array is converted to a string of chars and imported.
      data_string = data_array.tostring()
      dataImporter.CopyImportVoidPointer(data_string, len(data_string))
      # The type of the newly imported data is set to unsigned char (uint8)
      dataImporter.SetDataScalarTypeToDouble()
      # Because the data that is imported only contains an intensity value (it isnt RGB-coded or someting similar), the importer
      # must be told this is the case.
      dataImporter.SetNumberOfScalarComponents(1)
      # The following two functions describe how the data is stored and the dimensions of the array it is stored in. For this
      # simple case, all axes are of length 75 and begins with the first element. For other data, this is probably not the case.
      # I have to admit however, that I honestly dont know the difference between SetDataExtent() and SetWholeExtent() although
      # VTK complains if not both are used.
      print dimensions
      dataImporter.SetDataExtent( 0, dimensions[0]-1, 0, dimensions[1]-1, 0, dimensions[2]-1)
      dataImporter.SetWholeExtent(0, dimensions[0]-1, 0, dimensions[1]-1, 0, dimensions[2]-1)
      dataImporter.SetDataSpacing( spacing )
      dataImporter.SetDataOrigin(  origin  )

      # cast to float
      vtkModelPrediction = vtk.vtkImageCast() 
      vtkModelPrediction.SetInput(dataImporter.GetOutput())
      vtkModelPrediction.SetOutputScalarTypeToFloat()
      vtkModelPrediction.Update()

      # write output
      print "writing ", timeID
      vtkTemperatureWriter = vtk.vtkDataSetWriter()
      vtkTemperatureWriter.SetFileTypeToBinary()
      vtkTemperatureWriter.SetFileName("modeltemperature.%04d.vtk" % timeID )
      vtkTemperatureWriter.SetInput(vtkModelPrediction.GetOutput())
      vtkTemperatureWriter.Update()

