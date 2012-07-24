"""
background correction model 
"""

# import needed modules
import petsc4py, numpy, sys
import scipy.io as scipyio
PetscOptions = sys.argv
PetscOptions.append("-ksp_monitor") 
PetscOptions.append("-ksp_rtol") 
PetscOptions.append("1.e-10") 
# need to solve w/ constant null space for Neumann problem
#PetscOptions.append("-ksp_constant_null_space")
#PetscOptions.append("-help")
petsc4py.init(PetscOptions)

# get rank
from petsc4py import PETSc
petscRank = PETSc.COMM_WORLD.getRank()
petscSize = PETSc.COMM_WORLD.Get_size()
sys.stdout.write("petsc rank %d petsc nproc %d\n" % (petscRank, petscSize))

# set shell context
# TODO import vtk should be called after femLibrary ???? 
# FIXME WHY IS THIS????
import femLibrary
# initialize libMesh data structures
libMeshInit = femLibrary.PyLibMeshInit(PetscOptions,PETSc.COMM_WORLD) 
  
# load vtk modules to read imaging
import vtk 
import vtk.util.numpy_support as vtkNumPy 

FileNameTemplate = "/data/fuentes/biotex/MayoLiverDataJan2011/VTKData/226p-536/phase.%04d.vti"
FileNameTemplate = "/home/jyung/e137/Processed/s22000/phase.%04d.vtk"
# set the default reader based on extension
if( FileNameTemplate.split(".").pop() == "vtk"):
   vtkImageReader = vtk.vtkDataSetReader
elif( FileNameTemplate.split(".").pop() == "vti"):
   vtkImageReader = vtk.vtkXMLImageDataReader
else:
   raise RuntimeError("uknown file")

# get dimension info from header
vtkSetupReader = vtkImageReader() 
vtkSetupReader.SetFileName(FileNameTemplate % 0 ) 
vtkSetupReader.Update() 
dimensions = vtkSetupReader.GetOutput().GetDimensions()
numberPointsImage =  vtkSetupReader.GetOutput().GetNumberOfPoints()
spacing = vtkSetupReader.GetOutput().GetSpacing()
origin  = vtkSetupReader.GetOutput().GetOrigin()
print "#points",numberPointsImage , "dimensions ",dimensions , "spacing ",spacing , "origin ",origin 
# temperature map factor
alpha = +0.0097    # FIXME should be neg
# FIXME need to extract automagically at some point
#(0018|0081) Echo Time = 9.648
#(0018|0082) Inversion Time = 0
#(0018|0083) Number of Averages = 2
#(0018|0084) Imaging Frequency = 63.869849
echoTime = 9.648
imagFreq = 63.869849
tmap_factor = 1.0 / (2.0 * numpy.pi * imagFreq * alpha * echoTime * 1.e-3)
print "tmap_factor = ", tmap_factor

# return image data from raw file names
def GetNumpyPhaseData(filename):
  vtkReader = vtkImageReader() 
  vtkReader.SetFileName(filename) 
  vtkReader.Update() 
  phase_image = vtkReader.GetOutput().GetPointData() 
  # convert to phase
  phase_array =  vtkNumPy.vtk_to_numpy(phase_image.GetArray(0)) 
  #return (2.0 * numpy.pi / 4095. ) *phase_array 
  return phase_array.reshape(dimensions)

# write a numpy data to disk in vtk format
def ConvertNumpyVTKImage(NumpyImageData):
  # Create initial image
  # imports raw data and stores it.
  dataImporter = vtk.vtkImageImport()
  # array is converted to a string of chars and imported.
  data_string = NumpyImageData.tostring()
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
  dataImporter.SetDataExtent( 0, dimensions[0]-1, 0, dimensions[1]-1, 0, dimensions[2]-1)
  dataImporter.SetWholeExtent(0, dimensions[0]-1, 0, dimensions[1]-1, 0, dimensions[2]-1)
  dataImporter.SetDataSpacing( spacing )
  dataImporter.Update()
  return dataImporter.GetOutput()

# store control variables
getpot = femLibrary.PylibMeshGetPot(PetscOptions) 

# initialize FEM Mesh
femMesh = femLibrary.PylibMeshMesh()
ROI = [[40,50],   # pixel # of xbounds
       [30,50],   # pixel # of ybounds
       [ 0, 0]]   # pixel # of zbounds
npixelROI = tuple( [ (pixel[1] - pixel[0] ) for pixel in ROI] )
nnodeROI  = tuple( [ (pixel[1] - pixel[0] + 1 ) for pixel in ROI] )
xbounds = [ origin[0], origin[0] + spacing[0] * npixelROI[0] ]
ybounds = [ origin[1], origin[1] + spacing[1] * npixelROI[1] ]
zbounds = [ origin[2], origin[2] + spacing[2] * npixelROI[2] ]
# setup structure Grid expect # of elements
femMesh.SetupStructuredGrid( npixelROI, 
                             xbounds,ybounds,zbounds,[1,1,1,1,1,1]) 

# setup imaging to interpolate onto FEM mesh
femImaging = femLibrary.PytttkImaging(getpot, dimensions ,origin,spacing) 

# add the data structures for the Background System Solve
eqnSystems =  femLibrary.PylibMeshEquationSystems(femMesh,getpot)
eqnSystems.AddBackgroundSystem( "Background" ) 
eqnSystems.AddExplicitSystem( "ImageMask" ,1,1 ) 
# initialize libMesh data structures
eqnSystems.init( ) 
# print info
eqnSystems.PrintSelf() 
  
# create data structure for cumulative sum w/o background correction
TempSumOrig = numpy.zeros(dimensions,dtype=numpy.float)
# create data structure for cumulative sum w/ background correction
TempSum     = numpy.zeros(dimensions,dtype=numpy.float)
# shift/scale the phase data set data type to float
phase_prev = GetNumpyPhaseData(FileNameTemplate % 0 ) 

# write IC
MeshOutputFile = "fem_data.e"
exodusII_IO = femLibrary.PylibMeshExodusII_IO(femMesh)
exodusII_IO.WriteTimeStep(MeshOutputFile,eqnSystems, 1, 0.0 )  

# loop over desired time instances
ntime=93
for timeID in range(0,ntime+1):
   print "working on time id %d " % timeID
   # read in data
   phase_curr = GetNumpyPhaseData(FileNameTemplate % timeID ) 
   
   # pixel wise subtract to create deltat image 
   delta_temp =tmap_factor * (phase_curr - phase_prev) 
   # store for next time
   phase_prev = phase_curr 

   # write original delta tmap 
   vtkTempImage = ConvertNumpyVTKImage(delta_temp)
   vtkTmpWriter = vtk.vtkDataSetWriter()
   vtkTmpWriter.SetFileName("origdeltatemp.%04d.vtk" % timeID )
   vtkTmpWriter.SetInput(vtkTempImage)
   vtkTmpWriter.Update()

   # write original tmap
   TempSumOrig = TempSumOrig + delta_temp
   vtkTempImage = ConvertNumpyVTKImage(TempSumOrig )
   vtkTmpWriter = vtk.vtkDataSetWriter()
   vtkTmpWriter.SetFileName("origtemp.%04d.vtk" % timeID )
   vtkTmpWriter.SetInput(vtkTempImage )
   vtkTmpWriter.Update()

   # get image mask
   vtkImageMask = vtk.vtkImageThreshold() 
   vtkImageMask.ReplaceOutOn()
   vtkImageMask.SetOutValue(0.0)
   vtkImageMask.ReplaceInOn()
   vtkImageMask.SetInValue(1.0)
   vtkImageMask.SetOutputScalarTypeToDouble()
   try:
     #default look for a user defined image mask
     maskFileName = "%s/SNRuncert.%04d.vtk" %  (
                   FileNameTemplate[:FileNameTemplate.rfind("/")],timeID)
     vtkMaskReader = vtk.vtkDataSetReader() 
     vtkMaskReader.SetFileName( maskFileName ) 
     vtkMaskReader.Update() 
     # set threshold
     vtkImageMask.ThresholdByLower( 100.0)
     vtkImageMask.SetInput(vtkMaskReader.GetOutput())
   except: # if nothing available threshold the phase image
     # take gradient
     vtkPhaseData = ConvertNumpyVTKImage( phase_curr )
     vtkGradientImage = vtk.vtkImageGradient() 
     vtkGradientImage.SetInput(vtkPhaseData) 
     vtkGradientImage.Update() 
     
     # take magnitude of gradient
     vtkImageNorm = vtk.vtkImageMagnitude() 
     vtkImageNorm.SetInput(vtkGradientImage.GetOutput())
     vtkImageNorm.Update() 
     # set threshold
     vtkImageMask.ThresholdByLower( 100.0* tmap_factor * 2.0*numpy.pi/4095.)
     vtkImageMask.SetInput(vtkImageNorm.GetOutput())
   vtkImageMask.Update( )
   # check output
   vtkWriter = vtk.vtkXMLImageDataWriter()
   vtkWriter.SetFileName("threshold.%04d.vti" % timeID )
   vtkWriter.SetInput( vtkImageMask.GetOutput() )
   vtkWriter.Update()

   # pass pointer to c++
   image_cells = vtkImageMask.GetOutput().GetPointData() 
   data_array = vtkNumPy.vtk_to_numpy( image_cells.GetArray(0) ) 
   v1 = PETSc.Vec().createWithArray(phase_curr, comm=PETSc.COMM_SELF)
   v2 = PETSc.Vec().createWithArray(data_array, comm=PETSc.COMM_SELF)

   # Project imaging onto libMesh data structures
   femImaging.ProjectImagingToFEMMesh("Background",0.0,v1,eqnSystems)  
   femImaging.ProjectImagingToFEMMesh("ImageMask" ,0.0,v2,eqnSystems)  
   eqnSystems.SystemSolve( "Background" ) 
   exodusII_IO.WriteTimeStep(MeshOutputFile,eqnSystems, timeID+1, timeID )  

   # get libMesh Background Solution as numpy data structure
   maxwell_array = eqnSystems.GetSolutionVector( "Background" )[...]
   maxwell_data  = phase_curr 
   maxwell_data[ROI[0][0]:(ROI[0][1]+1),
                ROI[1][0]:(ROI[1][1]+1),
                ROI[2][0]:(ROI[2][1]+1)] = maxwell_array.reshape(nnodeROI)
   # write numpy to disk in matlab
   scipyio.savemat("background.%04d.mat"%(timeID), {'maxwell':maxwell_array} )
   # check output
   vtkTempImage = ConvertNumpyVTKImage(maxwell_data)
   vtkWriterTmpTwo = vtk.vtkXMLImageDataWriter()
   vtkWriterTmpTwo.SetFileName("maxwell.%04d.vti" % timeID )
   vtkWriterTmpTwo.SetInput( vtkTempImage )
   vtkWriterTmpTwo.Update()

   # compute temperature difference
   delta_temp =tmap_factor * (phase_curr - maxwell_data) 

   # check output
   vtkTempImage = ConvertNumpyVTKImage(delta_temp)
   vtkWriterTmpTwo = vtk.vtkXMLImageDataWriter()
   vtkWriterTmpTwo.SetFileName("deltat.%04d.vti" % timeID )
   vtkWriterTmpTwo.SetInput( vtkTempImage )
   vtkWriterTmpTwo.Update()

   # write tmap
   TempSum = TempSum + delta_temp
   vtkTempImage = ConvertNumpyVTKImage( TempSum )
   # check output
   vtkTmpWriter3 = vtk.vtkXMLImageDataWriter()
   vtkTmpWriter3.SetFileName("temperature.%04d.vti" % timeID )
   vtkTmpWriter3.SetInput( vtkTempImage )
   vtkTmpWriter3.Update()
