#import vtk module
import vtk
import vtk.util.numpy_support as vtkNumPy 

# load FEM data
vtkExodusIIReader = vtk.vtkExodusIIReader()
vtkExodusIIReader.SetFileName("fem_data.e")
vtkExodusIIReader.Update()
ntime = vtkExodusIIReader.GetNumberOfTimeSteps()
femData = vtkExodusIIReader.GetOutput()

#get the bounding box
femBounds = None
if femData.IsA("vtkMultiBlockDataSet"):
    iter = femData.NewIterator()
    iter.UnRegister(None)
    iter.InitTraversal()
    while not iter.IsDoneWithTraversal():
        curInput = iter.GetCurrentDataObject()
        bounds = curInput.GetBounds()
        if femBounds == None:
           femBounds = bounds
        else:
           femBounds = (
                min(femBounds[0],bounds[0]),
                max(femBounds[1],bounds[1]),
                min(femBounds[2],bounds[2]),
                max(femBounds[3],bounds[3]),
                min(femBounds[4],bounds[4]),
                max(femBounds[5],bounds[5])
                       )
        #print bounds
        iter.GoToNextItem();


# set imaging dimensions and echo FEM Data
dimensions = (100,100,20)
origin = (femBounds[0], femBounds[2], femBounds[4])
spacing = ( (femBounds[1]-femBounds[0])/ dimensions[0] ,
            (femBounds[3]-femBounds[2])/ dimensions[1] ,
            (femBounds[5]-femBounds[4])/ dimensions[2]  
          )
print femBounds, origin, spacing
print ntime

# create imaging data structure
import numpy 
data_matrix = numpy.zeros(dimensions,dtype=numpy.float)
# For VTK to be able to use the data, it must be stored as a VTK-image. This can be done by the vtkImageImport-class which
# imports raw data and stores it.
dataImporter = vtk.vtkImageImport()
# The preaviusly created array is converted to a string of chars and imported.
data_string = data_matrix.tostring()
dataImporter.CopyImportVoidPointer(data_string, len(data_string))
# The type of the newly imported data is set to unsigned char (uint8)
dataImporter.SetDataScalarTypeToFloat()
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
dataImporter.SetDataOrigin(  origin  )
dataImporter.Update()
    
# write one image file for each timestep in mesh
for timeID in range(1,ntime):
  vtkExodusIIReader.SetPointResultArrayStatus("u0",1)
  vtkExodusIIReader.SetTimeStep(timeID-1) 
  vtkExodusIIReader.Update()
  # reuse ShiftScale Geometry
  vtkResample = vtk.vtkCompositeDataProbeFilter()
  vtkResample.SetInput( dataImporter.GetOutput() )
  vtkResample.SetSource( vtkExodusIIReader.GetOutput() ) 
  vtkResample.Update()
  
  # FIXME this is prob the longest round about way possible...
  # FIXME convert to binary first then read back in a single component
  # For VTK to be able to use the data, it must be stored as a VTK-image. This can be done by the vtkImageImport-class which
  # imports raw data and stores it.
  FEMdataImporter = vtk.vtkImageImport()
  imageData = vtkResample.GetOutput().GetPointData().GetArray('u0')
  data_array = vtkNumPy.vtk_to_numpy(imageData) 
  # write as binary also 
  #data_array.tofile("background.%04d.raw" % timeID )
  # The previously created array is converted to a string of chars and imported.
  data_string = data_array.tostring()
  FEMdataImporter.CopyImportVoidPointer(data_string, len(data_string))
  # The type of the newly imported data is set to unsigned char (uint8)
  FEMdataImporter.SetDataScalarTypeToDouble()
  # Because the data that is imported only contains an intensity value (it isnt RGB-coded or someting similar), the importer
  # must be told this is the case.
  FEMdataImporter.SetNumberOfScalarComponents(1)
  # The following two functions describe how the data is stored and the dimensions of the array it is stored in. For this
  # simple case, all axes are of length 75 and begins with the first element. For other data, this is probably not the case.
  # I have to admit however, that I honestly dont know the difference between SetDataExtent() and SetWholeExtent() although
  # VTK complains if not both are used.
  print dimensions
  FEMdataImporter.SetDataExtent( 0, dimensions[0]-1, 0, dimensions[1]-1, 0, dimensions[2]-1)
  FEMdataImporter.SetWholeExtent(0, dimensions[0]-1, 0, dimensions[1]-1, 0, dimensions[2]-1)
  FEMdataImporter.SetDataSpacing( spacing )
  FEMdataImporter.SetDataOrigin(  origin  )
  
  # cast to float
  vtkModelPrediction = vtk.vtkImageCast() 
  vtkModelPrediction.SetInput(FEMdataImporter.GetOutput())
  vtkModelPrediction.SetOutputScalarTypeToFloat()
  vtkModelPrediction.Update()
  
  # write output
  print "writing ", timeID
  vtkTemperatureWriter = vtk.vtkDataSetWriter()
  vtkTemperatureWriter.SetFileTypeToBinary()
  vtkTemperatureWriter.SetFileName("modeltemperature.%04d.vtk" % timeID )
  vtkTemperatureWriter.SetInput(vtkModelPrediction.GetOutput())
  vtkTemperatureWriter.Update()

# save as matlab file
import scipy.io as sio
vect = numpy.arange(10)
print vect.shape
sio.savemat('vector.mat', {'vect':vect})
