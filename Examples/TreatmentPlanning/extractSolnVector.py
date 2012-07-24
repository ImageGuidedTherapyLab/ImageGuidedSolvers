#import vtk module
import vtk
import vtk.util.numpy_support as vtkNumPy
# echo vtk version info
print "using vtk version", vtk.vtkVersion.GetVTKVersion()

# load FEM data
vtkExodusIIReader = vtk.vtkExodusIIReader()
vtkExodusIIReader.SetFileName("fem_data.e")
vtkExodusIIReader.SetPointResultArrayStatus("u0",1)
vtkExodusIIReader.Update()
nblock = vtkExodusIIReader.GetNumberOfElementBlockArrays()
ntime = vtkExodusIIReader.GetNumberOfTimeSteps()

#extract the data to a file 
maxData = None
for idTime in range(ntime):
  vtkExodusIIReader.SetTimeStep(idTime) 
  vtkExodusIIReader.Update()
  femData = vtkExodusIIReader.GetOutput()
  rawdata = []
  if femData.IsA("vtkMultiBlockDataSet"):
    iter = femData.NewIterator()
    iter.UnRegister(None)
    iter.InitTraversal()
    block = 0 
    while not iter.IsDoneWithTraversal():
      curInput = iter.GetCurrentDataObject()
      #print "number of points" , curInput.GetNumberOfPoints()
      curData = curInput.GetPointData().GetArray(0)
      rawdata.append(
                [curData.GetValue(i) for i in range(curData.GetSize())] ) 
      iter.GoToNextItem();
      block = block + 1 
  else:
    raise ValueError("unknown data type")
  if maxData == None:
     maxData = rawdata
  else:
     tmpData = maxData 
     for idblock in range(block):
        for i,(x,y) in enumerate(zip(tmpData[idblock],rawdata[idblock])):
            maxData[idblock][i] = max(x,y)
        #print max(maxData[idblock])

# write out max to a file as the only "timestep"
vtkExodusIIWriter = vtk.vtkExodusIIWriter()
vtkExodusIIWriter.SetFileName("max_fem_data.e" )
vtkExodusIIWriter.WriteAllTimeStepsOff()
vtkExodusIIWriter.WriteOutBlockIdArrayOn ();
vtkExodusIIWriter.WriteOutGlobalNodeIdArrayOn ();
vtkExodusIIWriter.WriteOutGlobalElementIdArrayOn ();
#vtkExodusIIWriter.WriteAllTimeStepsOn ();

# needs >= VTK 5.6.1 to get vtkexodusIIwriter to work for multiblock data
# overwrite w max data
if femData.IsA("vtkMultiBlockDataSet"):
  iter = femData.NewIterator()
  iter.UnRegister(None)
  iter.InitTraversal()
  block = 0 
  while not iter.IsDoneWithTraversal():
    curInput = iter.GetCurrentDataObject()
    #print "number of points" , curInput.GetNumberOfPoints()
    curData = curInput.GetPointData().GetArray(0)
    for i in range(curData.GetSize()):
        curData.SetValue(i,maxData[block][i]) 
    iter.GoToNextItem();
    block = block + 1 

#vtkExodusIIWriter.SetInput(curInput)
vtkExodusIIWriter.SetInput(femData)
vtkExodusIIWriter.Update()
