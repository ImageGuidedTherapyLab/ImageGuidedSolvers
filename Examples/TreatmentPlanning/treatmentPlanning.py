import os
import ConfigParser
####################################################################
def pennesModeling(**kwargs):
  """
  treatment planning model 
  """
  # import needed modules
  import petsc4py, numpy, sys
  PetscOptions =  sys.argv
  PetscOptions.append("-ksp_monitor")
  PetscOptions.append("-ksp_rtol")
  PetscOptions.append("1.0e-15")
  #PetscOptions.append("-idb")
  petsc4py.init(PetscOptions)
  
  from petsc4py import PETSc
  
  # break processors into separate communicators
  petscRank = PETSc.COMM_WORLD.getRank()
  petscSize = PETSc.COMM_WORLD.Get_size()
  sys.stdout.write("petsc rank %d petsc nproc %d\n" % (petscRank, petscSize))
  
  # set shell context
  import femLibrary
  # initialize libMesh data structures
  libMeshInit = femLibrary.PyLibMeshInit(PetscOptions,PETSc.COMM_WORLD) 
  
  # the original configuration ini file should be stored
  config = kwargs['config_parser']

  # store control variables
  getpot = femLibrary.PylibMeshGetPot(PetscOptions) 

  # copy all values from the input file
  for section in config.sections():
    for name,value in  config.items(section):
      #print "%s/%s" % (section,name) , value
      getpot.SetIniValue( "%s/%s" % (section,name) , value ) 

  # set tissue lookup tables
  k_0Table  = {"default":config.getfloat("thermal_conductivity","k_0_healthy")  ,
               "grey"   :config.getfloat("thermal_conductivity","k_0_grey"   )  ,
               "white"  :config.getfloat("thermal_conductivity","k_0_white"  )  ,
               "csf"    :config.getfloat("thermal_conductivity","k_0_csf"    )  ,
               "tumor"  :config.getfloat("thermal_conductivity","k_0_tumor"  )  }
  w_0Table  = {"default":config.getfloat("perfusion","w_0_healthy")  ,
               "grey"   :config.getfloat("perfusion","w_0_grey"   )  ,
               "white"  :config.getfloat("perfusion","w_0_white"  )  ,
               "csf"    :config.getfloat("perfusion","w_0_csf"    )  ,
               "tumor"  :config.getfloat("perfusion","w_0_tumor"  )  }
  mu_aTable = {"default":config.getfloat("optical","mu_a_healthy")  ,
               "grey"   :config.getfloat("optical","mu_a_grey"   )  ,
               "white"  :config.getfloat("optical","mu_a_white"  )  ,
               "csf"    :config.getfloat("optical","mu_a_csf"    )  ,
               "tumor"  :config.getfloat("optical","mu_a_tumor"  )  }
  mu_sTable = {"default":config.getfloat("optical","mu_s_healthy")  ,
               "grey"   :config.getfloat("optical","mu_s_grey"   )  ,
               "white"  :config.getfloat("optical","mu_s_white"  )  ,
               "csf"    :config.getfloat("optical","mu_s_csf"    )  ,
               "tumor"  :config.getfloat("optical","mu_s_tumor"  )  }
  labelTable= {config.get("labels","greymatter" ):"grey" , 
               config.get("labels","whitematter"):"white", 
               config.get("labels","csf"        ):"csf"  , 
               config.get("labels","tumor"      ):"tumor", 
               config.get("labels","vessel"     ):"vessel"}

  # echo lookup table
  if( petscRank ==0 ):
     print "lookup tables"
     print "labels"      , labelTable
     print "conductivity", k_0Table  
     print "perfusion"   , w_0Table  
     print "absorption"  , mu_aTable  
     print "scattering"  , mu_sTable  
  # imaging params
  import vtk 
  import vtk.util.numpy_support as vtkNumPy 
  SegmentFile=config.get("exec","segment_file")
  # set the default reader based on extension
  if( SegmentFile.split(".").pop() == "vtk"):
     vtkImageReader = vtk.vtkDataSetReader
  elif( SegmentFile.split(".").pop() == "vti"):
     vtkImageReader = vtk.vtkXMLImageDataReader
  else:
     raise RuntimeError("uknown file")
  
  # get dimension info from header
  vtkSetupReader = vtkImageReader() 
  vtkSetupReader.SetFileName(SegmentFile ) 
  vtkSetupReader.Update() 
  vtkImageMask = vtkSetupReader.GetOutput()
  dimensions = vtkSetupReader.GetOutput().GetDimensions()
  numberPointsImage =  vtkSetupReader.GetOutput().GetNumberOfPoints()
  spacing_mm = vtkSetupReader.GetOutput().GetSpacing()
  origin_mm  = vtkSetupReader.GetOutput().GetOrigin()
  # convert to meters
  spacing = [dx*.001 for dx in spacing_mm]
  origin  = [x0*.001 for x0 in  origin_mm]
  
  # pass pointer to c++
  image_cells = vtkImageMask.GetPointData() 
  data_array = vtkNumPy.vtk_to_numpy( image_cells.GetArray(0) ) 
  # need to pass numpy array's w/ Fortran storage... ie painful to debug
  v1 = PETSc.Vec().createWithArray(numpy.ravel(data_array,order='F'), comm=PETSc.COMM_SELF)

  # FIXME - center around quadrature in out-of-plane direction
  # FIXME - need better out of plane cooling model
  quadratureOffset  = 1./numpy.sqrt(3.0) * spacing[2]/2.0
  # expecting roi and subsample of the form:
  #     roi = [(40,210),(30,220),(6,78)]
  #     subsample = [3,3,2]
  ROI = eval(config.get('exec','roi'))
  subsample = eval(config.get('exec','subsample'))
  nelemROI  = [ (pixel[1] - pixel[0] - 1 )/sub for pixel,sub in zip(ROI,subsample)] 
  xbounds = [ origin[0]+spacing[0]*(ROI[0][0]+0.5),origin[0]+spacing[0]*(ROI[0][1]+0.5) ]
  ybounds = [ origin[1]+spacing[1]*(ROI[1][0]+0.5),origin[1]+spacing[1]*(ROI[1][1]+0.5) ]
  zbounds = [ origin[2]+spacing[2]*(ROI[2][0]+0.5),origin[2]+spacing[2]*(ROI[2][1]+0.5) ]
  if( petscRank ==0 ):
    print "#points",numberPointsImage , "dimensions ",dimensions , "spacing ",spacing , "origin ",origin 
    print "ROI",ROI , "nelemROI ",nelemROI , "bounds", xbounds, ybounds, zbounds

  #set to steady state solve 
  getpot.SetIniValue("steadystate/domain_0","true") 

  # initialize FEM Mesh
  femMesh = femLibrary.PylibMeshMesh()
  #femMesh.SetupUnStructuredGrid(kwargs['mesh_file'],0,RotationMatrix, Translation  ) 
  #femMesh.ReadFile(kwargs['mesh_file'])

  femMesh.SetupStructuredGrid(nelemROI,xbounds,ybounds,zbounds,
                              [2,2,2,2,2,2]) 
  MeshOutputFile = "fem.e"

  # add the data structures for the Background System Solve
  # set deltat, number of time steps, power profile, and add system
  eqnSystems =  femLibrary.PylibMeshEquationSystems(femMesh,getpot)
    
  #getpot.SetIniPower(1,[[19,119],[90.0,100.0]] )
  getpot.SetIniPower(1,kwargs['powerHistory'] )
  # AddPennesSDASystem
  # AddPennesRFSystem
  # AddPennesDeltaPSystem
  pennesSystem = eval("eqnSystems.%s('StateSystem',kwargs['deltat'])" %  kwargs['physics'])
  pennesSystem.AddStorageVectors( kwargs['ntime']+1 ) 
  
  # add system for labels 
  # FIXME: need both nodal, for dirichlet bc, and element version, for parameter masks
  maskElemSystem = eqnSystems.AddExplicitSystem( "ElemImageMask" ) 
  maskElemSystem.AddConstantMonomialVariable( "maskElem" ) 

  # initialize libMesh data structures
  eqnSystems.init( ) 
  
  # print info
  eqnSystems.PrintSelf() 
  
  # setup imaging to interpolate onto FEM mesh
  femImaging = femLibrary.PytttkImaging(getpot, dimensions ,origin,spacing) 
  # Project imaging onto libMesh data structures
  femImaging.ProjectImagingToFEMMesh("ElemImageMask" ,0.0,v1,eqnSystems)  
  femImaging.ProjectImagingToFEMMesh("StateSystem" ,0.0,v1,eqnSystems)  

  # create dirichlet nodes from this mask
  numNodes = pennesSystem.PetscFEMSystemCreateNodeSetFromMask(1.0,1)
  print "# of dirichlet nodes %d" %numNodes 


  # get image label as numpy array
  imageLabel = maskElemSystem.GetSolutionVector()[...]
  imageLabel = numpy.floor(10.0*imageLabel.copy())+1
  k_0Label  = numpy.floor(10.0*imageLabel.copy())+1
  w_0Label  = numpy.floor(10.0*imageLabel.copy())+1
  mu_aLabel = numpy.floor(10.0*imageLabel.copy())+1
  mu_sLabel = numpy.floor(10.0*imageLabel.copy())+1
  for (idpos,label) in enumerate(imageLabel):
     try:
       tissueType = labelTable[int(label)]
     except KeyError:
       tissueType = "default"
     k_0Label[ idpos] = k_0Table[ tissueType]
     w_0Label[ idpos] = w_0Table[ tissueType]
     mu_aLabel[idpos] = mu_aTable[tissueType]
     mu_sLabel[idpos] = mu_sTable[tissueType]
  # create array of imaging data as petsc vec
  k_0Vec  = PETSc.Vec().createWithArray(k_0Label , comm=PETSc.COMM_SELF)
  w_0Vec  = PETSc.Vec().createWithArray(w_0Label , comm=PETSc.COMM_SELF)
  mu_aVec = PETSc.Vec().createWithArray(mu_aLabel, comm=PETSc.COMM_SELF)
  mu_sVec = PETSc.Vec().createWithArray(mu_sLabel, comm=PETSc.COMM_SELF)
  # copy material properties to the system
  eqnSystems.GetSystem("k_0").SetSolutionVector( k_0Vec )
  eqnSystems.GetSystem("w_0").SetSolutionVector( w_0Vec )
  eqnSystems.GetSystem("mu_a").SetSolutionVector(mu_aVec)
  eqnSystems.GetSystem("mu_s").SetSolutionVector(mu_sVec)

  # write IC
  exodusII_IO = femLibrary.PylibMeshExodusII_IO(femMesh)
  exodusII_IO.WriteTimeStep(MeshOutputFile,eqnSystems, 1, 0.0 )  
  exodusII_IO.WriteParameterSystems(eqnSystems)  
  
  # setup IC 
  pennesSystem.PetscFEMSystemSetupInitialConditions( ) 

  # create list of laser positions
  laserPositionList = [ ((.015,.015,2*spacing[2]+quadratureOffset),
                         (.018,.018,2*spacing[2]+quadratureOffset) ),
                        ((.015,.015,2*spacing[2]+quadratureOffset),
                         (.016,.018,2*spacing[2]+quadratureOffset) ),
                        ((.015,.015,1*spacing[2]+quadratureOffset),
                         (.014,.014,1*spacing[2]+quadratureOffset) ),
                        ((.015,.015,1*spacing[2]+quadratureOffset),
                         (.015,.018,3*spacing[2]+quadratureOffset) )]
  # time stamp
  import pickle,time
  timeStamp =0 
  # set power id to turn laser on 
  pennesSystem.PetscFEMSystemUpdateTimeStep( 1 ) 
  # loop over laser position and solve steady state equations for each position
  #for (idpos,(pos0,pos1)) in enumerate(laserPositionList):
  # loop and read new laser parameters
  while(True):
    if(os.path.getmtime(fem_params['ini_filename']) > timeStamp):
      timeStamp = os.path.getmtime(fem_params['ini_filename'] ) 
      newIni = ConfigParser.SafeConfigParser({})
      newIni.read(fem_params['ini_filename'])

      laserParams = {}
      laserParams['position1'] = [newIni.getfloat("probe",varID ) for varID in ["x_0","y_0","z_0"] ]
      laserParams['position2'] = [newIni.getfloat("probe",varID ) for varID in ["x_1","y_1","z_1"] ]

      print "laser position = ", newIni.getfloat("timestep","power") , laserParams
      pennesSystem.PennesSDASystemUpdateLaserPower(newIni.getfloat("timestep","power"),1)
      pennesSystem.PennesSDASystemUpdateLaserPosition(laserParams['position1'],laserParams['position2'])
      pennesSystem.SystemSolve( ) 
      #fem.StoreTransientSystemTimeStep("StateSystem",timeID ) 
      #exodusII_IO.WriteTimeStep(MeshOutputFile ,eqnSystems, idpos+2, idpos*kwargs['deltat'])  
      exodusII_IO.WriteTimeStep(MeshOutputFile ,eqnSystems, 2, 1.0)  
      exodusII_IO.WriteParameterSystems(eqnSystems)  
      # write to txt file
      if( petscRank ==0 ):
        time.sleep(1)
        vtkExodusIIReader = vtk.vtkExodusIIReader()
        print "opening %s " % MeshOutputFile 
        vtkExodusIIReader.SetFileName( MeshOutputFile )
        vtkExodusIIReader.Update()
        ntime  = vtkExodusIIReader.GetNumberOfTimeSteps()
        variableID = "u0"
        print "ntime %d %s " % (ntime,variableID)
        vtkExodusIIReader.SetTimeStep(1) 
        vtkExodusIIReader.SetPointResultArrayStatus(variableID,1)
        vtkExodusIIReader.Update()
        curInput = None
        # multi block
        if vtkExodusIIReader.GetOutput().IsA("vtkMultiBlockDataSet"):
          iter = vtkExodusIIReader.GetOutput().NewIterator()
          iter.UnRegister(None)
          iter.InitTraversal()
          curInput = iter.GetCurrentDataObject()
        else: 
          curInput = vtkExodusIIReader.GetOutput()
        #fem_point_data= curInput.GetPointData().GetArray('u0') 
        #Soln=vtkNumPy.vtk_to_numpy(fem_point_data)
        #numpy.savetxt( "temperature.txt" ,Soln)
        vtkContour = vtk.vtkContourFilter()
        vtkContour.SetInput( curInput )
        # TODO: not sure why this works...
        # set the array to process at the temperature == u0
        vtkContour.SetInputArrayToProcess(0,0,0,0,'u0')
        contourValuesList  = eval(config.get('exec','contours'))
        vtkContour.SetNumberOfContours( len(contourValuesList ) )
        print "plotting array:", vtkContour.GetArrayComponent( )
        for idContour,contourValue in enumerate(contourValuesList):
           print "plotting contour:",idContour,contourValue
           vtkContour.SetValue( idContour,contourValue )
        vtkContour.Update( )
        stlWriter = vtk.vtkSTLWriter()
        stlWriter.SetInput(vtkContour.GetOutput( ))
        stlWriter.SetFileName("fem.stl")
        stlWriter.SetFileTypeToBinary()
        stlWriter.Write()
    else:
      print "waiting on user input.."
      time.sleep(2)
## end def pennesModeling(**kwargs)
####################################################################

# setup command line parser to control execution
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
parser.add_option( "--ini", 
                  action="store", dest="config_ini", default=None,
                  help="ini FILE containing setup info", metavar="FILE")
(options, args) = parser.parse_args()

# run planning solver w/ default options from ini file
if (options.config_ini != None):

  # read config file
  config = ConfigParser.SafeConfigParser({})
  config.read(options.config_ini)

  
  fem_params = {}
  fem_params['powerHistory']  = [[1,2],[0.0,config.getfloat('timestep','power')]]
  fem_params['deltat']        =  1.0
  fem_params['ntime']         =  fem_params['powerHistory'][0][-1]
  fem_params['nsubstep']      =  1
  fem_params['physics']       = "AddPennesSDASystem"
  fem_params['u_init']        = 37.0
  fem_params['fileID']        = "0"
  # store the entire configuration file for convienence
  fem_params['config_parser'] = config
  fem_params['ini_filename'] = options.config_ini

  pennesModeling(**fem_params)
else:
  parser.print_help()
  print options
