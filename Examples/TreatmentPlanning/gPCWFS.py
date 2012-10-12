# Read DAKOTA parameters file (aprepro or standard format) and call a
# Read DAKOTA parameters file (aprepro or standard format) and call a
# Python module for fem analysis.

# DAKOTA will execute this script as
#   gPCWFS.py params.in results.out
# so sys.argv[1] will be the parameters file and
#    sys.argv[2] will be the results file to return to DAKOTA

# necessary python modules
import sys
import re
import os
import pickle
import ConfigParser

SolnOutputTemplate = "variable.%s.%04d"
MeshOutputTemplate = "fem_data.%04d.e" 
####################################################################
def RunSubProcessQueue(CodeExeCmds,ErrorFile):
  import subprocess
  import time
  # max number of local jobs
  numlocalJob = 100 
  process = []; 
  # create an error log file
  errorFile=open(ErrorFile ,"w")
  while ( len(CodeExeCmds) > 0 or len(process) > 0 ):
      # only run numlocalJob at a time
      if (len(process) > numlocalJob):
        raise RuntimeError("\n\n running too many jobs at a time??")
      elif (len(process) == numlocalJob):
        print len(CodeExeCmds), " jobs remaining..."
        time.sleep(30) # pause wait for jobs to finish
      elif( len(CodeExeCmds) > 0 ):
        cmd = CodeExeCmds.pop(0)
        print "running " , cmd
        process.append( [subprocess.Popen(cmd,shell=True),cmd] )
      if( len(process) > 0 ):
        runningJob = process.pop(0)
        if ( runningJob[0].poll() == None ):
          # job not done put it back in the list
          # not that we pop from the front of the list and pushback at the 
          # end to cycle through
          print " pid ",runningJob[0].pid, " still running"
          process.append( runningJob )
          time.sleep(2) # pause 
        elif ( runningJob[0].poll() == 0 ):
          pass # job is done 
        else:
          print "job exiting with ", runningJob[0].poll() 
          errorFile.write("error in  %s   \n" % runningJob[1] )
          #raise RuntimeError("\n\n unknown exit code ")
  errorFile.close; errorFile.flush() 
####################################################################
def GetMeshNodes(file_name):
  """
  return the number of DOF for a exodus file
  """
  import vtk
  vtkExodusIIReader = vtk.vtkExodusIIReader()
  vtkExodusIIReader.SetFileName(file_name)
  #vtkExodusIIReader.SetPointResultArrayStatus("u0",1)
  vtkExodusIIReader.Update()
  # multi block
  if vtkExodusIIReader.GetOutput().IsA("vtkMultiBlockDataSet"):
    iter = vtkExodusIIReader.GetOutput().NewIterator()
    iter.UnRegister(None)
    iter.InitTraversal()
    # initialize list for storage
    NumberOfNodes = 0
    # loop over blocks...
    while not iter.IsDoneWithTraversal():
      curInput = iter.GetCurrentDataObject()
      NumberOfNodes = NumberOfNodes + curInput.GetNumberOfPoints()
      iter.GoToNextItem();
    # return total number of nodes
    return NumberOfNodes 
  # single block
  else:
    return vtkExodusIIReader.GetNumberOfNodes()
####################################################################
def AssembleStatistics(data_dir):
  """
  collect all statistics into one exodus file
  """
  # post process stats on FEM mesh
  # import petsc and numpy
  import vtk, numpy
  import vtk.util.numpy_support as vtkNumPy 

  # initialize FEM Mesh
  pkl_file = open('%s/CaseInfo.pkl' % ".", 'rb')
  fem_params = pickle.load(pkl_file)
  pkl_file.close()
   
  vtkExodusIIReader = vtk.vtkExodusIIReader()
  print "opening %s " % fem_params['mesh_file'] 
  vtkExodusIIReader.SetFileName( fem_params['mesh_file'] )
  #vtkExodusIIReader.SetFileName( "/data/fuentes/utsa/vasculature_july10/vessel_0/realization.1//fem_data.0001.e")
  vtkExodusIIReader.ExodusModelMetadataOn ()
  vtkExodusIIReader.Update()
  exodusObject = vtkExodusIIReader.GetOutput()
  responseLevelVarList = fem_params['responseLevelVarList'] 
  probabilityLevelList = fem_params['probabilityLevelList'] 
  reliabilityLevelList = fem_params['reliabilityLevelList'] 

  # loop over time steps and import data
  # vtkTemporalDataSet = vtk.vtkTemporalDataSet()
  # multiBlockData = {}
  #for timeID in [69]:
  #TODO how can we assemble all time steps at once ? 
  #for timeID in range(0,fem_params['ntime']+1):
  for timeID in [int(data_dir.split(".").pop())]:
    for variable,responseList in responseLevelVarList:
      #(variable,responseList) = responseLevelVarList[0]
      basePath = "%s/%s" % (data_dir,variable)
      # mean
      meanFile=open("%s/meanFile.txt" % (basePath) ,"r")
      meandataList = [float(line.strip()) for line in meanFile]
      meanFile.close()
      # std dev
      stdFile=open("%s/stddFile.txt" % (basePath),"r")
      stdddataList = [float(line.strip()) for line in stdFile]
      stdFile.close()
      # skewness
      skewFile=open("%s/skewFile.txt" % (basePath),"r")
      skewdataList = [float(line.strip()) for line in skewFile]
      skewFile.close()
      # kurtosis
      kurtFile=open("%s/kurtFile.txt" % (basePath),"r")
      kurtdataList = [float(line.strip()) for line in kurtFile]
      kurtFile.close()
      # response data
      responsedataList = []
      for iii,response in enumerate(responseList):
        responseFile=open("%s/response.%d.txt" %(basePath,iii),"r")
        singleresponseList = [float(line.strip()) for line in responseFile]
        responseFile.close()
        responsedataList.append(singleresponseList)
      # probability data
      probabilitydataList = []
      for iii,probability in enumerate(probabilityLevelList):
        probFile=open("%s/probability.%d.txt" %(basePath,iii),"r")
        probList = [float(line.strip()) for line in probFile]
        probFile.close()
        probabilitydataList.append(probList)
      # reliability data
      reliabilitydataList = []
      for iii,reliability in enumerate(reliabilityLevelList):
        reliabFile=open("%s/reliability.%d.txt" %(basePath,iii),"r")
        reliabList = [float(line.strip()) for line in reliabFile]
        reliabFile.close()
        reliabilitydataList.append(reliabList) 
      # multi block
      if exodusObject.IsA("vtkMultiBlockDataSet"):
        iter = exodusObject.NewIterator()
        iter.UnRegister(None)
        iter.InitTraversal()
        # iter.GoToNextItem();
        # iter.GoToNextItem();
        metadata = exodusObject.GetMetaData(iter)
        # initialize list for storage
        listSliceInit = 0
        # loop over blocks...
        while not iter.IsDoneWithTraversal():
          curInput = iter.GetCurrentDataObject()
          curNumberPoints = curInput.GetNumberOfPoints()
          fem_point_data= curInput.GetPointData() 
          DeepCopy = 1
          #print timeID,listSliceInit,curNumberPoints 
          vtkMean= vtkNumPy.numpy_to_vtk( meandataList[listSliceInit:listSliceInit+curNumberPoints ] , DeepCopy) 
          vtkStdD= vtkNumPy.numpy_to_vtk( stdddataList[listSliceInit:listSliceInit+curNumberPoints ] , DeepCopy) 
          vtkSkew= vtkNumPy.numpy_to_vtk( skewdataList[listSliceInit:listSliceInit+curNumberPoints ] , DeepCopy) 
          vtkKurt= vtkNumPy.numpy_to_vtk( kurtdataList[listSliceInit:listSliceInit+curNumberPoints ] , DeepCopy) 
          vtkMean.SetName(       "Var%sMean"        % variable  )
          vtkStdD.SetName(       "Var%sStdDev"      % variable  )
          vtkSkew.SetName(       "Var%sSkew"        % variable  )
          vtkKurt.SetName(       "Var%sKurt"        % variable  )
          fem_point_data.AddArray( vtkMean )
          fem_point_data.AddArray( vtkStdD )
          fem_point_data.AddArray( vtkSkew )
          fem_point_data.AddArray( vtkKurt )
          fem_point_data.Update()
          for iii,response in enumerate(responsedataList):
            vtkResponse= vtkNumPy.numpy_to_vtk( response[listSliceInit:listSliceInit+curNumberPoints ] , DeepCopy) 
            vtkResponse.SetName(   "Var%sresponse%d" %(variable,iii) )
            fem_point_data.AddArray( vtkResponse )
            fem_point_data.Update()
          for iii,probability in enumerate(probabilitydataList):
            vtkProb= vtkNumPy.numpy_to_vtk( probability[listSliceInit:listSliceInit+curNumberPoints ] , DeepCopy) 
            vtkProb.SetName(   "Var%sprobability%d" %(variable,iii) )
            fem_point_data.AddArray( vtkProb )
            fem_point_data.Update()
          for iii,reliability in enumerate(reliabilitydataList):
            vtkReliab= vtkNumPy.numpy_to_vtk( reliability[listSliceInit:listSliceInit+curNumberPoints ] , DeepCopy) 
            vtkReliab.SetName(   "Var%sreliability%d" %(variable,iii) )
            fem_point_data.AddArray( vtkReliab )
            fem_point_data.Update()
          #vtkSoln = vtkNumPy.numpy_to_vtk( listSliceInit * numpy.ones(curNumberPoints ), DeepCopy ) 
          #vtkSoln.SetName("%d" % listSliceInit) 
          curInput.Update()
          # multiBlockData['%d' % timeID] = vtk.vtkUnstructuredGrid()
          # multiBlockData['%d' % timeID].DeepCopy( curInput )
          listSliceInit = listSliceInit + curNumberPoints 
          iter.GoToNextItem();
      # single block
      else:
        raise RuntimeError("not implemented yet... " )
      vtkExodusIIWriter = vtk.vtkExodusIIWriter()
      vtkExodusIIWriter.SetFileName( '%s/fem_stats.%04d.e' % (".",timeID) )
      vtkExodusIIWriter.SetInput( exodusObject)
      vtkExodusIIWriter.Update()
      #vtkTemporalDataSet.SetTimeStep( timeID, multiBlockData['%d' % timeID] )
      #vtkTemporalDataSet.Update()

  ## exodusObject.Update()
  ## vtkExodusIIWriter = vtk.vtkExodusIIWriter()
  ## vtkExodusIIWriter.SetFileName( '%s/fem_stats.e' % data_dir )
  ## vtkExodusIIWriter.WriteOutBlockIdArrayOn ();
  ## vtkExodusIIWriter.WriteOutGlobalNodeIdArrayOn ();
  ## vtkExodusIIWriter.WriteOutGlobalElementIdArrayOn ();
  ## vtkExodusIIWriter.WriteAllTimeStepsOn ();
  ## print "#time", vtkTemporalDataSet.GetNumberOfTimeSteps()
  ## vtkExodusIIWriter.SetInput( exodusObject )
  ## vtkExodusIIWriter.Update()

  ## # hold imaging
  ## momentSystem      = {} 
  ## responseSystem    = {} 
  ## probabilitySystem = {} 
  ## reliabilitySystem = {} 
  ## for variable,responseList in responseLevelVarList:
  ##   varID = 'Var%d' % variable
  ##   momentSystem[      varID ] = eqnSystems.AddExplicitSystem( "Var%dMoments"     % variable  ) 
  ##   responseSystem   [ varID ] = eqnSystems.AddExplicitSystem( "Var%dResponse"    % variable  )
  ##   probabilitySystem[ varID ] = eqnSystems.AddExplicitSystem( "Var%dProbability" % variable  )
  ##   reliabilitySystem[ varID ] = eqnSystems.AddExplicitSystem( "Var%dReliability" % variable  )
  ##   momentSystem[      varID ].AddFirstLagrangeVariable(       "Var%dMean"        % variable  )
  ##   momentSystem[      varID ].AddFirstLagrangeVariable(       "Var%dStdDev"      % variable  )
  ##   momentSystem[      varID ].AddFirstLagrangeVariable(       "Var%dSkew"        % variable  )
  ##   momentSystem[      varID ].AddFirstLagrangeVariable(       "Var%dKurt"        % variable  )
  ##   
  ## # initialize libMesh data structures
  ## eqnSystems.init( ) 
  ## StatsOutputFile = "%s/fem_stats.e" % data_dir
  ## # write IC
  ## exodusII_IO = femLibrary.PylibMeshExodusII_IO(femMesh)
  ## exodusII_IO.WriteTimeStep(StatsOutputFile,eqnSystems, 1, 0.0 )  
  ## 
  ## # loop over time steps and import data
  ##   # write time
  ##   print timeID
  ##   exodusII_IO.WriteTimeStep(StatsOutputFile,eqnSystems, timeID+1, timeID*fem_params['deltat'])  
####################################################################
def ExtractSolutionFromExodus(inputExodusFile,variableID,timeID,fileID):
  """
  extract data as a vector from exodus file
  """
  import vtk
  import vtk.util.numpy_support as vtkNumPy 
  import numpy
  vtkExodusIIReader = vtk.vtkExodusIIReader()
  print "opening %s " % inputExodusFile 
  vtkExodusIIReader.SetFileName( inputExodusFile )
  vtkExodusIIReader.Update()
  ntime  = vtkExodusIIReader.GetNumberOfTimeSteps()
  print "ntime %d varID %s time %d " % (ntime,variableID,timeID)
  vtkExodusIIReader.SetTimeStep(timeID) 
  vtkExodusIIReader.SetPointResultArrayStatus(variableID,1)
  vtkExodusIIReader.Update()
  # multi block
  if vtkExodusIIReader.GetOutput().IsA("vtkMultiBlockDataSet"):
    iter = vtkExodusIIReader.GetOutput().NewIterator()
    iter.UnRegister(None)
    iter.InitTraversal()
    # initialize list for storage
    Soln = []
    # loop over blocks...
    while not iter.IsDoneWithTraversal():
      curInput = iter.GetCurrentDataObject()
      fem_point_data= curInput.GetPointData() 
      Soln.append( vtkNumPy.vtk_to_numpy(fem_point_data.GetArray(variableID)) ) 
      iter.GoToNextItem();
    # concatenate blocks
    numpy.savetxt( ("data/"+SolnOutputTemplate+".%d") % (variableID,timeID,fileID),numpy.concatenate(Soln))
  # single block
  else:
    curInput = vtkExodusIIReader.GetOutput()
    fem_point_data= curInput.GetPointData() 
    Soln = vtkNumPy.vtk_to_numpy(fem_point_data.GetArray(variableID)) 
    numpy.savetxt( ("data/"+SolnOutputTemplate+".%d") % (variableID,timeID,fileID),Soln)
# end def ExtractSolutionFromExodus(**kwargs):
####################################################################
def pennesModeling(**kwargs):
  """
  treatment planning model 
  """
  # import petsc and numpy
  import petsc4py, numpy
  # init petsc
  PetscOptions =  sys.argv
  PetscOptions.append("-ksp_monitor")
  PetscOptions.append("-ksp_rtol")
  PetscOptions.append("1.0e-15")
  #PetscOptions.append("-help")
  #PetscOptions.append("-idb")
  petsc4py.init(PetscOptions)
  
  # break processors into separate communicators
  from petsc4py import PETSc
  petscRank = PETSc.COMM_WORLD.getRank()
  petscSize = PETSc.COMM_WORLD.Get_size()
  sys.stdout.write("petsc rank %d petsc nproc %d\n" % (petscRank, petscSize))

  # the original configuration ini file should be stored
  config = kwargs['config_parser']

  # set shell context
  # TODO import vtk should be called after femLibrary ???? 
  # FIXME WHY IS THIS????
  import femLibrary
  # initialize libMesh data structures
  libMeshInit = femLibrary.PyLibMeshInit(PetscOptions,PETSc.COMM_WORLD) 
  
  # store control variables
  getpot = femLibrary.PylibMeshGetPot(PetscOptions) 

  # copy all values from the input file
  for section in config.sections():
    for name,value in  config.items(section):
      #print "%s/%s" % (section,name) , value
      getpot.SetIniValue( "%s/%s" % (section,name) , value ) 

  #####################################################
  # update values from dakota
  #####################################################
  # thermal conductivity from Duck/CRC Handbook
  try:
    getpot.SetIniValue( "thermal_conductivity/k_0_healthy",
                              kwargs['cv']['k_0_healthy'] ) 
  except KeyError:
    pass #use default value
  # nonlinear conductivity
  try:
    getpot.SetIniValue( "thermal_conductivity/k_1",
                                kwargs['cv']['k_1'] ) 
  except KeyError:
    pass #use default value
  # tumor conductivity
  try:
    getpot.SetIniValue( "thermal_conductivity/k_0_tumor",
                                kwargs['cv']['k_0_tumor'] ) 
  except KeyError:
    pass #use default value
  # perfusion from Duck/CRC Handbook
  try:
    getpot.SetIniValue( "perfusion/w_0_healthy",
                              kwargs['cv']['w_0_healthy'] ) 
  except KeyError:
    pass #use default value
  try:
    getpot.SetIniValue( "perfusion/w_1",
                     kwargs['cv']['w_1_coag'] ) 
  except KeyError:
    pass #use default value
  try:
    getpot.SetIniValue( "perfusion/w_0_tumor",
                     kwargs['cv']['w_0_tumor'] ) 
  except KeyError:
    pass #use default value
  # water properties from Duck
  # Adult white mater
  # 3.2 1/cm
  # Beek et al., 1993a
  # Adult grey mater
  # 5.0 1/cm
  try:
    getpot.SetIniValue( "optical/mu_a_healthy",
                   kwargs['cv']['mu_a_healthy'] ) 
  except KeyError:
    pass #use default value
  try:
    getpot.SetIniValue( "optical/mu_a_1",
                   kwargs['cv']['mu_a_coag'] ) 
  except KeyError:
    pass #use default value
  # 1-300
  try:
    getpot.SetIniValue( "optical/mu_a_tumor",
                   kwargs['cv']['mu_a_tumor'] ) 
  except KeyError:
    pass #use default value
  # FIXME large mu_s (> 30) in agar causing negative fluence to satisfy BC 
  try:
    getpot.SetIniValue( "optical/mu_s_healthy",
                   kwargs['cv']['mu_s_healthy'] ) 
  except KeyError:
    pass #use default value
  try:
    getpot.SetIniValue( "optical/mu_s_1",
                   kwargs['cv']['mu_s_coag'] ) 
  except KeyError:
    pass #use default value
  # 1-300
  try:
    getpot.SetIniValue( "optical/mu_s_tumor",
                   kwargs['cv']['mu_s_tumor'] ) 
  except KeyError:
    pass #use default value

  # from AE paper
  #http://scitation.aip.org/journals/doc/MPHYA6-ft/vol_36/iss_4/1351_1.html#F3
  # .9  - .99
  try:
    getpot.SetIniValue( "optical/anfact",
                   kwargs['cv']['anfact'] ) 
  except KeyError:
    pass #use default value
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
  
  # initialize FEM Mesh
  femMesh = femLibrary.PylibMeshMesh()
  #femMesh.SetupUnStructuredGrid(kwargs['mesh_file'],0,RotationMatrix, Translation  ) 
  femMesh.ReadFile(kwargs['mesh_file'])
  MeshOutputFile = MeshOutputTemplate % kwargs['fileID'] 
  #fem.SetupStructuredGrid( (10,10,4) ,[0.0,1.0],[0.0,1.0],[0.0,1.0]) 
  
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
  
  # initialize libMesh data structures
  eqnSystems.init( ) 
  
  # quick error check
  #errCheckSoln = fem.GetSolutionVector( "StateSystem" )[...]
  #if (  errCheckSoln.size + 1 != kwargs['functions'] ):
  #  print "ERROR!! number of response functions incorrect!!"
  #  raise RuntimeError("soln vec + 1 = %d .NE. num_response = %d"%(errCheckSoln.size+1,kwargs['functions']) )

  # print info
  eqnSystems.PrintSelf() 
  
  # setup IC 
  pennesSystem.PetscFEMSystemSetupInitialConditions( ) 

  # get responses for list of variable
  responseLevelVarList = kwargs['responseLevelVarList'] 

  # write IC
  exodusII_IO = femLibrary.PylibMeshExodusII_IO(femMesh)
  exodusII_IO.WriteTimeStep(MeshOutputFile,eqnSystems, 1, 0.0 )  
  
  ObjectiveFunction = 0.0
  # loop over time steps and solve
  print  "deltat = ",kwargs['deltat']
  #for timeID in range(1,2):
  for timeID in range(1,kwargs['ntime']+1):
     print "time step = " ,timeID
     #for subTimeID in range(1):
     for subTimeID in range(kwargs['nsubstep']):
       pennesSystem.PetscFEMSystemUpdateTimeStep( timeID ) 
       pennesSystem.SystemSolve( ) 
     #fem.StoreTransientSystemTimeStep("StateSystem",timeID ) 
     exodusII_IO.WriteTimeStep(MeshOutputFile ,eqnSystems, timeID+1, timeID*kwargs['acquisitionTime'])  
  retval = dict([])
  retval['fns'] = [ObjectiveFunction]
  retval['rank'] = petscRank 
  return(retval)
# end def pennesModeling(**kwargs):
####################################################################
def WriteWorstCastDakotaInputFile(filename,analysis_driver,
                         parameters_file,results_file,dakota_work_directory,
                         allow_existing_results,num_response,gPCVariables):
  dakotafile=open( filename,"w")
  # TODO better to pass None and capture ?
  if(allow_existing_results == None):
    allow_existing_results  = ""
  dakotafile.write(
"""
# see reference manual
#   http://dakota.sandia.gov/licensing/stable/html-ref/index.html

# run single method and output to a file
#  http://dakota.sandia.gov/licensing/stable/html-ref/StratCommands.html
strategy,
	single_method 

interface,
	system 
          asynchronous			#0,#p0
          evaluation_concurrency 16 	#3,#8,#18,#19
	  analysis_driver = '%s'
	  #analysis_driver = './ibrun_par_driver'
	  # this will guarantee that evaluations are replaced with
          # evaluations modulo the evaluation concurrency
 	  local_evaluation_static_scheduling
	  parameters_file = '%s'
	  results_file = '%s'
          work_directory named = "%s"
	  file_save file_tag directory_save 
	  deactivate restart_file
          %s

responses,
	num_response_functions = %d
	no_gradients
	no_hessians

method,
	multidim_parameter_study		
	  output verbose
""" % (analysis_driver,parameters_file,results_file,dakota_work_directory,allow_existing_results,num_response) )
  dakotafile.write( "          partitions    = " )
  for key,var in gPCVariables.iteritems():
     dakotafile.write( " 1  " )
  dakotafile.write( "\n")
  dakotafile.write( "variables,\n")
  dakotafile.write( "        discrete_design_set_real = %d\n" % len(gPCVariables) )
  dakotafile.write( "          set_values  = " )
  for key,var in gPCVariables.iteritems():
     dakotafile.write( " %f %f " % (var['lower_bound'],var['upper_bound']) )
  dakotafile.write( "\n")
  dakotafile.write( "          descriptors  = " )
  for key,var in gPCVariables.iteritems():
     dakotafile.write( " '%s' " % key )
  dakotafile.write( "\n")
  dakotafile.close; dakotafile.flush() 
# end def WriteWorstCastDakotaInputFile
####################################################################
def WritegPCDakotaInputFile(filename,evaluation_concurrency,analysis_driver,
                         parameters_file,results_file,dakota_work_directory,
                         file_save,directory_save,
                         optional_interface,num_response,gPCVariables,
                         responseLevelList,probabilityLevelList,reliabilityLevelList):
  dakotafile=open( filename,"w")
  # TODO better to pass None and capture ?
  if(optional_interface == None):
    optional_interface = ""
  if(file_save == None):
    file_save = ""
  if(directory_save == None):
    directory_save = ""
  dakotafile.write(
"""
# see reference manual
#   http://dakota.sandia.gov/licensing/stable/html-ref/index.html

# run single method and output to a file
#  http://dakota.sandia.gov/licensing/stable/html-ref/StratCommands.html
strategy,
	single_method 

interface,
	system 
          asynchronous			#0,#p0
          evaluation_concurrency %d 	#3,#8,#18,#19
	  analysis_driver = '%s'
	  #analysis_driver = './ibrun_par_driver'
	  # this will guarantee that evaluations are replaced with
          # evaluations modulo the evaluation concurrency
 	  local_evaluation_static_scheduling
	  parameters_file = '%s'
	  results_file = '%s'
          work_directory named = "%s"
	  file_tag   %s %s 
	  deactivate restart_file
          %s

responses,
	num_response_functions = %d
	no_gradients
	no_hessians

method,
	polynomial_chaos
	  samples = 10000		
	  seed = 12347 rng rnum2	
          # vector response input 
          # http://dakota.sandia.gov/licensing/stable/html-ref/IntroCommands.html#IntroCmdsInpSpec
	  variance_based_decomp #univariate_effects
	  #output verbose
	  output silent
""" % (evaluation_concurrency,analysis_driver,parameters_file,results_file,dakota_work_directory,file_save,directory_save, optional_interface,num_response) )
  dakotafile.write( "          response_levels    = " )
  for idresp in range(num_response):
    for response in responseLevelList:
       dakotafile.write( " %f " % response )
  dakotafile.write( "\n")
  dakotafile.write( "          probability_levels = " )
  for idresp in range(num_response):
    for probability in probabilityLevelList:
       dakotafile.write( " %f " % probability )
  dakotafile.write( "\n")
  dakotafile.write( "          reliability_levels = " )
  for idresp in range(num_response):
    for reliability in reliabilityLevelList:
       dakotafile.write( " %f " % reliability )
  dakotafile.write( "\n")
  dakotafile.write( "          quadrature_order  = ")
  for key,var in gPCVariables.iteritems():
     dakotafile.write( " %d " % var['quadrature_order'] )
  dakotafile.write( "\n")
  dakotafile.write( "variables,\n")
  # build dictionary of uniform variables
  uniformVariables   = dict((key,var) for key,var in gPCVariables.iteritems() if var['type'] == 'uniform_uncertain')
  logUniformVariables= dict((key,var) for key,var in gPCVariables.iteritems() if var['type'] == 'loguniform_uncertain')
  normalVariables    = dict((key,var) for key,var in gPCVariables.iteritems() if var['type'] == 'normal_uncertain' )
  logNormalVariables = dict((key,var) for key,var in gPCVariables.iteritems() if var['type'] == 'lognormal_uncertain' )
  # write uniform variables
  if len(uniformVariables) :
    try: 
      dakotafile.write( "        uniform_uncertain = %d\n" % len(uniformVariables) )
      dakotafile.write( "          lower_bounds  = " )
      for key,var in uniformVariables.iteritems():
         dakotafile.write( " %f " % var['lower_bound'] )
      dakotafile.write( "\n")
      dakotafile.write( "          upper_bounds  = " )
      for key,var in uniformVariables.iteritems():
         dakotafile.write( " %f " % var['upper_bound'] )
      dakotafile.write( "\n")
      dakotafile.write( "          descriptors  = " )
      for key,var in uniformVariables.iteritems():
         dakotafile.write( " '%s' " % key )
      dakotafile.write( "\n")
    except KeyError:
      print "\n\n error inputing uniform variable"
      raise RuntimeError("\n\n uniform variables should be of the form {'type':'uniform_uncertain', 'quadrature_order':4,'lower_bound':3.,'upper_bound':9.} ")
  # write log uniform variables
  if len(logUniformVariables) :
    try: 
      dakotafile.write( "        loguniform_uncertain = %d\n" % len(logUniformVariables) )
      dakotafile.write( "          lower_bounds  = " )
      for key,var in logUniformVariables.iteritems():
         dakotafile.write( " %f " % var['lower_bound'] )
      dakotafile.write( "\n")
      dakotafile.write( "          upper_bounds  = " )
      for key,var in logUniformVariables.iteritems():
         dakotafile.write( " %f " % var['upper_bound'] )
      dakotafile.write( "\n")
      dakotafile.write( "          descriptors  = " )
      for key,var in logUniformVariables.iteritems():
         dakotafile.write( " '%s' " % key )
      dakotafile.write( "\n")
    except KeyError:
      print "\n\n error inputing uniform variable"
      raise RuntimeError("\n\n uniform variables should be of the form {'type':'uniform_uncertain', 'quadrature_order':4,'lower_bound':3.,'upper_bound':9.} ")
  # write normal variables
  if len(normalVariables) :
    try:
      dakotafile.write( "        normal_uncertain  = %d\n" % len(normalVariables) )
      dakotafile.write( "          means  = " )
      for key,var in normalVariables.iteritems():
         dakotafile.write( " %f " % var['mean'] )
      dakotafile.write( "\n")
      dakotafile.write( "          std_deviations  = " )
      for key,var in normalVariables.iteritems():
         dakotafile.write( " %f " % var['std_dev'] )
      dakotafile.write( "\n")
      dakotafile.write( "          descriptors  = " )
      for key,var in normalVariables.iteritems():
         dakotafile.write( " '%s' " % key )
      dakotafile.write( "\n")
    except KeyError:
      print "\n\n error inputing normal variable"
      raise RuntimeError("\n\n normal variables should be of the form {'type':'normal_uncertain', 'quadrature_order':4,'mean':4.,'std_dev':2.,'lower_bound':3.,'upper_bound':9.} ")
  # write log normal variables
  if len(logNormalVariables) :
    try: 
      dakotafile.write( "        lognormal_uncertain  = %d\n" % len(logNormalVariables) )
      dakotafile.write( "          means  = " )
      for key,var in logNormalVariables.iteritems():
         dakotafile.write( " %f " % var['mean'] )
      dakotafile.write( "\n")
      dakotafile.write( "          std_deviations  = " )
      for key,var in logNormalVariables.iteritems():
         dakotafile.write( " %f " % var['std_dev'] )
      dakotafile.write( "\n")
      dakotafile.write( "          descriptors  = " )
      for key,var in logNormalVariables.iteritems():
         dakotafile.write( " '%s' " % key )
      dakotafile.write( "\n")
    except KeyError:
      print "\n\n error inputing log normal variable"
      raise RuntimeError("\n\n lognormal variables should be of the form {'type':'lognormal_uncertain', 'quadrature_order':4,'mean':4.,'std_dev':2.,'lower_bound':3.,'upper_bound':9.} ")
  # close the file
  dakotafile.close; dakotafile.flush() 
# end def WritegPCDakotaInputFile
##################################################################
def ParseInput(param_file):
  # ----------------------------
  # Parse DAKOTA parameters file
  # ----------------------------
  
  # setup regular expressions for parameter/label matching
  e = '-?(?:\\d+\\.?\\d*|\\.\\d+)[eEdD](?:\\+|-)?\\d+' # exponential notation
  f = '-?\\d+\\.\\d*|-?\\.\\d+'                        # floating point
  i = '-?\\d+'                                         # integer
  value = e+'|'+f+'|'+i                                # numeric field
  tag = '\\w+(?::\\w+)*'                               # text tag field
  
  # regular expression for aprepro parameters format
  aprepro_regex = re.compile('^\s*\{\s*(' + tag + ')\s*=\s*(' + value +')\s*\}$')
  # regular expression for standard parameters format
  standard_regex = re.compile('^\s*(' + value +')\s+(' + tag + ')$')
  
  # open DAKOTA parameters file for reading
  paramsfile = open(param_file, 'r')
  fileID = int(param_file.split(".").pop())
  #fileID = int(os.getcwd().split(".").pop())
  
  # extract the parameters from the file and store in a dictionary
  paramsdict = {}
  for line in paramsfile:
      m = aprepro_regex.match(line)
      if m:
          paramsdict[m.group(1)] = m.group(2)
      else:
          m = standard_regex.match(line)
          if m:
              paramsdict[m.group(2)] = m.group(1)
  
  paramsfile.close()
  
  # crude error checking; handle both standard and aprepro cases
  num_vars = 0
  if ('variables' in paramsdict):
      num_vars = int(paramsdict['variables'])
  elif ('DAKOTA_VARS' in paramsdict):
      num_vars = int(paramsdict['DAKOTA_VARS'])
  
  num_fns = 0
  if ('functions' in paramsdict):
      num_fns = int(paramsdict['functions'])
  elif ('DAKOTA_FNS' in paramsdict):
      num_fns = int(paramsdict['DAKOTA_FNS'])
  
  # set a dictionary for passing to fem via Python kwargs
  # read in params previously stored in dictionary and written
  pkl_file = open('../CaseInfo.pkl', 'rb')
  fem_params = pickle.load(pkl_file)
  pkl_file.close()

  # -------------------------------
  # Convert and send to application
  # -------------------------------
  
  # set up the data structures the rosenbrock analysis code expects
  # for this simple example, put all the variables into a single hardwired array
  continuous_vars = {} 

  try:
    continuous_vars['k_0_healthy' ] = paramsdict['k_0_healthy' ]
  except KeyError:
    pass

  try:
    continuous_vars['k_1'         ] = paramsdict['k_1'         ]
  except KeyError:
    pass

  try:
    continuous_vars['k_0_tumor'   ] = paramsdict['k_0_tumor'   ]
  except KeyError:
    pass

  try:
    continuous_vars['mu_a_healthy'] = paramsdict['mu_a_healthy']
  except KeyError:
    pass

  try:
    continuous_vars['mu_a_coag'  ]  = paramsdict['mu_a_coag'  ]
  except KeyError:
    pass

  try:
    continuous_vars['mu_a_tumor'  ] = paramsdict['mu_a_tumor'  ]
  except KeyError:
    pass

  try:
    continuous_vars['mu_s_healthy'] = paramsdict['mu_s_healthy']
  except KeyError:
    pass

  try:
    continuous_vars['mu_s_coag'  ] = paramsdict['mu_s_coag'  ]
  except KeyError:
    pass

  try:
    continuous_vars['mu_s_tumor'  ] = paramsdict['mu_s_tumor'  ]
  except KeyError:
    pass

  try:
    continuous_vars['w_0_healthy'] = paramsdict['w_0_healthy' ]  
  except KeyError:
    pass

  try:
    continuous_vars['w_1_coag'  ] = paramsdict['w_1_coag'   ] 
  except KeyError:
    pass
  
  try:
    continuous_vars['w_0_tumor'  ] = paramsdict['w_0_tumor'   ] 
  except KeyError:
    pass
  
  try:
    continuous_vars['anfact'] = paramsdict['anfact'   ] 
  except KeyError:
    pass
  
  try:
    active_set_vector = [ int(paramsdict['ASV_%d:response_fn_%d' % (i,i) ]) for i in range(1,num_fns+1)  ] 
  except KeyError:
    active_set_vector = [ int(paramsdict['ASV_%d:obj_fn' % (i) ]) for i in range(1,num_fns+1)  ] 
  
  # store dakota vars
  fem_params['cv']        = continuous_vars
  fem_params['asv']       = active_set_vector
  fem_params['functions'] = num_fns
  fem_params['fileID']    = fileID 

  return fem_params
  ## ----------------------------
  ## Return the results to DAKOTA
  ## ----------------------------
  #
  #if (fem_results['rank'] == 0 ):
  #  # write the results.out file for return to DAKOTA
  #  # this example only has a single function, so make some assumptions;
  #  # not processing DVV
  #  outfile = open('results.out.tmp.%d' % fileID, 'w')
  #  
  #  # write functions
  #  for func_ind in range(0, num_fns):
  #      if (active_set_vector[func_ind] & 1):
  #          functions = fem_results['fns']    
  #          outfile.write(str(functions[func_ind]) + ' f' + str(func_ind) + '\n')
  #  
  #  ## write gradients
  #  #for func_ind in range(0, num_fns):
  #  #    if (active_set_vector[func_ind] & 2):
  #  #        grad = rosen_results['fnGrads'][func_ind]
  #  #        outfile.write('[ ')
  #  #        for deriv in grad: 
  #  #            outfile.write(str(deriv) + ' ')
  #  #        outfile.write(']\n')
  #  #
  #  ## write Hessians
  #  #for func_ind in range(0, num_fns):
  #  #    if (active_set_vector[func_ind] & 4):
  #  #        hessian = rosen_results['fnHessians'][func_ind]
  #  #        outfile.write('[[ ')
  #  #        for hessrow in hessian:
  #  #            for hesscol in hessrow:
  #  #                outfile.write(str(hesscol) + ' ')
  #  #            outfile.write('\n')
  #  #        outfile.write(']]')
  #  #
  #  outfile.close();outfile.flush
  #  #
  #  ## move the temporary results file to the one DAKOTA expects
  #  #import shutil
  #  #shutil.move('results.out.tmp.%d' % fileID, sys.argv[2])
# end def ParseInput:
##################################################################

# setup command line parser to control execution
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--dakota_exe", 
                  action="store", dest="dakota_exe", default="dakota",
                  help="full path to dakota EXE", metavar="EXE")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
parser.add_option( "--pre_run", "--setup",
                  action="store", dest="pre_run", default=None,
                  help="path to driver setup FILE=[mpich2_setup_driver,ranger_setup_driver]", metavar="FILE")
parser.add_option( "--run_fem","--param_file", 
                  action="store", dest="param_file", default=None,
                  help="run code with parameter FILE", metavar="FILE")
parser.add_option( "--extract_exo",
                  action="store", dest="extract_exo", default=None,
                  help="extract data from exodus file for run code with parameter FILE", metavar="FILE")
parser.add_option( "--run_queue",
                  action="store", dest="run_queue", default=None,
                  help="run jobs in DIR", metavar="DIR")
parser.add_option( "--execution",
                  action="store", dest="execution", default="sh",
                  help="run jobs with EXE to run exe.qsub", metavar="EXE")
parser.add_option( "--assemble_stats",
                  action="store", dest="assemble_stats", default=None,
                  help="assemble stats in DIR to fem file", metavar = "DIR")
parser.add_option( "--config_ini",
                  action="store", dest="config_ini", default=None,
                  help="select input config file", metavar = "FILE")
parser.add_option("--k_0_healthy",
                  action="store_true", dest="k_0_healthy", default=False,
                  help="gPC for this variable")
parser.add_option("--k_0_tumor",
                  action="store_true", dest="k_0_tumor", default=False,
                  help="gPC for this variable")
parser.add_option("--s_0_healthy",
                  action="store_true", dest="s_0_healthy", default=False,
                  help="gPC for this variable")
parser.add_option("--s_0_tumor",
                  action="store_true", dest="s_0_tumor", default=False,
                  help="gPC for this variable")
parser.add_option("--k_1",
                  action="store_true", dest="k_1", default=False,
                  help="gPC for this variable")
parser.add_option("--w_0_healthy",
                  action="store_true", dest="w_0_healthy", default=False,
                  help="gPC for this variable")
parser.add_option("--w_0_tumor",
                  action="store_true", dest="w_0_tumor", default=False,
                  help="gPC for this variable")
parser.add_option("--w_1_coag",
                  action="store_true", dest="w_1_coag", default=False,
                  help="gPC for this variable")
parser.add_option("--mu_a_healthy",
                  action="store_true", dest="mu_a_healthy", default=False,
                  help="gPC for this variable")
parser.add_option("--mu_a_tumor",
                  action="store_true", dest="mu_a_tumor", default=False,
                  help="gPC for this variable")
parser.add_option("--mu_a_coag",
                  action="store_true", dest="mu_a_coag", default=False,
                  help="gPC for this variable")
parser.add_option("--mu_s_healthy",
                  action="store_true", dest="mu_s_healthy", default=False,
                  help="gPC for this variable")
parser.add_option("--mu_s_tumor",
                  action="store_true", dest="mu_s_tumor", default=False,
                  help="gPC for this variable")
parser.add_option("--mu_s_coag",
                  action="store_true", dest="mu_s_coag", default=False,
                  help="gPC for this variable")
(options, args) = parser.parse_args()

# write initial dakota input scripts to setup directories and files
if (options.pre_run != None):

  # read config file
  config = ConfigParser.SafeConfigParser({
                                          "nsubstep":'1'
                                         })
  config.read(options.config_ini)
  config.set("timestep","subtimestep","%(nsubstep)s") 

  # working directory
  work_dir = config.get('exec','work_dir')

  # mesh_file options should be of the form mesh_files = [mesh1.e,mesh2.e,mesh3.e,...]
  MeshFileList = config.get('pre_run','mesh_files')[1:-1].split(',')
  for MeshFile in MeshFileList: 
     # setup any data set cases and store in the directory
     DataSetParams              = {}
     DataSetParams['mesh_file'] = "%s/%s" % (os.getcwd(),MeshFile)

     # check that mesh is available and get # nodes
     if( not os.path.isfile(DataSetParams['mesh_file']) ):
      raise RuntimeError("mesh %s not found !!!!" % DataSetParams['mesh_file'] )

     DataSetParams['acquisitionTime'] = config.getfloat('timestep','acquisitiontime')
     DataSetParams['nsubstep'] = int(config.get('timestep','subtimestep'))
     DataSetParams['powerHistory'] = eval(config.get('timestep','powerhistory'))
     DataSetParams['default_perfusion'] = config.get('perfusion','w_0_healthy')
     DataSetParams['physics'] = config.get('exec','physics')
     #get last time point as ntime
     DataSetParams['ntime'] = DataSetParams['powerHistory'][0][-1]
     try:
       if( config.get('steadystate','domain_0') ):
         DataSetParams['ntime'] = 1
     except:
       pass
     DataSetParams['deltat'] = DataSetParams['acquisitionTime'] / DataSetParams['nsubstep']
     print "reading ini %s,nsubstep=%d ,ntime=%d" % (options.config_ini,DataSetParams['nsubstep'],DataSetParams['ntime'])

     # create work dir
     numjobid = 0
     workID   = "%s_%d" %(work_dir,numjobid)
     if(os.path.isdir(workID)):
       print "%s ALREADY EXISTS!!!!" % workID
       print ""
       while( os.path.isdir("%s_%d" %(work_dir,numjobid)) ):
          numjobid = numjobid + 1
       workID= "%s_%d" %(work_dir,numjobid)
     os.mkdir(workID)
     
     # get probability levels
     responseLevelVarList = eval(config.get('dakota','responselevels'   ))
     probabilityLevelList = eval(config.get('dakota','probabilitylevels'))
     reliabilityLevelList = eval(config.get('dakota','reliabilitylevels'))
     DataSetParams['responseLevelVarList'] = responseLevelVarList   
     DataSetParams['probabilityLevelList'] = probabilityLevelList
     DataSetParams['reliabilityLevelList'] = reliabilityLevelList

     # store the entire configuration file for convienence
     DataSetParams['config_parser'] = config

     # store data set params in the director
     pickleFile = open("%s/CaseInfo.pkl" % workID,'wb')
     pickle.dump(DataSetParams,pickleFile)
     pickleFile.close()

     #build list of gPC Vars
     # expecting an input of the form of a dictionary for each variable of interest, ie
     # w_0_healthy_gpc = {'type':'uniform_uncertain', 'quadrature_order':4,'lower_bound':3.,'upper_bound':9.}
     gPCVariables = {}
     if (    options.k_0_healthy ):
       gPCVariables['k_0_healthy' ]  =  eval(config.get('thermal_conductivity','k_0_healthy_gpc' ))
     if (    options.k_0_tumor ):
       gPCVariables['k_0_tumor' ]    =  eval(config.get('thermal_conductivity','k_0_tumor_gpc' ))
     if (    options.k_1 ):
       gPCVariables['k_1' ]          =  eval(config.get('thermal_conductivity','k_1_gpc' ))
     # perfusion parameters
     # TODO need to setup default values
     if (    options.w_0_healthy):
       gPCVariables['w_0_healthy' ]  = eval(config.get('perfusion','w_0_healthy_gpc' ))
     if (    options.w_0_tumor ):
       gPCVariables['w_0_tumor' ]    = eval(config.get('perfusion','w_0_tumor_gpc'   ))
     if (    options.w_1_coag ):
       gPCVariables['w_1_coag' ]     = eval(config.get('perfusion','w_1_coag_gpc'   ))
     #@article{yaroslavsky2002optical,
     #title={Optical properties of selected native and coagulated human brain tissues in vitro in the visible and near infrared spectral range},
     #author={Yaroslavsky, AN and Schulze, PC and Yaroslavsky, IV and Schober, R. and Ulrich, F. and Schwarzmaier, HJ},
     #journal={Physics in medicine and biology},
     #volume={47},
     #pages={2059},
     #year={2002},
     #publisher={IOP Publishing} }
     if (    options.mu_a_healthy ):
       gPCVariables['mu_a_healthy' ] = eval(config.get('optical','mu_a_healthy_gpc'   ))
     if (    options.mu_a_tumor ):
       gPCVariables['mu_a_tumor' ]   = eval(config.get('optical','mu_a_tumor_gpc'   ))
     if (    options.mu_a_coag ):
       gPCVariables['mu_a_coag' ]    = eval(config.get('optical','mu_a_coag_gpc'   ))
     if (    options.mu_s_healthy ):
       gPCVariables['mu_s_healthy' ] = eval(config.get('optical','mu_s_healthy_gpc'   ))
     if (    options.mu_s_tumor ):
       gPCVariables['mu_s_tumor' ]   = eval(config.get('optical','mu_s_tumor_gpc'   ))
     if (    options.mu_s_coag ):
       gPCVariables['mu_s_coag' ]    = eval(config.get('optical','mu_s_coag_gpc'   ))
     # RF conductivity parameters
     if (    options.s_0_healthy ):
       gPCVariables['s_0_healthy' ]  = eval(config.get('electric_conductivity','s_0_healthy_gpc'   ))
     if (    options.s_0_tumor ):
       gPCVariables['s_0_tumor' ]    = eval(config.get('electric_conductivity','s_0_tumor_gpc'   ))
     #error check
     if(len(gPCVariables) == 0 ):
       raise RuntimeError("\n\n no gPCVariables input ")

     # all in/out files should be different to avoid potential race condition
     WritegPCDakotaInputFile("%s/pce_setup_gPC.in" % workID,16,
                          "%s/%s" % (os.getcwd(),options.pre_run) , "pce.in", "pce.out" ,
                          "realization", "file_save","directory_save",None,1,gPCVariables,
                          [0.0,1.0],probabilityLevelList,reliabilityLevelList)

     # all in/out files should be different to avoid potential race condition
     WriteWorstCastDakotaInputFile("%s/pce_setup_worst.in" % workID,
                          "%s/%s" % (os.getcwd(),options.pre_run) , "worst.in", "worst.out" ,
                          "worst", None,1,gPCVariables)

     # write a second for actually running realizations
     WritegPCDakotaInputFile("%s/pce_run_gPC.in" % workID,16,
                          "%s/%s" % (os.getcwd(),"ibrun_par_driver") , "pce.in", "pce.out" ,
                          "realization","file_save","directory_save", None,1,gPCVariables,
                          [0.0,1.0],probabilityLevelList,reliabilityLevelList)
     # get # nodes
     num_func = GetMeshNodes( DataSetParams['mesh_file'] )
     # set the scratch directory
     scratchDir = config.get('exec','scratch_dir' )
     # write script for post processing stats
     postRunStatFile=open( "%s/paramlist" % workID ,"w")
     # create dir for storing all time files
     os.mkdir( "%s/timefiles" % (workID) )
     for idtime in range(0,DataSetParams['ntime']+1):
        for variable,responseList in responseLevelVarList:
           dakotaTimeFile = "%s/%s/timefiles/time.%s.%04d.in" % (os.getcwd(),workID,variable,idtime)
           fullpathJob =  "%s/%stime.%04d/%s" %(scratchDir,workID,idtime,variable)
           # all in/out files should be different to avoid potential race condition
           WritegPCDakotaInputFile(dakotaTimeFile ,8,
                               "echo" , "pce.%04d.in"%idtime ,SolnOutputTemplate %(variable,idtime),
                               fullpathJob+"/data",None,None,"allow_existing_results template_directory '%s/%s/realization'" %(os.getcwd(),workID),
                               num_func,gPCVariables, responseList ,probabilityLevelList,reliabilityLevelList)
           postRunStatFile.write("mkdir -p %s/data; cd %s; python $DDDAS_SRC/Examples/TreatmentPlanning/gPCWFS.py --extract_exo=%s/%s/realization;%s %s > /dev/null;" % ( fullpathJob,fullpathJob,os.getcwd(),workID,options.dakota_exe,dakotaTimeFile ))
        postRunStatFile.write( "cd %s/%s; python $DDDAS_SRC/Examples/TreatmentPlanning/gPCWFS.py --assemble_stats=%s/%stime.%04d;" % (os.getcwd(),workID,scratchDir,workID,idtime) )
        postRunStatFile.write( "rm -rf %s/%stime.%04d; " % ( scratchDir,workID,idtime ) )
        postRunStatFile.write("\n" )
     postRunStatFile.close; postRunStatFile.flush() 
     # set it up...
     os.system("cd %s; %s pce_setup_gPC.in" % (workID,options.dakota_exe) )
     os.system("cd %s; %s pce_setup_worst.in" % (workID,options.dakota_exe) )
#run all jobs in the queue
elif (options.run_queue != None):
  CODEEXEC=[]
  for job in os.listdir(options.run_queue):
    fullpathJob =  "%s/%s/%s" %(os.getcwd(),options.run_queue,job)
    # only consider directories
    if(os.path.isdir(fullpathJob)):
      # filter out directories not of interest
      if( job.find("realization") != -1 ) :
        CODEEXEC.append("cd %s; sleep 1; %s exe.qsub" % (fullpathJob,options.execution) )
  RunSubProcessQueue(CODEEXEC,"%s/error_run.log" % options.run_queue)
elif (options.param_file != None):
  fem_params = ParseInput(options.param_file)
  # execute the rosenbrock analysis as a separate Python module
  print "Running Pennes model..."
  fem_results = pennesModeling(**fem_params)
  print "Pennes Simulation complete."
elif (options.assemble_stats):
  AssembleStatistics(options.assemble_stats)
elif (options.extract_exo):
  directoryInfo = os.getcwd().split("/")
  variableID = directoryInfo.pop()
  timeID     = int(directoryInfo.pop().split(".").pop())
  nrealization = len(filter(lambda xfile: xfile.split(".")[1]=='in',os.listdir(options.extract_exo)))
  print variableID, nrealization,timeID
  for fileID in range(1,nrealization+1):
     inputExodusFile =  "/".join( [options.extract_exo] + [MeshOutputTemplate % fileID ] ) 
     #fem_params = ParseInput(inputFile)
     ExtractSolutionFromExodus(inputExodusFile,variableID,timeID,fileID)
else:
  parser.print_help()
  print options
