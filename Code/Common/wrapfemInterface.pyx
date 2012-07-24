include "wrapfemInterface.pxi"
# --------------------------------------------------------------------

cdef inline int CHKERR(int ierr) except -1:
    if ierr != 0: raise RuntimeError
    return 0  # no error

# -------------------- wrapper for LibMeshInit ----------------------- 
cdef class PyLibMeshInit:
    """Interface to libMesh::LibMeshInit class
        store a local pointer in cython
        to the the underlying C++ data structure
    """
    cdef LibMeshInit *thisptr  # hold a C++ instance which we're wrapping
    def __cinit__(self, object args,comm):
        cdef int    PyPetsc_Argc = 0       
        cdef char** PyPetsc_Argv = NULL   
        # default communicator to PETSC_COMM_SELF 
        # if comm is not set up properly...
        cdef MPI_Comm ccomm = GetComm(comm, PETSC_COMM_SELF)
        # get command line arguments                                    
        PetscGetArgs(&PyPetsc_Argc,&PyPetsc_Argv)                     
        # pop the error handler and stall for debugger if needed 
        CHKERR( PetscPopErrorHandlerAndDebug() )
        # initialize PETSc                                       
        self.thisptr = new LibMeshInit(PyPetsc_Argc, PyPetsc_Argv, ccomm)
# -------------------- wrapper for libMesh::Mesh --------------------- 
cdef class PylibMeshMesh:
    """Interface to libMesh::Mesh class
        store a local pointer in cython
        to the the underlying C++ data structure
    """
    cdef Mesh *thisptr  # hold a C++ instance which we're wrapping
    def __cinit__(self):
        cdef int    Dimension = 3       
        self.thisptr = new Mesh(Dimension);
    def SetupStructuredGrid(self,dimensions,xbounds,ybounds,zbounds,bcList):
        """"setup SetupStructuredGrid 
             input = dimensions,xmin,xmax,ymin,ymax,zmin,zmax
            dimensions[2] = 0 ==> 2D
            3D should default to HEX8  element type
            2D should default to QUAD4 element type
        """
        cdef vector[int] boundary
        for bc in bcList:
           boundary.push_back(bc)
        CHKERR(GenerateStructuredGrid( cython.operator.dereference(self.thisptr),
                       dimensions[0],dimensions[1],dimensions[2],
                         xbounds[0] ,  xbounds[1] ,
                         ybounds[0] ,  ybounds[1] ,
                         zbounds[0] ,  zbounds[1] , boundary ) )
    def ReadFile(self,PyStringMeshFile):
        cdef char* mesh_file= PyStringMeshFile
        # base class pointer should call the inherited type
        (<MeshBase *>self.thisptr).read( string(mesh_file) )
    def SetupUnStructuredGrid(self,PyStringMeshFile,Refinement,A,b):
        cdef char* mesh_file= PyStringMeshFile
        CHKERR(SetupUnStructuredGrid(self.thisptr,mesh_file,Refinement,
                                             A[0][0], A[0][1], A[0][2],
                                             A[1][0], A[1][1], A[1][2],
                                             A[2][0], A[2][1], A[2][2],
                                                b[0],    b[1],    b[2]))
    def PrintSelf(self):
        # base class pointer should call the inherited type
        (<MeshBase *>self.thisptr).print_info()
        CHKERR(FlushStdCoutCerr())
# -------------------- wrapper for libMesh::GetPot --------------------- 
cdef class PylibMeshGetPot:
    cdef GetPot *thisptr  # hold a C++ instance which we're wrapping
    def __cinit__(self,object args):
        cdef int    PyPetsc_Argc = 0       
        cdef char** PyPetsc_Argv = NULL   
        # get command line arguments                                    
        PetscGetArgs(&PyPetsc_Argc,&PyPetsc_Argv)                     
        self.thisptr = new GetPot(PyPetsc_Argc, PyPetsc_Argv)
    def SetIniValue(self,PyStringVarName,PyStringValue):
        cdef char* var_name = PyStringVarName
        cdef char* value    = PyStringValue
        self.thisptr.set_value(var_name,value)
    def SetIniPower(self,NsubStep,timePowerList):
        """
           the data structure is setup as follows and ASSUMES equispace time 
               distances, IDEAL_DT   \forall i
 
           Equally Divide the time steps by NsubStep
            
                                                     *NOTE* closed at beginning
           time = 0    ---------                        BUT open   at end
                           |                                |
                           |                               \|/
                        Power[0]    power between time = [0,1)  is Power[0]
                           |            
                           |
           time = 1    ---------
                           |
                           |          
                        Power[1]    power between time = [1,2)  is Power[1]   
                           |         
                           |
           time = 2    ---------
                           |
                           |        
                        Power[2]    power between time = [2,3)  is Power[2]
                           |       
                           |
           time = 3    ---------
                     .
                     .
                     .
                     .
        """
        # error check sorting
        assert timePowerList[0] == sorted(timePowerList[0])
        # variables to hold ini info
        cdef char* var_name 
        cdef char* value   
        # initialize and loop
        intervalID = 0
        for iBound in timePowerList[0]:
          while (intervalID < iBound):
            for istep in range(NsubStep):
              pyVarName= "power/power[%d]" %  (NsubStep * intervalID  + istep)
              var_name = pyVarName
              pyValue  = "%f" % timePowerList[1][timePowerList[0].index(iBound)]
              value    = pyValue  
              self.thisptr.set_value(var_name,value)
            intervalID = intervalID + 1
        # set size info
        pyValue  = "%d" % intervalID ; value = pyValue 
        self.thisptr.set_value("power/nsize",value)
# ---------------- wrapper for libMesh::System ---------------- 
cdef class PylibMeshSystem:
    cdef System *thisptr  # hold C++ instance we're wrapping
    cdef public object variables # hold variable dictionary
    def __cinit__(self):
        self.thisptr = NULL 
        self.variables = {} 
    def AddFirstLagrangeVariable( self,PyStringVariableName):
        cdef char*  VariableName=  PyStringVariableName
        self.variables[PyStringVariableName] = AddFirstLagrangeVariable(self.thisptr,VariableName) 
    def AddConstantMonomialVariable( self,PyStringVariableName):
        cdef char*  VariableName=  PyStringVariableName
        self.variables[PyStringVariableName] = AddConstantMonomialVariable(self.thisptr,VariableName) 
    def AddStorageVectors( self,Ntime):
        CHKERR( AddStorageVectors(self.thisptr,"stored_local_solution",Ntime) )
    def StoreSystemTimeStep( self,TimeID):
        " store time step "
        CHKERR(StoreSystemTimeStep( self.thisptr ,TimeID))
    def SystemSolve(self):
        " solve the system"
        self.thisptr.solve();
    def CopySolutionVector(self, PylibMeshSystem newSystem not None):
        "copy a solution vector from newSystem to THIS system "
        CHKERR (CopySolutionVector(
               cython.operator.dereference(     self.thisptr) ,
               cython.operator.dereference(newSystem.thisptr) ) )
        return 
    def GetSolutionVector(self):
        "get a copy of the solution data "
        #print type(imagedata) , dir(imagedata)
        cdef Vec solutiondata = Vec()
        CHKERR (GetSolutionVector(
               cython.operator.dereference(self.thisptr) ,
                &solutiondata.vec) )
        #cdef PetscViewer newvwr = NULL
        #CHKERR (VecView(solutiondata.vec,newvwr) )
        return solutiondata
    def SetSolutionVector(self,Vec vecdata not None):
        "set copy python vector to c++ "
        CHKERR (SetSolutionVector(
               cython.operator.dereference(self.thisptr) ,
                 vecdata.vec) )
        #cdef PetscViewer newvwr = NULL
        #CHKERR (VecView(solutiondata.vec,newvwr) )
        return 
    #FIXME : prob separater this into a derived class
    #FIXME : note that only passing the System pointer 
    #FIXME : and dynamic casting in C++ routine
    def PennesSDASystemUpdateLaserPower(self,double newPower,int timeID):
        "assume system of type LITTSystem<PennesStandardDiffusionApproximation>"
        #print type(imagedata) , dir(imagedata)
        CHKERR(PennesSDASystemUpdateLaserPower( self.thisptr , newPower,timeID))
    #FIXME : prob separater this into a derived class
    #FIXME : note that only passing the System pointer 
    #FIXME : and dynamic casting in C++ routine
    def PennesSDASystemUpdateLaserPosition(self,position0,position1):
        "assume system of type LITTSystem<PennesStandardDiffusionApproximation>"
        #print type(imagedata) , dir(imagedata)
        CHKERR(PennesSDASystemUpdateLaserPosition( self.thisptr , position0[0], position0[1], position0[2], position1[0], position1[1], position1[2]))
    #FIXME : prob separater this into a derived class
    #FIXME : note that only passing the System pointer 
    #FIXME : and dynamic casting in C++ routine
    def PetscFEMSystemCreateNodeSetFromMask(self, double labelValue, int nodeSetID):
        "assume system of type PetscFEMSystem"
        #print type(imagedata) , dir(imagedata)
        return PetscFEMSystemCreateNodeSetFromMask( self.thisptr ,labelValue,nodeSetID  )
    #FIXME : prob separater this into a derived class
    #FIXME : note that only passing the System pointer 
    #FIXME : and dynamic casting in C++ routine
    def PetscFEMSystemSetupInitialConditions(self):
        CHKERR(PetscFEMSystemSetupInitialConditions(self.thisptr))
    #FIXME : prob separater this into a derived class
    #FIXME : note that only passing the System pointer 
    #FIXME : and dynamic casting in C++ routine
    def PetscFEMSystemUpdateTimeStep(self,TimeID):
        CHKERR(PetscFEMSystemUpdateTimeStep(self.thisptr,TimeID))
    #FIXME : prob separater this into a derived class
    #FIXME : note that only passing the System pointer 
    #FIXME : and dynamic casting in C++ routine
    def PetscFEMSystemGetSolnSubVector(self,VarID):
        "get a copy of the solution data for the input variable ID"
        #print type(imagedata) , dir(imagedata)
        cdef Vec solutiondata = Vec()
        CHKERR(PetscFEMSystemGetSolnSubVector(self.thisptr,
                                              VarID,&solutiondata.vec))
        #cdef PetscViewer newvwr = NULL
        #CHKERR (VecView(solutiondata.vec,newvwr) )
        return solutiondata
# ---------------- wrapper for libMesh::EquationSystems ---------------- 
cdef class PylibMeshEquationSystems:
    cdef EquationSystems *thisptr  # hold C++ instance we're wrapping
    def __cinit__(self, PylibMeshMesh meshdata not None, PylibMeshGetPot controldata not None ):
        self.thisptr = new EquationSystems(
            cython.operator.dereference(<MeshBase *>meshdata.thisptr) )
        CHKERR(SetParameterIni(&self.thisptr.parameters,controldata.thisptr))
    def init(self):
        "Initialize the system"
        self.thisptr.init()
    def PrintSelf(self):
        # base class pointer should call the inherited type
        self.thisptr.print_info()
        self.thisptr.parameters.print_info()
        CHKERR(FlushStdCoutCerr())
    def GetSystem(self,PyStringSystemName):
        cdef char* system_name= PyStringSystemName
        # get pointer to the system
        cdef PylibMeshSystem CurrSystem = PylibMeshSystem()
        CurrSystem.thisptr = GetSystem(system_name,self.thisptr)
        return CurrSystem
    def AddExplicitSystem(self,PyStringSystemName):
        cdef char* system_name= PyStringSystemName
        # get pointer to the system
        cdef PylibMeshSystem ExplicitSystem = PylibMeshSystem()
        ExplicitSystem.thisptr = AddExplicitSystem(system_name,self.thisptr)
        return ExplicitSystem
    def AddBackgroundSystem(self,PyStringSystemName):
        " Add BackgroundSystem "
        cdef char* system_name= PyStringSystemName
        cdef PylibMeshSystem BackgroundSystem = PylibMeshSystem()
        BackgroundSystem.thisptr = AddBackgroundSystem(system_name,self.thisptr)
        return BackgroundSystem
    def AddPennesSDASystem(self,PyStringSystemName,double DeltaT):
        " Add Pennes SDA with time step = DeltaT "
        cdef char* system_name= PyStringSystemName
        cdef PylibMeshSystem PennesSDASystem = PylibMeshSystem()
        PennesSDASystem.thisptr = AddPennesSDASystem(system_name,self.thisptr,DeltaT)
        return PennesSDASystem
    def AddPennesDeltaPSystem(self,PyStringSystemName,double DeltaT):
        " Add Pennes DeltaP with time step = DeltaT "
        cdef char* system_name= PyStringSystemName
        cdef PylibMeshSystem PennesDeltaPSystem = PylibMeshSystem()
        PennesDeltaPSystem.thisptr=AddPennesDeltaPSystem(system_name,self.thisptr,DeltaT)
        return PennesDeltaPSystem
    def AddPennesRFSystem(self,PyStringSystemName, double DeltaT):
        " Add Coupled Pennes Radiofrequency System w timestep = DeltaT"
        cdef char* system_name= PyStringSystemName
        cdef PylibMeshSystem PennesRFSystem = PylibMeshSystem()
        PennesRFSystem.thisptr = AddPennesRFSystem(system_name,self.thisptr,DeltaT)
        return PennesRFSystem
# ----------------- wrapper for libMesh::ExodusII_IO -------------------- 
cdef class PylibMeshExodusII_IO:
    cdef ExodusII_IO *thisptr  # hold a C++ instance which we're wrapping
    def __cinit__(self, PylibMeshMesh meshdata not None):
        self.thisptr = new ExodusII_IO( cython.operator.dereference(<MeshBase *>meshdata.thisptr) );
    def WriteTimeStep(self,PyStringFileName, 
                      PylibMeshEquationSystems EqnSystem not None, 
                           int timestep,double time):
        cdef char* file_name = PyStringFileName
        if(timestep < 1):
          raise RuntimeError("time step in exodus file must be >= 1" )
        self.thisptr.write_timestep( string(file_name), 
                        cython.operator.dereference(EqnSystem.thisptr),
                                     timestep , time)
    def WriteParameterSystems(self,
                      PylibMeshEquationSystems EqnSystem not None):
        self.thisptr.write_element_data( cython.operator.dereference(EqnSystem.thisptr) )
    def WriteEquationSystems(self,PyStringFileName, 
                      PylibMeshEquationSystems EqnSystem not None):
        cdef char* file_name = PyStringFileName
        self.thisptr.write_equation_systems ( string(file_name), 
                        cython.operator.dereference(EqnSystem.thisptr) )
# ------------wrapper for Imaging Interface -------------------- 
cdef class PytttkImaging:
    cdef Imaging *thisptr  # hold a C++ instance which we're wrapping
    def __cinit__(self, PylibMeshGetPot controldata not None,
                                     dimensions,origin,spacing):
        "setup imaging = dimensions, origin, spacing "
        self.thisptr = new Imaging( cython.operator.dereference(controldata.thisptr) );
        CHKERR(self.thisptr.SetHeaderData(
                             dimensions[0],dimensions[1],dimensions[2],
                                 origin[0],    origin[1],    origin[2] ,
                                spacing[0],   spacing[1],   spacing[2] ))
    def ProjectImagingToFEMMesh(self,PyStringSystemName,DefaultValue,
                                           Vec imagedata not None,
                      PylibMeshEquationSystems EqnSystem not None ): 
        "project imaging data onto libMesh data structures"
        #print type(imagedata) , dir(imagedata)
        cdef char* SystemName = PyStringSystemName
        # get pointer to the system
        cdef System *currentSystem = &(EqnSystem.thisptr.get_system(SystemName))
        CHKERR(self.thisptr.ProjectImagingToFEMMesh(currentSystem,imagedata.vec ,DefaultValue))
# ------------wrapper for Kalman Filter Interface -------------------- 
cdef class PytttkKalmanFilterMRTI:
    cdef KalmanFilter *thisptr  # hold a C++ base class pointer to instance which we're wrapping
    cdef public object systems    # hold libmesh Systems
    cdef public object LinAlgebra # store linear algebra type
    def AddSystems(self, PylibMeshEquationSystems EqnSystem not None, deltat):
        """setup equation systems needed for Kalman filter for MRTI 
        """
        print "Kalman adding StateSystem StateStdDev MRTIMean MRTIStdDev Systems"
        self.systems = {}
        # add systems for state / measurement mean and stddev
        self.systems['StateSystem'] = EqnSystem.AddPennesSDASystem("StateSystem",deltat) 
        # add state StdDev for plotting variance
        self.systems['StateStdDev'] = EqnSystem.AddExplicitSystem( "StateStdDev" ) 
        self.systems['StateStdDev'].AddFirstLagrangeVariable( "du0" )
        # hold imaging and uncertainty
        self.systems['MRTIMean']    = EqnSystem.AddExplicitSystem( "MRTIMean"    ) 
        self.systems['MRTIMean'].AddFirstLagrangeVariable( "u0*" )
        self.systems['MRTIStdDev']  = EqnSystem.AddExplicitSystem( "MRTIStdDev"  ) 
        self.systems['MRTIStdDev'].AddFirstLagrangeVariable( "du0*" )
    def CreateIdentityMeasurementMap(self, int nodeSetID):
        "setup identity map for state to measurements map "
        return self.thisptr.CreateIdentityMeasurementMap(nodeSetID)
    def CreateMeasurementMapFromImaging(self,PyStringSystemName, int nodeSetID):
        "setup state to measurements map from Imaging"
        cdef char*  SystemName = PyStringSystemName        
        return self.thisptr.CreateMeasurementMapFromImaging(SystemName,nodeSetID)
    def Setup(self, int numMeasurements):
        "setup petsc data structures"
        self.thisptr.Setup(numMeasurements)
        self.thisptr.printSelf()
    def CreateROIAverageMeasurementMap(self, image_roi, int NSubDZ,
			               origin, spacing):
        "setup petsc measurement matrix from  average over image thickness"
        self.thisptr.CreateROIAverageMeasurementMap(
                                     image_roi[0][0],image_roi[0][1],
                                     image_roi[1][0],image_roi[1][1],
                                     image_roi[2][0],image_roi[2][1], NSubDZ, 
                                           origin[0],  origin[1],  origin[2],
                                          spacing[0], spacing[1], spacing[2])
    def CreateROINodeSetMeasurementMap(self, int nodeSetID):
        "setup petsc measurement matrix from  node set"
        self.thisptr.CreateROINodeSetMeasurementMap(nodeSetID)
    def StatePredict(self,istep):
        "prediction of state"
        self.thisptr.StatePredict(istep)
    def CovariancePredict(self):
        "prediction of covariance"
        self.thisptr.CovariancePredict()
    def StateUpdate(self,PyStringKalmanGainMethod):
        "state update"
        cdef char* KalmanGainMethod= PyStringKalmanGainMethod
        self.thisptr.StateUpdate(KalmanGainMethod)
    def CovarianceUpdate(self):
        "covariance update"
        self.thisptr.CovarianceUpdate()
    def ProjectMeasurementMatrix(self,Vec SolutionData not None):
        "project a system soln to imaging space & return petsc vec"
        #print type(imagedata) , dir(imagedata)
        cdef Vec imagedata = Vec()
        CHKERR (
        self.thisptr.ProjectMeasurementMatrix(SolutionData.vec,&imagedata.vec)
               )
        return imagedata
    def MeasurementMatrixTranspose(self,PyStringSystemName,double ROIValue):
        "extract measurement data"
        cdef char*  SystemName = PyStringSystemName        
        self.thisptr.MeasurementMatrixTranspose(SystemName,ROIValue)
    def ExtractMeasurementData(self,Vec MRTIdata not None,
                              PyStringUncertaintySystemName):
        "extract measurement data"
        cdef char*  UncertaintySystemName= PyStringUncertaintySystemName
        self.thisptr.ExtractMeasurementData(MRTIdata.vec, UncertaintySystemName)
    def ExtractVarianceForPlotting(self,PyStringUncertaintySystemName):
        "extract diagonal of covariance"
        cdef char*  UncertaintySystemName= PyStringUncertaintySystemName
        self.thisptr.ExtractVarianceForPlotting(UncertaintySystemName)
    def ExtractCoVarianceForPlotting(self,PyStringUncertaintySystemName, int column):
        "extract diagonal of covariance"
        cdef char*  UncertaintySystemName= PyStringUncertaintySystemName
        self.thisptr.ExtractCoVarianceForPlotting(UncertaintySystemName,column)
cdef class PytttkSparseKalmanFilterMRTI(PytttkKalmanFilterMRTI):
    def __cinit__(self, PylibMeshEquationSystems EqnSystem not None, deltat):
        """setup Kalman filter for MRTI with sparse linear algebra
        """
        print "setting up sparse "
        self.LinAlgebra = "Sparse" # store linear algebra type
        self.thisptr = <KalmanFilter*>(new SparseKalmanFilter(EqnSystem.thisptr))
        PytttkKalmanFilterMRTI.AddSystems(self, EqnSystem , deltat)
cdef class PytttkDenseKalmanFilterMRTI(PytttkKalmanFilterMRTI):
    def __cinit__(self, PylibMeshEquationSystems EqnSystem not None, deltat):
        """setup Kalman filter for MRTI with dense linear algebra
        """
        print "setting up dense "
        self.LinAlgebra = "Dense" # store linear algebra type
        self.thisptr = <KalmanFilter*>(new DenseKalmanFilter(EqnSystem.thisptr))
        PytttkKalmanFilterMRTI.AddSystems(self, EqnSystem , deltat)
cdef class PytttkUncorrelatedKalmanFilterMRTI(PytttkKalmanFilterMRTI):
    def __cinit__(self, PylibMeshEquationSystems EqnSystem not None, deltat):
        """setup Kalman filter for MRTI with uncorrelated covariance
        """
        print "setting up uncorrelated "
        self.LinAlgebra = "UnCorr" # store linear algebra type
        self.thisptr = <KalmanFilter*>(new UncorrelatedKalmanFilter(EqnSystem.thisptr))
        PytttkKalmanFilterMRTI.AddSystems(self, EqnSystem , deltat)
###################### additional utility functions ################
# FIXME note case difference in calling and wrapping function 
# FIXME  to avoid conflict with previous definition
def GetFactorMat(Mat OrigMatrix not None, 
                 object PyStringSolvertype, object PyStringOrdertype ):
        cdef char*  Solvertype= PyStringSolvertype
        cdef char*  Ordertype = PyStringOrdertype 
        cdef Mat FactorMatrix = Mat()
        FactorMatrix.mat = GETFactorMat(OrigMatrix.mat,
                                        Solvertype,Ordertype)
        return FactorMatrix 
# FIXME note case difference in calling and wrapping function 
# FIXME  to avoid conflict with previous definition
def WeightedL2Norm(PylibMeshSystem StateSystem not None, # state system
                    object PyStringStateVariableName,# state system variable name
                    PylibMeshSystem IdealSystem not None, # measurement system
                    object PyStringIdealVariableName,# variable name
                    PylibMeshSystem NormalSystem not None,# normalize 
                    object PyStringNormalVariableName): # var name
        """compute the weighted ||.||_2 norm
               (Computed-Measurement)/Normalization, 
           inputs are the System names and corresponding Variable names
           for the calculation
        """
        #variable names
        cdef char*  StateVariableName=  PyStringStateVariableName        
        cdef char*  IdealVariableName=  PyStringIdealVariableName        
        cdef char* NormalVariableName= PyStringNormalVariableName
        # compute the weighted l2 norm
        return WEIGHTEDL2Norm(
                      StateSystem.thisptr , StateVariableName,
                      IdealSystem.thisptr , IdealVariableName,
                     NormalSystem.thisptr ,NormalVariableName 
                             )
