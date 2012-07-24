from petsc4py.PETSc cimport MPI_Comm,GetComm
from petsc4py.PETSc cimport PetscVec,Vec

cdef extern from "petsc.h":
    MPI_Comm PETSC_COMM_SELF

include "wrapsignalmodel.pxi"
# --------------------------------------------------------------------

## Vile hack for raising a exception 
#
#cdef extern from *:
#    enum: PETSC_ERR_PYTHON "(-1)"
#
#cdef extern from *:
#    void pyx_raise"__Pyx_Raise"(object, object, void*)
#
#cdef extern from *:
#    void *PyExc_RuntimeError
#
#cdef object PetscError = <object>PyExc_RuntimeError
#
#cdef inline int SETERR(int ierr):
#    if (<void*>PetscError):
#        pyx_raise(PetscError, ierr, NULL)
#    else:
#        pyx_raise(<object>PyExc_RuntimeError, ierr, NULL)
#    return ierr
#
#cdef inline int CHKERR(int ierr) except -1:
#    if ierr == 0: return 0                 # no error
#    if ierr == PETSC_ERR_PYTHON: return -1 # error in Python call
#    if (<void*>PetscError):
#        pyx_raise(PetscError, ierr, NULL)
#    else:
#        pyx_raise(<object>PyExc_RuntimeError, ierr, NULL)
#    return -1
#
#if 0: raise RuntimeError # Do not remove this line !!!
cdef inline int CHKERR(int ierr) except -1:
    if ierr == 0: return 0                 # no error
    else:
       raise RuntimeError
    return -1

# --------------------------------------------------------------------
cdef class PySignalModel:
    cdef SignalModelInterface *thisptr  # hold a C++ instance which we're wrapping
    def __cinit__(self, double TE):
        self.thisptr = new SignalModelInterface()
        self.thisptr.m_echoTime = TE
    def __dealloc__(self):
        del self.thisptr
    def SetuplibMesh(self,comm):
        cdef MPI_Comm ccomm = GetComm(comm, PETSC_COMM_SELF)
        self.thisptr.SetuplibMesh(ccomm)
    def SetupStructuredGrid(self,dimensions,xbounds,ybounds,zbounds):
        "setup SetupStructuredGrid input = dimensions,xmin,xmax,ymin,ymax,zmin,zmax"
        self.thisptr.SetupStructuredGrid(dimensions[0],dimensions[1],dimensions[2],
                                         xbounds[0] , xbounds[1] ,
                                         ybounds[0] , ybounds[1] ,
                                         zbounds[0] , zbounds[1]  )
    def SetupImaging(self,dimensions,origin,spacing):
        "setup imaging = dimensions, origin, spacing "
        self.thisptr.SetupImaging(dimensions[0],dimensions[1],dimensions[2],
                                         origin[0] , origin[1] , origin[2] ,
                                         spacing[0], spacing[1], spacing[2] )
    def InitializeEquationSystems(self):
        CHKERR(self.thisptr.InitializeEquationSystems())
    def printSelf(self):
        "print class info"
        self.thisptr.printSelf()
    def WriteEquationSystems(self,PyStringFileName):
        cdef char* file_name = PyStringFileName
        return self.thisptr.WriteEquationSystems(file_name);
    def AddExplicitSystem(self,PyStringSystemName,VarOrder):
        cdef char* system_name= PyStringSystemName
        CHKERR(self.thisptr.AddExplicitSystem(system_name,VarOrder))
    def ProjectImagingToFEMMesh(self,PyStringSystemName, Vec t2starmap not None):
        "setup libMesh data structures"
        #print type(t2starmap) , dir(t2starmap)
        cdef char* system_name= PyStringSystemName
        return self.thisptr.ProjectImagingToFEMMesh(system_name,t2starmap.vec)
    def AssembleSignal(self,TimeValues):
        """
        perform matvec product and compute the real and imaginary 
        part of the signal at a given time instance. assume
        all data is preload and push the work onto Python
        TODO for higher compuational speed an ABC for the 
        gradients can be provided in C++. 
        This can be coded post development for efficiency.
        THe current design allows the development to be performed
        in python rather than C++ :) 
        """
        cdef vector[double] time_values
        cdef vector[double] RealSignal
        cdef vector[double] ImagSignal
        cdef double time
        for time in TimeValues:
           time_values.push_back(time)
        self.thisptr.AssembleSignal(&time_values,&RealSignal,&ImagSignal)
        return [(RealSignal[i],ImagSignal[i]) for i in range(len(TimeValues))] 
    def SetupLinearGradientModel(self,Gradient,Larmor):
        return self.thisptr.SetupLinearGradientModel(Gradient[0],
                                                     Gradient[1],
                                                     Gradient[2],
                                                     Larmor)
    def SetupSpiralGradientModel(self,Gradient,Larmor,AngularVel):
        return self.thisptr.SetupSpiralGradientModel(Gradient[0],
                                                     Gradient[1],
                                                     Gradient[2],
                                                     Larmor,
                                                     AngularVel)
    def SetupEPIGradientModel(self,Gradient,Larmor,NStepKspace):
        return self.thisptr.SetupEPIGradientModel(Gradient[0],
                                                     Gradient[1],
                                                     Gradient[2],
                                                     Larmor,
                                                     NStepKspace)
