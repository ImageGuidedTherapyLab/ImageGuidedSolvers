#cdef extern from "signalmodel.h" namespace "IterativeReconstruction": 
from libcpp.vector  cimport vector
cdef extern from "signalmodel.h": 
    cdef cppclass SignalModelInterface:
        SignalModelInterface()
        double m_echoTime
        void printSelf ()
        int SetuplibMesh (MPI_Comm)
        int SetupImaging(int,int,int,
                      double, double, 
                      double, double, 
                      double, double) 
        int SetupStructuredGrid(int,int,int,
                      double, double, double, 
                      double, double, double) 
        int InitializeEquationSystems()
        int ProjectImagingToFEMMesh (char *, PetscVec)
        int WriteEquationSystems(char *)
        int AddExplicitSystem (  char *,int)
        void AssembleSignal (vector[double]*,
                             vector[double]*,
                             vector[double]*)
        int SetupLinearGradientModel(
                      double, double, double, double) 
        int SetupSpiralGradientModel(
                      double, double, double, double, double) 
        int SetupEPIGradientModel(
                      double, double, double, double, int) 
