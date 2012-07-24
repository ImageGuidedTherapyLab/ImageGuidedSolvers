//prevent multiple inclusions of header file
#ifndef SOLVEHP3D_H
#define SOLVEHP3D_H


// forward declarations
class AppSolve;
class ExodusII_IO; 
namespace dddas
{
class Image; 
}
//----------------------------------------------------------------------
//   class UniversalCntrl  
//   latest revision  - Aug 07
//
//   purpose          - this base class data structure contains 
//                      general information
//----------------------------------------------------------------------
class UniversalCntrl {
public:
  //constructor (minimal variables setup before debugger setup)
  UniversalCntrl(GetPot &); 

  //setup the majority of the variables after the debugger setup
  void GenSetup(GetPot &); 
 
  //echo data
  void printSelf();

  std::string MRTIDATA_TYPE,FieldFile;

  PetscInt NRELEB,  /* NELEMS # of elements locally on this processor */
           NDOF_FORWARD[3],   
           NDOF_ADJOINT[3];  
           /* NDOF_???[0]: local # of dof OWNED by this processor in a Petsc Vec
              NDOF_???[1]: total # of dof shared amongst processors 
              NDOF_???[2]: total # of dof OWNED AND SHARED by this processor 
              the numbers of dofs are stored for the forward and adjoint problem
              since they may be different for goal oriented error estimation*/
  std::vector<PetscInt> locmap_for,locmap_adj; // data structures to store local to
                                          // global mappings
  PetscScalar scalefact;  /* scale factor for penalty term */
};
//----------------------------------------------------------------------
//   class DroneControl  
//   latest revision    - Aug 07
//
//   purpose            - this data structure contains the time step
//                        control information for each computational group
//----------------------------------------------------------------------
class DroneControl {
public:
  DroneControl(GetPot &,UniversalCntrl &, PetscInt ); //constructor
  // data structures to gather data from computations groups
  std::vector<PetscInt> recvcnts_elm, recvcnts_for, recvcnts_adj,
                   displs_elm,   displs_for,   displs_adj;
  std::string      pde,  compobj;
};

// signatures
PetscErrorCode InitControlGlobals(DroneControl &);
#endif
