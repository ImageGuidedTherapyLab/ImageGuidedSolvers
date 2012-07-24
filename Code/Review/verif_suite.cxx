// C++ include files that we need
#include <iostream>
#include <vector>

// libmesh
#include "libmesh.h"
#include "equation_systems.h"
#include "fe.h"
#include "fe_interface.h"
#include "quadrature_gauss.h"
#include "dof_map.h"
#include "sparse_matrix.h"
#include "elem.h"
#include "string_to_enum.h"
#include "getpot.h"
#include "boundary_info.h"

// The nonlinear solver and system we will be using
#include "nonlinear_implicit_system.h"
#include "linear_implicit_system.h"

// pdeopt
#include "applicationContext.h"
#include "src/tao_impl.h" // tao and petsc solvers
#include "baseInfo.h" // basic information
#include "transient_inverse_system.h"
#include "pennesSolver.h"
#include "pennesVerification.h" 
#include "quantityOfInterest.h"
#include "dddas.h"
#include "variable_map.h"

#include "mesh.h"

#include "fortrandftranslate.h" // header to call fortran routines
extern FORTRAN_FUNCTION
{ 
#include "pennes_model.h" // pennes model interface
// void FORTRAN_NAME(storefield)(PetscInt*,PetscScalar*,PetscInt*,PetscInt*,PetscInt*, void(*)( PetscInt*,PetscInt*,PetscScalar*), void(*)( PetscInt*,PetscInt*,PetscInt*,PetscScalar*));
// void FORTRAN_NAME(get__field)(PetscInt*,PetscScalar*,PetscInt*,PetscInt*,PetscInt*, PetscScalar(*)( PetscInt*,PetscInt*), PetscScalar(*)( PetscInt*,PetscInt*,PetscInt*));
// PetscScalar FORTRAN_NAME(get_vert_temp)( PetscInt*,PetscInt*);
// PetscScalar FORTRAN_NAME(get_node_temp)( PetscInt*,PetscInt*,PetscInt*);
// void        FORTRAN_NAME(put_vert_temp)( PetscInt*,PetscInt*,PetscScalar*);
// void        FORTRAN_NAME(put_node_temp)( PetscInt*,PetscInt*,PetscInt*,PetscScalar*);
// void FORTRAN_NAME(update_iterno)(PetscInt*,PetscInt*,PetscInt*,PetscInt*);
// void FORTRAN_NAME(soln2ideal)();
// void FORTRAN_NAME(getndiv)(PetscInt*,PetscInt*);
// void FORTRAN_NAME(print_error)(PetscInt*,PetscInt*);
}
/* ------------------------------------------------------------------- 
   computeideal_mc - compute a fake ideal temperature field
                            used for code verification

   Input Parameters:
   ctx  - user-defined context 
   
   Required Routines: inittempfield must be called before this routine
*/
//#undef __FUNCT__
//#define __FUNCT__ "computeideal_mc"
//PetscErrorCode computeideal_mc( AppSolve         *user  ){
//  UniversalCntrl *GenInfo = user->GenInfo; // global data
//  DroneControl *qoiOptimizer = user->qoiOptimizer; // global data
//  Vec  *Z;
//  PetscInt nsteps , iii, idstep;
//  PetscScalar *solnloc;
//  FormFcnGaussType savefncter ;
//  PetscErrorCode info;
//
//  // used for error handling
//  PetscFunctionBegin; 
//
//  nsteps = user->Nstephi - user->Nsteplo + 1 ;
//
//  // first store the data from hp3d data structures
//  VecDuplicateVecs(user->xloc,nsteps,&Z);
//  for(iii = 0;iii < nsteps ; iii++){
//    idstep = user->Nsteplo + iii;
//    info = VecGetArray(Z[iii],&solnloc); CHKERRQ(info);
////    for(mdleID=MdleData->begin() ; mdleID!=MdleData->end() ; mdleID++){
////       FORTRAN_NAME(get__field)(&mdleID->mdle,solnloc,
////                            &GenInfo->NDOF_FORWARD[2],&idstep,&idstate,
////                                               FORTRAN_NAME(get_vert_temp) , 
////                                               FORTRAN_NAME(get_node_temp) );
////    }
//    info = VecRestoreArray(Z[iii],&solnloc); CHKERRQ(info);
//  }
//
//  savefncter = qoiOptimizer->formfcngauss[0][0]; // save function pointer for later
//
//  // compute the ideal field from the ideal Monte Carlo Data
//  qoiOptimizer->formfcngauss[0][0] = FORTRAN_NAME(pennesidealmonte);
//  ForwardSolve(user);
//
//  qoiOptimizer->formfcngauss[0][0] = savefncter; // return fnc pointer to original state
//
//  //for verifying qoi calculation
//  //transfer data to qoi data structures
////  FORTRAN_NAME(soln2ideal)();
//  // return the original data to hp3d data structures
//  for(iii = 0;iii < nsteps ; iii++){
//     idstep = user->Nsteplo + iii;
//     info = VecGetArray(Z[iii],&solnloc); CHKERRQ(info);
////     FORTRAN_NAME(storefield)(&mdleID->mdle,solnloc,
////                          &GenInfo->NDOF_FORWARD[2],&idstep,&idstate,
////                                             FORTRAN_NAME(put_vert_temp) , 
////                                             FORTRAN_NAME(put_node_temp) );
//     info = VecRestoreArray(Z[iii],&solnloc); CHKERRQ(info);
//  }
//  VecDestroyVecs(Z,nsteps); // Free memory
//
//  PetscFunctionReturn(0);
//}
