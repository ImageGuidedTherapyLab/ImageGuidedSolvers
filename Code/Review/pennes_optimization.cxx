// C++ include files 
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <vector>

// Petsc/Tao include files
#include "tao.h" 
#include "src/tao_impl.h"    // tao solvers
#include "private/kspimpl.h" // access to the KSP data structures

// libMesh include files
#include "libmesh.h"
#include "mesh.h"
#include "mesh_function.h"
#include "equation_systems.h"
#include "numeric_vector.h"
#include "exact_solution.h"
#include "fe.h"
#include "quadrature_gauss.h"
#include "dof_map.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "elem.h"
#include "boundary_info.h"
#include "nonlinear_solver.h"
#include "nonlinear_implicit_system.h"
#include "linear_implicit_system.h"
#include "parallel.h" // file parser
#include "getpot.h" // file parser

// dddas include files
using namespace std;
#include "baseInfo.h" 
#include "applicationContext.h" 
#include "transient_inverse_system.h"
#include "fortrandftranslate.h" // header to call fortran routines
#include "parser_defaults.h" // ensure same default values between C and Fortran
#include "mri_params.h" // mri params
#include "solvehp3d.h" // main header
#include "variable_map.h" // ensure same default values between C and Fortran
#include "dddas.h" // for profiling
#include "pennesVerification.h" 

/* -------------- User-defined routines ---------- */
#include "pennes_optimization.h"  // interface to pennes model optimization routines
#include "pennes_solver.h"  // interface to pennes model optimization routines
extern FORTRAN_FUNCTION
{
#include "pennes_model.h"  // interface to pennes model
//  void FORTRAN_NAME(get__field)(PetscInt*,PetscScalar*,PetscInt*,PetscInt*,PetscInt*, PetscScalar(*)( PetscInt*,PetscInt*), PetscScalar(*)( PetscInt*,PetscInt*,PetscInt*));
//  PetscScalar FORTRAN_NAME(get_vert_temp)( PetscInt*,PetscInt*);
//  PetscScalar FORTRAN_NAME(get_node_temp)( PetscInt*,PetscInt*,PetscInt*);
//  PetscScalar FORTRAN_NAME(qoieval_spacetime)(PetscInt*, PetscInt*, PetscInt*);
//  PetscScalar FORTRAN_NAME(qoieval_lesbegue2)(PetscInt*, PetscInt*, PetscInt*);
//  PetscScalar FORTRAN_NAME(getdeltat)();
//  void FORTRAN_NAME(storefield)(PetscInt*,PetscScalar*,PetscInt*,PetscInt*, const PetscInt*, void(*)( PetscInt*,PetscInt*,PetscScalar*), void(*)( PetscInt*,PetscInt*,PetscInt*,PetscScalar*));
//  void        FORTRAN_NAME(put_vert_dual)( PetscInt*,PetscInt*,PetscScalar*);
//  void        FORTRAN_NAME(put_node_dual)( PetscInt*,PetscInt*,PetscInt*,PetscScalar*);
//  void FORTRAN_NAME(nelconb)(PetscInt *,PetscInt *);
//  void FORTRAN_NAME(getglobalmap)(PetscInt*,PetscInt*,PetscInt(*)(PetscInt *),PetscInt*,PetscInt*,PetscInt*);
//  void FORTRAN_NAME(initadjointfield)(PetscScalar *,PetscInt *,PetscInt *);
//  void FORTRAN_NAME(elemdualload)(PetscInt *,PetscInt *,PetscInt *,PetscScalar *,PetscInt *, PetscInt *,void (*)(PetscScalar*,PetscInt*,PetscInt*,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*));
//  void FORTRAN_NAME(elemjac)(PetscInt *,PetscInt *,PetscInt *,PetscScalar *,PetscInt *,PetscInt(*)(PetscInt *),PetscInt*, void(*)(PetscScalar*, PetscScalar*, PetscScalar*, PetscScalar* , PetscScalar*, PetscScalar*), void(*)(PetscInt* , PetscScalar* , PetscScalar* , PetscScalar*) );
//  PetscInt FORTRAN_NAME(get_order_adjnt)(PetscInt*) ;
//  void FORTRAN_NAME(adjoint_grad)(PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscScalar*, void(*)(const PetscScalar* , PetscScalar* , PetscScalar* , PetscScalar*) );
//  void FORTRAN_NAME(zero_gradbeta)(PetscInt*,PetscInt*);
//  PetscScalar FORTRAN_NAME(accumulate_error)(PetscInt*,PetscInt*);
//  PetscScalar FORTRAN_NAME(getdeltat)();
}

// Global Identifiers
// control variable id
static const PetscInt idstate=ISTATE , 
                      idadjnt=IADJNT , 
                      nzero  =0 ; 


/* -------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "Lesbegue2_QOI"
petscScalar Lesbegue2_QOI( AppSolve *user)
{
  DroneControl   &QOIInfo = *user->QOIInfo; // time step info
  PetscErrorCode info;

  /* compute the data structures for the QOI
     ComputeQoiData is a function pointer to  
                         verifcomputeidealfield: for verification problems
                         readMRTI_data    : for calibration problems
                         LoadIdealArr     : arrhenius base optimize
                         LoadIdeal_TS     : two-state base optimize
                         LoadIdealHSP     : HSP base optimize
                         LoadIdealTmp     : for temp base optimize */
  if(QOIInfo.ComputeQoiData){info=QOIInfo.ComputeQoiData(user);CHKERRQ(info);}

  PetscScalar qoi_loc = 0.0;
  //  use the last computed damage/hsp field to evaluate the qoi
//  for(mdleID=MdleData->begin() ; mdleID!=MdleData->end() ; mdleID++){
//     qoi_loc = qoi_loc + FORTRAN_NAME(qoieval_lesbegue2)(&mdleID->mdle,
//                                               &mdleID->idelem,&user->Nstephi);
//  }
  return qoi_loc;
}

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  read the ideal field into the zdofs_ideal data structure for the initial
  time instance and store the computation of the integral in the remaining
  time instances
*/
#undef __FUNCT__
#define __FUNCT__ "LoadIdealHSP"
PetscErrorCode LoadIdealHSP( AppSolve *user ){
  PetscInt idstep;
  PetscErrorCode info;

  // used for error handling
  PetscFunctionBegin; 

  // calculate the local the damage field for each time step from the beginning
  //  to account for the entire history
  for(idstep = user->Init_Nsteplo + 1 ; idstep<=user->Nstephi ; idstep++){
//     for(mdleID=MdleData->begin() ; mdleID!=MdleData->end() ; mdleID++){
//         FORTRAN_NAME(field2datastruct)(&mdleID->mdle,&mdleID->idelem, &idstep,
//                                                 &FORTRAN_NAME(getlochsp), 
//                                                 &FORTRAN_NAME(elemideal)   );
//     }
  }
  // before accumulation must zero the initial ideal field 
//  for(mdleID=MdleData->begin() ; mdleID!=MdleData->end() ; mdleID++){
//      FORTRAN_NAME(zeroideal)(&mdleID->mdle,&user->Init_Nsteplo);
//  } 
  // accumulate the time history of the damage field 
  for(idstep = user->Init_Nsteplo + 1 ; idstep<=user->Nstephi ; idstep++){
//      FORTRAN_NAME(unvisit)(); 
//      for(mdleID=MdleData->begin() ; mdleID!=MdleData->end() ; mdleID++){
//          FORTRAN_NAME(accumulateideal)(&mdleID->mdle,&idstep);
//      } 
  }
  // put ideal field into FEM ideal data structures at time instance 0
  idstep=0;
//  for(mdleID=MdleData->begin() ; mdleID!=MdleData->end() ; mdleID++){
//      FORTRAN_NAME(field2datastruct)(&mdleID->mdle,&mdleID->idelem, &idstep,
//                                               &FORTRAN_NAME(getidealfield), 
//                                               &FORTRAN_NAME(elemideal)   );
//  }

  PetscFunctionReturn(0);
}
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  read the ideal field into the zdofs_ideal data structure for the initial
  time instance and store the computation of the integral in the remaining
  time instances
*/
#undef __FUNCT__
#define __FUNCT__ "LoadIdeal_TS"
petscErrorCode LoadIdeal_TS( AppSolve *user){
  PetscInt idstep;
  PetscErrorCode info;

  // used for error handling
  PetscFunctionBegin; 

  // calculate the local the damage field for each time step from the beginning
  //  to account for the entire history
  for(idstep = user->Init_Nsteplo + 1 ; idstep<=user->Nstephi ; idstep++){
//     for(mdleID=MdleData->begin() ; mdleID!=MdleData->end() ; mdleID++){
//         FORTRAN_NAME(field2datastruct)(&mdleID->mdle,&mdleID->idelem, &idstep,
//                                                 &FORTRAN_NAME(getloctwodam), 
//                                                 &FORTRAN_NAME(elemideal)   );
//     }
  }
  // before accumulation must zero the initial ideal field 
//  for(mdleID=MdleData->begin() ; mdleID!=MdleData->end() ; mdleID++){
//      FORTRAN_NAME(zeroideal)(&mdleID->mdle,&user->Init_Nsteplo);
//  } 
  // accumulate the time history of the damage field 
  for(idstep = user->Init_Nsteplo + 1 ; idstep<=user->Nstephi ; idstep++){
//      FORTRAN_NAME(unvisit)(); 
//      for(mdleID=MdleData->begin() ; mdleID!=MdleData->end() ; mdleID++){
//          FORTRAN_NAME(accumulateideal)(&mdleID->mdle,&idstep);
//      } 
  }
  // put ideal field into FEM ideal data structures at time instance 0
  idstep=0;
//  for(mdleID=MdleData->begin() ; mdleID!=MdleData->end() ; mdleID++){
//      FORTRAN_NAME(field2datastruct)(&mdleID->mdle,&mdleID->idelem, &idstep,
//                                               &FORTRAN_NAME(getidealfield), 
//                                               &FORTRAN_NAME(elemideal)   );
//  }

  PetscFunctionReturn(0);
}

/* Function to compute stored solutions for plotting and QOI evaluation purposes.  */
number pennes_damage (const Point& p,
                      const Parameters& parameters,
                      const std::string& sys_name,
                      const std::string& unknown_name)
{
 MeshFunction *stored_values = parameters.get<MeshFunction*> ("MeshFunction");
 PetscScalar usln            = (*stored_values)(p);
 PetscScalar FEM_DT          = parameters.get<Real> ("FEM_DT");
 return FORTRAN_NAME(getlocarrdam)(&usln,&FEM_DT); 
}
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  read the ideal field into the zdofs_ideal data structure for the initial
  time instance and store the computation of the integral in the remaining
  time instances
*/
#undef __FUNCT__
#define __FUNCT__ "LoadIdealArr"
petscErrorCode LoadIdealArr( AppSolve *user ){
  PetscErrorCode info;

  // used for error handling
  PetscFunctionBegin; 

  // deltat
  PetscScalar &FEM_DT= baseInfo::FEM_DT;

  // get the forward system
  TransientInverseNonlinearImplicitSystem & forward_system =
        user->_main_es->get_system<TransientInverseNonlinearImplicitSystem>("StateSystem");

  // get the ideal system
  TransientInverseExplicitSystem & ideal_system =
        user->_main_es->get_system<TransientInverseExplicitSystem>("IdealSystem");

  // store deltat for the projection
  user->_main_es->parameters.set<Real> ("FEM_DT") = FEM_DT;

  // calculate the local the damage field for each time step from the beginning
  //  to account for the entire history
  for(PetscInt idstep = user->Init_Nsteplo + 1 ; idstep<=user->Nstephi ; idstep++)
   {
     // Prepare MeshFunction to project the stored solution for plotting
     MeshFunction state_values(*user->_main_es,
                        *forward_system.vector_solution.at(idstep),
       		       forward_system.get_dof_map(),
                       forward_system.variable_number("u1"));
     state_values.init();
     user->_main_es->parameters.set<MeshFunction*>("MeshFunction") = 
                                                           &state_values;
     //project stored solution for plotting
     ideal_system.project_solution(pennes_damage,NULL,
                                     user->_main_es->parameters);

     *ideal_system.vector_solution.at(idstep)=
                                        *ideal_system.current_local_solution;
   }

  // before accumulation must zero the initial ideal field 
  ideal_system.vector_solution.at(user->Init_Nsteplo)->zero();

  // accumulate the time history of the damage field 
  for(PetscInt idstep = user->Init_Nsteplo + 1 ; idstep<=user->Nstephi ; idstep++)
   {
     ideal_system.vector_solution.at(idstep)->add(
                   *ideal_system.vector_solution.at(idstep-1));
   }

  // put full damage into current local solution
  *ideal_system.current_local_solution = 
        *ideal_system.vector_solution.at(user->Nstephi);

  PetscFunctionReturn(0);
}
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
#undef __FUNCT__
#define __FUNCT__ "CurrentIdealArr"
petscErrorCode CurrentIdealArr( AppSolve *user ){
  PetscErrorCode info;

  // used for error handling
  PetscFunctionBegin; 

  // deltat
  PetscScalar &FEM_DT= baseInfo::FEM_DT;

  // get the forward system
  TransientInverseNonlinearImplicitSystem & forward_system =
        user->_main_es->get_system<TransientInverseNonlinearImplicitSystem>("StateSystem");

  // get the ideal system
  TransientInverseExplicitSystem & ideal_system =
        user->_main_es->get_system<TransientInverseExplicitSystem>("IdealSystem");

  // store deltat for the projection
  user->_main_es->parameters.set<Real> ("FEM_DT") = FEM_DT;

  // store the current damage field
  *ideal_system.old_vector_solution.at(0)=
                                     *ideal_system.current_local_solution;

  // calculate the local the damage field due to the current time step
  // Prepare MeshFunction to project the stored solution for plotting
  MeshFunction state_values(*user->_main_es,
                     *forward_system.current_local_solution,
    		       forward_system.get_dof_map(),
                    forward_system.variable_number("u1"));
  state_values.init();
  user->_main_es->parameters.set<MeshFunction*>("MeshFunction") = 
                                                        &state_values;
  //project stored solution for plotting
  ideal_system.project_solution(pennes_damage,NULL,
                                  user->_main_es->parameters);

  //Accumulate the Damage
  ideal_system.current_local_solution->add(
                   *ideal_system.old_vector_solution.at(0));


  PetscFunctionReturn(0);
}
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  read the ideal field into the zdofs_ideal data structure for all
  time instanced 
*/
#undef __FUNCT__
#define __FUNCT__ "LoadIdealTmp"
petscErrorCode LoadIdealTmp( AppSolve *user ){
  PetscInt idstep;
  PetscErrorCode info;

  // used for error handling
  PetscFunctionBegin; 

  // get the ideal system
  TransientInverseExplicitSystem & ideal_system =
        user->_main_es->get_system<TransientInverseExplicitSystem>("IdealSystem");

  for(idstep=user->Nsteplo ; idstep<=user->Nstephi ; idstep++){
     // put ideal field into FEM ideal data structures
     ideal_system.project_solution(pennes_ideal_field,NULL,
                                     user->_main_es->parameters);

     *ideal_system.vector_solution.at(idstep)=
                                        *ideal_system.current_local_solution;
  }

  PetscFunctionReturn(0);
}

//static CntrlVars *CurrentControlVars;
///* -------------------------------------------------------------------
//    set pointer for kernel compuation of adjoint used by Fortran element routine
//   ------------------------------------------------------------------- */
//#undef __FUNCT__
//#define __FUNCT__ "SetControlVariablePointer"
//PetscErrorCode SetControlVariablePointer(CntrlVars *OptParams){
//  // used for error handling
//  PetscFunctionBegin;
//  CurrentControlVars =  OptParams;
//  PetscFunctionReturn(0);
//}
///* -------------------------------------------------------------------
//    map local element contributions to the local gradient vector
//   ------------------------------------------------------------------- */
//void GradKernel(PetscInt *IdParam,
//          PetscScalar *Xpoint,PetscInt *Idstep     ,PetscScalar *Deltat,
//          PetscScalar *Usln  ,PetscScalar *Usln_prev, 
//          PetscScalar *Adiff ,PetscScalar *Creac    ,PetscScalar *Source){
//
//  CurrentControlVars->Parameters.at(*IdParam).evalgrad(Xpoint,Idstep,Deltat,
//                                                      Usln  , Usln_prev, 
//                                                      Adiff , Creac, Source);
//
//}
//
//#undef __FUNCT__
//#define __FUNCT__ "gradkernel_"
//extern FORTRAN_FUNCTION void FORTRAN_NAME(gradkernel)(PetscInt *IdParam,
//          PetscScalar *Xpoint,PetscInt *Idstep     ,PetscScalar *Deltat,
//          PetscScalar *Usln  ,PetscScalar *Usln_prev, 
//          PetscScalar *Adiff ,PetscScalar *Creac    ,PetscScalar *Source){
//
//  GradKernel(IdParam, Xpoint,Idstep,Deltat,Usln,Usln_prev,Adiff,Creac,Source);
//
//}

