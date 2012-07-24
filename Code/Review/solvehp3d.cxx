// C++ include files 
#include <iostream>
#include <fstream>
#include <vector>

// Petsc/Tao include files
#include "tao.h" // petsc solvers w/ Fortran Interface

// libMesh include files
#include "libmesh.h"
#include "mesh.h"
#include "equation_systems.h"
#include "exact_solution.h"
#include "getpot.h"

using namespace std;

// dddas include files
#include "baseInfo.h"
#include "fortrandftranslate.h" // header to call fortran routines
#include "parser_defaults.h" // ensure same default values between C and Fortran
#include "mri_params.h" // mri params
#include "solvehp3d.h" // main header
#include "applicationContext.h"
#include "read_mridata.h"
#include "dddas.h"
#include "pennes_optimization.h"
#include "pennesVerification.h"

extern FORTRAN_FUNCTION 
{
#include "pennes_model.h" // interface to pennes model routines 
}

FcnterHp3d W_0_fcns , W_N_fcns , W_I_fcns , W_D_fcns ,
           W_2_fcns , W_NI_fcns, W_ID_fcns,
           K_0_fcns , K_1_fcns , K_2_fcns , K_3_fcns ,
           X_0_fcns , Y_0_fcns , Z_0_fcns , 
           MU_A_fcns, MU_S_fcns, POW_fcns , ERR_fcns ; 

      // must map to the correct gradient position term for
      // power as a function of time optimization
PetscInt time_offset_pow(PetscInt Istep, PetscInt Isteps_per_ideal, 
                                                  PetscInt Ideal_nzero){
   // this value is returned to a C - data structure where the indexing 
   // begins w/ ZERO!!!! this is not used for the actual POW data structure!!
   // only used w/ local gradient in which indexing begins w/ zero
   if( Istep % Isteps_per_ideal )
      return Istep / Isteps_per_ideal - Ideal_nzero      ; 
   else
      return Istep / Isteps_per_ideal - Ideal_nzero - 1  ;
}

PetscInt Given3IntsReturnPetscIntZero( PetscInt Istep, 
                    PetscInt Isteps_per_ideal, PetscInt Ideal_nzero){return 0;}
PetscInt    GivenIntReturnPetscIntZero(PetscInt *Mdle){return 0;}
PetscInt    ReturnPetscIntOne(){return 1;}
PetscTruth  ReturnPETSC_TRUE(){          return PETSC_TRUE;}
PetscTruth  ReturnPETSC_FALSE(void *dum){return PETSC_FALSE;}
PetscScalar ReturnPetscScalarZero(PetscInt *Id,PetscScalar *Scalefact)
{
  return 0.0;
}
PetscScalar ReturnPetscScalarZero(OPTGAUSSARG)
{
  return 0.0; 
}
void NoParams(const PetscScalar *Params,PetscInt *SizeParam){ return; }
void NoParams(      PetscScalar *Params,PetscInt *SizeParam){ return; }
void getTAO_NINFINITY(PetscScalar *Params,PetscInt *SizeParam){ 
     *Params = TAO_NINFINITY; return;
}
void getTAO_INFINITY( PetscScalar *Params,PetscInt *SizeParam){ 
     *Params = TAO_INFINITY ; return;
}

#undef __FUNCT__
#define __FUNCT__ "InitControlGlobals"
PetscErrorCode InitControlGlobals(DroneControl  &QOIInfo)
{

  PetscFunctionBegin; 

     ERR_fcns.putparam =  NoParams                    ;
     W_0_fcns.putparam =  FORTRAN_NAME(putparam_w_0)  ;
     W_N_fcns.putparam =  FORTRAN_NAME(putparam_w_n)  ;
     W_I_fcns.putparam =  FORTRAN_NAME(putparam_w_i)  ;
     W_D_fcns.putparam =  FORTRAN_NAME(putparam_w_d)  ;
     W_2_fcns.putparam =  FORTRAN_NAME(putparam_w_2)  ;
     W_NI_fcns.putparam=  FORTRAN_NAME(putparam_w_ni)  ;
     W_ID_fcns.putparam=  FORTRAN_NAME(putparam_w_id)  ;
     K_0_fcns.putparam =  FORTRAN_NAME(putparam_k_0)  ;
     K_1_fcns.putparam =  FORTRAN_NAME(putparam_k_1)  ;
     K_2_fcns.putparam =  FORTRAN_NAME(putparam_k_2)  ;
     K_3_fcns.putparam =  FORTRAN_NAME(putparam_k_3)  ;
     MU_A_fcns.putparam=  FORTRAN_NAME(putparam_mu_a) ;
     MU_S_fcns.putparam=  FORTRAN_NAME(putparam_mu_s) ;
     X_0_fcns.putparam =  FORTRAN_NAME(putparam_x_0)  ;
     Y_0_fcns.putparam =  FORTRAN_NAME(putparam_y_0)  ;
     Z_0_fcns.putparam =  FORTRAN_NAME(putparam_z_0)  ;
     POW_fcns.putparam =  FORTRAN_NAME(putparam_pow)  ;
 
     ERR_fcns.getparam =  NoParams                    ;    
     W_0_fcns.getparam =  FORTRAN_NAME(getparam_w_0)  ;    
     W_N_fcns.getparam =  FORTRAN_NAME(getparam_w_n)  ;    
     W_I_fcns.getparam =  FORTRAN_NAME(getparam_w_i)  ;    
     W_D_fcns.getparam =  FORTRAN_NAME(getparam_w_d)  ;    
     W_2_fcns.getparam =  FORTRAN_NAME(getparam_w_2)  ;    
     W_NI_fcns.getparam=  FORTRAN_NAME(getparam_w_ni)  ;    
     W_ID_fcns.getparam=  FORTRAN_NAME(getparam_w_id)  ;    
     K_0_fcns.getparam =  FORTRAN_NAME(getparam_k_0)  ;    
     K_1_fcns.getparam =  FORTRAN_NAME(getparam_k_1)  ;    
     K_2_fcns.getparam =  FORTRAN_NAME(getparam_k_2)  ;    
     K_3_fcns.getparam =  FORTRAN_NAME(getparam_k_3)  ;    
     MU_A_fcns.getparam=  FORTRAN_NAME(getparam_mu_a) ;    
     MU_S_fcns.getparam=  FORTRAN_NAME(getparam_mu_s) ;    
     X_0_fcns.getparam =  FORTRAN_NAME(getparam_x_0)  ;    
     Y_0_fcns.getparam =  FORTRAN_NAME(getparam_y_0)  ;    
     Z_0_fcns.getparam =  FORTRAN_NAME(getparam_z_0)  ;    
     POW_fcns.getparam =  FORTRAN_NAME(getparam_pow)  ;    

     ERR_fcns.getparam_lb =  getTAO_NINFINITY               ;    
     W_0_fcns.getparam_lb =  FORTRAN_NAME(getparam_w_0_lb)  ;    
     W_N_fcns.getparam_lb =  FORTRAN_NAME(getparam_w_n_lb)  ;    
     W_I_fcns.getparam_lb =  FORTRAN_NAME(getparam_w_i_lb)  ;    
     W_D_fcns.getparam_lb =  FORTRAN_NAME(getparam_w_d_lb)  ;    
     W_2_fcns.getparam_lb =  FORTRAN_NAME(getparam_w_2_lb)  ;    
     W_NI_fcns.getparam_lb=  FORTRAN_NAME(getparam_w_ni_lb) ;    
     W_ID_fcns.getparam_lb=  FORTRAN_NAME(getparam_w_id_lb) ;    
     K_0_fcns.getparam_lb =  FORTRAN_NAME(getparam_k_0_lb)  ;    
     K_1_fcns.getparam_lb =  FORTRAN_NAME(getparam_k_1_lb)  ;    
     K_2_fcns.getparam_lb =  FORTRAN_NAME(getparam_k_2_lb)  ;    
     K_3_fcns.getparam_lb =  FORTRAN_NAME(getparam_k_3_lb)  ;    
     MU_A_fcns.getparam_lb=  FORTRAN_NAME(getparam_mu_a_lb) ;    
     MU_S_fcns.getparam_lb=  FORTRAN_NAME(getparam_mu_s_lb) ;    
     X_0_fcns.getparam_lb =  FORTRAN_NAME(getparam_x_0_lb)  ;    
     Y_0_fcns.getparam_lb =  FORTRAN_NAME(getparam_y_0_lb)  ;    
     Z_0_fcns.getparam_lb =  FORTRAN_NAME(getparam_z_0_lb)  ;    
     POW_fcns.getparam_lb =  FORTRAN_NAME(getparam_pow_lb)  ;    

     ERR_fcns.getparam_ub =  getTAO_INFINITY                ;    
     W_0_fcns.getparam_ub =  FORTRAN_NAME(getparam_w_0_ub)  ;    
     W_N_fcns.getparam_ub =  FORTRAN_NAME(getparam_w_n_ub)  ;    
     W_I_fcns.getparam_ub =  FORTRAN_NAME(getparam_w_i_ub)  ;    
     W_D_fcns.getparam_ub =  FORTRAN_NAME(getparam_w_d_ub)  ;    
     W_2_fcns.getparam_ub =  FORTRAN_NAME(getparam_w_2_ub)  ;    
     W_NI_fcns.getparam_ub=  FORTRAN_NAME(getparam_w_ni_ub) ;    
     W_ID_fcns.getparam_ub=  FORTRAN_NAME(getparam_w_id_ub) ;    
     K_0_fcns.getparam_ub =  FORTRAN_NAME(getparam_k_0_ub)  ;    
     K_1_fcns.getparam_ub =  FORTRAN_NAME(getparam_k_1_ub)  ;    
     K_2_fcns.getparam_ub =  FORTRAN_NAME(getparam_k_2_ub)  ;    
     K_3_fcns.getparam_ub =  FORTRAN_NAME(getparam_k_3_ub)  ;    
     MU_A_fcns.getparam_ub=  FORTRAN_NAME(getparam_mu_a_ub) ;    
     MU_S_fcns.getparam_ub=  FORTRAN_NAME(getparam_mu_s_ub) ;    
     X_0_fcns.getparam_ub =  FORTRAN_NAME(getparam_x_0_ub)  ;    
     Y_0_fcns.getparam_ub =  FORTRAN_NAME(getparam_y_0_ub)  ;    
     Z_0_fcns.getparam_ub =  FORTRAN_NAME(getparam_z_0_ub)  ;    
     POW_fcns.getparam_ub =  FORTRAN_NAME(getparam_pow_ub)  ;    

     if( QOIInfo.pde.find("verifprob")!=std::string::npos )
      ERR_fcns.dpde_dm = FORTRAN_NAME(error_estimate_verif)  ;
     else
      ERR_fcns.dpde_dm = FORTRAN_NAME(error_estimate)        ;
      W_0_fcns.dpde_dm = FORTRAN_NAME(dwdw0)                 ;
      W_N_fcns.dpde_dm = FORTRAN_NAME(dwdwn)                 ;
      W_I_fcns.dpde_dm = FORTRAN_NAME(dwdwi)                 ;
      W_D_fcns.dpde_dm = FORTRAN_NAME(dwdwd)                 ;
      W_2_fcns.dpde_dm = FORTRAN_NAME(dwdw2)                 ;
     W_NI_fcns.dpde_dm = FORTRAN_NAME(dwdwni)                ;
     W_ID_fcns.dpde_dm = FORTRAN_NAME(dwdwid)                ;
      K_0_fcns.dpde_dm = FORTRAN_NAME(dkdk0)                 ;
      K_1_fcns.dpde_dm = FORTRAN_NAME(dkdk1)                 ;
      K_2_fcns.dpde_dm = FORTRAN_NAME(dkdk2)                 ;
      K_3_fcns.dpde_dm = FORTRAN_NAME(dkdk3)                 ;
     if( QOIInfo.pde.find("verifprob")!=std::string::npos
                               ||
         QOIInfo.pde.find("pennesisolasercooling")!=std::string::npos
                               ||
         QOIInfo.pde.find("nonlinpennesrf")!=std::string::npos
                               ||
         QOIInfo.pde.find("nonlinpennesisolaser")!=std::string::npos) 
       { 
        MU_A_fcns.dpde_dm = FORTRAN_NAME(dqlaserdmu_a)          ;
        MU_S_fcns.dpde_dm = FORTRAN_NAME(dqlaserdmu_s)          ;
         X_0_fcns.dpde_dm = FORTRAN_NAME(dqlaserdx)             ;
         Y_0_fcns.dpde_dm = FORTRAN_NAME(dqlaserdy)             ;
         Z_0_fcns.dpde_dm = FORTRAN_NAME(dqlaserdz)             ;
         POW_fcns.dpde_dm = FORTRAN_NAME(dqlaserdpow)           ;
        cout << "set isotropic laser functers \n";
       } 
     else if( QOIInfo.pde.find("nonlinpennesmonte")!=std::string::npos) 
       { 
        MU_A_fcns.dpde_dm = ReturnPetscScalarZero               ;
        MU_S_fcns.dpde_dm = ReturnPetscScalarZero               ;
         X_0_fcns.dpde_dm = FORTRAN_NAME(dqmontedx)             ;
         Y_0_fcns.dpde_dm = FORTRAN_NAME(dqmontedy)             ;
         Z_0_fcns.dpde_dm = FORTRAN_NAME(dqmontedz)             ;
         POW_fcns.dpde_dm = FORTRAN_NAME(dqmontedpow)           ;
        cout << "set monte carlo grid functers \n";
       } 
     else 
     {
        cout << "unknown pde  " << QOIInfo.pde << endl; abort();
     }

      ERR_fcns.qoi = ReturnPetscScalarZero             ;
      W_0_fcns.qoi = FORTRAN_NAME(regularization_w_0)  ;
      W_N_fcns.qoi = FORTRAN_NAME(regularization_w_n)  ;
      W_I_fcns.qoi = FORTRAN_NAME(regularization_w_i)  ;
      W_D_fcns.qoi = FORTRAN_NAME(regularization_w_d)  ;
      W_2_fcns.qoi = FORTRAN_NAME(regularization_w_2)  ;
     W_NI_fcns.qoi = FORTRAN_NAME(regularization_w_ni) ;
     W_ID_fcns.qoi = FORTRAN_NAME(regularization_w_id) ;
      K_0_fcns.qoi = FORTRAN_NAME(regularization_k_0)  ;
      K_1_fcns.qoi = FORTRAN_NAME(regularization_k_1)  ;
      K_2_fcns.qoi = FORTRAN_NAME(regularization_k_2)  ;
      K_3_fcns.qoi = FORTRAN_NAME(regularization_k_3)  ;
     MU_A_fcns.qoi = FORTRAN_NAME(regularization_mu_a) ;
     MU_S_fcns.qoi = FORTRAN_NAME(regularization_mu_s) ;
      X_0_fcns.qoi = FORTRAN_NAME(regularization_x_0)  ;
      Y_0_fcns.qoi = FORTRAN_NAME(regularization_y_0)  ;
      Z_0_fcns.qoi = FORTRAN_NAME(regularization_z_0)  ;
      POW_fcns.qoi = FORTRAN_NAME(regularization_pow)  ;

      ERR_fcns.dqoi_dm = ReturnPetscScalarZero;
      W_0_fcns.dqoi_dm = ReturnPetscScalarZero;//FORTRAN_NAME(regularizderiv_w_0) ;
      W_N_fcns.dqoi_dm = ReturnPetscScalarZero;//FORTRAN_NAME(regularizderiv_w_n) ;
      W_I_fcns.dqoi_dm = ReturnPetscScalarZero;//FORTRAN_NAME(regularizderiv_w_i) ;
      W_D_fcns.dqoi_dm = ReturnPetscScalarZero;//FORTRAN_NAME(regularizderiv_w_d) ;
      W_2_fcns.dqoi_dm = ReturnPetscScalarZero;//FORTRAN_NAME(regularizderiv_w_2) ;
     W_NI_fcns.dqoi_dm = ReturnPetscScalarZero;//FORTRAN_NAME(regularizderiv_w_ni);
     W_ID_fcns.dqoi_dm = ReturnPetscScalarZero;//FORTRAN_NAME(regularizderiv_w_id);
      K_0_fcns.dqoi_dm = ReturnPetscScalarZero;//FORTRAN_NAME(regularizderiv_k_0) ;
      K_1_fcns.dqoi_dm = ReturnPetscScalarZero;//FORTRAN_NAME(regularizderiv_k_1) ;
      K_2_fcns.dqoi_dm = ReturnPetscScalarZero;//FORTRAN_NAME(regularizderiv_k_2) ;
      K_3_fcns.dqoi_dm = ReturnPetscScalarZero;//FORTRAN_NAME(regularizderiv_k_3) ;
     MU_A_fcns.dqoi_dm = ReturnPetscScalarZero;//FORTRAN_NAME(regularizderiv_mu_a);
     MU_S_fcns.dqoi_dm = ReturnPetscScalarZero;//FORTRAN_NAME(regularizderiv_mu_s);
      X_0_fcns.dqoi_dm = ReturnPetscScalarZero;//FORTRAN_NAME(regularizderiv_x_0) ;
      Y_0_fcns.dqoi_dm = ReturnPetscScalarZero;//FORTRAN_NAME(regularizderiv_y_0) ;
      Z_0_fcns.dqoi_dm = ReturnPetscScalarZero;//FORTRAN_NAME(regularizderiv_z_0) ;
      POW_fcns.dqoi_dm = ReturnPetscScalarZero;//FORTRAN_NAME(regularizderiv_pow) ;

      ERR_fcns.spatial_field = PETSC_FALSE                   ;
      W_0_fcns.spatial_field = FORTRAN_NAME(get_w_0_field)() ;
      W_N_fcns.spatial_field = PETSC_FALSE                   ;
      W_I_fcns.spatial_field = PETSC_FALSE                   ;
      W_D_fcns.spatial_field = PETSC_FALSE                   ;
      W_2_fcns.spatial_field = PETSC_FALSE                   ;
     W_NI_fcns.spatial_field = PETSC_FALSE                   ;
     W_ID_fcns.spatial_field = PETSC_FALSE                   ;
      K_0_fcns.spatial_field = FORTRAN_NAME(get_k_0_field)() ;
      K_1_fcns.spatial_field = PETSC_FALSE                   ;
      K_2_fcns.spatial_field = PETSC_FALSE                   ;
      K_3_fcns.spatial_field = PETSC_FALSE                   ;
     MU_A_fcns.spatial_field = PETSC_FALSE                   ;
     MU_S_fcns.spatial_field = PETSC_FALSE                   ;
      X_0_fcns.spatial_field = PETSC_FALSE                   ;
      Y_0_fcns.spatial_field = PETSC_FALSE                   ;
      Z_0_fcns.spatial_field = PETSC_FALSE                   ;
      POW_fcns.spatial_field = PETSC_FALSE                   ;

     ERR_fcns.verifparam      = 0.0e0                                  ;
     W_0_fcns.verifparam      = 1.0e+1                                 ;
     W_N_fcns.verifparam      = 5.0571139e0                            ;
     W_I_fcns.verifparam      = 69.571139e0                            ;
     W_D_fcns.verifparam      = 1.3571139e0                            ;
     W_2_fcns.verifparam      = 2.758939e0                             ;
     W_NI_fcns.verifparam     = 3.452626e2                             ;
     W_ID_fcns.verifparam     = 3.452626e2                             ;
     K_0_fcns.verifparam      = 0.5e0                                  ;
     K_1_fcns.verifparam      = 0.04061452e0                           ;
     K_2_fcns.verifparam      = 0.035438703e0                          ;
     K_3_fcns.verifparam      = 315.314579715e0                        ;
     MU_A_fcns.verifparam     = 10.0e0                                 ;
     MU_S_fcns.verifparam     = 1400.0e0                               ;
     X_0_fcns.verifparam      = 0.04e0                                 ;
     Y_0_fcns.verifparam      = 0.06e0                                 ;
     Z_0_fcns.verifparam      = 0.07e0                                 ;
     POW_fcns.verifparam      = 10.0e0                                 ;

//     ERR_fcns.time_offset     = Given3IntsReturnPetscIntZero           ;
//     W_0_fcns.time_offset     = Given3IntsReturnPetscIntZero           ;
//     W_N_fcns.time_offset     = Given3IntsReturnPetscIntZero           ;
//     W_I_fcns.time_offset     = Given3IntsReturnPetscIntZero           ;
//     W_D_fcns.time_offset     = Given3IntsReturnPetscIntZero           ;
//     W_2_fcns.time_offset     = Given3IntsReturnPetscIntZero           ;
//     W_NI_fcns.time_offset    = Given3IntsReturnPetscIntZero           ;
//     W_ID_fcns.time_offset    = Given3IntsReturnPetscIntZero           ;
//     K_0_fcns.time_offset     = Given3IntsReturnPetscIntZero           ;
//     K_1_fcns.time_offset     = Given3IntsReturnPetscIntZero           ;
//     K_2_fcns.time_offset     = Given3IntsReturnPetscIntZero           ;
//     K_3_fcns.time_offset     = Given3IntsReturnPetscIntZero           ;
//     MU_A_fcns.time_offset    = Given3IntsReturnPetscIntZero           ;
//     MU_S_fcns.time_offset    = Given3IntsReturnPetscIntZero           ;
//     X_0_fcns.time_offset     = Given3IntsReturnPetscIntZero           ;
//     Y_0_fcns.time_offset     = Given3IntsReturnPetscIntZero           ;
//     Z_0_fcns.time_offset     = Given3IntsReturnPetscIntZero           ;
//     POW_fcns.time_offset     = time_offset_pow                        ;

      ERR_fcns.time_vary = PETSC_FALSE           ;
      W_0_fcns.time_vary = PETSC_FALSE           ;
      W_N_fcns.time_vary = PETSC_FALSE           ;
      W_I_fcns.time_vary = PETSC_FALSE           ;
      W_D_fcns.time_vary = PETSC_FALSE           ;
      W_2_fcns.time_vary = PETSC_FALSE           ;
     W_NI_fcns.time_vary = PETSC_FALSE           ;
     W_ID_fcns.time_vary = PETSC_FALSE           ;
      K_0_fcns.time_vary = PETSC_FALSE           ;
      K_1_fcns.time_vary = PETSC_FALSE           ;
      K_2_fcns.time_vary = PETSC_FALSE           ;
      K_3_fcns.time_vary = PETSC_FALSE           ;
     MU_A_fcns.time_vary = PETSC_FALSE           ;
     MU_S_fcns.time_vary = PETSC_FALSE           ;
      X_0_fcns.time_vary = PETSC_FALSE           ;
      Y_0_fcns.time_vary = PETSC_FALSE           ;
      Z_0_fcns.time_vary = PETSC_FALSE           ;
      POW_fcns.time_vary = PETSC_TRUE            ;

  PetscFunctionReturn(0);
}
/* ------------------------------------------------------------------- 
     Default Constructor for UniversalCntrl class
      minimal variables setup before debugger setup
   ------------------------------------------------------------------- */
UniversalCntrl::UniversalCntrl(GetPot  &controlfile)
{
  Control_Task  = PETSC_FALSE;
  Total_Nopt=0;  // initialize total number of optimization steps
  NDOF_FORWARD[0]=0; NDOF_FORWARD[1]=0; NDOF_FORWARD[2]=0;
  NDOF_ADJOINT[0]=0; NDOF_ADJOINT[1]=0; NDOF_ADJOINT[2]=0;
  NfieldElem =0;
  IGNORE_BC = PETSC_FALSE;
}

/* ------------------------------------------------------------------- 
     Setup the majority of the variables
     for this class after the debugger is setup
   ------------------------------------------------------------------- */
void UniversalCntrl::GenSetup(GetPot  &controlfile)
{

 PetscFunctionBegin; 

  // control of Initial Conditions
  IC_MRTI= PETSC_FALSE; 
  if(controlfile("init_cond/ic_mrti",false))
                                         IC_MRTI=PETSC_TRUE;// load IC from MRTI
  TIMELAG_IC=controlfile("init_cond/timelag_ic",0 );

  // get the restart ID. IDrestart == 0  ==> this is an initial run
  IDrestart = controlfile("compexec/restartid",0);
  FieldFile = controlfile("compexec/fieldfile","files/field.e");

  // default is not to compute the final objective function
  COMPUTEFINALOBJECTIVE = PETSC_FALSE;
  /* Control whether the final time qoi is computed. This is used mainly
     in parameter studies. This cannot be used in real-time b/c all thermal
     images must be available                                              */
  if(controlfile("output/computefinalobjective",false))
                         COMPUTEFINALOBJECTIVE = PETSC_TRUE;

  // control plot file write at checkpoint iterations
  PLOTITER = PETSC_FALSE; 
  if(controlfile("output/plotiter",false)) PLOTITER = PETSC_TRUE;

  // get the scaling for the Penalty term!!!!!!!!!!!!
  scalefact= controlfile("compexec/scalefact",1000.0);

  // check point parameters
  checkpointiter = controlfile("compexec/checkpointiter",500);
  /* set the check point counter */
  checkpoint = checkpointiter ;

  PetscFunctionReturnVoid(); 
}
void UniversalCntrl::printSelf()
{
 PetscFunctionBegin; 
     std::vector<PetscScalar>::iterator TimeIter;
  // for(TimeIter=GenInfo.TIME.begin();TimeIter!=GenInfo.TIME.end(); TimeIter++)
  //      printf( "dddas: TIME[%d]=%f\n", 
  //               distance(GenInfo.TIME.begin(),TimeIter),*TimeIter);
 std::cout << "\ndddas: COMPUTEFINALOBJECTIVE ="<<COMPUTEFINALOBJECTIVE
           << "\ndddas: PLOTITER              ="<<PLOTITER
           << "\ndddas: [compexec] scalefact  ="<<scalefact
           << "\ndddas: [compexec] checkpoint ="<<checkpoint
           << "\ndddas: [init_cond ]IC_MRTI   ="<<IC_MRTI
           << "\ndddas: [init_cond ]TIMELAG_IC="<<TIMELAG_IC
           << "\ndddas: checkpointiter        ="<<checkpointiter
           << "\ndddas: IDrestart             ="<<IDrestart 
           << "\ndddas: IGNORE_BC             ="<<IGNORE_BC
           << std::endl;
 PetscFunctionReturnVoid(); 
}

