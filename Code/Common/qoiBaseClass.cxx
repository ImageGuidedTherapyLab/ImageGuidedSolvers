// C++ include files 
#include <iostream>
#include <fstream>
#include <vector>

// libMesh include files
#include "libmesh.h"
#include "libmesh_common.h"
#include "mesh.h"
#include "equation_systems.h"
#include "fe.h"
#include "fe_interface.h"
#include "dof_map.h"
#include "quadrature_gauss.h"
#include "getpot.h"
#include "parallel.h" // file parser
#include "petsc_vector.h"
#include "time_solver.h"
#include "petsc_matrix.h"
#include "petsc_diff_solver.h"

// tao interface
#include "src/tao_impl.h" 

// dddas include files
#include "applicationContext.h"
#include "pdeBaseClass.h"
#include "thermal_therapy_system.h"
#include "qoiBaseClass.h"
#include "BackgroundPhase.h"
#include "dddas.h"
#include "parser_defaults.h"

/* ------------------------------------------------------------------- 
     Constructor for qoiBaseClass class
      (constructor is minimal to use petsc error codes in the main setup)
   ------------------------------------------------------------------- */
qoiBaseClass::qoiBaseClass(AppSolve *user,GetPot &controlfile,PetscInt iii)
{
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  fncEvalCnt = 0;   // counter
  Penalty=PETSC_TRUE; // default is to use penalty function 
  NDOF_CONTROL[0] = 0;
  NDOF_CONTROL[1] = 0;
  NDOF_CONTROL[2] = 0;
  ComputeQoiData   = NULL;
  PostSolveQoiData = NULL;
  Init_Nsteplo = 0;
  tao    = NULL;  
  taoapp = NULL;  
  this->m_BackgroundCorrection = false; 
  QCONTRL       = NULL;
  CNTL_LOC      = NULL;
  GLOCSCAT_CNTL = NULL;

  // resize
  TIMEERR.resize(std::max( user->get_max_ideal() +1,0),0.0);
  MRTI_MEM.resize( user->get_num_MRTI() +1 ,0);

  std::ostringstream qoiID;
  qoiID << "qoi_"<< iii << "/" ;
  controlfile.set_prefix(qoiID.str().c_str());
  // time lag for the computations to finish
  TIMELAG= controlfile("timelag",0);

  TAO_MONITOR  = PETSC_FALSE; 

  // plot file write control
  PLOTOPTSOLVE = PETSC_FALSE; 
  // if need to compute the final value of the objective function
  // or want to plot the result of the optimization solve set to true
  if(user->ComputeFinalObjective() || controlfile("plotoptsolve",false))
                                                      PLOTOPTSOLVE = PETSC_TRUE;
  /* Get optimization problem information for this Computational Group*/

  ////@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ////@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ////@@ Quantity of Interest                                             
  ////@@                                                                  
  ////@@     || * || is understood as the space-time norm                 
  ////@@                                                                  
  ////@@   QOI(u) .eq.    1/2  *||\chi_MRI(x) {u-u_ideal}||^2             
  ////@@                +   1/2  *||\chi_MRI(x) {HSP(u)-HSP_ideal}||^2    
  ////@@                  +   1/2  *||\chi_MRI(x) {D(u)-D_ideal}||^2      
  ////@@                    +   1/2  *|| w_0(x) ||^2                      
  ////@@                      +   1/2  *|| \nabla  w_0(x) ||^2            
  ////@@                        +   1/2  *|| \chi_goal(x)  u ||^2         
  ////@@                          +   1/2  *||\chi_MRI(x) {u-u_MRTI}||^2  
  ////@@                                                                  
  ////@@     compobj .eq. qoitempbased -->  temp based optimization        
  ////@@     compobj .eq. qoihspbased  -->  hsp based optimization         
  ////@@     compobj .eq. qoidambased  -->  damage based optimization      
  ////@@     compobj .eq.      n/a     -->  L_2 regularization             
  ////@@     compobj .eq.      n/a     -->  H_1 regularization             
  ////@@     compobj .eq.      n/a     -->  goal-oriented error estimation 
  ////@@     compobj .eq.      n/a     -->  calibration                    
  ////@@                                                                  
  ////@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ////@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  //// default is not to check for thermal image availability
  //UpdateQoiData= ReturnPETSC_FALSE; 
  //// default is not to post process
  //PostSolveQoiData= NULL;

  //// if the images will be used as IC's then the compute task will
  //// check when the next image is ready... this is only useful for more
  //// than one optimization step as the 1st optimization step does not
  //// use and thermal images unless for calibration
  //if( NoptSteps > 1  && GenInfo.IC_MRTI) UpdateQoiData=checkMRTI_data4_IC; 
  //if (compobj.find("hsp_control")!=std::string::npos)
  //  {  // hsp based optimization 
  //     ComputeObjective= Lesbegue2_QOI;
  //     PriorLoadQoiData= NULL;
  //     ComputeQoiData= LoadIdealHSP;
  //  }
  //else if (compobj.find("dam_control_arr")!=std::string::npos)
  //  {
  //     //Arrhenius damage based optimization
  //     ComputeObjective= Lesbegue2_QOI;
  //     PriorLoadQoiData= NULL;
  //     ComputeQoiData  = LoadIdealArr;
  //     PostSolveQoiData= CurrentIdealArr;
  //  } 
  //else if (compobj.find("dam_control_two")!=std::string::npos)
  //  {
  //     //two state model damage based optimization
  //     ComputeObjective= Lesbegue2_QOI;
  //     PriorLoadQoiData= NULL;
  //     ComputeQoiData= LoadIdeal_TS;
  //  }  
  //else if (compobj.find("calibration")!=std::string::npos)
  //  {
  //     //calibration to MRTI data 
  //     ComputeObjective= SpaceTime_QOI;
  //     PriorLoadQoiData= readMRTI_data ;
  //     ComputeQoiData= NULL;
  //     // the calibration requires a check of images further ahead in time
  //     UpdateQoiData= checkMRTI_data4OPT; 
  //  }
  //else if (compobj.find("computeideal_mc")!=std::string::npos)
  //  {
  //     //compute ideal field from Monte Carlo data
  //     ComputeObjective= SpaceTime_QOI;
  //     // loading an ideal field of zero will be used to compute the space-time
  //     // norm of the temperature field only
  //     //PriorLoadQoiData= computeideal_mc;
  //     PriorLoadQoiData= LoadIdealTmp; // 
  //     ComputeQoiData= NULL;
  //  } 
  //else if (compobj.find("verifprob")!=std::string::npos)
  //  { 
  //     ComputeObjective= SpaceTime_QOI;
  //     switch(AppSolve::NEXACT)
  //       {
  //         case 15://Verification problems that do require MRTI data transfer
  //         PriorLoadQoiData  = readMRTI_data;
  //         break;
  //         default://Most Verification problems DO NOT require MRTI data transfer
  //         PriorLoadQoiData  = verifcomputeidealfield;
  //         //Cannot precompute ideal field if nothing setup
  //         GenInfo.COMPUTEFINALOBJECTIVE = PETSC_FALSE;
  //         break;
  //       }
  //     ComputeQoiData= NULL;
  //  } 
  //else if (compobj.find("temp_control")!=std::string::npos)
  //  { 
  //     ComputeObjective= SpaceTime_QOI;
  //     PriorLoadQoiData= LoadIdealTmp;
  //     ComputeQoiData= NULL;
  //     // NOTE: for goal oriented space-time L2 error estimate set the ideal
  //     //       field to zero
  //  } 
  //else if (compobj.find("notfound")!=std::string::npos)
  //  { 
  //     cout << "an objective function was not found"<< endl; abort();
  //  } 
  //else 
  //  {
  //     cout << "unknown objective function "<<compobj << endl; abort();
  //  } 


  // reset default search path
  controlfile.set_prefix("");
  //// get the function pointers to the QOI data
  //pde = controlfile("hp3d/pde","notfound");

  //// perturbations for Finite difference verification
  //ERR_fcns.perturbation    =          1.0e0                                  ;
  //K_0_fcns.perturbation    = controlfile("thermalcond/dk_0"  ,1.0e-2);
  //K_1_fcns.perturbation    = controlfile("thermalcond/dk_1"  ,2.5e-2);
  //K_2_fcns.perturbation    = controlfile("thermalcond/dk_2"  ,2.5e-2);
  //K_3_fcns.perturbation    = controlfile("thermalcond/dk_3"  ,2.5e-2);
  //W_0_fcns.perturbation    = controlfile("perfusivity/dw_0"  ,1.0e-5);
  //W_N_fcns.perturbation    = controlfile("perfusivity/dw_n"  ,1.0e-5);
  //W_I_fcns.perturbation    = controlfile("perfusivity/dw_i"  ,1.0e-5);
  //W_D_fcns.perturbation    = controlfile("perfusivity/dw_d"  ,1.0e-5);
  //W_2_fcns.perturbation    = controlfile("perfusivity/dw_2"  ,1.0e-5);
  //W_NI_fcns.perturbation   = controlfile("perfusivity/dw_ni" ,1.0e-5);
  //W_ID_fcns.perturbation   = controlfile("perfusivity/dw_id" ,1.0e-5);
  //MU_A_fcns.perturbation   = controlfile("source_laser/dmu_a",1.0e-2);
  //MU_S_fcns.perturbation   = controlfile("source_laser/dmu_s",1.0e-2);
  //POW_fcns.perturbation    = controlfile("source_laser/dpow" ,1.0e-2);
  //if( pde.find("verifprob")!=std::string::npos ||
  //    pde.find("pennesisolasercooling")!=std::string::npos ||
  //    pde.find("nonlinpennesisolaser")!=std::string::npos) { 
  //   X_0_fcns.perturbation = controlfile("source_laser/dx_0" ,1.0e-4);
  //   Y_0_fcns.perturbation = controlfile("source_laser/dy_0" ,1.0e-4);
  //   Z_0_fcns.perturbation = controlfile("source_laser/dz_0" ,1.0e-4);
  //} else if(pde.find("nonlinpennesmonte")!=std::string::npos) { 
  //   X_0_fcns.perturbation = controlfile("source_laser/dx_0" ,8.58e-2);
  //   Y_0_fcns.perturbation = controlfile("source_laser/dy_0" ,1.98e-2);
  //   Z_0_fcns.perturbation = controlfile("source_laser/dz_0" ,1.98e-2);
  //} else if(pde.find("nonlinpennesrf")!=std::string::npos) { 
  //   X_0_fcns.perturbation = controlfile("source_laser/dx_0" ,1.0e-4);
  //   Y_0_fcns.perturbation = controlfile("source_laser/dy_0" ,1.0e-4);
  //   Z_0_fcns.perturbation = controlfile("source_laser/dz_0" ,1.0e-4);
  //} else if(pde.find("notfound")!=std::string::npos) { 
  //   cout << "a pde type was not found"<< endl; abort();
  //} else {
  //   cout << "unknown pde "<<pde << endl; abort();
  //}

  minHessianColumn = 0; // hack for computing hessian subset
  maxHessianColumn =-1; // hack for computing hessian subset

}
/* -------------------------------------------------------------------
     initialize imaging
   ------------------------------------------------------------------- */
void qoiBaseClass::init_imaging(PetscInt qoiMethod,GetPot &controlfile)
{
  PetscFunctionBegin; 
  // input control
  ImageAcquisitionType inputmethod =  ImageAcquisitionType(
                       controlfile("mrti/inputmethod",(int)NO_IMAGING)
                                                          );
  switch(inputmethod)
   {
    case NO_IMAGING: // no imaging
      this->Images = NULL;
      // error check
      if( qoiMethod == 1 || qoiMethod == 4 ) 
        {
          std::cout << "expecting some form of imaging data..."
                    << std::endl << std::flush; abort();  
        }
      break;
    case DICOM2D: // 2d dicom series return snr and temperature as raw data
      this->Images = new ImagingComplex(controlfile);
      std::cout << "conventional complex difference snr temperature maps..."
                << std::endl ;  
      break;
    case BACKGROUNDPHASE: // return snr and phase as raw data
     {
      this->Images = new ImagingBackgroundPhase(controlfile);
      AppSolve *user = AppSolve::getInstance();
      EquationSystems   &EqnSystem = user->get_equation_systems();
      BackgroundPhaseSystem< BackgroundPhaseModel > & background_system =
       EqnSystem.add_system< BackgroundPhaseSystem < BackgroundPhaseModel > >("BackgroundSystem");
      std::cout << "setting up imaging to run w/ " << background_system.name()
                << std::endl ;  

      this->m_BackgroundCorrection = true; 
     }
      break;
    case VOLUME3D:  // 3d Volume image
      this->Images = new ImagingPreProcessed(controlfile);
      std::cout << "expecting preprocess temperature/snr maps..."
                << std::endl ;  
      break;
    case EIGHT_ECHO_CSI:   // CSI w/ 8 echos
      this->Images = new ImagingChemicalShift< 8 >(controlfile);    
      break;
    case SIXTEEN_ECHO_CSI: // CSI w/ 16 echos
      this->Images = new ImagingChemicalShift< 16 >(controlfile);
      break;
    default: 
      std::cout << "unknown inputmethod "<< inputmethod 
                << std::endl << std::flush; abort();
   }
 // get basic header info
 if(this->Images) this->Images->GetHeaderData(inputmethod);
 PetscFunctionReturnVoid();
}
/* Function to get MRTI uncertainty */
Number project_one (const Point& ,
                    const Parameters& ,
                    const std::string& ,
                    const std::string& )
{
  return 1.0e0;
}
/* -------------------------------------------------------------------
     get iteration number
   ------------------------------------------------------------------- */
PetscInt qoiBaseClass::getIteration() { return tao->iter; }
/* -------------------------------------------------------------------
     echo data
   ------------------------------------------------------------------- */
void qoiBaseClass::printSelf(std::ostream& os)
{
  os << "qoiBaseClass: NDOF_CONTROL[0] =" <<  NDOF_CONTROL[0]  << std::endl;
  os << "qoiBaseClass: NDOF_CONTROL[1] =" <<  NDOF_CONTROL[1]  << std::endl;
  os << "qoiBaseClass: NDOF_CONTROL[2] =" <<  NDOF_CONTROL[2]  << std::endl;
  os << "qoiBaseClass: TIMELAG         =" <<  TIMELAG          << std::endl;
  os << "qoiBaseClass: TAO_MONITOR     =" <<  TAO_MONITOR      << std::endl;
  os << "qoiBaseClass: PLOTOPTSOLVE    =" <<  PLOTOPTSOLVE     << std::endl;
  os << "qoiBaseClass: Penalty         =" <<  Penalty          << std::endl;
  os << "qoiBaseClass: GetParamSize()  =" <<  GetParamSize()   << std::endl;
  os << "qoiBaseClass: minHessianColumn=" <<  minHessianColumn << std::endl;
  os << "qoiBaseClass: maxHessianColumn=" <<  maxHessianColumn << std::endl;
}
void qoiBaseClass::echoOptWindow(std::ostream& os,
                                 PetscInt Nsteplo, PetscInt Nstephi)
{
  if(TAO_MONITOR && !rank )
    {
     AppSolve *user = AppSolve::getInstance();
     os << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << std::endl
        << " IDEAL_NZERO= " <<  user->IdealNzero()
        << " IDEAL_NTIME= " <<  user->IdealNtime() <<  std::endl       
        << " Nsteplo= "     <<  Nsteplo
        << " Nstephi= "     <<  Nstephi     << std::endl << std::endl ;
    }
}

/* ------------------------------------------------------------------- 
     Get optimization parameters from Hp3d 
   ------------------------------------------------------------------- */

#undef __FUNCT__
#define __FUNCT__ "qoiBaseClass::physicalBounds"
PetscErrorCode qoiBaseClass::physicalBounds(libMesh::MeshBase &mesh,Vec XL,Vec XU)
{

  AppSolve *user = AppSolve::getInstance();
  PetscErrorCode  info;
  // used for error handling
  PetscFunctionBegin; 

  std::cout << "physicalBounds not ready " << std::endl;
  libmesh_error();
  //// get pointer to variance data
  //PetscScalar    *lowerBoundData,*upperBoundData;
  //info = VecGetArray(XL,&lowerBoundData);CHKERRQ(info);
  //info = VecGetArray(XU,&upperBoundData);CHKERRQ(info);

  //// should have same owner range
  //PetscInt ilo,ihi;
  //info = VecGetOwnershipRange(XL,&ilo,&ihi); CHKERRQ(info);

  //// counter used for storage
  //for(IterParam=Parameters.begin(); IterParam!=Parameters.end(); IterParam++)
  //  {
  //   optimizationParameter* optParam = *IterParam;
  //   const int paramID   = distance(Parameters.begin(),IterParam);

  //   // loop over the elements 
  //   libMesh::MeshBase::const_element_iterator       it  = mesh.active_local_elements_begin();
  //   const libMesh::MeshBase::const_element_iterator end = mesh.active_local_elements_end(); 
  //   for ( ; it != end; ++it)
  //     {// element iterator
  //      const Elem* elem = *it;

  //      const int globalDof = elem->_param_opt_map.at(paramID);

  //      // parameter bounds are stored by block number
  //      int idBlock   = elem->subdomain_id() - 1;

  //      // not all procs may own part of the global vec
  //      if (globalDof >= ilo && globalDof < ihi) 
  //       {// plot all params
  //        optParam->getParam_lb(
  //           lowerBoundData[ globalDof - ilo ],idBlock);
  //        optParam->getParam_ub(
  //           upperBoundData[ globalDof - ilo ],idBlock);
  //       }
  //     }
  //   //  Possibility exists that a constant parameter is locally 
  //   //  owned in the control vector but the corresponding elements within the
  //   //  the entire block are on another procs domain. 
  //   if( !libMesh::processor_id() )
  //    for( PetscInt Ii = 0 ; Ii < user->get_num_elem_blk() ; Ii++ )
  //     {
  //      optParam->getParam_lb(lowerBoundData[user->get_num_elem_blk()*paramID+Ii],Ii);
  //      optParam->getParam_ub(upperBoundData[user->get_num_elem_blk()*paramID+Ii],Ii);
  //     }
  //  }	 

  //// restore pointers
  //info = VecRestoreArray(XL,&lowerBoundData);CHKERRQ(info);
  //info = VecRestoreArray(XU,&upperBoundData);CHKERRQ(info);

  PetscFunctionReturn(0);
}
/* ------------------------------------------------------------------- 
     Get optimization parameters 
   ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "qoiBaseClass::GetCntrlVars"
PetscErrorCode qoiBaseClass::GetCntrlVars(libMesh::MeshBase &mesh,Vec Xinit)
{
  PetscErrorCode  info; // used for error handling
  AppSolve *user = AppSolve::getInstance();
  PetscFunctionBegin; 

  std::cout << "GetCntrlVars not ready " << std::endl;
  libmesh_error();
  //// get pointer to variance data
  //PetscScalar    *initialData;
  //info = VecGetArray(Xinit,&initialData);CHKERRQ(info);

  //// get global range 
  //PetscInt ilo,ihi;
  //info = VecGetOwnershipRange(Xinit,&ilo,&ihi); CHKERRQ(info);

  //// do not loop over the error parameter
  //IterParam=Parameters.begin(); IterParam++;
  //for( ; IterParam!=Parameters.end(); IterParam++)
  //  {// get paramter pointer
  //   optimizationParameter* optParam = *IterParam;
  //   const int paramID   = distance(Parameters.begin(),IterParam);

  //   // loop over local elements to initialize all field parameters
  //   libMesh::MeshBase::const_element_iterator       it  = mesh.active_local_elements_begin();
  //   const libMesh::MeshBase::const_element_iterator end = mesh.active_local_elements_end(); 
  //   for ( ; it != end; ++it)
  //     {// element iterator
  //      const Elem* elem = *it;

  //      const int globalDof = elem->_param_opt_map.at(paramID);

  //      // parameter bounds are stored by block number
  //      int idBlock   = elem->subdomain_id() - 1;

  //      // not all procs may own part of the global vec
  //      if (globalDof >= ilo && globalDof < ihi) 
  //       {// get params
  //        optParam->getParam( initialData[ globalDof - ilo ] , idBlock );
  //       }
  //     }
  //   //  Possibility exists that a constant parameter is locally 
  //   //  owned in the control vector but the corresponding elements within the
  //   //  the entire block are on another procs domain. 
  //   if( !libMesh::processor_id() )
  //    for( PetscInt Ii = 0 ; Ii < user->get_num_elem_blk() ; Ii++ )
  //       optParam->getParam( initialData[ user->get_num_elem_blk()*paramID+Ii ] , Ii );
  //  }	 
  //
  //// restore pointers
  //info = VecRestoreArray(Xinit,&initialData);CHKERRQ(info);
  
  PetscFunctionReturn(0);
}
/* ------------------------------------------------------------------- 
     Get optimization parameters 
   ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "qoiBaseClass::GetLocalCntrlVars"
PetscErrorCode qoiBaseClass::GetLocalCntrlVars(Vec Xloc)
{

  PetscScalar *parambuffer;
  PetscErrorCode  info;

  // used for error handling
  PetscFunctionBegin; 

  // get pointer to local data buffer
  info = VecGetArray(Xloc,&parambuffer); CHKERRQ(info);

  // initialize counter
  PetscInt dofcount= 0;

  for(IterParam=Parameters.begin(); IterParam!=Parameters.end(); IterParam++)
    {
     optimizationParameter* optParam = *IterParam;

     // get the control variables 
     for(DofIter = optParam->dofs.begin(); DofIter != optParam->dofs.end();
                                                                     DofIter++)
      {
         unsigned int  dofID = distance(optParam->dofs.begin(),DofIter);
         optParam->getParam(parambuffer[dofcount+dofID],DofIter[0]);
      }

     // update counter
     dofcount = dofcount + optParam->dofs.size();
    }

  // restore pointer
  info = VecRestoreArray(Xloc,&parambuffer); CHKERRQ(info);
  
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- 
     Get optimization parameters from Hp3d 
   ------------------------------------------------------------------- */

#undef __FUNCT__
#define __FUNCT__ "qoiBaseClass::GetLocalVars_LB"
PetscErrorCode qoiBaseClass::GetLocalVars_LB(Vec Xloc) // deprecated
{

  PetscScalar *parambuffer;
  PetscErrorCode  info;

  // used for error handling
  PetscFunctionBegin; 

  // get pointer to local data buffer
  info = VecGetArray(Xloc,&parambuffer); CHKERRQ(info);

  // initialize counter
  PetscInt dofcount= 0;

  for(IterParam=Parameters.begin(); IterParam!=Parameters.end(); IterParam++)
   {

     optimizationParameter* optParam = *IterParam;
     // get the control variables from hp3d
     for(DofIter = optParam->dofs.begin(); DofIter != optParam->dofs.end();
                                                                     DofIter++)
      {
         unsigned int dofID = distance(optParam->dofs.begin(),DofIter);
         optParam->getParam_lb(parambuffer[dofcount+dofID],DofIter[0]);
      }

     // update counter
     dofcount = dofcount + optParam->dofs.size();
   }

  // restore pointer
  info = VecRestoreArray(Xloc,&parambuffer); CHKERRQ(info);
  
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- 
     Get optimization parameters from Hp3d 
   ------------------------------------------------------------------- */

#undef __FUNCT__
#define __FUNCT__ "qoiBaseClass::GetLocalVars_UB"
PetscErrorCode qoiBaseClass::GetLocalVars_UB(Vec Xloc) // deprecated
{

  PetscScalar *parambuffer;
  PetscErrorCode  info;

  // used for error handling
  PetscFunctionBegin; 

  // get pointer to local data buffer
  info = VecGetArray(Xloc,&parambuffer); CHKERRQ(info);

  // initialize counter
  PetscInt dofcount= 0;

  for(IterParam=Parameters.begin(); IterParam!=Parameters.end(); IterParam++)
   {

     optimizationParameter* optParam = *IterParam;
     // get the control variables from hp3d
     for(DofIter = optParam->dofs.begin(); DofIter != optParam->dofs.end();
                                                                     DofIter++)
      {
         unsigned int dofID = distance(optParam->dofs.begin(),DofIter);
         optParam->getParam_ub(parambuffer[dofcount+dofID],DofIter[0]);
      }

     // update counter
     dofcount = dofcount + optParam->dofs.size();
   }

  // restore pointer
  info = VecRestoreArray(Xloc,&parambuffer); CHKERRQ(info);
  
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- 
     Put optimization parameters into Hp3d 
   ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "qoiBaseClass::PutLocalCntrlVars"
PetscErrorCode qoiBaseClass::PutLocalCntrlVars(Vec Xloc)
{
  PetscScalar *parambuffer;
  AppSolve *user = AppSolve::getInstance();
  PetscErrorCode  info;

  // used for error handling
  PetscFunctionBegin; 

  // get pointer to local data buffer
  info = VecGetArray(Xloc,&parambuffer); CHKERRQ(info);

  // initialize counter
  PetscInt dofcount= 0;

  for(IterParam=Parameters.begin(); IterParam!=Parameters.end(); IterParam++)
   {
     optimizationParameter* optParam = *IterParam;
     // put the control variables into data structures
     for(DofIter = optParam->dofs.begin(); DofIter != optParam->dofs.end();
                                                                    DofIter++)
      {
         unsigned int dofID = distance(optParam->dofs.begin(),DofIter);
         optParam->putParam(parambuffer[dofcount+dofID],DofIter[0]);
      }
     // domain constant field parameter need to have the same paramters in 
     // the field section as the constant section 
     if( !optParam->spatial_field && !optParam->time_vary )
      { 
         for( PetscInt iii = user->get_num_elem_blk(); iii < optParam->size() ; iii++ )
           optParam->putParam(parambuffer[dofcount],iii);
      }
    
     // update counter
     dofcount = dofcount + optParam->dofs.size();
   }

  // restore pointer
  info = VecRestoreArray(Xloc,&parambuffer); CHKERRQ(info);
  
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- 
     Add regularization terms to objective function 
   ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "qoiBaseClass::AddRegularization"
PetscErrorCode qoiBaseClass::AddRegularization(PetscScalar &Qoi_loc,
                                                     PetscScalar &ScaleFact)
{

  // used for error handling
  PetscFunctionBegin; 

  for(IterParam=Parameters.begin(); IterParam!=Parameters.end(); IterParam++)
   {
     optimizationParameter* optParam = *IterParam;
     if(optParam->time_vary)
      { // varies in time
        if(!rank) {
           for( DofIter = optParam->dofs.begin(); 
                          DofIter != optParam->dofs.end() ; DofIter++ ){
             Qoi_loc=Qoi_loc+optParam->qoi(DofIter[0],ScaleFact);
           }
        }
      }
     else
      { //may vary spatially
        // to ensure that regularization terms and not added by every task 
        //    for the constant part of the fields. The constant portion of
        //    the regularization term is added only by rank 0 
        DofIter = optParam->dofs.begin(); 
        if(!rank) 
          Qoi_loc = Qoi_loc+optParam->qoi(DofIter[0],ScaleFact);
        DofIter++;

        while(DofIter != optParam->dofs.end()){
            // add the regularization term for this dof
            Qoi_loc = Qoi_loc + 
                        optParam->qoi(DofIter[0],ScaleFact);
            DofIter++;
        }
      }

   }

  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- 
     add derivatives of the regularization term to the gradient
#undef __FUNCT__
#define __FUNCT__ "qoiBaseClass::RegularizationDeriv"
PetscErrorCode qoiBaseClass::RegularizationDeriv(PetscScalar *Grad_loc,
                                              PetscScalar *ScaleFact)
{

  // used for error handling
  PetscFunctionBegin; 

  // initialize counter
  PetscInt dofcount= 0;

  for(IterParam=Parameters.begin(); IterParam!=Parameters.end(); IterParam++)
   {
    if(IterParam->time_vary)
     { //varies in time
       if(!rank) {
          for( DofIter = IterParam->dofs.begin(); 
                         DofIter != IterParam->dofs.end() ; DofIter++ ){
             dofID = distance(IterParam->dofs.begin(),DofIter);
             Grad_loc[dofcount+dofID] = Grad_loc[dofcount+dofID] + 
                              IterParam->regularizderiv(&DofIter[0],ScaleFact);
          }
       }
     }
    else
     { //may vary spatially
       // only rank 0 adds to the constant term
       DofIter = IterParam->dofs.begin(); 
       if(!rank) Grad_loc[dofcount]= Grad_loc[dofcount] + 
                          IterParam->regularizderiv(&DofIter[0],ScaleFact);
       DofIter++;

       // all tasks add to the field regularizations terms 
       while(DofIter != IterParam->dofs.end()){
          dofID = distance(IterParam->dofs.begin(),DofIter);
          Grad_loc[dofcount+dofID] = Grad_loc[dofcount+dofID] +
                          IterParam->regularizderiv(&DofIter[0],ScaleFact);
          DofIter++;
       }
     }

    // update counter
    dofcount = dofcount + IterParam->dofs.size();
   }

  PetscFunctionReturn(0);
}
   ------------------------------------------------------------------- */

/* ------------------------------------------------------------------- 
     Used to put parameters into Control vector for  Verification problems
   ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "qoiBaseClass::PutVerifCntrlVars"
PetscErrorCode qoiBaseClass::PutVerifCntrlVars()
{
  PetscScalar param;
  // used for error handling
  PetscFunctionBegin; 

  for(IterParam=Parameters.begin(); IterParam!=Parameters.end(); IterParam++)
   {

     optimizationParameter* optParam = *IterParam;
     // put the control variables into hp3d
     for(DofIter = optParam->dofs.begin(); DofIter != optParam->dofs.end();
                                                                     DofIter++)
      {
         unsigned int dofID = distance(optParam->dofs.begin(),DofIter);
         param = optParam->_verifparam * (dofID+1) ; 
         optParam->putParam(param,DofIter[0]);
      }

   }

  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- 
     Used to perturb parameters of the control vector for 
     FD computations Verification problems
   ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "qoiBaseClass::GetVerifDeltaCntrl"
PetscErrorCode qoiBaseClass::GetVerifDeltaCntrl(Vec Xloc,PetscInt IDdof)
{
  PetscInt dofcount,  locdof = -1; 
  PetscErrorCode info;

  // used for error handling
  PetscFunctionBegin; 


  // find the corresponding local dof for the global dof IDdof
  for( PetscInt iii = 0 ; iii < NDOF_CONTROL[2]; iii++){
      // note that multiple processors will OWN the same dof
      // HOWEVER using INSERT_VALUES for the VecScatter
      // will make the final global vector have the correct value
      if(locmap.at(iii) == IDdof) locdof = iii;
  }

  if(locdof != -1){ // this task owns this global dof
     // initialize counter
     dofcount=0;

     // find the range where IDdof belongs
     IterParam=Parameters.begin();
     optimizationParameter* optParam = *IterParam;

     dofcount=optParam->dofs.size(); 
     while(locdof >= dofcount)
      {
       IterParam++; 

       optParam = *IterParam; // FIXME is this needed?

       dofcount+=optParam->dofs.size(); 
      }

     // put the perturbation into the local vector
     VecSetValue(Xloc,locdof,optParam->_perturbation,INSERT_VALUES);

     info = VecAssemblyBegin(Xloc); CHKERRQ(info);
     info = VecAssemblyEnd(Xloc); CHKERRQ(info);
  }

  PetscFunctionReturn(0);
}

/* -------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "qoiBaseClass::getGlobalMap"
void qoiBaseClass::getGlobalMap(const Elem* elem,
                                std::vector<PetscInt> &globmap)
{
 AppSolve *user = AppSolve::getInstance();
 PetscFunctionBegin;

  std::cout << "getGlobalMap not ready " << std::endl;
  libmesh_error();
 ////get the element wise map to the global dofs 
 //for( IterParam=Parameters.begin(); IterParam!=Parameters.end(); IterParam++)
 // {
 //  int iii = distance(Parameters.begin(),IterParam);
 //  optimizationParameter* optParam = *IterParam;
 //  // global map ONLY altered if the parameter varies with time
 //  // the mapping uses the floor property of integer division which is
 //  // the reason for the " ISTEP - 1 "
 //  globmap.at(iii)  =elem->_param_opt_map.at(iii) ;
 //  globmap.at(iii) +=  (optParam->time_vary ? 1 : 0) * ( (AppSolve::ISTEP - 1) /
 //                       user->IstepsPerIdeal() - user->IdealNzero() );
 // }

 PetscFunctionReturnVoid();
}
/* -------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "qoiBaseClass::getGlobalParamValues"
void qoiBaseClass::getGlobalParamValues(const Elem* elem,
                                std::vector<PetscScalar> &elemParameter)
{
 PetscFunctionBegin;

 //get the element wise map to the global dofs 
 std::vector<PetscInt> globmap(this->GetParamSize(),0);
 this->getGlobalMap(elem,globmap);

 //get local map 
 std::vector<PetscInt> locdof;
 for( IterParam=Parameters.begin(); IterParam!=Parameters.end(); IterParam++)
  {
   int iii = distance(Parameters.begin(),IterParam);
   for( PetscInt jjj = 0 ; jjj < this->NDOF_CONTROL[2]; jjj++)
    {
       if(locmap.at(jjj) == globmap.at(iii)) locdof.push_back(jjj);
    }
  }

 //get values 
 PetscErrorCode info = VecGetValues(this->CNTL_LOC, locdof.size(), 
                                &locdof[0], &elemParameter[0] ); CHKERRV(info);

 PetscFunctionReturnVoid();
}
/* -------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "qoiBaseClass::getParameterPointer"
std::vector<optimizationParameter*>::iterator 
                       qoiBaseClass::getParameterPointer(const int &globalDof)
{
 PetscFunctionBegin;

 // find the corresponding local dof for the global dof IDdof
 PetscInt locdof = -1; 
 for( PetscInt jjj = 0 ; jjj < this->NDOF_CONTROL[2]; jjj++)
  {
     if(locmap.at(jjj) == globalDof) locdof = jjj;
  }

 // find the range where IDdof belongs
 if(locdof != -1) // this task owns this global dof
  {
   IterParam=Parameters.begin(); 
   optimizationParameter* optParam = *IterParam;

   // initialize counter and loop through parameters
   PetscInt dofcount = optParam->dofs.size(); 
   while(locdof >= dofcount)
    {
     IterParam++; 
     libmesh_assert( IterParam != Parameters.end() );  
     optParam = *IterParam; // FIXME is this needed?
     dofcount+=optParam->dofs.size(); // update counter
    }
   // return parameter pointer corresponding to the global dof
   PetscFunctionReturn( IterParam );
  } 
 else
  { 
   // parameter not found
   PetscFunctionReturn( Parameters.end() );
  } 
}
/* -------------------------------------------------------------------- 
   plot variances from hessian
*/
#undef __FUNCT__
#define __FUNCT__ "qoiBaseClass::plotElemData"
void qoiBaseClass::plotElemData(OStringStream &file_name, 
                                libMesh::MeshBase &mesh,Vec elemData )
{
  PetscFunctionBegin;

  // This function must be run on all processors at once
  parallel_only();

  std::cout << "plotElemData not ready " << std::endl;
  libmesh_error();
  ////const unsigned int dim = mesh.mesh_dimension();
  //const unsigned int ne  = mesh.n_elem();

  ////get number of field variables and name of each variable
  //unsigned int nv  = 0;
  //std::vector < std::string > names;
  //for(IterParam=Parameters.begin(); IterParam!=Parameters.end(); IterParam++)
  //  {
  //   optimizationParameter* optParam = *IterParam;
  //   names.push_back(optParam->name()); 
  //   nv++; 
  //  }

  //// We'd better have a contiguous node numbering
  //libmesh_assert (ne == mesh.max_elem_id());

  //// allocate storage to hold
  //// (number_of_nodes)*(number_of_variables) entries.
  //std::vector<Number> soln(ne*nv);

  //// TODO must get entire vector locally to plot
  //VecScatter allToAllctx;
  //Vec fullElemData;
  //PetscErrorCode info;
  //info = VecScatterCreateToAll(elemData,&allToAllctx,&fullElemData);
  //CHKERRV(info);
  //
  //// scatter as many times as you need 
  //info = VecScatterBegin(allToAllctx,elemData,fullElemData,INSERT_VALUES,SCATTER_FORWARD);
  //info = VecScatterEnd(  allToAllctx,elemData,fullElemData,INSERT_VALUES,SCATTER_FORWARD);
  //CHKERRV(info);
  //CHKERRV(info);
  //
  //// get pointer to elemData  data
  //PetscScalar    *varData;
  //info = VecGetArray(fullElemData,&varData);CHKERRV(info);

  //// Zero out the soln vector
  //std::fill (soln.begin(),       soln.end(),       libMesh::zero);

  //// For each system in this EquationSystems object,
  //// update the global solution and if we are on processor 0,
  //// loop over the elements and build the nodal solution
  //// from the element solution.  Then insert this nodal solution
  //// into the vector passed to build_solution_vector.

  //libMesh::MeshBase::const_element_iterator       it  = mesh.active_local_elements_begin();
  //const libMesh::MeshBase::const_element_iterator end = mesh.active_local_elements_end(); 
  //for ( ; it != end; ++it)
  //  {// element iterator
  //   const Elem* elem = *it;

  //   // counter used for storage
  //   for(IterParam=Parameters.begin(); IterParam!=Parameters.end(); IterParam++)
  //     {
  //      //optimizationParameter* optParam = *IterParam;
  //      const int paramID = distance(Parameters.begin(),IterParam);
  //      const int globalDof = elem->_param_opt_map.at(paramID);
  //      // plot all params
  //      soln[ nv*elem->id() + paramID ] = varData[ globalDof ];
  //     }
  //  }	 

  //// restore pointer
  //info = VecRestoreArray(elemData ,&varData);CHKERRV(info);

  //// destroy scatter context and local vector when no longer needed
  //info = VecScatterDestroy(allToAllctx); CHKERRV(info);
  //info = VecDestroy(fullElemData); CHKERRV(info);

  //CHKMEMA; // check for memory corruption use -malloc_debug to enable
  //// Now each processor has computed contriburions to the
  //// soln vector.  Gather them all up.
  //Parallel::sum(soln);

  //// initialize file output of main data structures
  //ExodusII_IO output_var(mesh);
  //// Write out the initialize vis of the field
  //output_var.write_elem_data(file_name.str(),soln,names);

  PetscFunctionReturnVoid();
}
/* -------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "qoiBaseClass::computeSensitivity"
void qoiBaseClass::computeSensitivity( AppSolve *user, Vec )
{
 PetscErrorCode info; /* used to check for functions returning nonzeros */
 PetscFunctionBegin;

 // libmesh eqn systems
 EquationSystems &es = user->get_equation_systems();  
 
 // Get a reference to the NonlinearImplicitSystem we have solved
 TransientFEMSystem& state_system = 
      es.get_system<TransientFEMSystem>("StateSystem");

 // Get a reference to the LinearImplicitSystem for sensitivity problem
 TransientFEMSystem & sensitivity_system =
      es.get_system<TransientFEMSystem>("SensitivitySystem");

 // compute sensitivity 
 // due to memory limitations need to keep the sensitivity solution in
 // memory and must recompute jacobian at each time step 
 /** 
  * one sensitivity solve for each parameter
  *  FIXME: The same matrix may be used for each sensivity solve but
  *  due to memory limitations need to keep the sensitivity solution in
  *  memory and must recompute jacobian at each time step. 
  *  TODO: test if it is faster to read/write senstivity solution from disk
  *  and not recompute jacobian
  */

 // ensure the same memory is used for the sensitivity matrix and state system
 // jacobian
 libmesh_assert( state_system.matrix == sensitivity_system.matrix);

 // assemble
 state_system.assembly(false,true);

 /* the input vector determines the sensitivity load */
 state_system.sensitivityLoad( *sensitivity_system.rhs);

 // debugging
 PetscVector<Number>* rhsVec = libmesh_cast_ptr<PetscVector<Number>*>( 
                                            &(*sensitivity_system.rhs) );
 //std::cout << "writeSensitivityLoad" << globalDof 
 //          << "istep" << AppSolve::ISTEP << std::endl;
 //info = VecView(rhsVec->vec(),0);CHKERRV(info);

 // Perform Nonlinear solve with PetscNonlinearSolver.
 PetscLogEventBegin(AppSolve::logevents[5],0,0,0,0); // libMesh Solve
 PetscDiffSolver* sensitivitySolver = 
                     libmesh_cast_ptr<PetscDiffSolver*>(
                      & (*sensitivity_system.time_solver->diff_solver()) );
 // use zero as default initial guess instead of previous sensitivity
 // TODO - is this the best way? previous sensitivity as the initial guess
 // caused  senstitivity solve to diverge do to KSP_DIVERGED_DTOL 
 // info = KSPSetInitialGuessNonzero (sensitivityLinearSolver->ksp(), 
 //                                   PETSC_FALSE); CHKERRV(info);
 sensitivity_system.solve(); //solve

 // debugging
 PetscVector<Number>* solnVec = libmesh_cast_ptr<PetscVector<Number>*>( 
                                       &(*sensitivity_system.solution) );
 //std::cout << "writeSensitivitySoln" << globalDof 
 //          << "istep" << AppSolve::ISTEP << std::endl;
 //info = VecView(solnVec->vec(),0);CHKERRV(info);

 //check convergence reason
 KSP  snesksp;
 SNESGetKSP(sensitivitySolver->snes() , &snesksp);
 KSPConvergedReason kspReason;
 KSPGetConvergedReason(snesksp , &kspReason);
 if(kspReason<0)
   {// echo debuggin info and write out diverging system to study
    if( !libMesh::processor_id() ) std::cerr
                                <<"Sensitivity System Diverged "<< kspReason
                                <<" at time step "<< AppSolve::ISTEP
                                << std::endl << std::flush ; 

    PetscViewer viewHandle; // file handle for petsc viewer
    // write initial Guess. the prev solution  should have been used 
    // as the initial guess and should be unchanged in the diverged system
    OStringStream sensitivitySystemFileName;
    sensitivitySystemFileName <<"files/sensitivityIG" 
                              <<"istep"    << AppSolve::ISTEP
                              <<"Diverged" << kspReason
                              << AppSolve::profileID << ".dat";
    info = PetscViewerBinaryOpen(PETSC_COMM_WORLD,
                                 sensitivitySystemFileName.str().c_str(),
                                 FILE_MODE_WRITE,&viewHandle);CHKERRV(info);
    info = VecView(solnVec->vec(),viewHandle);CHKERRV(info);
    info = PetscViewerDestroy(viewHandle);CHKERRV(info);

    // write rhs 
    sensitivitySystemFileName.str(""); // reset before reuse
    sensitivitySystemFileName <<"files/sensitivityRHS" 
                              <<"istep"    << AppSolve::ISTEP
                              <<"Diverged" << kspReason
                              << AppSolve::profileID << ".dat";
    info = PetscViewerBinaryOpen(PETSC_COMM_WORLD,
                                 sensitivitySystemFileName.str().c_str(),
                                 FILE_MODE_WRITE,&viewHandle);CHKERRV(info);
    info = VecView(rhsVec->vec(),viewHandle);CHKERRV(info);
    info = PetscViewerDestroy(viewHandle);CHKERRV(info);

    // write matrix 
    PetscMatrix<Number>* sysMatrix = libmesh_cast_ptr<PetscMatrix<Number>*>( 
                                             &(*sensitivity_system.matrix) );
    sensitivitySystemFileName.str(""); // reset before reuse
    sensitivitySystemFileName <<"files/sensitivityMatrix" 
                              <<"istep"    << AppSolve::ISTEP
                              <<"Diverged" << kspReason
                              << AppSolve::profileID << ".dat";
    info = PetscViewerBinaryOpen(PETSC_COMM_WORLD,
                                 sensitivitySystemFileName.str().c_str(),
                                 FILE_MODE_WRITE,&viewHandle);CHKERRV(info);
    info = MatView(sysMatrix->mat(),viewHandle);CHKERRV(info);
    info = PetscViewerDestroy(viewHandle);CHKERRV(info);
  
    // try zero solution as initial guess
    sensitivity_system.solution->zero() ;
    sensitivity_system.solve(); //solve

    //check convergence reason again
    KSPGetConvergedReason(snesksp , &kspReason);
    if(kspReason<0)
      {// echo debuggin info and write out diverging system to study
       if( !libMesh::processor_id() ) 
             std::cerr <<"Sensitivity System zero IG Diverged "<< kspReason
                       <<" at time step "<< AppSolve::ISTEP
                       << std::endl << std::flush ; 
       // no point in continuing
       //libmesh_error(); 
 
       // ensure sensitivity is zero for write out
       sensitivity_system.solution->zero() ;
      } 
   }

 PetscLogEventEnd(  AppSolve::logevents[5],0,0,0,0); // libMesh Solve

 // copy parallel data structures to local data structures
 sensitivity_system.solution->localize(
                               *sensitivity_system.current_local_solution);

 PetscLogEventBegin(AppSolve::logevents[16],0,0,0,0); // data write

 // store the sensitivity if necessary
 //this->storeSensitivity(user, globalDof, 
 //                             *sensitivity_system.current_local_solution );

 PetscLogEventEnd(  AppSolve::logevents[16],0,0,0,0); // data write

 PetscFunctionReturnVoid();
}

// subroutine to extract image data to a NumericVector
#undef __FUNCT__
#define __FUNCT__ "qoiBaseClass::ExtractImageData"
PetscErrorCode qoiBaseClass::
ExtractImageData( InputImageType::Pointer tempImage, // image to extract from
                  TransientFEMSystem & local_system ) // system to put soln
{
  PetscFunctionBegin;

  // set pointer
  this->Images->interpolator->SetInputImage( tempImage );

  // put MRTI thermal image into FEM data structures
  AppSolve *user = AppSolve::getInstance();
  EquationSystems    &EqnSystem = user->get_equation_systems();
  local_system.project_solution(GetITKImageData,NULL,
                                EqnSystem.parameters);

  PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "qoiBaseClass::PlotStoredData"
void qoiBaseClass::PlotStoredData( AppSolve *user)
{ // plot what is already stored
  TransientFEMSystem    *pdeSolver    = user->pdeSolver();
  libMesh::MeshBase                &mesh = user->get_mesh();
  EquationSystems    &EqnSystem = user->get_equation_systems();
  PetscFunctionBegin;

  std::cout << "PlotStoredData not ready " << std::endl;
  libmesh_error();
  //PetscLogEventBegin(AppSolve::logevents[16],0,0,0,0); // data write

  //// set the file output name
  //OStringStream file_name;

  //// initialize data buffers
  //std::vector<Number>      soln;
  //
  //PetscInt TaoIteration = 0 ;
  //if(this->tao) TaoIteration = this->tao->iter ;

  //// Write out the element fields at this iteration for visualization
  //file_name << "femvis/femqoi_"<<user->GroupID<<"nopt_"<<user->IDOPT;
  //file_name << "mat" << AppSolve::profileID  << ".e";
  //// FIXME: iter is being updated in TaoSolve_BLMVM before the write
  //// add one to avoid overwrite of initial time
  //Real FakeTime = TaoIteration + 1;
  ////pdeSolver->writeElementData(file_name,mesh,soln,FakeTime);

  //// initialize file output of main data structures
  //ExodusII_IO output_vis(mesh);

  //// set main solution filename
  //file_name.str(""); // reset before reuse
  //file_name << "femvis/femqoi_"<<user->GroupID<<"nopt_"<<user->IDOPT;
  //// OSSRealzeroright(file_name,4,0, idplot);
  //// FIXME: iter is being updated in TaoSolve_BLMVM before the write
  //file_name << "iter_" << TaoIteration 
  //          << AppSolve::profileID  << ".e";

  //// get names form the system
  //std::vector<std::string> names;
  //for (unsigned int i_sys = 0 ; i_sys < 
  //                              EqnSystem.n_systems() ; i_sys++)
  //  {  
  //    System& system  = EqnSystem.get_system(i_sys) ; 
  //    const unsigned int nv_sys = system.n_vars();
  //    for (unsigned int var=0; var<nv_sys; var++)
  //     names.push_back( system.variable_name(var) );
  //  }  

  ////plot the stored history...
  //for(AppSolve::ISTEP  = this->Nsteplo(); 
  //    AppSolve::ISTEP <= this->Nstephi(); 
  //    AppSolve::ISTEP++)
  // {
  //  // write every nth interval
  //  if( user->modWrite(AppSolve::ISTEP) ) //write every nth step
  //   {
  //    //PetscInt idplot = AppSolve::ISTEP / user->WriteInterval() ; 
  //    // Main_Visualization(user->ControlComm,user, &user->GroupID,
  //    //                    &user->IDOPT, &idplot, &AppSolve::ISTEP, 
  //    //                    GenInfo,QOIInfo,user->x, user->p);
  //    
  //    // set time
  //    Real Time = user->getTime(AppSolve::ISTEP);
  //    
  //    // A pretty update message
  //    std::cout << "\n\n*** Writing time step "  << AppSolve::ISTEP 
  //              << ", time = " << Time << " ***" << std::endl;

  //    // restore any stored data to be plotted in current_local_solution
  //    for (unsigned int i_sys = 0 ; i_sys < EqnSystem.n_systems() ; i_sys++)
  //     {
  //      System & local_system = EqnSystem.get_system(i_sys);
  //      TransientFEMSystem *transient_local_system = 
  //                         dynamic_cast<TransientFEMSystem*>(&local_system);
  //      if( transient_local_system && 
  //          AppSolve::ISTEP< (int)transient_local_system->vector_solution.size() )
  //        *transient_local_system->current_local_solution=
  //            *transient_local_system->vector_solution.at(AppSolve::ISTEP);
  //     }

  //    // build data buffer from current_local_solution and write
  //    build_equation_system_solution_vector(EqnSystem, soln);

  //    // Write out every ideal timestep to file.
  //    output_vis.write_nodal_data(file_name.str(),soln,names);
  //    // write timestep info to file 
  //    output_vis.write_timestep(Time);
  //    // increment timestep (1-based) 
  //    output_vis.increment_timestep();

  //    /* run code verification suite */ 
  //    if( user->get_num_exact() ) pdeSolver->verifySolution(EqnSystem);

  //   } // end if( AppSolve::modWrite(AppSolve::ISTEP) )
  // } // end loop over time steps

  //PetscLogEventEnd(  AppSolve::logevents[16],0,0,0,0); // data write


  PetscFunctionReturnVoid();
} // end plot

#undef __FUNCT__
#define __FUNCT__ "qoiBaseClass::PlotInterleaveCompute"
void qoiBaseClass::PlotInterleaveCompute( AppSolve *user)
{ // plot what is already stored
  libMesh::MeshBase                &mesh = user->get_mesh();
  EquationSystems    &EqnSystem = user->get_equation_systems();
  PetscErrorCode info; /* used to check for functions returning nonzeros */
  PetscFunctionBegin;

  std::cout << "PlotInterleaveCompute not ready " << std::endl;
  libmesh_error();
  //// initialize data buffers
  //std::vector<Number>      soln;
  //
  //// initialize file output of main data structures
  //ExodusII_IO output_vis(mesh);

  //// set the file output name
  //OStringStream file_name;
  //file_name << "femvis/fempostqoi_"<<user->GroupID 
  //          << "idopt_"<<user->IDOPT
  //          << AppSolve::profileID  << ".e";

  //// get names form the system
  //std::vector<std::string> names;
  //for (unsigned int i_sys = 0 ; i_sys < 
  //                              EqnSystem.n_systems() ; i_sys++)
  //  {  
  //    System& system  = EqnSystem.get_system(i_sys) ; 
  //    const unsigned int nv_sys = system.n_vars();
  //    for (unsigned int var=0; var<nv_sys; var++)
  //     names.push_back( system.variable_name(var) );
  //  }  
  //// compute all the way to the end of MRTI_ntime
  //for(AppSolve::ISTEP  = this->Nstephi()+1 ; 
  //    AppSolve::ISTEP <= user->get_num_MRTI() * user->IstepsPerIdeal();
  //    AppSolve::ISTEP++)
  // {
  //   // set time
  //   Real Time = user->getTime(AppSolve::ISTEP);

  //   // update the solution vector from the previous time step
  //   TransientFEMSystem & state_system =
  //      EqnSystem.get_system<TransientFEMSystem>("StateSystem");
  //   *state_system.old_vector_solution[0]= *state_system.current_local_solution;

  //   /* !!!NOTE!!! Do not need to store the soln in data structures 
  //                 after the solve b/c this was already done by the 
  //                 last function evaluation */
  //   // Perform NonlinearImplicitSystem::solve () 
  //   //     with PetscNonlinearSolver. 
  //   //      compute_residual  
  //   //      compute_jacobian 
  //   // This will put a local copy of solution into 
  //   // current_local_solution.
  //   PetscLogStagePush(AppSolve::logstages[2]); // fnc evaluation
  //   state_system.solve();
  //   PetscLogStagePop();     // fnc evaluation
  // 
  //   /* compute the data structures for the QOI
  //      PostSolveQoiData is a function pointer to  
  //                 CurrentIdealArr     : arrhenius base optimize
  //                                                               */

  //   if(this->PostSolveQoiData)
  //     {
  //      info=this->PostSolveQoiData(user);CHKERRV(info);
  //     }
  //   // write every nth interval
  //   if( user->modWrite(AppSolve::ISTEP) ) //write every nth step
  //    {
  //     // A pretty update message
  //     std::cout << "\n\n*** Writing time step " << AppSolve::ISTEP 
  //               << ", time = " << Time
  //               << " ***" << std::endl;
  //     // build data buffer from current_local_solution and write
  //     build_equation_system_solution_vector(EqnSystem, soln);
  //     output_vis.write_nodal_data(file_name.str(),soln,names);
  //     // write timestep info to file 
  //     output_vis.write_timestep(Time);
  //     // increment timestep (1-based) 
  //     output_vis.increment_timestep();
  //    }
  // } // end loop over time steps
  PetscFunctionReturnVoid();
} // end plot
