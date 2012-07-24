/*
       This file contains any fortran calls from C++ for all classes
*/
// libMesh 
#include "libmesh.h"
#include "equation_systems.h"
#include "exact_solution.h"
#include "fe.h"
#include "fe_interface.h"
#include "quadrature_gauss.h"
#include "dof_map.h"
#include "sparse_matrix.h"
#include "elem.h"
#include "string_to_enum.h"
#include "boundary_info.h"
#include "getpot.h" // file parser
#include "parallel.h"
#include "dense_submatrix.h"
// The nonlinear solver we will be using
#include "steady_solver.h"
#include "petsc_diff_solver.h"
#include "petsc_vector.h"
#include "petsc_matrix.h"

// tao interface
#include "src/tao_impl.h" 

// local
#include "applicationContext.h"
#include "pennesVerification.h"
#include "pennesInverseSystem.h"
#include "quantityOfInterest.h"
#include "dddas.h"
#include "variable_map.h"
#include "parser_defaults.h"
#include "mesh.h"
#include "thermal_therapy_system.txx"
#include "pennesSystem.txx"
#include "pennesInverseSystem.txx"
#include "quantityOfInterest.txx"

#include "fortrandftranslate.h" // header to call fortran routines
static PetscInt nzero  = 0 ,
                negone =-1 , 
                none   = 1 ;
/* --------------External Fortran Routines---------- */
extern FORTRAN_FUNCTION 
{
#include "global_params.h" // global data
#include "pennes_model.h" // interface to pennes model
  PetscTruth FORTRAN_NAME(islittleendian)()
  {
  	assert(sizeof(short int)==2);
  
          short int word = 0x0001;
          char *byte = (char *) &word;
          if(byte[0]) return PETSC_TRUE;
          else        return PETSC_FALSE;
  
  }
  /* read visualase variables from the ini file 
      to dynamically change it during runtime*/
  void FORTRAN_NAME(dynamic_laser_control)( PetscTruth  *Override, 
                                                   PetscScalar *Override_Power )
  {
    GetPot controlfile(dfltINIFile);
       
    //median filter neighborhood size
    try
    {  // note if file not open will return default
      *Override_Power=(PetscScalar)controlfile("visualase/override_power",0.0);
      if(controlfile("visualase/override",false)){
         *Override = PETSC_TRUE ; 
      } else { 
         *Override = PETSC_FALSE; 
      }
    }
    catch(...) //under any expection return defaults
    {
      std::cout << "dynamic_laser_control: input error, using defaults\n";
      *Override_Power = 0.0;
      *Override       = PETSC_FALSE; 
    }
    return;
  }
  /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  PetscTruth FORTRAN_NAME(stringcomp)(char *value,char *testvalue){ 
    if( !strcasecmp( value ,testvalue) ) return PETSC_TRUE ;
    else return PETSC_FALSE ;
  }
  /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  void FORTRAN_NAME(printpetscscalar)( char *message , PetscScalar *value ) 
  {
    std::cout << message << *value << std::endl ; 
    return;
  }
  /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  void FORTRAN_NAME(printpetscint)( char *message , PetscInt *value ) 
  {
    std::cout << message << *value << std::endl ; 
    return;
  }
  /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  // print out address of fortran arrays (doubles)
  void FORTRAN_NAME(printdble)(int *tid,char *label,double *a){ 
              printf("thread %d %s %p\n",*tid,label,(void *)a);
              return; 
  }
  /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  // print out address of fortran arrays (int)
  void FORTRAN_NAME(printint)(int *tid,char *label,int *a){ 
              printf("thread %d %s %p\n",*tid,label,(void *)a);
              return; 
  }

}
/* -------------------------------------------------------------------- 
   Verifications problems using fortran modules
   -------------------------------------------------------------------- */ 
void PennesVerification::printSelf(std::ostream& os)
{
  PennesStandardDiffusionApproximation::printSelf(os);
  //os << "PennesVerification: gs=" << gs <<" f= " << f << std::endl;
}
PennesVerification::PennesVerification(GetPot &controlfile,libMesh::MeshBase &mesh):
                        PennesInverseModel(controlfile,mesh)
{
  //FORTRAN_NAME(meshcheck)();
  //PetscPrintf(GenInfo.DDDAS_COMM_WORLD,"@@@@@@@@@@@@@initial mesh generated\n" );
 // scaling 
 m_TimeDerivativeScalingFactor = 1.0;

 // Get a constant reference to the mesh object.
 AppSolve *user = AppSolve::getInstance();
 //if( user->get_num_exact() == 1 )
 // {
 //  m_theta = 1.0; 
 // }

 if(user->get_num_exact() )
  {
   // default initial condition is domain wise constant
   InitValues.at(0) = &PDEModelBaseClass::getExactTemperature;
  }

  PetscErrorCode info =  PetscInitializeFortran(); 
  PetscTruth dftrue  = PETSC_TRUE;
  PetscTruth dffalse = PETSC_FALSE;
  FORTRAN_NAME(initialize_visualase)(&dftrue,&dftrue,&dffalse,&nzero);
  PetscInt numField = user->get_num_field_elem(); 
  FORTRAN_NAME(setup_pennes_model)(&numField); 
  PetscScalar spacing[3] = {1.0,1.0,1.0};
  FORTRAN_NAME(read_control_fortran)(&negone,spacing);
}
/* ---------------------------------------------------------- */
#undef __FUNCT__ 
#define __FUNCT__ "PennesVerification::WriteControlFile"
PetscErrorCode 
PennesVerification::WriteControlFile(PetscInt IDFileOut, PetscInt GroupID)
{
  PetscScalar beta; //scratch variable
  PetscFunctionBegin; 

  std::ofstream Inifile;
  std::ostringstream fileID ; // filename
  fileID << "files/control"<< IDFileOut <<".ini";
  Inifile.open(fileID.str().c_str(), std::ios::out);

  // perfusion
  Inifile <<"[perfusion]" <<std::endl;
  FORTRAN_NAME(getparam_w_n)(&beta,&none);
  Inifile << "w_n=" << beta << std::endl;
  FORTRAN_NAME(getparam_w_i)(&beta,&none);
  Inifile << "w_i=" << beta << std::endl;
  FORTRAN_NAME(getparam_w_d)(&beta,&none);
  Inifile << "w_d=" << beta << std::endl;
  FORTRAN_NAME(getparam_w_2)(&beta,&none);
  Inifile << "w_2=" << beta << std::endl;
  FORTRAN_NAME(getparam_w_ni)(&beta,&none);
  Inifile << "w_nI=" << beta << std::endl;
  FORTRAN_NAME(getparam_w_id)(&beta,&none);
  Inifile << "w_id=" << beta << std::endl;
  // perfusion lb
  FORTRAN_NAME(getparam_w_0_lb)(&beta,&none);
  Inifile << "w_0_lb=" << beta << std::endl;
  FORTRAN_NAME(getparam_w_n_lb)(&beta,&none);
  Inifile << "w_n_lb=" << beta << std::endl;
  FORTRAN_NAME(getparam_w_i_lb)(&beta,&none);
  Inifile << "w_i_lb=" << beta << std::endl;
  FORTRAN_NAME(getparam_w_d_lb)(&beta,&none);
  Inifile << "w_d_lb=" << beta << std::endl;
  FORTRAN_NAME(getparam_w_2_lb)(&beta,&none);
  Inifile << "w_2_lb=" << beta << std::endl;
  FORTRAN_NAME(getparam_w_ni_lb)(&beta,&none);
  Inifile << "w_ni_lb=" << beta << std::endl;
  FORTRAN_NAME(getparam_w_id_lb)(&beta,&none);
  Inifile << "w_id_lb=" << beta << std::endl;
  // perfusion ub
  FORTRAN_NAME(getparam_w_0_ub)(&beta,&none);
  Inifile << "w_0_ub=" << beta << std::endl;
  FORTRAN_NAME(getparam_w_n_ub)(&beta,&none);
  Inifile << "w_n_ub=" << beta << std::endl;
  FORTRAN_NAME(getparam_w_i_ub)(&beta,&none);
  Inifile << "w_i_ub=" << beta << std::endl;
  FORTRAN_NAME(getparam_w_d_ub)(&beta,&none);
  Inifile << "w_d_ub=" << beta << std::endl;
  FORTRAN_NAME(getparam_w_2_ub)(&beta,&none);
  Inifile << "w_2_ub=" << beta << std::endl;
  FORTRAN_NAME(getparam_w_ni_ub)(&beta,&none);
  Inifile << "w_ni_ub=" << beta << std::endl;
  FORTRAN_NAME(getparam_w_id_ub)(&beta,&none);
  Inifile << "w_id_ub=" << beta << std::endl;

  // thermal conductivity
  Inifile <<"[thermal_conductivity]" <<std::endl;
  FORTRAN_NAME(getparam_k_1)(&beta,&none);
  Inifile << "k_1=" << beta << std::endl;
  FORTRAN_NAME(getparam_k_2)(&beta,&none);
  Inifile << "k_2=" << beta << std::endl;
  FORTRAN_NAME(getparam_k_3)(&beta,&none);
  Inifile << "k_3=" << beta << std::endl;
  // thermal conductivity lb
  FORTRAN_NAME(getparam_k_0_lb)(&beta,&none);
  Inifile << "k_0_lb=" << beta << std::endl;
  FORTRAN_NAME(getparam_k_1_lb)(&beta,&none);
  Inifile << "k_1_lb=" << beta << std::endl;
  FORTRAN_NAME(getparam_k_2_lb)(&beta,&none);
  Inifile << "k_2_lb=" << beta << std::endl;
  FORTRAN_NAME(getparam_k_3_lb)(&beta,&none);
  Inifile << "k_3_lb=" << beta << std::endl;
  // thermal conductivity ub
  FORTRAN_NAME(getparam_k_0_ub)(&beta,&none);
  Inifile << "k_0_ub=" << beta << std::endl;
  FORTRAN_NAME(getparam_k_1_ub)(&beta,&none);
  Inifile << "k_1_ub=" << beta << std::endl;
  FORTRAN_NAME(getparam_k_2_ub)(&beta,&none);
  Inifile << "k_2_ub=" << beta << std::endl;
  FORTRAN_NAME(getparam_k_3_ub)(&beta,&none);
  Inifile << "k_3_ub=" << beta << std::endl;

  // probe
  Inifile <<"[probe]" <<std::endl;
  FORTRAN_NAME(getparam_x_0)(&beta,&none);
  Inifile << "x_0="  << beta << std::endl;
  FORTRAN_NAME(getparam_y_0)(&beta,&none);
  Inifile << "y_0="  << beta << std::endl;
  FORTRAN_NAME(getparam_z_0)(&beta,&none);
  Inifile << "z_0="  << beta << std::endl;
  // probe lb
  FORTRAN_NAME(getparam_x_0_lb)(&beta,&none);
  Inifile << "x_0_lb="  << beta << std::endl;
  FORTRAN_NAME(getparam_y_0_lb)(&beta,&none);
  Inifile << "y_0_lb="  << beta << std::endl;
  FORTRAN_NAME(getparam_z_0_lb)(&beta,&none);
  Inifile << "z_0_lb="  << beta << std::endl;
  // probe ub
  FORTRAN_NAME(getparam_x_0_ub)(&beta,&none);
  Inifile << "x_0_ub="  << beta << std::endl;
  FORTRAN_NAME(getparam_y_0_ub)(&beta,&none);
  Inifile << "y_0_ub="  << beta << std::endl;
  FORTRAN_NAME(getparam_z_0_ub)(&beta,&none);
  Inifile << "z_0_ub="  << beta << std::endl;

  // optical 
  Inifile <<"[optical]" <<std::endl;
  FORTRAN_NAME(getparam_mu_a)(&beta,&none);
  Inifile << "mu_a=" << beta << std::endl;
  FORTRAN_NAME(getparam_mu_s)(&beta,&none);
  Inifile << "mu_s=" << beta << std::endl;
  FORTRAN_NAME(getparam_mu_a_lb)(&beta,&none);
  Inifile << "mu_a_lb=" << beta << std::endl;
  FORTRAN_NAME(getparam_mu_s_lb)(&beta,&none);
  Inifile << "mu_s_lb=" << beta << std::endl;
  FORTRAN_NAME(getparam_mu_a_ub)(&beta,&none);
  Inifile << "mu_a_ub=" << beta << std::endl;
  FORTRAN_NAME(getparam_mu_s_ub)(&beta,&none);
  Inifile << "mu_s_ub=" << beta << std::endl;

  // write out the restart ID to restart from this file
  Inifile <<"[compexec]" << std::endl;
  Inifile <<"restartID="    << IDFileOut << std::endl;
  // tag the file with the group number that wrote it
  Inifile <<"id_qoi_group=" << GroupID   << std::endl;

  // write solution to file for next solve
  Inifile.close();

  /* write_power will write to "files/power%d.dat" % IDFileOut */
  FORTRAN_NAME(write_power)(&IDFileOut  );

//
//  /* gather field info */
//  vector<MdleInfo> *MdleData = &GenInfo.MdleData; 
//  /* allocate buffer for inter groups field transmission */
//  PetscInt iii,nelems = GenInfo.MdleData.size() ;
//  FORTRAN_NAME(alloc_hpelemfield)( &nelems );
//  /* pack the  buffer for inter groups field transmission */
//  for(mdleID=MdleData->begin() ; mdleID!=MdleData->end() ; mdleID++){
//     iii = distance(MdleData->begin(),mdleID );
//     FORTRAN_NAME(pack_hpelemfield)(&iii,&mdleID->idelem,&mdleID->mdle);
//  }
//  /* setup for varying perfusivity and conductivity fields */
//  PetscMPIInt fortran_comm = PetscFromPointerComm(user->ControlComm);
//  FORTRAN_NAME(gatherdealloc_hpelemfield)(&fortran_comm,
//                                          &QOIInfo.optimize_W_0,
//                                          &QOIInfo.optimize_K_0,
//                                          &QOIInfo.ntask,
//                                          &QOIInfo.recvcnts_elm[0],
//                                          &QOIInfo.displs_elm[0] );
  PetscFunctionReturn(0);
}
/* ---------------------------------------------------------- */
PetscScalar PennesVerification::getMass()
{
 PetscFunctionBegin; 
 const PetscScalar pennes_mass = FORTRAN_NAME(get_pennes_mass)() ;
 PetscFunctionReturn(pennes_mass); 
}
/* ---------------------------------------------------------- */
void PennesVerification::XferData()
{
 PetscFunctionBegin; 
   FORTRAN_NAME(putparam_w_0)(  &w_0[  0] , &nzero );
   FORTRAN_NAME(putparam_k_0)(  &k_0[  0] , &nzero );
   FORTRAN_NAME(putparam_mu_a)(&mu_a[  0] , &nzero );
   FORTRAN_NAME(putparam_mu_s)(&mu_s[  0] , &nzero );
   FORTRAN_NAME(putparam_x_0)(  &x_0[  0] , &nzero );
   FORTRAN_NAME(putparam_y_0)(  &y_0[  0] , &nzero );
   FORTRAN_NAME(putparam_z_0)(  &z_0[  0] , &nzero );
   FORTRAN_NAME(putparam_pow)(  &Power[0] , &nzero );
   //FORTRAN_NAME(putparam_w_n)(  );
   //FORTRAN_NAME(putparam_w_i)(  );
   //FORTRAN_NAME(putparam_w_d)(  );
   //FORTRAN_NAME(putparam_w_2)(  );
   //FORTRAN_NAME(putparam_w_ni)( );
   //FORTRAN_NAME(putparam_w_id)( );
   //FORTRAN_NAME(putparam_k_1)(  );
   //FORTRAN_NAME(putparam_k_2)(  );
   //FORTRAN_NAME(putparam_k_3)(  );
 PetscFunctionReturnVoid(); 
}

Number PennesVerification::getExactTemperature(unsigned int , 
                                         const Point& p,
                                         const Parameters& )
{
 AppSolve *user = AppSolve::getInstance();
 PetscScalar time = user->getTime(AppSolve::ISTEP);
 PetscScalar xpoint[3] = { p(0), p(1), p(2) };
 PetscScalar texact = 0.0;
 PetscScalar gradexact[3] = { 0.0, 0.0, 0.0 };

 FORTRAN_NAME(forwardexact)(xpoint,&time,&texact,gradexact);

 // :( special case for verification problem one
 if( !AppSolve::ISTEP && user->get_num_exact()  == 1 )  texact = 0.0;

 return texact;
}
// check verification problems for pennes solver
#undef __FUNCT__
#define __FUNCT__ "PennesVerification::verifySolution"
void PennesVerification::verifySolution( EquationSystems &EqnSystem)
{
  PetscFunctionBegin;

  // Compute the error 
  // **note** if no functions attached then the error will return the norm

  // forward problem.
  libMesh::ExactSolution state_exact(EqnSystem);
  state_exact.attach_exact_value(
           pdeExactSolution< PennesInverseSystem < PennesVerification > >
                                );
  state_exact.compute_error("StateSystem","u0");

  // adjoint problem.
  libMesh::ExactSolution adjoint_exact(EqnSystem);
  adjoint_exact.attach_exact_value(pennes_adjoint_exact);
  //project ideal solution for plotting
  AppSolve *user = AppSolve::getInstance();
  if( user->get_num_exact() )
    {
      EqnSystem.get_system<ExplicitSystem>
      ("ExactAdjoint").project_solution(pennes_adjoint_exact,
                              NULL,EqnSystem.parameters);
    }
  adjoint_exact.compute_error("AdjointSystem","p0");
   
  // get forward problem error
  PetscScalar forwarderror = state_exact.l2_error("StateSystem", "u0");
  // get adjoint problem error
  PetscScalar adjointerror = adjoint_exact.l2_error("AdjointSystem", "p0");

  // forward problem verification
  std::cout << "Forward:NEXACTVERIFNUMBER=" << user->get_num_exact()  << std::endl;
  // tolerances stored for verification suite 
  PetscScalar errorTol;
  switch (user->get_num_exact() ) 
   {
    case 1:
      errorTol=8.0e-1;
      break;
    case 2:
      errorTol=2.5e-1;
      break;
    case 3:
      errorTol=2.0e0;
      break;
    case 4:
      errorTol=4.0e0;
      break;
    case 5:
      errorTol=1.0e-06;
      break;
    case 6:
      errorTol=0.35e0;
      break;
    case 7:
      errorTol=5.0e-2;
      break;
    case 8:
      errorTol=2.0e0   ;
      break;
    default:
      errorTol=1.0e0;
      forwarderror=errorTol;
      break;
   }
  if(AppSolve::ISTEP) AppSolve::indicateError(forwarderror,errorTol);

  // adjoint problem verification
  std::cout << "Adjoint:NEXACTVERIFNUMBER=" << user->get_num_exact()  << std::endl;
  // errorTolerances stored for verification suite 
  switch (user->get_num_exact() ) 
   {
    case 1:
      errorTol=2.0e-1;
      break;
    case 3:
      errorTol=2.0e0;
      break;
    case 4:
      errorTol=4.0e0;
      break;
    case 5:
      errorTol=1.0e-06;
      break;
    case 6:
      errorTol=0.35e0;
      break;
    case 7:
      errorTol=5.0e-2;
      break;
    case 8:
      errorTol=2.0e0   ;
      break;
    default:
      errorTol=1.0e0;
      adjointerror=errorTol;
      break;
   }
  if(AppSolve::ISTEP) AppSolve::indicateError(adjointerror,errorTol);

  PetscFunctionReturnVoid();
}


void PennesVerification::accumulateResidual(unsigned int idVar, 
                                            DiffContext &context,
                                            TransientFEMSystem &system)
{

 FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

 std::cout <<"accumulateResidual not ready"<< std::endl;
 libmesh_error();
// // set the field id in the spatially varying data structures
// const unsigned int field_id =  c.elem->_mat_data_field;
//
// // First we get some references to cell-specific data that
// // will be used to assemble the linear system.
//
// // Element Jacobian * quadrature weights for interior integration
// const std::vector<Real> &JxW = 
//   c.element_fe_var[idVar]->get_JxW();
//
// // The velocity shape functions at interior quadrature points.
// const std::vector<std::vector<Real> >& phi = 
//   c.element_fe_var[idVar]->get_phi();
// this->phi[idVar] = &phi;
//
// // The velocity shape function gradients at interior
// // quadrature points.
// const std::vector<std::vector<RealGradient> >& dphi =
//   c.element_fe_var[idVar]->get_dphi();
// this->dphi[idVar] = &dphi;
//
// // Physical location of the quadrature points
// const std::vector<Point>& qpoint = 
//   c.element_fe_var[idVar]->get_xyz();
// 
// // The number of local degrees of freedom in each variable
// const unsigned int n_u_dofs = c.dof_indices_var[idVar].size(); 
// this->n_u_dofs[idVar] = n_u_dofs ;
// const DofMap& dof_map = system.get_dof_map();
// // get the number of variable in the state system
// const unsigned int n_vars = system.n_vars();
// for( unsigned int i_var = 0 ; i_var < n_vars ; i_var++)
//   dof_map.dof_indices (c.elem, this->dof_indices_u[i_var],i_var);
//
// // The subvectors and submatrices we need to fill:
// DenseSubVector<Number>  &Fu = *c.elem_subresiduals[idVar];
//      
// PetscFunctionBegin; 
// AppSolve *user = AppSolve::getInstance();
//
// this->XferData();// ensure same data
// for (unsigned int qp=0; qp< c.element_qrule->n_points(); qp++)
//  {
//   // evaluate Solution
//   this->evaluateSolutionForState(qp,system);
//
//   // Solution buffer 
//   std::vector<PetscScalar>  solnbuffer; 
//
//   // Timestep info 
//   const PetscInt    &ISTEP = AppSolve::ISTEP;
//   const PetscScalar &FEM_DT= user->get_fem_dt();
//
//   // notice that the current time gradient and solution
//   // are computed with the "solution" vector 
//   // "vector_solution" has not been updated yet
//   const unsigned n_vars =  system.n_vars();
//   for( unsigned int i_var = 0 ; i_var < n_vars ; i_var++)
//    {
//      // order is important in routines called by function pointers
//      solnbuffer.push_back(  u__current[i_var]   ); // i_var*5+0
//      solnbuffer.push_back(  u_previous[i_var]   ); // i_var*5+1
//      solnbuffer.push_back(grad_u_mtheta[i_var](0)); // i_var*5+2
//      solnbuffer.push_back(grad_u_mtheta[i_var](1)); // i_var*5+3
//      solnbuffer.push_back(grad_u_mtheta[i_var](2)); // i_var*5+4
//    }                                               
//   // initiailize buffers                          
//   const PetscScalar xpoint[3] = {qpoint[qp](0), 
//                                  qpoint[qp](1), 
//                                  qpoint[qp](2)};
//   const PetscInt Nsize = solnbuffer.size(); 
//   Real Adiff  = 0.0, Creac  = 0.0, Source = 0.0;
//   /* loop over the temperature degrees of freedom for each variable
//      and Evaluate the elliptic coefficients at the gauss points 
//      formfcngauss[] is a function pointer to: 
//                        - getverifformfcn
//                        - nonlinpennesisolaser
//                        - nonlinpennesmonte   
//                        - pennesidealmonte 
//                        - nonlinpennesrf
//                        - sourceemrf 
//                        - penaltytemperature
//                        - penaltyvoltage */
//   for( unsigned int i_var = 0 ; i_var < n_vars ; i_var++)
//    {
//      formfcngauss[field_id][i_var]( xpoint,&ISTEP,&FEM_DT,
//                                                &solnbuffer[0],&Nsize,
//                                                &Adiff,&Creac,&Source);
//      for (unsigned int i=0; i < n_u_dofs; i++)
//        {
//          Fu(i) += JxW[qp]*(
//                   Adiff*( dphi[i][qp]*grad_u_mtheta[i_var] ) 
//                            +
//                   (Creac-Source)*phi[i][qp]
//            	      );
//        }
//    }
//  } // end of the quadrature point (qp-loop) over element interior


 PetscFunctionReturnVoid(); 

}
void PennesVerification::accumulateJacobian(unsigned int idVar, 
                                            DiffContext &context,
                                            TransientFEMSystem &system)
{
 PetscFunctionBegin; 
 AppSolve *user = AppSolve::getInstance();

 std::cout <<"accumulateJacobian not ready"<< std::endl;
 libmesh_error();
// // Timestep info 
// const PetscScalar &FEM_DT= user->get_fem_dt();
//
// FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
//
// // set the field id in the spatially varying data structures
// const unsigned int field_id =  c.elem->_mat_data_field;
//
// // First we get some references to cell-specific data that
// // will be used to assemble the linear system.
//
// // Element Jacobian * quadrature weights for interior integration
// const std::vector<Real> &JxW = 
//   c.element_fe_var[idVar]->get_JxW();
//
// // The velocity shape functions at interior quadrature points.
// const std::vector<std::vector<Real> >& phi = 
//   c.element_fe_var[idVar]->get_phi();
// this->phi[idVar] = &phi;
//
// // The velocity shape function gradients at interior
// // quadrature points.
// const std::vector<std::vector<RealGradient> >& dphi =
//   c.element_fe_var[idVar]->get_dphi();
// this->dphi[idVar] = &dphi;
//
// // Physical location of the quadrature points
// const std::vector<Point>& qpoint = 
//   c.element_fe_var[idVar]->get_xyz();
// 
// // The number of local degrees of freedom in each variable
// const unsigned int n_u_dofs = c.dof_indices_var[idVar].size(); 
// this->n_u_dofs[idVar] = n_u_dofs ;
// const DofMap& dof_map = system.get_dof_map();
// // get the number of variable in the state system
// const unsigned int n_vars = system.n_vars();
// for( unsigned int i_var = 0 ; i_var < n_vars ; i_var++)
//   dof_map.dof_indices (c.elem, this->dof_indices_u[i_var],i_var);
//
// // The subvectors and submatrices we need to fill:
// DenseSubMatrix<Number> &Kuu = *c.elem_subjacobians[idVar][idVar];
//
// for (unsigned int qp=0; qp< c.element_qrule->n_points(); qp++)
//  {
//
//    // evaluate Solution
//    this->evaluateSolutionForState(qp,system);
//
//    /* function pointer for nonlinear function jacobian for SNESsolve 
//       assuming diffusion, reaction, convection pde form main and coupled terms
//        JacPennesLaser      : jacobian of pennes with laser source term
//        JacPennesRF         : jacobian of coupled pennes/RF model 
//        JacPenalty          : jacobian of penalty method
//        VerifJacPennesLaser : jacobian for verification problem of pennes laser
//    */
//    jacgauss[field_id](FEM_DT,qp,n_vars,
//                       this->n_u_dofs,this->phi,this->dphi,JxW,
//                        u__current,u_previous,grad_u_mtheta,Kuu);
//
//  } // end of the quadrature point (qp-loop) over element interior

 PetscFunctionReturnVoid(); 
}
#undef __FUNCT__
#define __FUNCT__ "JacPennesLaser"
/* --------------jacobian of pennes with laser source term---------- */
void JacPennesLaser(const PetscScalar &FEM_DT,const unsigned int &qp, const unsigned int &n_vars,
                    std::vector< unsigned int>  &n_u_dofs,
                    std::vector< const std::vector<std::vector<Real> >* > &phi,
                    std::vector< const std::vector<std::vector<RealGradient> >* > &dphi,
                    const std::vector<Real>& JxW,
                    std::vector<Real>   &u__current,
                    std::vector<Real>   &u_previous,
                    std::vector<RealGradient> &grad_u_theta,
                    DenseSubMatrix<Number> &Kij)
{

 PetscFunctionBegin; 

 libmesh_assert(n_vars == 1); // make sure solving the right system
 Real Adiff  = 0.0, Bconv = 0.0, Creac  = 0.0; 
 FORTRAN_NAME(nonlinpennesjac)(&FEM_DT,&u__current[0],&u_previous[0],
                                               &Adiff,&Bconv,&Creac);
 // Matrix contributions 
 for (unsigned int i=0; i<n_u_dofs[0]; i++)
   for (unsigned int j=0; j<n_u_dofs[0]; j++)
     Kij(i,j) += JxW[qp]*(
            Adiff * ( (*dphi[0])[i][qp] * (*dphi[0])[j][qp] ) 
                                   +
            Bconv * ( (*dphi[0])[i][qp] * grad_u_theta[0] ) * 
                                           (*phi[0])[j][qp]
                                   +
            Creac *    (*phi[0])[i][qp] *  (*phi[0])[j][qp] 
      	                          );
 PetscFunctionReturnVoid(); 
}
/* --------------jacobian of pennes with laser source term---------- */
#undef __FUNCT__
#define __FUNCT__ "VerifJacPennesLaser"
void VerifJacPennesLaser(const PetscScalar &FEM_DT,const unsigned int &qp, const unsigned int &n_vars,
                    std::vector< unsigned int>  &n_u_dofs,
                    std::vector< const std::vector<std::vector<Real> >* > &phi,
                    std::vector< const std::vector<std::vector<RealGradient> >* > &dphi,
                    const std::vector<Real>& JxW,
                    std::vector<Real>   &u__current,
                    std::vector<Real>   &u_previous,
                    std::vector<RealGradient> &grad_u_theta,
                    DenseSubMatrix<Number> &Kij)
{

 PetscFunctionBegin;

 libmesh_assert(n_vars == 1); // make sure solving the right system
 Real Adiff  = 0.0, Bconv = 0.0, Creac  = 0.0; 
 FORTRAN_NAME(getverifjac)(&FEM_DT,&u__current[0],&u_previous[0],
                                           &Adiff,&Bconv,&Creac);
 // Matrix contributions 
 for (unsigned int i=0; i<n_u_dofs[0]; i++)
   for (unsigned int j=0; j<n_u_dofs[0]; j++)
     Kij(i,j) += JxW[qp]*(
            Adiff * ( (*dphi[0])[i][qp] * (*dphi[0])[j][qp] ) 
                                   +
            Bconv * ( (*dphi[0])[i][qp] * grad_u_theta[0] ) * 
                                           (*phi[0])[j][qp]
                                   +
            Creac *    (*phi[0])[i][qp] *  (*phi[0])[j][qp] 
      	                          );
 PetscFunctionReturnVoid();
}
/* --jacobian of Penalty method for pennes with coupled RF source term------- */
#undef __FUNCT__
#define __FUNCT__ "JacPenalty"
void JacPenalty(const PetscScalar &,const unsigned int &qp, const unsigned int &n_vars,
                    std::vector< unsigned int>  &n_u_dofs,
                    std::vector< const std::vector<std::vector<Real> >* > &phi,
                    std::vector< const std::vector<std::vector<RealGradient> >* > &,
                    const std::vector<Real>& JxW,
                    std::vector<Real>   &,
                    std::vector<Real>   &,
                    std::vector<RealGradient> &,
                    std::vector< SubMatrixVector > &Kij)
{

 PetscFunctionBegin;

 Real Penalty  = FORTRAN_NAME(getpenalty)();
  
 // Loop over governing equations
 // diagonal submatrix contributions (no coupling terms)
 for (unsigned int i_var=0; i_var<n_vars; i_var++)
   for (unsigned int i=0; i<n_u_dofs[i_var]; i++)
    for (unsigned int j=0; j<n_u_dofs[i_var]; j++)
      Kij[i_var][i_var](i,j) += JxW[qp] * Penalty * 
                                   (*phi[i_var])[i][qp] * (*phi[i_var])[j][qp] ;
 PetscFunctionReturnVoid();
}
/* -----adjoint matrix verification of pennes with laser source term-------- */
#undef __FUNCT__
#define __FUNCT__ "VerifAdjMatPennesLaser"
void VerifAdjMatPennesLaser(const PetscScalar &FEM_DT,const unsigned int &qp, const unsigned int &n_vars,
                    std::vector< unsigned int>  &n_u_dofs,
                    std::vector< const std::vector<std::vector<Real> >* > &psi,
                    std::vector< const std::vector<std::vector<RealGradient> >* > &dpsi,
                    const std::vector<Real>& JxW,
                    std::vector<Real>   &u__current,
                    std::vector<Real>   &u_previous,
                    std::vector<RealGradient> &grad_u_theta,
                    DenseSubMatrix<Number> &Kij)
{
 PetscFunctionBegin;

 libmesh_assert(n_vars == 1); // make sure solving the right system
 Real Adiff  = 0.0, Bconv = 0.0, Creac  = 0.0; 
 FORTRAN_NAME(getverifdualjac)(&FEM_DT,&u__current[0],&u_previous[0],
                                           &Adiff,&Bconv,&Creac);
 // Matrix contributions 
 for (unsigned int i=0; i<n_u_dofs[0]; i++)
   for (unsigned int j=0; j<n_u_dofs[0]; j++)
     Kij(i,j) += JxW[qp]*(
            Adiff * ( (*dpsi[0])[i][qp] * (*dpsi[0])[j][qp] ) 
                                   +
            Bconv * ( (*dpsi[0])[j][qp] * grad_u_theta[0] ) *  // note indicies 
                                           (*psi[0])[i][qp]    // are switched
                                   +
            Creac *    (*psi[0])[i][qp] *  (*psi[0])[j][qp] 
      	                          );
 PetscFunctionReturnVoid();
}
void PennesVerification::accumulateAdjointPDE( const QGauss &qrule,
                                const unsigned int &, 
                                const std::vector<Real>& JxW,
                                std::vector< SubMatrixVector > &Kij,
                                std::vector< DenseSubVector<Number> > &Fi,  
                        TransientFEMSystem &state_system,
                        TransientFEMSystem &adjoint_system) 
{
 PetscFunctionBegin; 
 AppSolve *user = AppSolve::getInstance();

 PetscScalar neg_FEM_DT = -user->get_fem_dt(),dum;
 const PetscInt constdum=0;
 // The numeric integration is combined in the same loop
 for (unsigned int qp=0; qp<qrule.n_points(); qp++)
  {
   // evaluate Solution
   this->evaluateSolutionForAdjoint(qp,AppSolve::ISTEP,state_system,adjoint_system);

   // Solution buffer 
   std::vector<PetscScalar>  solnbuffer; 

   // order is important in routines called by function pointers
   solnbuffer.push_back(  u__current[0]   ); // i_var*8+0
   solnbuffer.push_back(  u_previous[0]   ); // i_var*8+1
   solnbuffer.push_back(  u___future[0]   ); // i_var*8+2
   solnbuffer.push_back(  s__current[0]   ); // i_var*8+3
   solnbuffer.push_back(  s_previous[0]   ); // i_var*8+4
   solnbuffer.push_back(  s___future[0]   ); // i_var*8+5
   solnbuffer.push_back(  s_zero[    0]   ); // i_var*8+6
   solnbuffer.push_back(  s_full[    0]   ); // i_var*8+7

   // initialize buffers
   const PetscScalar xpoint[3] = { 0.0,0.0,0.0 };
   //const PetscScalar xpoint[3] = { q_point[qp](0), 
   //                                q_point[qp](1), 
   //                                q_point[qp](2)};
   const PetscInt Nsize = solnbuffer.size(); 
   Real AdjLoadDiff = 0.0, 
        AdjLoadConv = 0.0, 
        AdjLoadReac = 0.0, 
        Dualsource  = 0.0;
   for( unsigned int i_var = 0 ; i_var < numSystemVars ; i_var++)
    {
      /*function pointer for evaluations at the gauss points that forms 
        the load vector for the dual/adjoint problem inherent to the QOI and PDE
                               adjloadtemp:    temp-based optimize
                               getdualsource:  verification problem
                               adjloaddam__ts: two state model
                               adjloaddam_arr: arrhenius model
                               adjloadhsp:     hsp model */
      /* Negative sign nightmare:
          the jacobian time the adjoint is subtracted from the load but 
          the timestepping term in the jacobian also has a negative sign*/
      FORTRAN_NAME(getdualsource)(xpoint,&neg_FEM_DT,&AppSolve::ISTEP,&constdum,
                                  &solnbuffer[0],&Nsize,&AdjLoadDiff,
                                  &AdjLoadConv,&AdjLoadReac,&dum);
      // Sub Vector contributions  
      //    Dualsource should be zero for verif1
      for (unsigned int i=0; i<n_p_dofs[i_var]; i++)
        Fi[i_var](i) += JxW[qp]*
          ( -AdjLoadDiff*(   grad_p___future[i_var]*(*dpsi[i_var])[i][qp] )
             +(-AdjLoadConv*(grad_p___future[i_var]*grad_u_ptheta[i_var])
             + -AdjLoadReac * p___future[i_var] 
             +  Dualsource  ) * (*psi[i_var])[i][qp]               );
    }
   /*  evaluate the Adjoint matrix at the guass points 
        AdjMatPennesLaser      : adjoint matrix of pennes with laser source term
        VerifAdjMatPennesLaser : verification problem of pennes laser (adjoint) */
   VerifAdjMatPennesLaser(user->get_fem_dt(),qp,numSystemVars,n_p_dofs,
                           psi,dpsi,JxW, u__current,u_previous,
                           grad_u_mtheta,Kij[0][0]);
  } // end loop over gauss points
 PetscFunctionReturnVoid(); 
}

template< typename PDEModelSolver >
void VerificationQOI< PDEModelSolver >::
accumulateAdjointQOI( AppSolve *user,
                      const QGauss &qrule,
                      std::vector<PetscScalar> &elemParameter,
                      const unsigned int &field_id, 
                      const std::vector<Real>& JxW,
                      const std::vector<Point>& q_point,
                      std::vector< DenseSubVector<Number> > &Fi)
{
 PetscFunctionBegin; 
 EquationSystems    &EqnSystem = user->get_equation_systems();
 switch(user->get_num_exact() )
   {
    /* Various Cases for Verification problems*/
    case 1: 
     {
       // Get a reference to the NonlinearImplicitSystem we are solving
       TransientFEMSystem& state_system = 
         EqnSystem.get_system<TransientFEMSystem>("StateSystem");
     
       // Get a reference to the ExplicitSystem for ideal data
       TransientFEMSystem & ideal_system =
         EqnSystem.get_system<TransientFEMSystem>("IdealSystem");
     
       /* the following system will hold uncertainty data */
       TransientFEMSystem & ideal_uncertainty_system =
         EqnSystem.get_system<TransientFEMSystem>("IdealUncertaintySystem");

       PetscScalar neg_FEM_DT = -user->get_fem_dt();
       const PetscInt NstephiPI = user->qoiOptimizer()->Nstephi();
       // The numeric integration is combined in the same loop
       for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
         // evaluate Solution
         this->m_pdeSolver->evaluateIdealSolution(qp,AppSolve::ISTEP,NstephiPI,
                           state_system, ideal_system,ideal_uncertainty_system);
     
         // initialize buffers
         const PetscScalar xpoint[3] = { q_point[qp](0), 
                                         q_point[qp](1), 
                                         q_point[qp](2)};
         // Solution buffer 
         std::vector<PetscScalar>  solnbuffer; 
     
         Real AdjLoadDiff, AdjLoadConv, AdjLoadReac, Dualsource;
         Real Adiff      , Bconv      , Creac                  ;
         for( unsigned int i_var = 0 ; i_var < this->m_pdeSolver->n_vars() ; i_var++)
          {
            // fill solution  buffer 
            this->m_pdeSolver->fillSolutionBuffer(solnbuffer,i_var);
            const PetscInt Nsize = solnbuffer.size(); 
	    /*function pointer for evaluations at the gauss points that forms 
	      the load vector for the dual/adjoint problem inherent to the QOI
              and PDE
                                     adjloadtemp:    temp-based optimize
                                     getdualsource:  verification problem
                                     adjloaddam__ts: two state model
                                     adjloaddam_arr: arrhenius model
                                     adjloadhsp:     hsp model */
            /* Negative sign nightmare:
                the jacobian time the adjoint is subtracted from the load but 
                the timestepping term in the jacobian also has a negative sign*/
           FORTRAN_NAME(getdualsource)(xpoint,&neg_FEM_DT,&AppSolve::ISTEP,
                                       &NstephiPI,
                                       &solnbuffer[0],&Nsize,&AdjLoadDiff,
                                       &AdjLoadConv,&AdjLoadReac,&Dualsource);
           // Sub Vector contributions 
           for (unsigned int i=0; i < this->m_pdeSolver->n_p_dofs(i_var); i++)
             Fi[i_var](i) += JxW[qp]* Dualsource   *
                      this->m_pdeSolver->psi(i_var,i,qp); 
          }
        } // end loop over gauss points
     }
    break;
    default:
      this->spaceTimeQOI< PDEModelSolver >::accumulateAdjointQOI(
                            user,qrule,elemParameter,field_id,JxW,q_point,Fi);
    break;
   }
 PetscFunctionReturnVoid(); 
}

// pointer to application data stored in src/pdeSolver.cxx

/* -------------------------------------------------------------------- 
   setup
   -------------------------------------------------------------------- */ 
void PennesVerification::SetupState(AppSolve *user)
{
 PetscErrorCode info;
 PetscFunctionBegin;
 
 PDEModelBaseClass::SetupState(user);
 
 //const DroneControl &QOIInfo = *_user_app->QOIInfo;  // QOI info

 // if (QOIInfo.compobj.find("hsp_control")!=std::string::npos)
 //   {  // hsp based optimization 
 //      adjload.push_back( FORTRAN_NAME(adjloadhsp));
 //   }
 // else if (QOIInfo.compobj.find("dam_control_arr")!=std::string::npos)
 //   {
 //      //Arrhenius damage based optimization
 //      adjload.push_back( FORTRAN_NAME(adjloaddam_arr));
 //   } 
 // else if (QOIInfo.compobj.find("dam_control_two")!=std::string::npos)
 //   {
 //      //two state model damage based optimization
 //      adjload.push_back( FORTRAN_NAME(adjloaddam__ts));
 //   }  
 // else if (QOIInfo.compobj.find("calibration")!=std::string::npos)
 //   {
 //      //calibration to MRTI data 
 //      adjload.push_back( FORTRAN_NAME(adjloadtemp));
 //   }
 // else if (QOIInfo.compobj.find("computeideal_mc")!=std::string::npos)
 //   {
 //      //compute ideal field from Monte Carlo data
 //      adjload.push_back( FORTRAN_NAME(adjloadtemp));
 //   } 
 // else if (QOIInfo.compobj.find("verifprob")!=std::string::npos)
 //   { 
 //      adjload.push_back( FORTRAN_NAME(adjloadtemp));
 //   } 
 // else if (QOIInfo.compobj.find("temp_control")!=std::string::npos)
 //   { 
 //      adjload.push_back( FORTRAN_NAME(adjloadtemp));
 //   } 
 // else if (QOIInfo.compobj.find("notfound")!=std::string::npos)
 //   { 
 //      std::cout << "an objective function was not found"<<  std::endl; abort();
 //   } 
 // else 
 //   {
 //      std::cout << "unknown objective function "<<QOIInfo.compobj << std::endl; abort();
 //   } 

 //if( QOIInfo.pde.find("nonlinpennesrf")!=std::string::npos )
 //  {
 //    info = SetupPennesRF(user);
 //  }
 //else
 //  {
     info = SetupPennesLaser(user);
 //  }

  PetscScalar MaxTime = user->MaxTime();
  FORTRAN_NAME(setup_power)(&MaxTime); 
  PetscScalar dum[3],time=0.0;
  m_bodyTemp = FORTRAN_NAME(getinittemp)(dum,&time);
}
/* ------------------------------------------------------------------- 
   SetupPennesLaser  - Setup Equations Systems for a single equation
   solve of the Pennes BioHeat Transfer Equation with a laser heat source */
#undef __FUNCT__
#define __FUNCT__ "PennesVerification::SetupPennesLaser"
PetscErrorCode PennesVerification::SetupPennesLaser(AppSolve *)
{

 //set funtion pointers to routines that represent the pde at the gauss points
 //if (qoiOptimizer->pde.find("verifprob")!=std::string::npos)
 //  { 
    /* verification problem */
    // nonlinear function evaluation and bc's
    formfcngauss.push_back( 
       std::vector<FormFcnGaussType>::vector(1,FORTRAN_NAME(getverifformfcn)) );
    fcnbcgauss.push_back(         FORTRAN_NAME(getverifbc)       );
    // jacobian
    jacgauss.push_back(VerifJacPennesLaser);
    jacbcgauss.push_back( FORTRAN_NAME(getverifbc) );
    // dual problem matrix and bc's
    adjmat= VerifAdjMatPennesLaser;
    adjbcgauss.push_back( FORTRAN_NAME(getverifbcdual) );
    // special verification problem cases for dual problem load vector
    //switch(user->get_num_exact() )
    // {
    //  case 1: case 2: case 8: 
    //  adjload.at(0)=FORTRAN_NAME(getdualsource);
    //  break;
    // }
 //  } 
 //else if (qoiOptimizer->pde.find("nonlinpennesmonte")!=std::string::npos)
 //  { 
 //  /* nonlinear pennes model with monte carlo heat source */
 //    // nonlinear function evaluation and bc's
 //    formfcngauss.push_back(     
 //      std::vector<FormFcnGaussType>::vector(1,FORTRAN_NAME(nonlinpennesmonte)) );
 //    fcnbcgauss.push_back(        FORTRAN_NAME(heattransbc)        );
 //    // open file containg monte carlo source for the forward solve
 //    // jacobian and bc's
 //    jacgauss.push_back(JacPennesLaser);
 //    jacbcgauss.push_back(            FORTRAN_NAME(heattransbc)      );
 //    // dual problem matrix and bc's
 //    adjmat= AdjMatPennesLaser;
 //    adjbcgauss.push_back(            FORTRAN_NAME(heattransbc)      );
 //  } 
 //else if (qoiOptimizer->pde.find("nonlinpennesisolaser")!=std::string::npos)
 //  { 
 //  /* default to nonlinear pennes model with isotropic laser heat source for
 //   * entire domain */
 //    // nonlinear function evaluation and bc's
 //    formfcngauss.push_back(
 //      std::vector<FormFcnGaussType>::vector(1,FORTRAN_NAME(nonlinpennesisolaser)) );
 //    formfcngauss.push_back( // same residual over entire domain
 //      std::vector<FormFcnGaussType>::vector(1,FORTRAN_NAME(nonlinpennesisolaser)) );
 //    fcnbcgauss.push_back(        FORTRAN_NAME(heattransbc)           );
 //    // jacobian and bc's
 //    jacgauss.push_back(JacPennesLaser);
 //    jacgauss.push_back(JacPennesLaser); // same jac over entire domain
 //    jacbcgauss.push_back(     FORTRAN_NAME(heattransbc)      );
 //    // dual problem matrix and bc's
 //    adjmat= AdjMatPennesLaser;
 //    adjbcgauss.push_back(     FORTRAN_NAME(heattransbc)      );
 //  } 
 //else if (qoiOptimizer->pde.find("pennesisolasercooling")!=std::string::npos)
 //  { 
 //  /* default to nonlinear pennes model with isotropic laser heat source for
 //   * entire domain */
 //    // nonlinear function evaluation and bc's
 //    formfcngauss.push_back(
 //      std::vector<FormFcnGaussType>::vector(1,FORTRAN_NAME(nonlinpennesisolaser)) );
 //    formfcngauss.push_back(  // penalty method for cooling
 //      std::vector<FormFcnGaussType>::vector(1,FORTRAN_NAME(penaltytemperature)) );
 //    fcnbcgauss.push_back(        FORTRAN_NAME(heattransbc)           );
 //    // jacobian and bc's
 //    jacgauss.push_back(JacPennesLaser);
 //    jacgauss.push_back(JacPenalty); // penalty method for cooling
 //    jacbcgauss.push_back(     FORTRAN_NAME(heattransbc)      );
 //    // dual problem matrix and bc's
 //    adjmat= AdjMatPennesLaser;
 //    adjbcgauss.push_back(     FORTRAN_NAME(heattransbc)      );

 //    //overwrite for initial condition with the probe
 //    state_system.attach_init_function(initial_condition);
 //  } 

 PetscFunctionReturn(0);
}
/* Function to get ideal field */
Number pennes_ideal_field (const Point& p,
                           const Parameters& ,
                           const std::string& ,
                           const std::string& )
{
 PetscScalar xpoint[3] = { p(0), p(1), p(2) };
 return FORTRAN_NAME(getidealfield)(xpoint);
}
/* Function to compute exact solution for verification purposes.  */
Number pennes_adjoint_exact (const Point& p,
                             const Parameters& ,
                             const std::string& ,
                             const std::string& )
{

 PetscScalar xpoint[3] = { p(0), p(1), p(2) };
 PetscScalar texact = 0.0;
 PetscScalar gradexact[3] = { 0.0, 0.0, 0.0 };

 AppSolve *user = AppSolve::getInstance();
 PetscScalar FEM_DT =  user->get_fem_dt(); 
 FORTRAN_NAME(adjointexact)(xpoint,&AppSolve::ISTEP,
                            &texact,gradexact,&FEM_DT );

 return texact;

}

PetscScalar PennesVerification::residualCauchyBC(const unsigned int i_var,
                                      const Real &temperature)
{
 PetscFunctionBegin; 

 // initiliaze buffers
 AppSolve *user = AppSolve::getInstance();
 const Real time = 0.5*(2*AppSolve::ISTEP-1) * user->get_fem_dt();
 Real coeffNewton,uInfty,uFlux,uValue,penalty;

 /* get the boundary data, boundary condition 
    function pointers to one of the following routines
                  getverifbc     : verification problem bc
                  rfvoltagebc    : bc for rf voltage
                  heattransbc    : bc for heat transfer */
 //fcnbcgauss[i_var](&time,&coeffNewton,&uInfty,
 //                              &uFlux,&uValue,&penalty);
 FORTRAN_NAME(getverifbc)(&time,&coeffNewton,&uInfty,
                               &uFlux,&uValue,&penalty);

 PetscFunctionReturn( coeffNewton*(temperature-uInfty) ); 
}
PetscScalar PennesVerification::residualNeumannBC(
                 const unsigned int ,const Real &)
{
 PetscFunctionBegin; 

 // initiliaze buffers
 AppSolve *user = AppSolve::getInstance();
 const Real time = 0.5*(2*AppSolve::ISTEP-1) * user->get_fem_dt();
 Real coeffNewton,uInfty,uFlux,uValue,penalty;

 /* get the boundary data, boundary condition 
    function pointers to one of the following routines
                  getverifbc     : verification problem bc
                  rfvoltagebc    : bc for rf voltage
                  heattransbc    : bc for heat transfer */
 //fcnbcgauss[i_var](&time,&coeffNewton,&uInfty,
 //                              &uFlux,&uValue,&penalty);
 FORTRAN_NAME(getverifbc)(&time,&coeffNewton,&uInfty,
                               &uFlux,&uValue,&penalty);

 PetscFunctionReturn(uFlux); 
}
PetscScalar PennesVerification::jacobianCauchyBC(const unsigned int )
{
 PetscFunctionBegin;

 // initiliaze buffers
 AppSolve *user = AppSolve::getInstance();
 const Real time = 0.5*(2*AppSolve::ISTEP-1) * user->get_fem_dt();
 Real coeffNewton,uInfty,uFlux,uValue,penalty;

 /* get the boundary data, boundary condition 
    function pointers to one of the following routines
                  getverifbc     : verification problem bc
                  rfvoltagebc    : bc for rf voltage
                  heattransbc    : bc for heat transfer */
 //fcnbcgauss[i_var](&time,&coeffNewton,&uInfty,
 //                              &uFlux,&uValue,&penalty);
 FORTRAN_NAME(getverifbc)(&time,&coeffNewton,&uInfty,
                               &uFlux,&uValue,&penalty);

 // special case for steady state offset m_theta
 if( user->get_num_exact() == 1 )
  {
   coeffNewton = 2.0 *  coeffNewton;
  }
 PetscFunctionReturn( coeffNewton );
}

void PennesVerification::gradSolutionBufferSetup(
                                   std::vector<PetscScalar>  &solnbuffer, 
                                   const unsigned int &field_id)
{
 PetscFunctionBegin; 
 // the state solution computed during the state solve is used
 // order is important in routines called by function pointers
 solnbuffer.push_back(   u__current[0]    ); // i_var * 0
 solnbuffer.push_back(   u_previous[0]    ); // i_var * 1
 solnbuffer.push_back( grad_u_mtheta[0](0) ); // i_var * 2
 solnbuffer.push_back( grad_u_mtheta[0](1) ); // i_var * 3
 solnbuffer.push_back( grad_u_mtheta[0](2) ); // i_var * 4

 // order is important in routines called by function pointers
 solnbuffer.push_back(      p__now[0]    );//n_vars_state*5 + i_var * 0
 solnbuffer.push_back( grad_p__now[0](0) );//n_vars_state*5 + i_var * 1
 solnbuffer.push_back( grad_p__now[0](1) );//n_vars_state*5 + i_var * 2
 solnbuffer.push_back( grad_p__now[0](2) );//n_vars_state*5 + i_var * 3

 // set field
 const PetscInt FieldID = field_id;
 FORTRAN_NAME(setfieldid)( &FieldID );
 PetscFunctionReturnVoid();
}
/*---------------------------------------------------------------------------*/
PetscScalar PennesVerification::dpde_dw_0( const unsigned int &field_id, const Point &q_point)
{
 PetscFunctionBegin; 

 // solution buffer
 std::vector<PetscScalar>  solnbuffer; 
 // gauss point buffer
 const PetscScalar xpoint[3] = { q_point(0), q_point(1), q_point(2)};
 this->gradSolutionBufferSetup(solnbuffer,field_id);
 const PetscInt Nsize=solnbuffer.size();
 AppSolve *user = AppSolve::getInstance();
 PetscScalar FEM_DT =  user->get_fem_dt(); 

 PetscFunctionReturn(
       FORTRAN_NAME(dwdw0)(xpoint,&AppSolve::ISTEP,&FEM_DT,
                           &solnbuffer[0],&Nsize));
}
/*--------------------------------------------------------------------------------*/
PetscScalar PennesVerification::dpde_dx_0( const unsigned int &field_id, const Point &q_point)
{
 PetscFunctionBegin; 

 // solution buffer
 std::vector<PetscScalar>  solnbuffer; 
 // gauss point buffer
 const PetscScalar xpoint[3] = { q_point(0), q_point(1), q_point(2)};
 this->gradSolutionBufferSetup(solnbuffer,field_id);
 const PetscInt Nsize=solnbuffer.size();
 AppSolve *user = AppSolve::getInstance();
 PetscScalar FEM_DT =  user->get_fem_dt(); 

 PetscFunctionReturn(
       FORTRAN_NAME(dqlaserdx)(xpoint,&AppSolve::ISTEP,&FEM_DT,
                           &solnbuffer[0],&Nsize));
}
/*--------------------------------------------------------------------------------*/
PetscScalar PennesVerification::dpde_dy_0( const unsigned int &field_id, const Point &q_point)
{
 PetscFunctionBegin; 

 // solution buffer
 std::vector<PetscScalar>  solnbuffer; 
 // gauss point buffer
 const PetscScalar xpoint[3] = { q_point(0), q_point(1), q_point(2)};
 this->gradSolutionBufferSetup(solnbuffer,field_id);
 const PetscInt Nsize=solnbuffer.size();
 AppSolve *user = AppSolve::getInstance();
 PetscScalar FEM_DT =  user->get_fem_dt(); 

 PetscFunctionReturn(
       FORTRAN_NAME(dqlaserdy)(xpoint,&AppSolve::ISTEP,&FEM_DT,
                           &solnbuffer[0],&Nsize));
}
/*--------------------------------------------------------------------------------*/
PetscScalar PennesVerification::dpde_dz_0( const unsigned int &field_id, const Point &q_point)
{
 PetscFunctionBegin; 

 // solution buffer
 std::vector<PetscScalar>  solnbuffer; 
 // gauss point buffer
 const PetscScalar xpoint[3] = { q_point(0), q_point(1), q_point(2)};
 this->gradSolutionBufferSetup(solnbuffer,field_id);
 const PetscInt Nsize=solnbuffer.size();
 AppSolve *user = AppSolve::getInstance();
 PetscScalar FEM_DT =  user->get_fem_dt(); 

 PetscFunctionReturn(
       FORTRAN_NAME(dqlaserdz)(xpoint,&AppSolve::ISTEP,&FEM_DT,
                           &solnbuffer[0],&Nsize));
}
/*--------------------------------------------------------------------------------*/
PetscScalar PennesVerification::dpde_dw_2( const unsigned int &field_id, const Point &q_point)
{
 PetscFunctionBegin; 

 // solution buffer
 std::vector<PetscScalar>  solnbuffer; 
 // gauss point buffer
 const PetscScalar xpoint[3] = { q_point(0), q_point(1), q_point(2)};
 this->gradSolutionBufferSetup(solnbuffer,field_id);
 const PetscInt Nsize=solnbuffer.size();
 AppSolve *user = AppSolve::getInstance();
 PetscScalar FEM_DT =  user->get_fem_dt(); 

 PetscFunctionReturn(
       FORTRAN_NAME(dwdw2)(xpoint,&AppSolve::ISTEP,&FEM_DT,
                           &solnbuffer[0],&Nsize));
}
/*--------------------------------------------------------------------------------*/
PetscScalar PennesVerification::dpde_dw_N( const unsigned int &field_id, const Point &q_point)
{
 PetscFunctionBegin; 

 // solution buffer
 std::vector<PetscScalar>  solnbuffer; 
 // gauss point buffer
 const PetscScalar xpoint[3] = { q_point(0), q_point(1), q_point(2)};
 this->gradSolutionBufferSetup(solnbuffer,field_id);
 const PetscInt Nsize=solnbuffer.size();
 AppSolve *user = AppSolve::getInstance();
 PetscScalar FEM_DT =  user->get_fem_dt(); 

 PetscFunctionReturn(
       FORTRAN_NAME(dwdwn)(xpoint,&AppSolve::ISTEP,&FEM_DT,
                           &solnbuffer[0],&Nsize));
}
/*--------------------------------------------------------------------------------*/
PetscScalar PennesVerification::dpde_dw_I( const unsigned int &field_id, const Point &q_point)
{
 PetscFunctionBegin; 

 // solution buffer
 std::vector<PetscScalar>  solnbuffer; 
 // gauss point buffer
 const PetscScalar xpoint[3] = { q_point(0), q_point(1), q_point(2)};
 this->gradSolutionBufferSetup(solnbuffer,field_id);
 const PetscInt Nsize=solnbuffer.size();
 AppSolve *user = AppSolve::getInstance();
 PetscScalar FEM_DT =  user->get_fem_dt(); 

 PetscFunctionReturn(
       FORTRAN_NAME(dwdwi)(xpoint,&AppSolve::ISTEP,&FEM_DT,
                           &solnbuffer[0],&Nsize));
}
/*--------------------------------------------------------------------------------*/
PetscScalar PennesVerification::dpde_dw_D( const unsigned int &field_id, const Point &q_point)
{
 PetscFunctionBegin; 

 // solution buffer
 std::vector<PetscScalar>  solnbuffer; 
 // gauss point buffer
 const PetscScalar xpoint[3] = { q_point(0), q_point(1), q_point(2)};
 this->gradSolutionBufferSetup(solnbuffer,field_id);
 const PetscInt Nsize=solnbuffer.size();
 AppSolve *user = AppSolve::getInstance();
 PetscScalar FEM_DT =  user->get_fem_dt(); 

 PetscFunctionReturn(
       FORTRAN_NAME(dwdwd)(xpoint,&AppSolve::ISTEP,&FEM_DT,
                           &solnbuffer[0],&Nsize));
}
/*--------------------------------------------------------------------------*/
PetscScalar PennesVerification::dpde_dk_0( const unsigned int &field_id, const Point &q_point)
{
 PetscFunctionBegin; 

 // solution buffer
 std::vector<PetscScalar>  solnbuffer; 
 // gauss point buffer
 const PetscScalar xpoint[3] = { q_point(0), q_point(1), q_point(2)};
 this->gradSolutionBufferSetup(solnbuffer,field_id);
 const PetscInt Nsize=solnbuffer.size();
 AppSolve *user = AppSolve::getInstance();
 PetscScalar FEM_DT =  user->get_fem_dt(); 

 PetscFunctionReturn(
       FORTRAN_NAME(dkdk0)(xpoint,&AppSolve::ISTEP,&FEM_DT,
                           &solnbuffer[0],&Nsize));
}
/*--------------------------------------------------------------------------------*/
PetscScalar PennesVerification::dpde_dk_1( const unsigned int &field_id, const Point &q_point)
{
 PetscFunctionBegin; 

 // solution buffer
 std::vector<PetscScalar>  solnbuffer; 
 // gauss point buffer
 const PetscScalar xpoint[3] = { q_point(0), q_point(1), q_point(2)};
 this->gradSolutionBufferSetup(solnbuffer,field_id);
 const PetscInt Nsize=solnbuffer.size();
 AppSolve *user = AppSolve::getInstance();
 PetscScalar FEM_DT =  user->get_fem_dt(); 

 PetscFunctionReturn(
       FORTRAN_NAME(dkdk1)(xpoint,&AppSolve::ISTEP,&FEM_DT,
                           &solnbuffer[0],&Nsize));
}
/*--------------------------------------------------------------------------------*/
PetscScalar PennesVerification::dpde_dmu_a(const unsigned int &field_id, const Point &q_point)
{
 PetscFunctionBegin; 

 // solution buffer
 std::vector<PetscScalar>  solnbuffer; 
 // gauss point buffer
 const PetscScalar xpoint[3] = { q_point(0), q_point(1), q_point(2)};
 this->gradSolutionBufferSetup(solnbuffer,field_id);
 const PetscInt Nsize=solnbuffer.size();
 AppSolve *user = AppSolve::getInstance();
 PetscScalar FEM_DT =  user->get_fem_dt(); 

 PetscFunctionReturn(
       FORTRAN_NAME(dqlaserdmu_a)(xpoint,&AppSolve::ISTEP,&FEM_DT,
                           &solnbuffer[0],&Nsize));
}
/*--------------------------------------------------------------------------------*/
PetscScalar PennesVerification::dpde_dmu_s(const unsigned int &field_id, const Point &q_point)
{
 PetscFunctionBegin; 

 // solution buffer
 std::vector<PetscScalar>  solnbuffer; 
 // gauss point buffer
 const PetscScalar xpoint[3] = { q_point(0), q_point(1), q_point(2)};
 this->gradSolutionBufferSetup(solnbuffer,field_id);
 const PetscInt Nsize=solnbuffer.size();
 AppSolve *user = AppSolve::getInstance();
 PetscScalar FEM_DT =  user->get_fem_dt(); 

 PetscFunctionReturn(
       FORTRAN_NAME(dqlaserdmu_s)(xpoint,&AppSolve::ISTEP,&FEM_DT,
                           &solnbuffer[0],&Nsize));
}
// ------------------------------------------------------------
// template instantiations
template class PennesInverseSystem<PennesVerification>;
template class VerificationQOI< PennesInverseSystem < PennesVerification > >;
template class spaceTimeFDVerification< PennesInverseSystem < PennesVerification > >;
