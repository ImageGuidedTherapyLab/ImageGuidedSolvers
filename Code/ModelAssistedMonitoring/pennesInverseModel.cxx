// libmesh includes
#include "libmesh.h"
#include "equation_systems.h"
#include "fe.h"
#include "fe_interface.h"
#include "quadrature_gauss.h"
#include "getpot.h"
// libmesh verification framework
#include "exact_solution.h"


// The nonlinear solver and system we will be using
#include "nonlinear_implicit_system.h"
#include "linear_implicit_system.h"
#include "mesh.h"

// dddas includes
#include "petsc.h"
#include "femInterface.h"
#include "transient_fem_system.h"
#include "pennesInverseModel.h"
#include "tttkUtilities.h"

/** 
 * @file pennesSolver.cxx
 * 
 * this file contains the pre-FEMSystem routines
 * base clase for pde constitutitve data for an arbitrary number of variables.
 * contains basic methods for solution evaluation and boundary condition
 * evaluations assuming that systems with multiple variables have uncoupled
 * boundary conditions
 */ 

PennesInverseModel::PennesInverseModel(GetPot &controlfile,libMesh::MeshBase &mesh):
          PennesStandardDiffusionApproximation(controlfile,          mesh)
{
}

/*----------------------------------------------------------------------*/
void PennesInverseModel::SetupOptimizationVariables( libMesh::MeshBase &mesh,
                 std::vector<optimizationParameter*> &Parameters)
{
 PetscFunctionBegin; 
                            
 PDEModelBaseClass::SetupOptimizationVariables( mesh,Parameters);

 if(  w_0.Optimize  ){  w_0.OneVarSetup(mesh); Parameters.push_back( &w_0  ); } 
 if(  k_0.Optimize  ){  k_0.OneVarSetup(mesh); Parameters.push_back( &k_0  ); } 
 if(  x_0.Optimize  ){  x_0.OneVarSetup(mesh); Parameters.push_back( &x_0  ); } 
 if(  y_0.Optimize  ){  y_0.OneVarSetup(mesh); Parameters.push_back( &y_0  ); } 
 if(  z_0.Optimize  ){  z_0.OneVarSetup(mesh); Parameters.push_back( &z_0  ); } 
 if( mu_a.Optimize  ){ mu_a.OneVarSetup(mesh); Parameters.push_back( &mu_a ); } 
 if( mu_s.Optimize  ){ mu_s.OneVarSetup(mesh); Parameters.push_back( &mu_s ); } 
 if(Power.Optimize  ){Power.OneVarSetup(mesh); Parameters.push_back( &Power); } 

 // allocate space for second derivs functers, most should return zero
 const unsigned int nParam = Parameters.size(); 
 d2pde_dmi_dmj.resize(nParam, std::vector< dpde_dmMemFn >::vector (nParam, NULL) );
 
 // accumulate element hessian sensitivities 
 std::vector<optimizationParameter*>::iterator rowIter,columnIter;
 for( rowIter=Parameters.begin(); rowIter!=Parameters.end(); rowIter++)
  {
   int iii = distance(Parameters.begin(),rowIter);
   optimizationParameter* iiiParam = *rowIter;
   // exploit symmetry  
   for( columnIter=rowIter; columnIter!=Parameters.end(); columnIter++)
    {
     int jjj = distance(Parameters.begin(),columnIter);
     optimizationParameter* jjjParam = *columnIter;
     // default returns zero
     dpde_dmMemFn tmpMemFcn = &PDEModelBaseClass::d2pde_d2m;
     if(      iiiParam->name().find("mu_a") !=std::string::npos
                              &&
              jjjParam->name().find("mu_a") !=std::string::npos )
                                   tmpMemFcn = &PDEModelBaseClass::d2pde_d2mu_a;
     else if( iiiParam->name().find("mu_s") !=std::string::npos
                              &&
              jjjParam->name().find("mu_s") !=std::string::npos )
                                   tmpMemFcn = &PDEModelBaseClass::d2pde_d2mu_s;
     // due to order of push back mu_a should be the iiiParam
     else if( iiiParam->name().find("mu_a") !=std::string::npos
                              &&
              jjjParam->name().find("mu_s") !=std::string::npos )
                                   tmpMemFcn = &PDEModelBaseClass::d2pde_dmu_a_dmu_s;
     this->d2pde_dmi_dmj[iii][jjj] = tmpMemFcn ;
     this->d2pde_dmi_dmj[jjj][iii] = tmpMemFcn ; // exploit symmetry  
    }
  }

 PetscFunctionReturnVoid(); 
}

void PennesInverseModel::accumulateAdjointPDE( const QGauss &qrule,
                                const unsigned int &field_id, 
                                const std::vector<Real>& JxW,
                                std::vector< SubMatrixVector > &Kij,
                                std::vector< DenseSubVector<Number> > &Fi,  
                        TransientFEMSystem &state_system,
                        TransientFEMSystem &adjoint_system) 
{
 PetscFunctionBegin; 
 FiniteElementInterface *user = FiniteElementInterface::getInstance();

 // The numeric integration is combined in the same loop
 for (unsigned int qp=0; qp<qrule.n_points(); qp++)
   {
    // evaluate Solution
    this->evaluateSolutionForAdjoint(qp,FiniteElementInterface::ISTEP,state_system,adjoint_system);

    const Real u_mtheta = (1.0-m_theta)*u_previous[0]+m_theta*u__current[0];
    const Real u_ptheta = (1.0-m_theta)*u__current[0]+m_theta*u___future[0];
    // FIXME - Negative sign nightmare... pulling pde terms to the RHS
    const Real AdjLoadDiff = thermal_cond(field_id,u_ptheta)*m_theta;
    const Real AdjLoadConv = dthermal_conddu(u_ptheta)      *m_theta;
    const Real AdjLoadReac =-rho*specific_heat/user->get_fem_dt() //transient term
                                   +    // perfusion term (check the SIGN)
                             perfusion(field_id,u_ptheta) 
                                           * bloodspecificheat*m_theta
                                   +    // perfusion derivative    
                             dperfusiondu(u_ptheta)
                                        * bloodspecificheat*(u_ptheta-u_artery);
    // Sub Vector contributions 
    for (unsigned int i=0; i<n_p_dofs[0]; i++)
     {
        Fi[0](i) += JxW[qp]*
          (  -AdjLoadDiff*(   grad_p___future[0]*(*dpsi[0])[i][qp] )
           +(-AdjLoadConv*(grad_p___future[0]*grad_u_ptheta[0])
           + -AdjLoadReac * p___future[0] ) * (*psi[0])[i][qp]  );
      // Sub Matrix contributions 
      for (unsigned int j=0; j<n_p_dofs[0]; j++)
        Kij[0][0](i,j) += JxW[qp]*(
            thermal_cond(field_id,u_mtheta)*m_theta 
                   * ( (*dpsi[0])[i][qp] * (*dpsi[0])[j][qp] ) 
                                   +
            dthermal_conddu(u_mtheta) * m_theta
                   * ( (*dpsi[0])[j][qp] * grad_u_mtheta[0] ) * // note indicies
                                           (*psi[0])[i][qp]     // are switched
                                   +
                         (rho*specific_heat/user->get_fem_dt() //transient term
                                   +    // perfusion term (check the SIGN)
                          perfusion(field_id,u_mtheta)
                                          *bloodspecificheat*m_theta
                                   +    // perfusion derivative    
                          dperfusiondu(u_mtheta) * m_theta
                                          *bloodspecificheat*(u_mtheta-u_artery)
                         ) *(*psi[0])[i][qp] *  (*psi[0])[j][qp] 
      	                          );
     }
   } // end loop over gauss points
 PetscFunctionReturnVoid(); 
}
//accumulate the load for the sensitivity solve  useing previous sensitivity
void PennesInverseModel::accumulatePreviousSensitivity(const QGauss &qrule, 
                                const unsigned int &field_id, 
                                const std::vector<Real>& JxW,
                                const std::vector<Point>& ,
                                std::vector< DenseSubVector<Number> > &Fi,
                     TransientFEMSystem &state_system,
                     TransientFEMSystem  &sensitivity_system)
{
 PetscFunctionBegin; 

 // floor fem time step by ISTEPS_PER_IDEAL  to get power step
 FiniteElementInterface *user = FiniteElementInterface::getInstance();
 //PetscInt idpow= user->get_id_power();

 for (unsigned int qp=0; qp<qrule.n_points(); qp++)
  {
   // evaluate Solution used in the dpde_dm functer
   this->evaluateSolutionForState(qp,state_system);
   this->evaluateSolutionForSensitivity(qp,sensitivity_system);

   const Real u_theta = (1.0-m_theta)*u_previous[0]+m_theta*u__current[0];
   /*   evaluate the sensitivity load at the guass points */
   // loop over shape functions
   for (unsigned int i=0; i<n_u_dofs[0]; i++)
    {
      Fi[0](i) += JxW[qp] * (
             //transient term
             ( rho*specific_heat/user->get_fem_dt() * 
               du_previous[0] 
                       -    // perfusion term (check the SIGN)
               perfusion(field_id,u_theta) * (1.0 - m_theta) * du_previous[0]
                  * bloodspecificheat  
                       -    // perfusion derivative    
               dperfusiondu(u_theta) * (1.0 - m_theta) * du_previous[0]
                  * bloodspecificheat * (u_theta-u_artery) ) * (*phi[0])[i][qp]
                       -    // diffusion term
             thermal_cond(field_id,u_theta)*(1.0-m_theta)
                             * grad_du_previous[0] * (*dphi[0])[i][qp]
                       -    // diffusion derivative (non-symmetric)
             dthermal_conddu(u_theta) * (1.0-m_theta) * du_previous[0]
                       * grad_u_mtheta[0] * (*dphi[0])[i][qp] 
                            );
    }
  } // end of the quadrature point (qp-loop) over element interior

 PetscFunctionReturnVoid(); 
}
/*----------------------------------------------------------------------------*/
PetscScalar PennesInverseModel::dpde_derr( const unsigned int &field_id, const Point &qpoint)
{
 PetscFunctionBegin; 
 FiniteElementInterface *user = FiniteElementInterface::getInstance();
 PetscInt idpow= user->get_id_power();

 const Real u_theta = (1.0-m_theta)*u_previous[0]+m_theta*u__current[0];
 Real residual = -( //MINUS
                    (rho*specific_heat*(u__current[0]-u_previous[0])/user->get_fem_dt() //transient term
                           +    // perfusion term (check the SIGN)
                  perfusion(field_id,u_theta) *bloodspecificheat*(u_theta-u_artery)
                           -    // source term 
                  mu_a[field_id] * laserFluence(field_id,qpoint, idpow )
                  )* p__now[0]
                           +    // diffusion term
                 thermal_cond(field_id,u_theta)*grad_u_mtheta[0]*grad_p__now[0]
                  );
 PetscFunctionReturn(residual);
}
/*----------------------------------------------------------------------------*/
PetscScalar PennesInverseModel::dpde_dw_0( const unsigned int &, const Point &)
{
 PetscFunctionBegin; 
 const Real u_mtheta = (1.0-m_theta)*u_previous[0]+m_theta*u__current[0];
 PetscFunctionReturn(1.0*bloodspecificheat*(u_mtheta-u_artery)* p__now[0]);
}
/*-------------------------------------------------------------------------*/
PetscScalar PennesInverseModel::dpde_dk_0( const unsigned int &, const Point &)
{
 PetscFunctionBegin; 

 PetscFunctionReturn( 1.0*(grad_u_mtheta[0]* grad_p__now[0]) );
}
/*-------------------------------------------------------------------------*/
PetscScalar PennesInverseModel::dpde_dmu_a( const unsigned int &field_id, const Point &q_point)
{
 PetscFunctionBegin; 

 // floor fem time step by ISTEPS_PER_IDEAL  to get power step
 FiniteElementInterface *user = FiniteElementInterface::getInstance();
 PetscInt idpow= user->get_id_power();

 PetscFunctionReturn(-dqlaserdmu_a(field_id,q_point,idpow ) * p__now[0] );
                  
}
/*-------------------------------------------------------------------------*/
PetscScalar PennesInverseModel::dpde_dmu_s( const unsigned int &field_id, const Point &q_point)
{
 PetscFunctionBegin; 

 // floor fem time step by ISTEPS_PER_IDEAL  to get power step
 FiniteElementInterface *user = FiniteElementInterface::getInstance();
 PetscInt idpow= user->get_id_power();

 PetscFunctionReturn(-dqlaserdmu_s(field_id,q_point,idpow ) * p__now[0] );
}
/*-------------------------------------------------------------------------*/
PetscScalar PennesInverseModel::dpde_dprobeTemp( const unsigned int &, const Point &)
{
 PetscFunctionBegin; 

 Real dpdeddirichlet = 0.0;

 //const Real u_mtheta = (1.0-m_theta)*u_previous[0]+m_theta*u__current[0];

 //// shape functions should be zeroed for non dirichelet nodes
 //for (unsigned int i=0; i<n_u_dofs[0]; i++)
 //   {
 //   dpdeddirichlet = 
 //     thermal_cond(field_id,u_mtheta) * (*dphi[0])[i][qp] * grad_p__now[0]
 //               + 
 //     perfusion(field_id,u_mtheta) * bloodspecificheat
 //                                 * (*phi[0])[i][qp] * p__now[0] ;
 //   }

 PetscFunctionReturn( dpdeddirichlet  );
}
/*-------------------------------------------------------------------------*/
PetscScalar PennesInverseModel::d2pde_du_dk_0( const unsigned int &, const Point &)
{
 PetscFunctionBegin; 

 PetscFunctionReturn( 1.0*(grad_du_mtheta[0]* grad_p__now[0]) );
}
/*-------------------------------------------------------------------------*/
PetscScalar PennesInverseModel::d2pde_du_dw_0( const unsigned int &, const Point &)
{
 PetscFunctionBegin; 

 const Real du_mtheta = (1.0-m_theta)*du_previous[0]
                            +m_theta *du__current[0];
 PetscFunctionReturn(1.0*bloodspecificheat*du_mtheta*p__now[0]);
}
/*-------------------------------------------------------------------------*/
PetscScalar PennesInverseModel::dqoidu_dudm( const unsigned int &)
{
 PetscFunctionBegin; 

 //const Real du_mtheta = (1.0-m_theta)*du_previous[0]
 //                           +m_theta *du__current[0];
 //const Real  u_theta  = (1.0-m_theta)* u_previous[0]
 //                           +m_theta * u__current[0];
 PetscScalar grad = this->getMass() * (u__current[0]-s__current[0]) /
                    ds_current[0] / ds_current[0] * du__current[0] ;
 //PetscScalar grad = this->getMass() * p__now[0]  / user->get_fem_dt()  
 //                    * (du__current[idParam][0] - du_previous[idParam][0]) 
 //                       + 
 //                       perfusion(field_id,u_theta) * 
 //                         bloodspecificheat * du_mtheta * p__now[0] 
 //                       + 
 //                       thermal_cond(field_id,u_theta)  * 
 //                         grad_du_mtheta[idParam][0] * grad_p__now[0];
 PetscFunctionReturn( grad );
}
/*-------------------------------------------------------------------------*/
PetscScalar PennesInverseModel::d2pde_d2mu_a( const unsigned int &field_id, const Point &q_point)
{
 PetscFunctionBegin; 

 // floor fem time step by ISTEPS_PER_IDEAL  to get power step
 FiniteElementInterface *user = FiniteElementInterface::getInstance();
 PetscInt idpow= user->get_id_power();

 PetscFunctionReturn(-d2qlaserd2mu_a(field_id,q_point,idpow ) * p__now[0] );
                  
}
/*-------------------------------------------------------------------------*/
PetscScalar PennesInverseModel::d2pde_d2mu_s( const unsigned int &field_id, const Point &q_point)
{
 PetscFunctionBegin; 

 // floor fem time step by ISTEPS_PER_IDEAL  to get power step
 FiniteElementInterface *user = FiniteElementInterface::getInstance();
 PetscInt idpow= user->get_id_power();

 PetscFunctionReturn(-d2qlaserd2mu_s(field_id,q_point,idpow ) * p__now[0] );
                  
}
/*-------------------------------------------------------------------------*/
PetscScalar PennesInverseModel::d2pde_dmu_a_dmu_s( const unsigned int &field_id, const Point &q_point)
{
 PetscFunctionBegin; 

 // floor fem time step by ISTEPS_PER_IDEAL  to get power step
 FiniteElementInterface *user = FiniteElementInterface::getInstance();
 PetscInt idpow= user->get_id_power();

 PetscFunctionReturn(-d2qlaserdmu_a_dmu_s(field_id,q_point,idpow ) * p__now[0] );
                  
}
/* -------------------------------------------------------------------- 
   Direct coupled fluence solver
   -------------------------------------------------------------------- */ 

void PennesFluencePDE::printSelf(std::ostream& os)
{
  PennesInverseModel::printSelf(os);
  printStdVector<Real>(os, "         alpha[",      alpha);
  printStdVector<Real>(os, "    speedLight[", speedLight);
  printStdVector<Real>(os, "initialFluence[", initialFluence);
}
PennesFluencePDE::PennesFluencePDE(GetPot &controlfile,libMesh::MeshBase &mesh):
                                 PennesInverseModel(controlfile,mesh)
{
 // one temperature variable and one variable for the fluence
 numSystemVars = 2;
 FiniteElementInterface *user = FiniteElementInterface::getInstance();

 const Real speedLightVacuum = 2.98e8; // m/s
 // setup domain parameters
 for(int subdomain_id = 0 ; 
         subdomain_id < user->get_num_elem_blk() ; subdomain_id ++)
   {
    // speed of light in tissue is speed of light 
    // in vacuum divided by refraction index
    speedLight.push_back(
           speedLightVacuum/controlfile("optical/refraction",1.37) ); 
    alpha.push_back( speedLight[subdomain_id]/
             ( 3.0*(mu_a[subdomain_id]+(1.0-anfact)*mu_s[subdomain_id]) ) );
   }

 // neuman bc data
 m_NeumannFlux.push_back( controlfile("bc/phi_flux",0.0) ) ; 

 // cauchy bc data
 m_newton_coeff.push_back( controlfile("bc/phi_newton_coeff",1.0)); 
 m_u_infty.push_back(  controlfile("bc/phi_infty",20.0)); 

 // initial condition for temperature should be set
 //   get IC for fluence
 for ( unsigned int Ii = 0 ; Ii < user->get_num_elem_blk()  ; Ii++ ) 
                              initialFluence.push_back(0.0);
 initialFluence.at(user->get_num_elem_blk()  - 1)  
           = controlfile("initial_condition/z_init",100.0); 

 // default initial condition is domain wise constant
 //std::vector<InitialConditionFncterType>  
 //    fluence_IC(user->get_num_elem_blk() ,PennesFluencePDEInitialFluence);
 //InitValues.push_back(fluence_IC);

 // domains above the dirichletID will have dirichlet BC's applied
 std::vector<PetscTruth>  fluenceDirichletBC(user->get_num_elem_blk() ,PETSC_FALSE);

 // default is dirichlet data on final domain
 //PetscInt fluenceDirichletID 
 //             = controlfile("bc/phi_dirichletid",user->get_num_elem_blk() -1); 
 //libmesh_assert(fluenceDirichletID > 0 ); 
 //for ( unsigned int Ii = fluenceDirichletID ; Ii < user->get_num_elem_blk()  ; Ii++ ) 
 //                                       fluenceDirichletBC[Ii] = PETSC_TRUE;

 //dirichletBC.push_back(fluenceDirichletBC);
}
void PennesFluencePDE::verifySolution( EquationSystems &eqnSystem )
{
  PetscFunctionBegin; 

  // call base class
  //this->PennesInverseModel::verifySolution(user);

  // compute error
  libMesh::ExactSolution exactState( eqnSystem );
  //exactState.attach_exact_value(pdeExactSolution);
  exactState.compute_error("StateSystem","u1");
  // get forward problem error
  PetscScalar solutionError = 
              exactState.l2_error("StateSystem", "u1");
  std::cout << "error: " << solutionError << std::endl;

  // printout error result
  //if(FiniteElementInterface::ISTEP) FiniteElementInterface::indicateError(solutionError,errorTol);

  PetscFunctionReturnVoid(); 
}

// Wrapper to evaluate initial condition
Number PennesFluencePDEInitialFluence(const Point& ,
                                      const Parameters& parameters,
                                      const std::string& ,
                                      const std::string& )
{ // get the pointer to the solver
  PennesInverseModel* pdeSolver = parameters.get<PennesInverseModel*>("pdeSolver"); 
  PennesFluencePDE* PennesFluenceSolver= dynamic_cast<PennesFluencePDE*>(pdeSolver);
  unsigned int subdomain_id =  parameters.get<unsigned int>("subdomain_id"); 
  return PennesFluenceSolver->getInitFluence(subdomain_id);
}  	 	 

Number PennesFluencePDE::exactSolution(const Point& p,
                                       const Parameters& ,
                                       const std::string& )
{
  FiniteElementInterface *user = FiniteElementInterface::getInstance();
  PetscScalar time = user->getTime(FiniteElementInterface::ISTEP);
  int n = 1;
  Number L = 10.;
  Number x = p(0);

  Number solutionValue = 
      100. * std::exp(- (alpha[0] *std::pow((n*libMesh::pi/L),2) 
                      + speedLight[0] * mu_a[0] ) * time ) 
           * std::sin( n * libMesh::pi / L * x);
  return solutionValue;  	 	 
}
//void PennesFluencePDE::accumulateResidual(const QGauss &qrule, 
//                                          const unsigned int &field_id, 
//                                          const std::vector<Real>& JxW,
//                                          const std::vector<Point>& ,
//                                          std::vector< DenseSubVector<Number> > &Ri,  
//                                          TransientFEMSystem &system)
//{
// PetscFunctionBegin;
// FiniteElementInterface *user = FiniteElementInterface::getInstance();
//
// for (unsigned int qp=0; qp<qrule.n_points(); qp++)
//  {
//   // evaluate Solution
//   this->evaluateSolutionForState(qp,system);
//
//   const Real u_theta = (1.0-m_theta)*u_previous[0]+m_theta*u__current[0];
//   const Real z_theta = (1.0-m_theta)*u_previous[1]+m_theta*u__current[1];
//   // accumulate the residual for temperature
//   for (unsigned int i=0; i<n_u_dofs[0]; i++)
//     {
//       Ri[0](i) += JxW[qp]*(
//         (rho*specific_heat*(u__current[0]-u_previous[0])/user->get_fem_dt() //transient term
//                   +    // perfusion term (check the SIGN)
//          perfusion(field_id,u_theta) *bloodspecificheat*(u_theta-u_artery)
//                   -    // source term 
//          mu_a[field_id] * z_theta 
//         )* (*phi[0])[i][qp]
//                   +    // diffusion term
//         thermal_cond(field_id,u_theta)*grad_u_mtheta[0]*(*dphi[0])[i][qp]
//                           );
//     }
//
//   // accumulate the residual for the fluence
//   for (unsigned int i=0; i<n_u_dofs[1]; i++)
//     {
//       Ri[1](i) += JxW[qp]*(
//                 ( (u__current[1]-u_previous[1])/user->get_fem_dt() //transient term
//                           +    // source term
//                  speedLight[field_id]* mu_a[field_id] * z_theta 
//                 )* (*phi[1])[i][qp]
//                           +    // diffusion term
//                 alpha[field_id]*grad_u_mtheta[1]*(*dphi[1])[i][qp]
//                           );
//     }
//  } // end of the quadrature point (qp-loop) over element interior
//
// PetscFunctionReturnVoid();
//}
//
//void PennesFluencePDE::accumulateJacobian( const QGauss &qrule,
//            const unsigned int &field_id, 
//            const std::vector<Real>& JxW,
//            std::vector< SubMatrixVector > &Kij,
//            TransientFEMSystem &system)  
//{
// PetscFunctionBegin;
// FiniteElementInterface *user = FiniteElementInterface::getInstance();
//
// for (unsigned int qp=0; qp<qrule.n_points(); qp++)
//  {
//   // evaluate Solution
//   this->evaluateSolutionForState(qp,system);
//
//   const Real u_theta = (1.0-m_theta)*u_previous[0]+m_theta*u__current[0];
//   // Juu and Juz contributions
//   for (unsigned int i=0; i<n_u_dofs[0]; i++)
//    {
//     for (unsigned int j=0; j<n_u_dofs[0]; j++)
//      {
//       Kij[0][0](i,j) += JxW[qp]*(
//              thermal_cond(field_id,u_theta)*m_theta  // diffusion term
//                    * ( (*dphi[0])[i][qp] * (*dphi[0])[j][qp] ) 
//                                     +    // diffusion derivative (non-symmetric)
//              dthermal_conddu(u_theta)
//                    * ( (*dphi[0])[i][qp] * grad_u_mtheta[0] ) * 
//                                                   (*phi[0])[j][qp]
//                                     +
//              (rho*specific_heat/user->get_fem_dt() //transient term
//                        +    // perfusion term (check the SIGN)
//               perfusion(field_id,u_theta)
//                               *bloodspecificheat*m_theta
//                        +    // perfusion derivative    
//               dperfusiondu(u_theta)
//                               *bloodspecificheat*(u_theta-u_artery)
//              )     * (*phi[0])[i][qp] * (*phi[0])[j][qp] 
//        	                        );
//      }
//     for (unsigned int j=0; j<n_u_dofs[1]; j++)
//      {
//       Kij[0][1](i,j) += JxW[qp]*(
//           -mu_a[field_id] * m_theta * (*phi[1])[j][qp] * (*phi[0])[i][qp] 
//        	                        );
//      }
//    }
//   // Jzu and Jzz contributions
//   for (unsigned int i=0; i< n_u_dofs[1]; i++)
//    {
//     // no coupled term contribution
//     for (unsigned int j=0; j<n_u_dofs[1];j++)
//      {
//       Kij[1][1](i,j) += JxW[qp]*(
//               alpha[field_id]*m_theta  // diffusion term
//                    * ( (*dphi[1])[i][qp] * (*dphi[1])[j][qp] ) 
//                                     +
//              (1.0/user->get_fem_dt() //transient term
//                        +    // source term 
//               speedLight[field_id]*mu_a[field_id]*m_theta
//              )     * (*phi[1])[i][qp] * (*phi[1])[j][qp] 
//        	                        );
//      }
//    }
//  } // end of the quadrature point (qp-loop) over element interior
//
// PetscFunctionReturnVoid();
//}
/* -------------------------------------------------------------------- 
   Steady State coupled fluence solver
   -------------------------------------------------------------------- */ 

PennesFluenceSS::PennesFluenceSS(GetPot &controlfile, libMesh::MeshBase &mesh)
                :PennesFluencePDE(       controlfile,       mesh)
{ 
 /*
 if(FiniteElementInterface::NEXACT)
   {
    // set exact solution as IC
    std::vector<InitialConditionFncterType>  
                   fluence_IC(user->get_num_elem_blk() ,pdeExactSolution);
    InitValues.pop_back();
    InitValues.push_back(fluence_IC);
   }
 switch(FiniteElementInterface::NEXACT)
   {
    case 1:      
       //errorTol = 0.85; // get the stored tolerance for the verification problem
       break;
    default:
       //errorTol = 0.0;
   }
 */
}

Number PennesFluenceSS::exactSolution(const Point& p,
                                      const Parameters& parameters,
                                      const std::string& unknown_name)
{
 Number solutionValue ; 
 FiniteElementInterface *user = FiniteElementInterface::getInstance();
 PetscScalar time = user->getTime(FiniteElementInterface::ISTEP);
 Number L = 10.;
 Number x = p(0);
 switch( user->get_num_exact() )
   {
    case 1:      
      {
         /* 
           phi_xx + speedLight * mu_a / alpha phi = 0 

           general solution of the form
       
           phi = A sin( gamma x ) + B cos( gamma x )

                  where gamma == sqrt(speedLight * mu_a / alpha)

           BC phi(0) = phi(L) = C 
                              ==> B = C  
                                  A = C*(1 -cos( gamma * L ))/(sin( gamma * L ))
         */
         Number gamma = std::sqrt( speedLight[0] * mu_a[0] / alpha[0] );
         //Number Bcoeff = initialFluence[1];
         //Number Acoeff = initialFluence[1]*(1.0 - std::cos( gamma * L ))/
         //                                         std::sin( gamma * L )  ; 
         Number C1 = initialFluence[1];
         Number C2 = initialFluence[1];
         if ( x > 0 && x < L ) 
           {
               //Acoeff * std::sin( gamma * x ) + Bcoeff * std::cos( gamma * x );
            solutionValue = 
              (C2-C1*std::exp(-gamma*L))/(std::exp(gamma*L)-std::exp(-gamma*L))*std::exp(gamma*x)+(std::exp(gamma*L)*C1-C2)/(std::exp(gamma*L)-std::exp(-gamma*L))*std::exp(-gamma*x);

           }
         else
           { 
            solutionValue = initialFluence[1];
           }
      }
       break;
    case 2:     
       {
         const unsigned int field_id = 0 ; 
         //solutionValue = initialFluence.at(user->get_num_elem_blk()  - 1)  ;
         //if( p.size() > .003 ) 
         FiniteElementInterface *user = FiniteElementInterface::getInstance();
         PetscInt idpow= user->get_id_power();
         solutionValue = laserFluence(field_id ,p, idpow );
       }
       break;
    //default:
    //   //errorTol = 0.0;
   }


  return solutionValue;  	 	 
}
//void PennesFluenceSS::accumulateResidual(const QGauss &qrule, 
//                                          const unsigned int &field_id, 
//                                          const std::vector<Real>& JxW,
//                                          const std::vector<Point>& qpoint,
//                                          std::vector< DenseSubVector<Number> > &Ri,  
//                                          TransientFEMSystem &system)
//{
// PetscFunctionBegin;
// FiniteElementInterface *user = FiniteElementInterface::getInstance();
//
// for (unsigned int qp=0; qp<qrule.n_points(); qp++)
//  {
//   // evaluate Solution
//   this->evaluateSolutionForState(qp,system);
//
//   const Real u_theta = (1.0-m_theta)*u_previous[0]+m_theta*u__current[0];
//   const Real z_theta = (1.0-m_theta)*u_previous[1]+m_theta*u__current[1];
//   // accumulate the residual for temperature
//   for (unsigned int i=0; i<n_u_dofs[0]; i++)
//     {
//       Ri[0](i) += JxW[qp]*(
//         (rho*specific_heat*(u__current[0]-u_previous[0])/user->get_fem_dt() //transient term
//                   +    // perfusion term (check the SIGN)
//          perfusion(field_id,u_theta) *bloodspecificheat*(u_theta-u_artery)
//                   -    // source term 
//          mu_a[field_id] * z_theta 
//         )* (*phi[0])[i][qp]
//                   +    // diffusion term
//         thermal_cond(field_id,u_theta)*grad_u_mtheta[0]*(*dphi[0])[i][qp]
//                           );
//     }
//
//   // accumulate the residual for the fluence
//   for (unsigned int i=0; i<n_u_dofs[1]; i++)
//     {
//       Ri[1](i) += JxW[qp]*(
//                   // source term
//                  speedLight[field_id]* mu_a[field_id] * z_theta 
//                  * (*phi[1])[i][qp]
//                           +    // diffusion term
//                 alpha[field_id]*grad_u_mtheta[1]*(*dphi[1])[i][qp]
//                           );
//     }
//  } // end of the quadrature point (qp-loop) over element interior
//
// PetscFunctionReturnVoid();
//}
//
//void PennesFluenceSS::accumulateJacobian( const QGauss &qrule,
//            const unsigned int &field_id, 
//            const std::vector<Real>& JxW,
//            std::vector< SubMatrixVector > &Kij,
//            TransientFEMSystem &system)  
//{
// PetscFunctionBegin;
// FiniteElementInterface *user = FiniteElementInterface::getInstance();
//
// for (unsigned int qp=0; qp<qrule.n_points(); qp++)
//  {
//   // evaluate Solution
//   this->evaluateSolutionForState(qp,system);
//
//   const Real u_theta = (1.0-m_theta)*u_previous[0]+m_theta*u__current[0];
//   // Juu and Juz contributions
//   for (unsigned int i=0; i<n_u_dofs[0]; i++)
//    {
//     for (unsigned int j=0; j<n_u_dofs[0]; j++)
//      {
//       Kij[0][0](i,j) += JxW[qp]*(
//              thermal_cond(field_id,u_theta)*m_theta  // diffusion term
//                    * ( (*dphi[0])[i][qp] * (*dphi[0])[j][qp] ) 
//                                     +    // diffusion derivative (non-symmetric)
//              dthermal_conddu(u_theta)
//                    * ( (*dphi[0])[i][qp] * grad_u_mtheta[0] ) * 
//                                                   (*phi[0])[j][qp]
//                                     +
//              (rho*specific_heat/user->get_fem_dt() //transient term
//                        +    // perfusion term (check the SIGN)
//               perfusion(field_id,u_theta)
//                               *bloodspecificheat*m_theta
//                        +    // perfusion derivative    
//               dperfusiondu(u_theta)
//                               *bloodspecificheat*(u_theta-u_artery)
//              )     * (*phi[0])[i][qp] * (*phi[0])[j][qp] 
//        	                        );
//      }
//     for (unsigned int j=0; j<n_u_dofs[1]; j++)
//      {
//       Kij[0][1](i,j) += JxW[qp]*(
//           -mu_a[field_id] * m_theta * (*phi[1])[j][qp] * (*phi[0])[i][qp] 
//        	                        );
//      }
//    }
//   // Jzu and Jzz contributions
//   for (unsigned int i=0; i< n_u_dofs[1]; i++)
//    {
//     // no coupled term contribution
//     for (unsigned int j=0; j<n_u_dofs[1];j++)
//      {
//       Kij[1][1](i,j) += JxW[qp]*(
//               alpha[field_id]*m_theta  // diffusion term
//                    * ( (*dphi[1])[i][qp] * (*dphi[1])[j][qp] ) 
//                        +    // source term 
//               speedLight[field_id]*mu_a[field_id]*m_theta
//                    * (*phi[1])[i][qp] * (*phi[1])[j][qp] 
//        	                        );
//      }
//    }
//  } // end of the quadrature point (qp-loop) over element interior
//
// PetscFunctionReturnVoid();
//}

/* ------------------------------------------------------------------- 
     FIXME: Fortran Hack. mainly for fortran compatibility with old verification
     suite to pass data buffer back and forth
   
   ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "PennesInverseModel::fillSolutionBuffer"
void PennesInverseModel::fillSolutionBuffer( 
          std::vector<PetscScalar> &solnBuffer, const unsigned int idVar  )
{
 solnBuffer.clear( ); // reset
 // order is important in routines called by function pointers
 solnBuffer.push_back(  u__current[idVar]   ); // i_var*8+0
 solnBuffer.push_back(  u_previous[idVar]   ); // i_var*8+1
 solnBuffer.push_back(  u___future[idVar]   ); // i_var*8+2
 solnBuffer.push_back(  s__current[idVar]   ); // i_var*8+3
 solnBuffer.push_back(  s_previous[idVar]   ); // i_var*8+4
 solnBuffer.push_back(  s___future[idVar]   ); // i_var*8+5
 solnBuffer.push_back(  s_zero[    idVar]   ); // i_var*8+6
 solnBuffer.push_back(  s_full[    idVar]   ); // i_var*8+7
 return ;
}
