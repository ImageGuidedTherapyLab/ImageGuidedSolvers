// distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include "petsc_fem_context.h"
#include "dof_map.h"
#include "elem.h"
#include "fe_base.h"
#include "fe_interface.h"

//local
#include "optimizationParameter.h" 
#include "petsc_fem_system.h"

namespace libMesh
{

PetscFEMContext::PetscFEMContext (const System& sys) :
  FEMContext(sys)
{
  unsigned int n_vars = sys.n_vars();
  elem_old_subsolutions.reserve(n_vars);
  m_theta.resize(n_vars,0.5);// assume crank nicolson

  for (unsigned int i=0; i != n_vars; ++i)
  {
    elem_old_subsolutions.push_back(new DenseSubVector<Number>(elem_old_solution));
  }
}

PetscFEMContext::~PetscFEMContext ()
{
  for (unsigned int i=0; i != elem_old_subsolutions.size(); ++i)
  {
    delete elem_old_subsolutions[i];
  }
}

Number PetscFEMContext::interior_theta_value(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_old_subsolutions.size() > var);
  libmesh_assert (elem_old_subsolutions[var] != NULL);
  libmesh_assert (elem_subsolutions.size() > var);
  libmesh_assert (elem_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef     = *elem_subsolutions[var];
  DenseSubVector<Number> &coef_old = *elem_old_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<Real> > &phi =
    element_fe_var[var]->get_phi();

  // Accumulate solution value
  libmesh_assert (m_theta.size() > var);
  Number u_theta = 0.;

  for (unsigned int l=0; l != n_dofs; l++)
    u_theta += phi[l][qp] * (
                   m_theta[var]   * coef(l)
                            +
          ( 1.0 -  m_theta[var] ) * coef_old(l)
                            );
  return u_theta;
}

Gradient PetscFEMContext::interior_theta_gradient(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_old_subsolutions.size() > var);
  libmesh_assert (elem_old_subsolutions[var] != NULL);
  libmesh_assert (elem_subsolutions.size() > var);
  libmesh_assert (elem_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef     = *elem_subsolutions[var];
  DenseSubVector<Number> &coef_old = *elem_old_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<RealGradient> > &dphi =
    element_fe_var[var]->get_dphi();

  // Accumulate solution derivatives
  libmesh_assert (m_theta.size() > var);
  Gradient theta_du;

  for (unsigned int l=0; l != n_dofs; l++)
    theta_du.add_scaled(dphi[l][qp], 
                   m_theta[var]   * coef(l)
                            +
          ( 1.0 -  m_theta[var] ) * coef_old(l)
                       );

  return theta_du;
}

Number PetscFEMContext::side_theta_value(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_old_subsolutions.size() > var);
  libmesh_assert (elem_old_subsolutions[var] != NULL);
  libmesh_assert (elem_subsolutions.size() > var);
  libmesh_assert (elem_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef     = *elem_subsolutions[var];
  DenseSubVector<Number> &coef_old = *elem_old_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<Real> > &phi =
    side_fe_var[var]->get_phi();

  // Accumulate solution value
  libmesh_assert (m_theta.size() > var);
  Number u_theta = 0.;

  for (unsigned int l=0; l != n_dofs; l++)
    u_theta += phi[l][qp] *   (
                   m_theta[var]   * coef(l)
                            +
          ( 1.0 -  m_theta[var] ) * coef_old(l)
                            );

  return u_theta;
}

Gradient PetscFEMContext::side_theta_gradient(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_old_subsolutions.size() > var);
  libmesh_assert (elem_old_subsolutions[var] != NULL);
  libmesh_assert (elem_subsolutions.size() > var);
  libmesh_assert (elem_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef     = *elem_subsolutions[var];
  DenseSubVector<Number> &coef_old = *elem_old_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<RealGradient> > &dphi =
    side_fe_var[var]->get_dphi();

  // Accumulate solution derivatives
  libmesh_assert (m_theta.size() > var);
  Gradient du_theta;

  for (unsigned int l=0; l != n_dofs; l++)
    du_theta.add_scaled(dphi[l][qp], 
                   m_theta[var]   * coef(l)
                            +
          ( 1.0 -  m_theta[var] ) * coef_old(l)
                       );

  return du_theta;
}

Number PetscFEMContext::interior_diff_value(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_old_subsolutions.size() > var);
  libmesh_assert (elem_old_subsolutions[var] != NULL);
  libmesh_assert (elem_subsolutions.size() > var);
  libmesh_assert (elem_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef     = *elem_subsolutions[var];
  DenseSubVector<Number> &coef_old = *elem_old_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<Real> > &phi =
    element_fe_var[var]->get_phi();

  // Accumulate solution value
  libmesh_assert (m_theta.size() > var);
  Number u_diff = 0.;

  for (unsigned int l=0; l != n_dofs; l++)
    u_diff += phi[l][qp] * ( coef(l) - coef_old(l));
  return u_diff;
}

void PetscFEMContext::pre_fe_reinit(const System &sys, Elem *e)
{
  FEMContext::pre_fe_reinit(sys, e);
        
  // sys must be a QNTransientRBSystem
  const PetscFEMSystem& tfem_system = libmesh_cast_ref<const PetscFEMSystem&>(sys);

  // Initialize the per-element data for elem.
  unsigned int n_dofs = dof_indices.size();
  elem_old_solution.resize(n_dofs);

  for (unsigned int i=0; i != n_dofs; ++i)
    elem_old_solution(i) = tfem_system.prev_solution(dof_indices[i]);

  // Initialize the per-variable data for elem.
  unsigned int sub_dofs = 0;
  for (unsigned int i=0; i != tfem_system.n_vars(); ++i)
    {
      tfem_system.get_dof_map().dof_indices (elem, dof_indices_var[i], i);

      elem_old_subsolutions[i]->reposition
        (sub_dofs, dof_indices_var[i].size());

      sub_dofs += dof_indices_var[i].size();
    }
  libmesh_assert(sub_dofs == n_dofs);
}

///*----------------------------------------------------------------------*/
//void PDEModelBaseClass::SetupState(EquationSystems &)
//{
// PetscFunctionBegin;
//
// // Values to hold the solution & its gradient at the previous 
// // timestep and the current "position".
// u__current.resize(numSystemVars,0.0);
// u_previous.resize(numSystemVars,0.0);
// grad_u_mtheta.resize(numSystemVars,0.0);
//
// // variable dofs
// n_u_dofs.resize(numSystemVars,0); 
// dof_indices_u.resize(numSystemVars);
//
// // The element shape functions evaluated at the quadrature points.
// phi.resize(numSystemVars,NULL);
//
// // The element shape function gradients evaluated at the quadrature points.
// dphi.resize(numSystemVars,NULL);
//
// // error check
// libmesh_assert ( m_NeumannFlux.size() == numSystemVars       );
// libmesh_assert ( m_newton_coeff.size() == numSystemVars );
// libmesh_assert ( m_u_infty.size() == numSystemVars      );
//
// PetscFunctionReturnVoid();
//}
//
///*----------------------------------------------------------------------*/
//void PDEModelBaseClass::SetupAdjoint( const unsigned int )
//{
// PetscFunctionBegin;
//
// // Values to hold the solution & its gradient at the previous 
// // timestep and the current "position".
//         u___future.resize(numSystemVars,0.0);
//         p___future.resize(numSystemVars,0.0);
//             p__now.resize(numSystemVars,0.0);
//         s_previous.resize(numSystemVars,0.0);
//         s__current.resize(numSystemVars,0.0);
//         s___future.resize(numSystemVars,0.0);
//             s_zero.resize(numSystemVars,0.0);
//             s_full.resize(numSystemVars,0.0);
//         ds_current.resize(numSystemVars,0.0);
//      grad_u_ptheta.resize(numSystemVars,0.0);
//    grad_p___future.resize(numSystemVars,0.0);
//        grad_p__now.resize(numSystemVars,0.0);
//        du_previous.resize(numSystemVars,0.0);
//        du__current.resize(numSystemVars,0.0);
//     grad_du_mtheta.resize(numSystemVars,0.0);
//   grad_du_previous.resize(numSystemVars,0.0);
//
// // variable dofs
// n_p_dofs.resize(numSystemVars,0); 
// dof_indices_p.resize(numSystemVars);
//
// // The element shape functions evaluated at the quadrature points.
// psi.resize(numSystemVars,NULL);
//
// // The element shape function gradients evaluated at the quadrature points.
// dpsi.resize(numSystemVars,NULL);
//
// CHKMEMA; // check for memory corruption use -malloc_debug to enable
//  
// PetscFunctionReturnVoid();
//}

// /*----------------------------------------------------------------------*/
// void PDEModelBaseClass::evaluateSolutionForState(const unsigned int &qp, 
//                        PetscFEMSystem &system)
// {
//   PetscFunctionBegin;
//  
//   // notice that the current time gradient and solution
//   // are computed with the "solution" vector 
//   // "vector_solution" has not been updated yet
//   for( unsigned int i_var = 0 ; i_var < numSystemVars ; i_var++)
//    {
//      // re initialize
//      u__current[i_var]=0.0;
//      u_previous[i_var]=0.0;
//      grad_u_mtheta[i_var].zero() ; 
//      for (unsigned int i=0; i<n_u_dofs[i_var]; i++)
//       {
//         u__current[i_var]+=(*phi[i_var])[i][qp]*system.current_solution(
//                                                 dof_indices_u[i_var][i]);
//         u_previous[i_var]+=(*phi[i_var])[i][qp]*
//                                  system.prev_solution(
//                                                 dof_indices_u[i_var][i]);
//         grad_u_mtheta[i_var].add_scaled ((*dphi[i_var])[i][qp],
//             (1.0-m_theta)*system.prev_solution(   dof_indices_u[i_var][i])
//                                  +
//                  m_theta *system.current_solution(dof_indices_u[i_var][i])
//                                         );
//       }
//    }                                               
//   PetscFunctionReturnVoid();
// }
// /*----------------------------------------------------------------------*/
// void PDEModelBaseClass::evaluateSolutionForAdjoint(const unsigned int &qp,
//                                               const int istep, 
//                        PetscFEMSystem &state_system,
//                        PetscFEMSystem  &adjoint_system)
// {
//  PetscFunctionBegin; 
// 
//  // the state solution computed during the state solve is used
//  for( unsigned int i_var = 0 ; i_var < numSystemVars; i_var++)
//   {
//     // re initialize
//     u__current[i_var]=0.0;
//     u_previous[i_var]=0.0;
//     u___future[i_var]=0.0;
//     grad_u_mtheta[i_var].zero() ; 
//     grad_u_ptheta[i_var].zero() ; 
//     for (unsigned int i=0; i<n_u_dofs[i_var]; i++)
//      {
//        //state
//        u__current[i_var]+=(*phi[i_var])[i][qp]* state_system.stored_solution( 
//                                              dof_indices_u[i_var][i],istep  );
//        u_previous[i_var]+=(*phi[i_var])[i][qp]* state_system.stored_solution( 
//                                              dof_indices_u[i_var][i],istep-1);
//        u___future[i_var]+=(*phi[i_var])[i][qp]* state_system.stored_solution( 
//                                              dof_indices_u[i_var][i],istep+1);
//        grad_u_mtheta[i_var].add_scaled ((*dphi[i_var])[i][qp],
//             (1.0- m_theta)*state_system.stored_solution(
//                                        dof_indices_u[i_var][i],istep-1)
//                                   +
//                   m_theta *state_system.stored_solution(
//                                        dof_indices_u[i_var][i],istep  )
//                                        );
//        grad_u_ptheta[i_var].add_scaled ((*dphi[i_var])[i][qp],
//                   m_theta *state_system.stored_solution(
//                                        dof_indices_u[i_var][i],istep+1)
//                                   +
//             (1.0- m_theta)*state_system.stored_solution(
//                                        dof_indices_u[i_var][i],istep  )
//                                        );
//      }
//   }
// 
//  // variables to hold the adjoint and  adjoint gradient
//  for( unsigned int i_var = 0 ; i_var < numSystemVars; i_var++)
//   {
//     // re initialize
//     p___future[i_var]=0.0;
//     grad_p___future[i_var].zero() ; 
//     for (unsigned int i=0; i<n_p_dofs[i_var]; i++)
//      {
//       p___future[i_var]+=(*psi[i_var])[i][qp]*
//         adjoint_system.prev_solution(dof_indices_p[i_var][i])  ;
//       grad_p___future[i_var].add_scaled ( (*dpsi[i_var])[i][qp],
//         adjoint_system.prev_solution(dof_indices_p[i_var][i]) );
//      }
//   }
//  PetscFunctionReturnVoid(); 
// }
// /*----------------------------------------------------------------------*/
// void PDEModelBaseClass::evaluateSolutionForAdjoint(const unsigned int &qp,
//                        PetscFEMSystem  &adjoint_system)
// {
//  PetscFunctionBegin; 
// 
//  // variables to hold the adjoint and  adjoint gradient
//  for( unsigned int i_var = 0 ; i_var < numSystemVars; i_var++)
//   {
//     // re initialize
//     p__now[i_var]=0.0;
//     grad_p__now[i_var].zero() ; 
//     for (unsigned int i=0; i<n_p_dofs[i_var]; i++)
//      {
//       p__now[i_var]+=(*psi[i_var])[i][qp]*
//            adjoint_system.current_solution( dof_indices_p[i_var][i]  );
//       grad_p__now[i_var].add_scaled ( (*dpsi[i_var])[i][qp],
//            adjoint_system.current_solution(dof_indices_p[i_var][i] ) );
//      }
//   }
//  PetscFunctionReturnVoid(); 
// }
// /*----------------------------------------------------------------------*/
// void PDEModelBaseClass::evaluateIdealSolution(const unsigned int &qp,
//                                          const int istep, const int Nstephi,
//                      PetscFEMSystem    &state_system,
//                      PetscFEMSystem             &ideal_system, 
//                      PetscFEMSystem &ideal_uncertainty_system) 
// {
//  PetscFunctionBegin; 
//  // the state solution computed during the state solve is used
//  for( unsigned int i_var = 0 ; i_var < numSystemVars; i_var++)
//   {
//     // re initialize
//     u__current[i_var]=0.0;
//     u_previous[i_var]=0.0;
//     s_previous[i_var]=0.0;
//     s__current[i_var]=0.0;
//     s___future[i_var]=0.0;
//     s_zero[i_var]=0.0;
//     s_full[i_var]=0.0;
//     ds_current[i_var]=0.0;
//     for (unsigned int i=0; i<n_u_dofs[i_var]; i++)
//      {
//        //state
//        u__current[i_var]+=(*phi[i_var])[i][qp]* state_system.stored_solution( 
//                                              dof_indices_u[i_var][i],istep  );
//        u_previous[i_var]+=(*phi[i_var])[i][qp]* state_system.stored_solution( 
//                                              dof_indices_u[i_var][i],istep-1);
//        // ideal
//        s_zero[i_var]     += (*phi[i_var])[i][qp]*ideal_system.stored_solution(
//                                           dof_indices_u[i_var][i],0           );
//        s_previous[i_var] += (*phi[i_var])[i][qp]*ideal_system.stored_solution(
//                                           dof_indices_u[i_var][i],istep-1     );
//        s__current[i_var] += (*phi[i_var])[i][qp]*ideal_system.stored_solution(
//                                           dof_indices_u[i_var][i],istep       );
//        s___future[i_var] += (*phi[i_var])[i][qp]*ideal_system.stored_solution(
//                                           dof_indices_u[i_var][i],istep+1     );
//        s_full[i_var]     += (*phi[i_var])[i][qp]*ideal_system.stored_solution(
//                                           dof_indices_u[i_var][i],Nstephi);
//        // ideal uncertainty
//        ds_current[i_var] += (*phi[i_var])[i][qp]*
//                                      ideal_uncertainty_system.stored_solution(
//                                           dof_indices_u[i_var][i],istep       );
//      }
//   }
//  PetscFunctionReturnVoid(); 
// }
// /*----------------------------------------------------------------------*/
// void PDEModelBaseClass::evaluateSolutionForGradient(const unsigned int &qp,
//                                                const int istep, 
//                        PetscFEMSystem &state_system,
//                        PetscFEMSystem  &adjoint_system)
// {
//  PetscFunctionBegin; 
// 
//  // the state solution computed during the state solve is used
//  for( unsigned int i_var = 0 ; i_var < numSystemVars; i_var++)
//   {
//     // re initialize
//     u__current[i_var]=0.0;
//     u_previous[i_var]=0.0;
//     grad_u_mtheta[i_var].zero() ; 
//     for (unsigned int i=0; i<n_u_dofs[i_var]; i++)
//      {
//        //state
//        u__current[i_var]+=(*phi[i_var])[i][qp]* state_system.stored_solution( 
//                                              dof_indices_u[i_var][i],istep  );
//        u_previous[i_var]+=(*phi[i_var])[i][qp]* state_system.stored_solution( 
//                                              dof_indices_u[i_var][i],istep-1);
//        grad_u_mtheta[i_var].add_scaled ((*dphi[i_var])[i][qp],
//             (1.0- m_theta)*state_system.stored_solution(
//                                        dof_indices_u[i_var][i],istep-1)
//                                   +
//                   m_theta *state_system.stored_solution(
//                                        dof_indices_u[i_var][i],istep  )
//                                        );
//      }
//   }
// 
//  // variables to hold the adjoint and  adjoint gradient
//  for( unsigned int i_var = 0 ; i_var < numSystemVars; i_var++)
//   {
//     // re initialize
//     p__now[i_var]=0.0;
//     grad_p__now[i_var].zero() ; 
//     for (unsigned int i=0; i<n_p_dofs[i_var]; i++)
//      {
//       p__now[i_var]+=(*psi[i_var])[i][qp]*
//            adjoint_system.stored_solution(dof_indices_p[i_var][i],istep  );
//       grad_p__now[i_var].add_scaled ( (*dpsi[i_var])[i][qp],
//            adjoint_system.stored_solution(dof_indices_p[i_var][i],istep  ));
//      }
//   }
//  PetscFunctionReturnVoid(); 
// }
// /*----------------------------------------------------------------------*/
// void PDEModelBaseClass::evaluateSolutionForSensitivity(const unsigned int &qp,
//                    PetscFEMSystem  &sensitivity_system)
// {
//  PetscFunctionBegin; 
// 
//  // the state solution computed during the state solve is used
//  for( unsigned int i_var = 0 ; i_var < numSystemVars; i_var++)
//   {
//     // re initialize
//     du_previous[  i_var] = 0.0;
//     grad_du_previous[i_var].zero() ; 
//     for (unsigned int i=0; i<n_u_dofs[i_var]; i++)
//      { //sensitivities
//        du_previous[i_var] += (*phi[i_var])[i][qp]*
//                    sensitivity_system.prev_solution(dof_indices_u[i_var][i])
//                                              ;
//        grad_du_previous[i_var].add_scaled ( (*dphi[i_var])[i][qp],
//                    sensitivity_system.prev_solution(dof_indices_u[i_var][i])
//                                           );
//      }
//   }
// 
//  PetscFunctionReturnVoid(); 
// }
// /*----------------------------------------------------------------------*/
// void PDEModelBaseClass::evaluateSolutionForSensitivity(const unsigned int &qp,
//                                                   const int istep, 
//                    PetscFEMSystem  &sensitivity_system) 
// {
//  PetscFunctionBegin; 
// 
//  // the state solution computed during the state solve is used
//  for( unsigned int i_var = 0 ; i_var < numSystemVars; i_var++)
//   {
//     du__current[i_var] = 0.0;
//     du_previous[i_var] = 0.0;
//     grad_du_mtheta[i_var].zero() ; 
//     for (unsigned int i=0; i<n_u_dofs[i_var]; i++)
//      { //sensitivities
//        du_previous[i_var] += (*phi[i_var])[i][qp]*(
//          sensitivity_system.stored_solution( dof_indices_u[i_var][i],istep-1 ));
//        du__current[i_var] += (*phi[i_var])[i][qp]*(
//          sensitivity_system.stored_solution( dof_indices_u[i_var][i],istep   ));
//        grad_du_mtheta[i_var].add_scaled ( (*dphi[i_var])[i][qp],
//            (1.0- m_theta) * sensitivity_system.stored_solution( 
//                                              dof_indices_u[i_var][i],istep-1 )
//                                      +
//                  m_theta  * sensitivity_system.stored_solution(
//                                              dof_indices_u[i_var][i],istep   )
//                                             );
//      }
//   }
// 
//  PetscFunctionReturnVoid(); 
// }



} // namespace libMesh

