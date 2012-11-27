// rbOOmit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#ifndef __petsc_fem_context_h__
#define __petsc_fem_context_h__

// Local Includes
#include "fem_context.h"

namespace libMesh
{

/**
 * provides extra context
 * information relevant to Transient problems
 *
 * FEMContext should allow to cache data locally 
 *   an example would be taking
 *   the full soln vector (O(10000) entries) and
 *   storing the entries needed in a
 *   local vector (8 entries for hex) during element 
 *   operations to keep everything close to the core
 *   during the soln evaluation. During recomputation of
 *   solution at quadrature points during jacobian/residual
 *   everything should be in L1-L2 cache to speed up computations
 *   by minimizing memory operations
 * 
 * @author David Fuentes, 2011
 */

// ------------------------------------------------------------
// PetscFEMContext class definition

class PetscFEMContext : public FEMContext
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  PetscFEMContext (const System &sys);

  /**
   * Destructor.
   */
  virtual ~PetscFEMContext ();

  /**
   * Returns the value of the solution variable \p var at the quadrature
   * point \p qp on the current element interior at theta time point
   */
  Number interior_theta_value(unsigned int var, unsigned int qp);

  /**
   * Returns the value of the solution variable \p var at the quadrature
   * point \p qp on the current element side at theta time point
   */
  Number side_theta_value(unsigned int var, unsigned int qp);

  /**
   * Returns the value of the diffsolution variable \p var at the quadrature
   * point \p qp on the current element interior at theta time point
   */
  Number interior_diff_value(unsigned int var, unsigned int qp);

  /**
   * Returns the gradient of the solution variable \p var at the quadrature
   * point \p qp on the current element interior at theta time point
   */
  Gradient interior_theta_gradient(unsigned int var, unsigned int qp);

  /**
   * Returns the gradient of the solution variable \p var at the quadrature
   * point \p qp on the current element side at theta time point
   */
  Gradient side_theta_gradient(unsigned int var, unsigned int qp);

  /**
   * Reinitialize all the context data on a given
   * element for the given system. Overload to
   * also reinitialize the solution data from
   * the previous time step.
   * cache the needed data locally
   */
  virtual void pre_fe_reinit(const System&, Elem*);

  void SetThetaValue(unsigned int idVar,Number theta)  
       { m_theta.at(idVar)= theta; return ; }

  Number ThetaValue(unsigned int idVar)  const
       { return m_theta[idVar] ; }

  /**
   * Element by element components of nonlinear_solution
   * as adjusted by a time_solver
   */
  DenseVector<Number> elem_old_solution;
  std::vector<DenseSubVector<Number> *> elem_old_subsolutions;

  //Real  interior_mtheta_value(unsigned int idVar) 
  //       {return (1.0-m_theta)* u_previous[idVar] 
  //                   +m_theta * u__current[idVar]; }
  //Gradient interior_mtheta_gradient(unsigned int idVar)
  //       {return grad_u_mtheta[idVar] ; }
  //Real  interior_diff_value(unsigned int idVar) 
  //       {return u__current[idVar] - u_previous[idVar]; }
  //Real  ThetaValue() {return m_theta; }
  // /**
  //  * @deprecated {evaluate the solution at the gauss points using
  //  *              pre-femsystem data structures}
  //  * @todo { use FEMContext to cache data locally 
  //  *           an example would be taking
  //  *           the full soln vector (O(10000) entries) and
  //  *           storing the entries needed in a
  //  *           local vector (8 entries for hex) during element 
  //  *           operations to keep everything close to the core
  //  *           during the soln evaluation.
  //  *        }
  //  */
  // void evaluateSolutionForState(const unsigned int &,
  //                           TransientFEMSystem &);
  // /**
  //  * @deprecated {evaluate the solution at the gauss points using
  //  *              pre-femsystem data structures}
  //  * @todo { use FEMContext to cache data locally} 
  //  */
  // void evaluateSolutionForAdjoint(const unsigned int &, const int,
  //                           TransientFEMSystem &,
  //                           TransientFEMSystem  &);
  // /**
  //  * @deprecated {evaluate the solution at the gauss points using
  //  *              pre-femsystem data structures}
  //  * @todo { use FEMContext to cache data locally} 
  //  */
  // void evaluateSolutionForAdjoint(const unsigned int &, 
  //                           TransientFEMSystem  &);
  // /**
  //  * @deprecated {evaluate the solution at the gauss points using
  //  *              pre-femsystem data structures}
  //  * @todo { use FEMContext to cache data locally} 
  //  */
  // void evaluateIdealSolution(const unsigned int &, const int, const int,
  //                        TransientFEMSystem &,
  //                        TransientFEMSystem          &, 
  //                        TransientFEMSystem          &); 
  // /**
  //  * @deprecated {evaluate the solution at the gauss points using
  //  *              pre-femsystem data structures}
  //  * @todo { use FEMContext to cache data locally} 
  //  */
  // void evaluateSolutionForGradient(const unsigned int &, const int, 
  //                           TransientFEMSystem &,
  //                           TransientFEMSystem  &);
  // /**
  //  * @deprecated {evaluate the solution at the gauss points using
  //  *              pre-femsystem data structures}
  //  * @todo { use FEMContext to cache data locally} 
  //  */
  // void evaluateSolutionForSensitivity(const unsigned int &,
  //                   TransientFEMSystem  &);
  // /**
  //  * @deprecated {evaluate the solution at the gauss points using
  //  *              pre-femsystem data structures}
  //  * @todo { use FEMContext to cache data locally} 
  //  */
  // void evaluateSolutionForSensitivity(const unsigned int &, const int,
  //                                 TransientFEMSystem    &);
  // /** adjoint variable used in dpde_dm functer and used as storage space to
  //  * pass in the shape function values 
  //  * @todo { altering the interface to
  //  * pass adjoint variable in gradient accumulation and shape functions in
  //  * sensitivity solve will be a lot of work }
  //  */
  // void evaluateAdjointSolutionAsShapeFunction(const unsigned int &i, 
  //                                             const unsigned int &qp)
  //   {
  //     for( unsigned int i_var = 0 ; i_var < numSystemVars ; i_var++)
  //      {
  //        p__now[i_var]      =  (*phi[i_var])[i][qp];
  //        grad_p__now[i_var] = (*dphi[i_var])[i][qp];
  //      }
  //     return;
  //   }
  // /** 
  //  * sensitivity variable used in d2pde_du_dm functer and used as storage space
  //  * to pass in the shape function values 
  //  * @todo { altering the interface to pass
  //  * shape functions in adjoint sensitivity load will be a lot of work}
  //  */ 
  // void evaluateSensitivitySolutionAsShapeFunction(const unsigned int &i, 
  //                                                   const unsigned int &qp)
  //     {
  //       for( unsigned int i_var = 0 ; i_var < numSystemVars ; i_var++)
  //        {
  //          du_previous[   i_var] = (*phi[ i_var])[i][qp];
  //          du__current[   i_var] = (*phi[ i_var])[i][qp];
  //          grad_du_mtheta[i_var] = (*dphi[i_var])[i][qp];
  //        }
  //       return;
  //     }

  // /**
  //  * difference of solution at guass points
  //  * @deprecated {evaluate the solution at the gauss points using
  //  *              pre-femsystem data structures}
  //  * @todo { use FEMContext to cache data locally} 
  //  */
  // PetscScalar getWeightedDifference(const unsigned int &i)
  //      {
  //        return  (u__current[i]-s__current[i])/ds_current[i]/ds_current[i];
  //      }

  // /**
  //  * weighted sensitivity
  //  * @deprecated {evaluate the solution at the gauss points using
  //  *              pre-femsystem data structures}
  //  * @todo { use FEMContext to cache data locally} 
  //  */
  // PetscScalar getWeightedSensitivity(const unsigned int &i)
  //      {
  //        return  du__current[i]/ds_current[i]/ds_current[i];
  //      }
 
  // // variable dofs for state/ideal solution
  // std::vector< unsigned int > n_u_dofs;
  // std::vector< std::vector< unsigned int > >  dof_indices_u;

  // // variable dofs for adjoint solution (possibly different order)
  // std::vector< unsigned int > n_p_dofs;
  // std::vector< std::vector< unsigned int > >  dof_indices_p;

  // // The element shape functions and function gradients
  // // evaluated at the quadrature points for the state/ideal solution
  // std::vector< const std::vector<std::vector<Real> >* > phi;
  // std::vector< const std::vector<std::vector<RealGradient> >* > dphi;

  // // The element shape functions and function gradients
  // // evaluated at the quadrature points for the adjoint solution 
  // // (possibly different order)
  // std::vector< const std::vector<std::vector<Real> >* > psi;
  // std::vector< const std::vector<std::vector<RealGradient> >* > dpsi;
  // /**
  //  * Values to hold the solution & its gradient at the previous 
  //  * timestep and the current "position". 
  //  * Declared static to allow access to the solution at the 
  //  * quadrature points by the QOI classes
  //  * @deprecated {evaluate the solution at the gauss points using
  //  *              pre-femsystem data structures}
  //  * @todo { use FEMContext to cache data locally} 
  //  */
  // std::vector<Real> u__current,u_previous,u___future;
  // std::vector<RealGradient> grad_u_mtheta,grad_u_ptheta;
  // std::vector<Real>              p___future,p__now;
  // std::vector<RealGradient> grad_p___future,grad_p__now;
  // std::vector<Real>  s_previous, s__current, s___future, s_zero, s_full;
  // std::vector<Real>  ds_current; // for uncertainty

  // // sensitivity
  // std::vector< Real         >  du__current,du_previous;
  // std::vector< RealGradient >  grad_du_mtheta,grad_du_previous;

  // unsigned int numSystemVars ; 

private:
  /**
   * Theta values for each variable
   * time stepping control 
   *            m_theta = 0.0 ==> backward euler
   *            m_theta = 0.5 ==> crank nicolson
   *            m_theta = 1.0 ==> forward  euler
   */
  std::vector<Number> m_theta;
};

} // namespace libMesh

#endif
