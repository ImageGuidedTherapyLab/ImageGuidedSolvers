#ifndef __thermal_therapy_system_h__
#define __thermal_therapy_system_h__

// libmesh includes
#include "petsc_fem_system.h"

/**
 * Common Routines for thermal therapy simulations.
 * additional vectors/matrices may be added,
 * as offered in the parent classes.
 */
template< typename MathematicalModel  >
class ThermalTherapySystem : public PetscFEMSystem 
{
public:

  /**
   * Constructor.  Initializes required
   * data structures.
   */
  ThermalTherapySystem (EquationSystems& es,
                          const std::string& name,
                          const unsigned int number ) ;

  /**
   * The type of the parent.
   */
  typedef PetscFEMSystem Parent;
  
  /**
   * setup intial conditions 
   */
  virtual void SetupInitialConditions ();
  
  /**
   * scatter global parameters to all locally
   */
  virtual void ScatterParametersLocally()
    { this->m_MathModel.ScatterParametersLocally(); }

  /**
   * Builds a FEMContext object with enough information to do
   * evaluations on each element.
   *
   * For most problems, the default FEMSystem implementation is correct; users
   * who subclass FEMContext will need to also reimplement this method to build
   * it.
   */
  virtual AutoPtr<DiffContext> build_context();

  /**
   * BC
   */
  virtual bool side_time_derivative(bool request_jacobian,
                                    DiffContext& context);

  /**
   * setup dirichlet data
   */
  virtual void SetupDirichlet(libMesh::MeshBase& );

  /** 
    * Print constitutive data parameters 
    * This is a wrapper to the templated constitutive data
    */
  virtual void printSelf(std::ostream& os)
    { 
      this->m_MathModel.printSelf(os); 
      return; 
    }

  /** 
    * return exact solution 
    */
  PetscScalar exactSolution(const Point& p,
                            const Parameters& parameters,
                            const std::string& unknown_name)
    {  return this->m_MathModel.exactSolution(p,parameters,unknown_name); }

  /** setup optimization variables */
  virtual void SetupOptimizationVariables( std::vector<optimizationParameter*> &Parameters)
   { this->m_MathModel.SetupOptimizationVariables(this->get_mesh(),Parameters);return; }
  
  /**
   *  setup Adjoint
   */
  virtual void SetupAdjoint( const unsigned int nParam)
   { libmesh_not_implemented();
     //this->m_MathModel.SetupAdjoint(nParam);
     return; }

 
  // /** wrapper should be inlined @ref PDEModelInlineStrategy */
  // virtual void evaluateAdjointSolutionAsShapeFunction(const unsigned int &i,
  //                                                        const unsigned int &qp)
  //  { 
  //    this->m_MathModel.evaluateAdjointSolutionAsShapeFunction(i,qp); 
  //    return;
  //  }
  // /** wrapper should be inlined @ref PDEModelInlineStrategy */
  // virtual void evaluateSensitivitySolutionAsShapeFunction(const unsigned int &i,
  //                                                        const unsigned int &qp)
  //  { 
  //    this->m_MathModel.evaluateSensitivitySolutionAsShapeFunction(i,qp); 
  //    return;
  //  }
  // /** wrapper should be inlined @ref PDEModelInlineStrategy */
  // virtual void evaluateSolutionForState(const unsigned int &qp, 
  //                                       PetscFEMSystem &system)
  //  { this->m_MathModel.evaluateSolutionForState(qp,system); return; }

  // /** wrapper should be inlined @ref PDEModelInlineStrategy */
  // virtual void evaluateSolutionForGradient(const unsigned int &qp,
  //                                              const int istep, 
  //                      PetscFEMSystem &state_system,
  //                      PetscFEMSystem  &adjoint_system)
  //  { this->m_MathModel.evaluateSolutionForGradient(qp,istep, 
  //                                     state_system,adjoint_system); return; }
  // /** wrapper should be inlined @ref PDEModelInlineStrategy */
  // virtual void evaluateSolutionForAdjoint(const unsigned int &qp,
  //                                         PetscFEMSystem  &adjoint_system)
  //  { this->m_MathModel.evaluateSolutionForAdjoint(qp,adjoint_system);return;}

  // /** wrapper should be inlined @ref PDEModelInlineStrategy */
  // virtual void evaluateIdealSolution(const unsigned int &qp,
  //                                    const int istep, const int Nstephi,
  //                    PetscFEMSystem    &state_system,
  //                    PetscFEMSystem             &ideal_system, 
  //                    PetscFEMSystem &ideal_uncertainty_system) 
  //  { this->m_MathModel.evaluateIdealSolution( qp, istep, Nstephi, 
  //               state_system, ideal_system, ideal_uncertainty_system);return;}
  // /** wrapper should be inlined @ref PDEModelInlineStrategy */
  // virtual void evaluateSolutionForSensitivity(const unsigned int &qp,
  //                                             const int istep, 
  //                                    PetscFEMSystem  &sensitivity_system) 
  //  { this->m_MathModel.evaluateSolutionForSensitivity(qp,
  //                                           istep,sensitivity_system);return;}

  /** wrapper should be inlined @ref PDEModelInlineStrategy */
  virtual void fillSolutionBuffer(
          std::vector<PetscScalar> &solnBuffer, const unsigned int idVar  )
   { this->m_MathModel.fillSolutionBuffer(solnBuffer,idVar);return;}

  /** setup optimization variables */
  virtual PetscErrorCode GetDataFromDisk(libMesh::MeshBase &mesh,std::string &FieldFile)
   { return this->m_MathModel.GetDataFromDisk(mesh,FieldFile); }

  /** wrapper should be inlined @ref PDEModelInlineStrategy */
  virtual PetscScalar dpde_dm(optimizationParameter* optParam ,
                             const unsigned int &field_id, const Point &q_point)
   { return CALL_MEMBER_FN_W_REF(this->m_MathModel,optParam->dpde_dm)
                                                    (field_id,q_point); }
  /** wrapper should be inlined @ref PDEModelInlineStrategy */
  virtual PetscScalar d2qoi_du_dm(optimizationParameter* optParam ,
                             const unsigned int &field_id, const Point &q_point)
   { return CALL_MEMBER_FN_W_REF(this->m_MathModel,optParam->d2qoi_du_dm)
                                                    (field_id,q_point); }
  /** wrapper should be inlined @ref PDEModelInlineStrategy */
  virtual PetscScalar d2pde_du_dm(optimizationParameter* optParam ,
                             const unsigned int &field_id, const Point &q_point)
   { return CALL_MEMBER_FN_W_REF(this->m_MathModel,optParam->d2pde_du_dm)
                                                    (field_id,q_point); }
  /** wrapper should be inlined @ref PDEModelInlineStrategy */
  virtual PetscScalar dqoidu_dudm( const unsigned int &field_id )
   { return this->m_MathModel.dqoidu_dudm(field_id); }

  /** wrapper should be inlined @ref PDEModelInlineStrategy */
  virtual PetscScalar d2pde_dmi_dmj(int iii, int jjj,
                             const unsigned int &field_id, const Point &q_point)
   { return CALL_MEMBER_FN_W_REF(this->m_MathModel,
                                 this->m_MathModel.d2pde_dmi_dmj[iii][jjj])
                                                    (field_id,q_point); }
  ///** wrapper should be inlined @ref PDEModelInlineStrategy */
  //virtual PetscScalar getWeightedSensitivity(const unsigned int &i)
  // { return this->m_MathModel.getWeightedSensitivity(i); }

  ///** wrapper should be inlined @ref PDEModelInlineStrategy */
  //virtual PetscScalar getWeightedDifference(const unsigned int &i)
  // { return this->m_MathModel.getWeightedDifference(i); }

  ///** wrapper should be inlined @ref PDEModelInlineStrategy */
  //virtual unsigned int  n_p_dofs(const unsigned int &i)
  // { return this->m_MathModel.n_p_dofs[i]; }

  ///** wrapper should be inlined @ref PDEModelInlineStrategy */
  //virtual unsigned int  n_u_dofs(const unsigned int &i)
  // { return this->m_MathModel.n_u_dofs[i]; }

  ///** wrapper should be inlined @ref PDEModelInlineStrategy */
  //virtual Real psi(const unsigned int &var, 
  //                 const unsigned int &dof,const unsigned int &qp)
  // { return (*this->m_MathModel.psi[var])[dof][qp] ; }

  /** this class holds all constitutive data.
    * should be an instaniation of 
    *     - PennesStandardDiffusionApproximation
    *     - PennesVoltage 
    *     - PennesLineSource 
    *     - VerifySourceTerm 
    */
  MathematicalModel m_MathModel; 

protected:
  

  /**
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void init_data ();

/*
  not needed??? why put solution into all vectors?

   * Re-update the local values when the mesh has changed.
   * This method takes the data updated by \p update() and
   * makes it up-to-date on the current mesh.
  virtual void re_update ();
   */
};
#endif
