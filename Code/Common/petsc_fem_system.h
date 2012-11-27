#ifndef __petsc_fem_system_h__
#define __petsc_fem_system_h__
// C++ includes
#include <vector>

// libmesh includes
#include "numeric_vector.h"
#include "fem_system.h"
#include "getpot.h"

// local includes
#include "tttkUtilities.h"

class optimizationParameter; //forward declaration
#define CALL_MEMBER_FN(object,ptrToMember)  ((object)->*(ptrToMember))
#define CALL_MEMBER_FN_W_REF(object,ptrToMember)  ((object).*(ptrToMember))

/**
 * Common Routines for libmesh with PETSc Solvers 
 * additional vectors/matrices may be added,
 * as offered in the parent classes.
 * 
 * use SNESLineSearchNoNorms to treat SNES as a linear solver
 *  
 *  -snes_ls basicnonorms -snes_no_convergence_test -snes_max_it 1
 *  
 *  -snes_max_it 1                   # Do a maximum of one linear solve
 *  -snes_ls basicnonorms            # Don't do a line search and don't even
 *                                   # compute the residual (3)
 *  -snes_convergence_test skip      # Skip the convergence test before
 *                                   # evaluating the Jacobian.  SNES will
 *                                   # normally bail out early if it starts
 *                                   # with a sufficiently small residual.
 */
class PetscFEMSystem : public FEMSystem
{
public:

  /**
   * Constructor.  Initializes required
   * data structures.
   */
  PetscFEMSystem (EquationSystems& es,
                      const std::string& name,
                      const unsigned int number ) ;

  // called by destructor
  void clear()
  {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    for (unsigned int var=0; var<this->n_vars(); var++)
       {
          ierr = VecScatterDestroy(m_subvector_scatter[var]); 
          ierr = ISDestroy(m_subvector_is[var]);       
       }
    Parent::clear();
    PetscFunctionReturnVoid();
  }

  /** regression */
  virtual void verifySolution( EquationSystems &){}  

  /**
   * FEMSystem::assembly is overwritten to add dirichlet BC and 
   * to enable special processing for real-time
   */
  virtual void assembly(bool,bool);

  /** setup scatter context to vector of subvariabls */
  void SetupSubVector (); 
  
  /** get/set subvector for given variable ID */
  void GetSubVector (Vec &,int ); 
  void SetSubVector (Vec &,int ); 
  IS & GetSubVectorIS (int var ) {return m_subvector_is[var];} 
  
  /** Dirichlet BC */
  virtual void ApplyDirichlet (){}; 
  
  /** pre FEMSystem sensitivity load */
  virtual void sensitivityLoad( NumericVector<Number>& )
    {libmesh_error();return;}

  /** setup optimization variables */
  virtual void SetupOptimizationVariables( std::vector<optimizationParameter*> &)
    {libmesh_error();return;}

  /** typedefs to make function pointer declarations more readable */
  typedef PetscScalar(PetscFEMSystem::*ICMemFn)(unsigned int,unsigned int,
                                           const Point&, const Parameters& ) ;

  /** initial condition function pointers */
  std::vector<ICMemFn>  InitValues;
 
  /** voltage initial condition. overwritten in derived class*/
  virtual PetscScalar getInitialVoltage(unsigned int,unsigned int,
                                        const Point&,const Parameters& )
       {
         libmesh_error();return 0.0;
       }
  /** temperature initial condition. overwritten in derived class*/
  virtual PetscScalar getInitialTemperature(unsigned int,unsigned int, 
                                            const Point&  , const Parameters& )
       {
         libmesh_error();return 0.0;
       }

  /** initial damage */
  virtual PetscScalar getInitialDamage(unsigned int,unsigned int,
                                       const Point&,const Parameters& )
       {
         return 0.0;
       }
 
  /** fluence initial condition. overwritten in derived class*/
  virtual PetscScalar getInitialFluence(unsigned int,unsigned int, 
                                            const Point&  , const Parameters& )
       {
         return 0.0;
       }

  /** flux initial condition. overwritten in derived class*/
  virtual PetscScalar getExternalFlux_X(unsigned int,unsigned int, 
                                      const Point&  , const Parameters& )
       { libmesh_error();return 0.0; }
  virtual PetscScalar getExternalFlux_Y(unsigned int,unsigned int, 
                                      const Point&  , const Parameters& )
       { libmesh_error();return 0.0; }
  virtual PetscScalar getExternalFlux_Z(unsigned int,unsigned int, 
                                      const Point&  , const Parameters& )
       { libmesh_error();return 0.0; }
  /** flux initial condition. overwritten in derived class*/
  virtual PetscScalar getInterstitialFlux_X(unsigned int,unsigned int, 
                                          const Point&,const Parameters& )
       { libmesh_error();return 0.0; }
  virtual PetscScalar getInterstitialFlux_Y(unsigned int,unsigned int, 
                                          const Point&,const Parameters& )
       { libmesh_error();return 0.0; }
  virtual PetscScalar getInterstitialFlux_Z(unsigned int,unsigned int, 
                                          const Point&,const Parameters& )
       { libmesh_error();return 0.0; }
  /** flux initial condition. overwritten in derived class*/
  virtual PetscScalar getExternalIrradiance(unsigned int,unsigned int, 
                                            const Point&  , const Parameters& )
       {
         libmesh_error();return 0.0;
       }
  /** flux initial condition. overwritten in derived class*/
  virtual PetscScalar getInterstitialIrradiance(unsigned int,unsigned int, 
                                            const Point&  , const Parameters& )
       {
         libmesh_error();return 0.0;
       }
  /**
   * setup dirichlet data
   */
  virtual void SetupDirichlet(libMesh::MeshBase& );

  /**
   * setup intial conditions 
   */
  virtual void SetupInitialConditions ();

  //-----------------------------------------------------------------
  // access to the solution data fields
  
  /**
   * @returns the solution for the given ID
   * for the specified global DOF.
   */
  Number stored_solution (const unsigned int &global_dof_number,
                          const unsigned int id_vec=0) const ;
  Number old_solution    (const unsigned int &,
                          const unsigned int ) const  
   {
     libmesh_error(); // old_local_solution not used
     return 0.0;
   }
  Number prev_solution   (const unsigned int &global_dof_number) const ;

  /**
   * wrapper around solver to call regression tests
   *  also used to call PetscLinearSolver in the form 
   *  of a nonlinear solve
   *   \f[
   *      f(x_{n+1}) = f(x_n) + J(x_n) (x_{n+1} - x_n) \approx 0
   *      \qquad \Rightarrow \qquad
   *      J(x_n) ( x_n - x_{n+1} )  \equiv J(x_n) s_n =  f(x_n) 
   *      \qquad \Rightarrow \qquad
   *      x_{n+1} = x_n - s_n 
   *   \f]
   */
  virtual void solve();

  /** store nodeset data */
  std::map<PetscInt, std::map<PetscInt, std::vector<PetscInt> > > m_NodeSets;
  std::map<PetscInt, std::vector<PetscInt> > m_DirichletNodeSets;

  /** for Kalman system we want system dynamics matrix only */
  bool mass_residual (bool request_jacobian, DiffContext &)
    {
      std::cout << "handling time stepping directly";
      std::cout << "should not be called...";
      libmesh_error();
      return request_jacobian;
    }
  // control jacobian assembly
  void AssembleMassMatrix(){ m_MassMatrix      = PETSC_TRUE ; 
                             m_StiffnessMatrix = PETSC_FALSE;}
  void AssembleStiffnessMatrix(){ m_MassMatrix      = PETSC_FALSE; 
                                  m_StiffnessMatrix = PETSC_TRUE ;}
  void AssembleFullJacobian(){ m_MassMatrix      = PETSC_TRUE ; 
                               m_StiffnessMatrix = PETSC_TRUE ;}
  void SetJacobianComputed(PetscTruth state){m_jacobianComputed=state; }
  void SetPowerID(int CurrentTimeID){m_PowerID=CurrentTimeID; }

  //default is to recompute jacobian matrices
  virtual PetscTruth linearStatePDE()
   {  
    if( !m_LinearSolve )
     { // this is a nonlinear jacobian must recompute
       return PETSC_FALSE;
     }
    else if(m_jacobianComputed)
     { // the jacobian is linear and the jacobian has already been computed
       return PETSC_TRUE;
     }
    else 
     { // the jacobian is linear but we still need to compute jacobian once
       m_jacobianComputed = PETSC_TRUE;
       return PETSC_FALSE;
     }
   } 
  virtual void printSelf(std::ostream& os)
    { 
      os << "PetscFEMSystem:  linearPDE          =" << m_LinearPDE        << std::endl;
      os << "PetscFEMSystem:  m_LinearSolve      =" << m_LinearSolve      << std::endl;
      os << "PetscFEMSystem:  m_MassMatrix       =" << m_MassMatrix       << std::endl;
      os << "PetscFEMSystem:  m_StiffnessMatrix  =" << m_StiffnessMatrix  << std::endl;
      os << "PetscFEMSystem:  m_jacobianComputed =" << m_jacobianComputed << std::endl;
      printStdVector< Real >(os, "PDE.m_newton_coeff[", m_newton_coeff );
      printStdVector< Real >(os, "PDE.m_u_infty["     , m_u_infty      );
      printStdVector< Real >(os, "PDE.m_NeumannFlux[" , m_NeumannFlux  );
      return; 
    }

  virtual void ScatterParametersLocally() = 0 ; 
  // member function pointers for boundary conditions
  typedef PetscScalar (PetscFEMSystem::*ResidualBCMemFn)(const unsigned int,  const Real & );
  typedef PetscScalar (PetscFEMSystem::*JacobianBCMemFn)(const unsigned int);

  /** 
   * boundary conditions member function pointers to 
   *     - jacobianNothingBC
   *     - jacobianCauchyBC
   *     - residualNothingBC
   *     - residualNeumannBC
   *     - residualCauchyBC
   */
  std::vector< std::vector<ResidualBCMemFn> > ResidualBC;
  std::vector< std::vector<JacobianBCMemFn>  >  JacobianBC;
  
  /// BC computes diagonal terms only, off diag terms computed in derived class
  bool side_time_derivative(bool request_jacobian,DiffContext &context);

  // return true if this is a linear solve
  PetscTruth  LinearPDE()
         {return m_LinearPDE;}

  //boundary conditions  can be overridden in derived class if necessary
  // the member function pointer will call the derived class function
  virtual PetscScalar residualNothingBC(const unsigned int ,const Real &,
                                        int , const int ,
                                        const Point &, const Parameters& )
    { return 0.0; }
  virtual PetscScalar residualNeumannBC(const unsigned int i_var,const Real &,
                                        int , const int ,
                                        const Point & , const Parameters&)
    { return m_NeumannFlux[i_var]; }
  /**
   * Cauchy boundary data for each variable
   * \f[
   *    h(u - u_\infty)
   * \f]
   */
  virtual PetscScalar residualCauchyBC( const unsigned int i_var, const Real &temperature,
                                        int , const int ,
                                        const Point & , const Parameters& )
    { return m_newton_coeff[i_var]*(temperature-m_u_infty[i_var]) ; }

  virtual PetscScalar jacobianNothingBC(const unsigned int , int,
                                        const Point & , const Parameters& )
    { return 0.0; }
  virtual PetscScalar jacobianCauchyBC(const unsigned int i_var, int,
                                       const Point &, const Parameters& )
    { return m_newton_coeff[i_var]; }
  PetscInt    get_num_elem_blk(){     return n_block ; }

protected:
  
  PetscInt n_block;     ///< # of mesh blocks

  /**
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void init_data ();

  std::vector< std::vector<unsigned int> > m_subvector_indices;
  std::vector<VecScatter>                  m_subvector_scatter;
  std::vector<IS>                          m_subvector_is;

  unsigned int m_PowerID; ///<  time step in power model to use in the assembly

  /** error tolerance for regression tests */
  PetscScalar errorTol;
  
  // recomputation flags
  PetscTruth  m_LinearSolve      ,    ///< use KSPSolve
              m_MassMatrix       ,    ///< assemble mass Matrix only
              m_StiffnessMatrix  ,    ///< assemble Stiffness Matrix only
              m_jacobianComputed; ///< jacobian has been computed already

  // dirichlet, neumann, cauchy BC Data
  std::vector<Real>  m_newton_coeff,   ///<  \f$ h  \f$ 
                     m_u_infty,        ///<  \f$ u_\infty  \f$ 
                     m_NeumannFlux;    ///<  \f$ \mathcal{G}  \f$ 

  PetscTruth  m_LinearPDE;

private:
  /// BC data does note need update by default
  virtual void UpdateBoundaryData(DiffContext &,const unsigned int,
                                  const int ,  const int ,
                                  const Point &, const Parameters& ){}

  /** The type of the parent. */
  typedef FEMSystem Parent;
  
  // petsc profile logging tools
  PetscInt  PetscFEMSysLogResidual,
            PetscFEMSysLogJacobian,
            PetscFEMSysLogDefault ;
};
#endif
