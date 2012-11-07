#ifndef __pdeBaseClass_h
#define __pdeBaseClass_h
// libmesh types
#include "exodusII_io.h" 
#include "numeric_vector.h"
#include "vector_value.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "o_string_stream.h"
#include "getpot.h"

// For systems of equations the \p DenseSubMatrix
// and \p DenseSubVector provide convenient ways for
// assembling the element matrix and vector on a
// component-by-component basis.
#include "dense_submatrix.h"
#include "dense_subvector.h"

// class for optimization paramters
#include "optimizationParameter.h"

class TransientFEMSystem; //forward declarations
namespace libMesh
{
class QGauss; //forward declarations
}

/* macros to define the arguement list of functions used as function pointers*/
// arguement list of functers for jacobian evaluation
typedef std::vector< DenseSubMatrix<Number> >  SubMatrixVector;

/**@ingroup TreatmentPlanning 
 * Base Class for PDE Models
 * Assumes similar boundary conditions for all derived classes
 * 
 * @section PDEModelInlineStrategy Inlining of Member Functions
 * 
 * http://www.parashift.com/c++-faq-lite/inline-functions.html#faq-9.9
 * 
 * [9.9] With inline member functions that are defined outside the class, is
 * it best to put the inline keyword next to the declaration within the class
 * body, next to the definition outside the class body, or both?
 * 
 * Best practice: only in the definition outside the class body.
 * 
 * \verbatim
 *  class Foo {
 *  public:
 *    void method();  <- best practice: don't put the inline keyword here
 *    ...
 *  };
 *  
 *  inline void Foo::method()  <- best practice: put the inline keyword here
 *  { ... }
 * \endverbatim
 * 
 * Here's the basic idea:
 * 
 *     - The public: part of the class body is where you describe the
 *     observable semantics of a class, its public member functions, its friend
 *     functions, and anything else exported by the class. Try not to provide
 *     any inklings of anything that can't be observed from the caller's code.
 *     - The other parts of the class, including non-public: part of the class
 *     body, the definitions of your member and friend functions, etc. are pure
 *     implementation. Try not to describe any observable semantics that were
 *     not already described in the class's public: part.
 * 
 * From a practical standpoint, this separation makes life easier and safer
 * for your users. Say Chuck wants to simply "use" your class. Because you
 * read this FAQ and used the above separation, Chuck can read your class's
 * public: part and see everything he needs to see and nothing he doesn't
 * need to see. His life is easier because he needs to look in only one spot,
 * and his life is safer because his pure mind isn't polluted by
 * implementation minutiae.
 * 
 * Back to inline-ness: the decision of whether a function is or is not
 * inline is an implementation detail that does not change the observable
 * semantics (the "meaning") of a call. Therefore the inline keyword should
 * go next to the function's definition, not within the class's public: part.
 * 
 * NOTE: most people use the terms "declaration" and "definition" to
 * differentiate the above two places. For example, they might say, "Should I
 * put the inline keyword next to the declaration or the definition?"
 * Unfortunately that usage is sloppy and somebody out there will eventually
 * gig you for it. The people who gig you are probably insecure, pathetic
 * wannabes who know they're not good enough to actually acomplish something
 * with their lives, nonetheless you might as well learn the correct
 * terminology to avoid getting gigged. Here it is: every definition is also
 * a declaration. This means using the two as if they are mutually exclusive
 * would be like asking which is heavier, steel or metal? Almost everybody
 * will know what you mean if you use "definition" as if it is the opposite
 * of "declaration," and only the worst of the techie weenies will gig you
 * for it, but at least you now know how to use the terms correctly. 
 * @todo { test if caching current_local_solution into 
 *         DiffContext::elem_solution,
 *         DiffContext::elem_subsolutions would make a differenc in
 *         evaluateSolutionForState,evaluateSolutionForState... 
 *         ie copy FEMContext::interior_value
 *         }
 */
class PDEModelBaseClass
{
public:
  PDEModelBaseClass(GetPot &,EquationSystems &); // constructor

  /**
   * @returns a reference to this system's parent EquationSystems object.
   */
  EquationSystems & get_equation_systems() { return _equation_systems; }

  /**
   *  Write the element data
   */
  void writeElementData(OStringStream&, libMesh::MeshBase&, std::vector<Number>&, Real);

  /**
   *  Default exact solution value returns 0
   */
  virtual Number exactSolution(const Point&,
                               const Parameters&,
                               const std::string& ) {return 0.0;} 	 

  /**
   *  Default does nothing for old fortran compatibility
   */
  virtual void fillSolutionBuffer( 
          std::vector<PetscScalar> &, const unsigned int ){}
  /**
   *  Source Term
   */
  virtual PetscScalar  PennesSource(   const unsigned int &, const Real&,
                                    const Real&, const Real&,
                                    const Point &, const int ) 
                                   {libmesh_error();return 0.0;}

  virtual PetscScalar dPennesSourcedu( const unsigned int &, const Real&,
                                    const Real&,const Real&, const Real&,
                                    const Point &, const int ) 
                                   {libmesh_error();return 0.0;}

  /** default regression doesnt nothing */
  virtual void verifySolution( EquationSystems &){}
  
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
  virtual PetscScalar residualFluenceBC(const unsigned int , const Real &,
                                       int , const int ,
                                       const Point &, const Parameters& )
    { libmesh_error(); return 0.0; }
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
  virtual PetscScalar jacobianFluenceBC(const unsigned int , int,
                                        const Point & , const Parameters& )
    { libmesh_error();return 0.0; }
  virtual PetscScalar jacobianCauchyBC(const unsigned int i_var, int,
                                       const Point &, const Parameters& )
    { return m_newton_coeff[i_var]; }

  //boundary conditions  can be overridden in derived class if necessary
  // the member function pointer will call the derived class function
  virtual void adjointLoadNothingBC(const unsigned int ,libMesh::QGauss &,
                                    const std::vector<std::vector<Real> >&  , 
                                    const std::vector<Real>& ,
                                    std::vector< DenseSubVector<Number> > & ,
                                    TransientFEMSystem &){};

  virtual void adjointLoadCauchyBC(const unsigned int ,libMesh::QGauss &,
                                   const std::vector<std::vector<Real> >&  , 
                                   const std::vector<Real>& ,
                                   std::vector< DenseSubVector<Number> > & ,
                                   TransientFEMSystem &);

  virtual void gradientNothingBC(const unsigned int ,libMesh::QGauss &,
                                 const std::vector<std::vector<Real> >&, 
                                 const std::vector<std::vector<Real> >&, 
                                 const std::vector<Real>&, DenseVector<Number> &,
                                 TransientFEMSystem &,
                                 TransientFEMSystem &){};

  virtual void gradientNeumannBC(const unsigned int ,libMesh::QGauss &,
                                 const std::vector<std::vector<Real> >&, 
                                 const std::vector<std::vector<Real> >&, 
                                 const std::vector<Real>&, DenseVector<Number> &,
                                 TransientFEMSystem &,
                                 TransientFEMSystem &);

  virtual void gradientCauchyBC(const unsigned int i_var,libMesh::QGauss &,
                                const std::vector<std::vector<Real> >&, 
                                const std::vector<std::vector<Real> >&, 
                                const std::vector<Real>&, DenseVector<Number> &,
                                TransientFEMSystem &,
                                TransientFEMSystem &);
  /** Print the Constitutive Data */
  virtual void printSelf(std::ostream& );

  // member function pointers for boundary conditions
  typedef PetscScalar (PDEModelBaseClass::*ResidualBCMemFn)(const unsigned int, const Real &,
                                                            int , const int ,
                                                            const Point & , const Parameters&) ;
  typedef PetscScalar (PDEModelBaseClass::*JacobianBCMemFn)(const unsigned int, int,
                                                            const Point &, const Parameters& );

 // member function pointers for boundary conditions
 typedef void (PDEModelBaseClass::*pdeJacobianBCMemFn)(const unsigned int , libMesh::QGauss &,
                                    const std::vector<Real>& ,
                                    const std::vector<std::vector<Real> >& ,
                                    std::vector< SubMatrixVector > &) ;
 typedef void (PDEModelBaseClass::*pdeAdjointLoadBCMemFn)(const unsigned int,libMesh::QGauss &,
                                    const std::vector<std::vector<Real> >&  , 
                                    const std::vector<Real>& ,
                                    std::vector< DenseSubVector<Number> > & ,
                                    TransientFEMSystem &);
 typedef void (PDEModelBaseClass::*pdeGradientBCMemFn)(const unsigned int i_var,
                               libMesh::QGauss &,
                               const std::vector<std::vector<Real> >&, 
                               const std::vector<std::vector<Real> >&, 
                               const std::vector<Real>&, DenseVector<Number> &,
                               TransientFEMSystem &,
                               TransientFEMSystem &);
#define CALL_MEMBER_FN(object,ptrToMember)  ((object)->*(ptrToMember))
#define CALL_MEMBER_FN_W_REF(object,ptrToMember)  ((object).*(ptrToMember))

  /** 
   * boundary conditions member function pointers to 
   *     - jacobianNothingBC
   *     - jacobianCauchyBC
   *     - residualNothingBC
   *     - residualNeumannBC
   *     - residualCauchyBC
   */
  ResidualBCMemFn    ResidualBC[5];
  JacobianBCMemFn    JacobianBC[5];
  //boundary conditions
  pdeJacobianBCMemFn    accumulateJacobianBC[5];
  pdeAdjointLoadBCMemFn accumulateAdjointLoadBC[5];
  pdeGradientBCMemFn    accumulateGradientBC[5];

  /** nodal dirichlet data */
  virtual bool dirichletNodalData( Point &) {return false;}
 
  /** typedefs to make function pointer declarations more readable */
  typedef PetscScalar(PDEModelBaseClass::*ICMemFn)(unsigned int,unsigned int,
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
  *  derivatives of pde wrt the parameters 
  *   these are implemented at virtual member function pointers
  *     the base class function pointer is overridden in the derived class
  *   NEED one base class function for each possible parameter
  */
  // mixed derivatives
  virtual PetscScalar d2qoi_du_dm(  OPTGAUSSARG ){return 0.0;}
  virtual PetscScalar d2pde_du_dm(  OPTGAUSSARG ){return 0.0;}
  virtual PetscScalar dpde_derr(       OPTGAUSSARG ){return 0.0;}
  virtual PetscScalar dpde_dx_0(       OPTGAUSSARG ){return 0.0;}
  virtual PetscScalar dpde_dy_0(       OPTGAUSSARG ){return 0.0;}
  virtual PetscScalar dpde_dz_0(       OPTGAUSSARG ){return 0.0;}
  virtual PetscScalar dpde_dw_0(       OPTGAUSSARG ){return 0.0;}
  virtual PetscScalar dpde_dw_2(       OPTGAUSSARG ){return 0.0;}
  virtual PetscScalar dpde_dw_N(       OPTGAUSSARG ){return 0.0;}
  virtual PetscScalar dpde_dw_I(       OPTGAUSSARG ){return 0.0;}
  virtual PetscScalar dpde_dw_D(       OPTGAUSSARG ){return 0.0;}
  virtual PetscScalar dpde_dk_0(       OPTGAUSSARG ){return 0.0;}
  virtual PetscScalar dpde_dk_1(       OPTGAUSSARG ){return 0.0;}
  virtual PetscScalar dpde_ds_0(       OPTGAUSSARG ){return 0.0;}
  virtual PetscScalar dpde_dmu_a(      OPTGAUSSARG ){return 0.0;}
  virtual PetscScalar dpde_dmu_s(      OPTGAUSSARG ){return 0.0;}
  virtual PetscScalar dpde_dpow(       OPTGAUSSARG ){return 0.0;}
  virtual PetscScalar dpde_dprobeTemp( OPTGAUSSARG ){return 0.0;}
  virtual PetscScalar dpde_dchi(       OPTGAUSSARG ){return 0.0;}

  // second derivative wrt to state
  // default is linear pde with vanishing second state derivatives
  virtual PetscScalar d2pde_d2u(  const unsigned int &,const unsigned int &, OPTGAUSSARG ){return 0.0;} 
 
  // sensitivity gradient for verification
  virtual PetscScalar dqoidu_dudm( const unsigned int &){return 0.0;}
 
  // mixed derivatives
  virtual PetscScalar d2pde_du_dk_0(OPTGAUSSARG ){return 0.0;}
  virtual PetscScalar d2pde_du_ds_0(OPTGAUSSARG ){return 0.0;}
  virtual PetscScalar d2pde_du_dw_0(OPTGAUSSARG ){return 0.0;}
  virtual PetscScalar d2pde_du_dchi(OPTGAUSSARG ){return 0.0;}
 
  //  second partials with respect to parameters
  virtual PetscScalar d2pde_d2m(         OPTGAUSSARG ){return 0.0;}
  virtual PetscScalar d2pde_d2mu_a(      OPTGAUSSARG ){return 0.0;}
  virtual PetscScalar d2pde_d2mu_s(      OPTGAUSSARG ){return 0.0;}
  virtual PetscScalar d2pde_dmu_a_dmu_s( OPTGAUSSARG ){return 0.0;}

   void SetupOptimizationVariables(libMesh::MeshBase &,
                 std::vector<optimizationParameter*> &) ;

   virtual PetscScalar getMass(){return 1.0;}

   // load constitutive data from disk
   virtual PetscErrorCode GetDataFromDisk(libMesh::MeshBase &,std::string &);

  /** hessian contributions of the PDE - all function pointers for
   *  optimization should have the SAME argument lists given by OPTGAUSSARG 
   *  
   *  symmetric matrix storage not used
   */ 
   std::vector< std::vector< dpde_dmMemFn > > d2pde_dmi_dmj; 
/*
When a symmetric matrix is stored in upper-packed storage mode, the upper
triangular part of the symmetric matrix is stored, including the diagonal, in a
one-dimensional array. The upper triangle is packed by columns. (This is
equivalent to packing the lower triangle by rows.) The matrix is packed
sequentially column by column in n(n+1)/2 elements of a one-dimensional array.
To calculate the location of each element aij of matrix A in an array AP using
the upper triangular packed technique, use the following formula:

AP(i+(j(j-1)/2)) = aij    where j >= i
This results in the following storage arrangement for the elements of a
symmetric matrix A in an array AP:

Following is an example of a symmetric matrix that uses the element values to
show the order in which the matrix elements are stored in the array. Given the
following matrix A:

                      ┌                    ┐
                      |  1   2   4   7  11 |
                      |  2   3   5   8  12 |
                      |  4   5   6   9  13 |
                      |  7   8   9  10  14 |
                      | 11  12  13  14  15 |
                      └                    ┘

the array is:

           AP = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)

Following is an example of how to transform your symmetric matrix to
upper-packed storage mode:

       K = 0
       DO 1 J=1,N
          DO 2 I=1,J
             K = K+1
             AP(K)=A(I,J)
   2      CONTINUE
   1   CONTINUE

*/
  Real  TimeDerivativeScalingFactor() 
         {return m_TimeDerivativeScalingFactor;}

  // return true if this is a linear solve
  PetscTruth  LinearPDE()
         {return m_LinearPDE;}

  /**
   *  return PETSC_TRUE if we are solving time dependet equations  
   *  on this domain
   */
  PetscTruth  TransientTerm(const unsigned int DomainID) 
                  {return m_TransientDomain[DomainID];}

  PetscInt    get_num_elem_blk(){     return n_block ; }
protected:

  /**
   * Constant reference to the \p EquationSystems object
   * used for the simulation.
   */
  EquationSystems& _equation_systems;

  PetscInt n_block;     ///< # of mesh blocks

  /**
   *  return PETSC_TRUE if solving time dependent equations  
   *  on this domain
   */
  std::vector<PetscTruth>  m_TransientDomain;  

  // dirichlet, neumann, cauchy BC Data
  std::vector<Real>  m_newton_coeff,   ///<  \f$ h  \f$ 
                     m_u_infty,        ///<  \f$ u_\infty  \f$ 
                     m_NeumannFlux;    ///<  \f$ \mathcal{G}  \f$ 

  PetscTruth  m_LinearPDE;

  /**
   *  mass term to scale to time derivative to 1, \f$ \dot{u} = F(u)\f$
   *  @todo {m_TimeDerivativeScalingFactor will break verification problems}
   */
  Real  m_TimeDerivativeScalingFactor;

  // store a pointer to all field parameters for plotting
  std::vector<optimizationParameter*> _fieldParameters;


};

#endif
