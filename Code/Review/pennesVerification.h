#ifndef __pennes_Verification_h__
#define __pennes_Verification_h__

#include "pennesInverseModel.h"  // interface to pennes solver 

/* -------------------------------------------------------------------- 
   Verifications problems using fortran modules
    TODO: This class is a hodgepodge of hacks w/ switch statements to
          verify solutions against a manufactured or analytical solution
   -------------------------------------------------------------------- */ 
class PennesVerification: public PennesInverseModel
{
public:

 //constructor
 PennesVerification(GetPot &, libMesh::MeshBase& );

 // transfer C++ to F90 data structures 
 virtual void XferData();

 //setup the equation system
 virtual void SetupState(AppSolve* ) ;

 // jacobian is same as in base class only need to change residual
 virtual void accumulateResidual(unsigned int , 
                                 DiffContext &,TransientFEMSystem &);

 //accumulate the jacobian 
 virtual void accumulateJacobian(unsigned int , 
                                 DiffContext &,TransientFEMSystem &);


 // mass
 virtual PetscScalar getMass();

 virtual PetscErrorCode WriteControlFile(PetscInt, PetscInt );

 //accumulate the adjoint matrix and rhs inherent to the pde
 virtual void accumulateAdjointPDE(const QGauss &, 
                                const unsigned int &, 
                                const std::vector<Real>&,
                                std::vector< SubMatrixVector > &,
                                std::vector< DenseSubVector<Number> > &,
                                TransientFEMSystem &,
                                TransientFEMSystem &) ;

 //print data 
 virtual void printSelf(std::ostream& os=std::cout) ; 

 // verify the solution
 virtual void verifySolution( EquationSystems &);  

 // overwrite exactsolution
 virtual PetscScalar exactSolution(const Point& p,
                                   const Parameters& parameters,
                                   const std::string& )
  { return this->getExactTemperature(0,p,parameters); }

 virtual PetscScalar getExactTemperature(unsigned int, 
                                         const Point&  , const Parameters& );
 // special boundary conditions
 virtual PetscScalar residualCauchyBC( const unsigned int i_var,
                                        const Real &temperature);
 virtual PetscScalar residualNeumannBC(const unsigned int i_var,const Real &dum);

 virtual PetscScalar jacobianCauchyBC(const unsigned int i_var);
 /** call base class to avoid "hidden error"
  * http://www.parashift.com/c++-faq-lite/strange-inheritance.html#faq-23.9
  */
 virtual void jacobianCauchyBC( 
                    const unsigned int i_var, QGauss &qface,
                    const std::vector<Real>& JxW_face ,
                    const std::vector<std::vector<Real> >& phi_face,
                    std::vector< SubMatrixVector > &Kij)
   {
    PDEModelBaseClass::jacobianCauchyBC(i_var,qface,JxW_face ,
                                        phi_face, Kij) ; return;
   }
 


 virtual PetscScalar dpde_dw_0( const unsigned int &,const Point &);
 virtual PetscScalar dpde_dx_0( const unsigned int &,const Point &);
 virtual PetscScalar dpde_dy_0( const unsigned int &,const Point &);
 virtual PetscScalar dpde_dz_0( const unsigned int &,const Point &);
 virtual PetscScalar dpde_dw_2( const unsigned int &,const Point &);
 virtual PetscScalar dpde_dw_N( const unsigned int &,const Point &);
 virtual PetscScalar dpde_dw_I( const unsigned int &,const Point &);
 virtual PetscScalar dpde_dw_D( const unsigned int &,const Point &);
 virtual PetscScalar dpde_dk_0( const unsigned int &,const Point &);
 virtual PetscScalar dpde_dk_1( const unsigned int &,const Point &);
 virtual PetscScalar dpde_dmu_a(const unsigned int &,const Point &);
 virtual PetscScalar dpde_dmu_s(const unsigned int &,const Point &);

protected:
 
// arguement list of functers for residual evaluation
#define FORMFCNGAUSSARG const PetscScalar*,const PetscInt*,const PetscScalar*,PetscScalar*,const PetscInt*,PetscScalar*,PetscScalar*,PetscScalar*
#define JACGAUSSARG const PetscScalar &,const unsigned int &, const unsigned int &, std::vector< unsigned int>  &, std::vector< const std::vector<std::vector<Real> >* > &, std::vector< const std::vector<std::vector<RealGradient> >* > &d, const std::vector<Real>& , std::vector<Real>   &, std::vector<Real>   &, std::vector<RealGradient> &, DenseSubMatrix<Number>&
// arguement list of functers for BC evaluation
#define BCGAUSSARG const PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*

  // typedefs to make function pointer declarations more readable
  typedef void (*FormFcnGaussType)(FORMFCNGAUSSARG) ;
  typedef void (*JacGaussType)(JACGAUSSARG) ;
  typedef void (*BcGaussType)(BCGAUSSARG) ;

  /*function pointers for evaluations at the gauss points for element routine
    that forms nonlinear function for SNESsolve 
                Formfcngauss is a function pointer to: getverifformfcn
                                                       nonlinpennesisolaser
                                                       nonlinpennesmonte   
                                                       pennesidealmonte */

  std::vector < std::vector < FormFcnGaussType > > formfcngauss;
  /*function pointers for evaluations at the gauss points for element routine
    that forms nonlinear function jacobian for SNESsolve 
                         nonlinpennesjac  : jacobian for pennes model
                         getverifjac      : jacobian for verification problem */
  std::vector < JacGaussType > jacgauss;

  /*function pointers for evaluations at the gauss points for element routine
    that forms the load for the dual problem inherent to the pde 
                          nonlinpennesadj  :  pennes adjoint
                          getverifdualjac  :  verification problems for adjoint */
  void (*adjmat)(JACGAUSSARG) ;
  /* boundary condition function pointers to one of the following routines
                                 getverifbc           : verification problem bc
                                 heattransbc          : bc for heat transfer
                                 getverifbcdual       : bc for adjoint prob
  */
  std::vector< BcGaussType > fcnbcgauss;
  std::vector< BcGaussType > jacbcgauss;
  std::vector< BcGaussType > adjbcgauss;
private:
  PetscErrorCode SetupPennesLaser(AppSolve *);
  void gradSolutionBufferSetup(std::vector<PetscScalar> &,const unsigned int &);
};
void JacPennesLaser(        JACGAUSSARG);
void JacPennesRF(           JACGAUSSARG);
void JacPenalty(            JACGAUSSARG);
void VerifJacPennesLaser(   JACGAUSSARG);
void AdjMatPennesLaser(     JACGAUSSARG);
void AdjMatPennesRF(        JACGAUSSARG);
void VerifAdjMatPennesLaser(JACGAUSSARG);

void pennes_verification(PetscInt, PetscInt);
void pennes_sensitivity_load(OptGaussType, NumericVector<Number>& );
Number pennes_adjoint_exact    (const Point& , const Parameters& ,
                                const std::string& , const std::string& );
Number pennes_ideal_field      (const Point& , const Parameters& ,
                                const std::string& , const std::string& );
#endif
