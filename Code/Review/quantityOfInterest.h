#ifndef __quantityOfInterest_h
#define __quantityOfInterest_h
#include "qoiBaseClass.h"  // interface to qoi 
#include "Imaging.h"       // itk imaging data strucs

/**
 * Class to load MRTI data onto mesh only 
 *   allows for Background Phase Correction Computions
 *    should do no temperature predictions
 * @todo{
 *        templated on the m_pdeSolver. can still instantiate 
 *        with TransientFEMSystem as the base class for
 *        generality but can also test with particular Solver 
 *        to test for potential performance benifit 
 *      }
 */
template< typename PDEModelSolver  >
class onlyMRTIQOI : public qoiBaseClass
{
public : 
  //constructor
  onlyMRTIQOI(AppSolve *user,GetPot &controlfile, PetscInt iii); 

  /** method to extract image data to a NumericVector */
  PetscErrorCode ExtractImageDataFromBackground(InputImageType::Pointer, PetscInt ,
                                                TransientFEMSystem &);

  /** load imaging data */
  virtual PetscErrorCode PriorLoadQoiData(AppSolve *);   

  /** setup qoi specific auxillary systems */
  virtual PetscErrorCode SetupAuxSystems( EquationSystems &) ;

  /** qoi specific solution */
  virtual PetscErrorCode solveqoi(AppSolve *user, GetPot &);

protected : 
  /** method to extract image data to a NumericVector */
  void ExtractThermalImageTimeInstance(AppSolve*,PetscInt);

  PDEModelSolver  *m_pdeSolver; 
  PetscScalar     m_sigma; 
};
/**
 * Class for forward solve only 
 *  should require no stored vectors
 */
template< typename PDEModelSolver  >
class noOptQOI : public onlyMRTIQOI < PDEModelSolver >
{
public : 
  noOptQOI(AppSolve *user,GetPot &controlfile, PetscInt iii):
   onlyMRTIQOI < PDEModelSolver >::onlyMRTIQOI(user,controlfile,iii) {}

  /** setup Petsc Parallel Vector */
  virtual PetscErrorCode SetupCntrl(AppSolve *,GetPot &); 

  /** clean up Petsc Parallel Vector */
  virtual PetscErrorCode DestroyPetsc();

  /** setup qoi specific auxillary systems */
  virtual PetscErrorCode SetupAuxSystems( EquationSystems &) ;

  /** qoi specific solution */
  virtual PetscErrorCode solveqoi(AppSolve *user, GetPot &);

  /** initial guess */
  virtual void GetInitialSolutionVec(AppSolve *user);

  /** forward stepping */
  virtual PetscErrorCode ForwardSolve(AppSolve *);

  /** mainly for verif */
  virtual PetscErrorCode PriorLoadQoiData(AppSolve *);   
};
// l2 spacetime norm
PetscScalar SpaceTimeNorm(AppSolve*, TransientFEMSystem & );
/**
 *  Main Class for inverse problem
 */
template< typename PDEModelSolver  >
class spaceTimeQOI : public noOptQOI < PDEModelSolver >
{
public : 
  /** constructor */
  spaceTimeQOI(AppSolve *user,GetPot &controlfile, PetscInt iii): 
    noOptQOI < PDEModelSolver >::noOptQOI (user,controlfile,iii) {}
  
  virtual PetscErrorCode SetupTAO(AppSolve *); ///< setup TAO
  virtual PetscErrorCode SolveTAO(AppSolve *); ///< run TAO solvers
  virtual PetscErrorCode DestroyTAO();   ///< clean up TAO

  /** setup qoi specific auxillary systems */
  virtual PetscErrorCode SetupAuxSystems( EquationSystems &) ;

  /** qoi specific solution */
  virtual PetscErrorCode solveqoi(AppSolve *user, GetPot &);

  //call previous base class
  virtual PetscErrorCode PriorLoadQoiData(AppSolve *user)
    { return this->onlyMRTIQOI< PDEModelSolver >::PriorLoadQoiData(user) ; } 

  /* accumulate the adjoint rhs inherent to the QOI */
  virtual void accumulateAdjointQOI(ADJOINTLOADARG);

  /* accumulate the adjoint rhs inherent to the Sensitivity */
  virtual void accumulateAdjointSensitivity(ADJOINTLOADARG);

  /** 
   * evaluation of the qoi of the form of a space time norm
   */
  virtual PetscScalar ComputeObjective(AppSolve* user,
               TransientFEMSystem & uncertainty_system)
    { return SpaceTimeNorm(user,uncertainty_system ); } 

  /** accumulate the load for the sensitivity solve  */
  virtual void accumulateSensitivityLoad(AppSolve *,const QGauss &, 
                                 std::vector<PetscScalar> &elemParameter,
                                 const unsigned int &, 
                                 const std::vector<Real>&,
                                 const std::vector<Point>&,
                                 std::vector< DenseSubVector<Number> > &,
                                 TransientFEMSystem &);

  /** element wise gradient computations */
  void accumulateGradient(AppSolve *, const QGauss &,
                               const unsigned int &,const std::vector<Real>&,
                               const std::vector<Point>&,DenseVector<Number> &);

  /** element wise hessian computations */
  void accumulateHessian(AppSolve *, const QGauss &,
                           std::vector<PetscScalar> &,
                           const unsigned int &,const std::vector<Real>&,
                           const std::vector<Point>&, DenseVector<Number> &);
  /** element wise sensitivity gradient */
  PetscScalar accumulateSensitivityGradient(AppSolve *, const QGauss &,
                           const unsigned int &,const std::vector<Real>&,
                           const std::vector<Point>&);

  /** setup adjoint pde solver */
  virtual void SetupAdjoint(AppSolve *); 
};
///*Read Sensitivities from Disk to accumulate Hessian*/ 
//class qoiSpaceTimeDiskRead: public spaceTimeQOI
//{
//public : 
//  qoiSpaceTimeDiskRead(AppSolve *user,GetPot &controlfile, PetscInt iii): 
//                          spaceTimeQOI(user,controlfile,iii) {};//constructor
//
//  /* initialize sensitivity for each control dof and open file from disk */
//  virtual void initializeSensitivity( TransientInverseLinearImplicitSystem & ,
//                                                               const PetscInt);
//  /* read the sensitivity  from disk*/
//  virtual void getSensitivity(  AppSolve * ,const int &globalDof,
//              const    std::vector<optimizationParameter*>::iterator );
//
//  /* store the sensitivity to disk*/
//  virtual void storeSensitivity( AppSolve *, const int &,
//                                             NumericVector<Number>& );
//  /* cleanup this sensitivity */
//  virtual void finalizeSensitivity();
//private : 
//  // file handle information for petsc viewer
//  PetscViewer localViewHandle; 
//  int         localfileHandle;
//  PetscInt    rows,type,tr[2];
//};
/*Recompute Sensitivities from Disk to accumulate Hessian*/ 
template< typename PDEModelSolver  >
class qoiSpaceTimeRecompute: public spaceTimeQOI < PDEModelSolver >
{
public : 
  //constructor
  qoiSpaceTimeRecompute(AppSolve *user,GetPot &controlfile, PetscInt iii): 
      spaceTimeQOI< PDEModelSolver >::spaceTimeQOI(user,controlfile,iii) {};

  /* recompute the sensitivity */
  virtual void getSensitivity(  AppSolve * ,const int &globalDof,
              const    std::vector<optimizationParameter*>::iterator );
};
/********** dddas class ***************/
template< typename PDEModelSolver  >
class dddasQOI : public noOptQOI < PDEModelSolver >
{
public : 
  dddasQOI(AppSolve *user,GetPot &controlfile, PetscInt iii);//constructor

  virtual PetscErrorCode PriorLoadQoiData(AppSolve *)
                               { PetscFunctionReturn(0); }
  virtual PetscScalar ComputeObjective(AppSolve*,
                                       TransientFEMSystem & ) {return 0.0;} 
  //virtual void SendToServer(PetscInt *commands, PetscInt Message)
  // { MPI_Send(commands,2,MPI_INT,0,Message,AppSolve::DDDAS_COMM_WORLD); }
  // do nothing
  virtual void accumulateAdjointQOI(ADJOINTLOADARG){};
  virtual void accumulateAdjointSensitivity(ADJOINTLOADARG){};

};
/********** class for FD verification of imaging data ***************/
template< typename PDEModelSolver  >
class spaceTimeFDVerification : public spaceTimeQOI < PDEModelSolver >
{
public : 
  spaceTimeFDVerification(AppSolve *,GetPot &, PetscInt );

  // overwrite the FD gradient verification routine
  virtual PetscErrorCode compareFDGrad(TAO_APPLICATION, Vec ,Vec,AppSolve *);

  /** dirichlet domains for optimization problem
   *  dirichlet domains should impose adjoint solution = 0 
   *  and zero the corresponding gradient anyways
   *  mainly used for verification problems and hessian vector product 
   */
  void DirichletDomains(AppSolve *, Vec , const PetscScalar );

};

template< typename PDEModelSolver  >
class VerificationQOI : public spaceTimeFDVerification < PDEModelSolver >
{
/* -------------------------------------------------------------------- 
    TODO: This class is a hodgepodge of hacks w/ switch statements to
          verify solutions against a manufactured or analytical solution
   -------------------------------------------------------------------- */ 
public : 
  VerificationQOI(AppSolve *,GetPot &, PetscInt );

  // overwrite the data load routine
  virtual PetscErrorCode PriorLoadQoiData(AppSolve *); 

  /* accumulate the adjoint rhs inherent to the QOI
     *FIXME* will need to implement a new derived class for each different QOI
                   - this may or may not work out well...  */
  virtual void accumulateAdjointQOI(ADJOINTLOADARG);

};

#endif
