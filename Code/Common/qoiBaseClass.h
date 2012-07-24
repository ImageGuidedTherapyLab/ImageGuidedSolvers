#ifndef __qoiBaseClass_h
#define __qoiBaseClass_h
#include "Imaging.h"
/*----------------------------------------------------------------------
     base class for quantities of interest
  
     bridge pattern implemented for qoi and pde heirarchies
  
     http://www.parashift.com/c++-faq-lite/multiple-inheritance.html

     for example see:

               spaceTimeQOI::accumulateAdjointQOI
  
     the pdeBaseClass pointer is used to call any simple operations on the
     solutions at the gauss points, ie

                         pdeSolver->getWeightedDifference(0); 
  
     latest revision    - june 09
  
     purpose            - this is the main data structure that allows 
                          communication between libMesh/PETSC and TAO.
                          this is used for the optimizations
  ---------------------------------------------------------------------- */
class qoiBaseClass 
{
  friend PetscErrorCode TaoConverged_CheckPoint(TAO_SOLVER ,void *);
public : 
  qoiBaseClass(AppSolve *,GetPot &, PetscInt);//constructor

  /** Destructor  */
  virtual ~qoiBaseClass(){ if(this->Images) delete this->Images; } 

  /** setup Imaging data structures*/
  void init_imaging(PetscInt ,GetPot &);

  /** setup qoi specific auxillary systems */
  virtual PetscErrorCode SetupAuxSystems( EquationSystems &) = 0;

  /** clean up */
  virtual PetscErrorCode DestroyPetsc(){return 0;}

  /** qoi specific solution */
  virtual PetscErrorCode solveqoi(AppSolve *user, GetPot &) = 0 ;

  /** setup adjoint pde solver */
  PetscErrorCode physicalBounds(libMesh::MeshBase& ,Vec ,Vec );

  /** plot data already stored */
  void PlotStoredData( AppSolve *); 

  /** compute and plot in interleaved fasion */
  void PlotInterleaveCompute( AppSolve *); 

  /** communicate with server */
  virtual void SendToServer(PetscInt*,PetscInt){};

  /** method to extract image data to a NumericVector */
  PetscErrorCode ExtractImageData(InputImageType::Pointer,TransientFEMSystem &);

  // echo data
  virtual void printSelf(std::ostream& os=std::cout) ; 
  void echoOptWindow(std::ostream&, PetscInt , PetscInt );

  // iteration number
  PetscInt getIteration(); 

  /** forward stepping */
  virtual PetscErrorCode ForwardSolve(AppSolve *) {return 0;}

  /** function pointer for evaluation of the qoi */
  virtual PetscScalar ComputeObjective(AppSolve*,
                                       TransientFEMSystem & ) {return 0.0;} 

  /*function pointers for setting up the qoi data structures
      PriorLoadQoiData is used to load the data BEFORE the solve so that the 
                       data will be available to compute the QOI at each time step
                       this avoid unnecessary communication with project_vector
      ComputeQoiData is used to compute the QOI from the solution AFTER the solve
        the available functers are-
                     verifcomputeidealfield: for verification problems
                     readMRTI_data         : for calibration problems
                     LoadIdealArr          : arrhenius base optimize
                     LoadIdeal_TS          : two-state base optimize
                     LoadIdealHSP          : HSP base optimize
                     LoadIdealTmp          : for temp base optimize */
  virtual PetscErrorCode PriorLoadQoiData(  AppSolve*) {return 0;} 
  PetscErrorCode (*ComputeQoiData)(    AppSolve*);
  PetscErrorCode (*PostSolveQoiData)(  AppSolve*);
  /*function pointer for killing the current optimization solve 
      and moving to the next optimization window may point to:
                                     ReturnPETSC_FALSE  
                                     checkMRTI_data4_IC  
                                     checkMRTI_data4OPT  */
  virtual PetscTruth  UpdateQoiData(AppSolve *){return PETSC_FALSE;};

  PetscInt NDOF_CONTROL[3]; 
          /* NDOF_???[0]: local # of dof OWNED by this processor in a Petsc Vec
             NDOF_???[1]: total # of dof shared amongst processors 
             NDOF_???[2]: total # of dof OWNED AND SHARED by this processor */
  Vec        CNTL_LOC;// local scratch storage for optimization solve
  VecScatter GLOCSCAT_CNTL;/*scattering context of local petsc vector 
       		    GRAD_LOC to/from petsc parallel vectors for
       		    control problem*/
  PetscInt        minHessianColumn , // hack for computing hessian subset
                  maxHessianColumn ; // hack for computing hessian subset
  Mat        Hessian,// storage for Hessian matrix 
             HessianShell;// shell for Hessian-vector product 

  PetscErrorCode PutLocalCntrlVars(Vec);

  // initialize solution
  virtual void GetInitialSolutionVec (AppSolve *) {}

  /* compute sensitivity for each control dof */
  void computeSensitivity(  AppSolve * , Vec ); 

  /* initialize sensitivity default does nothing*/
  virtual void initializeSensitivity(TransientFEMSystem &, 
                                     const PetscInt){}
  /* dflt does nothing and compute only variance */
  virtual void getSensitivity(  AppSolve * ,const int &,
          const    std::vector<optimizationParameter*>::iterator ){}

  /* dflt does nothing */
  virtual void storeSensitivity( AppSolve *, const int &,
                                 NumericVector<Number>& ){}

  /* cleanup this sensitivity... default cleanup is not necessary
  */
  virtual void finalizeSensitivity(){};

  /** accumulate the load for the sensitivity solve pde deriv */
  virtual void accumulateSensitivityLoad(AppSolve *,const libMesh::QGauss &,
                             std::vector<PetscScalar> &,
                             const unsigned int &, 
                             const std::vector<Real>&,
                             const std::vector<Point>&,
                             std::vector< DenseSubVector<Number> > &,
                             TransientFEMSystem &) {}

  virtual void accumulateGradient(AppSolve *, const libMesh::QGauss &,
                                   const unsigned int &, 
                                   const std::vector<Real>& ,
                                   const std::vector<Point>& ,
                                   DenseVector<Number> &) {}
  virtual void DirichletDomains(AppSolve*,Vec,const PetscScalar) {}

  /** element wise hessian computations */
  virtual void accumulateHessian(AppSolve *, const libMesh::QGauss &,
                         std::vector<PetscScalar> &,
                         const unsigned int &,const std::vector<Real>&,
                         const std::vector<Point>&, DenseVector<Number> &) {}
  /** element wise sensitivity gradient */
  virtual PetscScalar accumulateSensitivityGradient(AppSolve *, const libMesh::QGauss &,
                           const unsigned int &,const std::vector<Real>&,
                           const std::vector<Point>&) {return 0.0;}
#define ADJOINTLOADARG AppSolve *, const libMesh::QGauss &, std::vector<PetscScalar> &, const unsigned int &, const std::vector<Real>&, const std::vector<Point>&, std::vector< DenseSubVector<Number> > & 

  typedef void (qoiBaseClass::*adjointLoadMemFn)(ADJOINTLOADARG);
  adjointLoadMemFn accumulateAdjointLoad;

  /* accumulate the adjoint rhs inherent to the QOI */
  virtual void accumulateAdjointQOI(        ADJOINTLOADARG) {}
  virtual void accumulateAdjointSensitivity(ADJOINTLOADARG) {}

  /* return  the optimization paramter pertaining to this global dof */
  std::vector<optimizationParameter*>::iterator getParameterPointer(const int&);

  // for verification
  virtual PetscErrorCode compareFDGrad(TAO_APPLICATION, Vec ,Vec,
                                         AppSolve *){ return 0; };

  // get parameter mapping
  void getGlobalMap(const Elem* , std::vector<PetscInt> &);

  // get element parameter values
  void getGlobalParamValues(const Elem* , std::vector<PetscScalar> &);

  unsigned int GetParamSize(){return Parameters.size();}

  virtual PetscScalar dqoi_dm(OPTGAUSSARG){return 0.0;}

  PetscInt TIMELAG,     // time lag for optimial control computations
           fncEvalCnt;  // counter
  PetscTruth  PLOTOPTSOLVE;   // control over plotting

  // write variance data
  virtual void plotElemData(OStringStream &, libMesh::MeshBase &, Vec );

  /* reference to TAO application context */
  TAO_APPLICATION &  TAOAPP(){ return taoapp; }

  /* reference to TAO application context */
  TAO_SOLVER &  TAO(){ return tao; }

  /* get TAO solver method */
  std::string GetTaoSolverMethod()
    { 
      TaoMethod methodType;
      TaoGetMethod(tao,&methodType); 
      return std::string(methodType);
    }

  PetscInt Nsteplo(){return m_Nsteplo;} 
  PetscInt Nstephi(){return m_Nstephi;} 
           
  void SetNsteplo(PetscInt value){m_Nsteplo = value ;return;} 
  void SetNstephi(PetscInt value){m_Nstephi = value ;return;} 

  std::vector<PetscScalar> TIMEERR;
  /** use in debugging */
  PetscInt hessianColumn;
  /** PETSC PARALLEL VECTORS */
  Vec     hessianDiagonal;      

  /** get imaging pointer */
  Imaging *ImagingPointer(){return Images;}

protected : 
  /** Pointer to itk based imaging data*/
  Imaging *Images;

   // keep track of which mrti data is in memory
   std::vector<PetscInt> MRTI_MEM; 

  /** if(m_BackgroundCorrection) background correction needed */
  bool m_BackgroundCorrection ;
  /** TAO application context */
  TAO_APPLICATION taoapp; 
  /** TAO_SOLVER solver context */
  TAO_SOLVER tao;  

  PetscErrorCode GetLocalVars_LB(Vec);  // deprecated
  PetscErrorCode GetLocalVars_UB(Vec);  // deprecated
  PetscErrorCode GetLocalCntrlVars(Vec);
  PetscErrorCode GetCntrlVars(libMesh::MeshBase &,Vec);
  PetscErrorCode PutVerifCntrlVars();
  PetscErrorCode GetVerifDeltaCntrl(Vec,PetscInt); 
  PetscErrorCode AddRegularization(PetscScalar&,PetscScalar&);
  PetscErrorCode RegularizationDeriv(PetscScalar&,PetscScalar&);

  // element wise storage for gradient and error estimates
  std::vector<PetscScalar> GradBeta; 
  PetscInt Init_Nsteplo; ///< initial lower bound of the fem time steps
  PetscInt m_Nsteplo,   ///< lower bound of fem time step
           m_Nstephi;   ///< upper bound of fem time step

  PetscTruth TAO_MONITOR , // monitoring
             Penalty ; // Penalty  = PETSC_TRUE ==> apply penalty function 
  //PETSC Sequential VECTORS
  Vec        QCONTRL ;// main global vector for optimization solve
  std::vector<optimizationParameter*> Parameters;
  /* hessian contributions of the QOI - all function pointers for
     optimization should have the SAME argument lists given by OPTGAUSSARG */ 

  //boost::numeric::ublas::symmetric_matrix< OptGaussType,
  //                                          boost::numeric::ublas::lower > d2qoi_dmi_dmj;
  // FIXME - with regularization should be a matrix of functers
  virtual PetscScalar d2qoi_dmi_dmj(OPTGAUSSARG){return 0.0;}; 

  std::vector<optimizationParameter*>::iterator IterParam;
  std::vector<PetscInt>::iterator DofIter;
  std::vector<PetscInt> locmap; //store mapping from the local dof to global dof
  PetscInt rank,size;   // mpi info

};
Number GetITKImageData (const Point& ,
                        const Parameters& ,
                        const std::string& ,
                        const std::string& );
Number project_one (const Point& ,
                    const Parameters& ,
                    const std::string& ,
                    const std::string& );
PetscErrorCode ForwardSolve(AppSolve *);
PetscErrorCode FormObjectiveAndGradient(TAO_APPLICATION,Vec,double*,Vec,void*);
PetscErrorCode FormObjective(TAO_APPLICATION,Vec,double*,void*);
PetscErrorCode SkipObjective(TAO_APPLICATION,Vec,double*,void*);
PetscErrorCode ZeroGradient(TAO_APPLICATION,Vec,Vec,void*);
PetscErrorCode FormGradient(TAO_APPLICATION,Vec,Vec,void*);
PetscErrorCode PhysBnds(TAO_APPLICATION,Vec,Vec,void*);
PetscErrorCode FormHessian(TAO_APPLICATION,Vec,Mat*,Mat*,MatStructure*,void*);
PetscErrorCode MatrixFreeHessian(TAO_APPLICATION,Vec,Mat*,Mat*,MatStructure*,void*);
PetscErrorCode hessianVectorProduct(Mat , Vec , Vec );
PetscErrorCode hessianVectorDiagonal(Mat , Vec );
#endif
