// system includes
//#include <math.h>

// C++ include files 
#include <iostream>
#include <fstream>
//#include <string>
#include <vector>
#include <algorithm>

// libMesh include files
#include "libmesh.h"
#include "mesh.h"
#include "equation_systems.h"
#include "boundary_info.h"
#include "fe.h"
#include "fe_interface.h"
#include "dof_map.h"
#include "quadrature_gauss.h"
#include "getpot.h"
#include "petsc_vector.h"
#include "petsc_matrix.h"
#include "steady_solver.h"
#include "petsc_diff_solver.h"
#include "mesh_function.h"
#include "fe_compute_data.h"
#include "point_locator_base.h"
#include "parallel.h" // mpi utilities
#include "private/kspimpl.h" // mpi utilities
// petsc basic copy
PetscErrorCode MatCopy_Basic(Mat A,Mat B,MatStructure str);

// local includes
#include "pennesModel.h" // Constitutive data
#include "petsc_fem_context.h" 
#include "thermal_therapy_system.h"
#include "pennesSystem.h"
#include "pennesSystem.txx"
#include "KalmanFilter.h"

// global initial tolerance
const PetscScalar rtol0 = 1.e-9;

/* ------------------------------------------------------------------------

    FD pennes model
    the partial differential equation
  
            \rho dudt - k Laplacian u = q_laser
  
    with dirichlet boundary conditions
   
             u = 0  for  x = x_min, x = x_max, 
                         y = y_min, y = y_max, 
                         z = z_min, z = z_max
  
    A finite difference approximation with the usual 7-point stencil
    is used to discretize the boundary value problem to obtain a linear 
    system of equations.


  ------------------------------------------------------------------------- */

/* constructor */
KalmanFilter::KalmanFilter(EquationSystems *EqnSystemsPtr)
{
  /*
    Gaussian of the form
     f(x) = 1/\sqrt{2 \pi P} \exp(-1/2/P (x -\mu)^2 ) 

    sigma^2 = P

    (\mu- \sigma,\mu+ \sigma)  ---> 0.683 confidence
    (\mu-2\sigma,\mu+2\sigma)  ---> 0.997 confidence
 
  */
  m_eqnSystems = EqnSystemsPtr;
  modelcov = 2.5*2.5;     // [deg C]
  statecov = 4.0    ;     // [deg C]
  SumCov    = NULL; 
  tmpMatDenState = NULL;
  CovP      = NULL;
  CovR      = NULL;
  CovQ      = NULL;
  m_zerotol = 1.e-6;
  m_fill    = 1.5;

  // kalman profile
  PetscLogEventRegister("KFSetup"       ,PETSC_VIEWER_COOKIE,&kfLogSetup       );
  PetscLogEventRegister("KFStatePredict",PETSC_VIEWER_COOKIE,&kfLogStatePredict);
  PetscLogEventRegister("KFCovarPredict",PETSC_VIEWER_COOKIE,&kfLogCovarPredict);
  PetscLogEventRegister("KFStateUpdate" ,PETSC_VIEWER_COOKIE,&kfLogStateUpdate );
  PetscLogEventRegister("KFCovarUpdate" ,PETSC_VIEWER_COOKIE,&kfLogCovarUpdate );
  PetscLogEventRegister("KFInvSumCov"   ,PETSC_VIEWER_COOKIE,&kfLogInvSumCov   );

  // command line parameters
  PetscTruth  flg;
  PetscErrorCode ierr;
  ierr=PetscOptionsGetScalar(PETSC_NULL,"-statecov",&statecov,PETSC_NULL);
  ierr=PetscOptionsGetScalar(PETSC_NULL,"-modelcov",&modelcov,PETSC_NULL);
  ierr=PetscOptionsGetScalar(PETSC_NULL,"-sparsezerotol" ,&m_zerotol,PETSC_NULL);
  ierr=PetscOptionsGetScalar(PETSC_NULL,"-matmatmultfill",&m_fill   ,PETSC_NULL);

}

void KalmanFilter::printSelf(std::ostream& os)
{
 PetscErrorCode ierr;
  os << "KalmanFilter:   modelcov =" <<   modelcov   << std::endl;
  os << "KalmanFilter:   statecov =" <<   statecov   << std::endl;
  os << "KalmanFilter:   zerotol  =" << m_zerotol    << std::endl;
  os << "KalmanFilter:   fill     =" << m_fill       << std::endl;
 return;
}

// solve A X = B using dense solver
#undef __FUNCT__  
#define __FUNCT__ "MatMultipleRHSSolve"
PetscErrorCode KalmanFilter::MatMultipleRHSSolve(Mat A, Mat B, Mat X)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  Mat FactorMatrix;
  ierr = MatGetFactor(A,this->solvertype, 
                      MAT_FACTOR_LU,&FactorMatrix); 
  IS isrow,iscol;
  ierr = MatGetOrdering(A,this->ordertype,&isrow,&iscol);

  MatFactorInfo factorInfo;
  ierr = MatFactorInfoInitialize(&factorInfo); 

  // factor
  MatLUFactorSymbolic(FactorMatrix,A,isrow,iscol,&factorInfo);
  ierr = MatLUFactorNumeric(FactorMatrix,A,&factorInfo);
 	 
  // solve
  PetscPrintf(PETSC_COMM_WORLD,"MatMultipleRHSSolve...\n");
  MatMatSolve(FactorMatrix,B,X);

  ierr = MatDestroy(FactorMatrix);CHKERRQ(ierr);
  ierr = ISDestroy(isrow);CHKERRQ(ierr);
  ierr = ISDestroy(iscol);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
// solve A X = B using KSP
//   if B == NULL  ==>  B == I  ==>  X == A^-1
// for A a sparse matrix, A O(N) * N < O(N^3)
#undef __FUNCT__  
#define __FUNCT__ "MatKSPSolveInverse_Basic"
PetscErrorCode MatKSPSolveInverse_Basic(Mat A, Mat B, Mat X)
{
  PetscErrorCode ierr;
  Vec            x,canonicalBasis,columnVec;
  PetscInt       m,N,i;
  PetscScalar    *xx;

  PetscFunctionBegin;
  ierr = MatGetArray(X,&xx);CHKERRQ(ierr);
  ierr = MatGetLocalSize(X,&m,PETSC_NULL);CHKERRQ(ierr);  /* number local rows */
  ierr = MatGetSize(X,PETSC_NULL,&N);CHKERRQ(ierr);       /* total columns in dense matrix */
  ierr = VecCreateMPIWithArray(((PetscObject)A)->comm,m,PETSC_DETERMINE,PETSC_NULL,&x);CHKERRQ(ierr);

  /* Krylov subspace method context */
  KSP            localKSP;   
  ierr = KSPCreate(PETSC_COMM_WORLD,&localKSP);CHKERRQ(ierr);
  ierr = KSPSetOperators(localKSP,A,A,
         // same non zero pattern for the multiple solves
         SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  // set defaults
  PC localPC;
  ierr = KSPGetPC(localKSP,&localPC); CHKERRQ(ierr);
  // get PC, default to block jacobi
  PetscTruth  flg;
  char    ksp_inverse_pc[PETSC_MAX_PATH_LEN];
  sprintf(ksp_inverse_pc,"%s",PCBJACOBI);
  ierr = PetscOptionsGetString(PETSC_NULL,"-ksp_inverse_pc",ksp_inverse_pc,
                                           PETSC_MAX_PATH_LEN,&flg);
  ierr = PCSetType(localPC,ksp_inverse_pc); CHKERRQ(ierr);
  PetscInt maxiter;
  ierr = MatGetSize(A,&maxiter,PETSC_NULL);CHKERRQ(ierr);
  ierr = KSPSetTolerances(localKSP,rtol0,PETSC_DEFAULT,PETSC_DEFAULT,
                          maxiter); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(localKSP);CHKERRQ(ierr);
  ierr = KSPMonitorCancel(localKSP);CHKERRQ(ierr);
  localKSP->printreason = PETSC_FALSE;

  // create vector for multiplication
  ierr = MatGetVecs(A,&canonicalBasis,PETSC_NULL); CHKERRQ(ierr);
  ierr = VecDuplicate( canonicalBasis,&columnVec); CHKERRQ(ierr);

  // O(n) * n should still be O(n^2) < O(n^3) for dense LA
  PetscPrintf(PETSC_COMM_WORLD,"%d total solves using %s pc... \n",N,ksp_inverse_pc);
  for (i=0; i<N; i++) {
    ierr = VecZeroEntries(canonicalBasis);
    ierr = VecSetValue(canonicalBasis,i,1.0,INSERT_VALUES);
    ierr = VecAssemblyBegin(canonicalBasis); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(  canonicalBasis); CHKERRQ(ierr);
    if( B == NULL )
     {
      ierr = VecCopy(  canonicalBasis,columnVec);CHKERRQ(ierr);
     }
    else
     {
      ierr = MatMult(B,canonicalBasis,columnVec);CHKERRQ(ierr);
     }
    ierr = VecPlaceArray(x,xx + i*m);CHKERRQ(ierr);
    ierr = KSPSolve(localKSP,columnVec,x);CHKERRQ(ierr);
    //check convergence reason
    KSPConvergedReason kspReason;
    PetscInt kspIterNum;
    KSPGetConvergedReason(localKSP, &kspReason);
    if(kspReason<0)
      {
       std::cerr<<"MatKSPSolveInverse_Basic Diverged "<< kspReason
                << std::endl <<std::flush; 
       // write matrix to a file to examine in matlab
       PetscViewer matViewer; 
       ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"IllMatrix.m",
                                   &matViewer);  CHKERRQ(ierr);
       ierr = PetscViewerSetFormat(matViewer,
                     PETSC_VIEWER_ASCII_MATLAB); CHKERRQ(ierr);
       ierr = MatView(A,matViewer); CHKERRQ(ierr);
       ierr = PetscViewerFlush(matViewer); CHKERRQ(ierr);
       ierr = PetscViewerDestroy(matViewer); CHKERRQ(ierr);
       PetscPrintf(PETSC_COMM_WORLD,"wrote IllMatrix.m ...\n" );
       libmesh_error(); 
      } 
    else
      { 
       KSPGetIterationNumber(localKSP, &kspIterNum);
      } 
    ierr = VecResetArray(x);CHKERRQ(ierr);
    if(i%100 == 0)
     PetscPrintf(PETSC_COMM_WORLD,"(%d,%d) ",i,kspIterNum);// print progress
  }
  PetscPrintf(PETSC_COMM_WORLD,"\n");// print progress
  ierr = VecDestroy(x);CHKERRQ(ierr);
  ierr = VecDestroy(canonicalBasis); CHKERRQ(ierr);
  ierr = VecDestroy(     columnVec); CHKERRQ(ierr);
  ierr = MatRestoreArray(X,&xx);CHKERRQ(ierr);

  PetscTruth  assembled ; 
  ierr = MatAssembled(X,&assembled); CHKERRQ(ierr);
  if(!assembled)
    { 
     ierr = MatAssemblyBegin(X,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
     ierr = MatAssemblyEnd(  X,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    }

  // clean up
  ierr = KSPDestroy(localKSP); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
//setup kalman infrastructure
#undef __FUNCT__
#define __FUNCT__ "KalmanFilter::Setup"
PetscErrorCode KalmanFilter::Setup(PetscInt numMeasurements)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  PetscLogEventBegin(kfLogSetup,0,0,0,0); // setup

  // default solver to petsc
  sprintf(solvertype,"%s",MAT_SOLVER_PLAPACK);
  sprintf(solvertype,"%s",MAT_SOLVER_MUMPS);
  sprintf(solvertype,"%s",MAT_SOLVER_SUPERLU_DIST);
  sprintf(solvertype,"%s",MAT_SOLVER_PETSC);
  PetscTruth  flg;
  ierr = PetscOptionsGetString(PETSC_NULL,"-solver",solvertype,
                                           PETSC_MAX_PATH_LEN,&flg);

  // default ordering to nested dissection
  sprintf(ordertype,"%s",MATORDERING_RCM);
  sprintf(ordertype,"%s",MATORDERING_ND);
  ierr = PetscOptionsGetString(PETSC_NULL,"-ordering",ordertype,
                                           PETSC_MAX_PATH_LEN,&flg);

  // get state system
  typedef LITTSystem< PennesStandardDiffusionApproximation > LITTSDASystem;
  LITTSDASystem &state_system = m_eqnSystems->get_system< LITTSDASystem >("StateSystem");

  // get temperature vector
  Vec templateVector;
  state_system.GetSubVector(templateVector, state_system.u_var);

  ierr = VecScatterCreateToZero( templateVector  , &gather,&imageVec);
  ierr = VecDuplicate( templateVector , &globalVec); CHKERRQ(ierr);
  CHKERRQ(ierr);

  ierr = VecDestroy(templateVector);CHKERRQ(ierr);

  // get reference to the Sub Vector IS. DONT DESTROY
  IS &TemperatureIS = state_system.GetSubVectorIS(state_system.u_var);
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create and store system matricies
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  Mat stiffnessMatrix;
  state_system.AssembleMassMatrix(); 
  state_system.assembly(false,true); 
  state_system.matrix->close();

  ierr = MatGetSubMatrix(
  (dynamic_cast< PetscMatrix<double>* > (state_system.matrix) )->mat() ,
                         TemperatureIS , TemperatureIS ,
                         MAT_INITIAL_MATRIX,&m_MassMatrix );

  // notice that the copy is done w/ the deltat scaling already included
  ierr = MatDuplicate( m_MassMatrix, MAT_COPY_VALUES,&m_CovJacobian);
  ierr = MatDuplicate( m_MassMatrix, MAT_COPY_VALUES,&m_CovRHSMat  );
  
  // remove the deltat scaling in the mass matrix to be used in the covariance
  // propagation
  ierr = MatScale(m_MassMatrix, state_system.deltat);CHKERRQ(ierr);

  CHKERRQ(ierr);

  state_system.AssembleStiffnessMatrix(); 
  state_system.assembly(false,true); 
  state_system.matrix->close();
  ierr = MatGetSubMatrix(
  (dynamic_cast< PetscMatrix<double>* > (state_system.matrix) )->mat() ,
                         TemperatureIS , TemperatureIS ,
                         MAT_INITIAL_MATRIX,&stiffnessMatrix );

  CHKERRQ(ierr);


  AutoPtr<DiffContext> con = state_system.build_context();
  PetscFEMContext &_femcontext = libmesh_cast_ref<PetscFEMContext&>(*con);
  state_system.init_context(_femcontext);

  // remove the theta scaling and multiply by -1 b/c the stiffness matrix is
  // on the RHS... this is the opposite side from which the jacobian was
  // initially computed 
  ierr = MatScale(stiffnessMatrix, -1.0/ 
	             _femcontext.ThetaValue(0) );CHKERRQ(ierr);

  // build matrices needed for covariance propagation
  // m_CovJacobian p^i = m_CovRHSMat p^{i-1} + m_MassMatrix q
  ierr = MatAXPY(m_CovJacobian,-1.0,stiffnessMatrix,
                            SAME_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = MatAXPY(m_CovRHSMat  , 1.0,stiffnessMatrix,
                            SAME_NONZERO_PATTERN);CHKERRQ(ierr);

  // apply dirichlet data if any
  // TODO Dirichlet implementations seem to be causing numerical instabilities
  //  leave out for now...
  //if( state_system.m_dirichletNodes.size() )
  //  {
  //    ierr = MatZeroRows(m_CovJacobian, state_system.m_dirichletNodes.size(),
  //                                     &state_system.m_dirichletNodes[0], 1.0);
  //    ierr = MatZeroRows(m_CovRHSMat  , state_system.m_dirichletNodes.size(),
  //                                     &state_system.m_dirichletNodes[0], 1.0);
  //  }

  PetscTruth  writeSystemMatrix=PETSC_FALSE;
  ierr = PetscOptionsGetTruth(PETSC_NULL,"-write_system_matrix",&writeSystemMatrix,PETSC_NULL); 
  if(writeSystemMatrix)
    { 
     PetscViewer matViewer; 
     // write inverse
     ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"SystemMatrix.m",
                                 &matViewer);  CHKERRQ(ierr);
     ierr = PetscViewerSetFormat(matViewer,
                   PETSC_VIEWER_ASCII_MATLAB); CHKERRQ(ierr);
     ierr = MatView(m_CovJacobian,matViewer); CHKERRQ(ierr);
     ierr = PetscViewerFlush(matViewer); CHKERRQ(ierr);
     ierr = PetscViewerDestroy(matViewer); CHKERRQ(ierr);
     PetscPrintf(PETSC_COMM_WORLD,"wrote SystemMatrix.m ...\n" );

     // write initial matrix in binary to read in
     ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, "SystemMatrix.bin",
                  FILE_MODE_WRITE,&matViewer);CHKERRQ(ierr);
     ierr = MatView(m_CovJacobian ,matViewer);CHKERRQ(ierr);
     ierr = PetscViewerFlush(matViewer); CHKERRQ(ierr);
     ierr = PetscViewerDestroy(matViewer);CHKERRQ(ierr);

     PetscPrintf(PETSC_COMM_WORLD,"wrote SystemMatrix.bin ...\n" );
    } 
  // reset to the full jacobian matrix
  state_system.AssembleFullJacobian(); 

  // state matrix dimensions
  PetscInt M_state,N_state,m_state,n_state;
  ierr = MatGetSize(m_MassMatrix,&M_state,&N_state);CHKERRQ(ierr);
  ierr = MatGetLocalSize(m_MassMatrix,&m_state,&n_state);CHKERRQ(ierr);

  // create dense matrix used in state propagation
  ierr = MatCreateMPIDense(PETSC_COMM_WORLD,m_state,n_state,M_state,N_state,
                                   PETSC_NULL,&tmpMatDenState);CHKERRQ(ierr);
  //// F = -M^-1 * K
  //if(!this->m_consistentMassMatrix)
  //  {// use lumped mass approximation
  //   ierr=PetscPrintf(PETSC_COMM_WORLD,
  //        "  using Lumped Mass matrix\n"); CHKERRQ(ierr);
  //   // multiply consistent mass matrix by all ones and invert
  //   Vec OneVec, MassDiagonal;
  //   ierr = MatGetVecs(m_MassMatrix,&OneVec,&MassDiagonal); 
  //   CHKERRQ(ierr);
  //   // lumped mass matrix formed by row sum of consistent mass
  //   // matrix entries and storing in the diagonal
  //   ierr = VecSet(OneVec,1.0);CHKERRQ(ierr);
  //   ierr = MatMult(m_MassMatrix,OneVec,MassDiagonal);CHKERRQ(ierr);
  //   // invert the matrix
  //   ierr = VecReciprocal(MassDiagonal);CHKERRQ(ierr);
  //   // multiply stiffness matrix by the inverse
  //   ierr = MatDuplicate(stiffnessMatrix,MAT_COPY_VALUES,&m_SystemDynamics);
  //   ierr = MatDiagonalScale(m_SystemDynamics,MassDiagonal,PETSC_NULL);CHKERRQ(ierr);
  //   
  //   ierr = VecDestroy(OneVec      );CHKERRQ(ierr);
  //   ierr = VecDestroy(MassDiagonal);CHKERRQ(ierr);
  //  }
  //else
  //  {// use consistent mass matrix
  //   ierr=PetscPrintf(PETSC_COMM_WORLD,
  //        "  using Consistent Mass matrix\n"); CHKERRQ(ierr);
  //   // create temporary dense scratch matrix to hold solution System Dynamics
  //   // matrix
  //   Mat localMatDen;
  //   ierr = MatCreateMPIDense(PETSC_COMM_WORLD,m_state,n_state,M_state,N_state,
  //                                    PETSC_NULL,&localMatDen);CHKERRQ(ierr);
  //   // TODO:  **NOTE** 
  //   // TODO:  (-M^-1 * K)^T = - K * M^-1  != -M^-1 * K (non commutative mat mult )
  //   // TODO:   ==>  F != F^T  
  //   // TODO:   CHECK NO SYMMETRY ASSUMPTIONS REMAIN ON SYSTEM DYNAMICS MATRIX
  //   // TODO:  **NOTE** 
  //   ierr=MatKSPSolveInverse_Basic(m_MassMatrix,stiffnessMatrix,localMatDen);
  //   ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD,m_state,n_state,M_state,N_state,
  //             n_state,PETSC_NULL,N_state-n_state,PETSC_NULL,&m_SystemDynamics);
  //   CHKERRQ(ierr);
  //   ierr=this->GetSparseMatrixFromDense(localMatDen,m_SystemDynamics);
  //   CHKERRQ(ierr);
  //   ierr = MatDestroy(localMatDen           );CHKERRQ(ierr); 
  //  }
 
  // clean up
  ierr = MatDestroy(stiffnessMatrix       );CHKERRQ(ierr); 

  // create identity matrix to be used as a template
  // ensure same parallel format at libmesh
  ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD,m_state,n_state,M_state,N_state,
                         1,PETSC_NULL,1,PETSC_NULL,&eyeState);
  ierr = MatZeroEntries(eyeState);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(eyeState,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(  eyeState,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatShift(eyeState,1.0);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(eyeState,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(  eyeState,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  // create id matrix in measurement space
  ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,
         numMeasurements,numMeasurements,1,PETSC_NULL,1,PETSC_NULL,&eyeMeasure);
  ierr = MatZeroEntries(eyeMeasure);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(eyeMeasure,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(  eyeMeasure,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatShift(eyeMeasure,1.0);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(eyeMeasure,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(  eyeMeasure,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  // initialize covariance data
  this->InitializeStateCov();

  // echo params
  ierr = PetscPrintf(PETSC_COMM_WORLD,
     "Factor Matrix SOLVER %s ORDER %s size (%dx%d)...\n",
                  solvertype,ordertype,M_state,N_state);CHKERRQ(ierr);

  // setup data structure for covariance update
  ierr=PetscPrintf(PETSC_COMM_WORLD,
            "Setup Work Data Struc for Covariance Update...\n");CHKERRQ(ierr);

  PetscLogEventEnd(  kfLogSetup,0,0,0,0); // setup

  // Software Guide : EndCodeSnippet
  PetscFunctionReturn(0);
}

/* setup measurement matrix to transform the state to the roi */
#undef __FUNCT__
#define __FUNCT__ "KalmanFilter::CreateIdentityMeasurementMap"
PetscInt KalmanFilter::CreateIdentityMeasurementMap(PetscInt MeasNodeSetID)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  typedef LITTSystem< PennesStandardDiffusionApproximation > LITTSDASystem;
  LITTSDASystem &state_system = m_eqnSystems->get_system< LITTSDASystem >("StateSystem");

  // image map reserve room for entire state
  std::vector<PetscInt> &imageMap =  state_system.m_NodeSets[state_system.u_var][MeasNodeSetID];
  if(imageMap.size())
    {
     std::cout << MeasNodeSetID << "nodeset already in use" << std::endl;
     std::cout <<  std::flush; libmesh_error();
    }

  // get temperature vector
  Vec templateVector;
  state_system.GetSubVector(templateVector, state_system.u_var);

  // state matrix dimensions to create id matrix
  PetscInt N_state;
  ierr = VecGetSize( templateVector ,&N_state);CHKERRQ(ierr);

  imageMap.resize(N_state,0);
  for (PetscInt Ii=0; Ii<N_state; Ii++) imageMap[Ii] = Ii;

  ierr = VecDestroy(templateVector);CHKERRQ(ierr);
  PetscFunctionReturn(imageMap.size());
}
/* setup measurement matrix to transform the state to the roi */
#undef __FUNCT__
#define __FUNCT__ "KalmanFilter::CreateMeasurementMapFromImaging"
PetscInt KalmanFilter::CreateMeasurementMapFromImaging(  char*  SystemName,
                                                        PetscInt MeasNodeSetID)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  // the # of measurements can is limited by the ROI size or the number of nodes
  libMesh::MeshBase &mesh = m_eqnSystems->get_mesh();

  // get dof map
  typedef LITTSystem< PennesStandardDiffusionApproximation > LITTSDASystem;
  LITTSDASystem &state_system = m_eqnSystems->get_system< LITTSDASystem >("StateSystem");
  const DofMap & dof_map = state_system.get_dof_map();

  // ROI Image mask
  System &ROISystem = m_eqnSystems->get_system(SystemName);

  // image map reserve room for entire state
  std::vector<PetscInt> &imageMap =  state_system.m_NodeSets[state_system.u_var][MeasNodeSetID];
  if(imageMap.size())
    {
     std::cout << MeasNodeSetID << "nodeset already in use" << std::endl;
     std::cout <<  std::flush; libmesh_error();
    }

  // get temperature vector
  Vec templateVector;
  state_system.GetSubVector(templateVector, state_system.u_var);

  // state matrix dimensions
  PetscInt N_state;
  ierr = VecGetSize( templateVector ,&N_state);CHKERRQ(ierr);
  imageMap.reserve(N_state);

  // clean up
  ierr = VecDestroy(templateVector);CHKERRQ(ierr);

  //  hold variable dofs
  std::vector<unsigned int> dof_indices_var;

  //  loop over elements and create map
  libMesh::MeshBase::const_element_iterator       el    =mesh.active_elements_begin();
  const libMesh::MeshBase::const_element_iterator end_el_full=mesh.active_elements_end();
  for ( ; el != end_el_full; ++el)
   {
    // Store element subdomain id to allow different sets of equations to be
    // solved on different parts of the mesh
    // ASSUMES SUBDOMAIN ORDERING STARTS FROM ONE but we need a 0 based
    // numbering scheme for std::vector
    Elem* elem = *el;
    const unsigned int subdomain_id = elem->subdomain_id() - 1;

    // only add healthy tissue domain
    if( subdomain_id == 0 )
      for( unsigned int i_var = 0 ; i_var < state_system.n_vars() ; i_var++)
       {
        // degree of freedom indices for individual components
        dof_map.dof_indices (elem,   dof_indices_var ,i_var);
        // nodal based setup of dirichlet data
        // indices should be global
        for(unsigned int Ii = 0 ; Ii < dof_indices_var.size() ; Ii++)
         if( ROISystem.current_solution( dof_indices_var[Ii] ) ) 
             imageMap.push_back( dof_indices_var[Ii] );  
       }
   } // end element loop
  // sort then erase duplicates
  std::sort( imageMap.begin(), imageMap.end() );
  std::vector<PetscInt>::iterator pos;

  pos = std::unique(imageMap.begin(),imageMap.end());
  imageMap.erase( pos,imageMap.end() ); 

  PetscFunctionReturn(imageMap.size());
}

/* setup measurement matrix to transform the state to the roi */
#undef __FUNCT__
#define __FUNCT__ "KalmanFilter::CreateROINodeSetMeasurementMap"
PetscErrorCode KalmanFilter::CreateROINodeSetMeasurementMap(
                                          PetscInt MeasNodeSetID)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  // get dof map
  typedef LITTSystem< PennesStandardDiffusionApproximation > LITTSDASystem;
  LITTSDASystem &state_system = m_eqnSystems->get_system< LITTSDASystem >("StateSystem");
  std::vector<PetscInt> &imageMap= state_system.m_NodeSets[state_system.u_var][MeasNodeSetID];

  // create identity matrix to be used as a template
  PetscInt numMeasurements = imageMap.size();
  if(numMeasurements <= 0)
    {
     std::cout << "map from state to measurements not set up " << std::endl;
     std::cout <<  std::flush; libmesh_error();
    }

  // create measurement matrix
  std::cout << "# Measurement Points =  " << numMeasurements << std::endl;
  PetscInt m_meas;
  PetscInt M_state,N_state,m_state,n_state;
  ierr = MatGetLocalSize(eyeMeasure,&m_meas,PETSC_NULL);CHKERRQ(ierr);
  ierr = MatGetSize(m_MassMatrix,&M_state,&N_state);CHKERRQ(ierr);
  ierr = MatGetLocalSize(m_MassMatrix,&m_state,&n_state);CHKERRQ(ierr);
  ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD,m_meas,n_state,
        numMeasurements,N_state,1,PETSC_NULL,1,PETSC_NULL,&MeasurementMat);
  CHKERRQ(ierr);

  // local ROI ordering implicit in the sorting
  for (PetscInt Ii=0; Ii<imageMap.size(); Ii++) 
    {
       MatSetValue(MeasurementMat,Ii, imageMap[Ii],1.0,INSERT_VALUES);
    }
  ierr = MatAssemblyBegin(MeasurementMat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(  MeasurementMat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  // error check
  PetscInt   rstart,rend,ncol;
  const PetscInt    *cwork;
  const PetscScalar *vwork;
  ierr = MatGetOwnershipRange(MeasurementMat,&rstart,&rend);CHKERRQ(ierr);
  for (PetscInt Ii=rstart; Ii<rend; Ii++) {
    ierr = MatGetRow(MeasurementMat,Ii,&ncol,&cwork,&vwork);CHKERRQ(ierr);
    if(ncol != 1)
      {
       std::cout << "unexpected entries in measurment matrix" << std::endl;
       std::cout << "ncol = " << ncol << std::endl;
       std::cout << "should we just average vwork??? " << std::endl;
       std::cout <<  std::flush; libmesh_error();
      }
    ierr = MatRestoreRow(MeasurementMat,Ii,&ncol,&cwork,&vwork);CHKERRQ(ierr);
  }

  // vectors for measurement based state update
  ierr = MatGetVecs(MeasurementMat,PETSC_NULL,&covMaxVec     ); CHKERRQ(ierr);
  ierr = MatGetVecs(MeasurementMat,PETSC_NULL,&measurementVec); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/* setup measurement matrix to transform the state to the roi 
   average over the axial dimension */
#undef __FUNCT__
#define __FUNCT__ "KalmanFilter::CreateROIAverageMeasurementMap"
PetscErrorCode KalmanFilter::CreateROIAverageMeasurementMap( 
                             int nx_roi_lower, int nx_roi_upper,
                             int ny_roi_lower, int ny_roi_upper,
                             int nz_roi_lower, int nz_roi_upper, int NSubDZ,
			     double X0, double Y0, double Z0,
			     double DX, double DY, double DZ)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  // the # of measurements can is limited by the ROI size or the number of nodes
  libMesh::MeshBase &mesh = m_eqnSystems->get_mesh();

  // get dof map
  PetscFEMSystem &state_system = m_eqnSystems->get_system<PetscFEMSystem>("StateSystem");
  const DofMap & dof_map = state_system.get_dof_map();
  const unsigned int u0_var = state_system.variable_number("u0");

  // create identity matrix
  const int nsize_roi[3] = {nx_roi_upper - nx_roi_lower +1, 
                            ny_roi_upper - ny_roi_lower +1, 
                            nz_roi_upper - nz_roi_lower +1}; 
  PetscInt numMeasurements = nsize_roi[0] * nsize_roi[1] * nsize_roi[2]; 

  // setup origin for ROI
  const double origin_roi[3] = {X0 + nx_roi_lower * DX,
                                Y0 + ny_roi_lower * DY,
                                Z0 + nz_roi_lower * DZ};

  // echo input 
  std::cout << "setting Measurement matrix for ROI " 
            << nx_roi_lower << " " << nx_roi_upper  << " "   
            << ny_roi_lower << " " << ny_roi_upper  << " "   
            << nz_roi_lower << " " << nz_roi_upper  << " "   
                            << " " << NSubDZ << std::endl;

  std::cout << "ROI Origin " 
            << origin_roi[0] << " "<< origin_roi[1] << " " << origin_roi[2] << std::endl;

  std::cout << "zplanes " ;  
  for (int kkk=0; kkk<nsize_roi[2]; kkk++)
    for (int ksub=0; ksub<NSubDZ; ksub++)
      {
       double zloc = static_cast<double>(kkk) 
                        + static_cast<double>(ksub)
                        / static_cast<double>(NSubDZ);
       std::cout << origin_roi[2]+zloc*DZ << " " ;  
      }
  std::cout << std::endl;

  // create measurement matrix
  std::cout << "# Measurement Points =  " << numMeasurements << std::endl;
  PetscInt m_meas;
  PetscInt M_state,N_state,m_state,n_state;
  ierr = MatGetLocalSize(eyeMeasure,&m_meas,PETSC_NULL);CHKERRQ(ierr);
  ierr = MatGetSize(m_MassMatrix,&M_state,&N_state);CHKERRQ(ierr);
  ierr = MatGetLocalSize(m_MassMatrix,&m_state,&n_state);CHKERRQ(ierr);
  ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD,m_meas,n_state,
        numMeasurements,N_state,20,PETSC_NULL,20,PETSC_NULL,&MeasurementMat);
  CHKERRQ(ierr);

  // create mesh function to get shape function values at arbitrary location
  MeshFunction state_values(*this->m_eqnSystems,
                            *state_system.current_local_solution, 
    		            dof_map, u0_var  );
  state_values.init();

  for (int kkk=0; kkk<nsize_roi[2]; kkk++)
    for (int jjj=0; jjj<nsize_roi[1]; jjj++)
      for (int iii=0; iii<nsize_roi[0]; iii++)
        for (int ksub=0; ksub<NSubDZ; ksub++)
          {
            int imageID = kkk*nsize_roi[0]*nsize_roi[1]+jjj*nsize_roi[0]+iii  ; 
            double zloc = static_cast<double>(kkk) 
                        + static_cast<double>(ksub)
                        / static_cast<double>(NSubDZ);
            const Point image_location(origin_roi[0]+iii*DX,
                                       origin_roi[1]+jjj*DY,
                                       origin_roi[2]+zloc*DZ);
            // locate the point in the other mesh
            const Elem* element = state_values.get_point_locator().operator()(image_location);

            if(element==NULL)
              {
                libmesh_error();
              }
            else
              {
                const unsigned int dim = mesh.mesh_dimension();
                
                /*
                 * Get local coordinates to feed these into compute_data().  
                 * Note that the fe_type can safely be used from the 0-variable,
                 * since the inverse mapping is the same for all FEFamilies
                 */
                const Point mapped_point (FEInterface::inverse_map (dim, 
                						    dof_map.variable_type(0),
                						    element, 
                						    image_location));
                
                /*
                 * the data for this variable
                 */
                const FEType& fe_type = dof_map.variable_type(u0_var);
                
                /**
                 * Build an FEComputeData that contains both input and output data
                 * for the specific compute_data method.
                 */
                {
                  FEComputeData data (*this->m_eqnSystems, mapped_point);
                  
                  FEInterface::compute_data (dim, fe_type, element, data);
                  
                  // where the solution values for the var-th variable are stored
                  std::vector<unsigned int> dof_indices;
                  dof_map.dof_indices (element, dof_indices, u0_var);
                  
                  // interpolate the solution
                  //{
                  //Number value = 0.;
                  //
                  //for (unsigned int i=0; i<dof_indices.size(); i++)
                  //  value += this->_vector(dof_indices[i]) * data.shape[i];
                  //
                  //output(index) = value;
                  //}
                  // Matrix entry should be shape functions normalized by number
                  // of subplanes to average over
                  for (unsigned int Ii=0; Ii<dof_indices.size(); Ii++)
                  MatSetValue(MeasurementMat,imageID,dof_indices[Ii],
                              data.shape[Ii]/static_cast<PetscScalar>(NSubDZ),
                                                               ADD_VALUES);
                  
                }
              }
          }

  // assemble
  ierr = MatAssemblyBegin(MeasurementMat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(  MeasurementMat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  // error check
  PetscInt   rstart,rend,ncol;
  const PetscInt    *cwork;
  const PetscScalar *vwork;
  ierr = MatGetOwnershipRange(MeasurementMat,&rstart,&rend);CHKERRQ(ierr);
  for (PetscInt Ii=rstart; Ii<rend; Ii++) {
    ierr = MatGetRow(MeasurementMat,Ii,&ncol,&cwork,&vwork);CHKERRQ(ierr);
    if(!ncol )
      {
       std::cout << "no entries in measurment matrix row "<< Ii << std::endl;
       std::cout <<  std::flush; libmesh_error();
      }
    ierr = MatRestoreRow(MeasurementMat,Ii,&ncol,&cwork,&vwork);CHKERRQ(ierr);
  }

  // vectors for measurement based state update
  ierr = MatGetVecs(MeasurementMat,PETSC_NULL,&covMaxVec     ); CHKERRQ(ierr);
  ierr = MatGetVecs(MeasurementMat,PETSC_NULL,&measurementVec); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

// Free Petsc Data structures
#undef __FUNCT__
#define __FUNCT__ "KalmanFilter::FinalizeDA"
PetscErrorCode KalmanFilter::FinalizeDA()
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  // destroy scatter context, vector, when no longer needed
  if(SumCov        !=NULL) {ierr = MatDestroy(SumCov       );CHKERRQ(ierr); }
  if(m_CovJacobian !=NULL) {ierr = MatDestroy(m_CovJacobian);CHKERRQ(ierr); }
  if(m_CovRHSMat   !=NULL) {ierr = MatDestroy(m_CovRHSMat  );CHKERRQ(ierr); }
  if(m_MassMatrix  !=NULL) {ierr = MatDestroy(m_MassMatrix );CHKERRQ(ierr); }
  if(tmpMatDenState!=NULL) {ierr = MatDestroy(tmpMatDenState);CHKERRQ(ierr); }

     
  ierr = MatDestroy(MeasurementMat);CHKERRQ(ierr);
  ierr = VecScatterDestroy(gather);CHKERRQ(ierr);
  ierr = VecDestroy(measurementVec);CHKERRQ(ierr);
  ierr = VecDestroy(covMaxVec);CHKERRQ(ierr);
  ierr = VecDestroy(imageVec);CHKERRQ(ierr);
  ierr = VecDestroy(globalVec);CHKERRQ(ierr); 

  PetscFunctionReturn(0);
}

// Destructor
#undef __FUNCT__
#define __FUNCT__ "KalmanFilter::~KalmanFilter"
KalmanFilter::~KalmanFilter()
{
  PetscErrorCode ierr = FinalizeDA(); 
}   

// Prediction of state 
#undef __FUNCT__
#define __FUNCT__ "KalmanFilter::StatePredict"
PetscErrorCode KalmanFilter::StatePredict(int istep)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  PetscLogEventBegin(kfLogStatePredict,0,0,0,0); // state predict

  /* assemble the load vector */
  // print info
  PetscFEMSystem &system = m_eqnSystems->get_system<PetscFEMSystem>("StateSystem");
  ierr = VecDataInfo(
        (dynamic_cast< PetscVector<double>* > (&(*system.solution)) )->vec() ,
                              "(Prediction) State_n-1");CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve linear system overwrite current state
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr=PetscPrintf(PETSC_COMM_WORLD,
                 "Prediction for state vector...\n");CHKERRQ(ierr);
  
  // This will put a local copy of solution into current_local_solution.
  system.solve();
  
  // also need to store the solution history for the adjoint
  ierr = StoreSystemTimeStep(&system, istep);
  ierr = VecDataInfo(
        (dynamic_cast< PetscVector<double>* > (&(*system.solution)) )->vec() ,
                                  "(Prediction) State_n");CHKERRQ(ierr);

  PetscLogEventEnd(  kfLogStatePredict,0,0,0,0); // state predict
  PetscFunctionReturn(0);
}
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   P(t^-_i) = P(t^+_{i-1})
     +
   (F P(t^+_{i-1}) + (F P(t^+_{i-1}) )^T + Q(t_{i-1})) \cdot \Delta t
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#undef __FUNCT__
#define __FUNCT__ "KalmanFilter::CovariancePredict"
PetscErrorCode KalmanFilter::
KalmanFilter::CovariancePredict()
{

  PetscErrorCode ierr;
  PetscFunctionBegin;

  PetscLogEventBegin(kfLogCovarPredict,0,0,0,0); // covariance predict

  ierr=PetscPrintf(PETSC_COMM_WORLD,
                   "Prediction for covariance...\n");CHKERRQ(ierr);

  // changed pointer to CovP to avoid reallocation during MatAXPY 
  ierr=PetscPrintf(PETSC_COMM_WORLD,"MatMatMult...\n");CHKERRQ(ierr);

  // build RHS of: m_CovJacobian p^i = m_CovRHSMat p^{i-1} + m_MassMatrix q
  // NOTE: SystemDynamics matrix is NOT SYMMETRIC
  Mat matrixRHSLoad,modelErrorLoad;
  // contribution from previous time step
  ierr = MatMatMult(m_CovRHSMat,CovP,MAT_INITIAL_MATRIX,
                    this->m_fill,&matrixRHSLoad); CHKERRQ(ierr);
  // modeling error source
  ierr = MatMatMult(m_MassMatrix,CovQ,MAT_INITIAL_MATRIX,
                    this->m_fill,&modelErrorLoad); CHKERRQ(ierr);
  // add the load
  ierr = MatAXPY(matrixRHSLoad,1.0,modelErrorLoad,
                            SUBSET_NONZERO_PATTERN);CHKERRQ(ierr);

  // update state covariance
  ierr = MatDataInfo(CovQ,"(Prediction) CovQ");CHKERRQ(ierr);

  // perform the matrix matrix solve
  // build RHS of: m_CovJacobian p^i = m_CovRHSMat p^{i-1} + m_MassMatrix q
  //ierr=MatKSPSolveInverse_Basic(m_CovJacobian,matrixRHSLoad,tmpMatDenState);
  ierr=this->MatMultipleRHSSolve(m_CovJacobian,matrixRHSLoad,tmpMatDenState);

  // Ensure symmetry 
  // CovP = 1/2 * (CovPHat + CovPHatTranspose)
  Mat CovPHatTranspose;
  ierr = MatTranspose(tmpMatDenState,MAT_INITIAL_MATRIX,&CovPHatTranspose);
  ierr = MatAXPY(tmpMatDenState,1.0,CovPHatTranspose,
                            SAME_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr= MatScale(tmpMatDenState,0.5);CHKERRQ(ierr);

  // copy the matrix 
  //  - in the case of a sparse localized matrix drop non zero's that violate the sparsity 
  //  - dense matrix should not lose any info
  ierr = MatCopy_Basic(tmpMatDenState,CovP,DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = MatDataSparseLossInfo(tmpMatDenState,CovP,"(Prediction) CovP");CHKERRQ(ierr);

  // apply the dirichlet data
  typedef LITTSystem< PennesStandardDiffusionApproximation > LITTSDASystem;
  LITTSDASystem &state_system = m_eqnSystems->get_system< LITTSDASystem >("StateSystem");
  std::map<PetscInt, std::vector<PetscInt> > &uvar_NodeSets = 
                                             state_system.m_NodeSets[state_system.u_var];
  if( uvar_NodeSets.size() )
   for( std::map<PetscInt,std::vector<PetscInt> >::iterator 
                    nsIt =uvar_NodeSets.begin(); 
                    nsIt!=uvar_NodeSets.end(); ++nsIt)
    {
      // store reference 
      std::vector<PetscInt> &nodeSetData = (*nsIt).second;

      ierr = MatZeroRows(CovP,nodeSetData.size(),
                             &nodeSetData[0],this->statecov);
      // transpose in place to zero columns 
      // should not need to transpose back due to symmetry
      ierr = MatTranspose(CovP,MAT_REUSE_MATRIX,&CovP);
      ierr = MatZeroRows(CovP,nodeSetData.size(),
                             &nodeSetData[0],this->statecov);

      //// state matrix dimensions
      //PetscInt M_state,N_state,m_state,n_state;
      //ierr = MatGetSize(CovP,&M_state,&N_state);CHKERRQ(ierr);
      //ierr = MatGetLocalSize(CovP,&m_state,&n_state);CHKERRQ(ierr);

      //// column mapping
      //std::vector<PetscInt> rows(M_state,0);
      //for( unsigned int Jj = 0; Jj< M_state;Jj++) rows[Jj] = Jj;

      //// dirichlet data uncorrelated with every thing else 
      //std::vector<PetscScalar> vals(M_state,0.0);
      //for( unsigned int Ii = 0; Ii < state_system.m_dirichletNodes.size();Ii++) 
      //  vals.at( state_system.m_dirichletNodes[Ii] ) = this->statecov;

      //// zero columns of RHS
      //for( unsigned int Ii = 0; Ii < state_system.m_dirichletNodes.size();Ii++) 
      //  {
      //   ierr = MatSetValues( matrixRHSLoad, M_state, &rows[0], 
      //                                             1, &state_system.m_dirichletNodes[Ii],
      //                                      &vals[0],INSERT_VALUES);
      //  }
      ierr = MatAssemblyBegin(CovP,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      ierr = MatAssemblyEnd(  CovP,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    }

  // free memory
  ierr = MatDestroy(matrixRHSLoad );CHKERRQ(ierr);
  ierr = MatDestroy(modelErrorLoad);CHKERRQ(ierr);
  ierr = MatDestroy(CovPHatTranspose);CHKERRQ(ierr);
  
  PetscLogEventEnd(  kfLogCovarPredict,0,0,0,0); // covariance predict
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "KalmanFilter::InitializeDenseUncorrCov"
PetscErrorCode KalmanFilter::
InitializeDenseUncorrCov(Mat &CovMatrix,Mat &TemplateMatrix,PetscScalar covariance)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  PetscInt M,N,m,n;
  ierr = MatGetSize(TemplateMatrix,&M,&N);CHKERRQ(ierr);
  ierr = MatGetLocalSize(TemplateMatrix,&m,&n);CHKERRQ(ierr);
  ierr = MatCreateMPIDense(PETSC_COMM_WORLD,m,n,M,N,
                           PETSC_NULL,&CovMatrix);CHKERRQ(ierr);

  /* Insert Ones and assemble to ensure full matrix is allocated */
  PetscInt Istart,Iend;
  ierr = MatGetOwnershipRange(CovMatrix,&Istart,&Iend);CHKERRQ(ierr);
  for (PetscInt Ii=Istart; Ii<Iend; Ii++) 
        MatSetValue(CovMatrix,Ii,Ii,covariance,INSERT_VALUES);
  ierr = MatAssemblyBegin(CovMatrix,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(  CovMatrix,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr=PetscPrintf(PETSC_COMM_WORLD,
                   "(Dense) Initial Covariance Formed...\n");CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "KalmanFilter::InitializeSparseUncorrCov"
PetscErrorCode KalmanFilter::
InitializeSparseUncorrCov(Mat &CovMatrix,Mat &TemplateMatrix,PetscScalar covariance)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr=MatDuplicate(TemplateMatrix,MAT_DO_NOT_COPY_VALUES,&CovMatrix);
  ierr=MatZeroEntries(CovMatrix);CHKERRQ(ierr);
  ierr=MatAssemblyBegin(CovMatrix,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr=MatAssemblyEnd(  CovMatrix,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr=MatShift(CovMatrix,covariance);CHKERRQ(ierr);

  /* Assemble matrices  */
  ierr=PetscPrintf(PETSC_COMM_WORLD,
            "(Sparse) Initial Covariance Formed...\n");CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
             disp('Compute factorization needed for Kalman gain factor');
             %K = P*inv(P+R);

             disp('State Correction based on observation');
              
             define (P+R) q = z-x ==>  Pq = z-x-Rq
                  x = x + K*(z-x) 
                    = x + P*inv(P+R) *(z-x)
                    = x + P*q
                    = z - R*q 
             
             disp('Covariance Correction based on observation');
              
             define (P+R) S = P ==> -P*S = R*S - P
                  P = P - K*P  
                    = P - P*inv(P+R)*P
                    = P - P*S
                    = P + R*S - P
                    = R*S 

  Perform state update AND setup part of the update for the covariance update
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#undef __FUNCT__
#define __FUNCT__ "KalmanFilter::StateUpdate"
PetscErrorCode KalmanFilter::StateUpdate(char *KalmanGainMethod)
{

  PetscErrorCode ierr;
  PetscFunctionBegin;

  PetscLogEventBegin(kfLogStateUpdate ,0,0,0,0); // state update

  // get state system
  typedef LITTSystem< PennesStandardDiffusionApproximation > LITTSDASystem;
  LITTSDASystem &state_system = m_eqnSystems->get_system< LITTSDASystem >("StateSystem");

  Vec State;
  state_system.GetSubVector(State, state_system.u_var);

  Vec Measurement = this->measurementVec;
  ierr=PetscPrintf(PETSC_COMM_WORLD,
                   "Update state from measurement data\n");CHKERRQ(ierr);
  
  // create vec to hold mat vec product
  Vec tmpVec;
  ierr = MatGetVecs(this->MeasurementMat,PETSC_NULL,&tmpVec); CHKERRQ(ierr);
  
  //tmpVec = Measurement - H * State
  ierr = MatMult(this->MeasurementMat,State,tmpVec);CHKERRQ(ierr);
  ierr = VecAYPX(tmpVec,-1.0,Measurement);CHKERRQ(ierr);

  // print info
  ierr = VecDataInfo(Measurement,"(Update) Measurement");CHKERRQ(ierr);
  ierr = VecDataInfo(tmpVec,"(Update) tmpVec");CHKERRQ(ierr);

  // compute the kalman gain
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Kalman Gain...\n");CHKERRQ(ierr);

  // SumCov = H P H^T + R
  Mat MeasMatTrans, PtimesHT;
  ierr = MatTranspose(MeasurementMat,MAT_INITIAL_MATRIX,&MeasMatTrans);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"(Update) MatMatMult1...\n");CHKERRQ(ierr);
  ierr = MatMatMult(this->CovP,MeasMatTrans,
                    MAT_INITIAL_MATRIX,1.0,&PtimesHT);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"(Update) MatMatMult2...\n");CHKERRQ(ierr);
  ierr = MatMatMult(this->MeasurementMat,PtimesHT,
                    MAT_INITIAL_MATRIX,1.0,&SumCov);
  //   clean up 
  ierr = MatDestroy(PtimesHT);CHKERRQ(ierr);
  ierr = MatDestroy(MeasMatTrans);CHKERRQ(ierr);

  const MatType storagetype;
  ierr = MatGetType(SumCov,&storagetype);CHKERRQ(ierr);

  // compute covariance sum
  ierr = MatAYPX(SumCov,1.0,CovR,SUBSET_NONZERO_PATTERN);CHKERRQ(ierr);

  ierr = MatDataInfo(SumCov,  "(Update) SumCov");CHKERRQ(ierr);

  ierr=PetscPrintf(PETSC_COMM_WORLD,"Update state... \n");CHKERRQ(ierr);

  // update tmpVec, compute the state update vector
  // SumCov* tmpVec = globalVec
  ierr = this->MeasurementSolve(KalmanGainMethod, tmpVec );

  // State= State + P * H^T *  tmpVec  
  ierr = MatMultTranspose(MeasurementMat,tmpVec,globalVec);CHKERRQ(ierr);
  Vec tmpState;
  ierr = VecDuplicate(State,&tmpState);
  ierr = MatMult(CovP,globalVec,tmpState);CHKERRQ(ierr);
  // don't update any dirichlet nodes
  std::vector<PetscInt> &dirichletNodeSet= state_system.m_NodeSets[state_system.u_var][1];
  if( dirichletNodeSet.size() )
    {
      // assume measurement is the same as dirichlet data
      for( unsigned int Ii = 0; Ii<dirichletNodeSet.size();Ii++) 
        ierr = VecSetValue (  tmpState, dirichletNodeSet[Ii],
                                         0.0, INSERT_VALUES);
      ierr = VecAssemblyBegin(tmpState);
      ierr = VecAssemblyEnd(  tmpState);
    }
  ierr = VecAXPY(State,1.0,tmpState);CHKERRQ(ierr);

  ierr = VecDataInfo(globalVec,"(Update) globalVec");CHKERRQ(ierr);
  ierr = VecDataInfo(tmpVec   ,"(Update) tmpVec");CHKERRQ(ierr);
  ierr = VecDataInfo(State    ,"(Update) State");CHKERRQ(ierr);

  // replace with the updated values
  state_system.SetSubVector(State, state_system.u_var);

  // clean up
  ierr = VecDestroy(tmpVec);
  ierr = VecDestroy(tmpState);
  ierr = VecDestroy(State);

  PetscLogEventEnd(  kfLogStateUpdate ,0,0,0,0); // state update

  PetscFunctionReturn(0);
}

// subroutine to extract image data to a petsc vec
#undef __FUNCT__
#define __FUNCT__ "KalmanFilter::CovarianceUpdate"
PetscErrorCode KalmanFilter::CovarianceUpdate()
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  PetscLogEventBegin(kfLogCovarUpdate ,0,0,0,0); // covariance update

  ierr=PetscPrintf(PETSC_COMM_WORLD,
          "Update covariance from measurement data\n");CHKERRQ(ierr);

  PetscLogEventBegin(kfLogInvSumCov ,0,0,0,0); // covariance update
  // easiest to do an mxm inverse on the stored
  //   covariance sum
  //
  //   InverseSumCov = [ H P H^T + R ]^-1
  //
  //     then matmat mult
  Mat InverseSumCov = NULL;
  this->InvertSumCov(InverseSumCov);
  PetscLogEventEnd(  kfLogInvSumCov ,0,0,0,0); // covariance update

  // clean up for StateUpdate
  ierr = MatDestroy(this->SumCov);CHKERRQ(ierr);

  // print info
  ierr = MatDataInfo(InverseSumCov,"(Cov Update) InverseSumCov");CHKERRQ(ierr);

  // compute projection
  Mat HtimesP, PtimesHT;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"MatMatMult... \n");CHKERRQ(ierr);
  ierr = MatMatMult(this->MeasurementMat,this->CovP,
                    MAT_INITIAL_MATRIX,1.0,&HtimesP);
  ierr = MatTranspose(HtimesP,MAT_INITIAL_MATRIX,&PtimesHT);

  // P  = P - PH^T InverseSumCov H P   ... (HP)^T = P H^T
  Mat tmpMattimesHP,ProjTmpMat;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"2xMatMatMult... \n");CHKERRQ(ierr);
  ierr = MatMatMult(InverseSumCov,HtimesP,
                    MAT_INITIAL_MATRIX,1.0,&tmpMattimesHP);
  ierr = MatMatMult(PtimesHT,tmpMattimesHP,
                    MAT_INITIAL_MATRIX,this->m_fill,&ProjTmpMat);

  // need  to ensure that SAME_NONZERO_PATTERN works for MatAYPX
  Mat CovPCopy;
  ierr=MatDuplicate(this->CovP,MAT_COPY_VALUES,&CovPCopy );
  // this will be the maximum sparsity of the covariance state
  ierr=MatSetOption(CovPCopy,MAT_NEW_NONZERO_LOCATIONS,
                    PETSC_FALSE); CHKERRQ(ierr);
  CHKERRQ(ierr);
  // copy the matrix and drop non zero's that violate the sparsity
  // assumption (SAME_NONZERO_PATTERN shouldn't matter)
  ierr = MatCopy_Basic(ProjTmpMat,this->CovP,SAME_NONZERO_PATTERN);
  CHKERRQ(ierr);

  // P  = P - PH^T InverseSumCov H P   ... (HP)^T = P H^T
  ierr = MatAYPX(this->CovP,-1.0,CovPCopy ,
                 SAME_NONZERO_PATTERN); CHKERRQ(ierr);

  // print info
  ierr = MatDataInfo(ProjTmpMat,"(Sparse Update) ProjTmpMat");CHKERRQ(ierr);
  ierr = MatDataInfo(CovPCopy  ,"(Sparse Update) CovPCopy  ");CHKERRQ(ierr);
  ierr = MatDataInfo(this->CovP,"(Sparse Update) this->CovP");CHKERRQ(ierr);

  // free memory
  ierr = MatDestroy( InverseSumCov );CHKERRQ(ierr);
  ierr = MatDestroy( HtimesP );CHKERRQ(ierr);
  ierr = MatDestroy( tmpMattimesHP);CHKERRQ(ierr);
  ierr = MatDestroy( PtimesHT);CHKERRQ(ierr);
  ierr = MatDestroy( ProjTmpMat);CHKERRQ(ierr);
  ierr = MatDestroy( CovPCopy);CHKERRQ(ierr);

  PetscLogEventEnd(  kfLogCovarUpdate ,0,0,0,0); // covariance update

  PetscFunctionReturn(0);
}

// extract variance for plotting
#undef __FUNCT__
#define __FUNCT__ "KalmanFilter::ExtractVarianceForPlotting"
PetscErrorCode KalmanFilter::ExtractVarianceForPlotting(char*  SystemName)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  System &uncertainty_system= m_eqnSystems->get_system(SystemName);
  ierr = MatGetDiagonal(this->CovP,
   (dynamic_cast< PetscVector<double>* > (&(*uncertainty_system.solution)) )->vec() );
  uncertainty_system.solution->localize(
                   *uncertainty_system.current_local_solution);
  PetscFunctionReturn(0);
}

// extract a covariance column for plotting
#undef __FUNCT__
#define __FUNCT__ "KalmanFilter::ExtractCoVarianceForPlotting"
PetscErrorCode KalmanFilter::ExtractCoVarianceForPlotting(char* SystemName,
                                                          PetscInt ColumnID)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  System &uncertainty_system= m_eqnSystems->get_system(SystemName);

  // create vector for multiplication
  Vec  canonicalBasis;
  ierr = MatGetVecs(this->CovP,
                           &canonicalBasis,PETSC_NULL); CHKERRQ(ierr);
  ierr = VecZeroEntries(    canonicalBasis);
  ierr = VecSetValue(       canonicalBasis,ColumnID,1.0,INSERT_VALUES);
  ierr = VecAssemblyBegin(  canonicalBasis); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(    canonicalBasis); CHKERRQ(ierr);
  ierr = MatMult(this->CovP,canonicalBasis,
   (dynamic_cast< PetscVector<double>* > (&(*uncertainty_system.solution)) )->vec() );
         CHKERRQ(ierr);

  uncertainty_system.solution->localize(
                   *uncertainty_system.current_local_solution);
  PetscFunctionReturn(0);
}
// subroutine to extract image data to a petsc vec
#undef __FUNCT__
#define __FUNCT__ "KalmanFilter::ExtractMeasurementData"
PetscErrorCode KalmanFilter::
ExtractMeasurementData(Vec MRTIData, char* UncertaintySystemName)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  // copy input imaging data to measurement vector
  ierr = VecCopy( MRTIData,this->measurementVec);

  // Get a reference to the ExplicitSystem for measured data
  System & ideal_uncertainty_system =
     m_eqnSystems->get_system(UncertaintySystemName);

  // extract uncertainty data into mx1 vector
  PetscTruth  overrideMax=PETSC_FALSE;
  ierr = PetscOptionsGetTruth(PETSC_NULL,"-override_max",&overrideMax,PETSC_NULL); 
  CHKERRQ(ierr);

  // extract uncertainty data into mx1 vector
  Vec tmpMeasure;
  ierr = VecDuplicate(this->measurementVec,&tmpMeasure);
  ierr = MatMult(this->MeasurementMat,
   (dynamic_cast< PetscVector<double>* > 
                 (&(*ideal_uncertainty_system.solution)))->vec(), tmpMeasure);CHKERRQ(ierr);

  if(overrideMax)
    ierr = VecCopy( tmpMeasure ,covMaxVec);
  else // max to this time point
    ierr = VecPointwiseMax(covMaxVec, tmpMeasure , covMaxVec);

  ierr = VecDestroy(tmpMeasure);

  //if(this->Images) this->Images->WriteFEMImage(
  //                   ideal_uncertainty_system,istep,0.0,"tmapstd");

  ierr = MatDiagonalSet(CovR,covMaxVec,INSERT_VALUES);
  ierr = MatDataInfo(CovR,"(Extract Covariance) CovR");CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "KalmanFilter::GetSparseMatrixFromDense"
PetscErrorCode KalmanFilter::
GetSparseMatrixFromDense(Mat &DenseMat, Mat &SparseMat)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  // store nonzero structure in AIJ format for later
 
  /* 
  PetscInt       m,N;
  Vec              x;
  PetscScalar    *xx,*data;
  ierr = MatGetArray(DenseMat,&xx);CHKERRQ(ierr);
  ierr = MatGetLocalSize(DenseMat,&m,PETSC_NULL);CHKERRQ(ierr);  // number local rows 
  ierr = MatGetSize(DenseMat,PETSC_NULL,&N);CHKERRQ(ierr);       // total columns in dense matrix
  ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,m,
                               PETSC_DETERMINE,PETSC_NULL,&x);CHKERRQ(ierr);
  for (PetscInt Jj=0; Jj<N; Jj++) {
    ierr = VecPlaceArray(x,xx + Jj*m);CHKERRQ(ierr);
    ierr = VecGetArray(   x,&data);CHKERRQ(ierr);
    for (PetscInt Ii=0; Ii<m; Ii++) {
      if( std::abs(data[Ii]) > zerotol ) { 
        MatSetValue(SparseMat     ,Ii,Jj,data[Ii],INSERT_VALUES);
      }
    }
    ierr = VecRestoreArray(x,&data);CHKERRQ(ierr);
    ierr = VecResetArray(x);CHKERRQ(ierr);
  }
  ierr = VecDestroy(x);CHKERRQ(ierr);
  ierr = MatRestoreArray(DenseMat,&xx);CHKERRQ(ierr);
  */

  const PetscScalar  *vwork;
  const PetscInt     *cwork;
  PetscInt Istart,Iend, nz;

  ierr=PetscPrintf(PETSC_COMM_WORLD,
                   "  Converting...\n");CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(DenseMat,&Istart,&Iend);CHKERRQ(ierr);
  for (PetscInt Ii=Istart; Ii<Iend; Ii++) {
    ierr = MatGetRow(DenseMat,Ii,&nz,&cwork,&vwork);CHKERRQ(ierr);
    for (PetscInt Jj=0; Jj<nz; Jj++) {
      //ierr = MatSetValues(SparseMat,1,&Ii,nz,cwork,vwork,INSERT_VALUES);CHKERRQ(ierr);
      if( std::abs(vwork[Jj]) > m_zerotol ) { 
        MatSetValue(SparseMat ,Ii,cwork[Jj],vwork[Jj],INSERT_VALUES);
      }
    }
    ierr = MatRestoreRow(DenseMat,Ii,&nz,&cwork,&vwork);CHKERRQ(ierr);
  }

  ierr=PetscPrintf(PETSC_COMM_WORLD,
                   "  Assembling Sparse From Dense...\n");CHKERRQ(ierr);
  ierr = MatAssemblyBegin(SparseMat ,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(  SparseMat ,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
// Constructor
#undef __FUNCT__
#define __FUNCT__ "SparseKalmanFilter::SparseKalmanFilter"
SparseKalmanFilter::SparseKalmanFilter(EquationSystems *EqnSystemsPtr):
                    KalmanFilter(EqnSystemsPtr) 
{
  // linear algebra control params
  this->m_maxSpread = 3 ;
  if( this->m_maxSpread < 0 )
    {
     std::cout << "m_maxSpread >= 0 " 
               << std::endl << std::flush; abort();
    }
}
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "SparseKalmanFilter::InitializeStateCov"
/* 
   SparseInitializeState - Forms initial state and covariance

   Input Parameters:
   StateTemp - vector

 */
PetscErrorCode SparseKalmanFilter::InitializeStateCov()
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  // echo data
  ierr=PetscPrintf(PETSC_COMM_WORLD,"statecov=%f modelcov=%f \n",
                                  this->statecov,this->modelcov); CHKERRQ(ierr);

  // initialize covariance for measurements
  PetscScalar dummy = 1.0; // just need to initialize... will be overwritten
  ierr = this->InitializeSparseUncorrCov(this->CovR,this->eyeMeasure,dummy); 
  CHKERRQ(ierr);
  ierr = MatDataInfo(this->CovR,"Init CovR");CHKERRQ(ierr);

  // initialize covariance for model
  ierr = this->InitializeSparseUncorrCov(this->CovQ,this->eyeState,this->modelcov); 
  CHKERRQ(ierr);
  ierr = MatDataInfo(this->CovQ,"Init CovQ");CHKERRQ(ierr);

  // initialize covariance for state
  // VERY important for performance to limit the sparsity
  Mat CovPTemplate;
  ierr=MatDuplicate(this->eyeState,MAT_COPY_VALUES,&CovPTemplate);

  // control the maximum sparsity of the state covariance matrix
  for( PetscInt ispread=0 ; ispread < this->m_maxSpread; ispread++ )
   {
    ierr = MatMatMult(this->m_MassMatrix,CovPTemplate,
                      MAT_INITIAL_MATRIX,1.0,&this->CovP);
    // clean up and reset
    ierr = MatDestroy(CovPTemplate);CHKERRQ(ierr);
    CovPTemplate = this->CovP;
   }

  ierr = this->InitializeSparseUncorrCov(this->CovP,CovPTemplate,
                                     this->statecov); CHKERRQ(ierr);

  // this will be the maximum sparsity of the covariance state
  ierr=MatSetOption(this->CovP,MAT_NEW_NONZERO_LOCATIONS,
                    PETSC_FALSE); CHKERRQ(ierr);

  ierr=PetscPrintf(PETSC_COMM_WORLD,
                 "\n...max spread = %d\n", this->m_maxSpread);
  ierr = MatDataInfo(this->CovP,"Init CovP");CHKERRQ(ierr);

  //// simple model covariance update for regression testing
  //if ( user->get_num_exact () && !this->noDenseUpdate )
  //  {
  //    // assume system dynamics matrix is symmetric
  //    ierr = MatAXPY(this->CovQ,-2.0,this->m_MassMatrix,
  //                                  DIFFERENT_NONZERO_PATTERN);
  //    ierr = MatDataInfo(this->CovQ,"Init CovQ");CHKERRQ(ierr);
  //  }

  PetscFunctionReturn(0);
} 

// Constructor
#undef __FUNCT__
#define __FUNCT__ "DenseKalmanFilter::DenseKalmanFilter"
DenseKalmanFilter::DenseKalmanFilter(EquationSystems *EqnSystemsPtr):
                        KalmanFilter(EqnSystemsPtr) 
{
}
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "DenseKalmanFilter::InitializeStateCov"
/* InitializeStateCov - Forms initial covariance */
PetscErrorCode DenseKalmanFilter::InitializeStateCov()
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  // echo data
  ierr=PetscPrintf(PETSC_COMM_WORLD,"statecov=%f modelcov=%f \n",
                                     this->statecov,this->modelcov); CHKERRQ(ierr);

  // initialize covariance for measurements
  PetscScalar dummy = 1.0; // just need to initialize... will be overwritten
  ierr = InitializeDenseUncorrCov(this->CovR,this->eyeMeasure,dummy); CHKERRQ(ierr);
  ierr = MatDataInfo(this->CovR,"Init CovR");CHKERRQ(ierr);

  // initialize covariance for model
  ierr = InitializeDenseUncorrCov(this->CovQ,this->eyeState, this->modelcov); CHKERRQ(ierr);
  ierr = MatDataInfo(this->CovQ,"Init CovQ");CHKERRQ(ierr);

  // initialize covariance for state
  ierr = InitializeDenseUncorrCov(this->CovP,this->eyeState, this->statecov);CHKERRQ(ierr);
  ierr = MatDataInfo(this->CovP,"Init CovP");CHKERRQ(ierr);

  PetscFunctionReturn(0);
} 
// Constructor
#undef __FUNCT__
#define __FUNCT__ "UncorrelatedKalmanFilter::UncorrelatedKalmanFilter"
UncorrelatedKalmanFilter::UncorrelatedKalmanFilter(EquationSystems *EqnSystemsPtr):
                        DenseKalmanFilter(EqnSystemsPtr) 
{ 
}
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                  s.P     =   s.A * s.P * s.A' + s.Q
                      \approx       s.P        + s.Q
   This routine approximates the state transition matrix as the 
   identity, s.A \approx I, to avoid any complexities associated 
   with the dense linear algebra.
   TODO: determine if this is a valid approximation, other than 
         just making the code run fast...
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#undef __FUNCT__
#define __FUNCT__ "UncorrelatedKalmanFilter::CovariancePredict"
PetscErrorCode UncorrelatedKalmanFilter::CovariancePredict()
{

  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr=PetscPrintf(PETSC_COMM_WORLD,
                   "Prediction for covariance...\n");CHKERRQ(ierr);
  ierr = MatDataInfo(this->CovQ,"(Prediction) CovQ");CHKERRQ(ierr);
  PetscFEMSystem &state_system = m_eqnSystems->get_system<PetscFEMSystem>("StateSystem");
  ierr = MatAXPY(this->CovP,state_system.deltat,this->CovQ,
                 SUBSET_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = MatDataInfo(this->CovP,"(Prediction) CovP");CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

// compute the state update
#undef __FUNCT__
#define __FUNCT__ "DenseKalmanFilter::MeasurementSolve"
PetscErrorCode DenseKalmanFilter::MeasurementSolve( char *KalmanGainMethod, Vec UpdateVec )
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  // params for matrix factorization
  MatFactorInfo factorInfo;
  ierr = MatFactorInfoInitialize(&factorInfo); CHKERRQ(ierr);
  Mat PreFactSumCov = NULL;

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Converting... \n");CHKERRQ(ierr);
  if( !strcasecmp(KalmanGainMethod,"superlu_dist") ) 
    {// superlu only works with aij matrix have to convert first
    PetscInt M_meas,N_meas,m_meas,n_meas;
    ierr = MatGetSize(eyeMeasure,&M_meas,&N_meas);CHKERRQ(ierr);
    ierr = MatGetLocalSize(eyeMeasure,&m_meas,&n_meas);CHKERRQ(ierr);
    ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD,m_meas,n_meas,M_meas,N_meas,n_meas,
                           PETSC_NULL,N_meas-n_meas,PETSC_NULL,&PreFactSumCov );CHKERRQ(ierr);
    ierr = MatCopy(this->SumCov,PreFactSumCov,DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
    }
  else // create and use a copy to be consistent with possible AIJ format above
    ierr = MatConvert(this->SumCov,MATSAME,MAT_INITIAL_MATRIX,&PreFactSumCov );

  ierr = MatGetFactor(PreFactSumCov ,KalmanGainMethod, 
                      MAT_FACTOR_LU,&SumCovFact); 
  IS isrow,iscol;
  ierr = MatGetOrdering(PreFactSumCov ,this->ordertype,&isrow,&iscol); CHKERRQ(ierr);

  // factor 
  ierr = PetscPrintf(PETSC_COMM_WORLD,"MatFactor... \n");CHKERRQ(ierr);
  MatLUFactorSymbolic(SumCovFact,PreFactSumCov ,isrow,iscol,&factorInfo);
  ierr = MatLUFactorNumeric(SumCovFact,PreFactSumCov ,&factorInfo);

  // solve
  // Update vec should be pass in to this routine as the rhs
  // MatSolve cannot overwrite so we must duplicate and solve
  Vec MeasRHS;
  ierr = VecDuplicate(UpdateVec,&MeasRHS);CHKERRQ(ierr);
  ierr = VecCopy(UpdateVec,MeasRHS);CHKERRQ(ierr);
  ierr = MatSolve(SumCovFact,MeasRHS,UpdateVec);CHKERRQ(ierr);

  // clean up 
  ierr = MatDestroy(PreFactSumCov);
  ierr = VecDestroy(MeasRHS);
  ierr = ISDestroy(isrow);
  ierr = ISDestroy(iscol);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DenseKalmanFilter::InvertSumCov"
PetscErrorCode DenseKalmanFilter::InvertSumCov( Mat &CovarianceSumInverse )
{

  PetscErrorCode ierr;
  PetscFunctionBegin;

  // get measurement size
  PetscInt M_meas,N_meas,m_meas,n_meas;
  ierr = MatGetSize(eyeMeasure,&M_meas,&N_meas);CHKERRQ(ierr);
  ierr = MatGetLocalSize(eyeMeasure,&m_meas,&n_meas);CHKERRQ(ierr);

  // scratch matrix to hold solution for covariance update
  ierr = MatCreateMPIDense(PETSC_COMM_WORLD,m_meas,n_meas,M_meas,N_meas,
                           PETSC_NULL,&CovarianceSumInverse );CHKERRQ(ierr);

  // create dense identity to invert from
  Mat DenseEye;
  ierr = MatCreateMPIDense(PETSC_COMM_WORLD,m_meas,n_meas,M_meas,N_meas,
                           PETSC_NULL,&DenseEye);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(DenseEye,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(  DenseEye,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr= MatShift(DenseEye,1.0);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(DenseEye,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(  DenseEye,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr=PetscPrintf(PETSC_COMM_WORLD,"MatMatSolve... \n");CHKERRQ(ierr);
  ierr = MatMatSolve(SumCovFact,DenseEye,CovarianceSumInverse);CHKERRQ(ierr);
  PetscTruth  assembled ; 
  ierr = MatAssembled(CovarianceSumInverse,&assembled); CHKERRQ(ierr);
  if(!assembled)
    { 
     ierr = MatAssemblyBegin(CovarianceSumInverse,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
     ierr = MatAssemblyEnd(CovarianceSumInverse,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    }

  // free memory
  ierr = MatDestroy(SumCovFact);CHKERRQ(ierr);
  ierr = MatDestroy(DenseEye);CHKERRQ(ierr);
 
  // write matrix to a file to examine in matlab
  // write matlab file to .mat file to read in numpy
  PetscTruth  verifyInverse=PETSC_FALSE;
  ierr = PetscOptionsGetTruth(PETSC_NULL,"-verify_inverse",&verifyInverse,PETSC_NULL); 
  if(verifyInverse)
    { 
     PetscViewer matViewer; 
     // write inverse
     ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"SumCovInv.m",
                                 &matViewer);  CHKERRQ(ierr);
     ierr = PetscViewerSetFormat(matViewer,
                   PETSC_VIEWER_ASCII_MATLAB); CHKERRQ(ierr);
     ierr = MatView(CovarianceSumInverse,matViewer); CHKERRQ(ierr);
     ierr = PetscViewerFlush(matViewer); CHKERRQ(ierr);
     ierr = PetscViewerDestroy(matViewer); CHKERRQ(ierr);
     PetscPrintf(PETSC_COMM_WORLD,"wrote SumCovInv.m ...\n" );

     // write original
     ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"SumCov.m",
                                 &matViewer);  CHKERRQ(ierr);
     ierr = PetscViewerSetFormat(matViewer,
                   PETSC_VIEWER_ASCII_MATLAB); CHKERRQ(ierr);
     ierr = MatView(this->SumCov,matViewer); CHKERRQ(ierr);
     ierr = PetscViewerFlush(matViewer); CHKERRQ(ierr);
     ierr = PetscViewerDestroy(matViewer); CHKERRQ(ierr);
     PetscPrintf(PETSC_COMM_WORLD,"wrote SumCov.m ...\n" );

     // write initial matrix in binary to read in
     //ierr = PetscViewerMatlabOpen(PETSC_COMM_WORLD,
     //                             "SumCov.mat",
     //                             FILE_MODE_WRITE,&matViewer);CHKERRQ(ierr);
     //ierr = MatView(this->SumCov    ,matViewer);CHKERRQ(ierr);
     //ierr = PetscViewerDestroy(matViewer);CHKERRQ(ierr);
    }

  PetscFunctionReturn(0);
}

// compute the state update
#undef __FUNCT__
#define __FUNCT__ "SparseKalmanFilter::MeasurementSolve"
PetscErrorCode SparseKalmanFilter::MeasurementSolve( char *, Vec UpdateVec )
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  KSP            kspStateUpdate;         /* Krylov subspace method context */
  ierr = KSPCreate(PETSC_COMM_WORLD,&kspStateUpdate);CHKERRQ(ierr);
  ierr = KSPSetOperators(kspStateUpdate,SumCov,SumCov,
         // DIFFERENT_NONZERO_PATTERN is
         //  irrelavent for this ONE solve KSPDestroy is called
         DIFFERENT_NONZERO_PATTERN);
  CHKERRQ(ierr);
  // set defaults
  PC pcstate;
  ierr = KSPGetPC(kspStateUpdate,&pcstate); CHKERRQ(ierr);
  ierr = PCSetType(pcstate,PCBJACOBI); CHKERRQ(ierr);
  PetscInt maxiter;
  ierr = MatGetSize(this->eyeState,&maxiter,PETSC_NULL);CHKERRQ(ierr);
  ierr = KSPSetTolerances(kspStateUpdate,rtol0,PETSC_DEFAULT,PETSC_DEFAULT,
                          maxiter); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(kspStateUpdate);CHKERRQ(ierr);
  ierr = KSPMonitorCancel(kspStateUpdate);CHKERRQ(ierr);
  // KSPSolve overwrites rhs w/ the solution
  ierr = KSPSolve(kspStateUpdate,UpdateVec,UpdateVec);CHKERRQ(ierr);
  //check convergence reason
  KSPConvergedReason kspReason;
  KSPGetConvergedReason(kspStateUpdate, &kspReason);
  if(kspReason<0)
    {
     std::cerr<<"StateUpdate Diverged "<< kspReason
              << std::endl <<std::flush; libmesh_error(); 
    } 
  else
    { 
     PetscInt kspIterNum;
     KSPGetIterationNumber(kspStateUpdate, &kspIterNum);
     std::cout<<"StateUpdate Converged "<< kspReason
              <<" iterations "<< kspIterNum
              << std::endl <<std::flush;
    } 
  ierr = KSPDestroy(kspStateUpdate);CHKERRQ(ierr); 

  PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "SparseKalmanFilter::InvertSumCov"
PetscErrorCode SparseKalmanFilter::InvertSumCov( Mat &CovarianceSumInverse )
{

  PetscErrorCode ierr;
  PetscFunctionBegin;

  // get measurement size
  PetscInt M_meas,N_meas,m_meas,n_meas;
  ierr = MatGetSize(eyeMeasure,&M_meas,&N_meas);CHKERRQ(ierr);
  ierr = MatGetLocalSize(eyeMeasure,&m_meas,&n_meas);CHKERRQ(ierr);

  // scratch matrix to hold solution for covariance update
  Mat tmpMatDenMeas;
  ierr = MatCreateMPIDense(PETSC_COMM_WORLD,m_meas,n_meas,M_meas,N_meas,
                                   PETSC_NULL,&tmpMatDenMeas);CHKERRQ(ierr);

  ierr=PetscPrintf(PETSC_COMM_WORLD,"MatKSPSolveInverse... \n");CHKERRQ(ierr);
  ierr=MatKSPSolveInverse_Basic(this->SumCov,NULL,tmpMatDenMeas);
  CHKERRQ(ierr);
 
  // get sparse matrix from dense matrix
  ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD,m_meas,n_meas,M_meas,N_meas,n_meas,
                         PETSC_NULL,N_meas-n_meas,PETSC_NULL,&CovarianceSumInverse);CHKERRQ(ierr);
  ierr = this->GetSparseMatrixFromDense(tmpMatDenMeas,CovarianceSumInverse); CHKERRQ(ierr);
 
  // print info
  ierr = MatDataInfo(tmpMatDenMeas,"(Sparse Update) tmpMatDenMeas");CHKERRQ(ierr);
  ierr = MatDestroy(tmpMatDenMeas);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

//#undef __FUNCT__
//#define __FUNCT__ "UncorrelatedKalmanFilter::InvertSumCov"
//PetscErrorCode UncorrelatedKalmanFilter::InvertSumCov( Mat &CovarianceSumInverse )
//{
//  PetscErrorCode ierr;
//  PetscFunctionBegin;
//  
//  /**
//   * update covariance from measurement data this assumes that upon entering
//   * this routine the state covariance matrix is a diagonal matrix. The nice
//   * linear algebra properties of the diagonal matrix are exploited to compute
//   * the kalman gain and update the state covariance.  
//   */
//  this->InitializeSparseUncorrCov(CovarianceSumInverse,this->eyeMeasure,1.0);
//
//  // inverse should be diagonal of the reciprocal
//  Vec diagSumCov;
//  ierr = MatGetVecs(this->MeasurementMat,PETSC_NULL,&diagSumCov); 
//  CHKERRQ(ierr);
//  ierr = MatGetDiagonal(this->SumCov,diagSumCov); CHKERRQ(ierr);
//  ierr = VecReciprocal(diagSumCov); CHKERRQ(ierr);
//
//  // restore to matrix
//  ierr = MatZeroEntries(  CovarianceSumInverse);CHKERRQ(ierr);
//  ierr = MatAssemblyBegin(CovarianceSumInverse,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//  ierr = MatAssemblyEnd(  CovarianceSumInverse,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//  ierr = MatDiagonalSet(  CovarianceSumInverse,diagSumCov,INSERT_VALUES);
//
//  // SumCov = H P H^T + R should still be diagonal ? 
//  // error check that we have the correct inverse
//  Mat locIdentity;
//  ierr= MatMatMult(this->SumCov,CovarianceSumInverse,MAT_INITIAL_MATRIX,
//                    1.0,&locIdentity); CHKERRQ(ierr);
//  ierr= MatShift(locIdentity,-1.0);CHKERRQ(ierr);
//  PetscScalar norm_one, norm_inf ;
//  ierr = MatNorm(locIdentity,   NORM_1    ,&norm_one);CHKERRQ(ierr);
//  ierr = MatNorm(locIdentity,NORM_INFINITY,&norm_inf);CHKERRQ(ierr);
//  if ( std::max(norm_one,norm_inf) > 1.0e-6)
//    {
//      std::cout << "incorrect inverse..." << std::endl << std::flush;
//      libmesh_error(); 
//    }
//  // clean up
//  ierr = VecDestroy(diagSumCov);CHKERRQ(ierr);
//  ierr = MatDestroy(locIdentity);CHKERRQ(ierr); 
//  
//  PetscFunctionReturn(0);
//}

#undef __FUNCT__
#define __FUNCT__ "KalmanFilter::ProjectMeasurementMatrix"
PetscErrorCode KalmanFilter::ProjectMeasurementMatrix(Vec Soln ,
                                                      Vec *VecData)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = MatGetVecs(this->MeasurementMat,PETSC_NULL,VecData);CHKERRQ(ierr);
  ierr = MatMult(MeasurementMat,Soln,*VecData);CHKERRQ(ierr);

  // don't destroy VecData let python destroy
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "KalmanFilter::MeasurementMatrixTranspose"
PetscErrorCode KalmanFilter::MeasurementMatrixTranspose(char*  SystemName,
                                                   PetscScalar ROIValue)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  System & system = m_eqnSystems->get_system(SystemName);
  Vec Soln = (dynamic_cast< PetscVector<double>* > (&(*system.solution)) )->vec();
  Vec tmpVec;
  ierr = MatGetVecs(this->MeasurementMat,PETSC_NULL,&tmpVec); CHKERRQ(ierr);
  ierr = VecSet(tmpVec,ROIValue);CHKERRQ(ierr);
  ierr = MatMultTranspose(MeasurementMat,tmpVec,Soln);CHKERRQ(ierr);

  // clean up
  ierr = VecDestroy(tmpVec);

  PetscFunctionReturn(0);
}
