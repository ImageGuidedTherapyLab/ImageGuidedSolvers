// C++ include files that we need
#include <iostream>
#include <vector>

// libmesh
#include "libmesh.h"
#include "equation_systems.h"
#include "fe.h"
#include "fe_interface.h"
#include "quadrature_gauss.h"
#include "dof_map.h"
#include "sparse_matrix.h"
#include "elem.h"
#include "string_to_enum.h"
#include "getpot.h"
#include "parallel.h"
#include "boundary_info.h"

// tao interface
#include "src/tao_impl.h" 
#include "src/petsctao/vector/taovec_petsc.h" 
#include "src/bound/impls/blmvm/blmvm.h" 

// The nonlinear solver and system we will be using
#include "nonlinear_implicit_system.h"
#include "linear_implicit_system.h"
#include "petsc_linear_solver.h"

// pdeopt
#include "applicationContext.h" // basic information
#include "pdeBaseClass.h"
#include "pennesInverseSystem.h"
#include "quantityOfInterest.h"
#include "dddas.h"
#include "variable_map.h"

#include "mesh.h"
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "PhysBnds"
/* 
   PhysBnds -  Returns the physically acceptable bounds of the model parameters

   Input Parameter:
   user - user-defined application context

   Output Parameter:
   XL,XU array of lower,upper bounds on model parameters
*/
PetscErrorCode PhysBnds(TAO_APPLICATION , Vec XL,Vec XU,void *ctx)
{
  // used for error handling
  PetscFunctionBegin; 

  AppSolve          *user = (AppSolve *) ctx;  
  libMesh::MeshBase          &mesh = user->get_mesh();
  qoiBaseClass*      qoiOptimizer = user->qoiOptimizer();
  PetscErrorCode info = qoiOptimizer->physicalBounds( mesh,XL,XU);
  CHKERRQ(info);

  // plot parameter bounds to double check
  OStringStream bndsFileName;

  // lower bounds
  bndsFileName<< "femvis/femqoi_" << user->GroupID
              << "lb_nopt_"       << user->IDOPT 
              << AppSolve::profileID  << ".e";
  qoiOptimizer->plotElemData(bndsFileName,mesh,XL); 

  // upper bounds
  bndsFileName.str(""); // reset before reuse
  bndsFileName<< "femvis/femqoi_" << user->GroupID
              << "ub_nopt_"       << user->IDOPT 
              << AppSolve::profileID  << ".e";
  qoiOptimizer->plotElemData(bndsFileName,mesh,XU); 

  PetscFunctionReturn(0);
}
/*  
    FormObjectiveAndGradient - Evaluates the function, f(X), and its gradient
    Input Parameters:
.   taoapp  - the TAO_APPLICATION context
.   X    - input vector
.   ctx  - optional user-defined context, as set by TaoSetObjectiveRoutine()
    
    Output Parameters:
.   f - function value
*/
#undef __FUNCT__
#define __FUNCT__ "FormObjectiveAndGradient"
PetscErrorCode FormObjectiveAndGradient(TAO_APPLICATION taoapp,Vec X, double *qoi,Vec G,void *ctx)
{
 PetscErrorCode info;
 PetscFunctionBegin;
 info = FormObjective(taoapp,X,qoi,ctx); CHKERRQ(info);
 info = FormGradient(taoapp,X,G,ctx); CHKERRQ(info);
 PetscFunctionReturn(0);
}
/*
    FormObjective - Evaluates the function, f(X)
    Input Parameters:
.   taoapp  - the TAO_APPLICATION context
.   X    - input vector
.   ctx  - optional user-defined context, as set by TaoSetObjectiveRoutine()

    Output Parameters:
.   f - function value
*/
#undef __FUNCT__
#define __FUNCT__ "FormObjective"
PetscErrorCode FormObjective(TAO_APPLICATION ,Vec X,double *qoi,void *ctx)
{
  AppSolve          *user = (AppSolve *) ctx;
  qoiBaseClass*      qoiOptimizer = user->qoiOptimizer();

  PetscScalar qoi_loc=0.0;
  PetscErrorCode info;

  // used for error handling
  PetscFunctionBegin;

  PetscLogStagePush(AppSolve::logstages[2]);// fnc evaluation

  //FORTRAN_NAME(update_iterno)(&qoiOptimizer->NumOptVar,&user->OptParams->tao->iter,
  //                            &,&user->OptParams->NDOF_CONTROL[1]);
  std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl
            << "optimization # " <<  user->IDOPT
            << " iteration # " << qoiOptimizer->getIteration()
            << " # of control dofs=" << qoiOptimizer->NDOF_CONTROL[1]
            << " func_eval = " << qoiOptimizer->fncEvalCnt++
            << std::endl;

  //copy global vector x locally
  PetscLogEventBegin(AppSolve::logevents[0],0,0,0,0); // tao param xfer

  info =  VecScatterBegin(qoiOptimizer->GLOCSCAT_CNTL,X,qoiOptimizer->CNTL_LOC,
                          INSERT_VALUES, SCATTER_REVERSE); CHKERRQ(info);
  info =  VecScatterEnd(  qoiOptimizer->GLOCSCAT_CNTL,X,qoiOptimizer->CNTL_LOC,
                          INSERT_VALUES, SCATTER_REVERSE); CHKERRQ(info);

  /* Store control variable data in pdeSolver data structures*/
  info = qoiOptimizer->PutLocalCntrlVars(qoiOptimizer->CNTL_LOC);
  CHKERRQ(info);

  PetscLogEventEnd(AppSolve::logevents[0],0,0,0,0); // tao param xfer

  PetscLogEventBegin(AppSolve::logevents[5],0,0,0,0); // libMesh Solve

  // forward time stepping
  info = qoiOptimizer->ForwardSolve(user); CHKERRQ(info); 

  PetscLogEventEnd(  AppSolve::logevents[5],0,0,0,0); // libMesh Solve

  /* evaluate the QOI */
  EquationSystems &es = user->get_equation_systems();// libmesh solver
  TransientFEMSystem & ideal_uncertainty_system =
   es.get_system<TransientFEMSystem>("IdealUncertaintySystem");
  qoi_loc = qoiOptimizer->ComputeObjective(user,ideal_uncertainty_system);


  // NOTE- there is an MPI_AllReduce in ComputeObjective
  *qoi = qoi_loc;

  PetscLogStagePop(); // fnc evaluation

  //  info = PetscLogFlops(nn*15); CHKERRQ(info);
  PetscFunctionReturn(0);
}
 
/* -------------------------------------------------------------------- 
    FormGradient - Evaluates the gradient, G(X). 

    Input Parameters:
.   taoapp  - the TAO_APPLICATION context
.   X    - input vector
.   ctx  - optional user-defined context, as set by TaoSetGradientRoutine()
    
    Output Parameters:
.   G - vector containing the newly evaluated gradient
    !!!note!!! when using tao_monitor, ||G||_2 is printed as the residual

this routine computes the numerical gradient of the QOI. No plotting
is done here. However, the adjoint field is stored in the data structure
and may be plotted later.  
    -------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "ZeroGradient"
PetscErrorCode ZeroGradient(TAO_APPLICATION ,Vec , Vec G,void *)
{
 PetscFunctionBegin;
 PetscErrorCode info = VecSet(G,0.0); CHKERRQ(info);
 PetscFunctionReturn(0);
}
/* -------------------------------------------------------------------- 
   FormGradient - Evaluates Gradient vector.

   Input Parameters:
.  taoapp   - the TAO_APPLICATION context
.  x     - input vector
.  ctx   - optional user-defined context, as set by TaoSetHessianRoutine()

   Output Parameters:
.  G     - Gradient 
   wrapper
   -------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormGradient"
PetscErrorCode FormGradient(TAO_APPLICATION taoapp,Vec X, Vec G,void *ctx)
{
 AppSolve *user = (AppSolve *) ctx;  
 TransientFEMSystem    *pdeSolver  = user->pdeSolver();  // abstract pde solver
 PetscFunctionBegin;
 pdeSolver->FormGradient(taoapp,X,G);
 PetscFunctionReturn(0);
}
/* -------------------------------------------------------------------- 
   FormHessian - Evaluates Hessian matrix.

Set Hessian evaluation routine by calling:

int TaoAppSetHessianRoutine(TAO_APPLICATION, int (*)(TAO_APPLICATION,Vec,Mat*,Mat*,MatStructure*,void*), void *)


The first argument is the TAO application, the second argument is the function
that evaluates the Hessian, and the third argument is a pointer to a user
defined context, cast as a void* pointer.  For solvers that evaluate the
Hessian, the matrices used to store the Hessian should be set using

TaoAppSetHessianMat(TAO_APPLICATION,Mat,Mat);

The first argument is the TAO application, the second argument is the Hessian
matrix, and the third argument is the preconditioning matrix. In most
applications, these two matrices will be the same structure.

Finite differences approximations can be used to compute the gradient and the
Hessian of an objective function. These approximations will slow down the solve
considerably and are only recommended for checking the accuracy of hand-coded
gradients and Hessians.  These routines are

TaoAppDefaultComputeGradient(TAO_APPLICATION, Vec, Vec, void*);

TaoAppDefaultComputeHessian( TAO_APPLICATION, Vec, Mat*, Mat*, MatStructure*, void*); 

and 

TaoAppDefaultComputeHessianColor( TAO_APPLICATION, Vec, Mat*, Mat*, MatStructure*, void* ); 

These routines can be set using TaoAppSetGradientRoutine() and
TaoAppSetHessianRoutine() or through the options database. If finite
differencing is used with coloring, the routine
TaoAppSetColoring(TAO_APPLICATION, ISColoring); should be used to specify the
coloring.  It is also possible to use finite difference approximations to
directly check the correctness of an application’s gradient and/or Hessian
evaluation routines. This can be done using the special TAO solver 
  -tao_method tao_fd_test   (calls TaoSolve_FD)
together with the options 

  -tao_test_gradient  (calls TaoAppDefaultComputeGradient)

              and/or 

  -tao_test_hessian  (calls TaoAppDefaultComputeHessian ) 

        AND

  -tao_test_display    -tao_fd_delta 1.e-2

   Input Parameters:
.  taoapp   - the TAO_APPLICATION context
.  x     - input vector
.  ctx   - optional user-defined context, as set by TaoSetHessianRoutine()

   Output Parameters:
.  H     - Hessian matrix

   Note:  Providing the Hessian may not be necessary.  Only some solvers
   require this matrix.  
   -------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormHessian"
int FormHessian(TAO_APPLICATION ,Vec X,Mat *HH, Mat *, MatStructure *flag,void *ctx)
{
  AppSolve *user = (AppSolve *) ctx;  
  qoiBaseClass    *qoiOptimizer = user->qoiOptimizer();
  PetscErrorCode info;
  Mat H=*HH;
  
  // Define data structures to contain the hessian diagonal for plotting the
  // variance this is also use to hold the sensitivity gradient
  // for verification of sensitivities calculations via adjoint gradient
  // sensitivity gradient comparison
  info = VecDuplicate(X,&qoiOptimizer->hessianDiagonal); CHKERRQ(info);
  info = VecSet(qoiOptimizer->hessianDiagonal,0.0); CHKERRQ(info);

  // initialize hessian to zero
  info = MatZeroEntries(H); CHKERRQ(info);

  // used for error handling
  PetscFunctionBegin; 

  // create temp vecs
  Vec            canonicalBasis,matVecProduct;
  info = VecDuplicate(X,&canonicalBasis); CHKERRQ(info);

  // verify that the matrix has dense storage
  const MatType storagetype;
  info = MatGetType(H,&storagetype);CHKERRQ(info);
  if( strcasecmp(storagetype,MATMPIDENSE) &&
      strcasecmp(storagetype,MATSEQDENSE) )
   {
       std::cerr<<"Expected Dense Hessian Matrix" <<std::endl; libmesh_error(); 
   }

  // get matrix data and setup columns of the matrix as the output of the 
  // hessian-vector product
  PetscInt       localRow;
  PetscScalar    *data;
  info = MatGetArray(H,&data);CHKERRQ(info);
  info = MatGetLocalSize(H,&localRow,PETSC_NULL);CHKERRQ(info);  // number local rows 
  info = VecCreateMPIWithArray(PETSC_COMM_WORLD,localRow,
                               PETSC_DETERMINE,PETSC_NULL,&matVecProduct);CHKERRQ(info);

  // hessian vector product with e_i, i=0,1,2,... to form Hessian
  for (PetscInt Jj=qoiOptimizer->minHessianColumn; 
                Jj<qoiOptimizer->maxHessianColumn; Jj++) 
   {
    std::cout<< "computing column "<< Jj <<" ..."<<std::endl<<std::flush;
    qoiOptimizer->hessianColumn = Jj;  // keep track of column for debugging
    // setup canonical basis
    info = VecZeroEntries(canonicalBasis);
    info = VecSetValue(canonicalBasis,Jj,1.0,INSERT_VALUES);
    info = VecAssemblyBegin(canonicalBasis); CHKERRQ(info);
    info = VecAssemblyEnd(  canonicalBasis); CHKERRQ(info);

    // setup matrix data pointer for the vector
    info = VecPlaceArray(matVecProduct,data + Jj*localRow);CHKERRQ(info);

    // mat vec product
    info = MatMult(qoiOptimizer->HessianShell,canonicalBasis,matVecProduct);
    #if defined(PETSC_USE_DEBUG)// print sensitivity gradient
    PetscPrintf(PETSC_COMM_WORLD,"Hessian Vector Product %d...\n",Jj );
    info = VecView(matVecProduct,0); CHKERRQ(info);
    #endif
    //info = VecGetArray(  matVecProduct ,&data);CHKERRQ(info);
    //for (PetscInt Ii=0; Ii<m; Ii++) {
    //   MatSetValue( H,Ii,Jj,data[Ii],INSERT_VALUES);
    //}
    //info = VecRestoreArray(matVecProduct,&data);CHKERRQ(info);

    // restore vector data pointer
    info = VecResetArray(matVecProduct);CHKERRQ(info);
   }

  // restore matrix data structure
  info = VecDestroy(matVecProduct);CHKERRQ(info);
  info = MatRestoreArray(H,&data);CHKERRQ(info);

  // assemble 
  info = MatAssemblyBegin(H,MAT_FINAL_ASSEMBLY); CHKERRQ(info);
  info = MatAssemblyEnd(  H,MAT_FINAL_ASSEMBLY); CHKERRQ(info);

  #if defined(PETSC_USE_DEBUG)// print sensitivity gradient
  PetscPrintf(PETSC_COMM_WORLD,"Sensitivity gradient...\n" );
  info = VecAssemblyBegin(qoiOptimizer->hessianDiagonal); CHKERRQ(info);
  info = VecAssemblyEnd(  qoiOptimizer->hessianDiagonal); CHKERRQ(info);
  info = VecView(qoiOptimizer->hessianDiagonal,0); CHKERRQ(info);
  #endif
  //info = MatGetDiagonal(H,hessianDiagonal);
  //OStringStream varFileName;
  //varFileName << "femvis/femqoi_" << user->GroupID
  //            << "variance_"  
  //            << "nopt_"          << user->IDOPT 
  //            << AppSolve::profileID  << ".e";
  //qoiOptimizer->plotElemData(varFileName,*(user->mesh),hessianDiagonal); 

  // write matrix to a file
  OStringStream hessFileName;
  hessFileName << "files/hessianMatrix_" << user->GroupID
               << "min_" << qoiOptimizer->minHessianColumn
               << "max_" << qoiOptimizer->maxHessianColumn
               << "nopt_"                 << user->IDOPT 
               << AppSolve::profileID     << ".m";
  PetscViewer matViewer; 
  info = PetscViewerASCIIOpen(PETSC_COMM_WORLD,hessFileName.str().c_str(),
                                               &matViewer);  CHKERRQ(info);
  info = PetscViewerSetFormat(matViewer,PETSC_VIEWER_ASCII_MATLAB); 
                                                             CHKERRQ(info);
  info = MatView(H,matViewer); CHKERRQ(info);
  info = PetscViewerFlush(matViewer); CHKERRQ(info);
  info = PetscViewerDestroy(matViewer); CHKERRQ(info);
  PetscPrintf(PETSC_COMM_WORLD,"wrote hessian Matrix...\n" );

  // check that matrix multiply gives the same answer for a non canonical basis
  Vec nonCanonicalProductOne, nonCanonicalProductTwo;
  info = VecDuplicate(X,&nonCanonicalProductOne); CHKERRQ(info);
  info = VecDuplicate(X,&nonCanonicalProductTwo); CHKERRQ(info);
  // reuse canonicalBasis as work vector
  for (PetscInt Jj=0; Jj<qoiOptimizer->NDOF_CONTROL[1] ; Jj++) 
          info = VecSetValue(canonicalBasis,Jj,1.0+Jj,INSERT_VALUES);
  info = VecAssemblyBegin(canonicalBasis); CHKERRQ(info);
  info = VecAssemblyEnd(  canonicalBasis); CHKERRQ(info);
  info = MatMult(H,canonicalBasis,nonCanonicalProductOne);CHKERRQ(info);
  info = MatMult(qoiOptimizer->HessianShell,canonicalBasis,
                                             nonCanonicalProductTwo);
  info = VecAXPY(nonCanonicalProductTwo,-1.0,nonCanonicalProductOne);
  PetscScalar productNorm;
  info = VecNorm(nonCanonicalProductTwo,NORM_2,&productNorm); CHKERRQ(info)
  if(productNorm > 1.e-3) 
     PetscPrintf(PETSC_COMM_WORLD,"\n\n\n  NOT VERIFIED!!!!!! \n\n\n");
  PetscPrintf(PETSC_COMM_WORLD,"\n\n\n  ||H(x)-Hx||_2 = %e \n\n\n",productNorm);

  // clean up and exit
  info = VecDestroy(nonCanonicalProductOne);CHKERRQ(info);
  info = VecDestroy(nonCanonicalProductTwo);CHKERRQ(info);
  info = VecDestroy(canonicalBasis);CHKERRQ(info);
  info = VecDestroy(qoiOptimizer->hessianDiagonal); CHKERRQ(info);
  *flag=SAME_NONZERO_PATTERN;
  CHKMEMQ; // check for memory corruption use -malloc_debug to enable
  PetscFunctionReturn(0);

}
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "MatrixFreeHessian"
/* 
   MatrixFreeHessian - Sets a pointer for use in computing Hessian-vector
   products.
    
   Input Parameters:
.  taoapp - the TAO_APPLICATION context
.  X    - the input vector
.  ptr  - optional user-defined context, as set by TaoSetHessian()
   
   Output Parameters:
.  H     - Hessian matrix
.  PrecH - optionally different preconditioning Hessian
.  flag  - flag indicating matrix structure
*/
int MatrixFreeHessian(TAO_APPLICATION ,Vec ,Mat *,Mat *,
                      MatStructure *,void *)
{
  /* do nothing  until the hessian vector product*/
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "hessianVectorDiagonal"
PetscErrorCode hessianVectorDiagonal(Mat Hess, Vec diagOut)
{ 
  void *ctx;
  // used for error handling
  PetscInt info;
  PetscFunctionBegin; 

  // Zero Entries first
  info = VecSet(diagOut,0.0);CHKERRQ(info);
  info = VecAssemblyBegin(diagOut); CHKERRQ(info);
  info = VecAssemblyEnd(  diagOut); CHKERRQ(info);

  // get solver context 
  info = MatShellGetContext(Hess,&ctx);CHKERRQ(info);

  PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "hessianVectorProduct"
PetscErrorCode hessianVectorProduct(Mat Hess, Vec xIn, Vec yOut)
{ 
  void *ctx;
  // used for error handling
  PetscInt info;
  PetscFunctionBegin; 

 info = MatShellGetContext(Hess,&ctx);CHKERRQ(info);
 AppSolve *user = (AppSolve *) ctx;  
 TransientFEMSystem    *pdeSolver    = user->pdeSolver();  // abstract pde solver
 info = pdeSolver->hessianVectorProduct(Hess, xIn, yOut) ;  // abstract pde solver

 PetscFunctionReturn(0);
 //  info = PetscLogFlops(nn*9); CHKERRQ(info);
}
/* ---------------------------------------------------------- 
   FIXME 
   build_system_solution_vector should be merged with 
   build_equation_system_solution_vector 
            - or - 
   build_equation_system_solution_vector  should call
   build_system_solution_vector 
 */
void build_system_solution_vector(EquationSystems &es, 
                                  System& system ,
                                  std::vector<Number>& soln)
{
  PetscFunctionBegin; 
  // This function must be run on all processors at once
  parallel_only();

  libmesh_assert (es.n_systems());

  // Get a constant reference to the mesh object.
  const libMesh::MeshBase& mesh = es.get_mesh();

  const unsigned int dim = mesh.mesh_dimension();
  const unsigned int nn  = mesh.n_nodes();
  // FIXME note this differs from build_equation_system_solution_vector 
  const unsigned int nv  = system.n_vars(); 

  // We'd better have a contiguous node numbering
  libmesh_assert (nn == mesh.max_node_id());

  // allocate storage to hold
  // (number_of_nodes)*(number_of_variables) entries.
  soln.resize(nn*nv);

  std::vector<Number>             sys_soln;

  // (Note that we use an unsigned short int here even though an
  // unsigned char would be more that sufficient.  The MPI 1.1
  // standard does not require that MPI_SUM, MPI_PROD etc... be
  // implemented for char data types. 12/23/2003 - BSK)  
  std::vector<unsigned short int> node_conn(nn);

  
  // Zero out the soln vector
  std::fill (soln.begin(),       soln.end(),       libMesh::zero);

  
  // Get the number of elements that share each node.  We will
  // compute the average value at each node.  This is particularly
  // useful for plotting discontinuous data.
  libMesh::MeshBase::const_element_iterator       e_it  = mesh.active_local_elements_begin();
  const libMesh::MeshBase::const_element_iterator e_end = mesh.active_local_elements_end(); 

  for ( ; e_it != e_end; ++e_it)
    for (unsigned int n=0; n<(*e_it)->n_nodes(); n++)
      node_conn[(*e_it)->node(n)]++;

  // Gather the distributed node_conn arrays in the case of
  // multiple processors
  Parallel::sum(node_conn);

  unsigned int var_num=0;

  // For each system in this EquationSystems object,
  // update the global solution and if we are on processor 0,
  // loop over the elements and build the nodal solution
  // from the element solution.  Then insert this nodal solution
  // into the vector passed to build_solution_vector.

      const unsigned int nv_sys = system.n_vars();
      
      system.current_local_solution->localize_to_one(sys_soln,
                                                     libMesh::processor_id() );
      
      std::vector<Number>       elem_soln;   // The finite element solution
      std::vector<Number>       nodal_soln;  // The FE solution interpolated to the nodes
      std::vector<unsigned int> dof_indices; // The DOF indices for the finite element 
      
      for (unsigned int var=0; var<nv_sys; var++)
	{
	  const FEType& fe_type    = system.variable_type(var);
	  
	  libMesh::MeshBase::const_element_iterator       it  = mesh.active_local_elements_begin();
	  const libMesh::MeshBase::const_element_iterator end = mesh.active_local_elements_end(); 

	  for ( ; it != end; ++it)
	    {
	      const Elem* elem = *it;
	      system.get_dof_map().dof_indices (elem, dof_indices, var);
	      
	      elem_soln.resize(dof_indices.size());
	      
	      for (unsigned int i=0; i<dof_indices.size(); i++)
		elem_soln[i] = sys_soln[dof_indices[i]];
		  
	      FEInterface::nodal_soln (dim,
				       fe_type,
				       elem,
				       elem_soln,
				       nodal_soln);

#ifdef ENABLE_INFINITE_ELEMENTS
	      // infinite elements should be skipped...
	      if (!elem->infinite())
#endif
		{ 
		  libmesh_assert (nodal_soln.size() == elem->n_nodes());
		  
		  for (unsigned int n=0; n<elem->n_nodes(); n++)
		    {
		      libmesh_assert (node_conn[elem->node(n)] != 0);
		      soln[nv*(elem->node(n)) + (var + var_num)] +=
			nodal_soln[n]/static_cast<Real>(node_conn[elem->node(n)]);
		    }
		}
	    }	 
	}

      var_num += nv_sys;

  // Now each processor has computed contriburions to the
  // soln vector.  Gather them all up.
  Parallel::sum(soln);
  PetscFunctionReturnVoid(); 
}
/* ---------------------------------------------------------- */
void build_equation_system_solution_vector(EquationSystems &es, 
                                           std::vector<Number>& soln)
{
  PetscFunctionBegin; 
  // This function must be run on all processors at once
  parallel_only();

  libmesh_assert (es.n_systems());

  // Get a constant reference to the mesh object.
  const libMesh::MeshBase& mesh = es.get_mesh();

  const unsigned int dim = mesh.mesh_dimension();
  const unsigned int nn  = mesh.n_nodes();
  // FIXME note this differs from build_system_solution_vector 
  const unsigned int nv  = es.n_vars();

  // We'd better have a contiguous node numbering
  libmesh_assert (nn == mesh.max_node_id());

  // allocate storage to hold
  // (number_of_nodes)*(number_of_variables) entries.
  soln.resize(nn*nv);

  std::vector<Number>             sys_soln;

  // (Note that we use an unsigned short int here even though an
  // unsigned char would be more that sufficient.  The MPI 1.1
  // standard does not require that MPI_SUM, MPI_PROD etc... be
  // implemented for char data types. 12/23/2003 - BSK)  
  std::vector<unsigned short int> node_conn(nn);

  
  // Zero out the soln vector
  std::fill (soln.begin(),       soln.end(),       libMesh::zero);

  
  // Get the number of elements that share each node.  We will
  // compute the average value at each node.  This is particularly
  // useful for plotting discontinuous data.
  libMesh::MeshBase::const_element_iterator       e_it  = mesh.active_local_elements_begin();
  const libMesh::MeshBase::const_element_iterator e_end = mesh.active_local_elements_end(); 

  for ( ; e_it != e_end; ++e_it)
    for (unsigned int n=0; n<(*e_it)->n_nodes(); n++)
      node_conn[(*e_it)->node(n)]++;

  // Gather the distributed node_conn arrays in the case of
  // multiple processors
  Parallel::sum(node_conn);

  unsigned int var_num=0;

  // For each system in this EquationSystems object,
  // update the global solution and if we are on processor 0,
  // loop over the elements and build the nodal solution
  // from the element solution.  Then insert this nodal solution
  // into the vector passed to build_solution_vector.

  for (unsigned int i_sys = 0 ; i_sys < es.n_systems() ; i_sys++)
    {  
      System& system  = es.get_system(i_sys) ; 
      const unsigned int nv_sys = system.n_vars();
      
      system.current_local_solution->localize_to_one(sys_soln,
                                                     libMesh::processor_id() );
      
      std::vector<Number>       elem_soln;   // The finite element solution
      std::vector<Number>       nodal_soln;  // The FE solution interpolated to the nodes
      std::vector<unsigned int> dof_indices; // The DOF indices for the finite element 
      
      for (unsigned int var=0; var<nv_sys; var++)
	{
	  const FEType& fe_type    = system.variable_type(var);
	  
	  libMesh::MeshBase::const_element_iterator       it  = mesh.active_local_elements_begin();
	  const libMesh::MeshBase::const_element_iterator end = mesh.active_local_elements_end(); 

	  for ( ; it != end; ++it)
	    {
	      const Elem* elem = *it;
	      system.get_dof_map().dof_indices (elem, dof_indices, var);
	      
	      elem_soln.resize(dof_indices.size());
	      
	      for (unsigned int i=0; i<dof_indices.size(); i++)
		elem_soln[i] = sys_soln[dof_indices[i]];
		  
	      FEInterface::nodal_soln (dim,
				       fe_type,
				       elem,
				       elem_soln,
				       nodal_soln);

#ifdef ENABLE_INFINITE_ELEMENTS
	      // infinite elements should be skipped...
	      if (!elem->infinite())
#endif
		{ 
		  libmesh_assert (nodal_soln.size() == elem->n_nodes());
		  
		  for (unsigned int n=0; n<elem->n_nodes(); n++)
		    {
		      libmesh_assert (node_conn[elem->node(n)] != 0);
		      soln[nv*(elem->node(n)) + (var + var_num)] +=
			nodal_soln[n]/static_cast<Real>(node_conn[elem->node(n)]);
		    }
		}
	    }	 
	}

      var_num += nv_sys;
    }

  // Now each processor has computed contriburions to the
  // soln vector.  Gather them all up.
  Parallel::sum(soln);
  PetscFunctionReturnVoid(); 
}
/* ---------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "TaoConverged_CheckPoint"
/*@ 
   TaoConverged_CheckPoint- Determines whether the solver should continue 
                            iterating or terminate. It calls the default 
                            Tao Convergence routine but also checkpoints
                            the data of the solve periodically.

   Collective on TAO_SOLVER

   Input Parameters:
+  tao - the TAO_SOLVER context
-  dummy - unused dummy context

   Output Parameter:    reason - for terminating

   Notes:
   This routine checks the residual in the optimality conditions, the 
   relative residual in the optimity conditions, the number of function
   evaluations, and the function value to test convergence.  Some
   solvers may use different convergence routines.

   Level: developer

.seealso: TaoSetTolerances(),TaoGetTerminationReason(),TaoSetTerminationReason()
@*/
PetscErrorCode TaoConverged_CheckPoint(TAO_SOLVER tao,void *ctx)
{
  AppSolve *user = (AppSolve *) ctx;  
  qoiBaseClass    *qoiOptimizer = user->qoiOptimizer();
  libMesh::MeshBase                &mesh = user->get_mesh();
  TransientFEMSystem    *pdeSolver    = user->pdeSolver();
  PetscInt commands[2], // command buffer
           TaskID ; // task identification for control task
  TaoTerminateReason reason;
  PetscErrorCode info; /* used to check for functions returning nonzeros */
  PetscFunctionBegin; 

  // call the default Tao Convergence Tester
  TaoConverged_Default(tao,ctx); 
  // get the termination reason
  info = TaoGetTerminationReason(tao,&reason); CHKERRQ(info);
  // check if the ideal field requires that we kick out of this optimization
  //  solve and go to the next optimization solve
  if( qoiOptimizer->UpdateQoiData(user) ) reason=TAO_CONVERGED_USER;
  //if(reason!=TAO_CONTINUE_ITERATING) reason=TAO_CONTINUE_ITERATING;
  // store the possibly modified the termination reason
  info = TaoSetTerminationReason(tao, reason); CHKERRQ(info);

  // plot after AppSolve::checkpoint iterations 
  if( user->CheckPoint(tao->iter) || reason!=TAO_CONTINUE_ITERATING)
   {
     user->UpdateCheckPoint(tao->iter); 
     commands[0] = AppSolve::outputfileID;

     if(reason!=TAO_CONTINUE_ITERATING) TaskID =COMPLETE;//SOLVE COMPLETE
     else                               TaskID =CHKPOINT;//THIS IS A CHECKPOINT
     commands[1] = TaskID;
     // tell control task that the next control file has been written
     if( !libMesh::processor_id() ) 
      {
       qoiOptimizer->SendToServer(commands,user->GroupID);
       info = pdeSolver->WriteControlFile(AppSolve::outputfileID,user->GroupID);
       CHKERRQ(info);
       //FORTRAN_NAME(vis_power)(&user->IDOPT);
      }
     // plot
     if( (TaskID == COMPLETE && qoiOptimizer->PLOTOPTSOLVE) || user->PlotIteration() ) 
       {
         // tell control task that the vis solve is complete
         if( !libMesh::processor_id() ) 
            qoiOptimizer->SendToServer(commands,VIS_SOLVE_COMPLETE);

          qoiOptimizer->PlotStoredData(user);

          // plot the gradient
          TaoVecPetsc *taoGrad=NULL;
          TaoMethod methodType;
          if( qoiOptimizer->TAO() )
           {
            info = TaoGetMethod( qoiOptimizer->TAO() ,&methodType); 
            if(std::string(methodType).find("tao_blmvm")!=std::string::npos) 
             {
               TAO_BLMVM *blm = (TAO_BLMVM *)qoiOptimizer->TAO()->data;
               taoGrad = dynamic_cast <TaoVecPetsc *> (blm->G);
             }
           }
          OStringStream gradFileName;
          PetscInt TaoIteration = 0 ;
          if( qoiOptimizer->TAO() ) 
              TaoIteration = qoiOptimizer->TAO()->iter ;
          gradFileName << "femvis/femqoi_" << user->GroupID
                      << "gradient_iter_"  << TaoIteration 
                      << "nopt_"           << user->IDOPT 
                      << AppSolve::profileID  << ".e";
          if(taoGrad) qoiOptimizer->plotElemData(gradFileName,
                                                 mesh,taoGrad->GetVec());
          std::cout << "end plot..." << std::endl;
       }

     ////get the output file ID to write next from the control task
     //if(TaskID == CHKPOINT) //THIS IS A RESTART/CHECKPOINT
     //   MPI_Bcast(&GenInfo.outputfileID,1,MPI_INT,0,user->ControlComm);
   }

  if(reason!=TAO_CONTINUE_ITERATING)
   {
    MatStructure matrixFlag;
    Vec solution;
    info = TaoAppGetSolutionVec(qoiOptimizer->TAOAPP(),&solution);CHKERRQ(info);
    info = FormHessian(qoiOptimizer->TAOAPP(),solution,
                &qoiOptimizer->Hessian,&qoiOptimizer->Hessian,&matrixFlag,ctx);
   }

  PetscFunctionReturn(0);
}
