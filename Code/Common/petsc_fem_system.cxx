//libmesh
#include "libmesh.h"
#include "equation_systems.h"
#include "time_solver.h"
#include "petsc_vector.h"
#include "petsc_matrix.h"
#include "petsc_diff_solver.h"
#include "sparse_matrix.h"
#include "dof_map.h"
#include "exodusII_io.h"

//local
#include "optimizationParameter.h" 
#include "petsc_fem_system.h"
#include "petsc_fem_context.h" 

/* ------------------------------------------------------------ */
// PetscFEMSystem implementation
PetscFEMSystem::PetscFEMSystem (EquationSystems& es, 
                                         const std::string& name,
					 const unsigned int number) :
  FEMSystem                 (es, name, number)
{ 
  // default to the petsc linear solver
  m_LinearSolve     = PETSC_TRUE ; 

  // control jacobian assembly
  m_MassMatrix      = PETSC_TRUE ; 
  m_StiffnessMatrix = PETSC_TRUE ; 
  m_jacobianComputed= PETSC_FALSE ; 

  // power time step to use
  m_PowerID  = 0; 

  // system profile
  PetscLogEventRegister("PetscFEMSysResidual",PETSC_VIEWER_COOKIE,&PetscFEMSysLogResidual);
  PetscLogEventRegister("PetscFEMSysJacobian",PETSC_VIEWER_COOKIE,&PetscFEMSysLogJacobian);
  PetscLogEventRegister("PetscFEMSysDefault" ,PETSC_VIEWER_COOKIE,&PetscFEMSysLogDefault );
}


Number PetscFEMSystem::
stored_solution( const unsigned int &global_dof_number,
                 const unsigned int id_vec) const
{
  // Check the sizes
  libmesh_assert (global_dof_number < this->get_dof_map().n_dofs());

  OStringStream vector_name;
  vector_name << "stored_local_solution" << id_vec;
  return (this->get_vector( vector_name.str() ))(global_dof_number);
}


Number PetscFEMSystem::
prev_solution     ( const unsigned int &global_dof_number) const
{
  // Check the sizes
  libmesh_assert (global_dof_number < this->get_dof_map().n_dofs());
  return (this->get_vector( "old_local_solution" ))(global_dof_number);
}


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
void PetscFEMSystem:: solve()
{
  PetscFunctionBegin;

  // print recomputation flags
  this->PetscFEMSystem::printSelf(std::cout);

  // try a petsc linear solve if possible 
  // else default to the nonlinear solver
  if( this->m_LinearSolve )
    {
     // assemble 
     this->assembly(true,true);
     this->matrix->close();
     this->rhs->close();


     // setup linear solver
     PetscErrorCode info;
     // set SNES solver to the equivalent of
     //   -snes_ls basicnonorms -snes_no_convergence_test -snes_max_it 1
     //info = SNESSetTolerances(SystemSolver->snes(),PETSC_DEFAULT,
     //              PETSC_DEFAULT,PETSC_DEFAULT,1,PETSC_DEFAULT); CHKERRV(info);
     //info = SNESSetConvergenceTest(SystemSolver->snes(),SNESSkipConverged,
     //                              PETSC_NULL,PETSC_NULL); CHKERRV(info);
     //info = SNESLineSearchSet(SystemSolver->snes(),SNESLineSearchNoNorms,
     //                              PETSC_NULL); CHKERRV(info);

     /* set Krylov method */
     KSP PetscLinearKSP; ///< local context for linear solve if needed

     info = KSPCreate(libMesh::COMM_WORLD,&PetscLinearKSP);CHKERRV(info);
     // same non zero pattern for the multiple solves
     PetscMatrix<Number> &matrix =
       *(libmesh_cast_ptr<PetscMatrix<Number>*>(this->matrix));
     info = KSPSetOperators(PetscLinearKSP,matrix.mat(),matrix.mat(),
                            DIFFERENT_NONZERO_PATTERN ); CHKERRV(info);
     PC pcPre;
     info = KSPGetPC(PetscLinearKSP,&pcPre); CHKERRV(info);
     info = PCSetType(pcPre,PCBJACOBI); CHKERRV(info);

     // set defaults
     PetscInt maxiter = 1000;
     info = KSPSetTolerances(PetscLinearKSP,PETSC_DEFAULT,
                             PETSC_DEFAULT,PETSC_DEFAULT,
                             maxiter); CHKERRV(info);
     // initial guess should be zero
     info = KSPSetInitialGuessNonzero(PetscLinearKSP,PETSC_FALSE);CHKERRV(info);

     // allow last minute changes from options
     info = KSPSetFromOptions(PetscLinearKSP);CHKERRV(info);
     //info = KSPMonitorCancel(PetscLinearKSP);CHKERRV(info);
     //PetscLinearKSP->printreason = PETSC_FALSE;
   

     // get proper pointers
     PetscVector<Number> &solution =
       *(libmesh_cast_ptr<PetscVector<Number>*>(this->solution.get()));
     PetscVector<Number> &residual =
       *(libmesh_cast_ptr<PetscVector<Number>*>(this->rhs));

     //PetscDiffSolver* SystemSolver = 
     //                libmesh_cast_ptr<PetscDiffSolver*>(
     //                & (*this->time_solver->diff_solver()) );
     //KSP  snesksp;
     //info = SNESGetKSP(SystemSolver->snes(),&snesksp);CHKERRV(info);

     // solve
     // overwrite residual w/ Newton Step
     info = KSPSolve(PetscLinearKSP,residual.vec(),residual.vec() );CHKERRV(info);

     // Notes:
     // The Newton-like methods typically solve linear systems of the form
     //     f'(x) x = -f(x),
     // where f'(x) denotes the Jacobian matrix and f(x) is the function.
     // 
     // in this case  for the linear system (IGNORING MINUS SIGN ON RESIDUAL)
     //     A s = (A x0 - b) = r(x) 
     //     A x =         b 
     //         =>  x = x0 - s
     info = VecAXPY( solution.vec(), -1.0, residual.vec() );

     // localize 
     this->solution->localize( *this->current_local_solution );

     // clean up
     info = KSPDestroy(PetscLinearKSP);CHKERRV(info);
    }
  else
    {
     Parent::solve(); // default nonlinear solver
    }

  // localize
  this->solution->localize(*this->current_local_solution);
  PetscFunctionReturnVoid();
}

// setup scatter context to extract subvectors
void PetscFEMSystem::SetupSubVector()
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  //  setup dirichlet data. must be the same on all procs or will freeze
  const DofMap & dof_map    = this->get_dof_map();
  // Get a constant reference to the mesh object.
  libMesh::MeshBase& mesh = this->get_equation_systems().get_mesh();

  // get reference to petsc vector to build scatter contexts
  Vec Soln    = (dynamic_cast< PetscVector<double>* > (&(*this->solution)) )->vec();

  // create space
  m_subvector_indices.resize(this->n_vars()) ;
  m_subvector_scatter.resize(this->n_vars()) ;
  m_subvector_is.resize(     this->n_vars()) ;

  // loop over variables and build global indicies
  for (unsigned int var=0; var<this->n_vars(); var++)
    {
      m_subvector_indices[var].clear();
      std::vector<unsigned int> dof_indices_var;
      // element based setup of dirichlet data
      libMesh::MeshBase::const_element_iterator       el    =mesh.active_elements_begin();
      const libMesh::MeshBase::const_element_iterator end_el_full=mesh.active_elements_end();
      for ( ; el != end_el_full; ++el)
       {
         Elem* elem = *el;
         dof_map.dof_indices (elem, dof_indices_var,var);
         for( unsigned int Ii = 0 ; Ii < dof_indices_var.size(); Ii++)
               m_subvector_indices[var].push_back(dof_indices_var[Ii]);
       } // end element loop
      // sort then erase duplicates
      std::vector<unsigned int>::iterator pos;
      std::sort(        m_subvector_indices[var].begin(),m_subvector_indices[var].end());
      pos = std::unique(m_subvector_indices[var].begin(),m_subvector_indices[var].end());
      m_subvector_indices[var].erase( pos,               m_subvector_indices[var].end()); 

      // Construct index sets
      IS local_is ;
      ierr = ISCreateStride( libMesh::COMM_WORLD,
            		 m_subvector_indices[var].size(),0,1,&local_is);

      ierr = ISCreateGeneral(libMesh::COMM_WORLD,
            		 m_subvector_indices[var].size(),
                 (int*) &m_subvector_indices[var][0], 
                        &m_subvector_is[var]); 

      // create subvector to scatter to
      Vec subvector; 
      ierr = VecCreateMPI(libMesh::COMM_WORLD,
                          PETSC_DECIDE,                 // n_local
                          m_subvector_indices[var].size(),// n_global
                         &subvector); CHKERRV(ierr);

      // Construct the scatter objects
      //   SCATTER_FORWARD FROM SUBVECTOR TO FULL VECTOR
      //   SCATTER_REVERSE FROM FULL VECTOR TO SUBVECTOR 
      ierr = VecScatterCreate(  subvector, local_is, Soln, 
                              m_subvector_is[var],
            		     &m_subvector_scatter[var]);

      // Clean up 
      ierr = ISDestroy(local_is  );    
      ierr = VecDestroy(subvector);    
  
    }

  PetscFunctionReturnVoid();
}

// get sub vector 
// note the Vec is created here and 
// the calling function is in charge of memory clean up
void PetscFEMSystem::GetSubVector (Vec &SubVector,int VarID)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;

  Vec Soln    = (dynamic_cast< PetscVector<double>* > (&(*this->solution)) )->vec();
  ierr = VecCreateMPI(libMesh::COMM_WORLD,
   		  PETSC_DECIDE,          // n_local
            this->m_subvector_indices.at(VarID).size(),   // n_global
   		  &SubVector); CHKERRV(ierr);
  // get subvector solution data
  ierr=VecScatterBegin(this->m_subvector_scatter[VarID],Soln,SubVector,INSERT_VALUES,SCATTER_REVERSE);
  ierr=VecScatterEnd(  this->m_subvector_scatter[VarID],Soln,SubVector,INSERT_VALUES,SCATTER_REVERSE);

  PetscFunctionReturnVoid();
}
// get sub vector 
// note the Vec is created here and 
// the calling function is in charge of memory clean up
void PetscFEMSystem::SetSubVector (Vec &SubVector,int VarID)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;

  Vec Soln    = (dynamic_cast< PetscVector<double>* > (&(*this->solution)) )->vec();
  // set subvector solution data
  ierr=VecScatterBegin(this->m_subvector_scatter[VarID],SubVector,Soln,INSERT_VALUES,SCATTER_FORWARD);
  ierr=VecScatterEnd(  this->m_subvector_scatter[VarID],SubVector,Soln,INSERT_VALUES,SCATTER_FORWARD);

  PetscFunctionReturnVoid();
}
void PetscFEMSystem::assembly(bool get_residual, bool get_jacobian)
{
 PetscFunctionBegin;

 // apply BC for solution then assemble
 this->ApplyDirichlet();

 PetscPrintf(PETSC_COMM_WORLD,"PetscFEMSystem::assembly residual %d jacobian %d\n", 
                              get_residual,get_jacobian);
 // only assemble jacobian once...
 if( this->m_LinearSolve && this->m_MassMatrix && this->m_StiffnessMatrix)
   {// treat residual and jacobian separately 
    if( get_residual )
      {
      PetscPrintf(PETSC_COMM_WORLD,"linear assembly residual...\n");
      PetscLogEventBegin(PetscFEMSysLogResidual,0,0,0,0); // log residual
      this->FEMSystem::assembly(get_residual,false);
      PetscLogEventEnd(  PetscFEMSysLogResidual,0,0,0,0); // log residual
      }

    if( !m_jacobianComputed && get_jacobian )
      {
      PetscPrintf(PETSC_COMM_WORLD,"linear assembly jacobian...\n");
      PetscLogEventBegin(PetscFEMSysLogJacobian,0,0,0,0); // log jacobian
      this->FEMSystem::assembly(false,get_jacobian);
      m_jacobianComputed = PETSC_TRUE;
      PetscLogEventEnd(  PetscFEMSysLogJacobian,0,0,0,0); // log jacobian
      }
   }
 else
   {// default... no special assembly needed ... call base class
    PetscLogEventBegin(PetscFEMSysLogDefault,0,0,0,0); // log default
    PetscPrintf(PETSC_COMM_WORLD,"default assembly...\n");
    this->FEMSystem::assembly(get_residual,get_jacobian);
    PetscLogEventEnd(  PetscFEMSysLogDefault,0,0,0,0); // log default
   }

 // apply dirichlet data if any
 for( unsigned int i_var = 0 ; i_var < this->n_vars() ; i_var++)
  {
   std::vector<PetscInt>::iterator setIDIter ; 
   for(setIDIter  = this->m_DirichletNodeSets[i_var].begin(); 
       setIDIter != this->m_DirichletNodeSets[i_var].end();  setIDIter++)
    {
     std::vector<PetscInt> &dirichletNodeSet= this->m_NodeSets[i_var][*setIDIter];
     //  only apply dirichlet data for typical assembly
     if( dirichletNodeSet.size() && 
         this->m_MassMatrix && this->m_StiffnessMatrix )
       {
         this->rhs->close();
         for( unsigned int Ii = 0; Ii<dirichletNodeSet.size();Ii++) 
                           this->rhs->set(dirichletNodeSet[Ii],0.0);
         this->rhs->close();
         this->matrix->close();
         this->matrix->zero_rows(dirichletNodeSet,1.0);
       }
    }
  }
 PetscFunctionReturnVoid();
}

