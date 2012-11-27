//libmesh
#include "libmesh.h"
#include "equation_systems.h"
#include "fe_base.h"
#include "fe_interface.h"
#include "quadrature.h"
#include "time_solver.h"
#include "petsc_vector.h"
#include "petsc_matrix.h"
#include "petsc_diff_solver.h"
#include "steady_solver.h"
#include "sparse_matrix.h"
#include "dof_map.h"
#include "exodusII_io.h"
#include "boundary_info.h"

//local
#include "optimizationParameter.h" 
#include "petsc_fem_system.h"

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

  NumericVector<Number> & GlobalVec  = this->add_vector("old_global_solution",true, PARALLEL);
  #ifdef LIBMESH_ENABLE_GHOSTED
  NumericVector<Number> & StorageVec = this->add_vector( "old_local_solution",true, GHOSTED);
  #else
  NumericVector<Number> & StorageVec = this->add_vector( "old_local_solution",true, SERIAL);
  #endif
  // system profile
  PetscLogEventRegister("PetscFEMSysResidual",PETSC_VIEWER_COOKIE,&PetscFEMSysLogResidual);
  PetscLogEventRegister("PetscFEMSysJacobian",PETSC_VIEWER_COOKIE,&PetscFEMSysLogJacobian);
  PetscLogEventRegister("PetscFEMSysDefault" ,PETSC_VIEWER_COOKIE,&PetscFEMSysLogDefault );
}

/* ------------------------------------------------------------ */
void PetscFEMSystem::SetupDirichlet(libMesh::MeshBase& mesh)
{
  PetscFunctionBegin;

 //  setup dirichlet data. must be the same on all procs or will freeze
 const DofMap & dof_map    = this->get_dof_map();

 // error check
 if( !this->m_NodeSets.empty() ) 
   {
     std::cout << "NodeSets not empty " << std::endl;
     //libmesh_error();
   }

 //  hold variable dofs
 std::vector<unsigned int> dof_indices_var;

 // element based setup of dirichlet data
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

   for( unsigned int i_var = 0 ; i_var < this->n_vars() ; i_var++)
    {
     // degree of freedom indices for individual components
     dof_map.dof_indices (elem,   dof_indices_var ,i_var);
     //if( this->m_MathModel.dirichletData(i_var,subdomain_id) ) 
     // { 
       for(unsigned int Ii = 0 ; Ii < dof_indices_var.size() ; Ii++)
         {
          const Node* node = elem->get_node(Ii);
          std::vector<short int> bcdata = mesh.boundary_info->boundary_ids (node);
          // store all node set IDs
          for(unsigned int Jj = 0 ; Jj < bcdata.size() ; Jj++)
            {
             this->m_NodeSets[i_var][ static_cast<unsigned int>(bcdata[Jj] ) ].push_back(
                                                             dof_indices_var[Ii] );
            }
         }
     // }
     //// nodal based setup of dirichlet data
     //for(unsigned int Ii = 0 ; Ii < dof_indices_var.size() ; Ii++)
     //if( this->m_MathModel.dirichletNodalData( elem->point(Ii) ) ) 
     //    this->m_NodeSets[1].push_back( dof_indices_var[Ii] ); 
     // look for dirichlet side sets
     for (unsigned int s=0; s<elem->n_sides(); s++)
       if (elem->neighbor(s) == NULL)
        {
         short int bc_id = mesh.boundary_info->boundary_id (elem,s);
         AutoPtr<Elem> side (elem->build_side(s));
         // Loop over the nodes on the side.
         for (unsigned int ns=0; ns<side->n_nodes(); ns++)
           for(unsigned int Ii = 0 ; Ii < dof_indices_var.size() ; Ii++)
             if (elem->node(Ii) == side->node(ns))
               this->m_NodeSets[i_var][bc_id].push_back( dof_indices_var[Ii] ); 
        } // end if (elem->neighbor(side) == NULL)
    } // end loop over variables
  } // end element loop

  // remove all duplicates
  for( unsigned int i_var = 0 ; i_var < this->n_vars() ; i_var++)
   {
     std::map<PetscInt, std::vector<PetscInt> > &ivar_NodeSets = 
                                             this->m_NodeSets[i_var];
     for( std::map<PetscInt,std::vector<PetscInt> >::iterator 
                       nsIt= ivar_NodeSets.begin(); 
                       nsIt!=ivar_NodeSets.end(); ++nsIt)
       {
         // store reference 
         std::vector<PetscInt> &nodeSetData = (*nsIt).second;

         // sort then erase duplicates
         std::sort( nodeSetData.begin(), nodeSetData.end() );
         std::vector<PetscInt>::iterator pos;

         pos = std::unique(nodeSetData.begin(),nodeSetData.end());
         nodeSetData.erase( pos,nodeSetData.end() ); 
     
         std::cout << " var # "   << i_var
                   << " node set " << (*nsIt).first 
                   << " size "<< nodeSetData.size() << std::endl;
         for(unsigned int Ii=0;Ii < nodeSetData.size();Ii++)
            std::cout << nodeSetData[Ii] << " " ; 
         std::cout << std::endl;

       }
   }

  PetscFunctionReturnVoid();
}
/* -------------------------------------------------------------------- 
   Here we define the initialization routine for the
   Convection-Diffusion system.  This routine is
   responsible for applying the initial conditions to
   the system.
   -------------------------------------------------------------------- */ 
#undef __FUNCT__
#define __FUNCT__ "PetscFEMSystem::SetupInitialConditions"
void PetscFEMSystem::SetupInitialConditions ()
{
  PetscFunctionBegin;

  // Get a constant reference to the mesh object.
  EquationSystems & es = this->get_equation_systems();

  // Get a constant reference to the mesh object.
  libMesh::MeshBase& mesh = es.get_mesh();

  // The number of variables in this system
  const unsigned int n_variables = this->n_vars();

  // The dimensionality of the current mesh
  const unsigned int dim = mesh.mesh_dimension();

  // The DofMap for this system
  const DofMap& dof_map = this->get_dof_map();

  // The new element coefficients
  DenseVector<Number> Ue;

  // Get the initial condition data
  bool ICerror = false;
  if( this->InitValues.size() != this->n_vars() ) ICerror = true;

  // error check
  if(ICerror)
   {
     std::cout << "Initial Conditions Improperly Setup" << std::endl;
     std::cout << "found # IC = " << this->InitValues.size()  << std::endl;
     std::cout << "Expected # System Variables = " << this->n_vars()  
               << " Expected # Domains = "         
               << this->get_num_elem_blk()  << std::endl;
     libmesh_error();
   }

  // scatter global parameters to local copy
  this->ScatterParametersLocally();

  // loop over subdomains so that the last subdomain has the final
  // overwrite of the nodal temperature
  for ( unsigned int id_sub = 0 ; 
        id_sub < static_cast<unsigned int>(this->get_num_elem_blk()); ++id_sub)
   {    
     // Loop over all the variables in the system
     for (unsigned int var=0; var<n_variables; var++)
       {
         PetscPrintf(PETSC_COMM_WORLD,
                     "user_initialization subdomain %d var %d\n",id_sub,var);
         // Get FE objects of the appropriate type
         const FEType& fe_type = dof_map.variable_type(var);     
         AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));      

         // Prepare variables for projection
         AutoPtr<QBase> qrule     (fe_type.default_quadrature_rule(dim));

         // The values of the shape functions at the quadrature
         // points
         //const std::vector<std::vector<Real> >& phi = fe->get_phi();

         const FEContinuity cont = fe->get_continuity();

         // only setup for C0 first order elements
         libmesh_assert (cont          == C_ZERO);
         libmesh_assert (fe_type.order == FIRST);

         // The Jacobian * quadrature weight at the quadrature points
         //const std::vector<Real>& JxW = fe->get_JxW();
        
         // The global DOF indices
         std::vector<unsigned int> dof_indices;
         // Side/edge DOF indices
         std::vector<unsigned int> side_dofs;
   
         // Now we will loop over all the elements in the mesh that
         // live on the global MESH. TO have the proper initial conditions
         // MUST LOOP OVER ENTIRE MESH SO THAT THERE IS NO POSSIBILITY OF
         // OVERWRITE. using 
         //     active_local_elements_begin, active_local_elements_end 
         // will cause errors and overwrite the subdomain initial conditions
         // on parallel assembly. Since the mesh
         // may be refined we want to only consider the ACTIVE elements,
         // hence we use a variant of the \p active_elem_iterator.
         libMesh::MeshBase::const_element_iterator       el    =mesh.active_local_elements_begin();
         const libMesh::MeshBase::const_element_iterator end_el=mesh.active_local_elements_end();
         for ( ; el != end_el; ++el)
            {    
             // Store a pointer to the element we are currently
             // working on.  This allows for nicer syntax later.
             const Elem* elem = *el;
             
             // Store element subdomain id to allow different sets of equations
             //  to be solved on different parts of the mesh
             std::vector<unsigned int> param_dof_indices;
             libMesh::System &template_parameter_system = 
                                 this->get_equation_systems().get_system("k_0");
             template_parameter_system.get_dof_map().dof_indices (elem, param_dof_indices);

             // current_solution calls map_global_to_local_index that will map
             // this global index to the local current_local_solution
             const unsigned int field_id = param_dof_indices[0];
             const unsigned int subdomain_id = static_cast<int>( elem->subdomain_id() ) - 1;

             // only proceed if the element belongs to the current subdomain
             if(id_sub != subdomain_id) continue;

             // Update the DOF indices for this element based on
             // the current mesh
             dof_map.dof_indices (elem, dof_indices, var);

             // The number of DOFs on the element
             const unsigned int n_dofs = dof_indices.size();

             // Fixed vs. free DoFs on edge/face projections
             std::vector<char> dof_is_fixed(n_dofs, false); // bools
             std::vector<int> free_dof(n_dofs, 0);

             // The element type
             const ElemType elem_type = elem->type();

             // The number of nodes on the new element
             const unsigned int n_nodes = elem->n_nodes();

             // Zero the interpolated values
             Ue.resize (n_dofs); Ue.zero();

             // In general, we need a series of
             // projections to ensure a unique and continuous
             // soln.  We start by interpolating nodes, then
             // hold those fixed and project edges, then
             // hold those fixed and project faces, then
             // hold those fixed and project interiors

             // store pointer to element
             es.parameters.set<const Elem*>("elem") = elem ; 

             // Interpolate node values first
             unsigned int current_dof = 0;
             for (unsigned int n=0; n!= n_nodes; ++n)
               {
                 // FIXME: this should go through the DofMap,
                 // not duplicate dof_indices code badly!
                 const unsigned int nc =
           	 FEInterface::n_dofs_at_node (dim, fe_type, elem_type,n);
                 if (!elem->is_vertex(n))
                   {
                     current_dof += nc;
                     continue;
                   }
                 // Assume that C_ZERO elements have a single nodal
                 // value shape function
                 libmesh_assert(nc == 1);

                 // call function pointers to member functions to setup IC
                 Ue(current_dof) = CALL_MEMBER_FN(this, this->InitValues[var])
                                    (subdomain_id,field_id,elem->point(n),es.parameters);
                 dof_is_fixed[current_dof] = true;
                 current_dof++;
               }

             // Make sure every DoF got reached!
             for (unsigned int i=0; i != n_dofs; ++i)
               libmesh_assert(dof_is_fixed[i]);

             const unsigned int
               first = this->solution->first_local_index(),
               last  = this->solution->last_local_index();

             // put soln in parallel data structures...
             for (unsigned int i = 0; i < n_dofs; i++) 
               // We may be projecting a new zero value onto
               // an old nonzero approximation - RHS
               // if (Ue(i) != 0.)
               if ((dof_indices[i] >= first) &&
                   (dof_indices[i] <  last))
                 {
                   this->solution->set( dof_indices[i], Ue(i));
                 }
           }  // end elem loop
       } // end variables loop
    // close after each subdomain so that the last subdomain will have
    // the final write...
    this->solution->close();
  } // end loop over subdomains    

  // copy to old data structures
  this->get_vector("old_global_solution") = *this->solution;

  // copy parallel data structures to local data structures
  this->solution->localize(*this->current_local_solution);

  // setup scatter context to subvector of variables 
  this->SetupSubVector();

  // setup dirichlet nodes
  this->SetupDirichlet(mesh);
  
  // **FIX Non zero Pattern of the Jacobian**
  // build jacobian once and set nonzeros
  this->FEMSystem::assembly(false,true);
  // make sure matrix is closed
  this->matrix->close();
  PetscMatrix<Number> &matrix =
      *(libmesh_cast_ptr<PetscMatrix<Number>*>(this->matrix));
  /*
     Tell the matrix we will never add a new nonzero location to the
     matrix. If we do, it will generate an error.
  */
  PetscPrintf(PETSC_COMM_WORLD,"ThermalTherapySystem::SetupInitialConditions fixing nonzero jacobian structure on first solve...\n\n\n" );

  MatSetOption(matrix.mat(),MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);
  // indicates when MatZeroRows() is called the zeroed entries
  // are kept in the nonzero structure
  MatSetOption(matrix.mat(),MAT_KEEP_NONZERO_PATTERN    ,PETSC_TRUE);

  PetscFunctionReturnVoid();
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
/* ------------------------------------------------------------ */
#undef __FUNCT__
#define __FUNCT__ "PetscFEMSystem::build_context"
AutoPtr<DiffContext> PetscFEMSystem::build_context ()
{
  AutoPtr<DiffContext> ap(new PetscFEMContext(*this));

  ap->set_deltat_pointer( &deltat );
  
  return ap;
}

void PetscFEMSystem::init_data ()
{
  // Use Solver to handle time stepping directly
  this->time_solver =
    AutoPtr<TimeSolver>(new SteadySolver(*this));
 
  // set nonlinear solver
  this->time_solver->diff_solver() =
        AutoPtr<DiffSolver>(new PetscDiffSolver(*this));

  // initialize parent data
  Parent::init_data();

  // echo data
  this->printSelf(std::cout); 

  // set the file output name to write initial element variables
  OStringStream file_name;
  // initialize data buffers
  std::vector<Number>      soln;

}

// boundary conditions
bool PetscFEMSystem::side_time_derivative(bool request_jacobian,DiffContext &context)
{
  PetscFEMContext &c = libmesh_cast_ref<PetscFEMContext&>(context);

  // get the number of variable in the state system
  const unsigned int n_vars = this->n_vars();
  const DofMap& dof_map = this->get_dof_map();

  // reference to equation systems
  EquationSystems & es = this->get_equation_systems();

  // set the field id in the spatially varying data structures
  // Initialize the per-element data for elem.
  std::vector<unsigned int> param_dof_indices;
  libMesh::System &template_parameter_system = 
                      this->get_equation_systems().get_system("k_0");
  template_parameter_system.get_dof_map().dof_indices (c.elem, param_dof_indices);

  // current_solution calls map_global_to_local_index that will map
  // this global index to the local current_local_solution
  const unsigned int field_id = param_dof_indices[0];

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  unsigned int n_sidepoints = c.side_qrule->n_points();

  for (unsigned int qp=0; qp != n_sidepoints; qp++)
    {
      // Compute the solution at the theta timestep
      this->UpdateBoundaryData(c,qp);

      // loop over variables and compute diagonal terms only
      // derive classes call this for off diag terms 
      for( unsigned int var = 0 ; var < n_vars ; var++)
        { 
        // Compute the solution at the theta timestep
        Number u_theta  = c.side_theta_value(var, qp);

        // residual and jacobian for this variable
        DenseSubMatrix<Number> &Kuu = *c.elem_subjacobians[var][var];
        DenseSubVector<Number> &Fu = *c.elem_subresiduals[var];

        // Element Jacobian * quadrature weight for side integration
        const std::vector<Real> &JxW_side = c.side_fe_var[var]->get_JxW();
      
        // The velocity shape functions at side quadrature points.
        const std::vector<std::vector<Real> >& phi_side =
          c.side_fe_var[var]->get_phi();
      
        // Physical location of the quadrature points on the side
        const std::vector<Point>& qpoint = c.side_fe_var[var]->get_xyz();
      
        // The number of local degrees of freedom in u and v
        const unsigned int n_dofs = c.dof_indices_var[var].size(); 

        // bc_id = 3 ==> cauchy
        // bc_id = 2 ==> neumann
        short int bc_id =
          this->get_mesh().boundary_info->boundary_id(c.elem, c.side);
        libmesh_assert (bc_id != BoundaryInfo::invalid_id);

        for (unsigned int i=0; i != n_dofs; i++)
          {
            Fu(i) += JxW_side[qp] * phi_side[i][qp] *
                     // add bc to residual
                     CALL_MEMBER_FN(this, this->ResidualBC[bc_id][var])(var,u_theta);

            if (request_jacobian)
              for (unsigned int j=0; j != n_dofs; j++)
                {
                  Kuu(i,j) += JxW_side[qp] * phi_side[i][qp] * phi_side[j][qp] *
                              // add bc to jacobian
                     CALL_MEMBER_FN(this, this->JacobianBC[bc_id][var])(var) 
		          * c.ThetaValue(var) ;
                }
          }
      } // end loop over variables
    } // end quadrature loop

  return request_jacobian;
}
// wrapper for default assembly
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

