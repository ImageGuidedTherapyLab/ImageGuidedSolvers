/* ------------------------------------------------------------ */
// ThermalTherapySystem implementation
template< typename BioheatTransferModel  >
ThermalTherapySystem<BioheatTransferModel>::
ThermalTherapySystem (EquationSystems& es, const std::string& name,
					   const unsigned int number) :
  PetscFEMSystem                 (es, name, number),
  // initialize constitutive data
  m_MathModel        (*(es.parameters.get<GetPot*>("controlfile")), es) 
{
  NumericVector<Number> & GlobalVec  = this->add_vector("old_global_solution",true, PARALLEL);
  #ifdef LIBMESH_ENABLE_GHOSTED
  NumericVector<Number> & StorageVec = this->add_vector( "old_local_solution",true, GHOSTED);
  #else
  NumericVector<Number> & StorageVec = this->add_vector( "old_local_solution",true, SERIAL);
  #endif
  // make sure the same
  m_LinearSolve =  m_MathModel.LinearPDE();
}


/* ------------------------------------------------------------ */
#undef __FUNCT__
#define __FUNCT__ "ThermalTherapySystem::build_context"
template< typename BioheatTransferModel  >
AutoPtr<DiffContext> ThermalTherapySystem< BioheatTransferModel >::
build_context ()
{
  AutoPtr<DiffContext> ap(new PetscFEMContext(*this));

  ap->set_deltat_pointer( &deltat );
  
  return ap;
}



/* ------------------------------------------------------------ */
template< typename BioheatTransferModel  >
void ThermalTherapySystem< BioheatTransferModel >::
SetupDirichlet(libMesh::MeshBase& mesh)
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




template< typename BioheatTransferModel  >
void ThermalTherapySystem<BioheatTransferModel>::
init_data ()
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


/* -------------------------------------------------------------------- 
   Here we define the initialization routine for the
   Convection-Diffusion system.  This routine is
   responsible for applying the initial conditions to
   the system.
   -------------------------------------------------------------------- */ 
#undef __FUNCT__
#define __FUNCT__ "ThermalTherapySystem::user_initialization"
template< typename BioheatTransferModel  >
void ThermalTherapySystem< BioheatTransferModel >::
SetupInitialConditions ()
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
  if( this->m_MathModel.InitValues.size() != this->n_vars() ) ICerror = true;

  // error check
  if(ICerror)
   {
     std::cout << "Initial Conditions Improperly Setup" << std::endl;
     std::cout << "found # IC = " << this->m_MathModel.InitValues.size()  << std::endl;
     std::cout << "Expected # System Variables = " << this->n_vars()  
               << " Expected # Domains = "         
               << this->m_MathModel.get_num_elem_blk()  << std::endl;
     libmesh_error();
   }

  // scatter global parameters to local copy
  this->ScatterParametersLocally();

  // loop over subdomains so that the last subdomain has the final
  // overwrite of the nodal temperature
  for ( unsigned int id_sub = 0 ; 
        id_sub < static_cast<unsigned int>(this->m_MathModel.get_num_elem_blk()); ++id_sub)
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
             typename libMesh::System &template_parameter_system = 
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
                 Ue(current_dof) = CALL_MEMBER_FN_W_REF(this->m_MathModel,
                                         this->m_MathModel.InitValues[var])
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

// boundary conditions
template< typename BioheatTransferModel  >
bool ThermalTherapySystem< BioheatTransferModel >::
side_time_derivative(bool request_jacobian,DiffContext &context)
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
  typename libMesh::System &template_parameter_system = 
                      this->get_equation_systems().get_system("k_0");
  template_parameter_system.get_dof_map().dof_indices (c.elem, param_dof_indices);

  // current_solution calls map_global_to_local_index that will map
  // this global index to the local current_local_solution
  const unsigned int field_id = param_dof_indices[0];

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  unsigned int n_sidepoints = c.side_qrule->n_points();

  // loop over variables and assume uncoupled
  for( unsigned int var = 0 ; var < n_vars ; var++)
    { 
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

    for (unsigned int qp=0; qp != n_sidepoints; qp++)
      {
        // Compute the solution at the theta timestep
        Number u_theta  = c.side_theta_value(var, qp);

        // bc_id = 3 ==> cauchy
        // bc_id = 2 ==> neumann
        short int bc_id =
          this->get_mesh().boundary_info->boundary_id(c.elem, c.side);
        libmesh_assert (bc_id != BoundaryInfo::invalid_id);

        for (unsigned int i=0; i != n_dofs; i++)
          {
            Fu(i) += JxW_side[qp] * phi_side[i][qp] *
                     // add bc to residual
        CALL_MEMBER_FN_W_REF(this->m_MathModel,
                             this->m_MathModel.ResidualBC[bc_id])(var,u_theta, 
                                                                  field_id , this->m_PowerID,
                                                                  qpoint[qp], es.parameters) ;

            if (request_jacobian)
              for (unsigned int j=0; j != n_dofs; j++)
                {
                  Kuu(i,j) += JxW_side[qp] * phi_side[i][qp] * phi_side[j][qp] *
                              // add bc to jacobian
        CALL_MEMBER_FN_W_REF(this->m_MathModel,
                             this->m_MathModel.JacobianBC[bc_id])(var, field_id, 
                                                                  qpoint[qp], es.parameters) 
		          * c.ThetaValue(var) ;
                }
          }
      } // end quadrature loop
    } // end loop over variables

  return request_jacobian;
}

