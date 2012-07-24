// Constructor
template < typename MathematicalModel >
PennesInverseSystem< MathematicalModel > 
::PennesInverseSystem(EquationSystems& es,
             const std::string& name,
             const unsigned int number ) 
//base class constructor
:LITTSystem< MathematicalModel >::LITTSystem(es, name, number) 
{
}
/**
 * assemble adjoint matrix and rhs
 */
#undef __FUNCT__
#define __FUNCT__ "PennesInverseSystem::assemble_adjoint"
template< typename MathematicalModel  >
void PennesInverseSystem< MathematicalModel >::
assemble_adjoint(EquationSystems& es)
{
 // application context
 AppSolve *user    = es.parameters.get<AppSolve*>("AppSolve")  ;
 qoiBaseClass*  qoiOptimizer = user->qoiOptimizer();  // optimization context

 PetscFunctionBegin;

 std::cout << "assemble_adjoint not ready " << std::endl;
 libmesh_error();
//  PetscLogEventBegin(AppSolve::logevents[20],0,0,0,0); // init libMesh
// 
//  // Get a reference to the NonlinearImplicitSystem we are solving
//  TransientFEMSystem& state_system = 
//           es.get_system<TransientFEMSystem>("StateSystem");
// 
//  // Get a reference to the LinearImplicitSystem for adjoint problem
//  TransientFEMSystem& adjoint_system = 
//            es.get_system<TransientFEMSystem>("AdjointSystem");
// 
//  // get the number of variables in each system
//  const unsigned int n_vars_state = state_system.n_vars();
//  const unsigned int n_vars_adjnt = adjoint_system.n_vars();
//  libmesh_assert(n_vars_adjnt == n_vars_state ); // should be the same
// 
//  // Get a constant reference to the mesh object.
//  const libMesh::MeshBase& mesh = es.get_mesh();
// 
//  // The dimension that we are running
//  const unsigned int dim = mesh.mesh_dimension();
// 
//  // A reference to the \p DofMap object for this system.  The \p DofMap
//  // object handles the index translation from node and element numbers
//  // to degree of freedom numbers.  We will talk more about the \p DofMap
//  // in future examples.
//  const DofMap& dof_map_adjnt = adjoint_system.get_dof_map();
//  const DofMap& dof_map_state = state_system.get_dof_map();
// 
//  // create containers to hold vector of FEM data structures
//  std::vector<FEType > fe_type_state(n_vars_state); // FEM type
//  std::vector<FEType > fe_type_adjnt(n_vars_adjnt); // FEM type
//  std::vector<FEBase*> fe_state(     n_vars_state); // FEM object 
//  std::vector<FEBase*> fe_adjnt(     n_vars_adjnt); // FEM object 
//  std::vector<FEBase*> fe_face_state(n_vars_state); // FEM object for BC
//  std::vector<FEBase*> fe_face_adjnt(n_vars_adjnt); // FEM object for BC
// 
//  // Get the Finite Element type for each variable
//  for(unsigned int i_var = 0 ; i_var < n_vars_state ; i_var++)
//       fe_type_state.at(i_var) = state_system.variable_type(i_var);
//  for(unsigned int i_var = 0 ; i_var < n_vars_adjnt ; i_var++)
//       fe_type_adjnt.at(i_var) = adjoint_system.variable_type(i_var);
// 
//  // Build the Gauss quadrature rules for numerical integration.  Let the \p
//  // FEType object decide what order rule is appropriate.  Assumes that the
//  // first variable is the highest poly order element to choose the quadrature
//  // rule. Boundary integration requires one quadraure rule, with dimensionality
//  // one less than the dimensionality of the element.
//  QGauss qrule(dim  , fe_type_adjnt.at(0).default_quadrature_order() );
//  QGauss qface(dim-1, fe_type_adjnt.at(0).default_quadrature_order() );
// 
//  // state
//  for( unsigned int i_var = 0 ; i_var < n_vars_state ; i_var++)
//    {
//     // Build a Finite Element object of the specified type for ith variable
//     fe_state.at(i_var) = FEBase::build(dim, fe_type_state.at(i_var) ).release();
//    
//     // Tell the finite element objects to use our quadrature rule.
//     // Use the same quadrature rule for both fe types
//     fe_state.at(i_var)->attach_quadrature_rule (&qrule);
// 
//     // Finite Element object of the specified type for ith variable on boundary
//     fe_face_state.at(i_var)=FEBase::build(dim,fe_type_state.at(i_var)).release();
// 
//     // attach the boundary quadrature rule
//     fe_face_state.at(i_var)->attach_quadrature_rule (&qface); 
// 
//     // shape functions and derivatives
//     this->m_MathModel.phi.at(i_var)  = &fe_state.at(i_var)->get_phi();
//     this->m_MathModel.dphi.at(i_var) = &fe_state.at(i_var)->get_dphi();
//    }
// 
//  // adjoint
//  for( unsigned int i_var = 0 ; i_var < n_vars_adjnt ; i_var++)
//    {
//     // Build a Finite Element object of the specified type for ith variable
//     fe_adjnt.at(i_var) =  FEBase::build(dim, fe_type_adjnt.at(i_var) ).release();
//    
//     // Tell the finite element objects to use our quadrature rule.
//     // Use the same quadrature rule for both fe types
//     fe_adjnt.at(i_var)->attach_quadrature_rule (&qrule);
// 
//     // Finite Element object of the specified type for ith variable on boundary
//     fe_face_adjnt.at(i_var)=FEBase::build(dim,fe_type_adjnt.at(i_var)).release();
// 
//     // attach the boundary quadrature rule
//     fe_face_adjnt.at(i_var)->attach_quadrature_rule (&qface); 
// 
//     // shape functions and derivatives
//     this->m_MathModel.psi.at(i_var)  = &fe_adjnt.at(i_var)->get_phi();
//     this->m_MathModel.dpsi.at(i_var) = &fe_adjnt.at(i_var)->get_dphi();
//    }
// 
//  // Here we define some references to cell-specific data that
//  // will be used to assemble the linear system.
//  // The element Jacobian * quadrature weight at each integration point.   
//  const std::vector<Real>& JxW      = fe_adjnt.at(0)->get_JxW();
//  const std::vector<Real>& JxW_face = fe_face_adjnt.at(0)->get_JxW();
// 
//  // The physical XY locations of the quadrature points on the element.
//  // These are useful for evaluating spatially varying material
//  // properties at the quadrature points.
//  const std::vector<Point>& q_point     = fe_adjnt.at(0)->get_xyz();
//  const std::vector<Point>& qpoint_face = fe_face_adjnt.at(0)->get_xyz();
// 
//  // Define data structures to contain the element matrix
//  // and right-hand-side vector contribution.  Following
//  // basic finite element terminology we will denote these
//  // "Ke" and "Fe".
//  DenseMatrix<Number> Ke;
//  DenseVector<Number> Fe;
//  // n_vars_adjnt x n_vars_adjnt submatrices
//  std::vector< SubMatrixVector > Kij(n_vars_adjnt, 
//   SubMatrixVector::vector(n_vars_adjnt, DenseSubMatrix<Number>::DenseSubMatrix(Ke)) );
// 
//  //Define the vectors holding the residual sub vectors
//  std::vector< DenseSubVector<Number> > Fi(n_vars_adjnt,DenseSubVector<Number>::DenseSubVector(Fe));
// 
//  // This vector will hold the degree of freedom indices for
//  // the element.  These define where in the global system
//  // the element degrees of freedom get mapped.
//  std::vector< unsigned int >                 dof_indices_state;
//  std::vector< unsigned int >                 dof_indices_adjnt;
// 
//  // space for parameter mappings
//  std::vector<PetscScalar> elemParameter(qoiOptimizer->GetParamSize(),0.0);
// 
//  PetscLogEventEnd(  AppSolve::logevents[20],0,0,0,0); // init libMesh
// 
//  // Now we will loop over all the elements in the mesh.
//  // We will compute the element Jacobian contribution.
//  libMesh::MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
//  const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
// 
//  for ( ; el != end_el; ++el)
//    {
// 
//     PetscLogEventBegin(AppSolve::logevents[7],0,0,0,0); // elemfnc setup
// 
//     // Store a pointer to the element we are currently
//     // working on.  This allows for nicer syntax later.
//     const Elem* elem = *el;
// 
//     // set the field id in the spatially varying data structures
//     const unsigned int FieldId =  elem->_mat_data_field;
// 
//     // Get the degree of freedom indices for the
//     // current element.  These define where in the global
//     // matrix and right-hand-side this element will
//     // contribute to.
//     dof_map_adjnt.dof_indices (elem, dof_indices_adjnt);
//     dof_map_state.dof_indices (elem, dof_indices_state);
//     const unsigned int n_dofs_adjnt = dof_indices_adjnt.size();
//     const unsigned int n_dofs_state = dof_indices_state.size();
// 
//     // Zero the element matrix and right-hand side before
//     // summing them.  We use the resize member here because
//     // the number of degrees of freedom might have changed from
//     // the last element.  Note that this will be the case if the
//     // element type is different (i.e. the last element was a
//     // triangle, now we are on a quadrilateral).
//     Ke.resize (n_dofs_adjnt, n_dofs_adjnt);
//     Fe.resize (n_dofs_adjnt);
// 
//     // state
//     for( unsigned int i_var = 0 ; i_var < n_vars_state ; i_var++)
//       { 
//        // degree of freedom indices for individual components
//        dof_map_state.dof_indices (elem, this->m_MathModel.dof_indices_u[i_var],i_var);
//        this->m_MathModel.n_u_dofs[i_var] =     this->m_MathModel.dof_indices_u[i_var].size(); 
// 
//        // Compute the element-specific data for the current
//        // element.  This involves computing the location of the
//        // quadrature points (q_point) and the shape functions
//        // (phi, dphi), (psi,dpsi),... for the current element.
//        fe_state[i_var]->reinit (elem);
//       } 
// 
//     // adjoint
//     for( unsigned int i_var = 0 ; i_var < n_vars_adjnt ; i_var++)
//       { 
//        // degree of freedom indices for individual components
//        dof_map_adjnt.dof_indices (elem, this->m_MathModel.dof_indices_p[i_var],i_var);
//        this->m_MathModel.n_p_dofs[i_var] =     this->m_MathModel.dof_indices_p[i_var].size(); 
// 
//        // Compute the element-specific data for the current
//        // element.  This involves computing the location of the
//        // quadrature points (q_point) and the shape functions
//        // (phi, dphi), (psi,dpsi),... for the current element.
//        fe_adjnt[i_var]->reinit (elem);
// 
//        // Reposition the submatrices and subvectors...  The idea is this:
//        //
//        //         -   -   -                       -  -
//        //        | J11 J12  . |                  | F1 |
//        //   Ke = | J21 J22  . |             Fe = | F2 |
//        //        | J31 J32  . |                  |  . |
//        //        |  .   .   . |                  |  . |
//        //         -   -   -                      |  . |
//        //                                         -  -  
// 
//        // The \p DenseSubVector.reposition () member
//        // takes the (row_offset, row_size)
//        Fi[i_var].reposition (i_var*this->m_MathModel.n_p_dofs[0],
//        this->m_MathModel.n_p_dofs[i_var]);
//       }
// 
//        // The \p DenseSubMatrix.repostition () member takes the
//        // (row_offset, column_offset, row_size, column_size).
//     for( unsigned int i_var = 0 ; i_var < n_vars_adjnt ; i_var++)
//        for( unsigned int j_var = 0 ; j_var < n_vars_adjnt ; j_var++)
//          Kij[i_var][j_var].reposition(i_var*this->m_MathModel.n_p_dofs[i_var],j_var*this->m_MathModel.n_p_dofs[i_var],
//                                             this->m_MathModel.n_p_dofs[i_var],
//          this->m_MathModel.n_p_dofs[j_var]);
// 
//     PetscLogEventEnd(  AppSolve::logevents[7],0,0,0,0); // elemfnc setup
// 
//     PetscLogEventBegin(AppSolve::logevents[8],0,0,0,0); // elemfnc assble
// 
//     // Now we will build the element matrix and rhw.  This involves
//     // a double loop to integrate the test funcions (i) against
//     // the trial functions (j).
//     //
//     // accumulate the adjoint matrix and rhs inherent to the pde
//     this->m_MathModel.accumulateAdjointPDE(qrule,FieldId,JxW,Kij,Fi,
//                                                state_system,adjoint_system);
// 
//     // get mappings
//     qoiOptimizer->getGlobalParamValues(elem,elemParameter);
// 
//     // add QOI contributions to the rhs
//     CALL_MEMBER_FN(qoiOptimizer,qoiOptimizer->accumulateAdjointLoad)(user,qrule,
//                                          elemParameter, FieldId,JxW,q_point,Fi);
// 
//     // At this point the interior element integration has
//     // been completed.  However, we have not yet addressed
//     // boundary conditions.
//       
//     PetscLogEventEnd(  AppSolve::logevents[8],0,0,0,0); // elemfnc assble
// 
//     PetscLogEventBegin(AppSolve::logevents[9],0,0,0,0); // elemfnc bndry
// 
//     //The interior element integration has now been completed.
//     //Here we begin dealing with boundary conditions.
//     //The following loops over the sides of the element.
//     //If the element has no neighbor on a side then that
//     //side MUST live on a boundary of the domain.
//     //This assume boundary conditions are UNCOUPLED.
//     for (unsigned int side=0; side<elem->n_sides(); side++)
//      if (elem->neighbor(side) == NULL)
//       {
//        // Get the boundary ID for side 's'.
//        // These are set internally by build_square() or
//        // by the external mesh generation package, i.e. 
//        // this retrieves the side set id in the exodus file
//        short int bc_id = mesh.boundary_info->boundary_id (elem,side);
//        
//        // loop over variables assuming boundary condition variables are
//        // uncoupled
//        for( unsigned int i_var = 0 ; i_var < n_vars_adjnt ; i_var++)
//          {
//           // Compute the shape function values on the element face.
//           fe_face_adjnt[i_var]->reinit(elem, side);
// 
//           // The value of the shape functions at the quadrature points.
//           const std::vector<std::vector<Real> >&  psi_face = 
//                                                fe_face_adjnt[i_var]->get_phi();
// 
//           // add bc to adjoint matrix (utilizing symmetry in BC)
//           CALL_MEMBER_FN_W_REF(this->m_MathModel,this->m_MathModel.accumulateJacobianBC[bc_id])
//                                       (i_var,qface,JxW_face,psi_face,Kij);
//           // add bc to adjoint matrix
//           CALL_MEMBER_FN_W_REF(this->m_MathModel,this->m_MathModel.accumulateAdjointLoadBC[bc_id])
//                        (i_var,qface,psi_face,JxW_face, Fi,adjoint_system);
// 
//          } //  end loop over variables
//       } //end BC
//     
//     PetscLogEventEnd(  AppSolve::logevents[9],0,0,0,0); // elemfnc bndry
// 
//     PetscLogEventBegin(AppSolve::logevents[21],0,0,0,0); // setup libMesh
// 
//     dof_map_adjnt.constrain_element_matrix (Ke, dof_indices_adjnt);
//     dof_map_adjnt.constrain_element_vector (Fe, dof_indices_adjnt);
//     // The element matrix and right-hand-side are now built
//     // for this element.  Add them to the global matrix and
//     // right-hand-side vector.  The \p PetscMatrix::add_matrix()
//     // and \p PetscVector::add_vector() members do this for us.
//     adjoint_system.matrix->add_matrix (Ke, dof_indices_adjnt);
//     adjoint_system.rhs->add_vector    (Fe, dof_indices_adjnt);
// 
//     PetscLogEventEnd(  AppSolve::logevents[21],0,0,0,0); // setup libMesh
// 
//    }//end loop over elements
// 
//  // apply dirichlet data if any
//  if(this->m_dirichletNodes.size())
//    {
//      adjoint_system.matrix->close();
// #if PETSC_VERSION_LESS_THAN(3,0,0)
//      adjoint_system.matrix->zero_rows((std::vector<unsigned int>&)this->m_dirichletNodes,1.0);
// #else
//      adjoint_system.matrix->zero_rows(this->m_dirichletNodes,1.0);
// #endif
//      adjoint_system.rhs->close();
//      for( unsigned int Ii = 0; Ii<this->m_dirichletNodes.size();Ii++) 
//                     adjoint_system.rhs->set(this->m_dirichletNodes[Ii],0.0);
//    }
// 
//  // free up memory
//  for( unsigned int i_var = 0 ; i_var < n_vars_state ; i_var++)
//    {
//     delete fe_state.at(i_var) ;
//     delete fe_face_state.at(i_var) ;
//    }
//  for( unsigned int i_var = 0 ; i_var < n_vars_adjnt ; i_var++)
//    {
//     delete fe_adjnt.at(i_var) ;
//     delete fe_face_adjnt.at(i_var) ;
//    }
// 
//  // That's it.

 PetscFunctionReturnVoid();

}

#undef __FUNCT__
#define __FUNCT__ "FormGradient"
template< typename MathematicalModel  >
PetscErrorCode PennesInverseSystem< MathematicalModel >::
FormGradient(TAO_APPLICATION taoapp,Vec X, Vec G)
{
 EquationSystems &es = this->get_equation_systems();// equation system solver
 AppSolve* user    = es.parameters.get<AppSolve*>("AppSolve")  ;
 qoiBaseClass*    qoiOptimizer = user->qoiOptimizer();  // optimization context

 PetscScalar *varloc,*g,*gradloc,deltat,invdeltat,error_estimate;
 KSPConvergedReason reason;
 MPI_Status status;
 PetscErrorCode info;

 std::cout << "FormGradient not ready " << std::endl;
 libmesh_error();
 PetscFunctionBegin;
//
// PetscLogStagePush(AppSolve::logstages[3]);// grad evaluation
//
// PetscLogEventBegin(AppSolve::logevents[22],0,0,0,0); // init variable
//
// // Get a reference to the NonlinearImplicitSystem we are solving
// PennesInverseSystem& state_system = 
//          es.get_system<PennesInverseSystem>("StateSystem");
//
// // Get a reference to the LinearImplicitSystem for adjoint problem
// TransientFEMSystem& adjoint_system = 
//          es.get_system<TransientFEMSystem>("AdjointSystem");
//
// // get the number of variable in the state system
// const unsigned int n_vars_state = state_system.n_vars();
// const unsigned int n_vars_adjnt = adjoint_system.n_vars();
//
// // Get a constant reference to the mesh object.
// const libMesh::MeshBase& mesh = es.get_mesh();
//
// // The dimension that we are running
// const unsigned int dim = mesh.mesh_dimension();
//
// // A reference to the \p DofMap object for this system.  The \p DofMap
// // object handles the index translation from node and element numbers
// // to degree of freedom numbers.  We will talk more about the \p DofMap
// // in future examples.
// const DofMap& dof_map_adjnt = adjoint_system.get_dof_map();
// const DofMap& dof_map_state = state_system.get_dof_map();
//
// /* dEBUGGING
// info = VecGetArray(user->xloc_ADJ,&varloc); CHKERRQ(info);
// FORTRAN_NAME(initadjointfield)(varloc,&user->NDOF_ADJOINT[2],&nzero);
// info = VecRestoreArray(user->xloc_ADJ,&varloc); CHKERRQ(info);
// // local to global scatter
// info =  VecScatterBegin(user->xloc_ADJ,user->p,INSERT_VALUES,SCATTER_FORWARD,
//       	  locscat); CHKERRQ(info);
// info =  VecScatterEnd(user->xloc_ADJ,user->p,INSERT_VALUES,SCATTER_FORWARD,
//       	  locscat); CHKERRQ(info);
// //AppSolve::ISTEP = 1; 
//   FORTRAN_NAME(output_ascii)(&nadjoint,&user->IDOPT,&AppSolve::ISTEP); 
// // end DEBUGGING */
//
// // set local gradient to zero and get the pointer to the local gradient
// info = VecSet(G,0.0); CHKERRQ(info);
//
// // zero out the error estimate in time for this optimization step
// for(PetscInt iii= user->IdealNzero() +1; 
//              iii <= user->IdealNtime() ; iii++)
//                                       qoiOptimizer->TIMEERR.at(iii)= 0.0;
//
// // Define data structures to contain the element Gradient
// DenseVector<Number> Grad;
//
// // optimization parameter map to global data structures
// std::vector<PetscInt> globmap(qoiOptimizer->GetParamSize(),0);
//
// PetscLogEventEnd(AppSolve::logevents[22],0,0,0,0); // init variable
//
// //step backward in time, Set IC
// adjoint_system.solution->zero();                  // parallel
// adjoint_system.current_local_solution->zero();    // local
// for(AppSolve::ISTEP = qoiOptimizer->Nstephi() ; AppSolve::ISTEP > qoiOptimizer->Nsteplo() ;
//AppSolve::ISTEP--)
//  {
//    // Get a reference to the timestep to make code more readable 
//    const PetscInt    &ISTEP = AppSolve::ISTEP;
//
//    // update from previous time step
//    *adjoint_system.old_vector_solution[0] = 
//                          *adjoint_system.current_local_solution;
//
//    /* assemble_adjoint routines used w/ adjoint_system.solve()
//       !!!NOTE!!! Do not need to store the soln in data structures 
//                  after the solve b/c this was already done by the 
//                  last function evaluation */
//
//    PetscLogEventBegin(AppSolve::logevents[5],0,0,0,0); // libMesh Solve
//
//    // set QOI function pointer
//    qoiOptimizer->accumulateAdjointLoad = &qoiBaseClass::accumulateAdjointQOI;
//
//    // assemble adjoint system 
//    adjoint_system.matrix->zero();
//    adjoint_system.rhs->zero();
//    state_system.assemble_adjoint(es); //@todo {separate adjoint load & matrix}
//    PetscVector<Number> &solution =
//      *(libmesh_cast_ptr<PetscVector<Number>*>(adjoint_system.solution.get()));
//    PetscMatrix<Number> &matrix =
//      *(libmesh_cast_ptr<PetscMatrix<Number>*>(adjoint_system.matrix));
//    PetscVector<Number> &load =
//      *(libmesh_cast_ptr<PetscVector<Number>*>(adjoint_system.rhs));
//    adjoint_system.matrix->close();
//    adjoint_system.rhs->close();
//
//    // solve
//    PetscDiffSolver* adjointSolver = 
//                     libmesh_cast_ptr<PetscDiffSolver*>(
//                      & (*adjoint_system.time_solver->diff_solver()) );
//    KSP  snesksp;
//    SNESGetKSP(adjointSolver->snes() , &snesksp);
//    info = KSPSetOperators(snesksp,matrix.mat(),matrix.mat(),
//                           SAME_NONZERO_PATTERN ); CHKERRQ(info);
//    info = KSPSolve(snesksp ,load.vec(),solution.vec() );CHKERRQ(info);
//
//    // localize 
//    adjoint_system.solution->localize( *adjoint_system.current_local_solution );
//
//
//    //check convergence reason
//    KSPConvergedReason kspReason;
//    KSPGetConvergedReason(snesksp , &kspReason);
//    if(kspReason<0)
//      {
//       std::cerr<<"Adjoint System Diverged " <<std::endl; libmesh_error(); 
//      } 
//
//    // store the lagrange multiplier solution history in data structures 
//    *adjoint_system.vector_solution.at(AppSolve::ISTEP)=
//                                       *adjoint_system.current_local_solution;
//    PetscLogEventEnd(  AppSolve::logevents[5],0,0,0,0); // libMesh Solve
//
//
//    // This vector will hold the degree of freedom indices for
//    // the element.  These define where in the global system
//    // the element degrees of freedom get mapped.
//    std::vector< unsigned int >                 dof_indices_state;
//    std::vector< unsigned int >                 dof_indices_adjnt;
//
//    // create containers to hold vector of FEM data structures
//    std::vector<FEType > fe_type_state(n_vars_state); // FEM type
//    std::vector<FEType > fe_type_adjnt(n_vars_adjnt); // FEM type
//    std::vector<FEBase*> fe_state(     n_vars_state); // FEM object 
//    std::vector<FEBase*> fe_adjnt(     n_vars_adjnt); // FEM object 
//    std::vector<FEBase*> fe_face_state(n_vars_state); // FEM object for BC
//    std::vector<FEBase*> fe_face_adjnt(n_vars_adjnt); // FEM object for BC
//  
//    // Get the Finite Element type for each variable
//    for(unsigned int i_var = 0 ; i_var < n_vars_state ; i_var++)
//         fe_type_state.at(i_var) = state_system.variable_type(i_var);
//    for(unsigned int i_var = 0 ; i_var < n_vars_adjnt ; i_var++)
//         fe_type_adjnt.at(i_var) = adjoint_system.variable_type(i_var);
//  
//    // Build the Gauss quadrature rules for numerical integration.  Let the \p
//    // FEType object decide what order rule is appropriate.  Assumes that the
//    // first variable is the highest poly order element to choose the quadrature
//    // rule. Boundary integration requires one quadraure rule, with dimensionality
//    // one less than the dimensionality of the element.
//    QGauss qrule(dim  , fe_type_adjnt.at(0).default_quadrature_order() );
//    QGauss qface(dim-1, fe_type_adjnt.at(0).default_quadrature_order() );
//  
//    // state
//    for( unsigned int i_var = 0 ; i_var < n_vars_state ; i_var++)
//      {
//       // Build a Finite Element object of the specified type for ith variable
//       fe_state.at(i_var) = FEBase::build(dim, fe_type_state.at(i_var) ).release();
//      
//       // Tell the finite element objects to use our quadrature rule.
//       // Use the same quadrature rule for both fe types
//       fe_state.at(i_var)->attach_quadrature_rule (&qrule);
//  
//       // Finite Element object of the specified type for ith variable on boundary
//       fe_face_state.at(i_var)=FEBase::build(dim,fe_type_state.at(i_var)).release();
//  
//       // attach the boundary quadrature rule
//       fe_face_state.at(i_var)->attach_quadrature_rule (&qface); 
//  
//       // shape functions and derivatives
//       this->m_MathModel.phi.at(i_var)  = &fe_state.at(i_var)->get_phi();
//       this->m_MathModel.dphi.at(i_var) = &fe_state.at(i_var)->get_dphi();
//      }
//  
//    // adjoint
//    for( unsigned int i_var = 0 ; i_var < n_vars_adjnt ; i_var++)
//      {
//       // Build a Finite Element object of the specified type for ith variable
//       fe_adjnt.at(i_var) = FEBase::build(dim, fe_type_adjnt.at(i_var) ).release();
//      
//       // Tell the finite element objects to use our quadrature rule.
//       // Use the same quadrature rule for both fe types
//       fe_adjnt.at(i_var)->attach_quadrature_rule (&qrule);
//  
//       // Finite Element object of the specified type for ith variable on boundary
//       fe_face_adjnt.at(i_var)=FEBase::build(dim,fe_type_adjnt.at(i_var)).release();
//  
//       // attach the boundary quadrature rule
//       fe_face_adjnt.at(i_var)->attach_quadrature_rule (&qface); 
//  
//       // shape functions and derivatives
//       this->m_MathModel.psi.at(i_var)  = &fe_adjnt.at(i_var)->get_phi();
//       this->m_MathModel.dpsi.at(i_var) = &fe_adjnt.at(i_var)->get_dphi();
//      }
//  
//    // Here we define some references to cell-specific data that
//    // will be used to assemble the linear system.
//    // The element Jacobian * quadrature weight at each integration point.   
//    const std::vector<Real>& JxW      = fe_adjnt.at(0)->get_JxW();
//    const std::vector<Real>& JxW_face = fe_face_adjnt.at(0)->get_JxW();
//  
//    // The physical XY locations of the quadrature points on the element.
//    // These are useful for evaluating spatially varying material
//    // properties at the quadrature points.
//    const std::vector<Point>& q_point     = fe_adjnt.at(0)->get_xyz();
//    const std::vector<Point>& qpoint_face = fe_face_adjnt.at(0)->get_xyz();
//
//    // loop over all the elements in the mesh, compute the elements
//    // contribution to the Hessian, then map the elements contribution
//    // to the global hessian
//    libMesh::MeshBase::const_element_iterator       el    =mesh.active_local_elements_begin();
//    const libMesh::MeshBase::const_element_iterator end_el=mesh.active_local_elements_end();
//    for ( ; el != end_el; ++el)
//     {
//
//      PetscLogEventBegin(AppSolve::logevents[10],0,0,0,0); // elemjac setup
//
//      // Store a pointer to the element we are currently
//      // working on.  This allows for nicer syntax later.
//      const Elem* elem = *el;
//
//      // set the field id in the spatially varying data structures
//      const unsigned int FieldId =  elem->_mat_data_field;
//
//      // Get the degree of freedom indices for the
//      // current element.  These define where in the global
//      // matrix and right-hand-side this element will
//      // contribute to.
//      dof_map_adjnt.dof_indices (elem, dof_indices_adjnt);
//      dof_map_state.dof_indices (elem, dof_indices_state);
//      const unsigned int n_dofs_adjnt = dof_indices_adjnt.size();
//      const unsigned int n_dofs_state = dof_indices_state.size();
//
//      // Ensure the element gradient is the right size and zero out the entries
//      Grad.resize (qoiOptimizer->GetParamSize());
//
//      // state
//      for( unsigned int i_var = 0 ; i_var < n_vars_state ; i_var++)
//        { 
//         // degree of freedom indices for individual components
//         dof_map_state.dof_indices (elem,this->m_MathModel.dof_indices_u[i_var],i_var);
//         this->m_MathModel.n_u_dofs[i_var]  =   this->m_MathModel.dof_indices_u[i_var].size();
//
//         // Compute the element-specific data for the current
//         // element.  This involves computing the location of the
//         // quadrature points (q_point) and the shape functions
//         // (phi, dphi), (psi,dpsi),... for the current element.
//         fe_state[i_var]->reinit (elem);
//        } 
//
//      // adjoint
//      for( unsigned int i_var = 0 ; i_var < n_vars_adjnt ; i_var++)
//        { 
//         // degree of freedom indices for individual components
//         dof_map_adjnt.dof_indices (elem,this->m_MathModel.dof_indices_p[i_var],i_var);
//         this->m_MathModel.n_p_dofs[i_var] =    this->m_MathModel.dof_indices_p[i_var].size();
//
//         // Compute the element-specific data for the current
//         // element.  This involves computing the location of the
//         // quadrature points (q_point) and the shape functions
//         // (phi, dphi), (psi,dpsi),... for the current element.
//         fe_adjnt[i_var]->reinit (elem);
//        }
//
//      PetscLogEventEnd(  AppSolve::logevents[10],0,0,0,0); // elemjac setup
//
//      PetscLogEventBegin(AppSolve::logevents[11],0,0,0,0); // elemjac assble
//
//      // accumulate the gradient for each parameter
//      qoiOptimizer->accumulateGradient(user,qrule,FieldId,JxW,q_point,Grad);
//
//      PetscLogEventEnd(  AppSolve::logevents[11],0,0,0,0); // elemjac assble
//
//      PetscLogEventBegin(AppSolve::logevents[12],0,0,0,0); // elemjac bndry
//
//      /* The following gets the error estimate contritubution on the 
//         boundary. Optimization parameters are assumed to have 
//         NO BC contribution */
//      for (unsigned int side=0; side<elem->n_sides(); side++)
//       if (elem->neighbor(side) == NULL)
//        {
//	 // Get the boundary ID for side 's'.  These are set internally by
//	 // build_square() or by the external mesh generation package, i.e.
//	 // this retrieves the side set id in the exodus file
//         short int bc_id = mesh.boundary_info->boundary_id (elem,side);
//         
//         // loop over variables assuming boundary condition variables are
//         // uncoupled
//         for( unsigned int i_var = 0 ; i_var < n_vars_adjnt ; i_var++)
//           {
//            // Compute the shape function values on the element face.
//            fe_face_state[i_var]->reinit(elem, side);
//            fe_face_adjnt[i_var]->reinit(elem, side);
//
//            // The value of the shape functions at the quadrature points.
//            const std::vector<std::vector<Real> >&  phi_face = 
//                                              fe_face_state[i_var]->get_phi();
//            const std::vector<std::vector<Real> >&  psi_face = 
//                                              fe_face_adjnt[i_var]->get_phi();
//            // add bc to adjoint matrix
//            CALL_MEMBER_FN_W_REF(this->m_MathModel,this->m_MathModel.accumulateGradientBC[bc_id])
//                               (i_var,qface,phi_face,psi_face,JxW_face,
//                                            Grad,state_system,adjoint_system);
//           } //  end loop over variables
//        } //end BC
//
//       PetscLogEventEnd(  AppSolve::logevents[12],0,0,0,0); // elemjac bndry
//
//       // get parameter mapping
//       qoiOptimizer->getGlobalMap(elem,globmap);
//
//       // scale the gradient by deltat
//       Grad.scale(user->get_fem_dt());
//
//       // put the gradient in the tao data structures
//       info =  VecSetValues(G,qoiOptimizer->GetParamSize(),&globmap[0],
//       	              (PetscScalar*) &Grad.get_values()[0] ,ADD_VALUES);
//       CHKERRQ(info);
//
//       // store error estimate 
//       PetscInt  idpow;
//       if( AppSolve::ISTEP % user->IstepsPerIdeal() ) 
//          idpow = user->get_id_power() + 1;
//       else 
//          idpow = user->get_id_power() ;
//       #if defined(PETSC_USE_DEBUG)
//       qoiOptimizer->TIMEERR.at(idpow) = // accumulate error over time
//                qoiOptimizer->TIMEERR.at(idpow) + Grad(0);
//       #else
//       qoiOptimizer->TIMEERR[idpow] = // accumulate error over time
//                qoiOptimizer->TIMEERR[idpow] + Grad(0);
//       #endif
//     } // end loop over element gradient computation
//
//    // free up memory
//    for( unsigned int i_var = 0 ; i_var < n_vars_state ; i_var++)
//      {
//       delete fe_state.at(i_var) ;
//       delete fe_face_state.at(i_var) ;
//      }
//    for( unsigned int i_var = 0 ; i_var < n_vars_adjnt ; i_var++)
//      {
//       delete fe_adjnt.at(i_var) ;
//       delete fe_face_adjnt.at(i_var) ;
//      }
//
//    CHKMEMQ; // check for memory corruption use -malloc_debug to enable
//  } // end time step loop
//
//
// PetscLogEventBegin(AppSolve::logevents[0],0,0,0,0); // tao param xfer
//
// //Communicate Assembly of the Gradient
// info = VecAssemblyBegin(G); CHKERRQ(info);
// info = VecAssemblyEnd(G); CHKERRQ(info)
//
// // print gradient norm
// PetscScalar gradientNorm;
// info = VecNorm(G,NORM_2,&gradientNorm); CHKERRQ(info)
// PetscPrintf(PETSC_COMM_WORLD,"||Grad||_2 = %e \n",gradientNorm);
//
// // gather the error estimates from all processors
// PetscInt Nsteps = user->IdealNtime()  - user->IdealNzero() ;
// // create temporary memory buffer
// std::vector<PetscScalar> scratcherr(Nsteps);
// //all reduce over PETSC_COMM_WORLD on the number of elements
// if(Nsteps > 0) MPI_Allreduce(&qoiOptimizer->TIMEERR[user->IdealNzero() +1],
//                              &scratcherr[0],Nsteps, MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD);
//
// CHKMEMQ; // check for memory corruption use -malloc_debug to enable
// // store buffer in main data structures
// error_estimate= 0.0;
// for(PetscInt iii= user->IdealNzero()  + 1; iii <= user->IdealNtime() ; iii++)
//  {
//   qoiOptimizer->TIMEERR.at(iii) = scratcherr.at(iii-(user->IdealNzero() +1));
//   error_estimate = error_estimate + qoiOptimizer->TIMEERR.at(iii) ;
//  }
//
// PetscPrintf(PETSC_COMM_WORLD,"error estimate = %e \n", error_estimate);
//
// if( user->get_num_exact() ){info = qoiOptimizer->compareFDGrad(taoapp,X,G,user); CHKERRQ(info);}
//
// PetscLogEventEnd(AppSolve::logevents[0],0,0,0,0); // tao param xfer
//
// PetscLogStagePop();// grad evaluation

 PetscFunctionReturn(0);
}

/* --------------Sensitivity Load Computation Routine---------- */
#undef __FUNCT__
#define __FUNCT__ "PennesInverseSystem::sensitivityLoad"
template< typename MathematicalModel  >
void PennesInverseSystem< MathematicalModel >::
sensitivityLoad( NumericVector<Number>& rhs )
{

 PetscFunctionBegin; 

 std::cout << "sensitivityLoad not ready " << std::endl;
 libmesh_error();
// //PetscLogEventBegin(AppSolve::logevents[22],0,0,0,0); // init variable
//
// EquationSystems &es = this->get_equation_systems();// equation system solver
// AppSolve* user    = es.parameters.get<AppSolve*>("AppSolve")  ;
// qoiBaseClass*    qoiOptimizer = user->qoiOptimizer();  // optimization context
//
// // Get a reference to the NonlinearImplicitSystem we are solving
// TransientFEMSystem& state_system = 
//   es.get_system<TransientFEMSystem>("StateSystem");
//
// TransientFEMSystem& sensitivity_system = 
//   es.get_system<TransientFEMSystem>("SensitivitySystem");
//
// // get the number of variable in the state system
// const unsigned int n_vars = state_system.n_vars();
//
// // Get a constant reference to the mesh object.
// const libMesh::MeshBase& mesh = es.get_mesh();
//
// // The dimension that we are running
// const unsigned int dim = mesh.mesh_dimension();
// 
// // A reference to the \p DofMap object for this system.  The \p DofMap
// // object handles the index translation from node and element numbers
// // to degree of freedom numbers.  We will talk more about the \p DofMap
// // in future examples.
// const DofMap & dof_map = state_system.get_dof_map();
//
// // create containers to hold vector of FEM data structures
// std::vector<FEType > fe_type(n_vars); // FEM type
// std::vector<FEBase*> fe(     n_vars); // FEM object 
// std::vector<FEBase*> fe_face(n_vars); // FEM object for boundary integration
//
// // Get the Finite Element type for each variable
// for(unsigned int i_var = 0 ; i_var < n_vars ; i_var++)
//      fe_type.at(i_var) = state_system.variable_type(i_var);
//
// // Build the Gauss quadrature rules for numerical integration.  Let the \p
// // FEType object decide what order rule is appropriate.  Assumes that the
// // first variable is the highest poly order element to choose the quadrature
// // rule. Boundary integration requires one quadraure rule, with dimensionality
// // one less than the dimensionality of the element.
// QGauss qrule(dim  , fe_type.at(0).default_quadrature_order() );
// QGauss qface(dim-1, fe_type.at(0).default_quadrature_order() );
//
// for( unsigned int i_var = 0 ; i_var < n_vars ; i_var++)
//   {
//    // Build a Finite Element object of the specified type for ith variable
//    fe.at(i_var) =  FEBase::build(dim, fe_type.at(i_var) ).release();
//   
//    // Tell the finite element objects to use our quadrature rule.
//    // Use the same quadrature rule for both fe types
//    fe.at(i_var)->attach_quadrature_rule (&qrule);
//
//    // Finite Element object of the specified type for ith variable on boundary
//    fe_face.at(i_var) =  FEBase::build(dim, fe_type.at(i_var) ).release();
//
//    // attach the boundary quadrature rule
//    fe_face.at(i_var)->attach_quadrature_rule (&qface); 
//
//    // shape functions and derivatives
//    this->m_MathModel.phi.at(i_var)  = &fe.at(i_var)->get_phi();
//    this->m_MathModel.dphi.at(i_var) = &fe.at(i_var)->get_dphi();
//   }
//
// // Here we define some references to cell-specific data that
// // will be used to assemble the linear system.
// // The element Jacobian * quadrature weight at each integration point.   
// const std::vector<Real>& JxW = fe.at(0)->get_JxW();
// const std::vector<Real>& JxW_face = fe_face.at(0)->get_JxW();
//
// // The physical XY locations of the quadrature points on the element.
// // These are useful for evaluating spatially varying material
// // properties at the quadrature points.
// const std::vector<Point>& q_point = fe.at(0)->get_xyz();
// const std::vector<Point>& qpoint_face = fe_face.at(0)->get_xyz();
//
// //Define the vector holding the R values on an element.
// DenseVector<Number> Fe;
//
// //Define the vectors holding the sensitivity load sub vectors
// std::vector< DenseSubVector<Number> > Fi(n_vars,DenseSubVector<Number>::DenseSubVector(Fe));
//
// //tmp space for element storage
// std::vector<PetscScalar> elemParameter(qoiOptimizer->GetParamSize(),0.0);
//
// // This vector will hold the degree of freedom indices for
// // the element.  These define where in the global system
// // the element degrees of freedom get mapped.
// std::vector< unsigned int > dof_indices;
//
// //PetscLogEventEnd(  AppSolve::logevents[22],0,0,0,0); // init variable
//
// // Now we will loop over all the elements in the mesh that
// // live on the local processor. We will compute the element
// // matrix and right-hand-side contribution.  Since the mesh
// // will be refined we want to only consider the ACTIVE elements,
// // hence we use a variant of the \p active_elem_iterator.
// libMesh::MeshBase::const_element_iterator       el    =mesh.active_local_elements_begin();
// const libMesh::MeshBase::const_element_iterator end_el=mesh.active_local_elements_end();
// rhs.zero();
// for ( ; el != end_el; ++el)
//   {    
// 
//    //PetscLogEventBegin(AppSolve::logevents[7],0,0,0,0); // elemfnc setup
//
//    // Store a pointer to the element we are currently
//    // working on.  This allows for nicer syntax later.
//    Elem* elem = *el;
//    
//    // set the field id in the spatially varying data structures
//    const unsigned int FieldId =  elem->_mat_data_field;
//
//    // Get the degree of freedom indices for the
//    // current element.  These define where in the global
//    // matrix and right-hand-side this element will
//    // contribute to.
//    dof_map.dof_indices (elem, dof_indices);
//    const unsigned int n_dofs   = dof_indices.size();
//
//    // Zero rhs element vector before summing it. 
//    // We use the resize member here because
//    // the number of degrees of freedom might have changed from
//    // the last element.  Note that this will be the case if the
//    // element type is different (i.e. the last element was a
//    // triangle, now we are on a quadrilateral).
//    Fe.resize (n_dofs);
//
//    for( unsigned int i_var = 0 ; i_var < n_vars ; i_var++)
//      { 
//      // degree of freedom indices for individual components
//       dof_map.dof_indices (elem,   this->m_MathModel.dof_indices_u[i_var],i_var);
//       this->m_MathModel.n_u_dofs[i_var] = this->m_MathModel.dof_indices_u[i_var].size(); 
//
//       // Compute the element-specific data for the current
//       // element.  This involves computing the location of the
//       // quadrature points (q_point) and the shape functions
//       // (phi, dphi), (psi,dpsi), ... for the current element.
//       fe[i_var]->reinit (elem);
//
//       // Reposition the subvectors...  The idea is this:
//       //
//       //         -  -
//       //        | R1 |
//       //   Fe = | R2 |
//       //        |  . |
//       //        |  . |
//       //        |  . |
//       //         -  -
//       //
//       // The \p DenseSubVector.reposition () member
//       // takes the (row_offset, row_size)
//       Fi[i_var].reposition (i_var*this->m_MathModel.n_u_dofs[0],
//this->m_MathModel.n_u_dofs[i_var]);
//      }
//
//    //PetscLogEventEnd(  AppSolve::logevents[7],0,0,0,0); // elemfnc setup
//
//    PetscLogEventBegin(AppSolve::logevents[8],0,0,0,0); // elemfnc assble
//
//    // Now we will build the element sensivitity load vector.
//    // Constructing this requires the soln and its
//    // gradient from the previous timestep, and the sensitivity from the
//    // previous time step.  This must be
//    // calculated at each quadrature point by summing the
//    // soln degree-of-freedom values by the appropriate
//    // weight functions. Accumulate load from previous sensitivity
//    this->m_MathModel.accumulatePreviousSensitivity(qrule,FieldId,JxW,
//                                   q_point,Fi,state_system,sensitivity_system );
//
//    // get the parameter values associated with this element 
//    /* Store control variable data in this->m_MathModel.ata structures*/
//    qoiOptimizer->getGlobalParamValues(elem,elemParameter);
//    qoiOptimizer->accumulateSensitivityLoad(user,qrule,elemParameter,FieldId,
//                                            JxW,q_point,Fi,state_system);
//
//    PetscLogEventEnd(  AppSolve::logevents[8],0,0,0,0); // elemfnc assble
//
//    /* At this point the interior element integration has
//       been completed. No boundary condition optimization at this point.
//       Thus no BC computations are needed in the Sensitivity Solve. */ 
//
//    // The element sensitivity load is now built for this element.  
//    // We must now add them to the global rhs.  The 
//    // \p PetscVector::add_vector() members do this for us.
//    dof_map.constrain_element_vector (Fe, dof_indices);
//    rhs.add_vector (                  Fe, dof_indices);
//
//   } // end element loop...   for ( ; el != end_el; ++el)
//
// // apply dirichlet data if any
// if(this->m_dirichletNodes.size())
//   {
//     rhs.close();
//     for( unsigned int Ii = 0; Ii<this->m_dirichletNodes.size();Ii++) 
//                       rhs.set(this->m_dirichletNodes[Ii],0.0);
//   }
//
// // free up memory
// for( unsigned int i_var = 0 ; i_var < n_vars ; i_var++)
//   {
//    delete fe.at(i_var) ;
//    delete fe_face.at(i_var) ;
//   }
//
// // That's it.

 PetscFunctionReturnVoid(); 

}

#undef __FUNCT__
#define __FUNCT__ "PennesInverseSystem::element_time_derivative "
template< typename MathematicalModel  >
bool PennesInverseSystem< MathematicalModel >::
element_time_derivative (bool request_jacobian, DiffContext &context)
{
  this->m_MathModel.accumulateResidual(this->u_var,
                                       context,*this);
  // Matrix contributions for the uu and vv couplings.
  if (request_jacobian)
     this->m_MathModel.accumulateJacobian(this->u_var, 
                                          context,*this);
  
  return request_jacobian;
}

#undef __FUNCT__
#define __FUNCT__ "hessianVectorProduct"
template< typename MathematicalModel  >
PetscErrorCode PennesInverseSystem< MathematicalModel >::
hessianVectorProduct(Mat Hess, Vec xIn, Vec yOut)
{
  PetscFunctionBegin; 

 std::cout << "hessianVectorProduct not ready " << std::endl;
 libmesh_error();
//  PetscErrorCode info;
//  // Zero Entries first
//  info = VecSet(yOut,0.0);CHKERRQ(info);
//  info = VecAssemblyBegin(yOut); CHKERRQ(info);
//  info = VecAssemblyEnd(  yOut); CHKERRQ(info);
//
//  // get solver context 
//  EquationSystems &es = this->get_equation_systems();// equation system solver
//  AppSolve* user    = es.parameters.get<AppSolve*>("AppSolve")  ;
//  qoiBaseClass*    qoiOptimizer = user->qoiOptimizer();  // optimization context
//
//  PetscLogStagePush(AppSolve::logstages[4]);
//
//  // scatter parameter vector to local data structures
//  info =  VecScatterBegin(qoiOptimizer->GLOCSCAT_CNTL,xIn,
//                          qoiOptimizer->CNTL_LOC,
//                          INSERT_VALUES, SCATTER_REVERSE); CHKERRQ(info);
//  info =  VecScatterEnd(  qoiOptimizer->GLOCSCAT_CNTL,xIn,
//                          qoiOptimizer->CNTL_LOC,
//                          INSERT_VALUES, SCATTER_REVERSE); CHKERRQ(info);
//
//  // Get a reference to the NonlinearImplicitSystem we are solving
//  TransientFEMSystem& state_system = 
//          es.get_system<TransientFEMSystem>("StateSystem");
// 
//  // Get a reference to the LinearImplicitSystem for adjoint problem
//  TransientFEMSystem& adjoint_system = 
//           es.get_system<TransientFEMSystem>("AdjointSystem");
// 
//  // Get a reference to the LinearImplicitSystem for sensitivity problem
//  TransientFEMSystem & sensitivity_system =
//       es.get_system<TransientFEMSystem>("SensitivitySystem");
//
//  // reset IC
//  *state_system.current_local_solution = *state_system.vector_solution.at(qoiOptimizer->Nsteplo() ) ;
//  sensitivity_system.vector_solution.at(qoiOptimizer->Nsteplo() )->zero();  // local
//  sensitivity_system.current_local_solution->zero();    // local
//  sensitivity_system.solution->zero();                  // parallel
//  sensitivity_system.old_vector_solution.at(1)->zero(); // local
//  // compute time history for this linear combination of sensitivities 
//  for(AppSolve::ISTEP  = qoiOptimizer->Nsteplo() +1 ; 
//      AppSolve::ISTEP <= qoiOptimizer->Nstephi()    ; AppSolve::ISTEP++)
//   { // ensure state entries are correct
//     *state_system.old_vector_solution[0] =
//         *state_system.current_local_solution;
//     *state_system.current_local_solution=
//         *state_system.vector_solution.at(AppSolve::ISTEP) ;
//
//     // update previous solution
//     *sensitivity_system.old_vector_solution.at(0) = 
//         *sensitivity_system.current_local_solution;
//
//     /* sensitivity load obtained from linear combination of input vector */
//     qoiOptimizer->computeSensitivity( user , xIn );
//
//     //std::vector<Number>      solnVector;
//     //// get names from the system
//     //std::vector<std::string> names;
//     //const unsigned int nv_sys = sensitivity_system.n_vars();
//     //for (unsigned int var=0; var<nv_sys; var++)
//     //     names.push_back( sensitivity_system.variable_name(var) );
//     //
//     // // write every nth interval
//     // if( AppSolve::modWrite(AppSolve::ISTEP) 
//     //        && qoiOptimizer->PLOTOPTSOLVE     ) //write every nth step
//     //  {
//     //   /* build data buffer from current_local_solution */
//     //   build_system_solution_vector(es,sensitivity_system,solnVector);
//     //   sensitivity_vis->write_nodal_data(file_name.str(),
//     //                                    solnVector, names);
//     //   // write timestep info to file 
//     //   Real Time = AppSolve::getTime(AppSolve::ISTEP);
//     //   sensitivity_vis->write_timestep(Time);
//     //   // increment timestep (1-based) 
//     //   sensitivity_vis->increment_timestep();
//     //  }
//
//     // store sensitivity 
//     *sensitivity_system.vector_solution.at(AppSolve::ISTEP) =
//                        *sensitivity_system.current_local_solution;
//
//     // notice need to deference to call the copy constructor
//     //*sensitivity_system.old_vector_solution.at(0) = 
//     //    *sensitivity_system.vector_solution.at(AppSolve::ISTEP-1);
//     //*sensitivity_system.old_vector_solution.at(1) = 
//     //    *sensitivity_system.vector_solution.at(AppSolve::ISTEP  );
//
//     // debugging
//     //PetscVector<Number>* solnVec = 
//     //       libmesh_cast_ptr<PetscVector<Number>*>(
//     //       &(*sensitivity_system.old_vector_solution.at(1)) );
//     //std::cout << "SensitivityOldSoln" << jjj << std::endl;
//     //info = VecView(solnVec->vec(),0);CHKERRQ(info);
//   }
//
// // get the number of variable in the state system
// const unsigned int n_vars_state = state_system.n_vars();
// const unsigned int n_vars_adjnt = adjoint_system.n_vars();
//
// // Get a constant reference to the mesh object.
// const libMesh::MeshBase& mesh = es.get_mesh();
//
// // The dimension that we are running
// const unsigned int dim = mesh.mesh_dimension();
//
// // A reference to the \p DofMap object for this system.  The \p DofMap
// // object handles the index translation from node and element numbers
// // to degree of freedom numbers.  We will talk more about the \p DofMap
// // in future examples.
// const DofMap& dof_map_adjnt = adjoint_system.get_dof_map();
// const DofMap& dof_map_state = state_system.get_dof_map();
//
// // This vector will hold the degree of freedom indices for
// // the element.  These define where in the global system
// // the element degrees of freedom get mapped.
// std::vector<unsigned int> dof_indices;
//
// // optimization parameter map to global data structures
// std::vector<PetscInt>    globmap(qoiOptimizer->GetParamSize(),0);
// std::vector<PetscScalar> elemParameter(qoiOptimizer->GetParamSize(),0.0);
//
// // Define data structures to contain the element Hessian-Vector Product
// DenseVector<Number> HessVecProd;
//
// //accumulate sensitivity grad for verif
// #if defined(PETSC_USE_DEBUG)
// PetscScalar Gradiii =0.0;
// #endif
//
// /* FIXME: Pointer error when using nested loop through elements so must
//           separate sensitivites solves from Hessian accumulation        */
//
// //step backward in time, Set IC
// adjoint_system.solution->zero();                  // parallel
// adjoint_system.current_local_solution->zero();    // local
// for(AppSolve::ISTEP = qoiOptimizer->Nstephi() ; AppSolve::ISTEP > qoiOptimizer->Nsteplo() ;
//AppSolve::ISTEP--)
//  {
//    // Get a reference to the timestep to make code more readable 
//    const PetscInt    &ISTEP = AppSolve::ISTEP;
//
//    // update from previous time step
//    *adjoint_system.old_vector_solution[0] = 
//                          *adjoint_system.current_local_solution;
//
//    /* assemble_adjoint routines used w/ adjoint_system.solve()
//       !!!NOTE!!! Do not need to store the soln in data structures 
//                  after the solve b/c this was already done by the 
//                  last function evaluation */
//
//    PetscLogEventBegin(AppSolve::logevents[5],0,0,0,0); // libMesh Solve
//  
//    // set QOI function pointer
//    qoiOptimizer->accumulateAdjointLoad =
//                                    &qoiBaseClass::accumulateAdjointSensitivity;
//    adjoint_system.solve(); // solve
//    //print the reason of convergence
//    PetscDiffSolver* adjointSolver = 
//                     libmesh_cast_ptr<PetscDiffSolver*>(
//                      & (*adjoint_system.time_solver->diff_solver()) );
//    //check convergence reason
//    KSP  snesksp;
//    SNESGetKSP(adjointSolver->snes() , &snesksp);
//    KSPConvergedReason kspReason;
//    KSPGetConvergedReason(snesksp , &kspReason);
//    if(kspReason<0)
//      {
//       std::cerr<<"Adjoint Sensitivity Diverged " <<std::endl; libmesh_error(); 
//      } 
//
//    PetscLogEventEnd(  AppSolve::logevents[5],0,0,0,0); // libMesh Solve
//
//
//    // This vector will hold the degree of freedom indices for
//    // the element.  These define where in the global system
//    // the element degrees of freedom get mapped.
//    std::vector< unsigned int >                 dof_indices_state;
//    std::vector< unsigned int >                 dof_indices_adjnt;
//
//    // create containers to hold vector of FEM data structures
//    std::vector<FEType > fe_type_state(n_vars_state); // FEM type
//    std::vector<FEType > fe_type_adjnt(n_vars_adjnt); // FEM type
//    std::vector<FEBase*> fe_state(     n_vars_state); // FEM object 
//    std::vector<FEBase*> fe_adjnt(     n_vars_adjnt); // FEM object 
//    std::vector<FEBase*> fe_face_state(n_vars_state); // FEM object for BC
//    std::vector<FEBase*> fe_face_adjnt(n_vars_adjnt); // FEM object for BC
//  
//    // Get the Finite Element type for each variable
//    for(unsigned int i_var = 0 ; i_var < n_vars_state ; i_var++)
//         fe_type_state.at(i_var) = state_system.variable_type(i_var);
//    for(unsigned int i_var = 0 ; i_var < n_vars_adjnt ; i_var++)
//         fe_type_adjnt.at(i_var) = adjoint_system.variable_type(i_var);
//  
//    // Build the Gauss quadrature rules for numerical integration.  Let the \p
//    // FEType object decide what order rule is appropriate.  Assumes that the
//    // first variable is the highest poly order element to choose the quadrature
//    // rule. Boundary integration requires one quadraure rule, with dimensionality
//    // one less than the dimensionality of the element.
//    QGauss qrule(dim  , fe_type_adjnt.at(0).default_quadrature_order() );
//    QGauss qface(dim-1, fe_type_adjnt.at(0).default_quadrature_order() );
//  
//    // state
//    for( unsigned int i_var = 0 ; i_var < n_vars_state ; i_var++)
//      {
//       // Build a Finite Element object of the specified type for ith variable
//       fe_state.at(i_var) = FEBase::build(dim, fe_type_state.at(i_var) ).release();
//      
//       // Tell the finite element objects to use our quadrature rule.
//       // Use the same quadrature rule for both fe types
//       fe_state.at(i_var)->attach_quadrature_rule (&qrule);
//  
//       // Finite Element object of the specified type for ith variable on boundary
//       fe_face_state.at(i_var)=FEBase::build(dim,fe_type_state.at(i_var)).release();
//  
//       // attach the boundary quadrature rule
//       fe_face_state.at(i_var)->attach_quadrature_rule (&qface); 
//  
//       // shape functions and derivatives
//       this->m_MathModel.phi.at(i_var)  = &fe_state.at(i_var)->get_phi();
//       this->m_MathModel.dphi.at(i_var) = &fe_state.at(i_var)->get_dphi();
//      }
//  
//    // adjoint
//    for( unsigned int i_var = 0 ; i_var < n_vars_adjnt ; i_var++)
//      {
//       // Build a Finite Element object of the specified type for ith variable
//       fe_adjnt.at(i_var) = FEBase::build(dim, fe_type_adjnt.at(i_var) ).release();
//      
//       // Tell the finite element objects to use our quadrature rule.
//       // Use the same quadrature rule for both fe types
//       fe_adjnt.at(i_var)->attach_quadrature_rule (&qrule);
//  
//       // Finite Element object of the specified type for ith variable on boundary
//       fe_face_adjnt.at(i_var)=FEBase::build(dim,fe_type_adjnt.at(i_var)).release();
//  
//       // attach the boundary quadrature rule
//       fe_face_adjnt.at(i_var)->attach_quadrature_rule (&qface); 
//  
//       // shape functions and derivatives
//       this->m_MathModel.psi.at(i_var)  = &fe_adjnt.at(i_var)->get_phi();
//       this->m_MathModel.dpsi.at(i_var) = &fe_adjnt.at(i_var)->get_dphi();
//      }
//  
//    // Here we define some references to cell-specific data that
//    // will be used to assemble the linear system.
//    // The element Jacobian * quadrature weight at each integration point.   
//    const std::vector<Real>& JxW      = fe_adjnt.at(0)->get_JxW();
//    const std::vector<Real>& JxW_face = fe_face_adjnt.at(0)->get_JxW();
//  
//    // The physical XY locations of the quadrature points on the element.
//    // These are useful for evaluating spatially varying material
//    // properties at the quadrature points.
//    const std::vector<Point>& q_point     = fe_adjnt.at(0)->get_xyz();
//    const std::vector<Point>& qpoint_face = fe_face_adjnt.at(0)->get_xyz();
//
//    // loop over all the elements in the mesh, compute the elements
//    // contribution to the Hessian, then map the elements contribution
//    // to the global hessian
//    libMesh::MeshBase::const_element_iterator       el    =mesh.active_local_elements_begin();
//    const libMesh::MeshBase::const_element_iterator end_el=mesh.active_local_elements_end();
//    for ( ; el != end_el; ++el)
//     {
//
//      PetscLogEventBegin(AppSolve::logevents[10],0,0,0,0); // elemjac setup
//
//      // Store a pointer to the element we are currently
//      // working on.  This allows for nicer syntax later.
//      const Elem* elem = *el;
//
//      // set the field id in the spatially varying data structures
//      const unsigned int FieldId =  elem->_mat_data_field;
//
//      // Get the degree of freedom indices for the
//      // current element.  These define where in the global
//      // matrix and right-hand-side this element will
//      // contribute to.
//      dof_map_adjnt.dof_indices (elem, dof_indices_adjnt);
//      dof_map_state.dof_indices (elem, dof_indices_state);
//      const unsigned int n_dofs_adjnt = dof_indices_adjnt.size();
//      const unsigned int n_dofs_state = dof_indices_state.size();
//
//      // Ensure the element gradient is the right size and zero out the entries
//      HessVecProd.resize (qoiOptimizer->GetParamSize());
//
//      // state
//      for( unsigned int i_var = 0 ; i_var < n_vars_state ; i_var++)
//        { 
//         // degree of freedom indices for individual components
//         dof_map_state.dof_indices (elem,this->m_MathModel.dof_indices_u[i_var],i_var);
//         this->m_MathModel.n_u_dofs[i_var]  =   this->m_MathModel.dof_indices_u[i_var].size();
//
//         // Compute the element-specific data for the current
//         // element.  This involves computing the location of the
//         // quadrature points (q_point) and the shape functions
//         // (phi, dphi), (psi,dpsi),... for the current element.
//         fe_state[i_var]->reinit (elem);
//        } 
//
//      // adjoint
//      for( unsigned int i_var = 0 ; i_var < n_vars_adjnt ; i_var++)
//        { 
//         // degree of freedom indices for individual components
//         dof_map_adjnt.dof_indices (elem,this->m_MathModel.dof_indices_p[i_var],i_var);
//         this->m_MathModel.n_p_dofs[i_var] =    this->m_MathModel.dof_indices_p[i_var].size();
//
//         // Compute the element-specific data for the current
//         // element.  This involves computing the location of the
//         // quadrature points (q_point) and the shape functions
//         // (phi, dphi), (psi,dpsi),... for the current element.
//         fe_adjnt[i_var]->reinit (elem);
//        }
//
//      PetscLogEventEnd(  AppSolve::logevents[10],0,0,0,0); // elemjac setup
//
//      PetscLogEventBegin(AppSolve::logevents[11],0,0,0,0); // elemjac assble
//
//      /* accumulate hessian */
//      qoiOptimizer->getGlobalParamValues(elem,elemParameter);
//      qoiOptimizer->accumulateHessian(user,qrule,elemParameter,
//                                      FieldId,JxW,q_point,HessVecProd);
//
//      //accumulate sensitivity grad for verif
//      #if defined(PETSC_USE_DEBUG)
//      Gradiii += user->get_fem_dt() * // scale by deltat
//         qoiOptimizer->accumulateSensitivityGradient(user, qrule, 
//                                                     FieldId, JxW, q_point);
//      #endif
//      PetscLogEventEnd(  AppSolve::logevents[11],0,0,0,0); // elemjac assble
//
//      PetscLogEventBegin(AppSolve::logevents[12],0,0,0,0); // elemjac bndry
//
//      /* The following gets the error estimate contritubution on the 
//         boundary. Optimization parameters are assumed to have 
//         NO BC contribution */
//
//      PetscLogEventEnd(  AppSolve::logevents[12],0,0,0,0); // elemjac bndry
//
//      // get parameter mapping
//      qoiOptimizer->getGlobalMap(elem,globmap);
//
//      // scale the gradient by deltat
//      HessVecProd.scale(user->get_fem_dt());
//
//      // put the gradient in the tao data structures
//      info =  VecSetValues(yOut,qoiOptimizer->GetParamSize(),&globmap[0],
//      	              (PetscScalar*) &HessVecProd.get_values()[0] ,ADD_VALUES);
//      CHKERRQ(info);
//
//     } // end loop over element gradient computation
//
//    // free up memory
//    for( unsigned int i_var = 0 ; i_var < n_vars_state ; i_var++)
//      {
//       delete fe_state.at(i_var) ;
//       delete fe_face_state.at(i_var) ;
//      }
//    for( unsigned int i_var = 0 ; i_var < n_vars_adjnt ; i_var++)
//      {
//       delete fe_adjnt.at(i_var) ;
//       delete fe_face_adjnt.at(i_var) ;
//      }
//
//    CHKMEMQ; // check for memory corruption use -malloc_debug to enable
//  } // end time step loop
// //print sensitivity grad for verif
// #if defined(PETSC_USE_DEBUG)
// info =  VecSetValue(qoiOptimizer->hessianDiagonal,qoiOptimizer->hessianColumn,Gradiii,INSERT_VALUES);
// CHKERRQ(info);
// #endif
//
// //final assembly
// info = VecAssemblyBegin(yOut); CHKERRQ(info);
// info = VecAssemblyEnd(  yOut); CHKERRQ(info);
// // loop through all variance computations 1st to output to file
// // so don't have to wait for covariance plot
// //for( PetscInt ivariance = 0 ; ivariance < 2 ; ivariance++)
// // {
// //  if( ivariance && !libMesh::processor_id() ) 
// //                  std::cout << "Completed Hessian Row " << std::flush ; 
// //  for( PetscInt iii = 0 ; iii < qoiOptimizer->NDOF_CONTROL[1]; iii++)
// //   {
// //    /* first need to obtain the optimization paramter 
// //       pertaining to this global dof */
// //    const std::vector<optimizationParameter*>::iterator 
// //                  iiiParamIter = qoiOptimizer->getParameterPointer(iii) ;
//
// //    /* build upper triangular portion of hessian */
// //    //PetscInt jjj = iii ;
// //    for( PetscInt jjj = iii ; 
// //                  jjj < (ivariance?qoiOptimizer->NDOF_CONTROL[1]:iii+1) ;jjj++)
// //     {
// //      /* first need to obtain the optimization paramter 
// //         pertaining to this global dof */
// //      const std::vector<optimizationParameter*>::iterator 
// //                    jjjParamIter = qoiOptimizer->getParameterPointer(jjj) ;
//
// //      // for plotting sensitivities only instantiate during first pass through
// //      // sensitivities..
// //      // FIXME - there is a bug close the file otherwise
// //      ExodusII_IO *sensitivity_vis;
// //      if( !ivariance  ) sensitivity_vis = new ExodusII_IO::ExodusII_IO(user->mesh);
//
// //      // filename for plotting the sensitivitives
// //      OStringStream file_name;
// //      file_name << "femvis/femqoi_" << user->GroupID
// //                << "sensitivity_"   << jjj
// //                << "nopt_"          << user->IDOPT 
// //                << AppSolve::profileID  << ".e";
//
// //      if( ivariance ) qoiOptimizer->initializeSensitivity(sensitivity_system,jjj);
//
// //      for(AppSolve::ISTEP  = qoiOptimizer->Nsteplo() +1 ; 
// //          AppSolve::ISTEP <= qoiOptimizer->Nstephi()    ; AppSolve::ISTEP++)
// //       {
// //        if (ivariance) 
// //          {// read/recompute sensitivity into current_local_solution
// //           qoiOptimizer->getSensitivity( user , jjj,  jjjParamIter );
// //          }
// //        else
// //          {
// //          }
// //        else//variance entries 
// //          {
// //          }
//
// //        // loop over all the elements in the mesh, compute the elements
// //        // contribution to the Hessian, then map the elements contribution
// //        // to the global hessian
// //        for ( ; el != end_el; ++el)
// //         {
//
// //          /* Accumulate Hessian...  */
//
// //          // put the hessian in the tao data structures
// //          if ( iii != jjj ) 
// //           { // should be symmetric
// //            info =  MatSetValue(H,iii,jjj,Hessiiijjj,ADD_VALUES);
// //            info =  MatSetValue(H,jjj,iii,Hessiiijjj,ADD_VALUES);
// //            CHKERRQ(info);
// //           }
// //          else if ( !ivariance ) 
// //           { //accumulate variance
// //             info =  MatSetValue(H,iii,jjj,Hessiiijjj,ADD_VALUES);
// //           }
// //          else 
// //           {
// //             /* ivariance && iii == jjj ==> do nothing 
// //                just need to load data structures */
// //           }
// //         } // end loop over element hessian computation
//
// //       }//end time loop
// //      // clean up plotting (if necessary) 
// //      if( !ivariance  ) delete sensitivity_vis;
// //      // clean up this sensitivity (if necessary) 
// //      if( ivariance ) qoiOptimizer->finalizeSensitivity();
// //    }//end inner loop of Parameter loop
// //    if( ivariance && !libMesh::processor_id() )std::cout<<iii<<" "<<std::flush;
// //   }//end outer loop of Parameter loop
// //  // assemble after variance computation and plot
// //  info = MatAssemblyBegin(H,MAT_FINAL_ASSEMBLY); CHKERRQ(info);
// //  info = MatAssemblyEnd(  H,MAT_FINAL_ASSEMBLY); CHKERRQ(info);
// // }
// //if( !libMesh::processor_id() ) std::cout << std::endl << std::flush ; 
//
// CHKMEMQ; // check for memory corruption use -malloc_debug to enable
// PetscLogStagePop();
//
// // apply dirichlet domains if any
// std::string methodType = qoiOptimizer->GetTaoSolverMethod();
// if( this->m_dirichletNodes.size() 
//        && // skip if tao_fd_test
//     methodType.find("tao_fd_test")==std::string::npos
//   )
//   {
//    Vec dirichletVector;
//    info = VecDuplicate(yOut,&dirichletVector); CHKERRQ(info);
//    info = VecSet(dirichletVector,0.0); CHKERRQ(info);
//    qoiOptimizer->DirichletDomains( user, dirichletVector, 1.0 ); 
//    info = VecPointwiseMult( dirichletVector,dirichletVector, xIn) ; 
//    CHKERRQ(info);
//    info = VecAXPY(yOut,1.0,dirichletVector); CHKERRQ(info);
//    info = VecDestroy(dirichletVector); CHKERRQ(info);
//   }
 PetscFunctionReturn(0);
} // That's it.
