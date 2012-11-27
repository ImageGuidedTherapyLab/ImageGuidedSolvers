// Constructor
template < typename MathematicalModel >
BioHeatSystem< MathematicalModel > ::
BioHeatSystem(EquationSystems& es,
           const std::string& name,
           const unsigned int number ) 
//base class constructor
:PetscFEMSystem::PetscFEMSystem(es, name, number) ,
// initialize constitutive data
m_MathModel (*(es.parameters.get<GetPot*>("controlfile")), es) 
{
  GetPot &controlfile = *(es.parameters.get<GetPot*>("controlfile"));
  // Temperature and Damage
  // will be approximated using first-order approximation.
  // use ini file to possibly change variable name
  this->u_var = this->add_variable( "u0" , FIRST);
  this->a_var = this->add_variable( "d0" , FIRST);
  this->b_var = this->add_variable( "d1" , FIRST);
  // set for BC
  es.parameters.set<unsigned int>("u_var") = this->u_var;
  es.parameters.set<unsigned int>("a_var") = this->a_var;
  es.parameters.set<unsigned int>("b_var") = this->b_var;
  // set dirichlet BC
  int num_ids_u = controlfile.vector_variable_size("bc/u_dirichlet");
  for( int iii = 0; iii < num_ids_u ; iii++ )
    {
      int bc_id = controlfile("bc/u_dirichlet", -1, iii );
      this->m_DirichletNodeSets[this->u_var].push_back(bc_id);
    }
  // echo bc
  std::cout << " temperature dirichlet bc nodesets " ;
  std::vector<PetscInt>::iterator setIDIter ; 
  for(setIDIter  = this->m_DirichletNodeSets[this->u_var].begin(); 
      setIDIter != this->m_DirichletNodeSets[this->u_var].end();  setIDIter++)
   {
     std::cout << *setIDIter << " " ;
   }
  std::cout << std::endl << std::flush;
  // damage should default to all nodes a dirichlet data
  this->m_DirichletNodeSets[this->a_var].push_back(1);
  this->m_DirichletNodeSets[this->b_var].push_back(1);
  // make sure the same
  m_LinearSolve =  this->LinearPDE();

  // default initial condition is domain wise constant
  this->InitValues.push_back( &PetscFEMSystem::getInitialTemperature);

  // neuman bc data
  m_NeumannFlux.push_back( controlfile("bc/temp_flux",0.0) ) ; //[J/m^2/s]  
  // cauchy bc data
  m_newton_coeff.push_back( controlfile("bc/newton_coeff",1.0)) ; //[J/s/m^2/C] 
  m_u_infty.push_back(  controlfile("bc/u_infty",m_MathModel.GetBodyTemp())); // degC
  
}
// init the variables
template< typename MathematicalModel  >
void BioHeatSystem< MathematicalModel >::init_data ()
{
  PetscFunctionBegin;
  // setup the system variables
  this->m_MathModel.SetupState( this->get_equation_systems() );

  // Do the parent's initialization after variables are defined
  Parent::init_data();

  // Tell the system to march velocity forward in time, but 
  // leave p as a constraint only
  this->time_evolving(this->u_var);

  // set IC data
  this->InitValues.resize(3,NULL);
  this->InitValues[this->u_var]=&PetscFEMSystem::getInitialTemperature;
  this->InitValues[this->a_var]=&PetscFEMSystem::getInitialDamage;
  this->InitValues[this->b_var]=&PetscFEMSystem::getInitialDamage;

  m_LinearPDE = this->m_MathModel.CheckLinearity(this->get_equation_systems().parameters);

  PetscFunctionReturnVoid();
}

/* ------------------------------------------------------------ */
template< typename BioheatTransferModel  >
void BioHeatSystem< BioheatTransferModel >::
SetupDirichlet(libMesh::MeshBase& mesh)
{
  PetscFunctionBegin;

  // call base class function
  Parent::SetupDirichlet(mesh);

  // add damage and its derivative to the dirichlet nodesets
  for( unsigned int Ii = 0 ; Ii < this->m_subvector_indices.at(this->a_var).size(); ++Ii)
       this->m_NodeSets[this->a_var][1].push_back(this->m_subvector_indices[   this->a_var][Ii]);
  for( unsigned int Ii = 0 ; Ii < this->m_subvector_indices.at(this->b_var).size(); ++Ii)
       this->m_NodeSets[this->b_var][1].push_back(this->m_subvector_indices[   this->b_var][Ii]);

  PetscFunctionReturnVoid();
}

// init data structures
template< typename MathematicalModel  >
void BioHeatSystem< MathematicalModel >::
init_context(DiffContext &context)
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // We should prerequest all the data
  // we will need to build the linear system.
  c.element_fe_var[this->u_var]->get_JxW();
  c.element_fe_var[this->u_var]->get_phi();
  c.element_fe_var[this->u_var]->get_dphi();
  c.element_fe_var[this->u_var]->get_xyz();
  
  c.side_fe_var[this->u_var]->get_JxW();
  c.side_fe_var[this->u_var]->get_phi();
  c.side_fe_var[this->u_var]->get_xyz();
}


/**
 * @fn bool BioHeatSystem::element_time_derivative(bool request_jacobian, DiffContext &context)
 *
 * @section PennesGalkerinDiscretization Discretization of Equations
 * 
 * @htmlonly
 *   see pdf for in depth derivation of the 
 *   Element residual and jacobian calculations
 * @endhtmlonly
 *
 * @latexonly
 * The optimization problem in \eqn{mainoptprob} is solved using
 * an adjoint method to compute the gradient of the quantity of interest.
 * The following Galerkin representation of the temperature field and
 * adjoint variable is assumed.
 * 
 * \[
 * u(\mathbf{x},t) = \sum_{k=1}^{N_{step}} \sum_{j=1}^{N_{dof}}
 *                     \alpha_j^k(t) \phi_j(\textbf{x})
 * \qquad \qquad
 * p(\mathbf{x},t) = \sum_{k=1}^{N_{step}} \sum_{i=1}^{N_{dof}}
 *                        \lambda_i^k(t) \phi_i(\textbf{x})
 * \]
 * where $N_{step}$ is the number of time steps, $N_{dof}$ is the number
 * of Galerkin coefficients, and $\phi_i$'s are the finite element shape functions
 * of polynomial order p=1,2,3...
 * {\small
 * \[
 * \alpha_j^k(t)  = 
 * \left\{
 * \begin{split}
 *    \frac{t_k-t}{t_k-t_{k-1}} \alpha_j^{k-1} + 
 *    \frac{t-t_{k-1}}{t_k-t_{k-1}} \alpha_j^k & ,\; t \in [t_{k-1},t_k) \\
 *    0 \hspace{.8in} & , \; \text{ otherwise} \\
 * \end{split}
 * \right.
 * \quad 
 * \lambda_i^k(t)  = 
 * \left\{
 * \begin{split}
 *    \lambda_i^k & ,\; t \in [t_{k-1},t_k) \\
 *    0  & , \; \text{ otherwise} \\
 * \end{split}
 * \right.
 * \]
 * }
 * The time discretization of the power is assumed piecewise constant in time.
 * \[
 * P(t)  = 
 * \left\{
 * \begin{split}
 *    P_k & ,\qquad t \in [t_{k-1},t_k) \\
 *    0  & , \qquad \text{ otherwise} \\
 * \end{split}
 * \right.
 * \]
 * The spatial variation of the parameters fields is assumed
 * to have the following Galerkin representation
 * \[ \begin{split}
 * k_0(     \textbf{x})  & = \sum_{j} k_0^j \psi^j(\textbf{x}) \\
 * \omega_0(\textbf{x})  & = \sum_{j} \omega_0^j \psi^j(\textbf{x})
 * \end{split} \]
 * where $\psi(\textbf{x})$ are piecewise constant across elements.
 * 
 * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * \subsubsection{Time stepping.}
 * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * 
 * Assuming that the test function is piecewise constant in time
 * \[
 * v(\mathbf{x},t) = \sum_{k=1}^{N_{step}} \sum_{i=1}^{N_{dof}}
 *                        v_i^k(t) \phi_i(\textbf{x})
 * \qquad \qquad
 * v_i^k(t)  = 
 * \left\{
 * \begin{split}
 *    v_i^k & ,\qquad t \in [t_{k-1},t_k) \\
 *    0  & , \qquad \text{ otherwise} \\
 * \end{split}
 * \right.
 * \]
 * The governing equations \eqn{governpde} are solved with the
 * following Crank-Nicolson time stepping scheme.
 * \begin{equation} \label{diseqn}
 * \begin{split}
 *  \Delta t_k \; \int_{\Omega} & \rho  c_p \frac{u_k - u_{k-1}}{\Delta t_k} v_k
 *  +k(u_{\kmhalf},\textbf{x},\beta)  \nabla u_{\kmhalf}   \cdot \nabla  v_k 
 * \;  dx \\
 * \\
 *  & + \Delta t_k \; \int_{\Omega} \omega(u_{\kmhalf},\textbf{x},\beta)  c_{blood} ( u_{\kmhalf} - u_a ) \; v_k 
 * \;  dx \\
 * &+\Delta t_k \int_{\partial \Omega_C} h (u_{\kmhalf}-u_{\infty})\; v_k \; dA 
 *   =
 *  \Delta t_k \; \int_{\Omega}  Q_{laser}(\beta,\textbf{x},t_{\kmhalf}) v_k \; dx
 * \\ 
 * & \hspace{1.3in} - \int_{t_{k-1}}^{t_k} \int_{\partial \Omega_N} g  \; v_k  \; dA 
 * \quad \forall v_k  \quad k = 1,2, ..., N_{step} 
 * \end{split}
 * \end{equation}
 * where (using Einstein summation notation)
 * \[
 * u_k = \alpha_j^k \phi_j(\mathbf{x}) 
 * \qquad
 * u_{\kmhalf} = \frac{\alpha_j^{k-1}+\alpha_j^k}{2} \phi_j(\mathbf{x})
 * \qquad
 * v_k = v_i^k \phi_i(\mathbf{x}) 
 * \]
 * The discretization \eqn{diseqn} is of the form 
 * \[ \boxed{ \begin{split}
 *   \text{find } & \vec{\alpha}^k=(\alpha_1^k,\alpha_2^k,...) \text{ s.t. } \\
 *    & \vec{f}( \vec{\alpha}^k) = \vec{0}
 * \end{split} }  \]
 * The Jacobian of this system of equations is
 * \[
 * \begin{split}
 *   \frac{\partial f_i}{\partial \alpha_j^k} &= 
 *  \Delta t_k \; \int_{\Omega} \frac{\rho c_p}{\Delta t_k} \phi_j \phi_i \;dx
 *  +\Delta t_k \frac{1}{2} \int_{\partial \Omega_C} h \phi_j \; \phi_i \; dA \\
 * &+\Delta t_k \frac{1}{2} \int_{\Omega}  
 * \left(
 *   \frac{\partial k}{\partial u}(u_{\kmhalf},\textbf{x},\beta) \nabla u_{\kmhalf} \phi_j + 
 *   k(u_{\kmhalf},\textbf{x},\beta) \nabla \phi_j
 * \right) \cdot \nabla  \phi_i \; dx \\
 * &  +\Delta t_k \frac{1}{2} \int_{\Omega}  c_{blood}  
 * \left(
 * \frac{\partial \omega}{\partial u}(u_{\kmhalf},\textbf{x},\beta) 
 *              \left[ u_{\kmhalf} - u_a \right] +  \omega(u_{\kmhalf},\textbf{x},\beta) 
 * \right)   \phi_j \;  \phi_i \; dx \\
 * \end{split}
 * \]
 * @endlatexonly
 */
template< typename MathematicalModel  >
bool BioHeatSystem< MathematicalModel >::
element_time_derivative (bool request_jacobian, DiffContext &context)
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // set the field id in the spatially varying data structures

  // Initialize the per-element data for elem.
  std::vector<unsigned int> param_dof_indices;
  typename libMesh::System &template_parameter_system = 
                      this->get_equation_systems().get_system("k_0");
  template_parameter_system.get_dof_map().dof_indices (c.elem, param_dof_indices);

  // current_solution calls map_global_to_local_index that will map
  // this global index to the local current_local_solution
  const unsigned int field_id = param_dof_indices[0];
  const unsigned int subdomain_id = static_cast<int>( c.elem->subdomain_id() ) -1 ;

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> &JxW = 
    c.element_fe_var[this->u_var]->get_JxW();

  // The velocity shape functions at interior quadrature points.
  const std::vector<std::vector<Real> >& phi = 
    c.element_fe_var[this->u_var]->get_phi();

  // The velocity shape function gradients at interior
  // quadrature points.
  const std::vector<std::vector<RealGradient> >& dphi =
    c.element_fe_var[this->u_var]->get_dphi();

  // Physical location of the quadrature points
  const std::vector<Point>& qpoint = 
    c.element_fe_var[this->u_var]->get_xyz();
 
  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.dof_indices_var[this->u_var].size();

  // The subvectors and submatrices we need to fill:
  DenseSubMatrix<Number> &Kuu = *c.elem_subjacobians[this->u_var][this->u_var];
  DenseSubVector<Number>  &Fu = *c.elem_subresiduals[this->u_var];
      
  // Now we will build the element Jacobian and residual.
  // Constructing the residual requires the solution and its
  // gradient from the previous timestep.  This must be
  // calculated at each quadrature point by summing the
  // solution degree-of-freedom values by the appropriate
  // weight functions.
  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      // Compute the solution & its gradient at the old Newton iterate
      Number u_theta  = c.interior_value(   this->u_var,qp);
      Gradient grad_u = c.interior_gradient(this->u_var,qp);
      Number  z_value = 0.0; // not used

      // get damage values
      Number  damage  = c.interior_value(   this->a_var,qp);
      Number DdamageDu= c.interior_value(   this->b_var,qp);

      Gradient DiffusionDirection = this->m_MathModel.DiffusionDirection(subdomain_id) ; 
      Gradient TempDiffusionDirection( 
               grad_u(0)*DiffusionDirection(0)  ,
               grad_u(1)*DiffusionDirection(1)  ,
               grad_u(2)*DiffusionDirection(2)  
                                     ); 

      // First, an i-loop over the velocity degrees of freedom.
      // We know that n_u_dofs == n_v_dofs so we can compute contributions
      // for both at the same time.
      for (unsigned int i=0; i != n_u_dofs; i++)
        {
          Fu(i) += JxW[qp] * (
                phi[i][qp] * 
                 (              // perfusion term (check the SIGN)
                  this->m_MathModel.PennesReactionTerm(field_id,u_theta,damage)
                           -    // source term 
                  this->m_MathModel.PennesSource(field_id,u_theta,
                                                 damage,z_value, 
                                                 qpoint[qp],
                                                 this->m_PowerID)
                 )
                           +    // diffusion term
                this->m_MathModel.ThermalConductivity(field_id,u_theta,damage) *
                                         ( TempDiffusionDirection * dphi[i][qp] )
		) ;
          // convection term
          Fu(i) += JxW[qp] * phi[i][qp] * 
                ( this->m_MathModel.BulkFluidFlow(subdomain_id) * grad_u ) ; 
        }

      // Matrix contributions 
      if (request_jacobian)
        {
         if( this->m_StiffnessMatrix )
          for (unsigned int i=0; i != n_u_dofs; i++)
            for (unsigned int j=0; j != n_u_dofs; j++)
              {
                Kuu(i,j) += JxW[qp] * (
                       // perfusion term (check the SIGN)
                       (
                       this->m_MathModel.dPennesReactionTermdu(field_id,
                                           u_theta,damage,DdamageDu)
                       -
                       this->m_MathModel.dPennesSourcedu(field_id,
                                    u_theta,damage,DdamageDu,z_value,
                                    qpoint[qp],this->m_PowerID)
                       )
                                                   * phi[i][qp] * phi[j][qp]
                                  +    // diffusion term
                       this->m_MathModel.ThermalConductivity(field_id,
                                                             u_theta,damage) *
                                                  ( dphi[i][qp] * dphi[j][qp] )
                                  +    // diffusion derivative (non-symmetric)
                       this->m_MathModel.dThermalConductivitydu(field_id,
                                                     u_theta,damage,DdamageDu) *
                                          ( grad_u * dphi[i][qp] ) * phi[j][qp]
		) ;
                // convection derivative (non-symmetric)
                Kuu(i,j) += JxW[qp] * phi[i][qp] *
                       ( this->m_MathModel.BulkFluidFlow(subdomain_id) 
                         * dphi[j][qp] ) ;
              }
        }
    } // end of the quadrature point qp-loop
  
  return request_jacobian;
}
// boundary conditions
template< typename MathematicalModel  >
bool BioHeatSystem< MathematicalModel >::
side_time_derivative(bool request_jacobian, DiffContext &context)
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

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
      // loop over variables and compute diagonal terms only
      // derive classes call this for off diag terms 
      // Compute the solution at the theta timestep
      Number u_theta  = c.side_value(this->u_var, qp);

      // residual and jacobian for this variable
      DenseSubMatrix<Number> &Kuu = *c.elem_subjacobians[this->u_var][this->u_var];
      DenseSubVector<Number> &Fu = *c.elem_subresiduals[this->u_var];

      // Element Jacobian * quadrature weight for side integration
      const std::vector<Real> &JxW_side = c.side_fe_var[this->u_var]->get_JxW();
      
      // The velocity shape functions at side quadrature points.
      const std::vector<std::vector<Real> >& phi_side =
        c.side_fe_var[this->u_var]->get_phi();
      
      // Physical location of the quadrature points on the side
      const std::vector<Point>& qpoint = c.side_fe_var[this->u_var]->get_xyz();
      
      // The number of local degrees of freedom in u and v
      const unsigned int n_dofs = c.dof_indices_var[this->u_var].size(); 

      // bc_id = 3 ==> cauchy
      // bc_id = 2 ==> neumann
      short int bc_id =
        this->get_mesh().boundary_info->boundary_id(c.elem, c.side);
      libmesh_assert (bc_id != BoundaryInfo::invalid_id);

      for (unsigned int i=0; i != n_dofs; i++)
        {
          Fu(i) += JxW_side[qp] * phi_side[i][qp] *
                   // add bc to residual
                   CALL_MEMBER_FN(this, this->ResidualBC[bc_id][this->u_var])(this->u_var,u_theta);

          if (request_jacobian)
            for (unsigned int j=0; j != n_dofs; j++)
              {
                Kuu(i,j) += JxW_side[qp] * phi_side[i][qp] * phi_side[j][qp] *
                            // add bc to jacobian
                   CALL_MEMBER_FN(this, this->JacobianBC[bc_id][this->u_var])(this->u_var) ;
              }
        }
    } // end quadrature loop

  return request_jacobian;
}
/* --------------Pre Assembly for Dirichlet BC ---------- */
// called before nonlinear solve to setup diriclet BC
#undef __FUNCT__
#define __FUNCT__ "BioHeatSystem::ApplyDirichlet"
// assembly
template< typename MathematicalModel  >
void BioHeatSystem< MathematicalModel >::
ApplyDirichlet ()
{
  PetscFunctionBegin;
  PetscErrorCode ierr;

  // petsc Vecs
  Vec temperature, damageWork, Damage, DamageDeriv;
  // get subvector of temperature
  ierr = VecCreateMPI(libMesh::COMM_WORLD,
   		  PETSC_DECIDE,          // n_local
            this->m_subvector_indices.at(this->u_var).size(),   // n_global
   		  &temperature); CHKERRV(ierr);
  ierr = VecDuplicate( temperature,&damageWork ); CHKERRV(ierr);
  ierr = VecDuplicate( temperature,&Damage     ); CHKERRV(ierr);
  ierr = VecDuplicate( temperature,&DamageDeriv); CHKERRV(ierr);

  // PETSc data structures
  Vec Soln    = (dynamic_cast< PetscVector<double>* > (&(*this->solution)) )->vec();
  Vec OldSoln = (dynamic_cast< PetscVector<double>* > (&(this->get_vector("old_global_solution"))))->vec();

  // get temperature data
  ierr=VecScatterBegin(this->m_subvector_scatter[this->u_var],Soln,temperature,INSERT_VALUES,SCATTER_REVERSE);
  ierr=VecScatterEnd(  this->m_subvector_scatter[this->u_var],Soln,temperature,INSERT_VALUES,SCATTER_REVERSE);
  ierr = VecCopy(temperature,damageWork);
  
  // get subvectors  of old damage and derivative
  ierr=VecScatterBegin(this->m_subvector_scatter[this->a_var],OldSoln,Damage, INSERT_VALUES, SCATTER_REVERSE); 
  ierr=VecScatterEnd(  this->m_subvector_scatter[this->a_var],OldSoln,Damage, INSERT_VALUES, SCATTER_REVERSE); 
  ierr=VecScatterBegin(this->m_subvector_scatter[this->b_var],OldSoln,DamageDeriv,INSERT_VALUES,SCATTER_REVERSE); 
  ierr=VecScatterEnd(  this->m_subvector_scatter[this->b_var],OldSoln,DamageDeriv,INSERT_VALUES,SCATTER_REVERSE); 

  // compute arrhenius damage
  this->m_MathModel.ThermalDose(damageWork,this->deltat);
  ierr = VecAXPY (Damage,1.0,damageWork) ; 
  // compute arrhenius damage derivative wrt u
  this->m_MathModel.DoseDerivative(damageWork,temperature);
  ierr = VecAXPY (DamageDeriv,1.0,damageWork) ; 

  // scatter back to the global vector
  ierr=VecScatterBegin(this->m_subvector_scatter[this->a_var],Damage     ,Soln,INSERT_VALUES,SCATTER_FORWARD); 
  ierr=VecScatterEnd(  this->m_subvector_scatter[this->a_var],Damage     ,Soln,INSERT_VALUES,SCATTER_FORWARD); 
  ierr=VecScatterBegin(this->m_subvector_scatter[this->b_var],DamageDeriv,Soln,INSERT_VALUES,SCATTER_FORWARD); 
  ierr=VecScatterEnd(  this->m_subvector_scatter[this->b_var],DamageDeriv,Soln,INSERT_VALUES,SCATTER_FORWARD); 
  
  // localize for assembly
  this->solution->localize(*this->current_local_solution);

  // Clean up 
  ierr = VecDestroy(temperature );       
  ierr = VecDestroy(damageWork  );    
  ierr = VecDestroy(Damage      );    
  ierr = VecDestroy(DamageDeriv );    

  PetscFunctionReturnVoid();
}


//PetscErrorCode 
//BioHeatSystem::AssembleSystemDynamicsMatrix(Mat,int) 
//{
//  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     Build State Transition Matrix
//     Compute entries for the locally owned part of the Jacobian.
//      - Currently, all PETSc parallel matrix formats are partitioned by
//        contiguous chunks of rows across the processors. 
//      - Each processor needs to insert only elements that it owns
//        locally (but any non-local elements will be sent to the
//        appropriate processor during matrix assembly). 
//      - Here, we set all entries for a particular row at once.
//      - We can set matrix entries either using either
//        MatSetValuesLocal() or MatSetValues()
//   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//  PetscScalar    v[7],two = 2.0,hx,hy,hz;
//  MatStencil     col[7],row;
//
//  // spacing parameters
//  hx     = Setup::sp[0]; // 1.0/(PetscReal)(size[0]-1);
//  hy     = Setup::sp[1]; // 1.0/(PetscReal)(size[1]-1);
//  hz     = Setup::sp[2]; // 1.0/(PetscReal)(size[2]-1);
//
//  // fd parameters
//  PetscScalar sc      = hx*hz*hy;
//  PetscScalar hxhzdhy = hx*hz/hy;
//  PetscScalar hyhzdhx = hy*hz/hx;
//  PetscScalar hxhydhz = hx*hy/hz;
//  ierr=PetscPrintf(PETSC_COMM_WORLD,"sc=%12.5e hxhzdhy=%12.5e hyhzdhx=%12.5e hxhydhz=%12.5e \n",
//                                     sc,hxhzdhy,hyhzdhx,hxhydhz);CHKERRQ(ierr);
//  // build matrix entries
//  //  take FD equations and multiply by hx * hy * hz before assembling
//  if(ROIsize[2]==1){ // 2D
//     for (int j=ProcStart[1]; j<ProcStart[1]+ProcWidth[1]; j++) {
//       std::cout << " " << j << std::flush;
//       for (int i=ProcStart[0]; i<ProcStart[0]+ProcWidth[0]; i++) {
//         row.k = 0; row.j = j; row.i = i;
//         /* boundary points */
//         if (i==0 || j==0 || i==ROIsize[0]-1 || j==ROIsize[1]-1 ) {
//           v[0] = 1.0;
//           ierr = MatSetValuesStencil(StateA,1,&row,1,&row,v,INSERT_VALUES);CHKERRQ(ierr);
//           ierr = MatSetValuesStencil(StateB,1,&row,1,&row,v,INSERT_VALUES);CHKERRQ(ierr);
//         } else {
//         /* interior grid points */
//           v[0] = k_0 / rho / c_p * hxhzdhy; col[0].k=0;  col[0].j=j-1;col[0].i = i;
//           v[1] = k_0 / rho / c_p * hyhzdhx; col[1].k=0;  col[1].j=j;  col[1].i = i-1;
//           v[2] =  sc*( - w_0 * c_blood / rho / c_p ) 
//                      + 2.0 * k_0 / rho / c_p * (hyhzdhx+hxhzdhy);
//                            col[2].k=row.k;col[2].j=row.j;col[2].i = row.i;
//           v[3] = k_0 / rho / c_p * hyhzdhx; col[3].k=0;  col[3].j=j;  col[3].i = i+1;
//           v[4] = k_0 / rho / c_p * hxhzdhy; col[4].k=0;  col[4].j=j+1;col[4].i = i;
//           ierr = MatSetValuesStencil(m_SystemDynamics,1,&row,5,col,v,INSERT_VALUES);CHKERRQ(ierr);
//         }
//       } // end for (int i=ProcStart[0]; i<ProcStart[0]+ProcWidth[0]; i++) 
//     } // end for (int j=ProcStart[1]; j<ProcStart[1]+ProcWidth[1]; j++) 
//   }else{ // 3D
//     for (int k=ProcStart[2]; k<ProcStart[2]+ProcWidth[2]; k++) {
//       for (int j=ProcStart[1]; j<ProcStart[1]+ProcWidth[1]; j++) {
//         for (int i=ProcStart[0]; i<ProcStart[0]+ProcWidth[0]; i++) {
//           row.k = k; row.j = j; row.i = i;
//           /* boundary points */
//           if (i==0 || j==0 || i==ROIsize[0]-1 || j==ROIsize[1]-1 ) {
//             v[0] = 1.0;
//             ierr = MatSetValuesStencil(StateA,1,&row,1,&row,v,INSERT_VALUES);CHKERRQ(ierr);
//             ierr = MatSetValuesStencil(StateB,1,&row,1,&row,v,INSERT_VALUES);CHKERRQ(ierr);
//           } else if (k ==       0  ) {// out of plane zero flux bc
//             v[0] = -0.5 * k_0 / hz; col[0].k=k+1;col[0].j=j;  col[0].i = i;
//             v[1] =  0.5 * k_0 / hz; col[1].k=k  ;col[1].j=j;  col[1].i = i;
//             ierr = MatSetValuesStencil(StateA,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
//             ierr = MatSetValuesStencil(StateB,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
//           } else if (k == ROIsize[2]-1) {// out of plane zero flux bc
//             v[0] = -0.5 * k_0 / hz; col[0].k=k  ;col[0].j=j;  col[0].i = i;
//             v[1] =  0.5 * k_0 / hz; col[1].k=k-1;col[1].j=j;  col[1].i = i;
//             ierr = MatSetValuesStencil(StateA,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
//             ierr = MatSetValuesStencil(StateB,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
//           } else {
//           /* interior grid points */
//             v[0] =  k_0  / rho / c_p * hxhydhz; col[0].k=k-1;col[0].j=j;  col[0].i = i;
//             v[1] =  k_0  / rho / c_p * hxhzdhy; col[1].k=k;  col[1].j=j-1;col[1].i = i;
//             v[2] =  k_0  / rho / c_p * hyhzdhx; col[2].k=k;  col[2].j=j;  col[2].i = i-1;
//             v[3] =  sc*( - w_0 * c_blood / rho / c_p ) 
//                          + 2.0 * k_0 / rho / c_p * (hyhzdhx+hxhzdhy+hxhydhz);
//                              col[3].k=row.k;col[3].j=row.j;col[3].i = row.i;
//             v[4] =  k_0 / rho / c_p * hyhzdhx; col[4].k=k;  col[4].j=j;  col[4].i = i+1;
//             v[5] =  k_0 / rho / c_p * hxhzdhy; col[5].k=k;  col[5].j=j+1;col[5].i = i;
//             v[6] =  k_0 / rho / c_p * hxhydhz; col[6].k=k+1;col[6].j=j;  col[6].i = i;
//             ierr = MatSetValuesStencil(m_SystemDynamics,1,&row,7,col,v,INSERT_VALUES);CHKERRQ(ierr);
//           }
//         } // end for (int i=ProcStart[0]; i<ProcStart[0]+ProcWidth[0]; i++) 
//       } // end for (int j=ProcStart[1]; j<ProcStart[1]+ProcWidth[1]; j++) 
//     } // end for (int k=ProcStart[2]; k<ProcStart[2]+ProcWidth[2]; k++) 
//   }
//
//  /* Assemble matrix, using the 2-step process:
//       MatAssemblyBegin(), MatAssemblyEnd().  */
//  ierr = MatAssemblyBegin(m_SystemDynamics,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//  ierr = MatAssemblyEnd(  m_SystemDynamics,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//}
//// Assemble Load Vector for LITT Simulation
//PetscErrorCode 
//BioHeatSystem::AssembleSystemDynamicsMatrix(Mat,int) 
//{
//  // spacing parameters
//  PetscScalar    hx,hy,hz,sc;
//  hx     = Setup::sp[0]; // 1.0/(PetscReal)(size[0]-1);
//  hy     = Setup::sp[1]; // 1.0/(PetscReal)(size[1]-1);
//  hz     = Setup::sp[2]; // 1.0/(PetscReal)(size[2]-1);
//  sc     = hx*hz*hy;
//  //SOURCE TERM
//  PetscScalar dist=0.0;
//  std::vector< PetscScalar > xpoint(3,0.0);
//
//  /* Get pointers to vector data */
//  PetscScalar    ***source;
//  ierr = DAVecGetArray(dac,globalVec,&source);CHKERRQ(ierr);
//
//  /* Compute function over the locally owned part of the grid */
//  for (int k=ProcStart[2]; k<ProcStart[2]+ProcWidth[2]; k++) {
//    for (int j=ProcStart[1]; j<ProcStart[1]+ProcWidth[1]; j++) {
//      for (int i=ProcStart[0]; i<ProcStart[0]+ProcWidth[0]; i++) {
//        if (k ==  0 || k == ROIsize[2]-1 ) {// out of plane flux bc
//           source[k][j][i] = gflux; 
//        } else {
//
//        // laser model
//        xpoint[0] = ROIorgn[0] + (i+0.5) * Setup::sp[0] ;
//        xpoint[1] = ROIorgn[1] + (j+0.5) * Setup::sp[1] ;
//        xpoint[2] = ROIorgn[2] +    k    * Setup::sp[2] ;
//
//        source[k][j][i] = sc * m_Source.GetSource(xpoint,istep) ;
//
//        //if(source[k][j][i] != 0.0) 
//        //       std::cout<< source[k][j][i] << " i= " << i 
//        //                                   << " j= " << j 
//        //                                   << " k= " << k
//        //                << std::endl << std::flush;
//
//        // perfusion term
//        source[k][j][i] = source[k][j][i] + sc * w_0 * c_blood * u_a;
//
//        } 
//      }
//    }
//  }
//
//  /* Restore vectors */
//  ierr = DAVecRestoreArray(dac,globalVec,&source);CHKERRQ(ierr);
//
//}

// Constructor
template < typename MathematicalModel >
RFASystem< MathematicalModel > ::
RFASystem(EquationSystems& es,
           const std::string& name,
           const unsigned int number ) 
//base class constructor
:BioHeatSystem< MathematicalModel >::BioHeatSystem(es, 
                                             name,
                                             number) 
{
  // Temperature and Voltage
  // will be approximated using first-order approximation.
  this->z_var = this->add_variable ("u1", FIRST);
  // set for BC
  es.parameters.set<unsigned int>("z_var") = this->z_var;
  // get ref to control file
  GetPot &controlfile = *(es.parameters.get<GetPot*>("controlfile"));
  // set electrode domain
  this->m_ElectrodeNodeSet = controlfile("probe/electrodenodeset",5);
  std::cout << "Electrode Node Set " <<  this->m_ElectrodeNodeSet << std::endl;
  // set dirichlet BC
  int num_ids_z = controlfile.vector_variable_size("bc/z_dirichlet");
  for( int iii = 0; iii < num_ids_z ; iii++ )
    {
      int bc_id = controlfile("bc/z_dirichlet", -1, iii );
      this->m_DirichletNodeSets[this->z_var].push_back(bc_id);
    }
  // echo bc
  std::cout << " voltage dirichlet bc nodesets " ;
  std::vector<PetscInt>::iterator setIDIter ; 
  for(setIDIter  = this->m_DirichletNodeSets[this->z_var].begin(); 
      setIDIter != this->m_DirichletNodeSets[this->z_var].end();  setIDIter++)
   {
     std::cout << *setIDIter << " " ;
   }
  std::cout << std::endl << std::flush;

  // default initial condition is domain wise constant
  this->InitValues.push_back( &PetscFEMSystem::getInitialVoltage);

  // neuman bc data
  this->m_NeumannFlux.push_back( controlfile("bc/volt_flux",0.0) ) ;
  // cauchy bc data
  this->m_newton_coeff.push_back( controlfile("bc/v_newton_coeff",100.0));
  this->m_u_infty.push_back(  controlfile("bc/v_infty",0.0)); 

}
// init the variables
template< typename MathematicalModel  >
void RFASystem< MathematicalModel >::init_data ()
{

  // Do the parent's initialization after variables are defined
  Parent::init_data();

  // Tell the system to march temperature forward in time, but 
  // leave voltage as a constraint only
  this->time_evolving(this->u_var);

  // set IC data
  this->InitValues.resize(4,NULL);
  this->InitValues[this->u_var ]=&PetscFEMSystem::getInitialTemperature;
  this->InitValues[this->a_var]=&PetscFEMSystem::getInitialDamage;
  this->InitValues[this->b_var]=&PetscFEMSystem::getInitialDamage;
  this->InitValues[this->z_var ]=&PetscFEMSystem::getInitialVoltage;
}

/* --------------Pre Assembly for Dirichlet BC ---------- */
// called before nonlinear solve to setup diriclet BC
#undef __FUNCT__
#define __FUNCT__ "RFASystem::ApplyDirichlet"
// assembly
template< typename MathematicalModel  >
void RFASystem< MathematicalModel >::
ApplyDirichlet ()
{
 PetscFunctionBegin;

 // update the damage
 Parent::ApplyDirichlet();

 // apply the current power as the dirichlet data if any
 PetscScalar CurrentPower = this->m_MathModel.getPower(this->m_PowerID); 
 std::vector<PetscInt> &dirichletNodeSet= 
               this->m_NodeSets[this->z_var][this->m_ElectrodeNodeSet];
 //  only apply dirichlet data for typical assembly
 if( dirichletNodeSet.size() && 
     this->m_MassMatrix && this->m_StiffnessMatrix )
   {
     this->solution->close();
     for( unsigned int Ii = 0; Ii<dirichletNodeSet.size();Ii++) 
                       this->solution->set(dirichletNodeSet[Ii],CurrentPower);
     this->solution->close();
   }

 // localize for assembly
 this->solution->localize(*this->current_local_solution);

 PetscFunctionReturnVoid();

}

template< typename MathematicalModel  >
bool RFASystem< MathematicalModel >::
element_time_derivative (bool request_jacobian, DiffContext &context)
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // Initialize the per-element data for elem.
  std::vector<unsigned int> param_dof_indices;
  typename libMesh::System &template_parameter_system = 
                      this->get_equation_systems().get_system("s_0");
  template_parameter_system.get_dof_map().dof_indices (c.elem, param_dof_indices);

  // current_solution calls map_global_to_local_index that will map
  // this global index to the local current_local_solution
  const unsigned int field_id = param_dof_indices[0];
  const unsigned int subdomain_id = static_cast<int>( c.elem->subdomain_id() ) -1 ;

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> &JxW = 
    c.element_fe_var[this->u_var]->get_JxW();

  // shape function and shape function gradients at interior
  // quadrature points. generally phi = psi but used different shape functions
  // for illustrative purposes
  const std::vector<std::vector<Real> >& phi =
    c.element_fe_var[this->u_var]->get_phi();
  const std::vector<std::vector<Real> >& psi =
    c.element_fe_var[this->z_var]->get_phi();
  //const std::vector<std::vector<Real> >& psi =
  //  c.element_fe_var[this->z_var]->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi =
    c.element_fe_var[this->u_var]->get_dphi();
  const std::vector<std::vector<RealGradient> >& dpsi =
    c.element_fe_var[this->z_var]->get_dphi();

  // Physical location of the quadrature points
  const std::vector<Point>& qpoint = 
    c.element_fe_var[this->u_var]->get_xyz();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.dof_indices_var[this->u_var].size(); 
  const unsigned int n_z_dofs = c.dof_indices_var[this->z_var].size(); 
 
  // The subvectors and submatrices we need to fill:
  //const unsigned int dim = this->get_mesh().mesh_dimension();
  DenseSubMatrix<Number> &Kuu = *c.elem_subjacobians[this->u_var][this->u_var];
  DenseSubMatrix<Number> &Kzz = *c.elem_subjacobians[this->z_var][this->z_var];
  DenseSubMatrix<Number> &Kuz = *c.elem_subjacobians[this->u_var][this->z_var];
  DenseSubMatrix<Number> &Kzu = *c.elem_subjacobians[this->z_var][this->u_var];
  DenseSubVector<Number> &Fu = *c.elem_subresiduals[this->u_var];
  DenseSubVector<Number> &Fz = *c.elem_subresiduals[this->z_var];
      
  // Now we will build the element Jacobian and residual.
  // Constructing the residual requires the solution and its
  // gradient from the previous timestep.  This must be
  // calculated at each quadrature point by summing the
  // solution degree-of-freedom values by the appropriate
  // weight functions.
  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      // Compute the solution & its gradient at the old Newton iterate
      Number u_theta  = c.interior_value(   this->u_var,qp);
      Gradient grad_u = c.interior_gradient(this->u_var,qp);
      Number z_value  = c.interior_value(   this->z_var,qp);
      Gradient grad_z = c.interior_gradient(this->z_var,qp);
      // get damage values
      Number  damage  = c.interior_value(   this->a_var,qp);
      Number DdamageDu= c.interior_value(   this->b_var,qp);

      // First, an i-loop over the velocity degrees of freedom.
      // We know that n_u_dofs == n_v_dofs so we can compute contributions
      // for both at the same time.
      for (unsigned int i=0; i != n_u_dofs; i++)
        {
          Fu(i) += JxW[qp] * (
                phi[i][qp] * 
                 (              // perfusion term (check the SIGN)
                  this->m_MathModel.PennesReactionTerm(field_id,
                                                       u_theta,damage)
                           -    // source term 
                  this->m_MathModel.PennesSource(field_id,u_theta,
                                                 damage,(grad_z*grad_z),
                                                 qpoint[qp],
                                                 this->m_PowerID)
                 )
                           +    // diffusion term
                this->m_MathModel.ThermalConductivity(field_id,
                                                       u_theta,damage) *
                                                       (grad_u * dphi[i][qp])
		) ;
          // convection term
          Fu(i) += JxW[qp] * phi[i][qp] * 
                ( this->m_MathModel.BulkFluidFlow(subdomain_id) * grad_u ) ; 

          // Matrix contributions for the uu and vv couplings.
          if (request_jacobian)
            for (unsigned int j=0; j != n_u_dofs; j++)
              {
                Kuu(i,j) += JxW[qp] * (
                         ( // perfusion term (check the SIGN)
                       this->m_MathModel.dPennesReactionTermdu(field_id,
                                              u_theta,damage,DdamageDu)
                              -    // RF source term 
                       this->m_MathModel.dPennesSourcedu(field_id,
                                           u_theta,damage,DdamageDu,
                                           (grad_z*grad_z),
                                           qpoint[qp],this->m_PowerID)
                         ) * ( phi[i][qp] * phi[j][qp] )
                                  +    // diffusion term
                      this->m_MathModel.ThermalConductivity(field_id,
                                              u_theta,damage) *
                                                  ( dphi[i][qp] * dphi[j][qp] )
                                  +    // diffusion derivative (non-symmetric)
                      this->m_MathModel.dThermalConductivitydu(field_id,
                                                     u_theta,damage,DdamageDu) *
                                    ( grad_u * dphi[i][qp] ) * phi[j][qp]
		) ;
                // convection derivative (non-symmetric)
                Kuu(i,j) += JxW[qp] * phi[i][qp] *
                       ( this->m_MathModel.BulkFluidFlow(subdomain_id) *
dphi[j][qp] ) ;

                // coupling terms
                Kuz(i,j) += JxW[qp] * (
                  2.0 * this->m_MathModel.ElectricConductivity(field_id,u_theta,damage) *
                                   phi[i][qp] * ( dpsi[j][qp] * grad_z ) 
		) ;
              }

        }
    } // end of the quadrature point qp-loop
  
  return request_jacobian;
}
template< typename MathematicalModel  >
bool RFASystem< MathematicalModel >::
element_constraint (bool request_jacobian, DiffContext &context)
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // Initialize the per-element data for elem.
  std::vector<unsigned int> param_dof_indices;
  typename libMesh::System &template_parameter_system = 
                      this->get_equation_systems().get_system("s_0");
  template_parameter_system.get_dof_map().dof_indices (c.elem, param_dof_indices);

  // current_solution calls map_global_to_local_index that will map
  // this global index to the local current_local_solution
  const unsigned int field_id = param_dof_indices[0];
  const unsigned int subdomain_id = static_cast<int>( c.elem->subdomain_id() ) -1 ;

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> &JxW = 
    c.element_fe_var[this->u_var]->get_JxW();

  // shape function and shape function gradients at interior
  // quadrature points. generally phi = psi but used different shape functions
  // for illustrative purposes
  const std::vector<std::vector<Real> >& phi =
    c.element_fe_var[this->u_var]->get_phi();
  const std::vector<std::vector<Real> >& psi =
    c.element_fe_var[this->z_var]->get_phi();
  //const std::vector<std::vector<Real> >& psi =
  //  c.element_fe_var[this->z_var]->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi =
    c.element_fe_var[this->u_var]->get_dphi();
  const std::vector<std::vector<RealGradient> >& dpsi =
    c.element_fe_var[this->z_var]->get_dphi();

  // Physical location of the quadrature points
  const std::vector<Point>& qpoint = 
    c.element_fe_var[this->u_var]->get_xyz();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.dof_indices_var[this->u_var].size(); 
  const unsigned int n_z_dofs = c.dof_indices_var[this->z_var].size(); 
 
  // The subvectors and submatrices we need to fill:
  //const unsigned int dim = this->get_mesh().mesh_dimension();
  DenseSubMatrix<Number> &Kuu = *c.elem_subjacobians[this->u_var][this->u_var];
  DenseSubMatrix<Number> &Kzz = *c.elem_subjacobians[this->z_var][this->z_var];
  DenseSubMatrix<Number> &Kuz = *c.elem_subjacobians[this->u_var][this->z_var];
  DenseSubMatrix<Number> &Kzu = *c.elem_subjacobians[this->z_var][this->u_var];
  DenseSubVector<Number> &Fu = *c.elem_subresiduals[this->u_var];
  DenseSubVector<Number> &Fz = *c.elem_subresiduals[this->z_var];
      
  // Now we will build the element Jacobian and residual.
  // Constructing the residual requires the solution and its
  // gradient from the previous timestep.  This must be
  // calculated at each quadrature point by summing the
  // solution degree-of-freedom values by the appropriate
  // weight functions.
  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      // Compute the solution & its gradient at the old Newton iterate
      Number u_theta  = c.interior_value(   this->u_var,qp);
      Gradient grad_u = c.interior_gradient(this->u_var,qp);
      Number z_value  = c.interior_value(   this->z_var,qp);
      Gradient grad_z = c.interior_gradient(this->z_var,qp);
      // get damage values
      Number  damage  = c.interior_value(   this->a_var,qp);
      Number DdamageDu= c.interior_value(   this->b_var,qp);

      // Now a loop over the voltage degrees of freedom.  This
      // computes the contributions of the continuity equation.
      for (unsigned int i=0; i != n_z_dofs; i++)
        {
          Fz(i) += JxW[qp] * ( grad_z * dpsi[i][qp] ) *
                   this->m_MathModel.ElectricConductivity(field_id,u_theta,damage) ;

          if (request_jacobian)
            for (unsigned int j=0; j != n_u_dofs; j++)
              {
                Kzu(i,j) += JxW[qp] * ( grad_z * dpsi[i][qp]) * phi[j][qp]   *
                      this->m_MathModel.dElectricConductivitydu(field_id,u_theta,damage,DdamageDu) ;
                Kzz(i,j) += JxW[qp] * ( dpsi[i][qp] * dpsi[j][qp] ) *
                      this->m_MathModel.ElectricConductivity(   field_id,u_theta,damage) ;
              }
        }
    } // end of the quadrature point qp-loop
  
  return request_jacobian;
}
// boundary conditions  pick up off diagonal terms
template< typename MathematicalModel  >
bool RFASystem< MathematicalModel >::
side_constraint(bool request_jacobian,DiffContext &context)
{
  // call parent by default
  Parent::side_time_derivative(request_jacobian,context);

  // continue with off diagonal updates
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

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

  // The subvectors and submatrices we need to fill:
  //const unsigned int dim = this->get_mesh().mesh_dimension();
  DenseSubVector<Number> &Fz  = *c.elem_subresiduals[this->z_var];
  DenseSubMatrix<Number> &Kzz = *c.elem_subjacobians[this->z_var][this->z_var];
  DenseSubMatrix<Number> &Kuz = *c.elem_subjacobians[this->u_var][this->z_var];

  // Element Jacobian * quadrature weight for side integration
  const std::vector<Real> &JxW_side = c.side_fe_var[this->u_var]->get_JxW();
  
  // The shape functions at side quadrature points.
  const std::vector<std::vector<Real> >& phi_side = c.side_fe_var[this->u_var]->get_phi();
  const std::vector<std::vector<Real> >& psi_side = c.side_fe_var[this->z_var]->get_phi();
  
  // Physical location of the quadrature points on the side
  const std::vector<Point>& qpoint = c.side_fe_var[this->u_var]->get_xyz();
  
  // The number of local degrees of freedom in u and v
  const unsigned int n_u_dofs = c.dof_indices_var[this->u_var].size(); 
  const unsigned int n_z_dofs = c.dof_indices_var[this->z_var].size(); 

  // bc_id = 3 ==> cauchy
  // bc_id = 2 ==> neumann
  short int bc_id =
    this->get_mesh().boundary_info->boundary_id(c.elem, c.side);
  libmesh_assert (bc_id != BoundaryInfo::invalid_id);

    // residual
    for (unsigned int qp=0; qp != n_sidepoints; qp++)
      {
      Number u_theta  = c.side_value(   this->u_var,qp);
      Number z_value  = c.side_value(         this->z_var,qp);
        // loop over the fluence degrees of freedom.  
        if (request_jacobian)
        for (unsigned int i=0; i != n_z_dofs; i++)
          {
            Fz(i) += JxW_side[qp] * phi_side[i][qp] *
                     // add bc to residual
                     CALL_MEMBER_FN(this, this->ResidualBC[bc_id][this->z_var])(this->z_var,u_theta);
            if (request_jacobian)
              for (unsigned int j=0; j != n_z_dofs; j++)
                {
                  Kzz(i,j) += JxW_side[qp] * phi_side[i][qp] * phi_side[j][qp] *
                              // add bc to jacobian
                     CALL_MEMBER_FN(this, this->JacobianBC[bc_id][this->z_var])(this->z_var) ;
                }
          }
      } // end quadrature loop

  return request_jacobian;
}

// Constructor
template < typename MathematicalModel >
RHTESystem< MathematicalModel > ::
RHTESystem(EquationSystems& es,
           const std::string& name,
           const unsigned int number ) 
//base class constructor
:BioHeatSystem< MathematicalModel >::BioHeatSystem(es, 
                                             name,
                                             number) 
{
  // Temperature and fluence
  // will be approximated using first-order approximation.
   this->z_var = this->add_variable ("u1", FIRST);
   this->e_var = this->add_variable ("u2", FIRST);
  this->fx_var = this->add_variable ("fx", FIRST);
  this->fy_var = this->add_variable ("fy", FIRST);
  this->fz_var = this->add_variable ("fz", FIRST);
  // set for BC
  es.parameters.set<unsigned int>("z_var") = this->z_var;

  // determine whether should solve as interstitial or external
  GetPot &controlfile = *(es.parameters.get<GetPot*>("controlfile"));
  if (controlfile("probe/interstitial",false)) 
    {
     PetscPrintf(PETSC_COMM_WORLD,"RHTESystem setting up interstitial fiber...\n");
     this->m_ExternalFiber = PETSC_FALSE ; 
    }
  else
    {
     PetscPrintf(PETSC_COMM_WORLD,"RHTESystem setting up External fiber...\n");
     this->m_ExternalFiber = PETSC_TRUE ; 
    }
  int num_ids_z = controlfile.vector_variable_size("bc/z_dirichlet");
  for( int iii = 0; iii < num_ids_z ; iii++ )
    {
      int bc_id = controlfile("bc/z_dirichlet", -1, iii );
      this->m_DirichletNodeSets[this->z_var].push_back(bc_id);
    }
  // echo bc
  std::cout << " fluence dirichlet bc nodesets " ;
  std::vector<PetscInt>::iterator setIDIter ; 
  for(setIDIter  = this->m_DirichletNodeSets[this->z_var].begin(); 
      setIDIter != this->m_DirichletNodeSets[this->z_var].end();  setIDIter++)
   {
     std::cout << *setIDIter << " " ;
   }
  std::cout << std::endl << std::flush;

  // neuman bc data
  this->m_NeumannFlux.push_back( controlfile("bc/fluence_flux",0.0) ) ; 
  // cauchy bc data
  this->m_newton_coeff.push_back( controlfile("bc/fluence_newton_coeff",0.0)); 
  this->m_u_infty.push_back(  controlfile("bc/fluence_infty",0.0)); 

}

// init the variables
template< typename MathematicalModel  >
void RHTESystem < MathematicalModel >::init_data ()
{

  // Do the parent's initialization after variables are defined
  Parent::init_data();

  // Tell the system to march temperature forward in time, but 
  // leave voltage as a constraint only
  this->time_evolving(this->u_var);

  // set IC data
  this->InitValues.resize(8,NULL);
  this->InitValues[ this->u_var]=&PetscFEMSystem::getInitialTemperature;
  this->InitValues[ this->a_var]=&PetscFEMSystem::getInitialDamage;
  this->InitValues[ this->b_var]=&PetscFEMSystem::getInitialDamage;
  this->InitValues[ this->z_var]=&PetscFEMSystem::getInitialFluence;

  if(this->m_ExternalFiber)
    {// set pointers to external fiber routines
     this->InitValues[ this->e_var]=&PetscFEMSystem::getExternalIrradiance;
     this->InitValues[this->fx_var]=&PetscFEMSystem::getExternalFlux_X;
     this->InitValues[this->fy_var]=&PetscFEMSystem::getExternalFlux_Y;
     this->InitValues[this->fz_var]=&PetscFEMSystem::getExternalFlux_Z;
    }
  else
    {// set pointers to interstitial applicator
     this->InitValues[ this->e_var]=&PetscFEMSystem::getInterstitialIrradiance;
     this->InitValues[this->fx_var]=&PetscFEMSystem::getInterstitialFlux_X;
     this->InitValues[this->fy_var]=&PetscFEMSystem::getInterstitialFlux_Y;
     this->InitValues[this->fz_var]=&PetscFEMSystem::getInterstitialFlux_Z;
    }
}

/* ------------------------------------------------------------ */
template< typename BioheatTransferModel  >
void RHTESystem< BioheatTransferModel >::
SetupDirichlet(libMesh::MeshBase& mesh)
{
  PetscFunctionBegin;

  // call base class function
  Parent::SetupDirichlet(mesh);

  // scattered fluence should have not dirichlet data
  this->m_NodeSets[this->z_var][1].clear();
  // add damage and its derivative to the dirichlet nodesets
  for( unsigned int Ii = 0 ; Ii < this->m_subvector_indices.at(this->e_var).size(); ++Ii)
       this->m_NodeSets[this->e_var][1].push_back(this->m_subvector_indices[   this->e_var][Ii]);
  for( unsigned int Ii = 0 ; Ii < this->m_subvector_indices.at(this->fx_var).size(); ++Ii)
       this->m_NodeSets[this->fx_var][1].push_back(this->m_subvector_indices[   this->fx_var][Ii]);
  for( unsigned int Ii = 0 ; Ii < this->m_subvector_indices.at(this->fy_var).size(); ++Ii)
       this->m_NodeSets[this->fy_var][1].push_back(this->m_subvector_indices[   this->fy_var][Ii]);
  for( unsigned int Ii = 0 ; Ii < this->m_subvector_indices.at(this->fz_var).size(); ++Ii)
       this->m_NodeSets[this->fz_var][1].push_back(this->m_subvector_indices[   this->fz_var][Ii]);

  PetscFunctionReturnVoid();
}
// element assembly
template< typename MathematicalModel  >
bool RHTESystem< MathematicalModel >::
element_time_derivative (bool request_jacobian, DiffContext &context)
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // Initialize the per-element data for elem.
  std::vector<unsigned int> param_dof_indices;
  typename libMesh::System &template_parameter_system = 
                      this->get_equation_systems().get_system("k_0");
  template_parameter_system.get_dof_map().dof_indices (c.elem, param_dof_indices);

  // current_solution calls map_global_to_local_index that will map
  // this global index to the local current_local_solution
  const unsigned int field_id = param_dof_indices[0];
  const unsigned int subdomain_id = static_cast<int>( c.elem->subdomain_id() ) -1 ;

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> &JxW = 
    c.element_fe_var[this->u_var]->get_JxW();

  // shape function and shape function gradients at interior
  // quadrature points. generally phi = psi but used different shape functions
  // for illustrative purposes
  const std::vector<std::vector<Real> >& phi =
    c.element_fe_var[this->u_var]->get_phi();
  const std::vector<std::vector<Real> >& psi =
    c.element_fe_var[this->z_var]->get_phi();
  //const std::vector<std::vector<Real> >& psi =
  //  c.element_fe_var[this->z_var]->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi =
    c.element_fe_var[this->u_var]->get_dphi();
  const std::vector<std::vector<RealGradient> >& dpsi =
    c.element_fe_var[this->z_var]->get_dphi();

  // Physical location of the quadrature points
  const std::vector<Point>& qpoint = 
    c.element_fe_var[this->u_var]->get_xyz();
 
  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.dof_indices_var[this->u_var].size(); 
  const unsigned int n_z_dofs = c.dof_indices_var[this->z_var].size(); 

  // The subvectors and submatrices we need to fill:
  //const unsigned int dim = this->get_mesh().mesh_dimension();
  DenseSubMatrix<Number> &Kuu = *c.elem_subjacobians[this->u_var][this->u_var];
  DenseSubMatrix<Number> &Kzz = *c.elem_subjacobians[this->z_var][this->z_var];
  DenseSubMatrix<Number> &Kuz = *c.elem_subjacobians[this->u_var][this->z_var];
  //DenseSubMatrix<Number> &Kzu = *c.elem_subjacobians[this->z_var][this->u_var];
  DenseSubVector<Number> &Fu = *c.elem_subresiduals[this->u_var];
  DenseSubVector<Number> &Fz = *c.elem_subresiduals[this->z_var];
      
  // Now we will build the element Jacobian and residual.
  // Constructing the residual requires the solution and its
  // gradient from the previous timestep.  This must be
  // calculated at each quadrature point by summing the
  // solution degree-of-freedom values by the appropriate
  // weight functions.
  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      // Compute the solution & its gradient 
      Number u_theta  = c.interior_value(   this->u_var,qp);
      Gradient grad_u = c.interior_gradient(this->u_var,qp);
      Number z_value  = c.interior_value(         this->z_var,qp);
      Gradient grad_z = c.interior_gradient(      this->z_var,qp);
      // get irradiance values
      Number  e_value = c.interior_value(   this->e_var,qp); // irradiance value
      Number  e_s0_x  = c.interior_value(  this->fx_var,qp); // irradiance flux in x direction
      Number  e_s0_y  = c.interior_value(  this->fy_var,qp); // irradiance flux in y direction
      Number  e_s0_z  = c.interior_value(  this->fz_var,qp); // irradiance flux in z direction
      Gradient e_s0(e_s0_x,e_s0_y,e_s0_z); //irradiance flux 
      // get damage values
      Number  damage  = c.interior_value(   this->a_var,qp);
      Number DdamageDu= c.interior_value(   this->b_var,qp);

      // First, an i-loop over the velocity degrees of freedom.
      // We know that n_u_dofs == n_v_dofs so we can compute contributions
      // for both at the same time.
      for (unsigned int i=0; i != n_u_dofs; i++)
        {
          Fu(i) += JxW[qp] * (
                phi[i][qp] * 
                 (              // perfusion term (check the SIGN)
                  this->m_MathModel.PennesReactionTerm(field_id,
                                                       u_theta,damage)
                           -  
                  // source term is sum of scattered and collimated light
                  this->m_MathModel.PennesSource(field_id,u_theta,
                                                 damage,e_value , 
                                                 qpoint[qp],
                                                 this->m_PowerID)
                           -  
                  this->m_MathModel.ScatteredSource(field_id,u_theta,damage,z_value) 
                 )
                           +    // diffusion term
                this->m_MathModel.ThermalConductivity(field_id,u_theta,damage) *
                                                       (grad_u * dphi[i][qp])
		) ;

          // Matrix contributions for the uu and vv couplings.
          if (request_jacobian)
            for (unsigned int j=0; j != n_u_dofs; j++)
              {
                Kuu(i,j) += JxW[qp] * (
                         ( // perfusion term (check the SIGN)
                      this->m_MathModel.dPennesReactionTermdu(field_id,
                                               u_theta,damage,DdamageDu) 
                      -
                      this->m_MathModel.dPennesSourcedu(field_id,
                                    u_theta,damage,DdamageDu,z_value,
                                    qpoint[qp],this->m_PowerID)
                         ) * ( phi[i][qp] * phi[j][qp] )
                                  +    // diffusion term
                      this->m_MathModel.ThermalConductivity(field_id,u_theta,damage) *
                                                  ( dphi[i][qp] * dphi[j][qp] )
                                  +    // diffusion derivative (non-symmetric)
                      this->m_MathModel.dThermalConductivitydu(field_id,
                                                     u_theta,damage,DdamageDu) *
                                          ( grad_u * dphi[i][qp] ) * phi[j][qp]
		) ;
                // coupling terms
                Kuz(i,j) += JxW[qp] * (
                           -   phi[i][qp] *  // source term 
                  this->m_MathModel.ScatteredSource(field_id,u_theta,damage,psi[j][qp]) 

		) ;
              }

        }
    } // end of the quadrature point qp-loop
  
  return request_jacobian;
}

// element assembly
template< typename MathematicalModel  >
bool RHTESystem< MathematicalModel >::
element_constraint (bool request_jacobian, DiffContext &context)
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // Initialize the per-element data for elem.
  std::vector<unsigned int> param_dof_indices;
  typename libMesh::System &template_parameter_system = 
                      this->get_equation_systems().get_system("k_0");
  template_parameter_system.get_dof_map().dof_indices (c.elem, param_dof_indices);

  // current_solution calls map_global_to_local_index that will map
  // this global index to the local current_local_solution
  const unsigned int field_id = param_dof_indices[0];
  const unsigned int subdomain_id = static_cast<int>( c.elem->subdomain_id() ) -1 ;

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> &JxW = 
    c.element_fe_var[this->u_var]->get_JxW();

  // shape function and shape function gradients at interior
  // quadrature points. generally phi = psi but used different shape functions
  // for illustrative purposes
  const std::vector<std::vector<Real> >& phi =
    c.element_fe_var[this->u_var]->get_phi();
  const std::vector<std::vector<Real> >& psi =
    c.element_fe_var[this->z_var]->get_phi();
  //const std::vector<std::vector<Real> >& psi =
  //  c.element_fe_var[this->z_var]->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi =
    c.element_fe_var[this->u_var]->get_dphi();
  const std::vector<std::vector<RealGradient> >& dpsi =
    c.element_fe_var[this->z_var]->get_dphi();

  // Physical location of the quadrature points
  const std::vector<Point>& qpoint = 
    c.element_fe_var[this->u_var]->get_xyz();
 
  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.dof_indices_var[this->u_var].size(); 
  const unsigned int n_z_dofs = c.dof_indices_var[this->z_var].size(); 

  // The subvectors and submatrices we need to fill:
  //const unsigned int dim = this->get_mesh().mesh_dimension();
  DenseSubMatrix<Number> &Kuu = *c.elem_subjacobians[this->u_var][this->u_var];
  DenseSubMatrix<Number> &Kzz = *c.elem_subjacobians[this->z_var][this->z_var];
  DenseSubMatrix<Number> &Kuz = *c.elem_subjacobians[this->u_var][this->z_var];
  //DenseSubMatrix<Number> &Kzu = *c.elem_subjacobians[this->z_var][this->u_var];
  DenseSubVector<Number> &Fu = *c.elem_subresiduals[this->u_var];
  DenseSubVector<Number> &Fz = *c.elem_subresiduals[this->z_var];
      
  // Now we will build the element Jacobian and residual.
  // Constructing the residual requires the solution and its
  // gradient from the previous timestep.  This must be
  // calculated at each quadrature point by summing the
  // solution degree-of-freedom values by the appropriate
  // weight functions.
  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      // Compute the solution & its gradient 
      Number u_theta  = c.interior_value(   this->u_var,qp);
      Gradient grad_u = c.interior_gradient(this->u_var,qp);
      Number z_value  = c.interior_value(         this->z_var,qp);
      Gradient grad_z = c.interior_gradient(      this->z_var,qp);
      // get irradiance values
      Number  e_value = c.interior_value(   this->e_var,qp); // irradiance value
      Number  e_s0_x  = c.interior_value(  this->fx_var,qp); // irradiance flux in x direction
      Number  e_s0_y  = c.interior_value(  this->fy_var,qp); // irradiance flux in y direction
      Number  e_s0_z  = c.interior_value(  this->fz_var,qp); // irradiance flux in z direction
      Gradient e_s0(e_s0_x,e_s0_y,e_s0_z); //irradiance flux 
      // get damage values
      Number  damage  = c.interior_value(   this->a_var,qp);
      Number DdamageDu= c.interior_value(   this->b_var,qp);

      // Now a loop over the fluence degrees of freedom.  This
      // computes the contributions of the continuity equation.
      for (unsigned int i=0; i != n_z_dofs; i++)
        {
          Fz(i) += JxW[qp] * (
             ( (  grad_z -this->m_MathModel.FluxSource(field_id,u_theta,
                                           damage,e_s0,this->m_PowerID) )
                     * dpsi[i][qp] ) 
                    + 
                   psi[i][qp] *  (
                   this->m_MathModel.FluenceReactionTerm(field_id,u_theta,damage,z_value)
                  -this->m_MathModel.IrradianceSource(field_id,u_theta,damage,e_value,this->m_PowerID)
                                 )
                             );
          if (request_jacobian)
            for (unsigned int j=0; j != n_z_dofs; j++)
              {
                // no coupling terms for linear system:w
                Kzz(i,j) += JxW[qp] * (
                ( dpsi[i][qp] * dpsi[j][qp] ) +
                   psi[i][qp] *  
                   this->m_MathModel.FluenceReactionTerm(field_id,u_theta,damage,psi[j][qp])
		                      ) ;
              }
        }
    } // end of the quadrature point qp-loop
  
  return request_jacobian;
}
// boundary conditions  pick up off diagonal terms
template< typename MathematicalModel  >
bool RHTESystem< MathematicalModel >::
side_constraint(bool request_jacobian,DiffContext &context)
{
  // call parent by default
  Parent::side_time_derivative(request_jacobian,context);

  // continue with off diagonal updates
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

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

  // The subvectors and submatrices we need to fill:
  //const unsigned int dim = this->get_mesh().mesh_dimension();
  DenseSubVector<Number> &Fz  = *c.elem_subresiduals[this->z_var];
  DenseSubMatrix<Number> &Kzz = *c.elem_subjacobians[this->z_var][this->z_var];
  DenseSubMatrix<Number> &Kuz = *c.elem_subjacobians[this->u_var][this->z_var];

  // Element Jacobian * quadrature weight for side integration
  const std::vector<Real> &JxW_side = c.side_fe_var[this->u_var]->get_JxW();
  
  // The shape functions at side quadrature points.
  const std::vector<std::vector<Real> >& phi_side = c.side_fe_var[this->u_var]->get_phi();
  const std::vector<std::vector<Real> >& psi_side = c.side_fe_var[this->z_var]->get_phi();
  
  // Physical location of the quadrature points on the side
  const std::vector<Point>& qpoint = c.side_fe_var[this->u_var]->get_xyz();
  
  // The number of local degrees of freedom in u and v
  const unsigned int n_u_dofs = c.dof_indices_var[this->u_var].size(); 
  const unsigned int n_z_dofs = c.dof_indices_var[this->z_var].size(); 

  // bc_id = 3 ==> cauchy
  // bc_id = 2 ==> neumann
  short int bc_id =
    this->get_mesh().boundary_info->boundary_id(c.elem, c.side);
  libmesh_assert (bc_id != BoundaryInfo::invalid_id);

    // residual
    for (unsigned int qp=0; qp != n_sidepoints; qp++)
      {
      Number u_theta  = c.side_value(   this->u_var,qp);
      Number z_value  = c.side_value(         this->z_var,qp);
        // loop over the fluence degrees of freedom.  
        if (request_jacobian)
        for (unsigned int i=0; i != n_z_dofs; i++)
          {
            Fz(i) += JxW_side[qp] * phi_side[i][qp] *
                     // add bc to residual
                     CALL_MEMBER_FN(this, this->ResidualBC[bc_id][this->z_var])(this->z_var,u_theta);
            if (request_jacobian)
              for (unsigned int j=0; j != n_z_dofs; j++)
                {
                  Kzz(i,j) += JxW_side[qp] * phi_side[i][qp] * phi_side[j][qp] *
                              // add bc to jacobian
                     CALL_MEMBER_FN(this, this->JacobianBC[bc_id][this->z_var])(this->z_var) ;
                  Kuz(i,j) += JxW_side[qp] * psi_side[i][qp] * phi_side[j][qp] * 1.0;
                  libmesh_not_implemented();
                  //this->m_MathModel.FluenceOffDiagonalBoundary(c,qp, field_id, 
                  //                                               qpoint[qp], es.parameters) ;
                }
          }
      } // end quadrature loop

  return request_jacobian;
}
// element assembly
template< typename MathematicalModel  >
void RHTESystem< MathematicalModel >::
UpdateBoundaryData(DiffContext &context,const unsigned int qp,
                       const int domainId,  const int timeId,
                       const Point &qpoint, const Parameters& parameters )
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
  Number u_theta  = c.side_value(this->u_var, qp);
  Number  damage  = c.side_value(   this->a_var,qp);
  //this->m_MathModel.UpdateBoundaryData(u_theta,damage,
  //                                     domainId,  timeId,
  //                                     qpoint, parameters );
  return;
}
