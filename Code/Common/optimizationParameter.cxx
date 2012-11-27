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
#include "petsc_vector.h"
#include "sparse_matrix.h"
#include "elem.h"
#include "string_to_enum.h"
#include "getpot.h"
#include "boundary_info.h"
#include "mesh.h"

// libmesh includes
#include "numeric_vector.h"
#include "dense_vector.h"

// pdeopt
#include "optimizationParameter.h"


/* --------------------- base class constructor -------------------------- */
optimizationParameter::optimizationParameter(const char *parameterName,
                       bool optimize, 
                       //dpde_dmMemFn dpde_dmMemberFunctor,
                       //d2pde_dudmMemFn  d2pde_dudumMemberFunctor,
                       PetscScalar lb, PetscScalar ub,
                       PetscScalar perturbation, 
                       PetscScalar verifparam):_name(parameterName)
{
    time_vary     = false;
    spatial_field = false;
    //dqoi_dm       = &qoiBaseClass::dqoi_dm;
    dqoi_dm       = NULL; 
    dpde_dm       = NULL; //dpde_dmMemberFunctor;
    d2qoi_du_dm   = NULL; //&PDEModelBaseClass::d2qoi_du_dm; 
    d2pde_du_dm   = NULL; //d2pde_dudumMemberFunctor;
    Optimize      = optimize;
    _verifparam   = verifparam;
    _perturbation  = perturbation;
    // initialize with at least enough room for each domain
    // FIXME 10 domains should be plenty... nice hard code...
    _value_lb.resize(10, lb);
    _value_ub.resize(10, ub);
}
/* --------------------- constructor -------------------------- */
spatialParameter::spatialParameter( const char *parameterName,
           EquationSystems &,
           bool opt, 
           //dpde_dmMemFn dpde_dmMemberFunctor,
           //d2pde_dudmMemFn  d2pde_dudumMemberFunctor,
           PetscScalar lb, PetscScalar ub,
           PetscScalar perturbation, PetscScalar verifparam, bool field )
           : optimizationParameter(parameterName,opt,
                                   //dpde_dmMemberFunctor,
                                   //d2pde_dudumMemberFunctor,
                                   lb,ub, perturbation,verifparam)
{
    // initialize pointer to data
    m_data                = NULL; 
    m_globalScatter       = NULL ;
    m_LocalGlobalSolution = NULL ;
    spatial_field = field;
}
// clean up
spatialParameter::~spatialParameter()
{
  // destroy scatter context and local vector when no longer needed
  if(m_globalScatter) VecScatterDestroy(m_globalScatter);
  if(m_LocalGlobalSolution) VecDestroy(m_LocalGlobalSolution);
}
/* ------------------------------------------------------------------- */
discreteParameter::discreteParameter( const char *parameterName,
           bool opt, 
           //dpde_dmMemFn dpde_dmMemberFunctor,
           //d2pde_dudmMemFn  d2pde_dudumMemberFunctor,
           PetscScalar lb, PetscScalar ub,
           PetscScalar perturbation, PetscScalar verifparam)
           : optimizationParameter(parameterName,opt,
                                   //dpde_dmMemberFunctor,
                                   //d2pde_dudumMemberFunctor,
                                   lb,ub, perturbation,verifparam),
             _value()
{
  time_vary = true;
}
// The matrix assembly function to be called at each time step to
// prepare for the linear solve.
void initialize_field_parameters (EquationSystems& es,
                                 const std::string& field_name)
{
  PetscFunctionBegin; 

  // Get a constant reference to the mesh object.
  const libMesh::MeshBase& mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();
  
  // Get a reference to the systems objects.
  ExplicitSystem & FieldSystem =
    es.get_system<ExplicitSystem> (field_name);

  // Numeric ids corresponding to each variable in the system
  const unsigned int f_var = FieldSystem.variable_number (field_name);
  
  // Get the Finite Element type for "f".
  FEType fe_type = FieldSystem.variable_type(f_var);
  
  // Build a Finite Element object of the specified type for
  // the imaging variables.
  AutoPtr<FEBase> fe_base  (FEBase::build(dim, fe_type));
    
  // A Gauss quadrature rule for numerical integration.
  // Let the \p FEType object decide what order rule is appropriate.
  QGauss qrule (dim, fe_type.default_quadrature_order());

  // Tell the finite element objects to use our quadrature rule.
  fe_base->attach_quadrature_rule (&qrule);
  
  // A reference to the \p DofMap object for this system.  The \p DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the \p DofMap
  // in future examples.
  const DofMap & dof_map = FieldSystem.get_dof_map();

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<unsigned int> dof_indices;

  // dense vector for storing solution
  DenseVector<Number> Ue;

  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix and right-hand-side contribution.  Since the mesh
  // will be refined we want to only consider the ACTIVE elements,
  // hence we use a variant of the \p active_elem_iterator.
  libMesh::MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 
  
  // loop over elements and compute signal
  for ( ; el != end_el; ++el)
    {    
     // Store a pointer to the element we are currently
     // working on.  This allows for nicer syntax later.
     const Elem* elem = *el;
     
     const int subdomain_id = elem->subdomain_id();
     // Get the degree of freedom indices for the
     // current element.  These define where in the global
     // matrix and right-hand-side this element will
     // contribute to.
     dof_map.dof_indices (elem, dof_indices);

     const unsigned int n_dofs   = dof_indices.size();
     
     // Compute the element-specific data for the current
     // element.  This involves computing the location of the
     // quadrature points (q_point) and the shape functions
     // (phi, dphi) for the current element.
     fe_base->reinit  (elem);

     // Zero the element matrix and right-hand side before
     // summing them.  We use the resize member here because
     // the number of degrees of freedom might have changed from
     // the last element.  Note that this will be the case if the
     // element type is different (i.e. the last element was a
     // triangle, now we are on a quadrilateral).
     Ue.resize (n_dofs);

     // should be piece wise constant across elements
     libmesh_assert(qrule.n_points() == 1);
     libmesh_assert(n_dofs == 1);
     // calculate signal  at each quadrature point by summing the
     // solution degree-of-freedom values by the appropriate
     // weight functions.
     for (unsigned int qp=0; qp<qrule.n_points(); qp++)
       {
        
        // Compute the velocity & its gradient from the previous timestep
        // and the old Newton iterate.
        for (unsigned int l=0; l<n_dofs; l++)
          {
            std::ostringstream var_name;
            if(subdomain_id == 2)
               var_name<< field_name <<  "_tumor";
            else if(subdomain_id >= 3)
               var_name<< field_name <<  "_probe";
            else 
               var_name<< field_name <<  "_healthy";
            Ue(l) =  es.parameters.get<PetscScalar>( var_name.str() ) ;
          }

       } // end of the quadrature point qp-loop

     FieldSystem.solution->insert (Ue,dof_indices); 
    } // end of element loop

  FieldSystem.solution->close();
  // copy parallel data structures to local data structures
  FieldSystem.solution->localize(*FieldSystem.current_local_solution);

  PetscFunctionReturnVoid(); 
}
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "spatialParameter::OneVarSetup"
PetscErrorCode spatialParameter::OneVarSetup(EquationSystems &es)
{
  PetscFunctionBegin; 

  m_data = &(es.add_system<ExplicitSystem>(_name));
  m_data->add_variable (_name,CONSTANT,MONOMIAL );
  m_data->attach_init_function ( initialize_field_parameters  );
  //m_data->init ( ); 
  //FiniteElementInterface *user = FiniteElementInterface::getInstance();
  //// for simplicity ALWAYS count all constant parameters in the map
  //for(PetscInt iii=0; iii < user->get_num_elem_blk(); iii++)  dofs.push_back(iii);

  //libMesh::MeshBase::const_element_iterator el     = mesh.active_local_elements_begin();
  //const libMesh::MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
  //// this parameter's spatially and constant portions of the  field
  //// should be implicitly set in  InitMeshDomains
  //// a constant portion of the field should ALWAYS exist  in subdomain 0
  //if( this->spatial_field ) 
  // for ( ; el != el_end; ++el)
  //  {
  //     // Store a pointer to the element we are currently
  //     // working on.  This allows for nicer syntax later.
  //     Elem* elem = *el;
  //     //dofs.push_back( elem->_mat_data_field );
  //     std::cout << "fix numbering" <<std::endl;
  //     libmesh_error();
  //  }

  //// sort all element maps note that for a field/constant mix the first 
  //// baseInfo::n_block dofs corresponds to the constant portion
  //std::sort(dofs.begin(),dofs.end());
  //// remove all duplicates
  //std::vector<PetscInt>::iterator pos;
  //pos = std::unique(dofs.begin(),dofs.end());
  //// have to explicitly erase duplicates
  //dofs.erase(pos,dofs.end()); 

  //// should have at least one dof to optimize on multiple procs may be less than
  //// baseInfo::n_block bc of the partition
  //libmesh_assert( dofs.size() >= 
  //                static_cast<unsigned int>(user->get_num_elem_blk()) );  

  PetscFunctionReturn(0);
}
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "spatialParameter::ScatterGlobalToAll"
PetscErrorCode spatialParameter::ScatterGlobalToAll()
{
  PetscFunctionBegin; 

  Vec GlobalSoln = (dynamic_cast< PetscVector<double>* > 
                    (&(*m_data->solution)) )->vec();

  // need a local copy of the global solution vec
  // create the scater context if it hasn't been created yet
  if(m_globalScatter == NULL && m_LocalGlobalSolution == NULL) 
    { VecScatterCreateToAll(GlobalSoln,&m_globalScatter,&m_LocalGlobalSolution); }

  // scater global vector to local copy
  VecScatterBegin(m_globalScatter,GlobalSoln,m_LocalGlobalSolution,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(  m_globalScatter,GlobalSoln,m_LocalGlobalSolution,INSERT_VALUES,SCATTER_FORWARD);

  PetscFunctionReturn(0);
}
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "spatialParameter::GetGlobalSolution"
PetscScalar spatialParameter::GetGlobalSolution(const unsigned int globalID)
{
  PetscErrorCode ierr;
  PetscFunctionBegin; 

  PetscScalar  GlobalValue;
  const PetscInt globalIndex = globalID;
  ierr = VecGetValues(m_LocalGlobalSolution,1,&globalIndex,&GlobalValue);

  PetscFunctionReturn(GlobalValue);
}
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "discreteParameter::OneVarSetup"
PetscErrorCode discreteParameter::OneVarSetup(EquationSystems &)
{
  PetscFunctionBegin; 

  // create the map from the local contributions to the local gradient
  // parameter varies in time but NOT space

  PetscInt power_size = _value.size();
  for(PetscInt iii= 0 ; iii < power_size ; iii ++)  dofs.push_back( iii );

  PetscFunctionReturn(0);
}
