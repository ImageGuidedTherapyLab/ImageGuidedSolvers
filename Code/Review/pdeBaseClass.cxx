// Various include files needed for the mesh & solver functionality.
#include "libmesh.h"
#include "equation_systems.h"
#include "fe.h"
#include "fe_interface.h"
#include "dof_map.h"
#include "quadrature_gauss.h"
#include "getpot.h"
#include "parallel.h"

// The nonlinear solver and system we will be using
#include "nonlinear_solver.h"
#include "nonlinear_implicit_system.h"
#include "linear_implicit_system.h"

#include "pdeBaseClass.h"
#include "petsc_fem_system.h"
#include "tttkUtilities.h"
#include "mesh.h"

PDEModelBaseClass::PDEModelBaseClass(GetPot &controlfile,EquationSystems &es):
  _equation_systems(es) // store the pointer

{
  // default is a linear PDE 
  m_LinearPDE = PETSC_TRUE;

  // default residual boundary conditions
  ResidualBC[0] = &PDEModelBaseClass::residualNothingBC;
  ResidualBC[1] = &PDEModelBaseClass::residualNothingBC;
  ResidualBC[2] = &PDEModelBaseClass::residualNeumannBC;
  ResidualBC[3] = &PDEModelBaseClass::residualCauchyBC;
  ResidualBC[4] = &PDEModelBaseClass::residualNothingBC;

  // default jacobian boundary conditions
  JacobianBC[0] = &PDEModelBaseClass::jacobianNothingBC;
  JacobianBC[1] = &PDEModelBaseClass::jacobianNothingBC;
  JacobianBC[2] = &PDEModelBaseClass::jacobianNothingBC;
  JacobianBC[3] = &PDEModelBaseClass::jacobianCauchyBC;
  JacobianBC[4] = &PDEModelBaseClass::jacobianNothingBC;

  // default jacobian boundary conditions (mainly for adjoint)
  //accumulateJacobianBC[0] = &PDEModelBaseClass::jacobianNothingBC;
  //accumulateJacobianBC[1] = &PDEModelBaseClass::jacobianNothingBC;
  //accumulateJacobianBC[2] = &PDEModelBaseClass::jacobianNothingBC;
  //accumulateJacobianBC[3] = &PDEModelBaseClass::jacobianCauchyBC;
  //accumulateJacobianBC[4] = &PDEModelBaseClass::jacobianNothingBC;
 
  // default adjoint boundary conditions
  accumulateAdjointLoadBC[0] = &PDEModelBaseClass::adjointLoadNothingBC;
  accumulateAdjointLoadBC[1] = &PDEModelBaseClass::adjointLoadNothingBC;
  accumulateAdjointLoadBC[2] = &PDEModelBaseClass::adjointLoadNothingBC;
  accumulateAdjointLoadBC[3] = &PDEModelBaseClass::adjointLoadCauchyBC;
  accumulateAdjointLoadBC[4] = &PDEModelBaseClass::adjointLoadNothingBC;

  // default gradient boundary conditions
  accumulateGradientBC[0] = &PDEModelBaseClass::gradientNothingBC;
  accumulateGradientBC[1] = &PDEModelBaseClass::gradientNothingBC;
  accumulateGradientBC[2] = &PDEModelBaseClass::gradientNeumannBC;
  accumulateGradientBC[3] = &PDEModelBaseClass::gradientCauchyBC;
  accumulateGradientBC[4] = &PDEModelBaseClass::gradientNothingBC;

}
/* -------------------------------------------------------------------- 
   Print BC
   -------------------------------------------------------------------- */ 
void PDEModelBaseClass::printSelf(std::ostream& os)
{
  os << "PDE.m_TimeDerivativeScalingFactor = " 
          << m_TimeDerivativeScalingFactor << std::endl;
  printStdVector< PetscTruth >(os, "PDE.m_TransientDomain[" , m_TransientDomain);

}

#undef __FUNCT__
#define __FUNCT__ "PDEModelBaseClass::writeElementData"
void PDEModelBaseClass::writeElementData(OStringStream &file_name, 
                       libMesh::MeshBase &mesh, std::vector<Number>& soln,Real ElemTime)
{
  PetscFunctionBegin;
  // This function must be run on all processors at once
  parallel_only();

  libmesh_not_implemented();
  ////const unsigned int dim = mesh.mesh_dimension();
  //const unsigned int ne  = mesh.n_elem();

  ////return if no field variables to plot
  //const unsigned int nv  = _fieldParameters.size();
  //if(!nv) PetscFunctionReturnVoid();

  ////get number of field variables and name of each variable
  //std::vector<std::string> names;
  //std::vector<optimizationParameter*>::iterator locParamIter;
  //for(locParamIter  = _fieldParameters.begin(); 
  //    locParamIter != _fieldParameters.end()  ; locParamIter++)
  //  {
  //   optimizationParameter* optParam = *locParamIter;
  //   // spatially varying parameter
  //   names.push_back(optParam->name()); 
  //  }

  //// We'd better have a contiguous node numbering
  //libmesh_assert (ne == mesh.max_elem_id());

  //// allocate storage to hold
  //// (number_of_nodes)*(number_of_variables) entries.
  //soln.resize(ne*nv);

  //// Zero out the soln vector
  //std::fill (soln.begin(),       soln.end(),       libMesh::zero);
  //
  //// For each system in this EquationSystems object,
  //// update the global solution and if we are on processor 0,
  //// loop over the elements and build the nodal solution
  //// from the element solution.  Then insert this nodal solution
  //// into the vector passed to build_solution_vector.

  //libMesh::MeshBase::const_element_iterator       it  = mesh.active_local_elements_begin();
  //const libMesh::MeshBase::const_element_iterator end = mesh.active_local_elements_end(); 
  //for ( ; it != end; ++it)
  //  {
  //   // element iterator
  //   const Elem* elem = *it;

  //   // get the field id in the spatially varying data structures
  //   PetscInt FieldId =  elem->_mat_data_field;

  //   // counter used for storage
  //   for(locParamIter  = _fieldParameters.begin(); 
  //       locParamIter != _fieldParameters.end(); locParamIter++)
  //     {
  //      optimizationParameter* optParam = *locParamIter;
  //      const int paramID = distance(_fieldParameters.begin(),locParamIter);
  //      // get parameter to plot 
  //      optParam->getParam(soln[ nv*elem->id() + paramID ], FieldId);
  //     }
  //  }	 

  //// Now each processor has computed contriburions to the
  //// soln vector.  Gather them all up.
  //Parallel::sum(soln);

  //// Write out the initialize vis of the field
  //output_params.write_elem_data(file_name.str(),soln,names);
  //output_params.write_timestep(ElemTime);
  //// increment timestep counter
  //output_params.increment_timestep();

  PetscFunctionReturnVoid();
}
/*----------------------------------------------------------------------*/
void PDEModelBaseClass::SetupOptimizationVariables( libMesh::MeshBase &mesh,
                   std::vector<optimizationParameter*> &Parameters)
{
 PetscFunctionBegin;


 PetscFunctionReturnVoid();
}
/* ------------------------------------------------------------------- 
     Get consitutivie parameters from exodus file on disk
   ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "PDEModelBaseClass::GetDataFromDisk"
PetscErrorCode PDEModelBaseClass::GetDataFromDisk(libMesh::MeshBase &mesh,std::string &FieldFile)
{
  PetscFunctionBegin; 

  libmesh_not_implemented();
  //// application context
  //AppSolve *user = AppSolve::getInstance();

  //// Open the exodus file and get the header info 
  //ExodusII_IO_Helper field_data;
  //field_data.open( FieldFile.c_str() ); // Open the exodus file, if possible
  //field_data.read_header();        // Get header information from exodus file
  //field_data.read_block_info();    // Get information about all the blocks
  //field_data.verbose(true);        // Set verbose
  //// default to final time
  //const unsigned int idtime = field_data.num_time_steps;
  //if( user->get_num_elem_blk() != field_data.get_num_elem_blk() )
  //  {
  //    std::cout << "# of mesh block do not agree " << std::endl << std::flush; 
  //    libmesh_error();
  //  }

  //// storage for element solution
  //std::vector<double> elem_values;  

  ////build subdomain map to retrieve data (SERIAL)
  //std::map<unsigned int, std::vector< unsigned int >  > subdomain_map;
  //for(unsigned int i=0; i<static_cast<unsigned int>( mesh.n_elem() ); i++)
  // {
  //  Elem * elem = mesh.elem(i);
  //  unsigned int cur_subdomain = elem->subdomain_id();
  //  subdomain_map[cur_subdomain].push_back(elem->id());
  // }
  //// loop over subdomains
  //std::vector< int > elem_id_map( mesh.n_elem() ,-1);  
  //int n_DomainElemPrev =  0 ; 
  //for( int idBlk = 0 ; idBlk < user->get_num_elem_blk() ; idBlk++ )
  // { 
  //    std::vector< unsigned int > & tmp_vec = subdomain_map[idBlk+1];
  //    int n_DomainElem =  tmp_vec.size() ; 
  //    // error check the # of elements per domain
  //    field_data.read_elem_in_block (idBlk);
  //    if( field_data.get_num_elem_this_blk() != n_DomainElem )
  //      {
  //        std::cout << "# of elem in mesh block " << idBlk
  //                  << " do not agree " << std::endl << std::flush; 
  //        libmesh_error();
  //      }

  //    for( int Ii = 0 ; Ii < n_DomainElem ; Ii++ ) 
  //                  elem_id_map.at( tmp_vec[Ii] ) = n_DomainElemPrev + Ii ; 
  //    n_DomainElemPrev += n_DomainElem  ; 
  // }

  //// initialize iterator
  //std::vector<optimizationParameter*>::iterator locParamIter;
  //const libMesh::MeshBase::const_element_iterator el_end = mesh.elements_end();

  //// loop over all parameters
  //for( locParamIter  = _fieldParameters.begin();
  //     locParamIter != _fieldParameters.end()  ; locParamIter++)
  // {
  //   optimizationParameter* optParam = *locParamIter ;
  //   if( field_data.get_elem_var_values(optParam->name(),idtime,elem_values)  )
  //    {
  //     // error check
  //     if(mesh.n_elem() != elem_values.size()) 
  //       {
  //        std::cout << "field mesh doesn't have same number of elements" 
  //                  << std::endl << std::flush ;
  //        libmesh_error();
  //       }
  //     libMesh::MeshBase::const_element_iterator   el     = mesh.elements_begin();
  //     for ( ; el != el_end; ++el)
  //       {
  //         Elem* elem = *el;
  //         #if defined(PETSC_USE_DEBUG)
  //         optParam->at(elem->_mat_data_field) = elem_values.at(
  //                                               elem_id_map.at( elem->id() ) );
  //         #else
  //         (*optParam)[elem->_mat_data_field]  = elem_values[ 
  //                                               elem_id_map[ elem->id() ] ];
  //         //#endif
  //       }
  //    }
  // }

  PetscFunctionReturn(0);
}


/*----------------------------------------------------------------------*/
void PDEModelBaseClass::adjointLoadCauchyBC(const unsigned int i_var,QGauss &qface,
                             const std::vector<std::vector<Real> >&  psi_face, 
                             const std::vector<Real>& JxW_face ,
                             std::vector< DenseSubVector<Number> > & Fi,
                             TransientFEMSystem &system)
{
 PetscFunctionBegin; 

 libmesh_not_implemented();
 //// Loop over the face quadrature points for integration.
 //for (unsigned int qp=0; qp<qface.n_points(); qp++)
 //  {
 //   Real p_i_future=0.0; //initialize solution value (if needed)
 //   for (unsigned int i=0; i<psi_face.size(); i++)
 //      p_i_future += psi_face[i][qp]*system.stored_solution(dof_indices_p[i_var][i],AppSolve::ISTEP+1);
 //   for (unsigned int i=0; i<psi_face.size(); i++)
 //     {
 //         Fi[i_var](i) += JxW_face[qp] * m_theta * m_newton_coeff[i_var]
 //                               * p_i_future * psi_face[i][qp]; 
 //     }
 //  } // end loop over bc guass points 
 PetscFunctionReturnVoid(); 
}
/*----------------------------------------------------------------------*/
void PDEModelBaseClass::gradientNeumannBC(const unsigned int i_var,QGauss &qface,
                         const std::vector<std::vector<Real> >&  , 
                         const std::vector<std::vector<Real> >&  psi_face, 
                         const std::vector<Real>& JxW_face ,
                         DenseVector<Number> &Grad,
                         TransientFEMSystem &,
                         TransientFEMSystem &adjoint_system)
{
 PetscFunctionBegin; 
 /* get the boundary data, boundary condition 
    function pointers to one of the following routines
                  getverifbc     : verification problem bc
                  rfvoltagebc    : bc for rf voltage
                  heattransbc    : bc for heat transfer */
 //NOTE THE PLUS SIGN
 libmesh_not_implemented();
 //for (unsigned int qp=0; qp<qface.n_points(); qp++)
 //     for (unsigned int i=0; i<psi_face.size(); i++)
 //       Grad(0) += JxW_face[qp]*m_NeumannFlux[i_var]*psi_face[i][qp] * 
 //                 adjoint_system.current_solution(dof_indices_p[i_var][i]);
 PetscFunctionReturnVoid(); 
}
/*----------------------------------------------------------------------*/
//void PDEModelBaseClass::jacobianCauchyBC(
//                    const unsigned int i_var, QGauss &qface,
//                    const std::vector<Real>& JxW_face ,
//                    const std::vector<std::vector<Real> >& phi_face,
//                    std::vector< SubMatrixVector > &Kij)
//{
// PetscFunctionBegin;
//
// // Loop over the face quadrature points for integration.
// for (unsigned int qp=0; qp<qface.n_points(); qp++)
//    for (unsigned int i=0; i<phi_face.size(); i++)
//      for (unsigned int j=0; j<phi_face.size(); j++)
//        Kij[i_var][i_var](i,j) += JxW_face[qp] * m_theta *
//                m_newton_coeff[i_var] * phi_face[i][qp]*phi_face[j][qp];
//
// PetscFunctionReturnVoid();
//}
///*----------------------------------------------------------------------*/
void PDEModelBaseClass::gradientCauchyBC(const unsigned int i_var,QGauss &qface,
                          const std::vector<std::vector<Real> >&  phi_face, 
                          const std::vector<std::vector<Real> >&  psi_face, 
                          const std::vector<Real>& JxW_face ,
                          DenseVector<Number> &Grad,
                          TransientFEMSystem &state_system,
                          TransientFEMSystem &adjoint_system)
{
 PetscFunctionBegin; 

 libmesh_not_implemented();
 //// Loop over the face quadrature points for integration.
 //for (unsigned int qp=0; qp<qface.n_points(); qp++)
 // {
 //   // Compute the solution at the theta timestep
 //   Real usln_theta=0.0; //initialize solution value 
 //   for (unsigned int i=0; i<phi_face.size(); i++)
 //       usln_theta +=      m_theta  * phi_face[i][qp]*
 //                   state_system.stored_solution(dof_indices_u[i_var][i],
 //                                                AppSolve::ISTEP) 
 //                                  +
 //                   (1.0-m_theta) * phi_face[i][qp]*
 //                   state_system.stored_solution(dof_indices_u[i_var][i],
 //                                                AppSolve::ISTEP-1);
 //   for (unsigned int i=0; i<psi_face.size(); i++)
 //     Grad(0) += JxW_face[qp]*psi_face[i][qp]*
 //                adjoint_system.stored_solution(dof_indices_p[i_var][i],
 //                                               AppSolve::ISTEP) *
 //                     m_newton_coeff[i_var]*(usln_theta-m_u_infty[i_var]) ;
 // }
 PetscFunctionReturnVoid(); 
}
