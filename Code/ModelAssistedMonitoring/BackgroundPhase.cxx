/* $Id: naviersystem.C 3039 2008-09-17 01:26:31Z benkirk $ */

/* The Next Great Finite Element Library. */
/* Copyright (C) 2003  Benjamin S. Kirk */

/* This library is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU Lesser General Public */
/* License as published by the Free Software Foundation; either */
/* version 2.1 of the License, or (at your option) any later version. */

/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/* Lesser General Public License for more details. */

/* You should have received a copy of the GNU Lesser General Public */
/* License along with this library; if not, write to the Free Software */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */

#include "getpot.h"


//libmesh
#include "equation_systems.h"
#include "boundary_info.h"
#include "fe_base.h"
#include "fe_interface.h"
#include "mesh.h"
#include "quadrature.h"
#include "steady_solver.h"
#include "petsc_diff_solver.h"
#include "dof_map.h"
#include "parallel.h"

//local
#include "pdeBaseClass.h"
#include "petsc_fem_context.h"
#include "BackgroundPhase.h"

const  PetscReal    fdEpsilon = PETSC_SQRT_MACHINE_EPSILON;

// constructor
BackgroundPhaseModel::BackgroundPhaseModel(GetPot &controlfile,
                                           EquationSystems &es):
  PDEModelBaseClass(controlfile,es),
  // magnetic Susceptibility
  m_Susceptibility(  
       "\\chi", es,
       controlfile("mrti/chi_optimize",  false),
       &PDEModelBaseClass::dpde_dchi,
       &PDEModelBaseClass::d2pde_du_dchi,
       controlfile("mrti/chi_lb"   ,    0.10e0),
       controlfile("mrti/chi_ub"   ,    0.70e0),
       controlfile("mrti/chi_dbeta", std::sqrt(fdEpsilon)), 0.5e0,  // TODO - cannot determine if this is a bug in adjoint or fd calc
       controlfile("mrti/chi_vary" ,   false) )
{
  m_GyroMagRatio  = 42.58 * 2.0 * libMesh::pi; // MHz/T *  rad
  m_ChemicalShift = controlfile("mrti/chemshift",0.0);
  m_ElectricPermeability  =   ///< [T*m/A]
      controlfile("mrti/permeability",4.0 * libMesh::pi*10.e-7);
  m_StaticFieldB0 = 1.5; // Tesla
  m_dSusceptibilitydz = 0.0;  // ppm / m
  // reserve space
  m_NeumannFlux.push_back(  0.0 ) ; 
  m_newton_coeff.push_back( 0.0 ) ; 
  m_u_infty.push_back(      0.0 ) ; 
  m_EchoTime      = controlfile("mrti/echotime",20.0);   // should be in millsec

  // use zero BC
  ResidualBC[2] = ResidualBC[1];  
  ResidualBC[3] = ResidualBC[1]; 
  JacobianBC[3] = JacobianBC[2];

  //field params
  es.parameters.set<PetscScalar>(   "susceptibility_healthy") = 
  controlfile(                 "mrti/susceptibility_healthy",0.0e0) ; 
  es.parameters.set<PetscScalar>(   "susceptibility_probe") = 
  controlfile(                 "mrti/susceptibility_probe",0.0e0) ; 
  es.parameters.set<PetscScalar>(   "susceptibility_tumor") = 
  controlfile(                 "mrti/susceptibility_tumor",
  es.parameters.get<PetscScalar>(   "susceptibility_healthy") ) ; 

  _fieldParameters.push_back( &m_Susceptibility  );  
}

/* -------------------------------------------------------------------- 
   Pennes Bioheat Equation with SDA laser is treated as the base class
   for the variety of laser source models
   -------------------------------------------------------------------- */ 
void BackgroundPhaseModel::printSelf(std::ostream& os)
{
  PetscFunctionBegin; 

  // print base class info
  this->PDEModelBaseClass::printSelf(os); 

  os << "Background: m_EchoTime             =" << m_EchoTime             << std::endl;
  os << "Background: m_StaticFieldB0        =" << m_StaticFieldB0        << std::endl;
  os << "Background: m_dSusceptibilitydz    =" << m_dSusceptibilitydz    << std::endl;
  os << "Background: m_ChemicalShift        =" << m_ChemicalShift        << std::endl;
  os << "Background: m_ElectricPermeability =" << m_ElectricPermeability << std::endl;
  os << "Background: m_GyroMagRatio         =" << m_GyroMagRatio         << std::endl;
  m_Susceptibility.printStdVector(     os, "  m_Susceptibility[");
  PetscFunctionReturnVoid(); 
}
// Constructor
template < typename MathematicalModel >
BackgroundPhaseSystem< MathematicalModel > ::
BackgroundPhaseSystem(EquationSystems& es, 
                      const std::string& name,
                      const unsigned int number) :
//base class constructor
  PetscFEMSystem                 (es, name, number),
  m_MathModel            (*(es.parameters.get<GetPot*>("controlfile")),
                             es) // initialize constitutive data
{
 this->m_jacobianComputed      = PETSC_FALSE;
}
template < typename MathematicalModel >
void BackgroundPhaseSystem< MathematicalModel > :: init_data ()
{
  const unsigned int dim = this->get_mesh().mesh_dimension();

  // phase data is not time evolving and is quasi static at each time step
  b_var = this->add_variable ("b", FIRST);

  this->time_solver = 
      AutoPtr<TimeSolver>(new SteadySolver(*this));

  // set nonlinear solver
  this->time_solver->diff_solver() =
        AutoPtr<DiffSolver>(new PetscDiffSolver(*this));

  // Do the parent's initialization after variables are defined
  Parent::init_data();

  // Useful debugging options
  // Set verify_analytic_jacobians to 1e-6 to use
  this->verify_analytic_jacobians = 0 ;
  this->print_jacobians = false ; 
  this->print_element_jacobians = false ; 

  // echo data
  this->m_MathModel.printSelf(std::cout); 

  // set default preconditioner
  PetscDiffSolver* SystemSolver = 
                  libmesh_cast_ptr<PetscDiffSolver*>(
                  & (*this->time_solver->diff_solver()) );
  PetscErrorCode info;
  PC pcPre;
  KSP  snesksp;
  info = SNESGetKSP(SystemSolver->snes(),&snesksp);CHKERRV(info);
  info = KSPSetType(snesksp,KSPCG); CHKERRV(info);
  info = KSPGetPC(snesksp,&pcPre); CHKERRV(info);
  info = PCSetType(pcPre,PCBJACOBI); CHKERRV(info);
}

template < typename MathematicalModel >
void BackgroundPhaseSystem< MathematicalModel > ::
init_context(DiffContext &context)
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // We should prerequest all the data
  // we will need to build the linear system.
  c.element_fe_var[b_var]->get_JxW();
  c.element_fe_var[b_var]->get_phi();
  c.element_fe_var[b_var]->get_dphi();
  c.element_fe_var[b_var]->get_xyz();
  
  c.side_fe_var[b_var]->get_JxW();
  c.side_fe_var[b_var]->get_phi();
  c.side_fe_var[b_var]->get_xyz();
}

template < typename MathematicalModel >
bool BackgroundPhaseSystem< MathematicalModel > ::
element_constraint ( bool request_jacobian, DiffContext &context)
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weight for interior integration
  const std::vector<Real> &JxW = c.element_fe_var[b_var]->get_JxW();

  // The shape function and gradients at interior
  // quadrature points.
  const std::vector<std::vector<Real> >& phi = 
    c.element_fe_var[b_var]->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi =
    c.element_fe_var[b_var]->get_dphi();

  // The number of local degrees of freedom in each variable
  const unsigned int n_b_dofs = c.dof_indices_var[b_var].size();

  // The subvectors and submatrices we need to fill:
  const unsigned int dim = this->get_mesh().mesh_dimension();
  DenseSubMatrix<Number> &Kpp = *c.elem_subjacobians[b_var][b_var];
  DenseSubVector<Number> &Fp = *c.elem_subresiduals[b_var];

  // Add the constraint given by the continuity equation
  unsigned int n_qpoints = c.element_qrule->n_points();
  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      // Compute the jacobian and residual for this LINEAR system
      Gradient grad_b = c.interior_gradient(b_var, qp);

      // Now a loop over the pressure degrees of freedom.  This
      // computes the contributions of the continuity equation.
      for (unsigned int i=0; i != n_b_dofs; i++)
        {
          // form the residual of this linear system
          Fp(i) += JxW[qp] * ( 
                dphi[i][qp] * grad_b + 
                 phi[i][qp] * this->m_MathModel.load()
                            );
          if (request_jacobian)
            for (unsigned int j=0; j != n_b_dofs; j++)
              {
                Kpp(i,j) += JxW[qp]*dphi[i][qp]*dphi[j][qp];
              }
        }
    } // end of the quadrature point qp-loop

  return request_jacobian;
}

//// Assemble System Dynamics Matrix for LITT Simulation
template< typename MathematicalModel  >
void BackgroundPhaseSystem< MathematicalModel >::
assembly(bool get_residual, bool get_jacobian)
{
 // update dirichlet BC before the assemble
 std::vector<PetscInt> &dirichletNodeSet= this->m_NodeSets[this->b_var][1];
 dirichletNodeSet.clear();

 EquationSystems &es = this->get_equation_systems();
 System &MaskSystem = es.get_system ("ImageMask");

 //  setup dirichlet data. 
 const DofMap & dof_map    = this->get_dof_map();

 //  hold variable dofs
 std::vector<unsigned int> dof_indices_var;

 // element based setup of dirichlet data
 libMesh::MeshBase &mesh = this->get_mesh();
 libMesh::MeshBase::const_element_iterator       el    =mesh.active_local_elements_begin();
 const libMesh::MeshBase::const_element_iterator end_el_full=mesh.active_local_elements_end();
 for ( ; el != end_el_full; ++el)
  {
   // Store element subdomain id to allow different sets of equations to be
   // solved on different parts of the mesh
   // ASSUMES SUBDOMAIN ORDERING STARTS FROM ONE but we need a 0 based
   // numbering scheme for std::vector
   Elem* elem = *el;

   for( unsigned int i_var = 0 ; i_var < this->n_vars() ; i_var++)
    {
     // degree of freedom indices for individual components
     dof_map.dof_indices (elem,   dof_indices_var ,i_var);
     // nodal based setup of dirichlet data
     // indices should be global
     for(unsigned int Ii = 0 ; Ii < dof_indices_var.size() ; Ii++)
      if( MaskSystem.current_solution( dof_indices_var[Ii] ) ) 
          dirichletNodeSet.push_back( dof_indices_var[Ii] );  
    }
  } // end loop over elements

 //  ensure all procs have the same dirichlet data or will freeze
 Parallel::allgather(dirichletNodeSet,false);
 std::cout << dirichletNodeSet.size() 
           << " dirichlet nodes detected" << std::endl;

 // be sure to reassemble Jacobian b/c dirichlet BC changing...
 this->m_jacobianComputed = PETSC_FALSE;
 // assemble
 this->Parent::assembly(get_residual,get_jacobian);

 return;
}
/**
 *  @mainpage Thermal Therapy ToolKit
 *  
 *  Thermal Therapy Toolkit is a suite of finite element solvers and image
 *  processing algorithms designed for 
 *  @ref TreatmentPlanning "Prospective Treatment Planning"
 *  of Laser Induced Thermal Therapy and Radio-Frequency Ablation.
 *  
 *  Input FEM Meshes and processed images are obtained from
 *    - Amira  http://www.amiravis.com
 *    - itkSNAP  http://www.itksnap.org
 *    - Cubit  http://cubit.sandia.gov
 *  
 *  The main C++ FEM routines are derived from 
 *    - libMesh  http://libmesh.sourceforge.net
 *  
 *  The main C++ Imaging routines are derived from 
 *    - ITK  http://www.itk.org
 *  
 *  Results may be visualized in 
 *    - Paraview  http://www.paraview.org
 *
 *  Routines for 
 *  @ref InverseSolver "Inverse Parameter Recovery" and 
 *  @ref ModelAssistedMonitoring "Model Assisted Monitoring"
 *  based on MR Temperature Imaging data are also available
 *  
 *  The main C++ optimization routines are derived from 
 *    - Petsc  http://www.mcs.anl.gov/petsc
 *    - Tao  http://www.mcs.anl.gov/tao
 *
 */
// template instantiations
template class BackgroundPhaseSystem<BackgroundPhaseModel>;
// Add an BackgroundSystem 
System * AddBackgroundSystem (        char *SystemName, 
                         EquationSystems *m_eqn_systems)
{
  PetscFunctionBegin; 

  FEMSystem *m_pdeSolver = static_cast< FEMSystem* > 
                                         ( &m_eqn_systems->add_system< 
          BackgroundPhaseSystem< BackgroundPhaseModel > >(SystemName)
                                         ) ; 

  m_pdeSolver->deltat = 1.0;
  PetscFunctionReturn(m_pdeSolver);

}
