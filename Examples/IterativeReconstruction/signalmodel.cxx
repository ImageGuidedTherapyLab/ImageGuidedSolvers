/* $Id: ex18.C 3043 2008-09-17 15:18:58Z benkirk $ */

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

 // <h1>Example 18 - Unsteady Navier-Stokes Equations with DiffSystem</h1>
 //
 // This example shows how the transient nonlinear problem from
 // example 13 can be solved using the new (and experimental)
 // DiffSystem class framework

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>

// Basic include files
#include "error_vector.h"
#include "getpot.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "mesh_refinement.h"
#include "uniform_refinement_estimator.h"
#include "fe.h"
#include "quadrature_gauss.h"
#include "dof_map.h"
#include "explicit_system.h"
#include "parallel.h"
#include "gmv_io.h"
#include "exodusII_io.h"

// Some (older) compilers do not offer full stream
// functionality, OStringStream works around this.
#include "o_string_stream.h"

// local class definition
#include "signalmodel.h"

SignalModelInterface::SignalModelInterface()
{
    m_mesh = NULL;
    m_eqn_systems = NULL;
    m_init = NULL;
    m_PhaseModel = NULL;
    m_echoTime = 0.0;
    m_InputGradientModel = UNKNOWNGRADIENT  ;
}

SignalModelInterface::~SignalModelInterface()
{ // clean up pointers
  if(m_mesh)        delete m_mesh;
  if(m_eqn_systems) delete m_eqn_systems;
  if(m_init)        delete m_init;
  if(m_PhaseModel ) delete m_PhaseModel ;
}

// print local info
void SignalModelInterface::printSelf()
{
  PetscFunctionBegin; 
  std::cout << "echo time = " << m_echoTime << std::endl;
  // Print information about the mesh to the screen.
  if(m_mesh) m_mesh->print_info();

  // Print information about the system to the screen.
  if(m_eqn_systems) m_eqn_systems->print_info();
  PetscFunctionReturnVoid(); 
}

PetscErrorCode SignalModelInterface::SetuplibMesh(MPI_Comm libMesh_COMM_IN)
{
  PetscFunctionBegin; 
  // initialize libMesh
  int argc = 1;
  char *notused[] = {"not","used..."};
  char **argv = notused;
  std::cout << "setting up libMesh" << std::endl;
  m_init = new LibMeshInit(argc, argv,libMesh_COMM_IN);
  PetscFunctionReturn(0);
}

PetscErrorCode SignalModelInterface::SetupImaging( 
                                       int nx_image, int ny_image, int nz_image,
				       double X0, double Y0, double Z0, 
                                       double DX, double DY, double DZ )
{
  PetscFunctionBegin; 

  // setup import filters
  m_importFilter = ImportFilterType::New();      

  // pixel dimensions
  m_size[0]  = nx_image;  // size along X
  m_size[1]  = ny_image;  // size along Y
  m_size[2]  = nz_image;  // size along Z

  ImportFilterType::IndexType start;
  start.Fill( 0 );

  ImportFilterType::RegionType region;
  region.SetIndex( start );
  region.SetSize(  m_size  );

  m_importFilter->SetRegion( region );

  // origin 
  m_origin[0] = X0;    // X coordinate 
  m_origin[1] = Y0;    // Y coordinate
  m_origin[2] = Z0;    // Z coordinate

  m_importFilter->SetOrigin( m_origin );

  // spacing
  m_spacing[0] = DX;    // along X direction 
  m_spacing[1] = DY;    // along Y direction
  m_spacing[2] = DZ;    // along Z direction

  m_importFilter->SetSpacing( m_spacing );

  m_interpolator = InterpolatorType::New();
  if(!m_eqn_systems)
    {
      std::cerr << "must call SetupStructuredGrid first " << std::endl;
      libmesh_error();
    }
  m_eqn_systems->parameters.set< InterpolatorType::Pointer 
                               >("interpolator") = m_interpolator;

  PetscFunctionReturn(0);
}

PetscErrorCode SignalModelInterface::SetupLinearGradientModel(
                double x_Gradient, double y_Gradient, double z_Gradient,
                double larmor)
{
  PetscFunctionBegin; 

  if(m_PhaseModel)
    {
      std::cerr << "gradient model already setup" << std::endl;
      libmesh_error();
    }
  m_InputGradientModel = LINEARGRADIENT;
  Point linearGradient(x_Gradient, y_Gradient, z_Gradient);
  m_PhaseModel = new LinearGradientModel(linearGradient,larmor);
  PetscFunctionReturn(0);
}
PetscErrorCode SignalModelInterface::SetupSpiralGradientModel(
                double x_Gradient, double y_Gradient, double z_Gradient,
                double larmor,double angularVel)
{
  PetscFunctionBegin; 

  if(m_PhaseModel)
    {
      std::cerr << "gradient model already setup" << std::endl;
      libmesh_error();
    }
  m_InputGradientModel = SPIRALGRADIENT; 
  Point linearGradient(x_Gradient, y_Gradient, z_Gradient);
  m_PhaseModel = new SpiralGradientModel(linearGradient,larmor,angularVel);

  PetscFunctionReturn(0);
}
PetscErrorCode SignalModelInterface::SetupEPIGradientModel(
                double x_Gradient, double y_Gradient, double z_Gradient,
                double larmor,int NSkipKSpace)
{
  PetscFunctionBegin; 

  if(m_PhaseModel)
    {
      std::cerr << "gradient model already setup" << std::endl;
      libmesh_error();
    }
  m_InputGradientModel = EPIGRADIENT;
  Point linearGradient(x_Gradient, y_Gradient, z_Gradient);
  m_PhaseModel = new EPIGradientModel(linearGradient,larmor,NSkipKSpace);
  PetscFunctionReturn(0);
}

PetscErrorCode SignalModelInterface::SetupStructuredGrid( 
                                       int nx_mesh, int ny_mesh, int nz_mesh,
				       double xmin, double xmax, 
				       double ymin, double ymax, 
				       double zmin, double zmax )
{
  PetscFunctionBegin; 

  // possible refinements
  const unsigned int coarserefinements = 0;

  libmesh_assert (m_Dimension == 3);

  // Create a n-dimensional mesh.
  m_mesh = new Mesh(m_Dimension);
  
  // And an object to refine it
  MeshRefinement mesh_refinement(*m_mesh);
  mesh_refinement.coarsen_by_parents() = true;
  //mesh_refinement.absolute_global_tolerance() = global_tolerance;
  //mesh_refinement.nelem_target() = nelem_target;
  mesh_refinement.refine_fraction() = 0.3;
  mesh_refinement.coarsen_fraction() = 0.3;
  mesh_refinement.coarsen_threshold() = 0.1;

  // Use the MeshTools::Generation mesh generator to create a uniform
  // grid on the square [-1,1]^D.  We instruct the mesh generator
  // to build a mesh of 8x8 \p Quad9 elements in 2D, or \p Hex27
  // elements in 3D.  Building these higher-order elements allows
  // us to use higher-order approximation, as in example 3.
  MeshTools::Generation::build_cube (*m_mesh,
                                       nx_mesh,ny_mesh,nz_mesh,
				       xmin, xmax, 
				       ymin, ymax, 
				       zmin, zmax,
                                       HEX8);

  mesh_refinement.uniformly_refine(coarserefinements);

  // Create an equation systems object.
  m_eqn_systems = new EquationSystems(*m_mesh);

  PetscFunctionReturn(0);
}

// Add an ExplicitSystem to store data
PetscErrorCode SignalModelInterface::AddExplicitSystem (char *SystemName,int VarOrder)
{
  PetscFunctionBegin; 

  if(!m_eqn_systems) 
    {
      std::cerr << "must call SetupStructuredGrid first " << std::endl;
      libmesh_error();
    }
  ExplicitSystem & system = 
    m_eqn_systems->add_system<ExplicitSystem> (SystemName);

  // use first order variable
  std::ostringstream variable_name;
  variable_name << SystemName << 0  ;

  switch(VarOrder)
   {
    case 0: 
      system.add_variable (variable_name.str(), CONSTANT,MONOMIAL );
      break;
    case 1: 
      system.add_variable (variable_name.str(), FIRST,   LAGRANGE );
      break;
    default: 
      std::cerr << "unexpected variable order" << std::endl;
      PetscFunctionReturn(-9);
   }

  // use coil sensitivity system to hold the real and imaginary parts of the coil sensitivity
  // in the solution vector. use rhs to hold the forward projected signal.
  PetscFunctionReturn(0);
}

PetscErrorCode SignalModelInterface::InitializeEquationSystems()
{
  PetscFunctionBegin; 
  // Initialize the system
  m_eqn_systems->init ();
  PetscFunctionReturn(0);
}

// Write out the FEM data
PetscErrorCode SignalModelInterface::WriteEquationSystems(char *FileName)
{
  PetscFunctionBegin; 
  const std::string file_name(FileName);
  ExodusII_IO(*m_mesh).write_equation_systems (file_name,
                                              *m_eqn_systems);
  //GMVIO(*m_mesh).write_equation_systems (file_name,
  //                                            *m_eqn_systems);
  PetscFunctionReturn(0);
}


// template to allow functionoids
template< typename GradientModel  >
void SignalAssembly(EquationSystems &es,
                    GradientModel   &PhaseModel,
                    std::vector<double> &TimeArray,
                    std::vector<double> &RealSignal,
                    std::vector<double> &ImagSignal)
{
  PetscFunctionBegin; 
  
  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es.get_mesh();
  
  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();
  
  // Get a reference to the systems objects.
  ExplicitSystem & SpinSystem =
    es.get_system<ExplicitSystem> ("Spin");
  ExplicitSystem & T2StarSystem =
    es.get_system<ExplicitSystem> ("T2Star");
  ExplicitSystem & RealCoilSensitivitySystem =
    es.get_system<ExplicitSystem> ("RealCoilSensitivity");
  ExplicitSystem & ImagCoilSensitivitySystem =
    es.get_system<ExplicitSystem> ("ImagCoilSensitivity");

  // Numeric ids corresponding to each variable in the system
  const unsigned int f_var = SpinSystem.variable_number ("Spin0");
  
  // Get the Finite Element type for "f".
  FEType fe_type = SpinSystem.variable_type(f_var);
  
  // Build a Finite Element object of the specified type for
  // the imaging variables.
  AutoPtr<FEBase> fe_base  (FEBase::build(dim, fe_type));
    
  // A Gauss quadrature rule for numerical integration.
  // Let the \p FEType object decide what order rule is appropriate.
  QGauss qrule (dim, fe_type.default_quadrature_order());

  // Tell the finite element objects to use our quadrature rule.
  fe_base->attach_quadrature_rule (&qrule);
  
  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  //
  // The element Jacobian * quadrature weight at each integration point.   
  const std::vector<Real>& JxW = fe_base->get_JxW();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> >& phi = fe_base->get_phi();

  // The element shape function gradients for the velocity
  // variables evaluated at the quadrature points.
  const std::vector<std::vector<RealGradient> >& dphi = fe_base->get_dphi();

  // The physical XY locations of the quadrature points on the element.
  // These are useful for evaluating spatially varying material
  // properties at the quadrature points.
  const std::vector<Point>& q_point = fe_base->get_xyz();

  // A reference to the \p DofMap object for this system.  The \p DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the \p DofMap
  // in future examples.
  const DofMap & dof_map = SpinSystem.get_dof_map();

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<unsigned int> dof_indices;

  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix and right-hand-side contribution.  Since the mesh
  // will be refined we want to only consider the ACTIVE elements,
  // hence we use a variant of the \p active_elem_iterator.
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 
  
  // initialize the signal to be returned
  RealSignal.resize(TimeArray.size(),0.0);
  ImagSignal.resize(TimeArray.size(),0.0);
  
  std::cout << "assembling signal " << std::endl;
  // loop over elements and compute signal
  for ( ; el != end_el; ++el)
    {    
     // Store a pointer to the element we are currently
     // working on.  This allows for nicer syntax later.
     const Elem* elem = *el;
     
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

     // calculate signal  at each quadrature point by summing the
     // solution degree-of-freedom values by the appropriate
     // weight functions.
     for (unsigned int qp=0; qp<qrule.n_points(); qp++)
       {
        // Values to hold the solution & its gradient at the previous timestep.
        Number   T2StarValue         = 0.;
        Number   SpinValue           = 0.;
        Number   RealCoilSensitivity = 0.;
        Number   ImagCoilSensitivity = 0.;
        
        // Compute the velocity & its gradient from the previous timestep
        // and the old Newton iterate.
        for (unsigned int l=0; l<n_dofs; l++)
          {
            // From the previous Newton iterate:
            SpinValue           += phi[l][qp]*               SpinSystem.current_solution (dof_indices[l]); 
            T2StarValue         += phi[l][qp]*             T2StarSystem.current_solution (dof_indices[l]); 
            RealCoilSensitivity += phi[l][qp]*RealCoilSensitivitySystem.current_solution (dof_indices[l]); 
            ImagCoilSensitivity += phi[l][qp]*ImagCoilSensitivitySystem.current_solution (dof_indices[l]); 
          }

        for( int timeID = 0 ; timeID <  TimeArray.size(); timeID++)
          {
           Number   PhaseValue  = PhaseModel.GetPhase(q_point[qp],TimeArray[timeID]); 
           // only compute exponential once
           Number JxWphiDecay = JxW[qp] * SpinValue * 
                            std::exp(-TimeArray[timeID]/T2StarValue);
  
           // Compute the signal
           RealSignal[timeID]+=JxWphiDecay * (RealCoilSensitivity * std::cos(PhaseValue) 
                                         - ImagCoilSensitivity * std::sin(PhaseValue));
           ImagSignal[timeID]+=JxWphiDecay * (RealCoilSensitivity * std::sin(PhaseValue) 
                                         + ImagCoilSensitivity * std::cos(PhaseValue));

          } // end of time loop
     } // end of the quadrature point qp-loop

     // At this point the interior element integration has
     // been completed.  However, we have not yet addressed
     // boundary conditions.  
     
    } // end of element loop

  // Add up the signal on all processors and return the values
  Parallel::sum(RealSignal);
  Parallel::sum(ImagSignal);

  // That's it.
  PetscFunctionReturnVoid(); 
}

/** perform matvec product and compute the real and imaginary 
  * part of the signal at a given time instances
  */ 
void SignalModelInterface::AssembleSignal(std::vector<double> *TimeArray,
                                          std::vector<double> *RealSignal,
                                          std::vector<double> *ImagSignal)
{
  PetscFunctionBegin; 
  //error check
  if(!m_PhaseModel)
   {
      std::cerr << "must call SetupGradientModel" << std::endl;
      libmesh_error();
   }

  // dynamic cast and dereference to setup functionoids 
  switch(m_InputGradientModel)
   {
    case LINEARGRADIENT: 
      SignalAssembly(*m_eqn_systems,
                     *m_PhaseModel,
                     *TimeArray, *RealSignal, *ImagSignal);
      break;
    case SPIRALGRADIENT: 
      SignalAssembly(*m_eqn_systems,
                     *dynamic_cast<SpiralGradientModel*>(m_PhaseModel),
                     *TimeArray, *RealSignal, *ImagSignal);
      break;
    case EPIGRADIENT:
      SignalAssembly(*m_eqn_systems,
                     *dynamic_cast<EPIGradientModel*>(m_PhaseModel),
                     *TimeArray, *RealSignal, *ImagSignal);
      break;
    default: 
      std::cerr << "unknown Gradient Model" << std::endl;
      libmesh_error();
   }
  // That's it.
  PetscFunctionReturnVoid(); 
}



LinearGradientModel::LinearGradientModel(Point &gradients,
                     PetscScalar Larmor):m_linearGradient(gradients)
{
 m_larmorFrequency = Larmor;
}

PetscScalar LinearGradientModel::GetPhase(const Point& position,
                                          PetscScalar Time)
{
  return (position(0) * m_linearGradient(0) +
          position(1) * m_linearGradient(1) +
          position(2) * m_linearGradient(2) )* Time;
}

SpiralGradientModel::SpiralGradientModel(Point &gradients,
                     PetscScalar Larmor,PetscScalar angularVel):
                     LinearGradientModel(gradients,Larmor)
{
 m_angularVelocity = angularVel;
}
PetscScalar SpiralGradientModel::GetPhase(const Point& position,
                                          PetscScalar Time)
{
  return (position(0) * m_linearGradient(0) +
          position(1) * m_linearGradient(1) +
          position(2) * m_linearGradient(2) )* Time * m_angularVelocity;
}

EPIGradientModel::EPIGradientModel(Point &gradients,
                     PetscScalar Larmor,int nKspaceSkip):
                     LinearGradientModel(gradients,Larmor)
{
 m_kSpaceSkip = nKspaceSkip;
}
PetscScalar EPIGradientModel::GetPhase(const Point& position,
                                          PetscScalar Time)
{
  return (position(0) * m_linearGradient(0) +
          position(1) * m_linearGradient(1) +
          position(2) * m_linearGradient(2) )* Time * m_kSpaceSkip;
}

