/* $Id: ex4.C 2501 2007-11-20 02:33:29Z benkirk $ */

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



 // <h1>Example 4 - Solving a 1D, 2D or 3D Poisson Problem in Parallel</h1>
 //
 // This is the fourth example program.  It builds on
 // the third example program by showing how to formulate
 // the code in a dimension-independent way.  Very minor
 // changes to the example will allow the problem to be
 // solved in one, two or three dimensions.
 //
 // This example will also introduce the PerfLog class
 // as a way to monitor your code's performance.  We will
 // use it to instrument the matrix assembly code and look
 // for bottlenecks where we should focus optimization efforts.
 //
 // This example also shows how to extend example 3 to run in
 // parallel.  Notice how litte has changed!  The significant
 // differences are marked with "PARALLEL CHANGE".


// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>

// Various include files needed for the mesh & solver functionality.
#include "libmesh.h"
#include "mesh.h"
#include "mesh_function.h"
#include "gmv_io.h"
#include "equation_systems.h"
#include "exact_solution.h"
#include "fe.h"
#include "fe_interface.h"
#include "quadrature_gauss.h"
#include "dof_map.h"
#include "sparse_matrix.h"
#include "elem.h"
#include "string_to_enum.h"
#include "getpot.h"
#include "boundary_info.h"

// The nonlinear solver and system we will be using
#include "nonlinear_implicit_system.h"
#include "linear_implicit_system.h"

// dddas include files
using namespace std;
#include "transient_inverse_system.h"
#include "tao.h" // Petsc/Tao include files
#include "solvehp3d.h" // main header
#include "fortrandftranslate.h" // header to call fortran routines
#include "variable_map.h" // variables using in both f90 and C++

/* --------------Global typedef ---------- */
/* --------------External Fortran Routines---------- */
extern FORTRAN_FUNCTION 
{
#include "global_params.h" // global data
#include "pennes_model.h" // interface to pennes model
}

/* -----------------------The Initial Condition------------------------ */
Number initial_temperature (const Point& p,
                    const Parameters& parameters,
                    const std::string&,
                    const std::string&)
{
  PetscScalar time = parameters.get<Real> ("Time");
  PetscScalar xpoint[3] = { p(0), p(1), p(2) };
  return FORTRAN_NAME(getinittemp)(xpoint,&time);
}

// We now define the function which provides the
// initialization routines for the "StateSystem"
// system.  This handles things like setting initial
// conditions and boundary conditions.
#undef __FUNCT__
#define __FUNCT__ "pennes_init_cd "
void pennes_init_cd (EquationSystems& es,
              const std::string& system_name)
{
  PetscFunctionBegin;
  // It is a good idea to make sure we are initializing
  // the proper system.
  libmesh_assert (system_name == "StateSystem");

  // Get a reference to the Convection-Diffusion system object.
  TransientInverseNonlinearImplicitSystem & system =
    es.get_system<TransientInverseNonlinearImplicitSystem>("StateSystem");

  system.project_solution(initial_temperature, NULL, es.parameters);

  // store the initial condition for all time instances before Nsteplo
  const AppSolve     &user = *_user_app;           // application context
  for(PetscInt iii = 0 ; iii <= user.Nsteplo ; iii++)
      *system.vector_solution.at(iii) = *system.current_local_solution;

  PetscFunctionReturnVoid();
}
// We now define the function which provides the
// initialization routines for the "StateSystem"
// system.  This handles things like setting initial
// conditions and boundary conditions.
#undef __FUNCT__
#define __FUNCT__ "pennes_init_rf "
void pennes_init_rf (EquationSystems& es,
              const std::string& system_name)
{
  PetscFunctionBegin;

  // It is a good idea to make sure we are initializing
  // the proper system.
  libmesh_assert (system_name == "StateSystem");

  // Get a reference to the Convection-Diffusion system object.
  TransientInverseNonlinearImplicitSystem & system =
    es.get_system<TransientInverseNonlinearImplicitSystem>("StateSystem");

  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es.get_mesh();

  // store the initial condition
  const AppSolve     &user = *_user_app;           // application context

  // The number of variables in this system
  const unsigned int n_variables = system.n_vars();


  // The dimensionality of the current mesh
  const unsigned int dim = system.get_mesh().mesh_dimension();

  // The DofMap for this system
  const DofMap& dof_map = system.get_dof_map();

  // The new element coefficients
  DenseVector<Number> Ue;

  // the initial values are are assumed subdomain dependent
  libmesh_assert (n_variables <= 2);  
  // FIXME this could prob be a little cleaner ranther than hard code "2"
  std::vector< std::vector < Real > > InitValues(2, 
                                      std::vector <Real> (2, 0.0 ) );
  PetscScalar       dum[3]={0.0,0.0,0.0};
  PetscScalar time = 0.5*(2*user.Nsteplo-1) * user.GenInfo->FEM_DT;
  // regular pennes/rf coupled equations for the Zeroth Subdomain
  InitValues[0][0] = FORTRAN_NAME(getinittemp)(dum,&time) ;
  InitValues[0][1] = 0.0 ; // set initial voltage to zero
   
  // In the 1st (second) subdomain the voltage is held constant (spatially) and
  // the temperature is assumed to be the temperature of the cooling fluid
  InitValues[1][0] = FORTRAN_NAME(get_probe_temp)() ;
  PetscInt  idpow = (user.Nsteplo * user.GenInfo->FEM_DT) / 
                           user.GenInfo->ISTEPS_PER_IDEAL  + 1 ;
  FORTRAN_NAME(getparam_pow)(&InitValues[1][1],&idpow) ;

  // loop over subdomains so that the last subdomain has the final
  // overwrite of the nodal temperature
  for ( unsigned int id_sub= 0 ; id_sub < 2; ++id_sub)
   {    
     // Loop over all the variables in the system
     for (unsigned int var=0; var<n_variables; var++)
       {
         // Get FE objects of the appropriate type
         const FEType& fe_type = dof_map.variable_type(var);     
         AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));      

         // Prepare variables for projection
         AutoPtr<QBase> qrule     (fe_type.default_quadrature_rule(dim));

         // The values of the shape functions at the quadrature
         // points
         const std::vector<std::vector<Real> >& phi = fe->get_phi();

         const FEContinuity cont = fe->get_continuity();

         // only setup for C0 first order elements
         libmesh_assert (cont          == C_ZERO);
         libmesh_assert (fe_type.order == FIRST);

         // The Jacobian * quadrature weight at the quadrature points
         const std::vector<Real>& JxW =
           fe->get_JxW();
        
         // The XYZ locations of the quadrature points
         const std::vector<Point>& xyz_values =
           fe->get_xyz();

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
         MeshBase::const_element_iterator       el    =mesh.active_elements_begin();
         const MeshBase::const_element_iterator end_el=mesh.active_elements_end();
         for ( ; el != end_el; ++el)
            {    
             // Store a pointer to the element we are currently
             // working on.  This allows for nicer syntax later.
             const Elem* elem = *el;
             
             // Store element subdomain id to allow different sets of equations
             //  to be solved on different parts of the mesh
             // ASSUMES SUBDOMAIN ORDERING STARTS FROM ONE but we need a 0 based
             // numbering scheme for std::vector
             const unsigned int subdomain_id = elem->subdomain_id() - 1;

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
             // solution.  We start by interpolating nodes, then
             // hold those fixed and project edges, then
             // hold those fixed and project faces, then
             // hold those fixed and project interiors

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

                 Ue(current_dof) = InitValues[subdomain_id][var];

                 dof_is_fixed[current_dof] = true;
                 current_dof++;
               }

             // Make sure every DoF got reached!
             for (unsigned int i=0; i != n_dofs; ++i)
               libmesh_assert(dof_is_fixed[i]);

             const unsigned int
               first = system.solution->first_local_index(),
               last  = system.solution->last_local_index();

             // put solution in parallel data structures...
             for (unsigned int i = 0; i < n_dofs; i++) 
               // We may be projecting a new zero value onto
               // an old nonzero approximation - RHS
               // if (Ue(i) != 0.)
               if ((dof_indices[i] >= first) &&
                   (dof_indices[i] <  last))
                 {
                   system.solution->set( dof_indices[i], Ue(i));
                 }
           }  // end elem loop
       } // end variables loop
    // close after each subdomain so that the last subdomain will have
    // the final write...
    system.solution->close();
  } // end loop over subdomains    

  // copy parallel data structures to local data structures
  system.solution->localize(*system.current_local_solution);
  *system.vector_solution.at(user.Nsteplo) = *system.current_local_solution;

  PetscFunctionReturnVoid();
}
/* Function to compute stored solutions for plotting and QOI evaluation purposes.  */
Number pennes_stored_solution (const Point& p,
                                const Parameters& parameters,
                                const std::string& sys_name,
                                const std::string& unknown_name)
{
 PetscLogEventBegin(logevents[27],0,0,0,0); // eval mesh fcn
 MeshFunction *stored_values = parameters.get<MeshFunction*> ("MeshFunction");
 Number value =  (*stored_values)(p);
 PetscLogEventEnd(  logevents[27],0,0,0,0); // eval mesh fcn

 return value;
}



