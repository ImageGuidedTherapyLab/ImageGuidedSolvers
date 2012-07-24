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
#include "sparse_matrix.h"
#include "elem.h"
#include "string_to_enum.h"
#include "getpot.h"
#include "boundary_info.h"
#include "parallel.h" // mpi utilities
#include "petsc_linear_solver.h" // mpi utilities
#include "exact_solution.h"
#include "time_solver.h"
#include "petsc_matrix.h"
#include "petsc_diff_solver.h"
#include "steady_solver.h"

// pdeopt
#include "applicationContext.h"
#include "pennesInverseModel.h"
#include "pennesInverseSystem.h"
#include "thermal_therapy_system.txx"
#include "qoiBaseClass.h"  // interface to qoi 
#include "mesh.h"
#include "pennesSystem.txx"
#include "pennesInverseSystem.txx"

// ------------------------------------------------------------
/**
 * Wrapper for adjoint computation routine
 */
#undef __FUNCT__
#define __FUNCT__ "assemble_adjoint"
void assemble_adjoint(EquationSystems& es,
                      const std::string& )
{
 // application context
 AppSolve* user    = es.parameters.get<AppSolve*>("AppSolve")  ;
 TransientFEMSystem*  pdeSolver = user->pdeSolver();  // abstract pde solver

 pdeSolver->assemble_adjoint(es); // abstract pde solver

}
#undef __FUNCT__
#define __FUNCT__ "initial_sensitivity_matrix"
void initial_sensitivity_matrix (EquationSystems& es, const std::string& )
{
 PetscFunctionBegin;

 // Get a reference to the NonlinearImplicitSystem we are solving
 TransientFEMSystem& state_system = 
          es.get_system<TransientFEMSystem>("StateSystem");

 TransientFEMSystem& sensitivity_system = 
   es.get_system<TransientFEMSystem>("SensitivitySystem");

 // reuse the state system matrix
 if(state_system.matrix)
   {
    //delete sensitivity_system.matrix;  this should be delete by libmesh destructor
    sensitivity_system.matrix = state_system.matrix;
    sensitivity_system.SetAssembleJacobian(false);
   }
 else
    {std::cerr << "state system uninitialized?" << std::endl << std::flush; libmesh_error();}
 
 PetscFunctionReturnVoid();
}
// ------------------------------------------------------------
// template instantiations
template class PennesInverseSystem<PennesInverseModel>;

