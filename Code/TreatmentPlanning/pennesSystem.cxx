//libmesh
#include "equation_systems.h"
#include "fem_system.h"
#include "boundary_info.h"
#include "fe_base.h"
#include "mesh_base.h"
#include "fe_interface.h"
#include "quadrature.h"
#include "petsc_macro.h" // define petsc version PETSC_VERSION_LESS_THAN
#include "dof_map.h"
#include "sparse_matrix.h"
#include "getpot.h"
#include "exact_solution.h"
#include "steady_solver.h"
#include "petsc_diff_solver.h"
#include "petsc_vector.h"
#include "petsc_matrix.h"
#include "mesh_function.h"
#include "point_locator_base.h"
//local
#include "pennesModel.h" // Constitutive data
#include "thermal_therapy_system.h"
#include "thermal_therapy_system.txx"
#include "pennesSystem.h"
#include "pennesSystem.txx"

// ------------------------------------------------------------

/**
 *  Add SDA System
 */
#undef __FUNCT__
#define __FUNCT__ "AddPennesSDASystem"
System * AddPennesSDASystem( char *SystemName, 
                  EquationSystems *m_eqn_systems, double DeltaT)
{
  PetscFunctionBegin;

  // store
  std::cout <<"Setting up Bioheat with SDA Laser "<<std::endl;
  typedef  LITTSystem < PennesStandardDiffusionApproximation >  PennesSDASystem;
  FEMSystem *m_pdeSolver = static_cast< FEMSystem* >( 
            &m_eqn_systems->add_system< PennesSDASystem > (SystemName)
                                                    ) ; 
  m_pdeSolver->deltat = DeltaT;
  std::cout <<"deltat = "<< m_pdeSolver->deltat <<std::endl;
  PetscFunctionReturn(m_pdeSolver);

}

// update the laser position
#undef __FUNCT__
#define __FUNCT__ "PennesSDASystemUpdateLaserPosition"
PetscErrorCode PennesSDASystemUpdateLaserPosition(System * currentSystem, 
                           PetscScalar X0,PetscScalar Y0, PetscScalar Z0, 
                           PetscScalar X1,PetscScalar Y1, PetscScalar Z1) 
{
  PetscFunctionBegin; 
  PetscErrorCode ierr;

  // get pointer to petsc vector
  typedef LITTSystem<PennesStandardDiffusionApproximation>  PennesSDASystem;
  PennesSDASystem *penesSDASystem = 
                     dynamic_cast<PennesSDASystem *>(currentSystem);
  ierr = penesSDASystem->m_MathModel.UpdateLaserPosition(X0,Y0,Z0,X1,Y1,Z1);

  PetscFunctionReturn(ierr);
}

// update the laser position
#undef __FUNCT__
#define __FUNCT__ "PennesSDASystemUpdateLaserPower"
PetscErrorCode PennesSDASystemUpdateLaserPower(System * currentSystem, 
                           PetscScalar newPower,PetscInt timeID ) 
{
  PetscFunctionBegin; 
  PetscErrorCode ierr;

  // get pointer to petsc vector
  typedef LITTSystem<PennesStandardDiffusionApproximation>  PennesSDASystem;
  PennesSDASystem *penesSDASystem = 
                     dynamic_cast<PennesSDASystem *>(currentSystem);
  ierr = penesSDASystem->m_MathModel.setPower(timeID,newPower);

  PetscFunctionReturn(ierr);
}

/**
 *  Add RF System
 */
#undef __FUNCT__
#define __FUNCT__ "AddPennesRFSystem"
System * AddPennesRFSystem( char *SystemName, 
                  EquationSystems *m_eqn_systems, double DeltaT)
{
  PetscFunctionBegin;

  // store
  std::cout <<"Setting up Bioheat with coupled RF "<<std::endl;
  // get constitutive data from ini file and print to screen 
  FEMSystem *m_pdeSolver = static_cast< FEMSystem* > 
                                         ( &m_eqn_systems->add_system< 
                     RFASystem < PennesVoltage > > (SystemName)
                                         ) ; 
  m_pdeSolver->deltat = DeltaT;
  std::cout <<"deltat = "<< m_pdeSolver->deltat <<std::endl;
  PetscFunctionReturn(m_pdeSolver);

}
/**
 *  Add delta-p System
 */
#undef __FUNCT__
#define __FUNCT__ "AddPennesDeltaPSystem"
System * AddPennesDeltaPSystem( char *SystemName, 
                  EquationSystems *m_eqn_systems, double DeltaT)
{
  PetscFunctionBegin;

  // store
  std::cout <<"Setting up Bioheat with Delta P1"<<std::endl;
  // get constitutive data from ini file and print to screen 
  FEMSystem *m_pdeSolver = static_cast< FEMSystem* > 
                                         ( &m_eqn_systems->add_system< 
          RHTESystem < PennesDeltaP1 > > (SystemName)
                                         ) ; 
  m_pdeSolver->deltat = DeltaT;
  std::cout <<"deltat = "<< m_pdeSolver->deltat <<std::endl;
  PetscFunctionReturn(m_pdeSolver);

}


/* Function to compute stored solutions for plotting and QOI evaluation purposes.  */
Real get_fluence_data(const Point& point,
                      const Parameters& parameters,
                      const std::string& sys_name,
                      const std::string& unknown_name)
{
 //use mesh function to location element
 MeshFunction *solution_values= parameters.get<MeshFunction*> ("MeshFunction");
 const PointLocatorBase& PointLocate = solution_values->get_point_locator();
 const Elem* elem = PointLocate(point);

 //AppSolve *user = parameters.get<AppSolve*> ("AppSolve");
 //PennesSDA* PennesSDASolver= dynamic_cast<PennesSDA*>(user->pdeSolver);
 //const unsigned int subdomain_id = elem->subdomain_id() - 1;
 //return state_system.m_MathModel.PennesSource(subdomain_id,37.0,
 //                                                damage,0.0, 
 //                                                point, 0);
 return 0.0;
}
