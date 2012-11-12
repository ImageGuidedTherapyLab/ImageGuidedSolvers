// Copyright 2008 Google Inc.
// All Rights Reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following disclaimer
// in the documentation and/or other materials provided with the
// distribution.
//     * Neither the name of Google Inc. nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Author: wan@google.com (Zhanyong Wan)

// This sample shows how to test common properties of multiple
// implementations of the same interface (aka interface tests).

// The interface and its implementations are in this header.
//libmesh
#include "equation_systems.h"
#include "fem_system.h"
#include "boundary_info.h"
#include "fe_base.h"
#include "mesh_base.h"
#include "mesh_generation.h"
#include "fe_interface.h"
#include "quadrature.h"
#include "petsc_macro.h" // define petsc version PETSC_VERSION_LESS_THAN
#include "dof_map.h"
#include "sparse_matrix.h"
#include "mesh.h"
#include "getpot.h"
#include "exact_solution.h"
#include "steady_solver.h"
#include "petsc_diff_solver.h"
//local
#include "tttkUtilities.h"
#include "pennesModel.h" // Constitutive data
#include "thermal_therapy_system.h"
#include "petsc_fem_context.h" 
#include "thermal_therapy_system.txx"
#include "pennesSystem.h"
#include "pennesSystem.txx"
#include "BackgroundPhase.h"
#include "PennesRegression.h"
#include "gtest/gtest.h"

// template instantiations
template class LITTSystem<VerifyPennesConstantSourceTerm>;

//// Manufactured Solution w/ quadratic exact solution
//TEST(ManufacturedSolution, Verify1) {
//
//  PetscFunctionBegin; 
//
//  FiniteElementInterface* user=FiniteElementInterface::getInstance();
//  ASSERT_TRUE(user);
//
//  //user->SetupStructuredGrid();
//
//  //// compute error
//  //libMesh::ExactSolution exactState( eqnSystem );
//  //exactState.attach_exact_value( 
//  //           pdeExactSolution< LITTSystem < MathematicalModel > >  );
//  //exactState.compute_error("StateSystem","u0");
//  //// get forward problem error
//  PetscScalar solutionError =  
//                0.89;
//  //            exactState.l2_error("StateSystem", "u0");
//
//  // compare the norm
//  EXPECT_NEAR(0.0,solutionError,0.9);
//
//  PetscFunctionReturnVoid(); 
//}


// As a general rule, to prevent a test from affecting the tests that come
// after it, you should create and destroy the tested objects for each test
// instead of reusing them.  We will instantiate objects in test's
// SetUp() method and delete them in TearDown() method.
// Then we define a test fixture class template.
template <typename T>
class FiniteDifferenceJacobianTest : public testing::Test {
 protected:
  FiniteDifferenceJacobianTest(){}

  //virtual void SetUp() { }
  //virtual void TearDown() { }

  //virtual ~FiniteDifferenceJacobianTest() { delete table_; }

  typedef T TestSystemType;
  // Note that we test an implementation via the base interface
  // instead of the actual implementation class.  This is important
  // for keeping the tests close to the real world scenario, where the
  // implementation is invoked via the base interface.  It avoids
  // got-yas where the implementation class has a method that shadows
  // a method with the same name (but slightly different argument
  // types) in the base interface, for example.
  PetscFEMSystem *m_pdeSolver; 
};

#if GTEST_HAS_TYPED_TEST

// Google Test offers two ways for reusing tests for different types.
// The first is called "typed tests".  You should use it if you
// already know *all* the types you are gonna exercise when you write
// the tests.

// To write a typed test case, first use
//
//   TYPED_TEST_CASE(TestCaseName, TypeList);
//
// to declare it and specify the type parameters.  As with TEST_F,
// TestCaseName must match the test fixture name.

// The list of types we want to test.
typedef testing::Types< 
               LITTSystem < PennesStandardDiffusionApproximation >,
               RHTESystem < PennesDeltaP1 >,
               //RFASystem  < PennesVoltage > ,
    BackgroundPhaseSystem < BackgroundPhaseModel >
             > Implementations;

TYPED_TEST_CASE(FiniteDifferenceJacobianTest, Implementations);

// Then use TYPED_TEST(TestCaseName, TestName) to define a typed test,
// similar to TEST_F.
TYPED_TEST(FiniteDifferenceJacobianTest, CompareAnalyticFiniteDifferenceMatrixNorm) {

  // Parse the input file
  GetPot controlfile;

  // Read in parameters from the input file
  const Real deltat                    = controlfile("deltat", 0.005);
  const unsigned int coarsegridsize    = controlfile("coarsegridsize", 1);
  const unsigned int dim               = 3;
  
  controlfile.set("power/nsize", 1 );
  controlfile.set("power/power[0]", 0.1234 );

  // Create a mesh.
  Mesh mesh(dim);
  
  // generate mesh w/ BC
  std::vector<int> boundaryData(6,3);
  PetscErrorCode ierr = GenerateStructuredGrid(  mesh,
                                                 coarsegridsize,
                                                 coarsegridsize,
                                                 coarsegridsize,
                                                 0., 1.,
                                                 0., 1.,
                                                 0., 1.,
                                                 boundaryData) ;

  // Print information about the mesh to the screen.
  mesh.print_info();

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // Declare the system and its variables.
  // To refer to typedefs in the fixture, add the 'typename TestFixture::'
  // prefix.  The 'typename' is required to satisfy the compiler.
  equation_systems.parameters.set<GetPot*>("controlfile") = &controlfile;

  this->m_pdeSolver=static_cast< PetscFEMSystem* > ( 
          &equation_systems.add_system< typename TestFixture::TestSystemType > ("StateSystem")
                                                 ) ; 

  // Solve this as a time-dependent or steady system
  this->m_pdeSolver->time_solver =
    AutoPtr<TimeSolver>(new SteadySolver(*this->m_pdeSolver));

  // set nonlinear solver
  this->m_pdeSolver->time_solver->diff_solver() =
        AutoPtr<DiffSolver>(new PetscDiffSolver(*this->m_pdeSolver));

  // Initialize the system
  equation_systems.init ();

  // Set the time stepping options
  this->m_pdeSolver->deltat = deltat;

  // And the nonlinear solver options
  DiffSolver &solver = *(this->m_pdeSolver->time_solver->diff_solver().get());
  solver.quiet = controlfile("solver_quiet", false);

  // Print information about the system to the screen.
  equation_systems.print_info();

  // Compare FD to Jacobian
  // set database options to call SNESSolve_Test
  /* PetscOptionsClear(); */
  ierr = PetscOptionsSetValue("-snes_type","test");
  try 
   {
  this->m_pdeSolver->FEMSystem::solve();
   }
  catch(...) 
   {
   std::cout << "hello" << std::endl; 
   }

  // set to FD Compare Jacobian

  //EXPECT_NEAR(0.0,solutionError,0.9);

  // All done.  
  return ;
}

#endif  // GTEST_HAS_TYPED_TEST

/** @file PennesRegression.cxx
 * @section Analytical Solutions for Pennes with Constant Souce for Use in Regression Testing
 *
 * Constant source with 0 neumann BC
 * Classical  Solutions for the Heat Equation may be may be obtained from a 
 * fundamental solutions approach, See Evans PDE book 1998, Section 2.3 Heat Equation.
 * @htmlonly
 *   see pdf 
 * @endhtmlonly
 * 
 * @latexonly
 * A 1D solution of
 *   \[
 *     u_t - k u_{xx} = q_0
 *   \]
 * is easily verified
 *   \[
 *     u(x,t) = q_0 t \qquad u_{xx} = 0 \qquad u_t = q_0
 *   \]
 * Neumann BC appear the easiest to implement
 *   \[
 *     u_x(0,t) = u_x(L,t) = 0  
 *   \]
 * Dirichlet are also possible but require more effort on the coding side
 *   \[
 *     u(0,t) = u(L,t) =  q_0 t  
 *   \]
 * @endlatexonly
 */
Number PennesConstantSourceExactSolution(const Point& ,
                                         const Parameters& parameters,
                                         const std::string& ,
                                         const std::string& )
{ // return the exact solution
  PetscScalar time = parameters.get<PetscScalar>("time")  ;
  PetscScalar ConstPowerSourceValue = parameters.get<PetscScalar>("ConstPowerSourceValue")  ;
  return ConstPowerSourceValue * time;
}  	 	 
TEST(PennesConstantSourceTerm, NeumannBC) {

  // Parse the input file
  GetPot controlfile;

  // Read in parameters from the input file
  const Real deltat                    = controlfile("deltat", 1.0);
  const unsigned int coarsegridsize    = controlfile("coarsegridsize", 1);
  const unsigned int dim               = 3;
  const unsigned int nstep             = 5;

  // initial temp should be zero 
  controlfile.set("initial_condition/u_init", 0.0 );
  controlfile.set("material/rho", 1.0 );
  controlfile.set("material/specific_heat", 1.0 );

  // set power data
  controlfile.set("power/nsize", nstep );
  PetscScalar ConstPowerSourceValue = 4.83;
  for (unsigned int timeID = 0 ; timeID < nstep ; timeID++)
     {
       std::ostringstream variableName; 
       variableName << "power/power["<< timeID << "]" ;
       controlfile.set(variableName.str(), ConstPowerSourceValue  );
     }

  // Create a mesh.
  Mesh mesh(dim);
  
  // generate mesh w/ BC
  std::vector<int> boundaryData(6,2);
  PetscErrorCode ierr = GenerateStructuredGrid(  mesh,
                                                 10.0*coarsegridsize,
                                                      coarsegridsize,
                                                      coarsegridsize,
                                                 0., 10.,
                                                 0., 1.,
                                                 0., 1.,
                                                 boundaryData) ;

  // Print information about the mesh to the screen.
  mesh.print_info();

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // Declare the system and its variables.
  // To refer to typedefs in the fixture, add the 'typename TestFixture::'
  // prefix.  The 'typename' is required to satisfy the compiler.
  equation_systems.parameters.set<GetPot*>("controlfile") = &controlfile;
  equation_systems.parameters.set<PetscScalar>("ConstPowerSourceValue")  = ConstPowerSourceValue;

  LITTSystem<VerifyPennesConstantSourceTerm>  &state_system = 
      equation_systems.add_system< LITTSystem<VerifyPennesConstantSourceTerm> > ("StateSystem") ; 

  // Solve this as a time-dependent or steady system
  state_system.time_solver =
    AutoPtr<TimeSolver>(new SteadySolver(state_system));

  // set nonlinear solver
  state_system.time_solver->diff_solver() =
        AutoPtr<DiffSolver>(new PetscDiffSolver(state_system));

  // Initialize the system
  equation_systems.init ();

  // Set the time stepping options
  state_system.deltat = deltat;

  // setup the initial conditions
  state_system.SetupInitialConditions();

  // And the nonlinear solver options
  DiffSolver &solver = *(state_system.time_solver->diff_solver().get());
  solver.quiet = controlfile("solver_quiet", false);

  // Print information about the system to the screen.
  equation_systems.print_info();
  equation_systems.parameters.print(std::cout);
  std::cout << std::flush ; 

  // write IC
  ExodusII_IO outputfile(mesh) ;
  std::string visfileName("femdataConstant.e");
  outputfile.write_timestep(visfileName,equation_systems,1,0.0);

  // Compare FD to Jacobian
  // set database options to call default SNESSSolve
  ierr = PetscOptionsClearValue("-snes_type");
  ierr = PetscOptionsSetValue("-ksp_monitor","stdout");
  for (unsigned int timeID = 1 ; timeID < nstep ; timeID++)
   {
     state_system.SetPowerID(timeID);
     state_system.get_vector("old_local_solution") = *state_system.current_local_solution;
     state_system.solve();
     outputfile.write_timestep(visfileName,equation_systems,timeID+1,timeID*deltat);

     // set to FD Compare Jacobian
     equation_systems.parameters.set<PetscScalar>("time")  = timeID * deltat ;
     libMesh::ExactSolution state_exact(equation_systems);
     state_exact.attach_exact_value(PennesConstantSourceExactSolution) ;
     state_exact.compute_error("StateSystem","u0");
     // get forward problem error
     PetscScalar stateerror = state_exact.l2_error("StateSystem", "u0");
     EXPECT_NEAR(0.0,stateerror,0.01);
   }

  // All done.  
  return ;
}
/** @file PennesRegression.cxx
 * @section Analytical Solutions for Delta P1 w/ Exp Source for Use in Regression Testing
 *
 * Using 1D test for regression of coupled delta-p1 / pennes system
 * Classical  Solutions for the Heat Equation may be may be obtained from a 
 * fundamental solutions approach, See Evans PDE book 1998, Section 2.3 Heat Equation.
 * @htmlonly
 *   see pdf 
 * @endhtmlonly
 * 
 * @latexonly
 * A 1D solution of 
 *   \[
 *     u_t - k u_{xx} = \mu_a \exp\left( -\mu_a x \right)
 *   \]
 * is available and is easily verified
 *   \[
 *     u(x,t) = \frac{\exp\left( -\mu_a x \right)}{\mu_a k} \left( \exp\left( k\mu_a^2 t \right) - 1 \right)
 *     \qquad
 *     u_t = \mu_a\exp\left( -\mu_a x \right) \exp\left( k\mu_a^2 t \right) 
 *     \qquad
 *     u_{xx} = \mu_a\frac{\exp\left( -\mu_a x \right)}{ k} \left( \exp\left( k\mu_a^2 t \right) - 1 \right)
 *   \]
 * Unfortuneatly, Neumann BC depends on time but are implementable (always a
 * pleasure with the minus signs) 
 *   \[
 *   - k \cdot u_x(0,t) = \exp\left( -\mu_a 0 \right) \left( \exp\left( k\mu_a^2 t \right) - 1 \right)
 *     \qquad
 *   - k \cdot u_x(L,t) = \exp\left( -\mu_a L \right) \left( \exp\left( k\mu_a^2 t \right) - 1 \right)
 *   \]
 * @endlatexonly
 */
Number PennesExponentialSourceExactSolution(const Point&p ,
                                         const Parameters& parameters,
                                         const std::string& ,
                                         const std::string& varName)
{ // return the exact solution
  PetscScalar time = parameters.get<PetscScalar>("time")  ;
  PetscScalar mu_a_healthy = parameters.get<PetscScalar>("mu_a_healthy") ; 
  PetscScalar k_0_healthy  = parameters.get<PetscScalar>("k_0_healthy" ) ; 

  PetscScalar exactvalue = 0.0;
  if (varName.find("u0")!=std::string::npos ||
      varName.find("u0*")!=std::string::npos )
    {
     exactvalue = std::exp( - mu_a_healthy * p(0) )/ mu_a_healthy/ k_0_healthy  
          * ( std::exp( k_0_healthy * mu_a_healthy * mu_a_healthy * time ) - 1.0 ); 
    }
  else if (varName.find("u1")!=std::string::npos ||
           varName.find("u2")!=std::string::npos )
    {
     exactvalue = std::exp( - mu_a_healthy * p(0) ) ;  // input irradiance
     exactvalue = 0.0; // scattered fluence
    }
  else if (varName.find("d0")!=std::string::npos ||
           varName.find("fx")!=std::string::npos ||
           varName.find("fy")!=std::string::npos ||
           varName.find("fz")!=std::string::npos ||
           varName.find("d1")!=std::string::npos )
    {
     exactvalue = 0.0; 
    }
  else 
    {
     std::cerr << "unknown variable"  << std::endl; libmesh_error();
    }
          
  return exactvalue ;
}  	 	 
TEST(PennesExponentialSourceTerm, NeumannBC) {

  // Parse the input file
  GetPot controlfile;

  // Read in parameters from the input file
  const Real deltat                    = controlfile("deltat", 0.05);
  const unsigned int coarsegridsize    = controlfile("coarsegridsize", 1);
  const unsigned int dim               = 3;
  const unsigned int nstep             = 5;

  // initial temp should be zero 
  controlfile.set("initial_condition/u_init", 0.0 );
  controlfile.set("material/rho", 1.0 );
  controlfile.set("material/specific_heat", 1.0 );
  controlfile.set("optical/mu_s_healthy", 0.0 );
  controlfile.set("optical/mu_a_healthy", 4.1 );
  controlfile.set("thermal_conductivity/k_0_healthy", 0.57 );
  controlfile.set("probe/x_orientation",1.0) ;
  controlfile.set("probe/z_orientation",0.0) ;

  // set power data
  controlfile.set("power/nsize", nstep );
  PetscScalar gaussRad = 1.0e4;
  controlfile.set("optical/guass_radius", gaussRad );  
  PetscScalar ConstPowerSourceValue = libMesh::pi*gaussRad *gaussRad /2.0; 
  for (unsigned int timeID = 0 ; timeID < nstep ; timeID++)
     {
       std::ostringstream variableName; 
       variableName << "power/power["<< timeID << "]" ;
       controlfile.set(variableName.str(), ConstPowerSourceValue );
     }

  // Create a mesh.
  Mesh mesh(dim);
  
  // generate mesh w/ BC
  std::vector<int> boundaryData(6,1);
  boundaryData[2] = 2;
  boundaryData[4] = 2;
  PetscScalar xmin = 0.0;
  PetscScalar xmax = 1.0; 
  PetscErrorCode ierr = GenerateStructuredGrid(  mesh,
                                                 50*coarsegridsize,
                                                    coarsegridsize,
                                                    coarsegridsize,
                                                 xmin, xmax,
                                                 -.05, .05,
                                                 -.05, .05,
                                                 boundaryData) ;

  // Print information about the mesh to the screen.
  mesh.print_info();

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // store some infor for bc evaluation
  equation_systems.parameters.set<PetscScalar>("xmin") = xmin;
  equation_systems.parameters.set<PetscScalar>("xmax") = xmax;

  // Declare the system and its variables.
  // To refer to typedefs in the fixture, add the 'typename TestFixture::'
  // prefix.  The 'typename' is required to satisfy the compiler.
  equation_systems.parameters.set<GetPot*>("controlfile") = &controlfile;

  RHTESystem<VerifyPennesExponentialSourceTerm> &state_system = 
      equation_systems.add_system< RHTESystem<VerifyPennesExponentialSourceTerm> > ("StateSystem") ; 

  // Solve this as a time-dependent or steady system
  state_system.time_solver =
    AutoPtr<TimeSolver>(new SteadySolver(state_system));

  // set nonlinear solver
  state_system.time_solver->diff_solver() =
        AutoPtr<DiffSolver>(new PetscDiffSolver(state_system));

  ExplicitSystem &exactState_system = 
      equation_systems.add_system< ExplicitSystem > ("ExactStateSystem") ; 
  exactState_system.add_variable ("u0*", FIRST);

  // Initialize the system
  equation_systems.init ();

  // Set the time stepping options
  state_system.deltat = deltat;

  // setup the initial conditions
  state_system.SetupInitialConditions();

  // And the nonlinear solver options
  DiffSolver &solver = *(state_system.time_solver->diff_solver().get());
  solver.quiet = controlfile("solver_quiet", false);

  // Print information about the system to the screen.
  equation_systems.print_info();
  equation_systems.parameters.print(std::cout);
  std::cout << std::flush ; 

  equation_systems.parameters.set<PetscScalar>("time")  = 0.0; 
  exactState_system.project_solution(PennesExponentialSourceExactSolution,NULL,
                                     equation_systems.parameters);
  // write IC
  ExodusII_IO outputfile(mesh) ;
  std::string visfileName("femdataExponential.e");
  outputfile.write_timestep(visfileName,equation_systems,1,0.0);

  // Compare FD to Jacobian
  // set database options to call default SNESSSolve
  ierr = PetscOptionsClearValue("-snes_type");
  ierr = PetscOptionsSetValue("-ksp_monitor","stdout");
  ierr = PetscOptionsSetValue("-ksp_converged_reason","TRUE");
  ierr = PetscOptionsSetValue("-ksp_rtol","1.e-9");
  for (unsigned int timeID = 1 ; timeID < nstep ; timeID++)
   {
     state_system.SetPowerID(timeID);
     state_system.get_vector("old_local_solution") = *state_system.current_local_solution;
     // bc at half time step
     equation_systems.parameters.set<PetscScalar>("timePrev")  = (timeID -1) * deltat ;
     equation_systems.parameters.set<PetscScalar>("time")  = timeID * deltat ;

     // check that can recover zero residual for a given exact solution
     state_system.project_solution(PennesExponentialSourceExactSolution,NULL,
                                        equation_systems.parameters);
     state_system.FEMSystem::assembly(true,false);
     state_system.rhs->close();
     EXPECT_NEAR(0.0, state_system.rhs->l2_norm(),2.e-4);
     
     // restore to previous time step and solve
     *state_system.current_local_solution = state_system.get_vector("old_local_solution");
     state_system.FEMSystem::solve();

     // plot at full time step 
     exactState_system.project_solution(PennesExponentialSourceExactSolution,NULL,
                                        equation_systems.parameters);
     outputfile.write_timestep(visfileName,equation_systems,timeID+1,timeID*deltat);

     // set to FD Compare Jacobian
     libMesh::ExactSolution state_exact(equation_systems);
     state_exact.attach_exact_value(PennesExponentialSourceExactSolution) ;
     state_exact.compute_error("StateSystem","u0");
     state_exact.compute_error("StateSystem","u1");
     // get forward problem error
     PetscScalar stateerror[2] = { state_exact.l2_error("StateSystem", "u0"),
                                   state_exact.l2_error("StateSystem", "u1") };
     EXPECT_NEAR(0.0,stateerror[0],0.002);
     EXPECT_NEAR(0.0,stateerror[1],0.001);
   }

  // All done.  
  return ;
}

/** 
 * @section Analytical Solutions for Delta P1 w/ Planar Source for Use in Regression Testing
 *
 * @htmlonly
 *   see pdf 
 * @endhtmlonly
 * 
 * @latexonly
 * A 1D solution for the diffus fluence is 
 *   \[
 *     \varphi_d = E_0 (1-R_s) \left[ \alpha \exp(-\mu_t^* z)  + \beta \exp (-\mu_eff z) \right] 
 *   \]
 * @endlatexonly
 */
Number PennesPlanarSourceExactSolution(const Point&p ,
                                         const Parameters& parameters,
                                         const std::string& ,
                                         const std::string& varName)
{ // return the exact solution
  PetscScalar alpha  = parameters.get<PetscScalar>("alpha" ) ; 
  PetscScalar beta   = parameters.get<PetscScalar>("beta"  ) ; 
  PetscScalar mu_t_star  = parameters.get<PetscScalar>("mu_t_star"  ) ; 
  PetscScalar mu_eff     = parameters.get<PetscScalar>("mu_eff"  ) ; 

  PetscScalar exactvalue = 0.0;
  if (varName.find("u1")!=std::string::npos ||
      varName.find("u1*")!=std::string::npos )
    {
     exactvalue = alpha * std::exp( - mu_t_star * p(0) ) + 
                  beta  * std::exp( - mu_eff    * p(0) ) ;
    }
          
  return exactvalue ;
}  	 	 
TEST(PennesPlanarSourceTerm, CauchyNeumannBC) {

  // Parse the input file
  GetPot controlfile;

  // Read in parameters from the input file
  const Real deltat                    = controlfile("deltat", 0.05);
  const unsigned int coarsegridsize    = controlfile("coarsegridsize", 1);
  const unsigned int dim               = 3;
  const unsigned int nstep             = 2;

  // initial temp should be zero 
  controlfile.set("initial_condition/u_init", 0.0 );
  controlfile.set("material/rho", 1.0 );
  controlfile.set("material/specific_heat", 1.0 );
  controlfile.set("optical/mu_s_healthy", 214.0 );
  controlfile.set("optical/mu_a_healthy", 4.1 );
  controlfile.set("thermal_conductivity/k_0_healthy", 0.57 );
  controlfile.set("probe/x_orientation",1.0) ;
  controlfile.set("probe/z_orientation",0.0) ;

  // set power data
  controlfile.set("power/nsize", nstep );
  PetscScalar ConstPowerSourceValue = 1.0;
  for (unsigned int timeID = 0 ; timeID < nstep ; timeID++)
     {
       std::ostringstream variableName; 
       variableName << "power/power["<< timeID << "]" ;
       controlfile.set(variableName.str(), ConstPowerSourceValue );
     }

  // Create a mesh.
  Mesh mesh(dim);
  
  // generate mesh w/ BC
  std::vector<int> boundaryData(6,0);
  boundaryData[2] = 2;
  boundaryData[4] = 4;
  PetscScalar xmin = 0.0;
  PetscScalar xmax = 0.08; 
  PetscErrorCode ierr = GenerateStructuredGrid(  mesh,
                                                 100*coarsegridsize,
                                                    coarsegridsize,
                                                    coarsegridsize,
                                                 xmin, xmax,
                                                 -.01, .01,
                                                 -.01, .01,
                                                 boundaryData) ;

  // Print information about the mesh to the screen.
  mesh.print_info();

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // store some infor for bc evaluation
  equation_systems.parameters.set<PetscScalar>("xmin") = xmin;
  equation_systems.parameters.set<PetscScalar>("xmax") = xmax;

  // Declare the system and its variables.
  // To refer to typedefs in the fixture, add the 'typename TestFixture::'
  // prefix.  The 'typename' is required to satisfy the compiler.
  equation_systems.parameters.set<GetPot*>("controlfile") = &controlfile;

  RHTESystem<VerifyPennesPlanarSourceTerm> &state_system = 
      equation_systems.add_system< RHTESystem<VerifyPennesPlanarSourceTerm> > ("StateSystem") ; 

  // Solve this as a time-dependent or steady system
  state_system.time_solver =
    AutoPtr<TimeSolver>(new SteadySolver(state_system));

  // set nonlinear solver
  state_system.time_solver->diff_solver() =
        AutoPtr<DiffSolver>(new PetscDiffSolver(state_system));

  ExplicitSystem &exactState_system = 
      equation_systems.add_system< ExplicitSystem > ("ExactStateSystem") ; 
  exactState_system.add_variable ("u1*", FIRST);

  // Initialize the system
  equation_systems.init ();

  // Set the time stepping options
  state_system.deltat = deltat;

  // setup the initial conditions
  state_system.SetupInitialConditions();

  // And the nonlinear solver options
  DiffSolver &solver = *(state_system.time_solver->diff_solver().get());
  solver.quiet = controlfile("solver_quiet", false);

  // Print information about the system to the screen.
  equation_systems.print_info();
  equation_systems.parameters.print(std::cout);
  std::cout << std::flush ; 

  equation_systems.parameters.set<PetscScalar>("time")  = 0.0; 
  exactState_system.project_solution(PennesPlanarSourceExactSolution,NULL,
                                     equation_systems.parameters);
  // write IC
  ExodusII_IO outputfile(mesh) ;
  std::string visfileName("femdataPlanar.e");
  outputfile.write_timestep(visfileName,equation_systems,1,0.0);

  // Compare FD to Jacobian
  // set database options to call default SNESSSolve
  ierr = PetscOptionsClearValue("-snes_type");
  ierr = PetscOptionsSetValue("-ksp_monitor","stdout");
  ierr = PetscOptionsSetValue("-ksp_converged_reason","TRUE");
  ierr = PetscOptionsSetValue("-ksp_rtol","1.e-9");
  unsigned int timeID = 1 ;
   
     state_system.SetPowerID(timeID);
     state_system.get_vector("old_local_solution") = *state_system.current_local_solution;
     // bc at half time step
     equation_systems.parameters.set<PetscScalar>("timePrev")  = (timeID -1) * deltat ;
     equation_systems.parameters.set<PetscScalar>("time")  = timeID * deltat ;

     // check that can recover zero residual for a given exact solution
     state_system.project_solution(PennesPlanarSourceExactSolution,NULL,
                                        equation_systems.parameters);
     state_system.PetscFEMSystem::assembly(true,false);
     state_system.rhs->close();
     EXPECT_NEAR(0.0, state_system.rhs->l2_norm(),2.e-4);

     // plot the residual
     *state_system.solution =  *state_system.rhs;
     outputfile.write_timestep(visfileName,equation_systems,timeID+1,timeID*deltat);

     // restore to previous time step 
     *state_system.current_local_solution = state_system.get_vector("old_local_solution");

     // solve
     state_system.FEMSystem::solve();

     // plot at full time step 
     exactState_system.project_solution(PennesPlanarSourceExactSolution,NULL,
                                        equation_systems.parameters);
     timeID++; 
     outputfile.write_timestep(visfileName,equation_systems,timeID+1,timeID*deltat);

     // set to FD Compare Jacobian
     libMesh::ExactSolution state_exact(equation_systems);
     state_exact.attach_exact_value(PennesPlanarSourceExactSolution) ;
     state_exact.compute_error("StateSystem","u1");
     // get forward problem error
     PetscScalar stateerror = state_exact.l2_error("StateSystem", "u1") ;
     EXPECT_NEAR(0.0,stateerror,0.001);

  // All done.  
  return ;
}
int main(int argc, char** argv) {
  // This allows the user to override the flag on the command line.
  ::testing::InitGoogleTest(&argc, argv);

  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  return RUN_ALL_TESTS();
}
