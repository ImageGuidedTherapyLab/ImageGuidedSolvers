/** @defgroup TreatmentPlanning Prospective Treatment Planning
 *  @example treatmentPlanning.cxx
 * 
 *  Routines for prospective treatment planning simulation. 
 *   - The @ref RFASystem "RFASystem" and @ref PennesVoltage "PennesVoltage" 
 *     classes are intend for simulating RFA
 *   - The @ref LITTSystem "LITTSystem" and @ref PennesStandardDiffusionApproximation "PennesStandardDiffusionApproximation " 
 *     classes are intend for simulating LITT
 *    
 *  Examples of usage of these classes are 
 *  provided in @ref treatmentPlanning.cxx "treatmentPlanning.cxx" 
 */
// petsc includes
#include "petsc.h"
// libmesh include files
#include "equation_systems.h"
#include "error_vector.h"
#include "exodusII_io.h"
#include "kelly_error_estimator.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "mesh_refinement.h"
#include "uniform_refinement_estimator.h"
#include "getpot.h"
#include "diff_solver.h"
#include "time_solver.h"

// Some (older) compilers do not offer full stream
// functionality, OStringStream works around this.
#include "o_string_stream.h"

//local
#include "applicationContext.h"
#include "optimizationParameter.h"
#include "thermal_therapy_system.h"
#include "Utilities.h"

/**@ingroup TreatmentPlanning 
 *  @file treatmentPlanning.cxx
 *  This is a modification of libMesh Example 18 
 *          http://libmesh.sourceforge.net/ex18.php
 *  and is intended as an example for for prospective modeling of thermal
 *  therapy procedures
 */
int main (int argc, char** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

#ifndef LIBMESH_ENABLE_AMR
  if (libMesh::processor_id() == 0)
    std::cerr << "ERROR: This example requires libMesh to be\n"
              << "compiled with AMR support!"
              << std::endl;
  return 0;
#else

  PetscErrorCode ierr;          PetscTruth  flg;

  // Parse the input file
  std::cout << "getting control file..." << std::endl;
  char controlFileName[PETSC_MAX_PATH_LEN];
  // default to "files/control.ini
  sprintf(controlFileName,"%s","files/control.ini");
  ierr = PetscOptionsGetString(PETSC_NULL,"-controlfile",controlFileName,
                                           PETSC_MAX_PATH_LEN,&flg);
  GetPot controlfile(controlFileName);

  // application context
  AppSolve *user = AppSolve::getInstance();

  // setup basic information
  user->Setup(controlfile); user->setupProfile();

  // setup time info
  PetscScalar maxTime = user->setupTime(user->get_num_MRTI(), 0 , 1 );
  std::cout << "max time ..." << maxTime << std::endl;

  // Read in parameters from the input file
  const Real global_tolerance          = controlfile("global_tolerance", 0.);
  const unsigned int nelem_target      = controlfile("n_elements", 400);
  const unsigned int write_interval    = controlfile("write_interval", 1);
  const unsigned int coarserefinements = controlfile("coarserefinements", 0);
  const unsigned int max_adaptivesteps = controlfile("max_adaptivesteps", 0);

  // Create the mesh.
  const unsigned int dim = 3;     // dimension
  Mesh mesh(dim); 
  std::cout << "reading mesh..." << std::endl;
  mesh.read(controlfile("compexec/meshdata","./mesh.e"));

  // And an object to refine it
  MeshRefinement mesh_refinement(mesh);
  mesh_refinement.coarsen_by_parents() = true;
  mesh_refinement.absolute_global_tolerance() = global_tolerance;
  mesh_refinement.nelem_target() = nelem_target;
  mesh_refinement.refine_fraction() = 0.3;
  mesh_refinement.coarsen_fraction() = 0.3;
  mesh_refinement.coarsen_threshold() = 0.1;
  mesh_refinement.uniformly_refine(coarserefinements);

  // Print information about the mesh to the screen.
  mesh.print_info();
  user->printSelf(std::cout);

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);
  equation_systems.parameters.set<GetPot*>("controlfile") = &controlfile ;

  // Setup the FEMSystem for the forward solve
  FEMSystem* system  = user->SetupFEMSystem(equation_systems);

  // Initialize the system
  equation_systems.init ();

  // And the nonlinear solver options
  DiffSolver &solver = *(system->time_solver->diff_solver().get());
  solver.quiet = controlfile("solver_quiet", true);
  solver.max_nonlinear_iterations =
    controlfile("max_nonlinear_iterations", 15);
  solver.relative_step_tolerance =
    controlfile("relative_step_tolerance", 1.e-3);
  solver.relative_residual_tolerance =
    controlfile("relative_residual_tolerance", 0.0);
  solver.absolute_residual_tolerance =
    controlfile("absolute_residual_tolerance", 0.0);
    
  // And the linear solver options
  solver.max_linear_iterations =
    controlfile("max_linear_iterations", 50000);
  solver.initial_linear_tolerance =
    controlfile("initial_linear_tolerance", 1.e-3);

  // Print information about the system to the screen.
  equation_systems.print_info();

  // initialize file output of main data structures
  ExodusII_IO output_vis(mesh);

  // set main solution filename
  OStringStream file_name;
  file_name << "femvis/fem_pdeMethod_"
            << equation_systems.parameters.get<PetscInt>("pdeMethod") << ".e";

  { // Write IC
    // Build the nodal solution values & get the variable
    // names from the EquationSystems object
    std::vector<Number>      soln; std::vector<std::string> names;
    equation_systems.build_variable_names  (names);
    equation_systems.build_solution_vector (soln);

    // Write out every ideal timestep to file.
    output_vis.write_nodal_data(file_name.str(),soln,names);

    // write timestep info to file 
    output_vis.write_timestep(system->time);
  }

  // Now we begin the timestep loop to compute the time-accurate
  // solution of the equations.
  for (AppSolve::ISTEP=0; AppSolve::ISTEP != user->get_num_MRTI(); 
                                                AppSolve::ISTEP++)
    {
      // A pretty update message
      std::cout << "\n\nSolving time step " << AppSolve::ISTEP<< ", time = "
                << system->time << std::endl;

      // Adaptively solve the timestep
      unsigned int a_step = 0;
      for (; a_step != max_adaptivesteps; ++a_step)
        {
          system->solve();

          ErrorVector error;

          AutoPtr<ErrorEstimator> error_estimator;

          // To solve to a tolerance in this problem we
          // need a better estimator than Kelly
          if (global_tolerance != 0.)
            {
              // We can't adapt to both a tolerance and a mesh
              // size at once
              libmesh_assert (nelem_target == 0);

              UniformRefinementEstimator *u =
                new UniformRefinementEstimator;

              // The lid-driven cavity problem isn't in H1, so
              // lets estimate H0 (i.e. L2) error
              u->sobolev_order() = 0;

              error_estimator.reset(u);
            }
          else
            {
              // If we aren't adapting to a tolerance we need a
              // target mesh size
              libmesh_assert (nelem_target > 0);

              // Kelly is a lousy estimator to use for a problem
              // not in H1 - if we were doing more than a few
              // timesteps we'd need to turn off or limit the
              // maximum level of our adaptivity eventually
              error_estimator.reset(new KellyErrorEstimator);
            }

          // Calculate error based on u and v but not p
          error_estimator->component_scale.push_back(1.0); // u
          error_estimator->component_scale.push_back(1.0); // v
          if (dim == 3)
            error_estimator->component_scale.push_back(1.0); // w
          error_estimator->component_scale.push_back(0.0); // p

          error_estimator->estimate_error(*system, error);

          // Print out status at each adaptive step.
          Real global_error = error.l2_norm();
          std::cout << "adaptive step " << a_step << ": ";
          if (global_tolerance != 0.)
            std::cout << "global_error = " << global_error
                      << " with ";
          std::cout << mesh.n_active_elem()
                    << " active elements and "
                    << equation_systems.n_active_dofs()
                    << " active dofs." << std::endl;
          if (global_tolerance != 0.)
            std::cout << "worst element error = " << error.maximum()
                      << ", mean = " << error.mean() << std::endl;

          if (global_tolerance != 0.)
            {
              // If we've reached our desired tolerance, we
              // don't need any more adaptive steps
              if (global_error < global_tolerance)
                break;
              mesh_refinement.flag_elements_by_error_tolerance(error);
            }
          else
            {
              // If flag_elements_by_nelem_target returns true, this
              // should be our last adaptive step.
              if (mesh_refinement.flag_elements_by_nelem_target(error))
                {
                  mesh_refinement.refine_and_coarsen_elements();
                  equation_systems.reinit();
                  a_step = max_adaptivesteps;
                  break;
                }
            }

          // Carry out the adaptive mesh refinement/coarsening
          mesh_refinement.refine_and_coarsen_elements();
          equation_systems.reinit();
        }
      // Do one last solve if necessary
      if (a_step == max_adaptivesteps)
        {
          system->solve();
        }

      // Advance to the next timestep in a transient problem
      system->time_solver->advance_timestep();

      // Write out this timestep if we're requested to
      if ((AppSolve::ISTEP+1)%write_interval == 0)
        {
          // Build the nodal solution values & get the variable
          // names from the EquationSystems object
          std::vector<Number>      soln; std::vector<std::string> names;
          equation_systems.build_variable_names  (names);
          equation_systems.build_solution_vector (soln);

          // Write out every ideal timestep to file.
          output_vis.write_nodal_data(file_name.str(),soln,names);

          // write timestep info to file 
          output_vis.write_timestep(system->time);
        }
      // increment timestep (1-based) 
      output_vis.increment_timestep();
    }
#endif // #ifndef LIBMESH_ENABLE_AMR
  
  // All done.  
  return 0;
}
