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

 // <h1>Background Correction Example - Modeled after libMesh ex18</h1>
 //
 // This is an attempt to provide the fewest lines possible to interact
 // with classes in main code for background phase correction

#include <vector>
#include "tao.h" // petsc solvers w/ Fortran Interface
#include "libmesh.h" // petsc solvers w/ Fortran Interface
#include "getpot.h" // file parser
#include "equation_systems.h"
#include "nonlinear_implicit_system.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "boundary_info.h"
#include "diff_solver.h"
#include "euler_solver.h"
#include "steady_solver.h"
#include "petsc_diff_solver.h"
#include "petsc_matrix.h"

// local
#include "applicationContext.h"
#include "pdeBaseClass.h"
#include "pennesInverseSystem.h"
#include "BackgroundPhase.h"
#include "Utilities.h"
#include "Imaging.h"
#include "parser_defaults.h"

/* Function to get MRTI data */
Number InterpolateITKImageData (const Point& p,
                                const Parameters& parameters,
                                const std::string& ,
                                const std::string& )
{
  InterpolatorType::Pointer interpolator = parameters.get<
                        InterpolatorType::Pointer>("interpolator") ;

  InterpolatorType::PointType point;
  point[0] = p(0);
  point[1] = p(1);
  point[2] = p(2);

  Number temp = parameters.get<PetscScalar>("DefaultImageValue");
  if( interpolator->IsInsideBuffer(point) )
     temp =  interpolator->Evaluate( point ); 

  return temp;
}

// The main program.
// TODO write python/matlab wrapper to this routine to interface with
// TODO matlab code and python code
int main (int argc, char** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);
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

  // Read in parameters from the input file
  const Real deltat                    = controlfile("deltat", 0.005);
  const unsigned int dim               = 3;
  libmesh_assert ( dim == 3 );

  // read phase image
  typedef itk::ImageFileReader< InputImageType  >  ReaderType;
  // setup data structures
  ReaderType::Pointer reader = ReaderType::New();
  // set filename 
  char phaseFileName[PETSC_MAX_PATH_LEN];
  // default to "files/control.ini
  sprintf(phaseFileName,"%s","phase.vtk");
  ierr = PetscOptionsGetString(PETSC_NULL,"-phasefile",phaseFileName,
                                           PETSC_MAX_PATH_LEN,&flg);
  std::string FileName( phaseFileName );
  reader->SetFileName( FileName );
  // Software Guide : BeginCodeSnippet
  reader->Update();
  InputImageType::Pointer phaseImage =  reader->GetOutput();

  // get image size information
  InputImageType::RegionType::SizeType 
      imageSize=phaseImage->GetRequestedRegion().GetSize();
  const InputImageType::SpacingType& imageSpacing = phaseImage->GetSpacing();
  const InputImageType::PointType& orgnImage = phaseImage->GetOrigin();

  // Create a n-dimensional mesh.
  Mesh mesh (dim);
  
  // Use the MeshTools::Generation mesh generator to create a uniform
  // grid on the square [-1,1]^D.  We instruct the mesh generator
  // to build a mesh of 8x8 \p Quad9 elements in 2D, or \p Hex27
  // elements in 3D.  Building these higher-order elements allows
  // us to use higher-order approximation, as in example 3.
  MeshTools::Generation::build_cube (mesh,
                                     imageSize[0],
                                     imageSize[1],
                                     imageSize[2],
       orgnImage[0] , orgnImage[0] + imageSize[0] * imageSpacing[0],
       orgnImage[1] , orgnImage[1] + imageSize[1] * imageSpacing[1],
       orgnImage[2] , orgnImage[2] + imageSize[2] * imageSpacing[2],
                                     HEX8);
    // setup exodus domain and remove existing BC
    libMesh::MeshBase::const_element_iterator       global_el   = mesh.elements_begin();
    const libMesh::MeshBase::const_element_iterator global_el_end = mesh.elements_end();
    for (  ; global_el != global_el_end ; ++global_el  )
      {
        // Store a pointer to the element we are currently
        // working on.  This allows for nicer syntax later.
        Elem* elem = *global_el;

        // exodus files start w/ domain one
        elem->subdomain_id() = 1;

        // remove any boundary conditions
        mesh.boundary_info->remove(elem);
      }
    // set custom boundary info
    // loops should be same as MeshTools::Generation::build_cube
    for (unsigned int k=0; k<imageSize[2] ; k++)
     for (unsigned int j=0; j<imageSize[1] ; j++)
      for (unsigned int i=0; i<imageSize[0] ; i++)
       {
        unsigned int elemID =  i + j*imageSize[0] + k*imageSize[0]*imageSize[1]; 
        Elem* elem = mesh.elem( elemID );
        if (k == 0)       
          {
          if ( j==0 || j==(imageSize[1]-1) || i==0  || i==(imageSize[0]-1) )
                              mesh.boundary_info->add_side(elem, 0, 2);
          else 
                              mesh.boundary_info->add_side(elem, 0, 3);
          }
        if (k == (imageSize[2]-1)) 
          {
          if ( j==0 || j==(imageSize[1]-1) || i==0  || i==(imageSize[0]-1) )
                              mesh.boundary_info->add_side(elem, 5, 2);
          else 
                              mesh.boundary_info->add_side(elem, 5, 3);
          }
        if (j == 0)           mesh.boundary_info->add_side(elem, 1, 3);
        if (j == (imageSize[1]-1)) mesh.boundary_info->add_side(elem, 3, 3);
        if (i == 0)           mesh.boundary_info->add_side(elem, 4, 3);
        if (i == (imageSize[0]-1)) mesh.boundary_info->add_side(elem, 2, 3);
       }

  // Print information about the mesh to the screen.
  mesh.print_info();
  user->printSelf(std::cout);

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // Declare the system "BackgroundSystem" and its variables.
  // the controlfile will be passed in when models are initialized
  equation_systems.parameters.set<GetPot*>("controlfile") = &controlfile;
  BackgroundPhaseSystem< BackgroundPhaseModel > & background_system =
       equation_systems.add_system< BackgroundPhaseSystem < BackgroundPhaseModel > >("BackgroundSystem");

  // Solve this as a steady system
  background_system.time_solver = AutoPtr<TimeSolver>(new SteadySolver(background_system));

  // Initialize the system
  equation_systems.init ();

  // Set the time stepping options
  background_system.deltat = deltat;

  // And the nonlinear solver options
  DiffSolver &solver = *(background_system.time_solver->diff_solver().get());
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
  // write out input parameters
  equation_systems.parameters.print();

  // A pretty update message
  std::cout << "\n\nSolving..." << std::endl;

  // put MRTI thermal image into FEM data structures
  equation_systems.parameters.set<PetscScalar>("DefaultImageValue") =  0.0;
  
  // initialize interpolator and set pointer
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetInputImage( phaseImage );

  // put MRTI thermal image into FEM data structures
  equation_systems.parameters.set<InterpolatorType::Pointer>(
                                 "interpolator") = interpolator;
  background_system.project_solution(InterpolateITKImageData,NULL,
                                     equation_systems.parameters);

  // solve Neumann Problem must remove NULL Space
  PetscDiffSolver* backgroundSolver = 
                  libmesh_cast_ptr<PetscDiffSolver*>(
                  & (*background_system.time_solver->diff_solver()) );
  PetscErrorCode info;
  KSP  snesksp;
  MatNullSpace NullSpace;
  info = SNESGetKSP(backgroundSolver->snes(),&snesksp);CHKERRQ(info);
  info = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,
                            PETSC_NULL,&NullSpace); CHKERRQ(info);
  info = KSPSetNullSpace(snesksp,NullSpace);CHKERRQ(info);

  // assemble and solve 
  background_system.solve();

  // clean up
  info = MatNullSpaceDestroy(NullSpace);CHKERRQ(info);

  // We write the file name in the gmv auto-read format.
  OStringStream file_name;
  file_name << "phase.e";
  ExodusII_IO(mesh).write_equation_systems (file_name.str(),
                                            equation_systems);
  
  // write matrix
  OStringStream matrixFileName;
  matrixFileName << "stiffnessMatrix.m" ;
  PetscViewer matViewer; 
  info = PetscViewerASCIIOpen(PETSC_COMM_WORLD,matrixFileName.str().c_str(),
                                               &matViewer);  CHKERRQ(info);
  info = PetscViewerSetFormat(matViewer,PETSC_VIEWER_ASCII_MATLAB); 
                                                             CHKERRQ(info);
  PetscMatrix<Number>* sysMatrix = libmesh_cast_ptr<PetscMatrix<Number>*>( 
                                            &(*background_system.matrix) );
  info = MatView(sysMatrix->mat(),matViewer); CHKERRQ(info);
  info = PetscViewerFlush(matViewer); CHKERRQ(info);
  info = PetscViewerDestroy(matViewer); CHKERRQ(info);
  PetscPrintf(PETSC_COMM_WORLD,"wrote stiffness Matrix...\n" );

  // All done.  
  return 0;
}
