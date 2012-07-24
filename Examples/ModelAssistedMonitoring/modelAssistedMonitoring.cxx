/** @defgroup ModelAssistedMonitoring Model Assisted Monitoring 
 *  Routines for Model Assisted Monitoring 
 */
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

// C include files 
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>

// Basic include file needed for the mesh functionality.
#include "libmesh.h"
#include "mesh.h"
#include "mesh_refinement.h"
#include "mesh_generation.h"
#include "gmv_io.h"
#include "equation_systems.h"
#include "fe.h"
#include "fe_interface.h"
#include "quadrature_gauss.h"
#include "dof_map.h"
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "petsc_nonlinear_solver.h"
#include "nonlinear_implicit_system.h"
#include "perf_log.h"
#include "boundary_info.h"
#include "utility.h"
#include "getpot.h" 

// The definition of a geometric element
#include "elem.h"

// local pde optimization includes
#include "applicationContext.h" 
#include "optimizationParameter.h"
#include "thermal_therapy_system.h"
#include "Utilities.h"
#include "quantityOfInterest.h"
#include "dddas.h"

static  char help[] = "Model Assisted Monitoring Driver\n";
/** -------------- For Debugging  ------------ */
volatile int endstall=0;
/** 
 * Give time to attach a debugger
 */
void stallForDebugger(PetscInt rank)
{
     /* get the process id */
     pid_t pid = getpid();
     std::ofstream debugfile;
     std::ostringstream file_name;
     file_name << getenv("HOME") << "/debugrank"<< rank;
     debugfile.open(file_name.str().c_str(), std::ios::out);
     debugfile << "cd ${DDDAS_SRC}; python debug.py --pid="<<pid<<std::endl;
     debugfile.close();
     // make the file executable
     chmod(file_name.str().c_str(),S_IREAD|S_IWRITE|S_IEXEC);
     // echo info
     std::cout << getenv("SSH_CONNECTION") << " rank " << rank 
               << " process id is " << pid << std::endl << std::flush ;
     while(!endstall)(sleep(5));
}
 
/**@ingroup InverseSolver 
 *  @file inverseSolver.cxx
 *  Main Routine for Inverse Solver
 */
int main (int argc, char** argv)
{

  // Initialize libMesh.
  std::cout << "initializing libMesh..." << std::endl;
  LibMeshInit init (argc, argv);

  /* Initialize TAO */
  std::cout << "initializing TAO..." << std::endl;
  TaoInitialize(&argc,&argv,(char *)0,help);

  // Set the dimensionality of the mesh = 3
  const unsigned int dim = 3;     
  PetscErrorCode ierr;          PetscTruth  flg;
  
  /*Check for debugging flags*/
  PetscTruth  debug=PETSC_FALSE;// debug flag
  ierr = PetscOptionsGetTruth(PETSC_NULL,"-idb",&debug,&flg); CHKERRQ(ierr);
  PetscInt  debugid=-1; // default to all
  ierr = PetscOptionsGetInt(PETSC_NULL,"-idbrank",&debugid,&flg); CHKERRQ(ierr);
  if(debug){ // debug
      PetscInt    rank; // world communicator info
      ierr = MPI_Comm_rank(MPI_COMM_WORLD,&rank); CHKERRQ(ierr);
      if(debugid<0){stallForDebugger(rank);} 
      else if(rank==debugid){ stallForDebugger(rank); }
      ierr = MPI_Barrier(PETSC_COMM_WORLD);
  }

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

  PetscLogStagePush(AppSolve::logstages[0]); // Initialization 

  // Create the mesh.
  Mesh mesh(dim); 
  std::cout << "reading mesh..." << std::endl;
  mesh.read(controlfile("compexec/meshdata","./mesh.e"));

  // This class handles all the details of mesh refinement and coarsening.
  MeshRefinement mesh_refinement (mesh);
  // Possibly uniformly refine the mesh 
  mesh_refinement.uniformly_refine (user->Refine());

  // initialize mesh domains
  ierr = InitMeshDomains(mesh,controlfile); CHKERRQ(ierr);

  // Create an equation systems object 
  EquationSystems equation_systems (mesh); 

  // store controlfile to be passed in when models are initialized
  equation_systems.parameters.set<GetPot*>("controlfile") = &controlfile;

  // Setup the FEMSystem for the forward solve
  TransientFEMSystem* system=user->SetupFEMSystem(equation_systems);

  // setup the quantities of Interest
  user->SetupQOI(equation_systems);

  // Print information about the mesh to the screen.
  mesh.print_info();
  user->printSelf(std::cout);

  // Create a performance-logging object for this example
  PerfLog perf_log("PDE Solver");
  
  // We also set a standard linear solver flag in the EquationSystems object
  // which controls the maxiumum number of linear solver iterations allowed.
  equation_systems.parameters.set<unsigned int>("nonlinear solver maximum iterations") = 250;

  // At the beginning of each solve, reset the linear solver tolerance
  // to a "reasonable" starting value.
  const Real initial_solver_tol = 1.e-6;
  equation_systems.parameters.set<Real> ("nonlinear solver tolerance") = initial_solver_tol;

  // pass a pointer to the solvers
  equation_systems.parameters.set<AppSolve*>("AppSolve") = user;

  // memory used...
  PetscLogDouble memoryUsed;
  ierr = PetscMemoryGetCurrentUsage(&memoryUsed);
  std::cout << "Memory Used " << memoryUsed <<  std::endl << std::flush;
  // Initialize the data structures for the equation system.
  std::cout << "Initializing Equations Systems..." <<  std::endl << std::flush;
  equation_systems.init(); 

  // Prints information about the system to the screen.
  equation_systems.print_info(); 
  // write out input parameters
  equation_systems.parameters.print();

  PetscLogStagePop();     // Initialization 

  // optimize parameters and solve 
  /* begin multiple optimization loop
     COUNTING STARTS FROM ZERO!!!!!!!!!!!!!!!
     LEAVE THE COUNTING STARTING FROM ZERO!!!!!!!!!!!!!!!!!!!
     the logic of MANY of the routines that follow depend on it!!!!!!  
     IT TOOK TWO DAYS TO DEBUG THE FIRST TIME AROUND!!!!!!!!!!!!!!!! */
  std::cout << "begin solver..." << std::endl << std::flush;
  for(user->IDOPT = 0;user->IDOPT < user->NoptSteps() ; user->IDOPT++)
    {
         /* reset the check point counter before the optimization solve */
         user->UpdateCheckPoint(0);

         /* get the time step info for this optimization step
            NOTE this is called before the number of 
            optimization variables is determined       
            assuming constant time step for now        */ 
         qoiBaseClass *qoiOptimizer= user->qoiOptimizer(); 
         qoiOptimizer->SetNsteplo( user->getFEMTimeStep( user->IdealNzero() ) );
         qoiOptimizer->SetNstephi( user->getFEMTimeStep( user->IdealNtime() ) );

         // Prints information about the QOI to the screen.
         qoiOptimizer->printSelf();

         /* solve the optimization problem */
         ierr = qoiOptimizer->solveqoi(user,controlfile);CHKERRQ(ierr);

         //update time window data structures
         user->UpdateIdealWindow();
    }

  // clean up
  for(user->IDOPT = 0;user->IDOPT < user->NumQOI() ; user->IDOPT++)
    {
      qoiBaseClass *qoiOptimizer= user->qoiOptimizer(); 
      //FIXME need to reset data structures...
      ierr = qoiOptimizer->DestroyPetsc(); CHKERRQ(ierr);
    }
  std::ostringstream profileFile;// filename for profile
  profileFile << "profile."<<AppSolve::profileID<<".log";
  PetscLogPrintSummary(PETSC_COMM_WORLD,profileFile.str().c_str()); 

  // All done.  
  TaoFinalize(); 
  return 0;
}

// generate CSI thermal images
#undef __FUNCT__
#define __FUNCT__ "GenerateCSITmap"
PetscErrorCode GenerateCSITmap()
{
  PetscFunctionBegin;
//
//  // file names
//  std::vector< std::vector<std::string> > filenames( 2 * necho , 
//                             std::vector< std::string >::vector(nslice,"") );
//
//  // We construct one instance of the series reader object. Set the DICOM
//  // image IO object to be use with it, and set the list of filenames to read.
//  ReaderType::Pointer  reader = ReaderType::New();
//  ImageIOType::Pointer gdcmIO = ImageIOType::New();
//  reader->SetImageIO( gdcmIO );
//
//  // get header info
//  std::vector< PetscScalar > echotime(necho,0.0), imagfreq(necho,0.0);
//  for (int jjj = 0 ; jjj < necho ; jjj ++ )
//    {
//     // generate list of file names
//     RealTimeGenerateFileNames(0,jjj,filenames);
//
//     // only need to read in the first set to get the header info
//     reader->SetFileNames(filenames[0] );
//     InputImageType::Pointer tmpImage; // image pointer for header info
//     tmpImage= GetImageData(reader, rank, filenames[0]);
//
//     // scratch storage for header value
//     std::string value; 
//     // Get Echo Time 
//     std::string echotimekey = "0018|0081";
//     try
//     {  
//        gdcmIO->GetValueFromTag(echotimekey, value) ;
//        echotime[jjj]=boost::lexical_cast<double>(value);
//     }
//     catch(const std::exception& e) //catch bad lexical cast
//     {
//        std::cout<<"Error getting echo time " << " (" << echotimekey << ")\n";
//        std::cout<<"value returned is "<<value<< "\n";
//        std::cout<<e.what() << std::endl; abort();
//     }
//     std::string imagfreqkey = "0018|0084";
//     try
//     {  
//        gdcmIO->GetValueFromTag(imagfreqkey, value) ;
//        // trailing space on string causing cast error
//        value.erase(value.find(' ')) ; 
//        imagfreq[jjj]=boost::lexical_cast<double>(value);
//     }
//     catch(const std::exception& e) //catch bad lexical cast
//     {
//        std::cout << "Error getting Imaging Freq"<<" ("<<imagfreqkey << ")\n";
//        std::cout << "value returned is "<<value << "\n";
//        std::cout << e.what() << std::endl;
//        ierr = PetscOptionsGetScalar(PETSC_NULL,"-imagfreq",&imagfreq[jjj],PETSC_NULL);
//        CHKERRQ(ierr);
//        std::cout << "using value from command line "<< imagfreq[jjj] << "\n";
//     }
//     // echo data
//     std::string studyIdkey  = "0020|0010" ;
//     std::string seriesIdkey = "0020|0011" ;
//     gdcmIO->GetValueFromTag(studyIdkey  , value) ;
//     if(!rank) std::cout << "Study Id " << value  ;
//     gdcmIO->GetValueFromTag(seriesIdkey , value) ;
//     if(!rank) std::cout << " Series Number " << value 
//                         << " echo time " << echotime[jjj]<< " (ms)"<<std::endl;
//    }
//
//  //Define gamma*B0 and the echo spacing
//  PetscScalar      gB0=imagfreq[0], esp=echotime[1]-echotime[0];
//  if(!rank) std::cout << "echo spacing " << esp   << " (ms) "
//                      << "imaging freq " << gB0   << " (MHz) "
//                      << "ntime  "       << ntime << std::endl;
//
//  if(!rank) std::cout << "   ****Header Info****   " << std::endl << std::endl ;
//  // containers for real and imaginary data
//  VecFilterType::Pointer   realImageFilter=VecFilterType::New(),// real images
//                           imagImageFilter=VecFilterType::New();// imag images
//
//  // pointers to real and imaginary vector images for CSI computation
//  VecImageType::Pointer realImages, imagImages;
//
//  // loop over time instances
//  for( int iii = 0 ; iii <= ntime ; iii++)
//   {
//
//    // get image data
//    realImages = GetVecEchoImage(realImageFilter,iii,0); //real images
//    imagImages = GetVecEchoImage(imagImageFilter,iii,1); //imag images
//
//    // Software Guide : BeginLatex
//    //
//    // \index{Iterators!construction of} \index{Iterators!and image regions} The
//    // necessary images and region definitions are now in place.  All that is
//    // left to do is to create the iterators and perform the copy.  Note that
//    // image iterators are not accessed via smart pointers so they are
//    // light-weight objects that are instantiated on the stack.  Also notice how
//    // the input and output iterators are defined over the \emph{same
//    // corresponding region}.  Though the images are different sizes, they both
//    // contain the same target subregion.
//    //
//    // Software Guide : EndLatex
//    
//    // Software Guide : BeginCodeSnippet
//
//    // initialize ITK iterators
//    VecIteratorType realIt( realImages, realImages->GetRequestedRegion() );
//    VecIteratorType imagIt( imagImages, imagImages->GetRequestedRegion() );
//    realIt.SetFirstDirection(  0 );   realIt.SetSecondDirection( 1 );
//    imagIt.SetFirstDirection(  0 );   imagIt.SetSecondDirection( 1 );
//
//    realIt.GoToBegin(); imagIt.GoToBegin();
//
//    // get pointer to petsc data structures
//    PetscScalar   ****MapPixel;
//    DAGetGlobalVector(dac,&globalVec);
//    DAVecGetArrayDOF(dac,globalVec,&MapPixel);
//    //ierr = PetscObjectSetName((PetscObject)MapPixel,"MapPixel");CHKERRQ(ierr);
//
//    PetscInt PixelCounter = 0 ; 
//    /*
//       loop through parallel data structures
//    */
//    for (PetscInt k=ProcStart[2]; k<ProcStart[2]+ProcWidth[2]; k++) 
//    {
//      for (PetscInt j=ProcStart[1]; j<ProcStart[1]+ProcWidth[1]; j++) 
//      {
//        for (PetscInt i=ProcStart[0]; i<ProcStart[0]+ProcWidth[0]; i++) 
//        {
//          //std::cout << "pixel # " << PixelCounter++ << std::endl;
//          //for(PetscInt ivar = 0; ivar < nvarplot ; ivar++) 
//          //          MapPixel[k][j][i][ivar] = realIt.Get()[ivar] ;
//          ierr= PetscMatlabEnginePutArray(PETSC_MATLAB_ENGINE_WORLD,
//                                          necho,1,&realIt.Get()[0],
//                                                     "RealPixel");CHKERRQ(ierr);
//          ierr= PetscMatlabEnginePutArray(PETSC_MATLAB_ENGINE_WORLD,
//                                          necho,1,&imagIt.Get()[0],
//                                                     "ImagPixel");CHKERRQ(ierr);
//          ierr= PetscMatlabEnginePutArray(PETSC_MATLAB_ENGINE_WORLD,
//                                          nvarplot,1,MapPixel[k][j][i],
//                                                     "MapPixel");CHKERRQ(ierr);
//          ierr= PetscMatlabEngineEvaluate(PETSC_MATLAB_ENGINE_WORLD,
//                "MapPixel=MapCSI(RealPixel,ImagPixel,%18.16e,%18.16e)",gB0,esp);
//          CHKERRQ(ierr);
//          ierr= PetscMatlabEngineGetArray(PETSC_MATLAB_ENGINE_WORLD,
//                                          nvarplot,1,MapPixel[k][j][i],
//                                                     "MapPixel");CHKERRQ(ierr);
//          ++realIt; ++imagIt; // update iterators
//          CHKMEMQ; // check for memory corruption use -malloc_debug to enable
//        }
//        // get next line
//        realIt.NextLine(); imagIt.NextLine();
//      }
//      // get next slice
//      realIt.NextSlice(); imagIt.NextSlice();
//    }
//
//    /*
//       Restore vector
//    */
//    DAVecRestoreArrayDOF(dac,globalVec,&MapPixel);
//    DARestoreGlobalVector(dac,&globalVec);
//   
//    // gather the image buffer on processor zero
//    VecScatterBegin(gather,globalVec,imageVec,INSERT_VALUES,SCATTER_FORWARD);
//    VecScatterEnd(gather,globalVec,imageVec,INSERT_VALUES,SCATTER_FORWARD);
//  
//    // output various maps
//    if(!rank) // write image from root process
//     {
//       std::ostringstream unfiltered_filename;
//       unfiltered_filename << OutputDir <<"/unfiltered."<< rank << ".";
//       OSSRealzeroright(unfiltered_filename,4,0,iii);
//       unfiltered_filename << ".vtk" ;
//       ierr= WriteDAImage(unfiltered_filename); CHKERRQ(ierr);
//     }
//
//   } // end loop over time instances
//
  PetscFunctionReturn(0);
}
