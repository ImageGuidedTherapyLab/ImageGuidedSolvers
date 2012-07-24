// C++ include files 
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <vector>

// C include files 
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>

// Petsc/Tao include files
#include "tao.h" 
#include "src/tao_impl.h" // tao solvers

// libMesh include files
#include "libmesh.h"
#include "mesh.h"
#include "mesh_function.h"
#include "equation_systems.h"
#include "exact_solution.h"
// The nonlinear solver and system we will be using
#include "linear_implicit_system.h"
#include "nonlinear_implicit_system.h"
#include "nonlinear_solver.h"
//#include "string_to_enum.h" // polynomial order definitions
#include "getpot.h" // file parser

// dddas include files
#include "fortrandftranslate.h" // header to call fortran routines
#include "parser_defaults.h" // ensure same default values between C and Fortran
#include "mri_params.h" // mri params
#include "solvehp3d.h" // main header
#include "dddas.h"
#include "dddas_control.h"
#include "read_mridata.h"
#include "variable_map.h" // ensure same default values between C and Fortran


//#include <boost/lexical_cast.hpp>
//using boost::lexical_cast;

static  char help[] = "This is the driver routine for the\n\
Computations \n";
/* -------------- User-defined routines ---------- */
#include "baseInfo.h" // basic information
#include "transient_inverse_system.h"
#include "applicationContext.h" // interface to optimization routines for pennes eqn 
#include "pdeSolver.h"  // pde solver functions
#include "pennes_optimization.h" // interface to optimization routines for pennes eqn 
#include "pennesVerification.h"  // interface to pennes solver 

extern FORTRAN_FUNCTION 
{
#include "pennes_model.h" // interface to pennes model 
#include "global_params.h" // interface to fortran global functions
}
            
// global integers 
static PetscInt nzero = 0 ,
                none  = 1 ;
/**********************************************
    Main Routine for Computational Groups
 **********************************************/
#undef __FUNCT__
#define __FUNCT__ "Compute_Drone"
PetscErrorCode Compute_Drone(int argc,char **argv, GetPot &controlfile,
                             UniversalCntrl &GenInfo, DroneControl &QOIInfo,
                                              Mesh &mesh,  PetscInt GroupID )
{
  PetscErrorCode info; /* used to check for functions returning nonzeros */
  PetscFunctionBegin; // used for error handling

//  // allocate memory for element matrices
//  info = PetscMalloc(sizeof(PetscScalar)*MAXDOF*MAXDOF,&ELEMMATX);
//  info = PetscMalloc(sizeof(PetscScalar)*MAXDOF,&ELEMLOAD);
//  info = PetscMalloc(sizeof(PetscInt)*MAXDOF,&GLOBMAP);
//
//  /* reserve room to hold the middle maps owned by this task */
//  GenInfo.MdleData.reserve(FORTRAN_NAME(get_maxelib)()/GenInfo.petscsize);
//
//  // allocate the scratch storage used in through out the element routines
//  FORTRAN_NAME(allocateelemscratch)(&QOIInfo.NumOptVar);

  // setup imaging and power data structures 
  // Sequence is IMPORTANT!!!! Setup Power uses data setup in imaging
  dddas::Image Image(GenInfo);
  info = SetupPower(GenInfo);

  /* instantiate user-defined application context to faclitate Petsc Solve */
  AppSolve     user(GenInfo,QOIInfo,mesh,Image,GroupID); 

  // Set pointer that is global to the pde Solver
  setApplicationContext(&user);
      
  PetscPrintf(PETSC_COMM_WORLD, 
            "@@@@@@@@@@@@@@@@@group %d data structures initialized \n",GroupID);


  // establish intercommunicator between this group and the Control Task
  // this intercommunicator will be used by the data thread
  MPI_Intercomm_create(PETSC_COMM_WORLD, 0, user.GenInfo->DDDAS_COMM_WORLD, 0, 
                                             user.GroupID,&user.DataComm);
  // merge the intercommunicators to create the communicators for
  // the control/vis 
  MPI_Intercomm_merge(user.DataComm,1, &user.ControlComm); 

  PetscPrintf(PETSC_COMM_WORLD, 
              "@@@@@@@@@@@@@@@@@group %d communicators initialized \n",GroupID);
  /*Get the # of threads per tasks on an SMP machine*/
  PetscTruth  flg;
  PetscInt    ncompthrd=1; // omp threads
  info=PetscOptionsGetInt(PETSC_NULL,"-ncompthrd",&ncompthrd,&flg);
  CHKERRQ(info);

  /* Initialize TAO */
  TaoInitialize(&argc,&argv,(char *)0,help);

  // initialize function pointers for Optimization solve;
  info  = InitControlGlobals(QOIInfo); CHKERRQ(info);

  PetscLogStagePush(logstages[1]); //Begin Parallel Infrastructure Stage
  /**
   * Create two equation systems objects: 
   *       main_es- used for the forward solve and handling the ideal 
   *                solution. this is the main object that will be plotted.
   *       work_es- used for the adjoint solve, the sensitivity solve and
   *                any other data structures (exact adjoint) needed for 
   *                gradient/hessian computation and verification problems
   *                this object should only need to be plotted during debugging.
   */
  EquationSystems main_es(mesh), 
                  work_es(mesh);

  // set pointer to libmesh EquationSystems within application context
  // to call libmesh solvers from optimization routines
  user._main_es = &main_es;
  user._work_es = &work_es;


  PetscLogEventBegin(logevents[21],0,0,0,0); // setup libMesh

  if( QOIInfo.pde.find("verifprob"            )!=std::string::npos ||
      QOIInfo.pde.find("nonlinpennesmonte"    )!=std::string::npos ||
      QOIInfo.pde.find("pennesisolasercooling")!=std::string::npos ||
      QOIInfo.pde.find("nonlinpennesisolaser" )!=std::string::npos)
    {
      std::cout <<"Setting up Verification Problems with SDA"<<std::endl;
      user.pdeSolver = new PennesVerification(controlfile,baseInfo::n_block);
    }
  else if( QOIInfo.pde.find("nonlinpennessda")!=std::string::npos )
    {
      std::cout <<"Setting up Bioheat with SDA"<<std::endl;
      user.pdeSolver = new PennesSDA(controlfile,baseInfo::n_block);
    } 
  else if( QOIInfo.pde.find("nonlinpennesrf")!=std::string::npos )
    { 
      std::cout <<"Setting up Bioheat with RF"<<std::endl;
      user.pdeSolver = new PennesRadioFrequency(controlfile,baseInfo::n_block);
    } 
  else if( QOIInfo.pde.find("tgm")!=std::string::npos )
    { 
       std::cout << "not implemented yet! \n"; abort();
    } 
  else if( QOIInfo.pde.find("pennesdeltap1")!=std::string::npos )
    { 
      std::cout <<"Setting up Bioheat with Delta P1"<<std::endl;
      user.pdeSolver = new PennesDeltaP1(controlfile,baseInfo::n_block);
    } 
  else if( QOIInfo.pde.find("pennesdirectfluence")!=std::string::npos )
    { 
      std::cout <<"Setting up Coupled Bioheat Direct Fluence solver"<<std::endl;
      user.pdeSolver = new PennesFluencePDE(controlfile,baseInfo::n_block);
    } 
  else 
    { 
       std::cout << "unknown pde solver \n"; abort();
    } 

  // setup the system variables
  user.pdeSolver->SetupState(user);


  CHKMEMQ; // check for memory corruption use -malloc_debug to enable

  PetscPrintf(PETSC_COMM_WORLD, 
              "@@@@@@@@@@@@@@@@@group %d data initialized\n",GroupID);


  // Declare the exact solution 
  ExactSolution exact_sol(main_es);
  user._exact_sol = &exact_sol;

  // set initial time
  main_es.parameters.set<Real> ("Time") = 0.0 ; 

  // Set the user-specifiied nonlinear solver tolerance
  main_es.parameters.set<Real>("nonlinear solver tolerance") =  1.e-9;

  // Set the user-specified maximum # of nonlinear solver iterations
  main_es.parameters.set<unsigned int>("nonlinear solver maximum iterations") = 25;

  PetscLogEventEnd(logevents[21],0,0,0,0);   // setup libMesh
  			 
  /*  Initialize the data structures for the equation system, 
      setup initial conditions, etc... 
                               pennes_init_cd  
                               pennes_init_rf  */
  // user.Nsteplo is used to set the IC time point (not necessarily at idstep=0)
  user.Nsteplo = QOIInfo.IDEAL_NZERO * baseInfo::ISTEPS_PER_IDEAL;

  PetscLogEventBegin(logevents[20],0,0,0,0); // init libMesh
  main_es.init();
  PetscLogEventEnd(logevents[20],0,0,0,0);   // init libMesh

  // Print information about the system to the screen.
  main_es.print_info();
  PetscLogStagePop(); // End Parallel Infrastructure Stage

  /* Initial data structures ready, now get message from 
              control task on which ini file to read in and execute */

//  // get initial domain decomposition
//  MdleInfo mdletmp;
//  PetscInt    mdle=0;
//  for(iii = 0 ; iii < FORTRAN_NAME(get_nreleb)() ; iii++){
//     FORTRAN_NAME(nelconb)(&mdle,&mdle); 
//     if(GenInfo.petscrank == FORTRAN_NAME(irankmdlenode)(&mdle)){
//        mdletmp.mdle   = mdle;
//        //return element id correspoding to the mdle node
//        mdletmp.idelem = FORTRAN_NAME(get_mdle_elem)(&mdle);
//        GenInfo.MdleData.push_back(mdletmp); 
//        FORTRAN_NAME(computeglobalmaps)(&mdletmp.mdle,&mdletmp.idelem,&idstate);
//     }
//  }
//  vector<MdleInfo> *MdleData = &GenInfo.MdleData; 
//
//  /* setup the characteristic function. mainly used to reduce the
//      domain of interest for the qoi in goal-oriented error estimation */
//  for(mdleID=MdleData->begin() ; mdleID!=MdleData->end() ; mdleID++){
//     FORTRAN_NAME(field2datastruct)(&mdleID->mdle,&mdleID->idelem,&nzero,
//                   &FORTRAN_NAME(getcharfield),&FORTRAN_NAME(elemchar)); 
//  }
//
//  /* before the multiple optimization loop setup the initial conditions 
//     to start propagating temperature forward in time */
//  // get the initial FEM time instance
//  for(mdleID=MdleData->begin() ; mdleID!=MdleData->end() ; mdleID++){
//     FORTRAN_NAME(field2datastruct)(&mdleID->mdle,&mdleID->idelem,
//                    &user.Init_Nsteplo,&FORTRAN_NAME(getinittemp),
//                                       &FORTRAN_NAME(elemzdofs)); 
//  }
//  PetscPrintf(PETSC_COMM_WORLD,
//                     "@@@@@@@@group %d setup initial conditions \n",GroupID);

  // need to setup imaging data structures in this is a calibration problem
  //  OR if initial conditions are taken from thermal imaging data
  if(QOIInfo.compobj.find("calibration")!=std::string::npos || GenInfo.IC_MRTI)
    {
      // read dicom header info if necessary
      Image.ImageSetup(GenInfo);
      Image.SetImagePointers(main_es.parameters);
    }
  // to compute the objective function over the entire time history while 
  //  optimizing over a smaller time window, preload all data to the ideal data 
  //  structures to begin with. this is used only for parameter studies
  if(GenInfo.COMPUTEFINALOBJECTIVE){
     PetscInt  save_IDEAL_NZERO = QOIInfo.IDEAL_NZERO; // save IDEAL_NZERO
     PetscInt  save_IDEAL_NTIME = QOIInfo.IDEAL_NTIME; // save IDEAL_NTIME

     QOIInfo.IDEAL_NZERO = baseInfo::MRTI_nzero; 
     QOIInfo.IDEAL_NTIME = baseInfo::MRTI_ntime;
     user.Nsteplo = QOIInfo.IDEAL_NZERO * baseInfo::ISTEPS_PER_IDEAL;
     user.Nstephi = QOIInfo.IDEAL_NTIME * baseInfo::ISTEPS_PER_IDEAL;

     /* load the data structures for the QOI
        LoadQoiData is a function pointer to  
                            verifcomputeidealfield: for verification problems
                            readMRTI_data         : for calibration problems
                            verifcomputeidealfield: for verification problems
                            readMRTI_data    : for calibration problems
                            LoadIdealArr     : arrhenius base optimize
                            LoadIdeal_TS     : two-state base optimize
                            LoadIdealHSP     : HSP base optimize
                            LoadIdealTmp     : for temp base optimize */
     info=QOIInfo.PriorLoadQoiData(&user);CHKERRQ(info);
     QOIInfo.IDEAL_NZERO = save_IDEAL_NZERO ; // return orig IDEAL_NZERO
     QOIInfo.IDEAL_NTIME = save_IDEAL_NTIME ; // return orig IDEAL_NTIME

     // if the times don't align then set thermal image as IC
     if ( QOIInfo.IDEAL_NZERO > baseInfo::MRTI_nzero )
                                   Image.loadMRTItempfield(&user);
  }
  if(baseInfo::NEXACT>=0)
    {
     /* begin multiple optimization loop
        COUNTING STARTS FROM ZERO!!!!!!!!!!!!!!!
        LEAVE THE COUNTING STARTING FROM ZERO!!!!!!!!!!!!!!!!!!!
        the logic of MANY of the routines that follow depend on it!!!!!!  
        IT TOOK TWO DAYS TO DEBUG THE FIRST TIME AROUND!!!!!!!!!!!!!!!! */
     for(user.IDOPT = 0;user.IDOPT < QOIInfo.NoptSteps ; user.IDOPT++)
       {

         /* THE SAME IC is used for all optimization iterations.
               The first optimization step ALWAYS uses the
               given initial condition. On subsequent optimization 
               steps the IC_MRTI flag can be used to load the
               initial conditions from the thermal images 
           !@!NOTE!@!  that this is called before the communication
           with the control task. This makes the Compute Drone 
           wait to recieve its model parameters until the thermal image 
           is available. 
           At the point in time when the thermal image is available
           the control task should have had time to accumulate at least
           one checkpoint of the solutions from other Compute_Drones.
         */
         if(user.IDOPT && GenInfo.IC_MRTI) Image.loadMRTItempfield(&user);

         // wait until full set of calibration data is available
         //  so that can wait for real-time treatment user input
         //  power update commands
         if(QOIInfo.compobj.find("calibration")!=std::string::npos)
                                           Image.checkMRTI_data4Pow(&user);

         // tell control task that we are ready for the parameters in the
         // next control file
         PetscInt commands[2] = { 0 , SETUP };
         if( !libMesh::processor_id() ) 
           MPI_Send(commands,2,MPI_INT,0,user.GroupID,
                                         GenInfo.DDDAS_COMM_WORLD); 
         PetscPrintf(PETSC_COMM_WORLD,
              "@@@@@@@@group %d sent control task setup command\n",GroupID);

         // get file id number of global elements from control server
         MPI_Bcast(&commands,2,MPI_INT,0,user.ControlComm); 
         GenInfo.controlfileID = commands[0];
         GenInfo.outputfileID  = commands[1];
         PetscPrintf(PETSC_COMM_WORLD,
              "@@@@@@@@group %d received commands from control task\n",GroupID);

         // to pass to fortran
         PetscScalar spacing[3]; Image.GetSpacing(spacing);

         // read in control parameters for fortran code
         FORTRAN_NAME(read_control_fortran)(&GenInfo.controlfileID,spacing);

         //read in entire power history 
         FORTRAN_NAME(read_power)(&GenInfo.controlfileID,&nzero);
         // read in field info
         //FORTRAN_NAME(read_w_0field)(&GenInfo.controlfileID);
         //FORTRAN_NAME(read_k_0field)(&GenInfo.controlfileID);
//
//
//         // Setup the Hp3d Data structures
//         info = SetupHp3d(GenInfo,user.ControlComm, &QOIInfo.recvcnts_elm,
//                                                    &QOIInfo.recvcnts_for,
//                                                    &QOIInfo.recvcnts_adj,
//                                                    &QOIInfo.displs_elm,
//                                                    &QOIInfo.displs_for,
//                                                    &QOIInfo.displs_adj);
//         CHKERRQ(info);
//         PetscPrintf(PETSC_COMM_WORLD,
//              "@@@@@@@@group %d setup hp3d data structures done \n",GroupID);
//
//         /*finally, to setup the communication between hp3d and 
//           PETSC DATA STRUCTURES create a mapping from local dof
//           to global dof which will then be used to setup global 
//           matrices and vectors also Setup Solver Contexts */
//         info = SetupPetscDstruc(&user); CHKERRQ(info);
//
//
//         PetscPrintf(PETSC_COMM_WORLD,
//              "@@@@@@@@group %d setup petsc data structures done\n",GroupID);
//
//

         /* reset the check point counter before the optimization solve */
         GenInfo.checkpoint = GenInfo.checkpointiter ;

         /* solve the optimization problem */
         info=SolveOptimizationProb(user,QOIInfo); CHKERRQ(info);

         //update time window data structures
         QOIInfo.IDEAL_NZERO = QOIInfo.IDEAL_NZERO + QOIInfo.NOFFSET ;
         QOIInfo.IDEAL_NTIME = QOIInfo.IDEAL_NTIME + QOIInfo.NUMIDEAL;

//         PetscLogStagePush(logstages[1]); //Parallel Infrastructure 
//         //update time window data structures
//         info = CleanPetscDstruc(&user); CHKERRQ(info);
//         PetscLogStagePop();
       }
    }
  else
    { // only checking domain decomposition
//     if(!GroupID){ //only send from group 0
//        // Setup the Hp3d Data structures
//        info = SetupHp3d(GenInfo,user.ControlComm, &QOIInfo.recvcnts_elm,
//                                                   &QOIInfo.recvcnts_for,
//                                                   &QOIInfo.recvcnts_adj,
//                                                   &QOIInfo.displs_elm,
//                                                   &QOIInfo.displs_for,
//                                                   &QOIInfo.ispls_adj);
//        CHKERRQ(info);
//        PetscPrintf(PETSC_COMM_WORLD,
//             "@@@@@@@@group %d setup hp3d data structures done \n",GroupID);
//        /* gather field info */
//        /* allocate buffer for inter groups field transmission */
//        PetscInt iii,nelems = GenInfo.MdleData.size() ;
//        FORTRAN_NAME(alloc_hpelemfield)( &nelems );
//        /* pack the  buffer for inter groups field transmission */
//        for(mdleID=MdleData.begin() ; mdleID!=MdleData.end() ; mdleID++){
//           iii = distance(MdleData.begin(),mdleID );
//           FORTRAN_NAME(pack_hpelemfield)(&iii,&mdleID.idelem,&mdleID.mdle);
//        }
//        /* setup for varying perfusivity and conductivity fields */
//        PetscMPIInt fortran_comm = PetscFromPointerComm(user.ControlComm);
//        FORTRAN_NAME(gatherdealloc_hpelemfield)(&fortran_comm,
//                                                &QOIInfo.optimize_W_0,
//                                                &QOIInfo.optimize_K_0,
//                                                &QOIInfo.ntask,
//                                                &QOIInfo.recvcnts_elm[0],
//                                                &QOIInfo.displs_elm[0] );
//     }
    }

  if(GenInfo.COMPUTEFINALOBJECTIVE){
     PetscInt save_Nsteplo = user.Nsteplo; // save Nsteplo
     PetscInt save_Nstephi = user.Nstephi; // save Nstephi

     /* evaluate the error ||u_ideal-u|| over the entire time history */
     user.Nsteplo = baseInfo::getMinFEMTimeStep();
     user.Nstephi = baseInfo::getMaxFEMTimeStep();

     if ( QOIInfo.IDEAL_NZERO > baseInfo::MRTI_nzero )
      { // recomputation is needed w/ a separate plot
        info = ForwardSolve(&user); CHKERRQ(info); // forward time stepping

        // initialize file output of main data structures
        ExodusII_IO output_vis(*user.mesh);

        // solution buffer
        std::vector<Number>      soln;

        // set main solution filename
        OStringStream file_name;
        file_name.str(""); // reset before reuse
        file_name << "femvis/femqoi_"<<user.GroupID<<"nopt_"<<user.IDOPT;
        // OSSRealzeroright(file_name,4,0, idplot);
        file_name << ".e";

        // get names form the system
        std::vector<std::string> names;
        for (unsigned int i_sys = 0 ; i_sys < 
                                      user._main_es->n_systems() ; i_sys++)
          {  
            System& system  = user._main_es->get_system(i_sys) ; 
            const unsigned int nv_sys = system.n_vars();
            for (unsigned int var=0; var<nv_sys; var++)
             names.push_back( system.variable_name(var) );
          }  

        //plot from beginning of optimization solve to the final time
        for(baseInfo::ISTEP  = user.Nsteplo;
            baseInfo::ISTEP <= baseInfo::getMaxFEMTimeStep() ; 
            baseInfo::ISTEP++)
         if( baseInfo::modWrite(baseInfo::ISTEP) ) //write every nth step
          {
           PetscInt idplot = baseInfo::ISTEP / baseInfo::WriteInterval ; 
           // Main_Visualization(user.ControlComm,user, &user.GroupID,
           //                    &user.IDOPT, &idplot, &baseInfo::ISTEP, 
           //                    GenInfo,QOIInfo,user.x, user.p);
           
           // set time
           Real Time = baseInfo::getTime(baseInfo::ISTEP);
           user._main_es->parameters.set<Real> ("Time") = Time;
           
           // create data buffer
           build_transient_system_solution_vector(*user._main_es, soln);

           // Write out every ideal timestep to file.
           output_vis.write_nodal_data(file_name.str(),soln,names);
           // write timestep info to file 
           output_vis.write_timestep(Time);
           // increment timestep (1-based) 
           output_vis.increment_timestep();

          } // end loop over time steps
      }

     /* evaluate the error */
     PetscScalar qoi = QOIInfo.ComputeObjective(&user);

     PetscPrintf(PETSC_COMM_WORLD,
      "@@@@group %d full error: ||u_ideal-u||_(%d,%d)=%e\n",
                                   GroupID, user.Nsteplo, user.Nstephi, qoi);

     ///* ensure that the temperature field stored is equal to ZERO
     //     i.e.    ComputeObjective = ||u_ideal - u || 
     //                              = ||u_ideal ||   (u=0) */
     //for(PetscInt idstep= user.Nsteplo ; idstep <= user.Nstephi ; idstep++){
     //   for(mdleID  = GenInfo.MdleData.begin() ; 
     //       mdleID != GenInfo.MdleData.end()   ; mdleID++){
     //      FORTRAN_NAME(field2datastruct)(&mdleID.mdle,&mdleID.idelem,
     //                                &idstep,&FORTRAN_NAME(getzerotemp),
     //                                        &FORTRAN_NAME(elemzdofs)); 
     //   }
     //}
     ///* evaluate ||u_ideal|| over the full time history 
     //      used to scale space time norm                 */
     //qoi = SpaceTime_QOI(&user);

     //PetscPrintf(PETSC_COMM_WORLD,
     // "@@@@group %d full norm : ||u_ideal||_(%d,%d) = %e\n",
     //                              GroupID, user.Nsteplo, user.Nstephi, qoi);
     ///* evaluate ||u_ideal|| over the time history in the solve
     //      used to scale space time norm                 */
     user.Nsteplo = save_Nsteplo;
     user.Nstephi = save_Nstephi;
     //qoi = SpaceTime_QOI(&user);

     //PetscPrintf(PETSC_COMM_WORLD,
     // "@@@@group %d solve norm: ||u_ideal||_(%d,%d)  =%e\n",
     //                              GroupID, user.Nsteplo, user.Nstephi, qoi);
     ///* evaluate ||u_ideal|| at 0th time step
     //             used to scale L2 norm                 */
     //// zero last ideal field data entry
     //idstep  = GenInfo.MAXIDEAL * GenInfo.ISTEPS_PER_IDEAL;
     //for(mdleID  = GenInfo.MdleData.begin() ; 
     //    mdleID != GenInfo.MdleData.end()   ; mdleID++){
     //   FORTRAN_NAME(field2datastruct)(&mdleID.mdle,&mdleID.idelem,
     //                             &idstep,&FORTRAN_NAME(getzerotemp),
     //                                     &FORTRAN_NAME(elemideal)); 
     //}
     //Petsc Scalar qoi_loc = 0.0;
     ///* i.e.    ComputeObjective = ||u_ideal - u || 
     //                            = ||u_ideal ||   (u=0) */
     //for(mdleID=MdleData.begin() ; mdleID!=MdleData.end() ; mdleID++){
     //   qoi_loc = qoi_loc + FORTRAN_NAME(qoieval_lesbegue2)(&mdleID.mdle,
     //                                          &mdleID.idelem,&idstep);
     //}
     //MPI_Allreduce(&qoi_loc,&qoi, 1, MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD);
     //PetscPrintf(PETSC_COMM_WORLD,
     // "@@@@group %d solve norm: ||u_ideal||_(0,0)  =%e\n", GroupID, qoi);
  }
//
//  
//  info = PetscFree(ELEMMATX); CHKERRQ(info);
//  info = PetscFree(ELEMLOAD); CHKERRQ(info);
//  info = PetscFree(GLOBMAP); CHKERRQ(info);

  /* Finalize Tao */
  TaoFinalize();
  PetscFunctionReturn(0);
}
