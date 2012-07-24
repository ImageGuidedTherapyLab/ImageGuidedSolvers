/*$Id$*/

/* Program usage: mpiexec -np ? $EXEC [-help] [all TAO options] */

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

// libMesh include files
#include "libmesh.h"
#include "mesh.h"
#include "mesh_refinement.h"
#include "equation_systems.h"
#include "exact_solution.h"
#include "getpot.h" // file parser

// dddas include files
#include "fortrandftranslate.h" // header to call fortran routines
#include "parser_defaults.h" // ensure same default values between C and Fortran
#include "mri_params.h" // mri params
#include "baseInfo.h"  // global params header
#include "solvehp3d.h" // main header

#include <boost/lexical_cast.hpp>
using boost::lexical_cast;

/*T 
   Concepts: TAO - Solving an unconstrained minimization problem
   Routines: TaoInitialize(); TaoFinalize();
   Routines: TaoApplicationCreate(); TaoAppDestroy();
   Routines: TaoCreate(); TaoDestroy(); 
   Routines: TaoAppSetObjectiveAndGradientRoutine();
   Routines: TaoAppSetHessianMat(); TaoAppSetHessianRoutine();
   Routines: TaoSetOptions();
   Routines: TaoAppSetInitialSolutionVec();
   Routines: TaoSolveApplication();
   Routines: TaoGetTerminationReason();
T*/ 

/* -------------- User-defined routines ---------- */
  
#include "dddas.h" // main interface to dddas code 
#include "dddas_control.h" // interface to dddas control routines

extern FORTRAN_FUNCTION 
{
#include "global_params.h" // interface to fortran global functions
#include "plot_params.h" // interface to fortran pennes_model
//  void FORTRAN_NAME(main_par_initialize)(PetscInt*);
//  void FORTRAN_NAME(init_control)(PetscMPIInt*,PetscInt*,PetscScalar *);
//  void FORTRAN_NAME(init_compute_drone)(PetscScalar*);
}
///* -------------- For Debugging  ------------ */
volatile int endstall=0;
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


//static  char help[] = "This is the driver routine for the\n\ DDDAS project \n";

/* 
  The main routine breaks up the tasks into separate MPI communicators
   and reads in all data from the control.ini file
*/
#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  PetscErrorCode info; /* used to check for functions returning nonzeros */

  /* Initialize MPI*/
  MPI_Init(&argc,&argv);

  PetscInt    size,rank; // world communicator info
  info = MPI_Comm_size(MPI_COMM_WORLD,&size); CHKERRQ(info);
  info = MPI_Comm_rank(MPI_COMM_WORLD,&rank); CHKERRQ(info);

  // all processors read in control parameters from ini file 
  if(rank == 0) std::cout << "opening control file..."<< std::endl << std::flush; 
  GetPot controlfile("files/control.ini");
  if(rank == 0) std::cout << "control file read!"<< std::endl << std::flush; 

  // get the number of computational groups 
  //     ASSUME one computational group per quantity of interest
  PetscInt number_comp_group = controlfile("compexec/num_qoi",1);
  std::vector<PetscInt> leadcomp(number_comp_group,0);

  // starting rank of computations
  PetscInt  comp_rank_begin= controlfile("compexec/comp_rank_begin",1);


  // account for case that user inputs a value that is smaller than 1
  comp_rank_begin=(comp_rank_begin > 1) ? comp_rank_begin : 1;

  if(rank == 0){
   printf("comp_rank_begin= %d number_comp_group= %d\n",
           comp_rank_begin,    number_comp_group);
  }

  /* available tasks = number of world ranks MINUS 
                       the rank of the lead computation task on 1st group */
  PetscInt available_tasks = size - comp_rank_begin ; 
  if(available_tasks < number_comp_group){ 
    std::cout << available_tasks << " available tasks,must run with at least "
         << comp_rank_begin+number_comp_group 
         << " mpi tasks\n one task is set aside for control\n"
         << " one task is set aside as a data server\n "
         << comp_rank_begin-1<< " tasks are unused \n "
         << " the rest are for computation\n "
         << "    if only have one processor available,\n "
         << "    setup machine files accordingly \n";
    return 0; 
  }
  // get # of tasks per group
  std::vector<PetscInt> tasks_per_group(number_comp_group,0), 
                   tasks_accum(number_comp_group,0);
  PetscInt  tasks_dflt=available_tasks/number_comp_group; // default # of tasks per group
  if(rank == 0){
     std::cout <<"available_tasks="<< available_tasks
          <<" tasks_dflt= "   << tasks_dflt      << std::endl;
  }
  PetscInt  sum = 0 ; // keep track of sum for error checking
  // store ID's of lead ranks. rank ordering is with respect to DDDAS_COMM_WORLD
  PetscInt    leadcompid = 1; // task splitting info
  for (unsigned int iii = 0 ; iii < tasks_per_group.size() ; iii++) {

     std::ostringstream qoiID ;
     qoiID << "qoi_"<< iii << "/comptasks";
     tasks_per_group.at(iii) = controlfile(qoiID.str().c_str(),tasks_dflt);
     sum = sum + tasks_per_group.at(iii);
     tasks_accum.at(iii) = sum; // keep track of how many tasks have been used
     leadcomp.at(iii) = leadcompid;
     leadcompid = leadcompid + tasks_per_group.at(iii);
  }
  if(sum < available_tasks){ // put remaining tasks in last group
    PetscInt  ilast = tasks_per_group.size()-1;
    tasks_per_group.at(ilast) = available_tasks - tasks_accum.at(ilast-1);
    tasks_accum.at(ilast) = available_tasks; 
  } else if(sum > available_tasks){ // error 
    std::cout << "error distributing # of tasks per computational group\n"
         << "requested " << sum << "tasks only had "
         <<  available_tasks << "available\n";
    return 0; 
  }

  if(rank == 0){
     for(unsigned int iii = 0 ; iii < tasks_per_group.size() ; iii++) {
        std::cout <<"tasks_per_group["<<iii<<"]= "<< tasks_per_group[iii]
             <<" tasks_accum["   <<iii<<"]= "<< tasks_accum[iii]    <<" \n";
     }
  }

  /* instantiate the class that has general info used by all tasks */
  UniversalCntrl GenInfo(controlfile);

  PetscInt    dddasID, groupID;
  /* setup groups for subcommunicators of MPI_COMM_WORLD to dedicate a entire
     compute groups to thermal images, fem computation, and task control */
  if(rank == 0){ // task for control and writing
       GenInfo.Control_Task = PETSC_TRUE; //this is the control task;
       groupID = number_comp_group;
       dddasID = 1;
  }
  else if(rank >= comp_rank_begin){ // tasks for computations
       PetscInt  rel_rank =  rank-comp_rank_begin;
       groupID = 0;
       while( rel_rank >= tasks_accum.at(groupID) ) groupID++;
       dddasID = 1;
  }
  else { // rest of tasks exit and may be used in openmp
       groupID = MPI_UNDEFINED;
       dddasID = MPI_UNDEFINED;
  }

  if(groupID == number_comp_group){
     for(unsigned int iii = 0 ; iii < tasks_per_group.size() ; iii++){
        std::cout <<"group "<<iii<<" has "<<tasks_per_group.at(iii)<<" tasks\n"
             <<"group "<<iii
             <<" comp leader = dddas_rank "<<leadcomp.at(iii)<< std::endl;
     }
  }

  // datatype error check 
  if(sizeof(MPI_FLOAT) != sizeof(AVSDATA_TYPE)){
    std::cout <<"error in MPI AVS datatypes file \n"; abort();
  }
  
  /* create PETSC_COMM_WORLD as a subcommunicator of MPI_COMM_WORLD */ 
  info = MPI_Comm_split(MPI_COMM_WORLD,dddasID,rank,&GenInfo.DDDAS_COMM_WORLD);
  info = MPI_Comm_split(MPI_COMM_WORLD,groupID,rank,&PETSC_COMM_WORLD);

  if(PETSC_COMM_WORLD != MPI_COMM_NULL) 
   {
    PetscInt dddas_rank;
    info = MPI_Comm_rank(GenInfo.DDDAS_COMM_WORLD,&dddas_rank); CHKERRQ(info);
    std::cout<<"rank "    << rank    <<" dddas_rank "<< dddas_rank
        <<" groupID "<< groupID 
        <<" dddasID "<< dddasID << std::endl << std::flush ;

    // run DDDAS system
    DDDAS(argc,argv,GenInfo,controlfile,tasks_per_group,leadcomp,groupID);
   }
  else
   {
    std::cout<<"mpi world rank "<<rank<<" exiting..."<<std::endl<<std::flush;
   }

  MPI_Finalize(); 
  return 0; //  
}
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   DDDAS system, notice that the LibMeshInit class constructor/destructor
   is called as a subroutine of main, this allows initialization of 
   libMesh and petsc, separately from MPI
  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
int DDDAS(int argc,char **argv, UniversalCntrl &GenInfo,
          GetPot &controlfile, std::vector<PetscInt> &tasks_per_group, 
          std::vector<PetscInt> &leadcomp, PetscInt &groupID)
{

  LibMeshInit init (argc, argv,PETSC_COMM_WORLD);

  /* Initialize PETSC for Fortran USE */
  PetscErrorCode info; /* used to check for functions returning nonzeros */
  info =  PetscInitializeFortran(); CHKERRQ(info);
  info =  PetscLogBegin(); CHKERRQ(info);

  /*Flag for monitoring the process of the optimization*/
  PetscTruth  flg,tao_monitor=PETSC_FALSE;//monitoring flag
  info = PetscOptionsGetTruth(PETSC_NULL,"-tao_monitor",&tao_monitor,&flg);
  CHKERRQ(info);

  /*Check for debugging flags*/
  PetscTruth  debug=PETSC_FALSE;// debug flag
  info = PetscOptionsGetTruth(PETSC_NULL,"-idb",&debug,&flg); CHKERRQ(info);
  PetscInt  debugid=-1; // default to all
  info = PetscOptionsGetInt(PETSC_NULL,"-idbrank",&debugid,&flg); CHKERRQ(info);
  if(debug){ // debug
      PetscInt    rank; // world communicator info
      info = MPI_Comm_rank(MPI_COMM_WORLD,&rank); CHKERRQ(info);
      if(debugid<0){stallForDebugger(rank);} 
      else if(rank==debugid){ stallForDebugger(rank); }
      info = MPI_Barrier(GenInfo.DDDAS_COMM_WORLD);
  }

  info = PetscOptionsGetTruth(PETSC_NULL,"-ignore_bc",&GenInfo.IGNORE_BC,&flg); 
  CHKERRQ(info);

  // getInfo
  GenInfo.GenSetup(controlfile);
  baseInfo::GenSetup(controlfile);

  // Read in a 3D mesh. 
  Mesh mesh(3);    
  mesh.read(controlfile("compexec/meshdata","files/mesh.ucd"));

  // This class handles all the details of mesh refinement and coarsening.
  MeshRefinement mesh_refinement (mesh);

  // Possibly uniformly refine the mesh 
  mesh_refinement.uniformly_refine (baseInfo::Refine);

  // Use the MeshTools::Generation mesh generator to create a uniform
  // grid on the square [-1,1]^D.  We instruct the mesh generator
  // to build a mesh of 8x8 \p Quad9 elements in 2D, or \p Hex27
  // elements in 3D.  Building these higher-order elements allows
  // us to use higher-order approximation, as in example 3.
  //if (order != "FIRST")
  //  mesh.all_second_order();
  //

  mesh.print_info();

  PetscPrintf(GenInfo.DDDAS_COMM_WORLD,"@@@@@@@@@@@@@@@Mesh read\n");

  // variables for array bounds error checking on # of thermal images
  //               and error checking the optimization variables
  PetscInt numopt_K_0  =0 ,  numopt_K_1 =0 , 
           numopt_K_2  =0 ,  numopt_K_3 =0 , numopt_W_0  =0, 
           numopt_W_N  =0 ,  numopt_W_I =0 , numopt_W_D  =0,
           numopt_W_2  =0 ,  numopt_W_NI=0 , numopt_W_ID =0,
           numopt_X_0  =0 ,  numopt_MU_A=0 , 
           numopt_Y_0  =0 ,  numopt_MU_S=0 , 
           numopt_Z_0  =0 ,  numopt_POW =0 ;
                                             
  /* get the number of optimization steps that will be done for each qoi */
  std::vector<DroneControl> QOIinfo;
  for (unsigned int iii = 0 ; iii < tasks_per_group.size() ; iii++) 
    {

       DroneControl qoitmp(controlfile,GenInfo,iii);

       /* get group-wise maximum and minimum of the last ideal step */
       baseInfo::MAXIDEAL=baseInfo::MAXIDEAL < qoitmp.FinalOpt_Ntime ? 
                                           qoitmp.FinalOpt_Ntime : baseInfo::MAXIDEAL ;
       baseInfo::MINIDEAL=baseInfo::MINIDEAL > qoitmp.IDEAL_NZERO ?
                                           qoitmp.IDEAL_NZERO : baseInfo::MINIDEAL ;

       /* get optimization variables for each group 
          NOTE: if more than one group tries to optimize the same variable
                this will cause a race condition when collecting the
                optimized variable                                           */

       std::ostringstream qoiID ;
       qoiID << "qoi_"<< iii << "/";
       controlfile.set_prefix(qoiID.str());
       // thermal conductivities
       if( controlfile(  "optimize_k_0", false) ){
          qoitmp.NumOptVar++; qoitmp.optimize_K_0 = PETSC_TRUE; numopt_K_0++;
       }
       if( controlfile(  "optimize_k_1", false) ){
          qoitmp.NumOptVar++; qoitmp.optimize_K_1 = PETSC_TRUE; numopt_K_1++;
       }
       if( controlfile(  "optimize_k_2", false) ){
          qoitmp.NumOptVar++; qoitmp.optimize_K_2 = PETSC_TRUE; numopt_K_2++;
       }
       if( controlfile(  "optimize_k_3", false) ){
          qoitmp.NumOptVar++; qoitmp.optimize_K_3 = PETSC_TRUE; numopt_K_3++;
       }

       // blood perfusion
       if( controlfile(  "optimize_w_0", false) ){
          qoitmp.NumOptVar++; qoitmp.optimize_W_0 = PETSC_TRUE; numopt_W_0++;
       }
       if( controlfile(  "optimize_w_n", false) ){
          qoitmp.NumOptVar++; qoitmp.optimize_W_N = PETSC_TRUE; numopt_W_N++;
       }
       if( controlfile(  "optimize_w_i", false) ){
          qoitmp.NumOptVar++; qoitmp.optimize_W_I = PETSC_TRUE; numopt_W_I++;
       }
       if( controlfile(  "optimize_w_d", false) ){
          qoitmp.NumOptVar++; qoitmp.optimize_W_D = PETSC_TRUE; numopt_W_D++;
       }
       if( controlfile(  "optimize_w_2", false) ){
          qoitmp.NumOptVar++; qoitmp.optimize_W_2 = PETSC_TRUE; numopt_W_2++;
       }
       if( controlfile(  "optimize_w_ni", false) ){
          qoitmp.NumOptVar++; qoitmp.optimize_W_NI = PETSC_TRUE; numopt_W_NI++;
       }
       if( controlfile(  "optimize_w_id", false) ){
          qoitmp.NumOptVar++; qoitmp.optimize_W_ID = PETSC_TRUE; numopt_W_ID++;
       }

       // laser parameters
       if( controlfile(  "optimize_x_0" , false) ){
          qoitmp.NumOptVar++; qoitmp.optimize_X_0 = PETSC_TRUE; numopt_X_0++;
       }
       if( controlfile(  "optimize_y_0" , false) ){
          qoitmp.NumOptVar++; qoitmp.optimize_Y_0 = PETSC_TRUE; numopt_Y_0++;
       }
       if( controlfile(  "optimize_z_0" , false) ){
          qoitmp.NumOptVar++; qoitmp.optimize_Z_0 = PETSC_TRUE; numopt_Z_0++;
       }
       if( controlfile(  "optimize_pow" , false) ){
          qoitmp.NumOptVar++; qoitmp.optimize_POW = PETSC_TRUE; numopt_POW++;
       }
       if( controlfile(  "optimize_mu_a", false) ){
          qoitmp.NumOptVar++; qoitmp.optimize_MU_A= PETSC_TRUE; numopt_MU_A++;
       }
       if( controlfile(  "optimize_mu_s", false) ){
          qoitmp.NumOptVar++; qoitmp.optimize_MU_S= PETSC_TRUE; numopt_MU_S++;
       }
       /*NOTE: b/c we are ALWAYS computing error estimates the logic of
               the mesh optimization logic is slightly different  we DO NOT
               increment the qoitmp.NumOptVar counter for this variable*/
       if( controlfile(  "optimize_mesh", false) ) 
                                            GenInfo.optimize_MESH= PETSC_TRUE;

       qoitmp.ntask = tasks_per_group[iii];
       //compute total # of optimization solves that need be done.
       GenInfo.Total_Nopt += qoitmp.NoptSteps;
       QOIinfo.push_back(qoitmp);
    }
  controlfile.set_prefix(""); // reset back to the default prefix

  // error check to ensure that each optimization varirable is optimized by
  // only ONE group of processors
  if( numopt_K_0  > 1 ||  numopt_W_0  > 1 ||
      numopt_K_1  > 1 ||  numopt_W_N  > 1 ||
      numopt_W_I  > 1 ||  numopt_W_D  > 1 ||
      numopt_K_2  > 1 ||  numopt_W_2  > 1 ||
      numopt_K_3  > 1 ||  numopt_W_NI > 1 || numopt_W_ID > 1 ||
      numopt_X_0  > 1 ||  numopt_MU_A > 1 ||
      numopt_Y_0  > 1 ||  numopt_MU_S > 1 ||
      numopt_Z_0  > 1 ||  numopt_POW  > 1 ) { 
      std::cout << "A variable cannot be optimized by more than one group \n"; 
      std::cout << "Race Condition created \n" ; abort(); 
  }

  // initialize booleans to read in monte carlo grid
  PetscTruth MCData=PETSC_FALSE, MCIdl=PETSC_FALSE;

  for (unsigned int iii = 0 ; iii < tasks_per_group.size() ; iii++) 
   {
   if(GenInfo.Control_Task)
     { // echo parameters
       std::cout<<"\ndddas:[qoi_"<<iii<<"]ntask       = " << QOIinfo.at(iii).ntask ;
       std::cout<<"\ndddas:[qoi_"<<iii<<"]NumOptVar   = " << QOIinfo.at(iii).NumOptVar ;
       std::cout<<"\ndddas:[qoi_"<<iii<<"]OPTIMIZE_K_0= " << QOIinfo.at(iii).optimize_K_0 ;
       std::cout<<"\ndddas:[qoi_"<<iii<<"]OPTIMIZE_K_1= " << QOIinfo.at(iii).optimize_K_1 ;
       std::cout<<"\ndddas:[qoi_"<<iii<<"]OPTIMIZE_K_2= " << QOIinfo.at(iii).optimize_K_2 ;
       std::cout<<"\ndddas:[qoi_"<<iii<<"]OPTIMIZE_K_3= " << QOIinfo.at(iii).optimize_K_3 ;
       std::cout<<"\ndddas:[qoi_"<<iii<<"]OPTIMIZE_W_0= " << QOIinfo.at(iii).optimize_W_0 ;
       std::cout<<"\ndddas:[qoi_"<<iii<<"]OPTIMIZE_W_N= " << QOIinfo.at(iii).optimize_W_N ;
       std::cout<<"\ndddas:[qoi_"<<iii<<"]OPTIMIZE_W_I= " << QOIinfo.at(iii).optimize_W_I ;
       std::cout<<"\ndddas:[qoi_"<<iii<<"]OPTIMIZE_W_D= " << QOIinfo.at(iii).optimize_W_D ;
       std::cout<<"\ndddas:[qoi_"<<iii<<"]OPTIMIZE_W_2= " << QOIinfo.at(iii).optimize_W_2 ;
       std::cout<<"\ndddas:[qoi_"<<iii<<"]OPTIMIZE_W_NI=" << QOIinfo.at(iii).optimize_W_NI;
       std::cout<<"\ndddas:[qoi_"<<iii<<"]OPTIMIZE_W_ID=" << QOIinfo.at(iii).optimize_W_ID;
       std::cout<<"\ndddas:[qoi_"<<iii<<"]OPTIMIZE_X_0= " << QOIinfo.at(iii).optimize_X_0 ;
       std::cout<<"\ndddas:[qoi_"<<iii<<"]OPTIMIZE_Y_0= " << QOIinfo.at(iii).optimize_Y_0 ;
       std::cout<<"\ndddas:[qoi_"<<iii<<"]OPTIMIZE_Z_0= " << QOIinfo.at(iii).optimize_Z_0 ;
       std::cout<<"\ndddas:[qoi_"<<iii<<"]OPTIMIZE_MU_A=" << QOIinfo.at(iii).optimize_MU_A;
       std::cout<<"\ndddas:[qoi_"<<iii<<"]OPTIMIZE_MU_S=" << QOIinfo.at(iii).optimize_MU_S;
       std::cout<<"\ndddas:[qoi_"<<iii<<"]OPTIMIZE_POW= " << QOIinfo.at(iii).optimize_POW ;
       std::cout<<"\ndddas:[qoi_"<<iii<<"]NoptSteps   = " << QOIinfo.at(iii).NoptSteps ;
       std::cout<<"\ndddas:[qoi_"<<iii<<"]NUMIDEAL    = " << QOIinfo.at(iii).NUMIDEAL ;
       std::cout<<"\ndddas:[qoi_"<<iii<<"]NOFFSET     = " << QOIinfo.at(iii).NOFFSET ;
       std::cout<<"\ndddas:[qoi_"<<iii<<"]IDEAL_NZERO = " << QOIinfo.at(iii).IDEAL_NZERO ;
       std::cout<<"\ndddas:[qoi_"<<iii<<"]IDEAL_NTIME = " << QOIinfo.at(iii).IDEAL_NTIME ;
       std::cout<<"\ndddas:[qoi_"<<iii<<"]TIMELAG     = " << QOIinfo.at(iii).TIMELAG ;
       std::cout<<"\ndddas:[qoi_"<<iii<<"]PLOTOPTSOLVE= " << QOIinfo.at(iii).PLOTOPTSOLVE;
       std::cout<<"\ndddas:[qoi_"<<iii<<"]using pde     " << QOIinfo.at(iii).pde;
       std::cout<<"\ndddas:[qoi_"<<iii<<"]using obj fnc " << QOIinfo.at(iii).compobj;
       std::cout<<"\ndddas:[qoi_"<<iii<<"]FinalOpt_Ntime="<< QOIinfo.at(iii).FinalOpt_Ntime;
       std::cout<<"\ndddas:[qoi_"<<iii<<"]FinalOpt_Nzero="<< QOIinfo.at(iii).FinalOpt_Nzero;
     }
   // determine if need to read in MC data
   if(QOIinfo.at(iii).pde.find("nonlinpennesmonte")  !=std::string::npos) MCData=PETSC_TRUE;
   if(QOIinfo.at(iii).compobj.find("computeideal_mc")!=std::string::npos) MCIdl=PETSC_TRUE;
  }

  //error check make sure that read in the lower # thermal image that is needed
  if( baseInfo::MRTI_nzero > baseInfo::MINIDEAL ) 
                baseInfo::MRTI_nzero = baseInfo::MINIDEAL;

  if(GenInfo.Control_Task)
   { // echo parameters
     GenInfo.printSelf();
     baseInfo::printSelf();
   }

  /* the sequence of initialization is IMPORTANT!!!!!!!!!!!!!
        initialization of the visualase data is used in dddas::Image  class  */
  FORTRAN_NAME(initialize_visualase)(&tao_monitor,&GenInfo.Control_Task,
                                     &GenInfo.optimize_MESH,&groupID);
  PetscMPIInt dddascommpointer = PetscFromPointerComm(GenInfo.DDDAS_COMM_WORLD);
  FORTRAN_NAME(init_plot_params)(&dddascommpointer);
  if(GenInfo.Control_Task) printf("@@@@@@@@@@@@@@@@@Fortran Data Structures Initialized\n");

  //// initialize Monte Carlo data if necessary
  // 
  //if(MCData ){info=Load_MCdataNopt(GenInfo.Control_Task);  CHKERRQ(info); }
  //if(MCIdl  ){info=Load_MCdataIdeal(GenInfo.Control_Task); CHKERRQ(info); }
  //if(MCData || MCIdl){info=Init_MCdata(GenInfo.Control_Task); CHKERRQ(info); }

  //CHKMEMQ; // check for memory corruption use -malloc_debug to enable

  //char vistransfer_template[MAXLEN], vistransfer[MAXLEN], transferfile[32]; 

  // setup the mesh
  PetscLogStagePush(logstages[0]); //Begin Mesh Generation Stage

  // setup all model parameters and setup the 
  // mesh to handle spatially varying material parameters 
  info = SetupModel(controlfile,GenInfo,QOIinfo,mesh);

  PetscLogStagePop(); // End Mesh generation Stage

  std::ostringstream profileFile;// filename for profile
  ///* Specific Processors have specific jobs */ 
  if ( (unsigned int) groupID == tasks_per_group.size() ) { 
  //   // write file to HOME directory to find the name of the Fem Write Node
  //   for(iii = 0 ; iii < number_comp_group ; iii++ ){
  //      sprintf(vistransfer_template,"%s",
  //                 ini.GetValue("compexec","vistransferfem","",&bHasMulti));
  //      sprintf(vistransfer , vistransfer_template      ,iii);
  //      sprintf(transferfile,"datatransferfemqoi_%d.txt",iii);
  //      info = TheRedPill(vistransfer,transferfile); 
  //   }
  //   // write file to HOME directory to find the name of the Data Server
  //   sprintf(vistransfer,"%s",
  //             ini.GetValue("compexec", "vistransfermri","",&bHasMulti));
  //   sprintf(transferfile,"%s","datatransfermri.txt");
  //   info = TheRedPill(vistransfer,transferfile); 

     // main code for control
     info = DDDAS_Control(leadcomp,GenInfo,QOIinfo,mesh);
     /* Print Profile */
     profileFile << "control.log" << GenInfo.profileID;
  } else {  

     info = Compute_Drone(argc,argv,controlfile,GenInfo,
                          QOIinfo.at(groupID),mesh,groupID);

     /* Print Profile */
     profileFile << "comp"<< groupID << ".log" << GenInfo.profileID;
  } 

  PetscLogPrintSummary(PETSC_COMM_WORLD,profileFile.str().c_str()); 
  
  //  libMesh::close(); // PetscFinalize();  called by libMesh::close();
  return 0; //  
}
///* ------------------------------------------------------------------- 
//   TheRedPill - write a file to $HOME to find the Data Server
//*/
//#undef __FUNCT__
//#define __FUNCT__ "TheRedPill"
//PetscErrorCode TheRedPill(const char *XferCmd,const char *filename){
//  char file_o[MAXLEN],DataXferCmd[MAXLEN];
//  ofstream xferfile;
//
//  PetscFunctionBegin; 
//
//  // write file in home dir
//  sprintf(file_o,"%s/%s",getenv("HOME"),filename); 
//  // ssh command
//  sprintf(DataXferCmd,"ssh `echo %s | cut -d ':' -f7 | cut -d ' ' -f1` %s",
//                                           getenv("SSH_CONNECTION"),XferCmd);
//  xferfile.open(file_o, ios::out);
//  xferfile << DataXferCmd << std::endl;
//  xferfile.close();
//  // make the file executable
//  chmod(file_o,S_IREAD|S_IWRITE|S_IEXEC);
//  PetscFunctionReturn(0);
//}
//
// 
