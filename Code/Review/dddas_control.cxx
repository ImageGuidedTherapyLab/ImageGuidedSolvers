// C++ include files 
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

// Petsc/Tao include files
#include "tao.h" // petsc solvers w/ Fortran Interface

// libMesh include files
#include "libmesh.h"
#include "mesh.h"
#include "equation_systems.h"
#include "exact_solution.h"
#include "getpot.h"
#include "parallel.h"
#include "mesh_base.h"
#include "dof_map.h"
#include "elem.h"
#include "fe_interface.h"
// The nonlinear solver and system we will be using
#include "nonlinear_implicit_system.h"
#include "linear_implicit_system.h"

// dddas include files
#include "transient_inverse_system.h"
#include "fortrandftranslate.h" // header to call fortran routines
#include "parser_defaults.h" // ensure same default values between C and Fortran
#include "mri_params.h" // mri params
#include "solvehp3d.h" // main header


#include "mri_params.h" // default values for mri data
#include "variable_map.h" // default values for mri data
#include "read_mridata.h" // imaging data structures


/* --------------------- User-defined routines ------------------------ */
#include "dddas.h"
#include "baseInfo.h"
#include "pennesVerification.h"
extern FORTRAN_FUNCTION 
{
#include "global_params.h"
#include "pennes_model.h"
}




/* -------------------------------------------------------------------- */

#undef __FUNCT__
#define __FUNCT__ "DDDAS_Control"
PetscErrorCode DDDAS_Control(std::vector<PetscInt> &LeadComp,UniversalCntrl &GenInfo,
                             std::vector<DroneControl> &QOIInfo,          Mesh &mesh)
{
  PetscInt nzero=0; //global to the file
  PetscErrorCode info; /* used to check for functions returning nonzeros */
  
  // used for error handling
  PetscFunctionBegin; 

  CHKMEMQ; // check for memory corruption use -malloc_debug to enable

  /* Establish communicators with all groups  
     data client-server communicator
     the root task is the lead in all groups */ 
  std::vector<MPI_Comm> ControlComm(LeadComp.size(),MPI_COMM_NULL), 
                      DataComm(LeadComp.size(),MPI_COMM_NULL); 
  for (unsigned int iii = 0 ; iii < LeadComp.size() ; iii++) {
       MPI_Intercomm_create(PETSC_COMM_WORLD, 0, GenInfo.DDDAS_COMM_WORLD, 
                                   LeadComp.at(iii), iii , &DataComm.at(iii) );
       MPI_Intercomm_merge(DataComm.at(iii),0, &ControlComm.at(iii) );
  }

  // allocate memory for visualization data structures
  //FORTRAN_NAME(visallocate)();
  PetscPrintf(PETSC_COMM_SELF,
              "@@@@@@@@@@@@@@@@@Control Task vis data structure allocated\n");

//  // only checking domain decomposition
//  if(GenInfo.NEXACT < 0){
//     groupID=0;
//     // Setup the Hp3d Data structures
//     info = SetupHp3d(GenInfo,ControlComm.at(groupID),
//                        &QOIInfo.at(groupID).recvcnts_elm,
//                        &QOIInfo.at(groupID).recvcnts_for,
//                        &QOIInfo.at(groupID).recvcnts_adj,
//                        &QOIInfo.at(groupID).displs_elm,
//                        &QOIInfo.at(groupID).displs_for,
//                        &QOIInfo.at(groupID).displs_adj); 
//     /* allocate buffer for inter group field transmission */
//     nreleb=FORTRAN_NAME(get_nreleb)();
//     FORTRAN_NAME(alloc_hpelemfield)(&nreleb);
//     /* gather perfusivity and conductivity fields */
//     PetscMPIInt fortran_comm=PetscFromPointerComm(
//                                               ControlComm.at(groupID));
//     FORTRAN_NAME(gatherdealloc_hpelemfield)(&fortran_comm, 
//                                   &QOIInfo.at(groupID).optimize_W_0,
//                                   &QOIInfo.at(groupID).optimize_K_0,
//                                   &QOIInfo.at(groupID).ntask,
//                                   &QOIInfo.at(groupID).recvcnts_elm[0],
//                                   &QOIInfo.at(groupID).displs_elm[0] );
//     //create data structure to collect temp field for plotting
//     info=SetupVisualization(&groupID,&nzero,
//                    &(QOIInfo.at(groupID)),&nzero );
//     PetscFunctionReturn(0);
//  }


  // initialize image server
  // Sequence is IMPORTANT!!!! Setup Power uses data setup in imaging
  dddas::ImageServer MainImageServer(GenInfo,QOIInfo);
  info = SetupPower(GenInfo);

  // check if thers is a calibration problem
  PetscTruth setupimage=PETSC_FALSE; 
  for (unsigned int iii = 0 ; iii < QOIInfo.size() ; iii++) 
     if(QOIInfo[iii].compobj.find("calibration")!=std::string::npos) setupimage=PETSC_TRUE;
  // setup imaging data structures
  if(setupimage || GenInfo.IC_MRTI) 
   {
     //wait for dicom header info to be transferred/available before proceeding
     MainImageServer.ServerSetup(GenInfo,QOIInfo);

     // setup image server buffers for writing
     MainImageServer.ImageSetupFile();
   }
  else
   {
     MainImageServer.NoWrite();
   }

  // to pass to fortran
  PetscScalar spacing[3]; MainImageServer.GetSpacing(spacing);

  // initialize parameters on server before writing
  PetscInt NegOne = -1; 
  FORTRAN_NAME(read_control_fortran)(&NegOne,spacing);

  //write initial ini, power, and field files to get control system started
  if(GenInfo.IDrestart ==0){// if this is an initial run write out the initial
                             // control file to get the process started        
     WriteControlFile( nzero , QOIInfo.size() );
  }
  FORTRAN_NAME(write_w_0field)(&nzero);
  FORTRAN_NAME(write_k_0field)(&nzero);
  // need an intial power visualization file at MRTI_nzero the power file 
  // written in write_visualase_file starts at MRTI_nzero+1
  FORTRAN_NAME(vis_power)(&baseInfo::MRTI_nzero);

  /* initialize counters and declare control variables */
  PetscInt inifileID =GenInfo.IDrestart,   // ini file ID
           inifileout=GenInfo.IDrestart+1, // id of next file to be written out
           IDimgfile =baseInfo::MRTI_nzero,   // image file ID
           IDpow    =IDimgfile+1,           // visualase file ID
           commands[2] , // buffer to send commands to comp groups
           icntopts=0,// optimization counter
           recv_complete=0,send_complete=0,
           istep,idplot,nsteplo,nstephi,
           TaskID; // TASKID == SETUP =>     this message sent before an
                   //                        optimization solve, setup the 
                   //                        data structures
                   // TASKID == CHKPOINT =>  this message sent during the 
                   //                        optimization solve an needs 
                   //                        to write a restart file
                   // TASKID == COMPLETE =>  this message sent after an
                   //                        optimization solve , update 
                   //                        the optimization counter and
                   //                        plot if necessary
                  
  // counters to keep track of which optimization step each group is on
  std::vector<PetscInt> noptcnt(LeadComp.size(),nzero), 
                   npltcnt(LeadComp.size(),nzero); 

  // get maximum visualase temperature & write initial visualase power file 
  PetscScalar MaxTemp = FORTRAN_NAME(get_visualase_max_temp)();
  FORTRAN_NAME(write_visualase_file)(&IDpow,&MaxTemp);

  CHKMEMQ; // check for memory corruption use -malloc_debug to enable

  // setup MPI data structures and begin initial recieve
  MPI_Status  status;  
  std::vector<MPI_Request> Send_Request;
  std::vector<MPI_Request>::iterator SendIter;
  MPI_Request recv_request;
  MPI_Irecv(commands,2,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,
                        GenInfo.DDDAS_COMM_WORLD,&recv_request);

  /*all data structures have been setup. now orchestrate
    the multiple optimization groups and reading/sending the 
    thermal imaging data */
  while(icntopts < GenInfo.Total_Nopt || IDimgfile<=baseInfo::MRTI_ntime
                                      ||    Send_Request.size() != 0) 
  {
     /*Notice that MRTIdata is initialized by default to the initial
       temperature value of the domain. If Thermal images not read in, 
       i.e. mrtiavs .eq. PETSC_FALSE but data still sent to computational
       groups then, the groups will just get an image with a constant value.
       read in all data processing all requests if needed */
     if(IDimgfile<=baseInfo::MRTI_ntime) 
      MainImageServer.LaserImageControl(IDimgfile,&GenInfo,&QOIInfo,&DataComm,&Send_Request);
      
     // destroy the list of send requests as they are recieved must do 
     //  an MPI_Test before the corresponding MPI_Recv will post and unblock
     for(SendIter =Send_Request.begin();
         SendIter!=Send_Request.end()  ; SendIter++)
     {
         send_complete=0;
         MPI_Test(&SendIter[0],&send_complete,&status);
         // if the data was successfully send break the loop to avoid 
         // any problems changing the loop bounds
         if(send_complete){ 
             Send_Request.erase(SendIter); break; 
         }
     }
     // main control of optimization problems
     if(icntopts < GenInfo.Total_Nopt)
     {
        MPI_Test(&recv_request,&recv_complete,&status);
        if(recv_complete)
        {
           TaskID  = commands[1];
           // the group's ID is sent as the tag
           PetscInt groupID=status.MPI_TAG;
  
           switch(TaskID){ 
             // need to setup data structures
             case SETUP:

               // OVERWRITE power w/ the most up to data power field 
               FORTRAN_NAME(write_power)(&inifileID);
 
               /* signal drone to execute computations from 
                          files/control%d.ini % inifileID
                           and write the next restart file to 
                          files/control%d.ini % inifileout */
               commands[0]=inifileID;
               commands[1]=inifileout; 
               inifileout++; //update file counter for next group
               MPI_Bcast(&commands,2,MPI_INT,0,ControlComm.at(groupID));
               // Setup the Hp3d Data structures
               //info = SetupHp3d(GenInfo,ControlComm.at(groupID),
               //                   &QOIInfo.at(groupID).recvcnts_elm,
               //                   &QOIInfo.at(groupID).recvcnts_for,
               //                   &QOIInfo.at(groupID).recvcnts_adj,
               //                   &QOIInfo.at(groupID).displs_elm,
               //                   &QOIInfo.at(groupID).displs_for,
               //                   &QOIInfo.at(groupID).displs_adj); 
               printf("@@@@@@@@@@@@@@Control Task Setup Group %d\n",groupID);
               break;
             default:
               // groupID has checkpointed/completed an optimization solve 
               //    and just wrote file files/control%d.ini % inifileID 
               inifileID  = commands[0];
               /* set the time window that was just optimized for
                  this is used in plotting */ 
               nsteplo= baseInfo::getFEMTimeStep( 
                            QOIInfo.at(groupID).IDEAL_NZERO );
               nstephi= baseInfo::getFEMTimeStep( 
                            QOIInfo.at(groupID).IDEAL_NTIME );
               PetscScalar Tau = baseInfo::getTime( nstephi);
               FORTRAN_NAME(setfemwindow)(&QOIInfo.at(groupID).IDEAL_NZERO,
                                          &QOIInfo.at(groupID).IDEAL_NTIME,
                                          &nsteplo,&nstephi,&Tau);
  
               /* allocate buffer for inter group field transmission */
               //nreleb=FORTRAN_NAME(get_nreleb)();
               //FORTRAN_NAME(alloc_hpelemfield)(&nreleb);
               /* gather perfusivity and conductivity fields */
               PetscMPIInt fortran_comm=PetscFromPointerComm(
                                                    ControlComm.at(groupID));
               //FORTRAN_NAME(gatherdealloc_hpelemfield)(&fortran_comm, 
               //                       &QOIInfo.at(groupID).optimize_W_0,
               //                       &QOIInfo.at(groupID).optimize_K_0,
               //                       &QOIInfo.at(groupID).ntask,
               //                       &QOIInfo.at(groupID).recvcnts_elm[0],
               //                       &QOIInfo.at(groupID).displs_elm[0] );
  
               /* write out field files for next group */
               FORTRAN_NAME(write_w_0field)(&inifileID);
               FORTRAN_NAME(write_k_0field)(&inifileID);
  
               //write out error field to drive adaptivity
               //FORTRAN_NAME(write_errfield)(&inifileID); 
               {//create scope so destructor called before file overwritten
                  // open the Ini file just written by groupID
                  std::ostringstream filename;
                  filename<< "files/control"<<inifileID<<".ini";
                  GetPot LatestIni(filename.str()); 
  
                  /* with multiple groups, must read only from the file 
                     written by the group optimizing the power. otherwise 
                     the possibility of overwritting exists. 
                     i.e. group 0 optimizes power group and group 1 DOES 
                     NOT both groups start computing at the same time
                     group 0 finishes first writes out the power
                     group 1 finishes and writes out the power and 
                     reading these in will be the latest ones in memory.
                     also have to check that the id of the
                     time instance is within the bounds of 
                     the power data structure.
                     to avoid overwriting previous history the power file
                     is read in from the current step of the power control
                     which is being up dated through shared memory by the
                     data thread */
                  IDpow =IDimgfile+1;   
                  if(QOIInfo.at(groupID).optimize_POW
                        &&  IDpow<=baseInfo::MRTI_ntime ){
                     FORTRAN_NAME(read_power)(&inifileID,&IDpow);
                  } 
                  // get latest solution to optimization problem
                  if(QOIInfo.at(groupID).optimize_W_N){
                     PetscScalar w_n = LatestIni("perfusivity/w_n",0.0) ;
                     FORTRAN_NAME(putparam_w_n)(&w_n,&nzero);
                  }
                  if(QOIInfo.at(groupID).optimize_W_I){
                     PetscScalar w_i = LatestIni("perfusivity/w_i",0.0) ;
                     FORTRAN_NAME(putparam_w_i)(&w_i,&nzero);
                  }
                  if(QOIInfo.at(groupID).optimize_W_D){
                     PetscScalar w_d = LatestIni("perfusivity/w_d",0.0) ;
                     FORTRAN_NAME(putparam_w_d)(&w_d,&nzero);
                  }
                  if(QOIInfo.at(groupID).optimize_W_2){
                     PetscScalar w_2 = LatestIni("perfusivity/w_2",0.0) ;
                     FORTRAN_NAME(putparam_w_2)(&w_2,&nzero);
                  }
                  if(QOIInfo.at(groupID).optimize_W_NI){
                     PetscScalar w_ni = LatestIni("perfusivity/w_ni",0.0) ;
                     FORTRAN_NAME(putparam_w_ni)(&w_ni,&nzero);
                  }
                  if(QOIInfo.at(groupID).optimize_W_ID){
                     PetscScalar w_id = LatestIni("perfusivity/w_id",0.0) ;
                     FORTRAN_NAME(putparam_w_id)(&w_id,&nzero);
                  }
                  if(QOIInfo.at(groupID).optimize_K_1){
                     PetscScalar k_1 = LatestIni("thermalcond/k_1",0.0) ;
                     FORTRAN_NAME(putparam_k_1)(&k_1,&nzero);
                  }
                  if(QOIInfo.at(groupID).optimize_K_2){
                     PetscScalar k_2 = LatestIni("thermalcond/k_2",0.0) ;
                     FORTRAN_NAME(putparam_k_2)(&k_2,&nzero);
                  }
                  if(QOIInfo.at(groupID).optimize_K_3){
                     PetscScalar k_3 = LatestIni("thermalcond/k_3",0.0) ;
                     FORTRAN_NAME(putparam_k_3)(&k_3,&nzero);
                  }
                  if(QOIInfo.at(groupID).optimize_X_0){
                     PetscScalar x_0 = LatestIni("source_laser/x_0",0.0) ;
                     FORTRAN_NAME(putparam_x_0)(&x_0,&nzero);
                  }
                  if(QOIInfo.at(groupID).optimize_Y_0){
                     PetscScalar y_0 = LatestIni("source_laser/y_0",0.0) ;
                     FORTRAN_NAME(putparam_y_0)(&y_0,&nzero);
                  }
                  if(QOIInfo.at(groupID).optimize_Z_0){
                     PetscScalar z_0 = LatestIni("source_laser/z_0",0.0) ;
                     FORTRAN_NAME(putparam_z_0)(&z_0,&nzero);
                  }
                  if(QOIInfo.at(groupID).optimize_MU_A){
                     PetscScalar mu_a = LatestIni("source_laser/mu_a",0.0) ;
                     FORTRAN_NAME(putparam_mu_a)(&mu_a,&nzero);
                  }
                  if(QOIInfo.at(groupID).optimize_MU_S){
                     PetscScalar mu_s = LatestIni("source_laser/mu_s",0.0) ;
                     FORTRAN_NAME(putparam_mu_s)(&mu_s,&nzero);
                  }
               } // ending the scope closes the ini file
  
               // over LatestIni with MainIni file to record all 
               // changes amongst computational groups
               WriteControlFile( inifileID, QOIInfo.size() );
  
               // see if need to plot
               if( (TaskID == COMPLETE && QOIInfo.at(groupID).PLOTOPTSOLVE)
                                             || GenInfo.PLOTITER ){
                  Vec var_for,var_adj;
                  // create local vector for forward problem 
                  //                         and adjoint solutions
                  info=VecCreateSeq(PETSC_COMM_SELF,GenInfo.NDOF_FORWARD[1],
                                    &var_for); 
                  info=VecCreateSeq(PETSC_COMM_SELF,GenInfo.NDOF_ADJOINT[1],
                                    &var_adj);
                  //create data structure to collect temp field for plotting
                  //info=SetupVisualization(&groupID,&noptcnt.at(groupID),
                  //            &(QOIInfo.at(groupID)),&npltcnt.at(groupID) );
                  printf("@@@@@@@@@@@@@@@@@Control Task setup visualization for group %d\n", groupID );
                  MPI_Irecv(commands,2,MPI_INT,LeadComp.at(groupID),
                            VIS_SOLVE_COMPLETE,GenInfo.DDDAS_COMM_WORLD,
                                                           &recv_request);
                  recv_complete=0;
                  while(!recv_complete){
                     // while waiting for computation to finish, 
                     // check if new image is ready and update if necessary
                     if(IDimgfile<=baseInfo::MRTI_ntime) 
                        MainImageServer.LaserImageControl(IDimgfile,&GenInfo,
                                                          &QOIInfo,&DataComm,&Send_Request);
                     MPI_Test(&recv_request,&recv_complete,&status);
                  }
                  for(idplot = QOIInfo.at(groupID).IDEAL_NZERO; 
                                     idplot <= baseInfo::MAXIDEAL ; idplot++){
                      istep = baseInfo::getFEMTimeStep( idplot );
                      //Main_Visualization(ControlComm.at(groupID),NULL,
                      //     &groupID,&noptcnt.at(groupID),&idplot,&istep,
                      //     GenInfo,&(QOIInfo.at(groupID)), var_for,var_adj);
                      // check if new image is ready and update if necessary
                      if(IDimgfile<=baseInfo::MRTI_ntime) 
                        MainImageServer.LaserImageControl(IDimgfile,&GenInfo,
                                                          &QOIInfo,&DataComm,&Send_Request);
                  }
                  info = VecDestroy(var_for); 
                  info = VecDestroy(var_adj); 
                  printf("@@@@@@@@@@@@@@@@@Control Task Completed Final Visualization for Group %d \n",groupID );
                  npltcnt.at(groupID)++; // update plot counter
               }
  
               switch(TaskID){ 
                  // this was a checkpoint 
                  case CHKPOINT:
                     MPI_Bcast(&inifileout,1,MPI_INT,0,
                                             ControlComm.at(groupID));
                     inifileout++; // update file counter for next group
                  break;
                  // the optimization solve completed 
                  case COMPLETE: 
                     QOIInfo.at(groupID).IDEAL_NZERO = 
                     QOIInfo.at(groupID).IDEAL_NZERO 
                            + QOIInfo.at(groupID).NOFFSET ;
                     QOIInfo.at(groupID).IDEAL_NTIME = 
                     QOIInfo.at(groupID).IDEAL_NTIME 
                            + QOIInfo.at(groupID).NUMIDEAL;
                     printf("@@@@@@@@@@@@@@@@@Control Task Completed Optimization Step %d for Group %d\n",noptcnt.at(groupID),groupID );
                     noptcnt.at(groupID)++; // update group counter
                     icntopts++; // update global optimization step counter
                  break;
               } 
           } 
           MPI_Irecv(commands,2,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,
                                 GenInfo.DDDAS_COMM_WORLD,&recv_request);
        }
     } // end main control of optimization problems
  }
//
//  // the computation has completed, save a final control file
//  //  this is really only used in parameters studies to have a single
//  //  control file with the same id to look at.
//  MainIni.SaveFile("files/control_final.ini");

  PetscFunctionReturn(0);
}


