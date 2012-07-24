//C++ include files that we need
#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>

// Basic include file needed for the mesh functionality.
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "equation_systems.h"
#include "exact_solution.h"
#include "dense_vector.h"
// The nonlinear solver and system we will be using
#include "linear_implicit_system.h"
#include "nonlinear_implicit_system.h"
#include "nonlinear_solver.h"

// The definition of a geometric element
#include "elem.h"

// file parser
#include "getpot.h"

//tao
#include "tao.h" 

//dddas includes
#include "fortrandftranslate.h" // header to call fortran routines
#include "mri_params.h" // mri data types
//#include "elem_matvec.h" // element variables
#include "variable_map.h" 

// dddas include files
#include "transient_inverse_system.h"
#include "baseInfo.h" 
#include "solvehp3d.h" 
#include "pennesVerification.h"

/* -------------- User-defined routines ---------- */
//PetscErrorCode FormFunction(SNES,Vec,Vec,void*);
//PetscErrorCode FormJacobian(SNES,Vec,Mat*,Mat*,MatStructure*,void*);
extern FORTRAN_FUNCTION
{
#include "pennes_model.h" // interface to pennes model
//  void FORTRAN_NAME(meshcheck)();
//  void FORTRAN_NAME(hp3gendf)(PetscInt *);
//  void FORTRAN_NAME(setup_fields)
//  void FORTRAN_NAME(set_sbs_flags)();
//  void FORTRAN_NAME(forward2adjoint)();
//  void FORTRAN_NAME(setupserial)(PetscInt *,PetscInt *);
//  void FORTRAN_NAME(dealloc_vis)();
//  void FORTRAN_NAME(sortandcondense)();
//  void FORTRAN_NAME(preprocess_plot)(PetscInt*);
//  void FORTRAN_NAME(compdamagefrac)( PetscInt*,PetscInt*);
//  void FORTRAN_NAME(storefield)(PetscInt*,PetscScalar*,PetscInt*,PetscInt*,PetscInt*, void(*)( PetscInt*,PetscInt*,PetscScalar*), void(*)( PetscInt*,PetscInt*,PetscInt*,PetscScalar*));
//  void FORTRAN_NAME(get__field)(PetscInt*,PetscScalar*,PetscInt*,PetscInt*,PetscInt*, PetscScalar(*)( PetscInt*,PetscInt*), PetscScalar(*)( PetscInt*,PetscInt*,PetscInt*));
//  void FORTRAN_NAME(output_fem)(PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*, void (*)(PetscScalar*,PetscInt*,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*));
//  void FORTRAN_NAME(createlocalmap)(const PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*);
//  void FORTRAN_NAME(countmatentries)(const PetscInt*,PetscInt*,PetscInt*);
//  void FORTRAN_NAME(nelconb)(PetscInt *,PetscInt *);
//  void FORTRAN_NAME(computeglobalmaps)(PetscInt *,PetscInt *,PetscInt *);
//  void FORTRAN_NAME(setupgausspoints)(PetscInt*,PetscInt*);
//  void FORTRAN_NAME(getglobalmap)(PetscInt*,PetscInt*,PetscInt(*)(PetscInt *),PetscInt*,PetscInt*,PetscInt*);
//  void FORTRAN_NAME(elemmass)(const PetscInt *,PetscInt *,PetscInt *,PetscInt *,PetscScalar *);
//  void FORTRAN_NAME(read_w_0field)(PetscInt *);
//  void FORTRAN_NAME(read_k_0field)(PetscInt *);
//  void     FORTRAN_NAME(put_vert_global)(PetscInt *,PetscInt *,PetscInt *);
//  void     FORTRAN_NAME(put_node_global)(PetscInt *,PetscInt *,PetscInt *);
//  void     FORTRAN_NAME(put_vert_local)(PetscInt *,PetscInt *,PetscInt *);
//  void     FORTRAN_NAME(put_node_local)(PetscInt *,PetscInt *,PetscInt *);
//  PetscInt FORTRAN_NAME(get_nreleb)() ;
//  PetscInt FORTRAN_NAME(get_maxelib)();
//  PetscInt FORTRAN_NAME(get_order_state)(PetscInt*) ;
//  PetscInt FORTRAN_NAME(get_order_adjnt)(PetscInt*) ;
//  PetscInt FORTRAN_NAME(get_mdle_elem)(PetscInt*);
//  void FORTRAN_NAME(setup_numbering)(PetscInt*,PetscInt*,PetscInt*, void(*)(PetscInt *,PetscInt *,PetscInt *), void(*)(PetscInt *,PetscInt *,PetscInt *));
//  void FORTRAN_NAME(unvisit)();
//  PetscScalar FORTRAN_NAME(get_vert_temp)( PetscInt*,PetscInt*);
//  PetscScalar FORTRAN_NAME(get_node_temp)( PetscInt*,PetscInt*,PetscInt*);
//  void        FORTRAN_NAME(put_vert_temp)( PetscInt*,PetscInt*,PetscScalar*);
//  void        FORTRAN_NAME(put_node_temp)( PetscInt*,PetscInt*,PetscInt*,PetscScalar*);
//  PetscScalar FORTRAN_NAME(get_vert_ideal)(PetscInt*,PetscInt*);
//  PetscScalar FORTRAN_NAME(get_node_ideal)(PetscInt*,PetscInt*,PetscInt*);
//  void        FORTRAN_NAME(put_vert_ideal)(PetscInt*,PetscInt*,PetscScalar*);
//  void        FORTRAN_NAME(put_node_ideal)(PetscInt*,PetscInt*,PetscInt*,PetscScalar*);
//  PetscScalar FORTRAN_NAME(get_vert_cost)( PetscInt*,PetscInt*);
//  PetscScalar FORTRAN_NAME(get_node_cost)( PetscInt*,PetscInt*,PetscInt*);
//  void        FORTRAN_NAME(put_vert_cost)( PetscInt*,PetscInt*,PetscScalar*);
//  void        FORTRAN_NAME(put_node_cost)( PetscInt*,PetscInt*,PetscInt*,PetscScalar*);
//  PetscScalar FORTRAN_NAME(get_vert_dual)( PetscInt*,PetscInt*);
//  PetscScalar FORTRAN_NAME(get_node_dual)( PetscInt*,PetscInt*,PetscInt*);
//  void        FORTRAN_NAME(put_vert_dual)( PetscInt*,PetscInt*,PetscScalar*);
//  void        FORTRAN_NAME(put_node_dual)( PetscInt*,PetscInt*,PetscInt*,PetscScalar*);
//  void        FORTRAN_NAME(allocateshapscratch)();
//
}
// Global Identifiers
// control variable id
//static vector<MdleInfo>::iterator mdleID; // mdle node iterator
//static PetscScalar zero=0.0 ; 
//static PetscInt idstate=ISTATE , 
//                idadjnt=IADJNT , 
//                idinit =IDINIT ,
//                idflds =IDFLDS ,
//                id_slv =ID_SLV ,
//                nzero  = 0 ,
//                none   = 1 ;
/* ------------------------------------------------------------------- 
 SetupModel  - Setup data structures for the particular
                 PDE we are solving

   Preprocess Data structures to set up for spatially varying 
   material parameter inversion

          check the field distribution amongst processors there are 2 cases 
          (1) the parameter field is allowed to vary on every
              element on every processor and rank 0 owns NO dofs 
              in the parallel vec
          (2) there exists at least one element on one processor
              that belongs to the "constant" field portion of the
              parameter field. Consequently, rank 0 owns ONE dof 
              in the parallel vec

   NOTICE the fact that the entire mesh is stored on each processor
   is exploited. The loop over elements is NOT local to the processor

*/
#undef __FUNCT__
#define __FUNCT__ "SetupModel"
PetscErrorCode SetupModel(GetPot &controlfile, UniversalCntrl &GenInfo,
                              std::vector<DroneControl> &QOIinfo, Mesh &mesh)
{
  PetscFunctionBegin; 


  PetscLogEventBegin(logevents[3],0,0,0,0); // read disk

  // setup the fields for various pde constitutive equations
  PetscTruth setuppennes=PETSC_FALSE,  // setup pennes bio heat transfer
                setuptgm=PETSC_FALSE,  // setup tumor growth equations
              setupchabl=PETSC_FALSE;  // setup chemical ablation
  for (unsigned int iii = 0 ; iii < QOIinfo.size() ; iii++) 
    {
      if( QOIinfo.at(iii).pde.find("verifprob"           ) !=std::string::npos||
          QOIinfo.at(iii).pde.find("nonlinpennesmonte"   ) !=std::string::npos||
          QOIinfo.at(iii).pde.find("nonlinpennesrf"      ) !=std::string::npos||
          QOIinfo.at(iii).pde.find("pennesisolasercooling")!=std::string::npos||
          QOIinfo.at(iii).pde.find("nonlinpennesisolaser") !=std::string::npos
        ) setuppennes=PETSC_TRUE;
    }
  if(setuppennes)
    {
      PetscMPIInt dddascommpointer = PetscFromPointerComm(GenInfo.DDDAS_COMM_WORLD);
      /* setup data structures for varying perfusion and conductivity fields */
      FORTRAN_NAME(setup_pennes_model)(&dddascommpointer,&GenInfo.NfieldElem); 
  

      /* if this is a restart run read in the power as a function of time
         from the restart file overwriting the values that are currently there. 
         NOTE that init_power_dstruc setup the data structures so that when the 
         file is read in the data_structures are there ready to be overwritten */
     if(GenInfo.IDrestart !=0)
        {// this is a restart run
         PetscInt nzero=0;
         FORTRAN_NAME(read_power)(&GenInfo.IDrestart,&nzero);
         /*read in perfusivity and conductivity fields overwrite the values
           that are currently there. NOTE that setup fields setup the
           data structures so that when the file is read in the data_structures
           are there ready to be overwritten */
         //FORTRAN_NAME(read_w_0field)(&GenInfo.IDrestart);
         //FORTRAN_NAME(read_k_0field)(&GenInfo.IDrestart);
         
        }
    }

  PetscPrintf(GenInfo.DDDAS_COMM_WORLD,"@@@@@@@@@@@@@fields setup \n" );

  PetscLogEventEnd(logevents[3],0,0,0,0); // read disk

  CHKMEMQ; // check for memory corruption use -malloc_debug to enable
  /* check that the mesh lay within the IMAGE bounds
                          AND 
     dof's for  nodes and vertices on the processor boundary of a 
     parallel computation are not allocated in the above initial
     mesh generation must make a pass through the initial elements
     to check if all required dof's are allocated                       */
  PetscFunctionReturn(0);
}
/* ------------------------------------------------------------------- 
   SetupPower  - Setup Power data structures 
   the sequence of initialization is IMPORTANT!!!!!!!!!!!!!
     This must be called after the imaging data structures have been initialized
     All tasks initialize data structure that holds power as a function of time
        the control task uses this to write out the vis files
        the compute tasks uses this to get the power as a function of time
        the data server uses this to control the laser                       */
#undef __FUNCT__
#define __FUNCT__ "SetupPower"
PetscErrorCode SetupPower(UniversalCntrl &GenInfo)
{
  PetscFunctionBegin; 
  // all processors read in control parameters from ini file
  GetPot controlfile("files/control.ini");
  controlfile.set_prefix("");

  PetscScalar MaxTime = baseInfo::setupTime(controlfile);

  // compute groups read in control parameters
  if( !GenInfo.Control_Task ) FORTRAN_NAME(init_compute_drone)(&MaxTime);

  /* setup data structures for power*/
  PetscMPIInt dddascommpointer = PetscFromPointerComm(GenInfo.DDDAS_COMM_WORLD);
  FORTRAN_NAME(setup_power)(&dddascommpointer,&MaxTime); 

  if(GenInfo.Control_Task)
   { // echo parameters
     std::cout << "\ndddas: MaxTime      ="<<MaxTime  << std::endl;
     GenInfo.printSelf();
     baseInfo::printSelf();
   }
  PetscFunctionReturn(0);
}
///* ------------------------------------------------------------------- 
//   SetupHp3d - Setup for Hp3d data structures
//
//   Input Parameters:
//   user  - user-defined context 
//   
//*/
//#undef __FUNCT__
//#define __FUNCT__ "SetupHp3d"
//PetscErrorCode SetupHp3d(UniversalCntrl *GenInfo, MPI_Comm ControlComm, 
//                                           vector<PetscInt> *Recvcnts_Elm,
//                                           vector<PetscInt> *Recvcnts_For,
//                                           vector<PetscInt> *Recvcnts_Adj,
//                                           vector<PetscInt> *Displs_Elm,
//                                           vector<PetscInt> *Displs_For,
//                                           vector<PetscInt> *Displs_Adj){
//  vector<MdleInfo> *MdleData = &GenInfo->MdleData; 
//  PetscInt iii,icnt,mdle;   // counters
//  PetscErrorCode info;/* used to check for functions returning nonzeros */
//  PetscFunctionBegin; 
//
//  /* scatter refinement flags 
//     first broadcast refinement level anh then how to refine
//  VecScatterBegin(ElemVecRoot,ElemDataVec,INSERT_VALUES,
//                                          SCATTER_REVERSE,&ScatElem);
//  VecScatterEnd(  ElemVecRoot,ElemDataVec,INSERT_VALUES,
//                                          SCATTER_REVERSE,&ScatElem);
//  info = VecGetArray(ElemDataVec,&elemdata); CHKERRQ(info);
//  FORTRAN_NAME(storerefinementflags)(&elemdata,&NelemLocal,);
//  info = VecRestoreArray(ElemDataVec,&elemdata); CHKERRQ(info); 
//  
//  refine ...
//
//  global element vectors and scattering contexts have changed due to refinement
//  must destroy then recreate
//
//  VecDestroy()
//  info = VecScatterDestroy(GenInfo->ScatElem); CHKERRQ(info);
//
//  */
//
//  PetscLogEventBegin(logevents[26],0,0,0,0); // setup serial 
//  // setup the global dofs orderings
//  // ordering for state solve
//  FORTRAN_NAME(unvisit)(); GenInfo->NDOF_FORWARD[1]=0; mdle=0;
//  for(iii = 0 ; iii < FORTRAN_NAME(get_nreleb)() ; iii++){
//     FORTRAN_NAME(nelconb)(&mdle,&mdle); 
//     FORTRAN_NAME(setup_numbering)(&mdle,&GenInfo->NDOF_FORWARD[1],&idstate,
//                                               &FORTRAN_NAME(put_vert_global),
//                                               &FORTRAN_NAME(put_node_global) );
//  }
//  // ordering for adjoint solve
//  FORTRAN_NAME(unvisit)(); GenInfo->NDOF_ADJOINT[1]=0; mdle=0;
//  for(iii = 0 ; iii < FORTRAN_NAME(get_nreleb)() ; iii++){
//     FORTRAN_NAME(nelconb)(&mdle,&mdle); 
//     FORTRAN_NAME(setup_numbering)(&mdle,&GenInfo->NDOF_ADJOINT[1],&idadjnt,
//                                               &FORTRAN_NAME(put_vert_global),
//                                               &FORTRAN_NAME(put_node_global) );
//  }
//
//  // setup local numbering on  Control Task
//  if(GenInfo->Control_Task){
//     // ordering for state solve
//     FORTRAN_NAME(unvisit)(); GenInfo->NDOF_FORWARD[2]=0; mdle=0;
//     for(iii = 0 ; iii < FORTRAN_NAME(get_nreleb)() ; iii++){
//        FORTRAN_NAME(nelconb)(&mdle,&mdle); 
//        FORTRAN_NAME(setup_numbering)(&mdle,&GenInfo->NDOF_FORWARD[2],&idstate,
//                                               &FORTRAN_NAME(put_vert_local),
//                                               &FORTRAN_NAME(put_node_local) );
//     }
//     // ordering for adjoint solve
//     FORTRAN_NAME(unvisit)(); GenInfo->NDOF_ADJOINT[2]=0; mdle=0;
//     for(iii = 0 ; iii < FORTRAN_NAME(get_nreleb)() ; iii++){
//        FORTRAN_NAME(nelconb)(&mdle,&mdle); 
//        FORTRAN_NAME(setup_numbering)(&mdle,&GenInfo->NDOF_ADJOINT[2],&idadjnt,
//                                               &FORTRAN_NAME(put_vert_local),
//                                               &FORTRAN_NAME(put_node_local) );
//     }
//  } else {
//     /* ----------------   forward problem   ------------------- */
//     // get number of local dof for state problem
//     FORTRAN_NAME(unvisit)(); GenInfo->NDOF_FORWARD[2]=0; mdle=0;
//     for(mdleID = MdleData->begin(); mdleID != MdleData->end() ; mdleID++){
//         FORTRAN_NAME(setup_numbering)(&mdleID->mdle,&GenInfo->NDOF_FORWARD[2],
//                                       &idstate,&FORTRAN_NAME(put_vert_local),
//                                                &FORTRAN_NAME(put_node_local) );
//     }
//     // initialized mapping
//     GenInfo->locmap_for.assign(GenInfo->NDOF_FORWARD[2],0);
//     // create mapping from hp3d data structures
//     FORTRAN_NAME(unvisit)(); GenInfo->NDOF_FORWARD[0] = 0  ;icnt = 0 ;
//     for(mdleID = MdleData->begin(); mdleID != MdleData->end() ; mdleID++){
//          FORTRAN_NAME(createlocalmap)(&mdleID->mdle,&idstate,
//                                       &icnt,&GenInfo->NDOF_FORWARD[0],
//                             &GenInfo->locmap_for[0],&GenInfo->NDOF_FORWARD[2]);
//     }
//     /* ----------------   adjoint problem   ------------------- */
//     // get number of local dof for adjoint problem
//     FORTRAN_NAME(unvisit)(); GenInfo->NDOF_ADJOINT[2]=0; mdle=0;
//     for(mdleID = MdleData->begin(); mdleID != MdleData->end() ; mdleID++){
//         FORTRAN_NAME(setup_numbering)(&mdleID->mdle,&GenInfo->NDOF_ADJOINT[2],
//                                       &idadjnt,&FORTRAN_NAME(put_vert_local),
//                                                &FORTRAN_NAME(put_node_local) );
//     }
//     // initialized mapping
//     GenInfo->locmap_adj.assign(GenInfo->NDOF_ADJOINT[2],0);
//     // create mapping from hp3d data structures
//     FORTRAN_NAME(unvisit)(); GenInfo->NDOF_ADJOINT[0] = 0  ;icnt = 0 ;
//     for(mdleID = MdleData->begin(); mdleID != MdleData->end() ; mdleID++){
//          FORTRAN_NAME(createlocalmap)(&mdleID->mdle,&idadjnt,
//                                       &icnt,&GenInfo->NDOF_ADJOINT[0],
//                             &GenInfo->locmap_adj[0],&GenInfo->NDOF_ADJOINT[2]);
//     }
//     /* alternate ownership of task dofs in the petsc vec 
//     GenInfo->NDOF_FORWARD[0] = GenInfo->NDOF_FORWARD[1]/GenInfo->petscsize;
//     GenInfo->NDOF_ADJOINT[0] = GenInfo->NDOF_ADJOINT[1]/GenInfo->petscsize;
//     if(GenInfo->petscrank == GenInfo->petscsize - 1 ){//using integer division
//        GenInfo->NDOF_FORWARD[0] = GenInfo->NDOF_FORWARD[1] - 
//         (GenInfo->NDOF_FORWARD[1]/GenInfo->petscsize) * (GenInfo->petscsize-1);
//        GenInfo->NDOF_ADJOINT[0] = GenInfo->NDOF_ADJOINT[1] - 
//         (GenInfo->NDOF_ADJOINT[1]/GenInfo->petscsize) * (GenInfo->petscsize-1);
//     }  */
//   
//  } 
//
//  PetscInt cntrlrank , cntrlsize;
//  info = MPI_Comm_rank(ControlComm,&cntrlrank); CHKERRQ(info);
//  info = MPI_Comm_size(ControlComm,&cntrlsize); CHKERRQ(info);
//  // construct global processor element number information
//  vector<PetscInt> cntsbuffelm(cntrlsize+1,0);
//  vector<PetscInt> cntsbufffor(cntrlsize+1,0);
//  vector<PetscInt> cntsbuffadj(cntrlsize+1,0);
//  // array who's ith entry contains the
//  // # of elements on rank i
//  Recvcnts_Elm->resize(cntrlsize,0);
//  Recvcnts_For->resize(cntrlsize,0);
//  Recvcnts_Adj->resize(cntrlsize,0);
//  // array who's ith entry contains the
//  // sum # of elements on rank 0 thru i-1
//  Displs_Elm->resize(cntrlsize+1,0);
//  Displs_For->resize(cntrlsize+1,0);
//  Displs_Adj->resize(cntrlsize+1,0);
//
//  for(iii = cntrlrank+1 ; iii < cntrlsize+1 ; iii++ ){
//      cntsbuffelm[iii] =  GenInfo->MdleData.size();
//      cntsbufffor[iii] =  GenInfo->NDOF_FORWARD[0];
//      cntsbuffadj[iii] =  GenInfo->NDOF_ADJOINT[0];
//  }
//
//  info = MPI_Allreduce(&cntsbuffelm[0],&(*Displs_Elm)[0],cntrlsize+1,
//                          MPI_INTEGER,MPI_SUM,ControlComm);CHKERRQ(info);
//  info = MPI_Allreduce(&cntsbufffor[0],&(*Displs_For)[0],cntrlsize+1,
//                          MPI_INTEGER,MPI_SUM,ControlComm);CHKERRQ(info);
//  info = MPI_Allreduce(&cntsbuffadj[0],&(*Displs_Adj)[0],cntrlsize+1,
//                          MPI_INTEGER,MPI_SUM,ControlComm);CHKERRQ(info);
//
//
//  for(iii = 0 ; iii < cntrlsize ; iii++ ){
//      Recvcnts_Elm->at(iii) = Displs_Elm->at(iii+1) - Displs_Elm->at(iii);
//      Recvcnts_For->at(iii) = Displs_For->at(iii+1) - Displs_For->at(iii);
//      Recvcnts_Adj->at(iii) = Displs_Adj->at(iii+1) - Displs_Adj->at(iii);
//  }
//
//  PetscLogEventEnd(  logevents[26],0,0,0,0); // setup serial 
//  
//  PetscFunctionReturn(0);
//}
///* ------------------------------------------------------------------- 
//   SetupPetscDstruc - Setup for Petsc Solvers
//
//   Input Parameters:
//   user  - user-defined context 
//   
//*/
//#undef __FUNCT__
//#define __FUNCT__ "SetupPetscDstruc"
//PetscErrorCode SetupPetscDstruc(AppSolve *user){
//  UniversalCntrl *GenInfo=user->GenInfo; // global data
//  vector<MdleInfo> *MdleData = &GenInfo->MdleData; 
//  DroneControl *QOIInfo=user->QOIInfo; // global data
//  PetscErrorCode info;
//  IS    locis_for,globis_for,locis_adj,globis_adj;
//  PC pcstate,pcadjoint;
//  PetscScalar rtol; 
//  PetscInt nzerodiag,nzerooffdiag, // # of non-zero entries in parallel matrix
//      mdle,    // mdle node number
//      nrmax,   // max elem dof assuming polynomial order same for mid-edge
//               // mid-face and middle node
//      maxiter; // maximum iterations for solver
//  // used for error handling
//  PetscFunctionBegin; 
//
//  // initialize shape fnc,derivatives,and jacobians at gauss points
//  for(mdleID = MdleData->begin(); mdleID != MdleData->end() ; mdleID++){
//      FORTRAN_NAME(setupgausspoints)(&mdleID->mdle,&mdleID->idelem);
//      FORTRAN_NAME(computeglobalmaps)(&mdleID->mdle,&mdleID->idelem,&idstate);
//      FORTRAN_NAME(computeglobalmaps)(&mdleID->mdle,&mdleID->idelem,&idadjnt);
//  }
//  PetscPrintf(PETSC_COMM_WORLD,
//                    "@@@@@@@@@@@@@@@@@@@@@@@@@gauss points initialized\n" );
//
//  /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//       Create parallel vectors, local vector, and 
//       scattering context to transfer petsc data structures to hp3d
//    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
//  /* ----------------   forward problem   ------------------- */
//  info =  VecCreateMPI(PETSC_COMM_WORLD,GenInfo->NDOF_FORWARD[0],
//                               GenInfo->NDOF_FORWARD[1],&user->x);CHKERRQ(info);
//  info = VecSetFromOptions(user->x);CHKERRQ(info);
//  info = VecDuplicate(user->x,&user->res       );CHKERRQ(info);
//  // create local vector for forward problem
//  info = VecCreateSeq(PETSC_COMM_SELF,GenInfo->NDOF_FORWARD[2],&user->xloc);
//  CHKERRQ(info);
//  info = VecSetFromOptions(user->xloc);CHKERRQ(info);
//  info = VecSet(user->xloc,0.0);CHKERRQ(info);
//  // index set of stride one
//  info=ISCreateStride(PETSC_COMM_SELF,GenInfo->NDOF_FORWARD[2],0,1,&locis_for);
//  CHKERRQ(info);
//  info=ISCreateGeneral(PETSC_COMM_SELF,GenInfo->NDOF_FORWARD[2],
//		            &GenInfo->locmap_for[0],&globis_for); CHKERRQ(info);
//  // create scattering context 
//  info = VecScatterCreate(user->xloc,locis_for,user->x,globis_for,
//		                                &user->GLOCSCAT); CHKERRQ(info);
//  // destroy index sets when done 
//  info=ISDestroy(locis_for ); CHKERRQ(info);
//  info=ISDestroy(globis_for); CHKERRQ(info);
//  /* ----------------   adjoint problem   ------------------- */
//  info =  VecCreateMPI(PETSC_COMM_WORLD,GenInfo->NDOF_ADJOINT[0],
//    	                     GenInfo->NDOF_ADJOINT[1],&user->p);CHKERRQ(info);
//  info = VecSetFromOptions(user->p);CHKERRQ(info);
//  info = VecDuplicate(user->p,&user->p_prev);CHKERRQ(info);
//  info = VecDuplicate(user->p,&user->load);CHKERRQ(info);
//  info = VecDuplicate(user->p,&user->load_prev);CHKERRQ(info);
//  info = VecDuplicate(user->p,&user->load_tmp);CHKERRQ(info);
//  info = VecDuplicate(user->p,&user->load_scratch);CHKERRQ(info);
//  // create local vector for forward problem
//  info = VecCreateSeq(PETSC_COMM_SELF,GenInfo->NDOF_ADJOINT[2],&user->xloc_ADJ);
//  CHKERRQ(info);
//  info = VecSetFromOptions(user->xloc_ADJ);CHKERRQ(info);
//  info = VecSet(user->xloc_ADJ,0.0);CHKERRQ(info);
//  // index set of stride one
//  info=ISCreateStride(PETSC_COMM_SELF,GenInfo->NDOF_ADJOINT[2],0,1,&locis_adj);
//  CHKERRQ(info);
//  // index set for global mapping
//  info=ISCreateGeneral(PETSC_COMM_SELF,GenInfo->NDOF_ADJOINT[2],
//                            &GenInfo->locmap_adj[0],&globis_adj); CHKERRQ(info);
//  // create scattering context 
//  info = VecScatterCreate(user->xloc_ADJ,locis_adj,user->p,globis_adj,
//		                            &user->GLOCSCAT_ADJ); CHKERRQ(info);
//  // destroy index sets when done 
//  info=ISDestroy(locis_adj ); CHKERRQ(info);
//  info=ISDestroy(globis_adj); CHKERRQ(info);
//  /* ----------------   initialize all to zero   ------------------- */
//  info = VecSet(user->x           ,0.0);CHKERRQ(info);
//  info = VecSet(user->res         ,0.0);CHKERRQ(info);
//  info = VecSet(user->p           ,0.0);CHKERRQ(info);
//  info = VecSet(user->p_prev      ,0.0);CHKERRQ(info);
//  info = VecSet(user->load        ,0.0);CHKERRQ(info);
//  info = VecSet(user->load_prev   ,0.0);CHKERRQ(info);
//  info = VecSet(user->load_tmp    ,0.0);CHKERRQ(info);
//  info = VecSet(user->load_scratch,0.0);CHKERRQ(info);
//  /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//    @@@@@@    Create Petsc parallel global matrices   @@@@@@@
//    Create parallel matrix, specifying only its global dimensions.
//    When using MatCreate(), the matrix format can be specified at
//    runtime. Also, the parallel partitioning of the matrix is
//    determined by PETSc at runtime.
//    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
//  // parallel matrix for forward problem
//  /* ----------------   forward problem   ------------------- */
//  FORTRAN_NAME(countmatentries)(&idstate,&nzerodiag,&nzerooffdiag);
//  info = MatCreateMPIAIJ(PETSC_COMM_WORLD,
//         GenInfo->NDOF_FORWARD[0],GenInfo->NDOF_FORWARD[0],
//         GenInfo->NDOF_FORWARD[1],GenInfo->NDOF_FORWARD[1], 
//	 nzerodiag,PETSC_NULL,nzerooffdiag,PETSC_NULL,&user->J); CHKERRQ(info);
//  //fortran store column major storage
//  info = MatSetOption(user->J,MAT_COLUMN_ORIENTED); CHKERRQ(info);
//  //for initialization purposes setup mass matrix for forward problem
//  for(mdleID = GenInfo->MdleData.begin(); mdleID != GenInfo->MdleData.end() ; mdleID++){
//    FORTRAN_NAME(getglobalmap)(&mdleID->mdle,&mdleID->idelem,
//                               &FORTRAN_NAME(get_order_state),&idstate,
//                                                              &nrmax,GLOBMAP);
//    FORTRAN_NAME(elemmass)(&idstate,&mdleID->mdle,&mdleID->idelem,
//                                                  &nrmax,ELEMMATX);
//    info =  MatSetValues(user->J,nrmax,GLOBMAP,nrmax,GLOBMAP,
//                                   ELEMMATX,ADD_VALUES);CHKERRQ(info); 
//  }
//  info = MatAssemblyBegin(user->J,MAT_FINAL_ASSEMBLY); CHKERRQ(info);
//  info = MatAssemblyEnd(user->J,MAT_FINAL_ASSEMBLY); CHKERRQ(info);
//  //let this be the template no new zeros will be added
//  info = MatSetOption(user->J,MAT_NO_NEW_NONZERO_LOCATIONS); CHKERRQ(info);
//  /* ----------------   adjoint problem   ------------------- */
//  FORTRAN_NAME(countmatentries)(&idadjnt,&nzerodiag,&nzerooffdiag);
//  info = MatCreateMPIAIJ(PETSC_COMM_WORLD,
//         GenInfo->NDOF_ADJOINT[0],GenInfo->NDOF_ADJOINT[0],
//         GenInfo->NDOF_ADJOINT[1],GenInfo->NDOF_ADJOINT[1], 
//	 nzerodiag,PETSC_NULL,nzerooffdiag,PETSC_NULL,&user->M); CHKERRQ(info);
//  //fortran store column major storage
//  info = MatSetOption(user->M,MAT_COLUMN_ORIENTED); CHKERRQ(info);
//  //for initialization purposes setup mass matrix for forward problem
//  for(mdleID = GenInfo->MdleData.begin(); mdleID != GenInfo->MdleData.end() ; mdleID++){
//    FORTRAN_NAME(getglobalmap)(&mdleID->mdle,&mdleID->idelem,
//                               &FORTRAN_NAME(get_order_adjnt),&idadjnt,
//                                                              &nrmax,GLOBMAP);
//    FORTRAN_NAME(elemmass)(&idadjnt,&mdleID->mdle,&mdleID->idelem,
//                                                  &nrmax,ELEMMATX);
//    info =  MatSetValues(user->M,nrmax,GLOBMAP,nrmax,GLOBMAP,
//                                   ELEMMATX,ADD_VALUES);CHKERRQ(info); 
//  }
//  info = MatAssemblyBegin(user->M,MAT_FINAL_ASSEMBLY); CHKERRQ(info);
//  info = MatAssemblyEnd(user->M,MAT_FINAL_ASSEMBLY); CHKERRQ(info);
//  //let this be the template no new zeros will be added
//  info = MatSetOption(user->M,MAT_NO_NEW_NONZERO_LOCATIONS); CHKERRQ(info);
//  info = MatDuplicate(user->M,MAT_COPY_VALUES,&user->A);CHKERRQ(info);
//  info = MatDuplicate(user->M,MAT_COPY_VALUES,&user->B);CHKERRQ(info);
//  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//  //   Create PETSC SOLVER contexts for forward and adjoint problem
//  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
//  /* Create nonlinear solver context */
//  info = SNESCreate(PETSC_COMM_WORLD,&user->SNESSTATE); CHKERRQ(info);
//  /* Set function evaluation routine and vector.  Whenever the nonlinear
//     solver needs to compute the nonlinear function, it will call this
//     routine.
//      - Note that the final routine argument is the user-defined
//        context that provides application-specific data for the
//        function evaluation routine.  */
//  info = SNESSetFunction(user->SNESSTATE,user->res,FormFunction,user);
//  CHKERRQ(info);
//
//  // Set Jacobian matrix data structure and Jacobian evaluation routine
//  // Here the matrix that defines the linear system
//  // also serves as the preconditioning matrix.
//
//  info = SNESSetJacobian(user->SNESSTATE,user->J,user->J,FormJacobian,user);
//  CHKERRQ(info);
//
//  // default tolerance for iterative solver
//  rtol = 1.e-9;
//  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  //  Customize nonlinear solver; set runtime options
//  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  //  Set linear solver defaults for this problem. By extracting the
//  //  KSP, KSP, and PC contexts from the SNES context, we can then
//  //  directly call any KSP, KSP, and PC routines to set various options.
//  info = SNESGetKSP(user->SNESSTATE,&user->KSPSTATE); CHKERRQ(info);
//  info = KSPGetPC(user->KSPSTATE,&pcstate); CHKERRQ(info);
//  info = PCSetType(pcstate,PCBJACOBI); CHKERRQ(info);
//  maxiter=GenInfo->NDOF_FORWARD[1];
//  info = KSPSetTolerances(user->KSPSTATE,rtol,PETSC_DEFAULT,PETSC_DEFAULT,
//		  maxiter); CHKERRQ(info);
//  /*  Set SNES/KSP/KSP/PC runtime options, e.g.,
//          -snes_view -snes_monitor -ksp_type <ksp> -pc_type <pc> -ksp_rtol
//          <rtol>
//      These options will override those specified above as long as
//      SNESSetFromOptions() is called _after_ any other customization
//      routines.                                                           */
//  info = SNESSetFromOptions(user->SNESSTATE); CHKERRQ(info);
//
//
//  /* Create linear solver context for dual problem */
//  info = KSPCreate(PETSC_COMM_WORLD,&user->KSPADJOINT); CHKERRQ(info);
//  info = KSPSetOperators(user->KSPADJOINT,user->A,user->A,
//    	      SAME_NONZERO_PATTERN); CHKERRQ(info);
//  // set default solver values
//  info = KSPGetPC(user->KSPADJOINT,&pcadjoint); CHKERRQ(info);
//  info = PCSetType(pcadjoint,PCBJACOBI); CHKERRQ(info);
//  maxiter=GenInfo->NDOF_ADJOINT[1];
//  info = KSPSetTolerances(user->KSPADJOINT,rtol, PETSC_DEFAULT, PETSC_DEFAULT,
//		  maxiter); CHKERRQ(info);
//
//  /*  Set KSP/PC runtime options, e.g.,
//          -ksp_type <ksp> -pc_type <pc> -ksp_rtol <rtol>
//      These options will override those specified above as long as
//      KSPSetFromOptions() is called _after_ any other customization
//      routines.                                                           */
//  info = KSPSetFromOptions(user->KSPADJOINT); CHKERRQ(info);
//  PetscFunctionReturn(0);
//}
///* ------------------------------------------------------------------- 
//   CleanPetscDstruc - Free memory for Petsc Solvers
//
//   Input Parameters:
//   user  - user-defined context 
//   
//*/
//#undef __FUNCT__
//#define __FUNCT__ "CleanPetscDstruc"
//PetscErrorCode CleanPetscDstruc(AppSolve *user){
//  PetscErrorCode info;
//  // used for error handling
//  PetscFunctionBegin; 
//  info = VecDestroy(user->x           ); CHKERRQ(info);
//  info = VecDestroy(user->res         ); CHKERRQ(info);
//  info = VecDestroy(user->p           ); CHKERRQ(info);
//  info = VecDestroy(user->p_prev      ); CHKERRQ(info);
//  info = VecDestroy(user->load        ); CHKERRQ(info);
//  info = VecDestroy(user->load_prev   ); CHKERRQ(info);
//  info = VecDestroy(user->load_tmp    ); CHKERRQ(info);
//  info = VecDestroy(user->load_scratch); CHKERRQ(info);
//  info = VecDestroy(user->xloc        ); CHKERRQ(info);
//  info = VecDestroy(user->xloc_ADJ    ); CHKERRQ(info);
//  info = VecScatterDestroy(user->GLOCSCAT    ); CHKERRQ(info);
//  info = VecScatterDestroy(user->GLOCSCAT_ADJ); CHKERRQ(info);
//  info = MatDestroy(user->J); CHKERRQ(info);
//  info = MatDestroy(user->M); CHKERRQ(info);
//  info = MatDestroy(user->A); CHKERRQ(info);
//  info = MatDestroy(user->B); CHKERRQ(info);
//  info = SNESDestroy(user->SNESSTATE ); CHKERRQ(info);
//  info = KSPDestroy( user->KSPADJOINT); CHKERRQ(info);
//  PetscFunctionReturn(0);
//}
//
///* ------------------------------------------------------------------- 
//   SetupVisualization - Setup for Creating Visualization files
//
//   !!!!!NOTE!!!!!!   all visualization done by interpolating hp mesh
//   !!!!!NOTE!!!!!!   w/ linear hexas
//
//   Input Parameters:
//   *user  - user-defined context 
//   
//*/
//#undef __FUNCT__
//#define __FUNCT__ "SetupVisualization"
//PetscErrorCode SetupVisualization(PetscInt *GroupID,PetscInt* Idopt,
//                                  DroneControl *QOIInfo,PetscInt *Idplot){
//  PetscInt mdle,idelem,iii;
//  PetscErrorCode info;
//  // used for error handling
//  PetscFunctionBegin; 
//
//
//  PetscLogEventBegin(logevents[25],0,0,0,0); // sort and merge
//  // sort and condense the plotted nodes/elements 
//  FORTRAN_NAME(sortandcondense)();
//  PetscLogEventEnd(  logevents[25],0,0,0,0); // sort and merge
//
//
//  // compute mappings to compute the solution
//  mdle = 0 ;
//  for(iii = 0 ; iii < FORTRAN_NAME(get_nreleb)() ; iii++){
//     FORTRAN_NAME(nelconb)(&mdle,&mdle); 
//     idelem = FORTRAN_NAME(get_mdle_elem)(&mdle);
//     FORTRAN_NAME(computeglobalmaps)(&mdle,&idelem,&idstate);
//  }
//  // preproprocess AVSBUFFER
//  FORTRAN_NAME(preprocess_plot)(GroupID);
//  printf("@@@@@@@@@@@@@@@@@plotting data preprocessed\n" );
//
//  // lets see what it looks like
//  // plot bc's and setup of dofs
//  FORTRAN_NAME(output_fem)(GroupID,&QOIInfo->ntask,&idinit,Idopt,Idplot,
//                                                 &nzero, QOIInfo->formfcngauss);
//  // perfusivity thermal conductivity and error fields
//  FORTRAN_NAME(output_fem)(GroupID,&QOIInfo->ntask,&idflds,Idopt,Idplot,
//                                                 &nzero, QOIInfo->formfcngauss);
//
//  PetscFunctionReturn(0);
//}
///* ------------------------------------------------------------------- 
//   Main_Visualization - Final Visualization files and deallocation
//
//   !!!!!NOTE!!!!!!   all visualization done by interpolating hp mesh
//   !!!!!NOTE!!!!!!   w/ linear hexas
//
//   Input Parameters:
//   *user  - user-defined context 
//   
//*/
//#undef __FUNCT__
//#define __FUNCT__ "Main_Visualization"
//PetscErrorCode Main_Visualization(MPI_Comm ControlComm,AppSolve *user,PetscInt *GroupID,PetscInt *Idopt,PetscInt *Idplot,PetscInt *Idstep,UniversalCntrl *GenInfo, DroneControl *QOIInfo,Vec PlotVar_for, Vec PlotVar_adj){
//  PetscScalar *varloc,dum=0.0;
//  PetscErrorCode info;
//  vector<MdleInfo> *MdleData = &GenInfo->MdleData; 
//  vector<PetscInt> recvcnts_for = QOIInfo->recvcnts_for,
//                   recvcnts_adj = QOIInfo->recvcnts_adj,
//                   displs_for   = QOIInfo->displs_for,
//                   displs_adj   = QOIInfo->displs_adj;
//  PetscFunctionBegin; 
//
//  if(GenInfo->Control_Task){
//     PetscInt iii,mdle;
//     VecGetArray(PlotVar_for,&varloc);
//
//     // gather and store solution field
//     info = MPI_Gatherv(&dum,0,MPIU_SCALAR,varloc,&recvcnts_for[0],
//                             &displs_for[0],MPIU_SCALAR,0,ControlComm);
//     mdle=0; // must loop over all elements for control task
//     for(iii = 0 ; iii < FORTRAN_NAME(get_nreleb)() ; iii++){
//        FORTRAN_NAME(nelconb)(&mdle,&mdle); 
//        FORTRAN_NAME(storefield)( &mdle,varloc,
//                               &GenInfo->NDOF_FORWARD[1], &none,&idstate,
//                                                  FORTRAN_NAME(put_vert_temp) , 
//                                                  FORTRAN_NAME(put_node_temp) );
//     }
//
//     // gather and store ideal field
//     info = MPI_Gatherv(&dum,0,MPIU_SCALAR,varloc,&recvcnts_for[0],
//                             &displs_for[0],MPIU_SCALAR,0,ControlComm);
//     mdle=0; // must loop over all elements for control task
//     for(iii = 0 ; iii < FORTRAN_NAME(get_nreleb)() ; iii++){
//        FORTRAN_NAME(nelconb)(&mdle,&mdle); 
//        FORTRAN_NAME(storefield)( &mdle,varloc,
//                               &GenInfo->NDOF_FORWARD[1], &none,&idstate,
//                                                  FORTRAN_NAME(put_vert_ideal), 
//                                                  FORTRAN_NAME(put_node_ideal));
//     }
//
//     // gather and store field representing objective function
//     info = MPI_Gatherv(&dum,0,MPIU_SCALAR,varloc,&recvcnts_for[0],
//                             &displs_for[0],MPIU_SCALAR,0,ControlComm);
//     mdle=0; // must loop over all elements for control task
//     for(iii = 0 ; iii < FORTRAN_NAME(get_nreleb)() ; iii++){
//        FORTRAN_NAME(nelconb)(&mdle,&mdle); 
//        FORTRAN_NAME(storefield)( &mdle,varloc,
//                               &GenInfo->NDOF_FORWARD[1], &none,&idstate,
//                                                  FORTRAN_NAME(put_vert_cost), 
//                                                  FORTRAN_NAME(put_node_cost));
//     }
//     VecRestoreArray(PlotVar_for,&varloc);
//     
//     // gather and store adjoint field
//     VecGetArray(PlotVar_adj,&varloc);
//     info = MPI_Gatherv(&dum,0,MPIU_SCALAR,varloc,&recvcnts_adj[0],
//                             &displs_adj[0],MPIU_SCALAR,0,ControlComm);
//     mdle=0; // must loop over all elements for control task
//     for(iii = 0 ; iii < FORTRAN_NAME(get_nreleb)() ; iii++){
//        FORTRAN_NAME(nelconb)(&mdle,&mdle); 
//        FORTRAN_NAME(storefield)( &mdle,varloc,
//                               &GenInfo->NDOF_ADJOINT[1], &none,&idadjnt,
//                                                  FORTRAN_NAME(put_vert_dual) , 
//                                                  FORTRAN_NAME(put_node_dual) );
//     }
//     VecRestoreArray(PlotVar_adj,&varloc);
//     FORTRAN_NAME(output_fem)(GroupID,&QOIInfo->ntask,&id_slv,Idopt,Idplot,
//                                                 Idstep,QOIInfo->formfcngauss);
//  }else{
//     VecScatter locscat = user->GLOCSCAT;
//
//     //temp field should already be there after the solve
//     // gather the field on the control task
//     info = VecGetArray(user->xloc,&varloc); CHKERRQ(info);
//     for(mdleID=MdleData->begin() ; mdleID!=MdleData->end() ; mdleID++){
//       FORTRAN_NAME(get__field)(&mdleID->mdle,varloc,
//                                &GenInfo->NDOF_FORWARD[2],&user->ISTEP,&idstate,
//                                                  FORTRAN_NAME(get_vert_temp) , 
//                                                  FORTRAN_NAME(get_node_temp) );
//     }
//     info = VecRestoreArray(user->xloc,&varloc); CHKERRQ(info);
//     info =  VecScatterBegin( locscat, user->xloc,PlotVar_for, // local
//                   INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(info);//  to
//     info =  VecScatterEnd( locscat, user->xloc,PlotVar_for, // global
//                   INSERT_VALUES,SCATTER_FORWARD ); CHKERRQ(info);//scatter
//     VecGetArray(PlotVar_for,&varloc);//get pntr to local portion of MPI Vec
//     info = MPI_Gatherv(varloc,GenInfo->NDOF_FORWARD[0],MPIU_SCALAR,
//                &dum,&recvcnts_for[0],&displs_for[0],MPIU_SCALAR,0,ControlComm);
//     VecRestoreArray(PlotVar_for,&varloc);
//
//     //pack and gather the ideal field on the control task
//     info = VecGetArray(user->xloc,&varloc); CHKERRQ(info);
//     for(mdleID=MdleData->begin() ; mdleID!=MdleData->end() ; mdleID++){
//       FORTRAN_NAME(get__field)(&mdleID->mdle,varloc,
//                                &GenInfo->NDOF_FORWARD[2],&user->ISTEP,&idstate,
//                                                  FORTRAN_NAME(get_vert_ideal), 
//                                                  FORTRAN_NAME(get_node_ideal));
//     }
//     info = VecRestoreArray(user->xloc,&varloc); CHKERRQ(info);
//     info =  VecScatterBegin( locscat, user->xloc,PlotVar_for, // local
//                   INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(info);//  to
//     info =  VecScatterEnd( locscat, user->xloc,PlotVar_for, // global
//                  INSERT_VALUES, SCATTER_FORWARD ); CHKERRQ(info);//scatter
//     VecGetArray(PlotVar_for,&varloc);//get pntr to local portion of MPI Vec
//     info = MPI_Gatherv(varloc,GenInfo->NDOF_FORWARD[0],MPIU_SCALAR,
//                &dum,&recvcnts_for[0],&displs_for[0],MPIU_SCALAR,0,ControlComm);
//     VecRestoreArray(PlotVar_for,&varloc);
//
//     //pack and gather the field representing the cost fcn on the control task
//     // compute damage fraction at this time step 
//     FORTRAN_NAME(unvisit)(); 
//     for(mdleID=MdleData->begin() ; mdleID!=MdleData->end() ; mdleID++){
//         FORTRAN_NAME(compdamagefrac)(&mdleID->mdle,&user->ISTEP);
//     } 
//     info = VecGetArray(user->xloc,&varloc); CHKERRQ(info);
//     for(mdleID=MdleData->begin() ; mdleID!=MdleData->end() ; mdleID++){
//       // data for plotting is stored at ..%zdofs_cost(0)
//       FORTRAN_NAME(get__field)(&mdleID->mdle,varloc,
//                                &GenInfo->NDOF_FORWARD[2],&nzero,&idstate,
//                                                  FORTRAN_NAME(get_vert_cost) , 
//                                                  FORTRAN_NAME(get_node_cost) );
//     }
//     info = VecRestoreArray(user->xloc,&varloc); CHKERRQ(info);
//     info =  VecScatterBegin(locscat, user->xloc,PlotVar_for,// local
//                   INSERT_VALUES,  SCATTER_FORWARD ); CHKERRQ(info);//  to
//     info =  VecScatterEnd( locscat, user->xloc,PlotVar_for,// global
//                   INSERT_VALUES,  SCATTER_FORWARD ); CHKERRQ(info);//scatter
//     VecGetArray(PlotVar_for,&varloc);//get pntr to local portion of MPI Vec
//     info = MPI_Gatherv(varloc,GenInfo->NDOF_FORWARD[0],MPIU_SCALAR,
//                &dum,&recvcnts_for[0],&displs_for[0],MPIU_SCALAR,0,ControlComm);
//     VecRestoreArray(PlotVar_for,&varloc);
//
//     //pack and gather the adjoint field on the control task
//     info = VecGetArray(user->xloc,&varloc); CHKERRQ(info);
//     for(mdleID=MdleData->begin() ; mdleID!=MdleData->end() ; mdleID++){
//       FORTRAN_NAME(get__field)(&mdleID->mdle,varloc,
//                                &GenInfo->NDOF_ADJOINT[2],&user->ISTEP,&idadjnt,
//                                                  FORTRAN_NAME(get_vert_dual) , 
//                                                  FORTRAN_NAME(get_node_dual) );
//     }
//     info = VecRestoreArray(user->xloc,&varloc); CHKERRQ(info);
//     info =  VecScatterBegin(locscat, user->xloc,PlotVar_adj,    // local
//                INSERT_VALUES, SCATTER_FORWARD ); CHKERRQ(info);//  to
//     info =  VecScatterEnd( locscat, user->xloc,PlotVar_adj,     // global
//                INSERT_VALUES, SCATTER_FORWARD ); CHKERRQ(info);//scatter
//     VecGetArray(PlotVar_adj,&varloc); //get pntr to local portion of MPI Vec
//     info = MPI_Gatherv(varloc,GenInfo->NDOF_ADJOINT[0],MPIU_SCALAR,
//                &dum,&recvcnts_adj[0],&displs_adj[0],MPIU_SCALAR,0,ControlComm);
//     VecRestoreArray(PlotVar_adj,&varloc);
//  }
//
//  PetscFunctionReturn(0);
//}
