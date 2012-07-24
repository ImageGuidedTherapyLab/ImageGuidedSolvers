/*
       Python C-style interface routines
 */ 
// system includes
#include <vector>

// C include files 
#include <sys/types.h>
#include <sys/stat.h>

//libmesh
#include "mesh_refinement.h"
#include "equation_systems.h"
#include "fem_system.h"
#include "boundary_info.h"
#include "mesh_generation.h"
#include "fe_base.h"
#include "mesh.h"
#include "fe_interface.h"
#include "quadrature_gauss.h"
#include "petsc_macro.h" // define petsc version PETSC_VERSION_LESS_THAN
#include "dof_map.h"
#include "sparse_matrix.h"
#include "getpot.h"
#include "parallel.h"

//local
#include "optimizationParameter.h"
#include "tttkUtilities.h"
#include "petsc_fem_system.h"

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
     debugfile << "cd ${DDDAS_SRC}; gdb python --pid="<<pid<<std::endl;
     debugfile.close();
     // make the file executable
     chmod(file_name.str().c_str(),S_IREAD|S_IWRITE|S_IEXEC);
     // echo info
     std::cout << getenv("SSH_CONNECTION") << " rank " << rank 
               << " process id is " << pid << std::endl << std::flush ;
     while(!endstall)(sleep(5));
}
 
#undef __FUNCT__
#define __FUNCT__ "PetscPopErrorHandlerAndDebug"
PetscErrorCode PetscPopErrorHandlerAndDebug()
{
  PetscFunctionBegin; 

  /*Check for debugging flags*/
  PetscErrorCode ierr;          PetscTruth  flg;

  // use default petsc error handling
  // TODO merge w/ petsc4py error handling...
  ierr = PetscPopErrorHandler();

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

  PetscFunctionReturn(0);
}

/*@@@@@@@@@@@@@@@call from debugger to view stl vectors@@@@@@@@@@@@@@@*/
PetscErrorCode ViewSTLPetscScalar( PetscScalar *X , PetscInt size ) {
  for(int iii=0;iii < size ;iii++) printf("(%d) %.16e\n",iii,X[iii]);
  fflush(stdout);
  PetscFunctionReturn(0);
}
PetscErrorCode ViewSTLPetscInt( PetscInt *X , PetscInt size  ) {
  for(int iii=0;iii < size ;iii++) printf("(%d) %d\n",iii,X[iii]);
  fflush(stdout);
  PetscFunctionReturn(0);
}
/*@@@@@@@@@@@@@@@call from debugger to view stl vectors@@@@@@@@@@@@@*/
PetscErrorCode ViewSTLPetscScalar( PetscScalar *X , PetscInt istart, PetscInt
istop ) {
  for(int iii=istart;iii < istop ;iii++) printf("(%d) %.16e\n",iii,X[iii]);
  fflush(stdout);
  PetscFunctionReturn(0);
}
PetscErrorCode ViewSTLPetscInt( PetscInt *X , PetscInt istart, PetscInt size ){
  for(int iii=istart;iii < size ;iii++) printf("(%d) %d\n",iii,X[iii]);
  fflush(stdout);
  PetscFunctionReturn(0);
}
// Add an ExplicitSystem to store data
#undef __FUNCT__
#define __FUNCT__ "AddStorageVectors"
PetscErrorCode AddStorageVectors(System* system, char * VectorID, int NumStorageVecs)
{
  PetscFunctionBegin; 
  // add any vectors needed for storage
  for(int iii = 0 ; iii < NumStorageVecs ; iii++)
    {
      // append time step to vector name
      OStringStream vector_name;
      vector_name << VectorID << iii ;
      // see  System::init_data ()
      #ifdef LIBMESH_ENABLE_GHOSTED
      NumericVector<Number> & StorageVec = system->add_vector( vector_name.str(),true, GHOSTED);
      #else
      NumericVector<Number> & StorageVec = system->add_vector( vector_name.str(),true, SERIAL);
      #endif
    }
  PetscFunctionReturn(0);
}
/*@@@@@@@@@@@@@@@@@@print info @@@@@@@@@@@@@@@@@@*/
#undef __FUNCT__
#define __FUNCT__ "MatDataInfo"
PetscErrorCode MatDataInfo(Mat MatOfInterest, const char MatID[])
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  PetscScalar norm_one, norm_frob, norm_infty ;
  ierr = MatNorm(MatOfInterest,   NORM_1     ,&norm_one);CHKERRQ(ierr);
  ierr = MatNorm(MatOfInterest,NORM_FROBENIUS,&norm_frob);CHKERRQ(ierr);
  ierr = MatNorm(MatOfInterest,NORM_INFINITY ,&norm_infty);CHKERRQ(ierr);
  ierr=PetscPrintf(PETSC_COMM_WORLD,
       "  %s |.|_1 =%12.5e |.|_frob =%12.5e |.|_infty =%12.5e\n",
       MatID, norm_one, norm_frob, norm_infty ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/**
 * print matrix norm info to determine the relative effect
 *   of matrix projection into the constrained sparsity pattern
 */
#undef __FUNCT__
#define __FUNCT__ "MatDataSparseLossInfo"
PetscErrorCode MatDataSparseLossInfo(Mat OldMat, 
                                     Mat NewMat,const char MatID[])
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  PetscScalar norm_one, norm_frob, norm_infty ;
  ierr = MatNorm(OldMat,   NORM_1     ,&norm_one);CHKERRQ(ierr);
  ierr = MatNorm(OldMat,NORM_FROBENIUS,&norm_frob);CHKERRQ(ierr);
  ierr = MatNorm(OldMat,NORM_INFINITY ,&norm_infty);CHKERRQ(ierr);
  ierr=PetscPrintf(PETSC_COMM_WORLD,
       "  %s *Orig* |.|_1 =%12.5e |.|_frob =%12.5e |.|_infty =%12.5e\n",
       MatID, norm_one, norm_frob, norm_infty ); CHKERRQ(ierr);
  ierr = MatNorm(NewMat,   NORM_1     ,&norm_one);CHKERRQ(ierr);
  ierr = MatNorm(NewMat,NORM_FROBENIUS,&norm_frob);CHKERRQ(ierr);
  ierr = MatNorm(NewMat,NORM_INFINITY ,&norm_infty);CHKERRQ(ierr);
  ierr=PetscPrintf(PETSC_COMM_WORLD,
       "  %s *New* |.|_1 =%12.5e |.|_frob =%12.5e |.|_infty =%12.5e\n",
       MatID, norm_one, norm_frob, norm_infty ); CHKERRQ(ierr);
  // compute the difference
  //ierr = MatAXPY(OldMat,-1.0,NewMat,DIFFERENT_NONZERO_PATTERN);
  //CHKERRQ(ierr);
  //ierr = MatNorm(OldMat,   NORM_1     ,&norm_one);CHKERRQ(ierr);
  //ierr = MatNorm(OldMat,NORM_FROBENIUS,&norm_frob);CHKERRQ(ierr);
  //ierr = MatNorm(OldMat,NORM_INFINITY ,&norm_infty);CHKERRQ(ierr);
  //ierr=PetscPrintf(PETSC_COMM_WORLD,
  //     "  %s *diff* |.|_1 =%12.5e |.|_frob =%12.5e |.|_infty =%12.5e\n",
  //     MatID, norm_one, norm_frob, norm_infty ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecDataInfo"
PetscErrorCode VecDataInfo(Vec VecOfInterest, const char VecID[])
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  PetscInt    idx;
  PetscReal   min,max,norm;

  ierr = VecMin( VecOfInterest, &idx  , &min ); CHKERRQ(ierr);
  ierr = VecMax( VecOfInterest, &idx  , &max ); CHKERRQ(ierr);
  ierr = VecNorm(VecOfInterest, NORM_2, &norm); CHKERRQ(ierr);

  ierr=PetscPrintf(PETSC_COMM_WORLD,"  %s Min=%12.5e Max=%12.5e Norm=%12.5e\n", 
                                   VecID, min, max, norm); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
/*@@@@@@@@@@ GenerateStructuredGrid w Exodus like info @@@@@@@@@@@@@@@*/
PetscErrorCode GenerateStructuredGrid( Mesh &mesh,
                                    int nx_mesh, int ny_mesh, int nz_mesh,
				    double xmin, double xmax, 
				    double ymin, double ymax, 
				    double zmin, double zmax, 
                                    std::vector<int> &boundaryData) 
{
  PetscFunctionBegin; 

  // Use the MeshTools::Generation mesh generator to create a uniform
  // grid on the square [-1,1]^D.  We instruct the mesh generator
  // to build a mesh of 8x8 \p Quad9 elements in 2D, or \p Hex27
  // elements in 3D.  Building these higher-order elements allows
  // us to use higher-order approximation, as in example 3.
  // 
  // Using ElemType type=INVALID_ELEM:
  //    nz_mesh > 0 ==> 3D  should give same as HEX8
  //    nz_mesh = 0 ==> 2D  should give same as QUAD4
  MeshTools::Generation::build_cube (mesh,
                                     nx_mesh,ny_mesh,nz_mesh,
				     xmin, xmax, 
				     ymin, ymax, 
				     zmin, zmax,
                                     INVALID_ELEM);

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
  libmesh_assert(boundaryData.size() == 6);
  // loops should be same as MeshTools::Generation::build_cube
  for (int k=0; k<nz_mesh ; k++)
   for (int j=0; j<ny_mesh ; j++)
    for (int i=0; i<nx_mesh ; i++)
     {
      unsigned int elemID =  i + j*nx_mesh + k*nx_mesh*ny_mesh; 
      Elem* elem = mesh.elem( elemID );
      if (k == 0)           mesh.boundary_info->add_side(elem, 0, boundaryData[0]);
      if (k == (nz_mesh-1)) mesh.boundary_info->add_side(elem, 5, boundaryData[5]);
      if (j == 0)           mesh.boundary_info->add_side(elem, 1, boundaryData[1]);
      if (j == (ny_mesh-1)) mesh.boundary_info->add_side(elem, 3, boundaryData[3]);
      if (i == 0)           mesh.boundary_info->add_side(elem, 4, boundaryData[4]);
      if (i == (nx_mesh-1)) mesh.boundary_info->add_side(elem, 2, boundaryData[2]);
     }
  PetscFunctionReturn(0);
}
/*@@@@@@@@@@ setup unStructuredGrid @@@@@@@@@@@@@@@*/
#undef __FUNCT__
#define __FUNCT__ "SetupUnStructuredGrid"
PetscErrorCode SetupUnStructuredGrid(Mesh *m_mesh, char *MeshFile,
                                                  int CoarseRefinements,
                                     double A00, double A01, double A02,
                                     double A10, double A11, double A12,
                                     double A20, double A21, double A22,
                                     double  b0, double  b1, double  b2 )
{
  PetscFunctionBegin; 

  // read the mesh.
  std::cout << "reading mesh..." << std::endl;
  m_mesh->read(MeshFile);

  // And an object to refine it
  MeshRefinement mesh_refinement(*m_mesh);
  mesh_refinement.coarsen_by_parents() = true;
  //mesh_refinement.absolute_global_tolerance() = global_tolerance;
  //mesh_refinement.nelem_target() = nelem_target;
  mesh_refinement.refine_fraction() = 0.3;
  mesh_refinement.coarsen_fraction() = 0.3;
  mesh_refinement.coarsen_threshold() = 0.1;

  // possible refinements = 0 ==> no refine
  mesh_refinement.uniformly_refine(CoarseRefinements);

  std::cout << std::endl ; 
  std::cout << A00 << " " << A01 << " " << A02 << " " << b0 <<  std::endl;
  std::cout << A10 << " " << A11 << " " << A12 << " " << b1 <<  std::endl;
  std::cout << A20 << " " << A21 << " " << A22 << " " << b2 <<  std::endl;
  std::cout << std::flush ; 

  // Loop over all the nodes in the mesh and apply a transformation
  libMesh::MeshBase::node_iterator nod     = m_mesh->nodes_begin();
  const libMesh::MeshBase::node_iterator nod_end = m_mesh->nodes_end();

  for ( ; nod != nod_end; ++nod)
    {
      // Store reference. This allows for nicer syntax later.
      Node &NODE = *(*nod);

      Point transformedNode;
      transformedNode(0)= A00*NODE(0)+A01*NODE(1)+A02*NODE(2)+b0;
      transformedNode(1)= A10*NODE(0)+A11*NODE(1)+A12*NODE(2)+b1;
      transformedNode(2)= A20*NODE(0)+A21*NODE(1)+A22*NODE(2)+b2;
      
      NODE = transformedNode;
      
    }
  CHKMEMQ; // check for memory corruption use -malloc_debug to enable

  PetscFunctionReturn(0);
}
/*@@@@@@@@@@@@@@@@@ Set parameter file @@@@@@@@@@@@@@@@@@@@@*/
PetscErrorCode SetParameterIni( Parameters *parameters ,
                                GetPot *controlfile)
{
  PetscFunctionBegin;
  parameters->set<GetPot*>("controlfile") = controlfile;
  PetscFunctionReturn(0);
}
/*@@@@@@@@@@@@@@@@@ flush output @@@@@@@@@@@@@@@@@@@@@*/
PetscErrorCode FlushStdCoutCerr( )
{
  PetscFunctionBegin;
  std::cout << std::flush ; fflush(stdout);
  std::cerr << std::flush ; fflush(stderr);
  PetscFunctionReturn(0);
}
// Store History
#undef __FUNCT__
#define __FUNCT__ "StoreSystemTimeStep"
PetscErrorCode StoreSystemTimeStep(System *system, int CurrentTimeID)
{
  PetscFunctionBegin;
  
  // store the solution history for the adjoint solve
  OStringStream vector_name;
  vector_name << "stored_local_solution" << CurrentTimeID ;
  NumericVector<Number> & StoredVec = system->get_vector( vector_name.str() );
  StoredVec = *system->current_local_solution;

  PetscFunctionReturn(0);
}

// Pass PetscVec to Python
#undef __FUNCT__
#define __FUNCT__ "CopySolutionVector"
PetscErrorCode CopySolutionVector(System & currentSystem, System & copySystem)
{
  PetscFunctionBegin; 

  // get pointer to petsc vector
  PetscVector<Number> &solution =
       *(libmesh_cast_ptr<PetscVector<Number>*>(currentSystem.solution.get()));
  PetscVector<Number> &copySoln =
       *(libmesh_cast_ptr<PetscVector<Number>*>(   copySystem.solution.get()));

  // VecCopy should be overloaded in =
  solution = copySoln ;

  PetscFunctionReturn(0);
}

// Pass PetscVec to Python
#undef __FUNCT__
#define __FUNCT__ "GetSolutionVector"
PetscErrorCode GetSolutionVector(System & currentSystem, Vec *VecData)
{
  PetscFunctionBegin; 

  // get pointer to petsc vector
  PetscVector<Number> &solution =
       *(libmesh_cast_ptr<PetscVector<Number>*>(currentSystem.solution.get()));

  PetscErrorCode ierr;
  VecScatter ctx;
  ierr = VecScatterCreateToAll(solution.vec(),&ctx,VecData);
  // scatter as many times as you need 
  ierr = VecScatterBegin(ctx,solution.vec(),*VecData,INSERT_VALUES,SCATTER_FORWARD);
  ierr = VecScatterEnd(  ctx,solution.vec(),*VecData,INSERT_VALUES,SCATTER_FORWARD);
  // destroy scatter context 
  // NOTE we do not destroy the vector passed back to python
  ierr = VecScatterDestroy(ctx);

  PetscFunctionReturn(ierr);
}

// Pass Python to C++
#undef __FUNCT__
#define __FUNCT__ "SetSolutionVector"
PetscErrorCode SetSolutionVector(System & currentSystem , Vec VecData)
{
  PetscFunctionBegin; 

  // get pointer to petsc vector
  PetscVector<Number> &solution =
       *(libmesh_cast_ptr<PetscVector<Number>*>(currentSystem.solution.get()));

  PetscErrorCode ierr;
  PetscInt pythonSize,libmeshSize;
  ierr = VecGetSize(VecData,        &pythonSize );
  ierr = VecGetSize(solution.vec(),&libmeshSize );
  if(libmeshSize != pythonSize ) // error check
    {
      std::cerr << "libmesh vector length "<< libmeshSize 
                << " not equal python length "<< pythonSize <<  std::endl;
      libmesh_error();
    }

  VecScatter ctx;
  Vec NotUsed;
  ierr = VecScatterCreateToAll(solution.vec(),&ctx,&NotUsed);
  // scatter as many times as you need 
  ierr = VecScatterBegin(ctx,VecData,solution.vec(),INSERT_VALUES,SCATTER_REVERSE);
  ierr = VecScatterEnd(  ctx,VecData,solution.vec(),INSERT_VALUES,SCATTER_REVERSE);
  // destroy scatter context and local vector when no longer needed
  ierr = VecScatterDestroy(ctx);
  ierr = VecDestroy(NotUsed);

  // localize the solution into the proc owned data structures
  currentSystem.solution->localize(*currentSystem.current_local_solution);
  PetscFunctionReturn(ierr);
}
/**
 *  get the petsc subvector for the given variable
 */
#undef __FUNCT__
#define __FUNCT__ "PetscFEMSystemGetSolnSubVector"
PetscErrorCode PetscFEMSystemGetSolnSubVector(System *currentSystem, 
                                       int VarID, Vec *VecData)
{
  PetscFunctionBegin;

  PetscErrorCode ierr;
  // get pointer to system 
  PetscFEMSystem *system = dynamic_cast<PetscFEMSystem *>(currentSystem);

  // get solution subvector
  Vec SubVector;
  system->GetSubVector(SubVector,VarID);

  // create and scatter to all procs
  VecScatter ctx;
  ierr = VecScatterCreateToAll(SubVector,&ctx,VecData);
  // scatter as many times as you need 
  ierr = VecScatterBegin(ctx,SubVector,*VecData,INSERT_VALUES,SCATTER_FORWARD);
  ierr = VecScatterEnd(  ctx,SubVector,*VecData,INSERT_VALUES,SCATTER_FORWARD);
  // destroy scatter context and local vector when no longer needed
  // NOTE we do not destroy the vector passed back to python
  ierr = VecScatterDestroy(ctx);
  ierr = VecDestroy(SubVector);

  PetscFunctionReturn(0);
}

/**
 *  setup intial conditions
 */
#undef __FUNCT__
#define __FUNCT__ "PetscFEMSystemSetupInitialConditions"
PetscErrorCode PetscFEMSystemSetupInitialConditions( System *currentSystem )
{
  PetscFunctionBegin;

  PetscFEMSystem *system = dynamic_cast<PetscFEMSystem *>(currentSystem);
  system->SetupInitialConditions();

  PetscFunctionReturn(0);
}
/**
 *  update the time step 
 */
#undef __FUNCT__
#define __FUNCT__ "PetscFEMSystemUpdateTimeStep"
PetscErrorCode PetscFEMSystemUpdateTimeStep(System *currentSystem ,
                                            int CurrentTimeID)
{
  PetscFunctionBegin;

  
  PetscFEMSystem *system = dynamic_cast<PetscFEMSystem *>(currentSystem);
  system->SetPowerID(CurrentTimeID);

  // store the global data for Thermal dose update
  // local to global map not available
  NumericVector<Number> &old_global_solution=system->get_vector("old_global_solution");
  system->solution->close();
  old_global_solution= *system->solution;

  // store the local data for assembly
  NumericVector<Number> &old_local_solution =system->get_vector("old_local_solution");
  system->solution->localize(*system->current_local_solution);
  old_local_solution     =   *system->current_local_solution;

  PetscFunctionReturn(0);
}

/* setup measurement matrix to transform the state to the roi */
#undef __FUNCT__
#define __FUNCT__ "PetscFEMSystemCreateNodeSetFromMask"
PetscInt PetscFEMSystemCreateNodeSetFromMask( System *currentSystem ,
                                               PetscScalar labelValue,
                                               PetscInt nodeSetID)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  // the # of measurements can is limited by the ROI size or the number of nodes
  libMesh::MeshBase &mesh = currentSystem->get_mesh();

  // get dof map
  PetscFEMSystem *state_system = 
                     dynamic_cast<PetscFEMSystem *>(currentSystem);
  const DofMap & dof_map = state_system->get_dof_map();

  //  hold variable dofs
  std::vector<unsigned int> dof_indices_var;

  for( unsigned int i_var = 0 ; i_var < state_system->n_vars() ; i_var++)
    {
     // image map reserve room for entire state
     std::vector<PetscInt> &dirichletNodeSet =  state_system->m_NodeSets[i_var][nodeSetID];
     if(dirichletNodeSet.size())
       {
        std::cout << nodeSetID << "nodeset already in use" << std::endl;
        std::cout <<  std::flush; libmesh_error();
       }
     //  loop over elements and create map
     libMesh::MeshBase::const_element_iterator       el    =mesh.active_local_elements_begin();
     const libMesh::MeshBase::const_element_iterator end_el_full=mesh.active_local_elements_end();
     for ( ; el != end_el_full; ++el)
      {
       // Store element subdomain id to allow different sets of equations to be
       // solved on different parts of the mesh
       // ASSUMES SUBDOMAIN ORDERING STARTS FROM ONE but we need a 0 based
       // numbering scheme for std::vector
       Elem* elem = *el;
       const unsigned int subdomain_id = elem->subdomain_id() - 1;

       // only add healthy tissue domain
       if( subdomain_id == 0 )
          {
           // degree of freedom indices for individual components
           dof_map.dof_indices (elem,   dof_indices_var ,i_var);
           // nodal based setup of dirichlet data
           // indices should be global
           for(unsigned int Ii = 0 ; Ii < dof_indices_var.size() ; Ii++)
            if( std::abs(state_system->current_solution( dof_indices_var[Ii] )-labelValue) < 1.e-6 ) 
                dirichletNodeSet.push_back( dof_indices_var[Ii] );  
          }
      } // end element loop
 
    // broadcast to all procs
    Parallel::allgather (dirichletNodeSet);
    // sort then erase duplicates
    std::sort( dirichletNodeSet.begin(), dirichletNodeSet.end() );
    std::vector<PetscInt>::iterator pos;

    pos = std::unique(dirichletNodeSet.begin(),dirichletNodeSet.end());
    dirichletNodeSet.erase( pos,dirichletNodeSet.end() ); 
  }

  PetscFunctionReturn(state_system->m_NodeSets[0][nodeSetID].size());
}
// Add an ExplicitSystem to store data
#undef __FUNCT__
#define __FUNCT__ "AddExplicitSystem"
System * AddExplicitSystem (char *SystemName, 
              EquationSystems *m_eqn_systems)
{
  PetscFunctionBegin; 

  ExplicitSystem & system = 
    m_eqn_systems->add_system<ExplicitSystem> (SystemName);

  PetscFunctionReturn( &system  );
  
}

// Get a pointer to SystemName
#undef __FUNCT__
#define __FUNCT__ "GetSystem"
System * GetSystem (char *SystemName, 
              EquationSystems *m_eqn_systems)
{
  PetscFunctionBegin; 

  System &system = m_eqn_systems->get_system (SystemName);

  PetscFunctionReturn( &system  );
  
}

// subroutine to apply dirichlet data to system
#undef __FUNCT__
#define __FUNCT__ "ApplyDirichletData"
PetscErrorCode 
ApplyDirichletData(EquationSystems *m_eqnSystems,
                   System *currentSystem, 
                   PetscInt varID, PetscInt nodeSetID,
                   PetscScalar DirichletValue )
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  //  get reference to dirichlet nodes
  PetscFEMSystem &state_system = m_eqnSystems->get_system<PetscFEMSystem>("StateSystem");
  std::vector<PetscInt> &dirichletNodeSet= state_system.m_NodeSets[varID][nodeSetID];

  // apply dirichlet data if any
  if( dirichletNodeSet.size() )
    {
      // assume measurement is the same as dirichlet data
      currentSystem->solution->close();
      for( unsigned int Ii = 0; Ii<dirichletNodeSet.size();Ii++) 
        currentSystem->solution->set(dirichletNodeSet[Ii],DirichletValue);
      currentSystem->solution->close();
    }
  PetscFunctionReturn(0);
}

// Add an ExplicitSystem to store data
#undef __FUNCT__
#define __FUNCT__ "AddConstantMonomialVariable"
unsigned int AddConstantMonomialVariable(System *currentSystem, 
                                             char *variable_name)
{
  PetscFunctionBegin; 
  PetscFunctionReturn(currentSystem->add_variable (variable_name,CONSTANT,MONOMIAL ));
}

// Add an ExplicitSystem to store data
#undef __FUNCT__
#define __FUNCT__ "AddFirstLagrangeVariable"
unsigned int AddFirstLagrangeVariable(System *currentSystem, 
                                          char *variable_name)
{
  PetscFunctionBegin; 
  PetscFunctionReturn(currentSystem->add_variable (variable_name, FIRST,   LAGRANGE ));
}

/* -------------------------------------------------------------------- */
//  compute the weighted ||.||_2 norm, (Computed-Measurement)/Normalization
#undef __FUNCT__
#define __FUNCT__ "WEIGHTEDL2Norm"
PetscScalar WEIGHTEDL2Norm(System * state_system, // the computed solution
                           char    *StateVariable , // name of state variable
                           System * ideal_system  , // measured data
                           char    *IdealVariable , // name of ideal variable 
                           System * normal_system ,// normalization
                           char    *NormalVariable ) // name of normalization variable 
{
  // used for error handling
  PetscFunctionBegin; 

  //PetscLogEventBegin(AppSolve::logevents[21],0,0,0,0); // Setup libMesh 

  // Get a constant reference to the mesh object.
  // should be the same mesh
  libmesh_assert( state_system->get_mesh() == ideal_system->get_mesh() );
  const libMesh::MeshBase& mesh = state_system->get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the dofmap and mesh for that system
  const DofMap& dof_map_state =  state_system->get_dof_map();
  const DofMap& dof_map_ideal =  ideal_system->get_dof_map();
  const DofMap& dof_map_normal= normal_system->get_dof_map();

  // Construct Quadrature rule based on default quadrature order
  const unsigned int  state_var=  state_system->variable_number( StateVariable);
  const unsigned int  ideal_var=  ideal_system->variable_number( IdealVariable);
  const unsigned int normal_var= normal_system->variable_number(NormalVariable);
  const FEType& fe_type  = dof_map_state.variable_type(state_var);

  // Build the Gauss quadrature rules for numerical integration. 
  QGauss qrule(dim  , fe_type.default_quadrature_order() );

  // Construct finite element object
  AutoPtr<FEBase> fe(FEBase::build(mesh.mesh_dimension(), fe_type));

  // Attach quadrature rule to FE object
  fe->attach_quadrature_rule (&qrule);
  
  // The Jacobian*weight at the quadrature points.
  const std::vector<Real>& JxW                               = fe->get_JxW();
  
  // The value of the shape functions at the quadrature points
  // i.e. phi(i) = phi_values[i][qp] 
  const std::vector<std::vector<Real> >&  phi_values         = fe->get_phi();
  
  // The global degree of freedom indices associated
  // with the local degrees of freedom.
  std::vector<unsigned int>  state_dof_indices;
  std::vector<unsigned int>  ideal_dof_indices;
  std::vector<unsigned int> normal_dof_indices;

  //PetscLogEventEnd(  AppSolve::logevents[21],0,0,0,0); // Setup libMesh 

  //PetscLogEventBegin(AppSolve::logevents[23],0,0,0,0); // qoi eval

  // This function must be run on all processors at once
  parallel_only();

  // intialize error buffer
  Number error_val = 0.; 

  // Begin the loop over the elements
  //
  libMesh::MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 
  for ( ; el != end_el; ++el)
    {
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el;

      // reinitialize the element-specific data
      // for the current element
      fe->reinit (elem);

      // Get the local to global degree of freedom maps
       dof_map_state.dof_indices (elem, state_dof_indices, state_var);
       dof_map_ideal.dof_indices (elem, ideal_dof_indices, ideal_var);
      dof_map_normal.dof_indices (elem,normal_dof_indices,normal_var);
      
      //
      // Begin the loop over the Quadrature points.
      //
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          Number u_h = 0., u_ideal = 0., u_uncert = 0.;
          // Compute solution values at the current
          // quadrature point.  This reqiures a sum
          // over all the shape functions evaluated
          // at the quadrature point.
          for (unsigned int i=0; i<state_dof_indices.size(); i++)
            {
              // Values from current solution.
              u_h      += phi_values[i][qp]* state_system->current_solution(
                                    state_dof_indices[i]);
              u_ideal  += phi_values[i][qp]* ideal_system->current_solution(
                                    ideal_dof_indices[i]);
              u_uncert += phi_values[i][qp]*normal_system->current_solution(
                                   normal_dof_indices[i]);
            }
          // compute normalized error 
          const Number val_error = (u_h - u_ideal)/u_uncert;

          // Add the squares of the error to each contribution
          error_val += JxW[qp]*val_error*val_error;
          
        } // end qp loop
    } // end element loop

  //PetscLogEventEnd(  AppSolve::logevents[23],0,0,0,0); // qoi eval

  // Add up the error values on all processors
  PetscScalar qoi_loc = error_val; 
  Parallel::sum(qoi_loc);

  //if( !libMesh::processor_id() ) 
  //  {
  //   for(AppSolve::ISTEP =user->qoiOptimizer()->Nsteplo() ; 
  //       AppSolve::ISTEP<=user->qoiOptimizer()->Nstephi() ; AppSolve::ISTEP++)
  //        std::cout << "qoi"<<normal_system.name()<<"["<<AppSolve::ISTEP
  //                  << "] = " << std::setprecision(10) 
  //                  << qoi_loc[AppSolve::ISTEP] << std::endl;
  //  }
  PetscFunctionReturn( qoi_loc ); 
}
/* -------------------------------------------------------------------- */
// factor matrix w/ arbitrary solver
#undef __FUNCT__
#define __FUNCT__ "GETFactorMat"
Mat GETFactorMat(Mat  OrigMatrix, // the unfactored matrix
                 char *Solvertype , // solver
                 char *Ordertype  ) //  order
{
  // used for error handling
  PetscFunctionBegin; 
  PetscErrorCode ierr;

  Mat FactorMatrix;
  ierr = MatGetFactor(OrigMatrix,Solvertype, 
                      MAT_FACTOR_LU,&FactorMatrix); 
  IS isrow,iscol;
  ierr = MatGetOrdering(OrigMatrix,Ordertype,&isrow,&iscol);

  MatFactorInfo factorInfo;
  ierr = MatFactorInfoInitialize(&factorInfo); 

  // factor
  MatLUFactorSymbolic(FactorMatrix,OrigMatrix,isrow,iscol,&factorInfo);
  ierr = MatLUFactorNumeric(FactorMatrix,OrigMatrix,&factorInfo);
 	 
  PetscFunctionReturn( FactorMatrix ); 
}
