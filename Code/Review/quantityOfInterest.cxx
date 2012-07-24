// C++ include files 
#include <iostream>
#include <fstream>
#include <vector>

// libMesh include files
#include "libmesh.h"
#include "mesh.h"
#include "equation_systems.h"
#include "fe.h"
#include "fe_interface.h"
#include "dof_map.h"
#include "quadrature_gauss.h"
#include "boundary_info.h"
#include "getpot.h"
#include "parallel.h" // mpi utilities
#include "steady_solver.h"
#include "petsc_diff_solver.h"
#include "petsc_vector.h"
#include "petsc_matrix.h"

// tao interface
#include "src/tao_impl.h" 

// dddas include files
#include "applicationContext.h"
#include "pennesInverseModel.h"
#include "pennesInverseSystem.h"
#include "thermal_therapy_system.txx"
#include "quantityOfInterest.h"
#include "dddas.h"
#include "parser_defaults.h"
#include "quantityOfInterest.txx"


///* -------------------------------------------------------------------- */
//#undef __FUNCT__
//#define __FUNCT__ "qoiSpaceTimeDiskRead::initializeSensitivity"
//void qoiSpaceTimeDiskRead::initializeSensitivity( 
//                  TransientFEMSystem &sensitivity_system, 
//                       const PetscInt globalDof )
//{
// PetscErrorCode info; /* used to check for functions returning nonzeros */
// PetscFunctionBegin;
//
// /* open files */
// // set the file output name for current time step
// OStringStream localFileName;
// localFileName << AppSolve::localDisk << "/sensitivity" << globalDof
//           << "rank"  << libMesh::processor_id()  
//           << AppSolve::profileID << ".dat";
//
// // get vec from local disk 
// PetscLogEventBegin(AppSolve::logevents[3],0,0,0,0); // read disk
// info = PetscViewerBinaryOpen(PETSC_COMM_SELF,localFileName.str().c_str(),
//                              FILE_MODE_READ,&localViewHandle);CHKERRV(info);
// info = PetscViewerBinaryRead(localViewHandle,tr,2,PETSC_INT);CHKERRV(info);
// type = tr[0];
// rows = tr[1];
// info = PetscViewerBinaryGetDescriptor(localViewHandle,&localfileHandle);
// CHKERRV(info);
// PetscLogEventEnd(  AppSolve::logevents[3],0,0,0,0); // read disk
//
// PetscFunctionReturnVoid();
//}
///* -------------------------------------------------------------------- */
//#undef __FUNCT__
//#define __FUNCT__ "qoiSpaceTimeDiskRead::finalizeSensitivity"
//void qoiSpaceTimeDiskRead::finalizeSensitivity( )
//{
// PetscErrorCode info; /* used to check for functions returning nonzeros */
// PetscFunctionBegin;
//
// // close file
// PetscLogEventBegin(AppSolve::logevents[15],0,0,0,0); // data transfer
// info = PetscViewerDestroy(localViewHandle);CHKERRV(info);
// PetscLogEventEnd(  AppSolve::logevents[15],0,0,0,0); // data transfer
//
// PetscFunctionReturnVoid();
//}
///* -------------------------------------------------------------------- */
//#undef __FUNCT__
//#define __FUNCT__ "qoiSpaceTimeDiskRead::storeSensitivity"
//void qoiSpaceTimeDiskRead::storeSensitivity( AppSolve *user, 
//                                const int &globalDof,
//                                NumericVector<Number>& current_local_solution )
//{
// PetscErrorCode info; /* used to check for functions returning nonzeros */
// PetscViewer viewHandle; // file handle for petsc viewer
// PetscFunctionBegin;
//
// // FIXME: * is overload in the AutoPtr class to return the raw pointer
// PetscVector<Number>* storeVec = libmesh_cast_ptr<PetscVector<Number>*>( 
//                                                   &current_local_solution );
// //
// // debugging
// //std::cout << "writeSensivity" << iii << std::endl;
// //info = VecView(storeVec->vec(),0);CHKERRV(info);
// //
//
// // set the file output name
// OStringStream file_name;
// file_name << AppSolve::localDisk << "/sensitivity" << globalDof
//           << "rank"  << libMesh::processor_id()  
//           //<< "istep" << AppSolve::ISTEP
//           << AppSolve::profileID << ".dat";
// // store all sensitivities at a given time step in one file on local disk 
// if( AppSolve::ISTEP != user->qoiOptimizer()->Nsteplo()+1 ) 
//  {
//   while( access( file_name.str().c_str() ,W_OK) )
//    {
//      std::cout<< "waiting to write "<<file_name.str()<<std::endl<<std::flush;
//    }
// info = PetscViewerBinaryOpen(PETSC_COMM_SELF,file_name.str().c_str(),
//                             FILE_MODE_APPEND,&viewHandle);CHKERRV(info);
//  }
// else
//  {
// info = PetscViewerBinaryOpen(PETSC_COMM_SELF,file_name.str().c_str(),
//                              FILE_MODE_WRITE,&viewHandle);CHKERRV(info);
//  }
// info = VecView(storeVec->vec(),viewHandle);CHKERRV(info);
// info = PetscViewerDestroy(viewHandle);CHKERRV(info);
//
// PetscFunctionReturnVoid();
//}
///* -------------------------------------------------------------------- */
//#undef __FUNCT__
//#define __FUNCT__ "qoiSpaceTimeDiskRead::getSensitivity"
//void qoiSpaceTimeDiskRead::getSensitivity( AppSolve *user, const int &globalDof,
//           const    std::vector<optimizationParameter*>::iterator idParamIter )
//{
// PetscErrorCode info; /* used to check for functions returning nonzeros */
// EquationSystems &es = *user->_eqnSystem;   // libmesh solver
// PetscFunctionBegin;
//
// // Get a reference to the LinearImplicitSystem for sensitivity problem
// TransientFEMSystem & sensitivity_system =
//      EqnSystem.get_system<TransientFEMSystem>("SensitivitySystem");
//
// // file position
// off_t   currentPos,
//         timeOffset = 
//          ( 2   * PETSC_BINARY_INT_SIZE  + 
//           rows * PETSC_BINARY_SCALAR_SIZE) *
//                              (AppSolve::ISTEP - user->qoiOptimizer()->Nsteplo() - 1)
//          + 2 *  PETSC_BINARY_INT_SIZE  ;
//
// /* get local vector */
// // FIXME: * is overload in the AutoPtr class to return the raw pointer
// PetscVector<Number>* localVec = libmesh_cast_ptr<PetscVector<Number>*>( 
//                                &(*sensitivity_system.current_local_solution) );
// info = PetscBinarySeek(localfileHandle, timeOffset,
//                        PETSC_BINARY_SEEK_SET,&currentPos);CHKERRV(info);
// PetscScalar    *avec;
// info = VecGetArray(localVec->vec(),&avec);CHKERRV(info);
// info = PetscBinaryRead(localfileHandle,avec,rows,PETSC_SCALAR);CHKERRV(info);
// info = VecRestoreArray(localVec->vec(),&avec);CHKERRV(info);
//
// PetscFunctionReturnVoid();
//}
// template instantiations
template class noOptQOI< TransientFEMSystem >;
template class onlyMRTIQOI< TransientFEMSystem >;
template class spaceTimeQOI< PennesInverseSystem<PennesInverseModel> >;
