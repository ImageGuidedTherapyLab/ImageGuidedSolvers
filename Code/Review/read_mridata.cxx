// C++ include files 
#include <iostream>
#include <fstream>
#include <vector>

// C include files 
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>

// ITK 
#include "itkVector.h"
#include "itkImage.h"
#include "itkImageRegionIterator.h" 
#include "itkImageSliceConstIteratorWithIndex.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"

// tao
#include "tao.h"

// libMesh include files
#include "libmesh.h"
#include "mesh.h"
#include "equation_systems.h"
#include "exact_solution.h"
#include "getpot.h" // file parser
#include "nonlinear_solver.h"
#include "nonlinear_implicit_system.h"
#include "linear_implicit_system.h"
#include "dense_submatrix.h"

// dddas include files
#include "applicationContext.h" 
#include "baseInfo.h" 
#include "pennesInverseModel.h" 
#include "quantityOfInterest.h"
#include "transient_fem_system.h"
#include "parser_defaults.h" // ensure same default values between C and Fortran
#include "byte_swapping.h"
#include "fortrandftranslate.h" // header to call fortran routines
#include "pennesVerification.h" 

// CImg
#include "CImg.h"
using namespace cimg_library;

#include "read_mridata.h" // mridata interface 
static CImgList<AVSDATA_TYPE> ImagingData;

namespace dddas
{
/* ------------------------------------------------------------------- 
   default constructor to for Image class
   setup params that do not depend on dicom header data */
Image::Image(GetPot &controlfile) 
{
  // create space and initialize
  MRTI_MEM.reserve(baseInfo::MRTI_ntime+1); 
  MRTI_MEM.resize( baseInfo::MRTI_ntime+1,0); 

  // default is not to write...
  WRITEAVS   = PETSC_FALSE; 
  WRITEMHA   = PETSC_FALSE; 
  WRITERAWIV = PETSC_FALSE; 

  // filename parameters
  std::string nofile("/no/file%d_im%d.bin");
  FILEIN = controlfile("mrti/filein",nofile);
  UNCERTAINTYIN = controlfile("mrti/uncertaintyin",nofile);

  // write an avs file
  if(controlfile("output/mrtiavs",false))     WRITEAVS = PETSC_TRUE ; 
  if(controlfile("mrti/transfer",false) )
    { 
     std::cout<<"Image overwriting writeavs file w/ TRUE!"<<std::endl;
     WRITEAVS = PETSC_TRUE ; 
    }
  // write a mha file
  if(controlfile("output/mrtimha",false))     WRITEMHA = PETSC_TRUE ; 
  // write a rawiv file
  if(controlfile("output/mrtirawiv",false)) WRITERAWIV = PETSC_TRUE ; 

  // echo input
  std::cout <<"\nImage WRITEAVS="  << WRITEAVS <<" WRITERAWIV="<<WRITERAWIV
            <<"\nImage WRITEMHA="  << WRITEMHA 
            <<" FILEIN = "         << FILEIN 
            <<" UNCERTAINTYIN = "  << UNCERTAINTYIN << std::endl;

  // initialize spacing in was to use w/o imaging data
  spacing[0]=(PetscScalar)controlfile("mrti/dx",1.0);
  spacing[1]=(PetscScalar)controlfile("mrti/dy",1.0);
  spacing[2]=(PetscScalar)controlfile("mrti/dz",1.0);

  // initialize conversion factors
  sigma = controlfile("mrti/sigma", 1.e9); // default uncertainty
  CONVERSIONFACT = controlfile("mrti/conversionfact", 0.0);
  AVSDATA_TYPE u_init = controlfile("initial_condition/u_init",37.0);//celciusk
  AVSDATA_TYPE probe_init = controlfile("initial_condition/probe_init",21.0);
  APPLICATORDIFF = controlfile("mrti/applicatordiff", u_init - probe_init);
} 
/* ------------------------------------------------------------------- 
   Turn off image write if image data structures not setup*/
void Image::NoWrite() 
{
  std::cout<< "turning off image file write "<< std::endl;
  WRITEAVS   = PETSC_FALSE; 
  WRITEMHA   = PETSC_FALSE; 
  WRITERAWIV = PETSC_FALSE; 
}
// wait for dicom header info to be transferred before proceeding..
void Image::ImageSetup(AppSolve &user)
{
  PetscInt izero=0,ithree=3,rank;

  std::string MRTIBase = FILEIN.substr( 0 , FILEIN.rfind('/') ); 

  std::ostringstream MRTIFile ; // filename
  MRTIFile << MRTIBase<< "/mrti.ini";

  // wait for dicom header info to be transferred before proceeding..
  while( access( MRTIFile.str().c_str() ,R_OK) )
   {
     std::cout<< "waiting for dicom header info... "<< std::endl;
     std::cout<< MRTIFile.str() <<   " not found..."<< std::endl;
     GetPot controlfile(dfltINIFile);
     std::string dfltFile= FILEIN;
     FILEIN= controlfile("mrti/filein",dfltFile);
     MRTIBase = FILEIN.substr( 0 , FILEIN.rfind('/') ); 
     MRTIFile.str(""); // reset before reuse
     MRTIFile << MRTIBase<< "/mrti.ini";
     sleep(5);
   }

  // instantiate ini file
  GetPot controlfile(MRTIFile.str());

  // physical space dimensions
  origin[0] =(PetscScalar)controlfile("mrti/x0",0.0);
  origin[1] =(PetscScalar)controlfile("mrti/y0",0.0);
  origin[2] =(PetscScalar)controlfile("mrti/z0",0.0);
  spacing[0]=(PetscScalar)controlfile("mrti/dx",0.0);
  spacing[1]=(PetscScalar)controlfile("mrti/dy",0.0);
  spacing[2]=(PetscScalar)controlfile("mrti/dz",0.0);

  // dimensions of MRTI data
  int xpixel=controlfile("mrti/xpixel"  ,256 );
  int ypixel=controlfile("mrti/ypixel"  ,256 );
  int nslice=controlfile("mrti/nslice"  , 1  );

  // error checking after input data 
  if(!nslice || !xpixel || !ypixel || !spacing[0] || !spacing[1] || !spacing[2] ){
     std::cout << "\nmrtiImage X0 = "<< origin[0] << " DX = "<< spacing[0] << " xpixel="<< xpixel
          << "\nmrtiImage Y0 = "<< origin[1] << " DY = "<< spacing[1] << " ypixel="<< ypixel
          << "\nmrtiImage Z0 = "<< origin[2] << " DZ = "<< spacing[2] << " nslice="<< nslice
          << "\n\nmrtiImage an image dimension = 0!!!\n\n"; abort();
  }

  // echo input
  std::cout <<"\nImage X0="<<origin[0]<<" DX="<<spacing[0]<<" xpixel="<< xpixel
            <<"\nImage Y0="<<origin[1]<<" DY="<<spacing[1]<<" ypixel="<< ypixel
            <<"\nImage Z0="<<origin[2]<<" DZ="<<spacing[2]<<" nslice="<< nslice 
            << std::endl;
  /*----------setup ITK data structure -------------------*/

  //  A filter object created using the \code{New()} method is then
  //  assigned to a \code{SmartPointer}.
  //  
  // \index{itk::ImportImageFilter!Pointer}
  // \index{itk::ImportImageFilter!New()}
  importFilter = ImportFilterType::New();      

  size[0]  = xpixel;  // size along X
  size[1]  = ypixel;  // size along Y
  size[2]  = nslice;  // size along Z

  ImportFilterType::IndexType start;
  start.Fill( 0 );

  ImportFilterType::RegionType region;
  region.SetIndex( start );
  region.SetSize(  size  );

  importFilter->SetRegion( region );

  // set physical spacing
  importFilter->SetOrigin( origin );
  importFilter->SetSpacing( spacing );

  // allocate base temperature image 
  ImageType::RegionType    fullRegion;
  fullRegion.SetSize(  size  );

  baseImage = ImageType::New();
  baseImage->SetRegions( fullRegion );
  baseImage->Allocate();
  baseImage->SetSpacing( spacing );
  baseImage->SetOrigin(  origin  );

  // dummy image to hold all ones
  sigmaImage = ImageType::New();
  sigmaImage->SetRegions( fullRegion );
  sigmaImage->Allocate();
  sigmaImage->SetSpacing( spacing );
  sigmaImage->SetOrigin(  origin  );
  sigmaImage->FillBuffer( sigma );

  // open using ITK library of file types
  typedef itk::ImageSeriesReader< ImageType >   LocalReaderType;
  LocalReaderType::Pointer  reader = LocalReaderType::New();

  reader->SetFileName( controlfile("mrti/baseimage", "baseimage.mha" ));

  try
    {
      reader->Update();
      // Define the iterators
      itk::ImageRegionIterator<ImageType> inputIt(reader->GetOutput() , 
                                     reader->GetOutput()->GetRequestedRegion()),
                                     outputIt(baseImage,
                                              baseImage->GetRequestedRegion()); 

      inputIt.GoToBegin();
      outputIt.GoToBegin();

      while( !inputIt.IsAtEnd() ) 
        {
        outputIt.Set( CONVERSIONFACT - APPLICATORDIFF * inputIt.Get() );
        ++inputIt;
        ++outputIt;
        }
      CHKMEMA; // check for memory corruption use -malloc_debug to enable
    }
  catch (itk::ExceptionObject &excp)
    {
    Point p; 
    PetscScalar u_fill =  CONVERSIONFACT +
       user.pdeSolver->InitValues(0,0,p, user._eqnSystem->parameters);
    std::cout << "Exception thrown reading baseline image" << excp << std::endl;
    std::cout << "fill baseline image w/  " << u_fill << std::endl;
    baseImage->FillBuffer( u_fill );
    }

  return;
}
/* --------- Set Image Pointers in libMesh data structures --------------------- */
void Image::SetImagePointers(AppSolve &user)
{
  /*  create a list of images for thermal data each 
      image initialized to  minimum temperature value 
      ASSUME THAT THE INITIAL IMAGE is labeled ZERO!!!   */
  Point dum;
  AVSDATA_TYPE inittemp = 
       user.pdeSolver->InitValues(0,0,dum, user._eqnSystem->parameters);

  // allocate CImgList Data structures for server
  ImagingData.assign(1,size[0],size[1],size[2],1,inittemp);

  //  The buffer is passed to the ImportImageFilter with the
  //  \code{SetImportPointer()}. Note that the last argument of this method
  //  specifies who will be responsible for deleting the memory block once it
  //  is no longer in use. A \code{false} value indicates that the
  //  ImportImageFilter will not try to delete the buffer when its
  //  destructor is called. A \code{true} value, on the other hand, will allow
  //  the
  //  filter to delete the memory block upon destruction of the import filter.
  //
  //  For the ImportImageFilter to appropriately delete the
  //  memory block, the memory must be allocated with the C++
  //  \code{new()} operator.  Memory allocated with other memory
  //  allocation mechanisms, such as C \code{malloc} or \code{calloc}, will not
  //  be deleted properly by the ImportImageFilter.  In
  //  other words, it is the application programmer's responsibility
  //  to ensure that ImportImageFilter is only given
  //  permission to delete the C++ \code{new} operator-allocated memory.
  bool importImageFilterWillOwnTheBuffer = false;
  importFilter->SetImportPointer( ImagingData[0].ptr(), ImagingData[0].size(),
                                  importImageFilterWillOwnTheBuffer );
  // must Update() before can user GetOutput()
  importFilter->Update();
  interpolator = InterpolatorType::New();
  interpolator->SetInputImage( importFilter->GetOutput() );
  return;
}
/* ------------------------------------------------------------------- 
  ImageServer Constructor - Initialize data Structures for C/C++ code
   setup default data to not run the image server... 
*/
ImageServer::ImageServer(GetPot &controlfile,
             std::vector<qoiBaseClass> &qoiOptimizer):Image(controlfile)  // base class constructor
{
  // Create array that contains precomputed groups that need data
  // initialize all to false
  // Need_to_send.resize(baseInfo::MRTI_ntime+1,
  //                               NeedData(qoiOptimizer.size(),PETSC_FALSE)); 
}
/* ------------------------------------------------------------------- 
  ServerSetup - Initialize data Structures for Image Server

*/
void ImageServer::ServerSetup(AppSolve &user,std::vector<qoiBaseClass>
&qoiOptimizer)
{

  PetscFunctionBegin; 

  this->ImageSetup(user);

  // allocate CImgList Data structures for server
  Point dum;
  AVSDATA_TYPE inittemp = 
       user.pdeSolver->InitValues(0,0,dum, user._eqnSystem->parameters);
  ImagingData.assign(baseInfo::MRTI_ntime+1,size[0],size[1],size[2],1,inittemp);

  // instantiate ini file
  GetPot controlfile(dfltINIFile);

  // strcasecmp != 0   ==> strings don't match (case insensitive)
  if(controlfile("output/byteswap",false)){ BYTESWAP = PETSC_TRUE ; }
          else                            { BYTESWAP = PETSC_FALSE; }
  if(BYTESWAP) printf("Byte SWAPPING!!!!!!!!!!!\n");
  std::string byteid;
  if(FORTRAN_NAME(islittleendian)()){
    PetscPrintf(PETSC_COMM_SELF,"native data type is LITTLE ENDIAN\n");
    byteid="litend"; if(BYTESWAP) byteid="bigend";
  }else{
    PetscPrintf(PETSC_COMM_SELF,"native data type is BIG    ENDIAN\n");
    byteid="bigend"; if(BYTESWAP) byteid="litend";
  }

  // cooridinate alignment for rawiv data 
  RAWORIGIN[0]= origin[0] ;
  RAWORIGIN[1]= origin[1] ;
  RAWORIGIN[2]= origin[2] ;
  
  // get the # of threads to use to read in mrti data slices
  NDATATHRD = controlfile("compexec/ndatathrd",1);

  std::string filebase     = controlfile("mrti/filebase","slice");
  std::string profileid    = controlfile("compexec/profileid","");
  std::string compwritemri = controlfile("output/compwritemri",dfltmriloc);
  std::string vis_sizefile = controlfile("avs/mrtivis_sizefile","mrtivis.size");
  FILEOUT      = compwritemri + "/" + profileid + filebase + byteid + "%04d.%s";
  VIS_SIZEFILE = compwritemri + "/" + vis_sizefile;
  if(baseInfo::Control_Task){
    std::cout  << "\nImageServer BYTESWAP= " << BYTESWAP
          << "\nImageServer FILEOUT = " << FILEOUT
          << "\nImageServer NDATATHRD=" << NDATATHRD
          << "\nImageServer VIS_SIZEFILE= "  <<VIS_SIZEFILE
          << "\nImageServer RAWORIGIN=" <<RAWORIGIN[0]<<RAWORIGIN[1]<<RAWORIGIN[2] << std::endl ;
  }

  // declate iterator
  //std::vector<qoiBaseClass>::iterator QoiIter;

  //for(QoiIter = qoiOptimizer.begin() ; QoiIter != qoiOptimizer.end() ; QoiIter++) {
  //   PetscInt iii = distance(qoiOptimizer.begin() , QoiIter);
  //   for(PetscInt IDopt = 0 ; IDopt < QoiIter->NoptSteps ; IDopt++) {
  //      //initial condition
  //      if(baseInfo::IC_MRTI) {
  //         // repeating logic from compute drone only send IC after first
  //         //   optimization step
  //         if(IDopt) Need_to_send.at(QoiIter->IDEAL_NZERO 
  //                                    + IDopt * QoiIter->NOFFSET 
  //                                    - baseInfo::TIMELAG_IC).at(iii)=PETSC_TRUE;
  //      } // end initial condition
  //      //optimizations need to send data for calibrations only
  //      //if (QoiIter->compobj.find("calibration")!=std::string::npos)
  //      //{
  //      //  // check if computing entire norm in time
  //      //  if(baseInfo::COMPUTEFINALOBJECTIVE){
  //      //     for(PetscInt jjj  = baseInfo::MRTI_nzero ;
  //      //                  jjj <= baseInfo::MRTI_ntime ; jjj++) 
  //      //                      Need_to_send.at(jjj     ).at(iii) = PETSC_TRUE;
  //      //  } else {
  //      //     for(PetscInt expectID  = QoiIter->IDEAL_NZERO 
  //      //                                + IDopt * QoiIter->NOFFSET ; 
  //      //                  expectID <= QoiIter->IDEAL_NTIME 
  //      //                                + IDopt * QoiIter->NUMIDEAL; 
  //      //                  expectID++)
  //      //                      Need_to_send.at(expectID).at(iii) = PETSC_TRUE;
  //      //  }
  //      //} // end check for calibration problem
  //   } // end loop over optimization steps
  //}

  return;
}
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
PetscTruth RealTimeImage::CheckImageData(PetscInt IDfile)
{
  char file_i[MAXLEN];
  const int isleep = 1;
  PetscTruth GoToNextOptStep=PETSC_TRUE;
  std::ifstream file;
  std::ifstream::pos_type size; // used to determine file size in bytes


  // if not writing any vis files no need to worry about reading it in
  // and return to the Data_Server that the file is found
  // **NOTE** for more than ONE optimization step this will kill the
  // **NOTE**  OPTIMIZATION STEP AS WELL
  if(!WRITEAVS && !WRITEMHA && !WRITERAWIV) return GoToNextOptStep;


  // check if user has changed the file name that is supposed to be read in
  GetPot controlfile(dfltINIFile);
  std::string Default=FILEIN;  // use old value as default
  try
  {  // note if file not open will return default
     FILEIN=controlfile("mrti/filein",Default);
     PIXEL_SIZE=controlfile("mrti/pixel_size",1);
  }
  catch(const std::exception& e) //catch bad lexical cast 
  {
     std::cout << e.what() << std::endl;
     PetscPrintf(PETSC_COMM_WORLD,"CheckImageData: setting defaults\n");
     PIXEL_SIZE=1; FILEIN = Default;  // use old value as default
  }
  PetscPrintf(PETSC_COMM_WORLD, "mrtiImage PIXEL_SIZE=%d\n"   , PIXEL_SIZE );

  // check file status
  for(unsigned int kkk=0 ; kkk < ImagingData[0].depth ; kkk++){
      sprintf(file_i,FILEIN.c_str(),kkk+1,IDfile); 
      // on lonestar must wait until file exists or will cause 
      // an infinite loop trying to open once the file is there
      if(access(file_i,R_OK) == 0 ){
          struct stat status;
          stat( file_i, &status );
          if ( status.st_mode & S_IFDIR ){
              std::cout << "The file is a directory. Enter a file name." <<
std::endl;
              GoToNextOptStep = PETSC_FALSE ; break;
          }
          else { // found the file and it is not a directory. proceed
              printf("found %s \n", file_i);
          }
      } else { //path not found not ready for next opt step
          GoToNextOptStep = PETSC_FALSE ; break;
      }
      if(FILEIN.substr(FILEIN.find('.')+1) == "bin") { // raw binary data
          file.open (file_i, std::ios::in|std::ios::binary|std::ios::ate); // open the file 
          size = file.tellg(); // file size in bytes
          unsigned int ipixfil = size/PIXEL_SIZE; // # of pixels in file
          /* check that have the full file if not close 
             the file b/c not ready for next step*/
          if(ipixfil!= ImagingData[0].width * ImagingData[0].height   ){
             printf("error opening %s  \n",file_i);
             printf("found %d pixels expected %d pixels \n",
                       ipixfil, ImagingData[0].width * ImagingData[0].height);
             GoToNextOptStep = PETSC_FALSE ; 
          }
          file.close(); 
      }
  }
  if(!GoToNextOptStep){
     PetscPrintf(PETSC_COMM_WORLD, 
                 "%s not available sleeping for %d secs \n", file_i,isleep);
     sleep(isleep);
  }
  return GoToNextOptStep;
}
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
             read MRTI data as initial condition 
  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
PetscErrorCode RealTimeImage::loadMRTItempfield(AppSolve *user)
{

  qoiBaseClass   &qoiOptimizer = *user->qoiOptimizer; // QOI data

  // used for error handling
  PetscFunctionBegin; 

  PetscInt idMRTI = qoiOptimizer.IDEAL_NZERO - baseInfo::TIMELAG_IC;
  if( idMRTI < 0 ){ std::cout << "loadMRTItempfield" << std::endl;  libmesh_error();}
  PetscInt idFEM  = idMRTI * baseInfo::ISTEPS_PER_IDEAL;


  // Get a reference to the ExplicitSystem for ideal data
  TransientFEMSystem & ideal_system =
     user->_eqnSystem->get_system<TransientFEMSystem>("IdealSystem");
  // Get a reference to the NonlinearImplicitSystem we are solving
  TransientFEMSystem& state_system = 
    user->_eqnSystem->get_system<TransientFEMSystem>("StateSystem");

  // determine if file already in memory or needs to be read in
  if(MRTI_MEM.at( idMRTI )){
      /* data already there and stored in qoi do not need to read anything in */
  }else{
      /* Data Server will read in the data from disk and broadcast */
      PetscLogEventBegin(baseInfo::logevents[3],0,0,0,0); // read disk

      /* Data Server is always rank 1 in DDDAS_COMM_WORLD 
         must do intra-communicator Isend then inter-communicator 
         broadcast... there is no MPI_Ibcast */
      if( !libMesh::processor_id() ) {
         MPI_Status status;
         MPI_Recv(ImagingData[0].data,ImagingData[0].size(),MPI_FLOAT, 0, idMRTI ,user->DataComm,&status);
         printf("group %d received mrti time instance %d\n",
                                        user->GroupID, idMRTI );
      }
      // Broadcast locally (MPI-1 implementation)
      MPI_Bcast(ImagingData[0].data,ImagingData[0].size(),MPI_FLOAT, 0, PETSC_COMM_WORLD);
      // put MRTI thermal image into FEM data structures for later
      //  ***NOTE*** the thermal image at the TIMELAG is put into the
      //  ***NOTE***   initial time instance of the optimization step
      //  ***NOTE***   for optimal control this assumes that there is no
      //  ***NOTE***   heating over the duration of the timelag
      // put MRTI thermal image into FEM data structures
      ideal_system.project_solution(get_mrti_data,NULL,
                                      user->_eqnSystem->parameters);
      // store the solution history
      *ideal_system.vector_solution.at(idFEM)=
                                         *ideal_system.current_local_solution;
     

      PetscLogEventEnd(baseInfo::logevents[3],0,0,0,0); // read disk
      MRTI_MEM.at( idMRTI  ) = 1;
  }
  // copy to solution field
  *state_system.vector_solution.at(idFEM)= *ideal_system.vector_solution.at(idFEM);

  PetscFunctionReturn(0);
}
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
             read MRTI data for QOI
  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
PetscErrorCode Image::ReadThermalData(AppSolve *user)
{
  PetscErrorCode info;

  PetscFunctionBegin; 
  qoiBaseClass   &qoiOptimizer = *user->qoiOptimizer; // global data
  TransientFEMSystem*  pdeSolver = user->pdeSolver;  // abstract pde solver


  // determine which files are already in memory and which need to be read in
  for(PetscInt idMRTIfile  = qoiOptimizer.IDEAL_NZERO  ; 
               idMRTIfile <= qoiOptimizer.IDEAL_NTIME  ; idMRTIfile++)
   {
     if(MRTI_MEM.at(idMRTIfile)){
         /* if data is already there 
                  DO NOTHING                */
     }else{
        /* Data Server will read in the data from disk and broadcast */
        PetscLogEventBegin(baseInfo::logevents[3],0,0,0,0); // read disk

        // read the image data
        this->ReadImageData(idMRTIfile);
        printf("group %d received mrti time instance %d\n",
                                 user->GroupID,idMRTIfile);

        // Get a reference to the ExplicitSystem for ideal data
        TransientFEMSystem & ideal_system =
           user->_eqnSystem->get_system<TransientFEMSystem>("IdealSystem");
        // put MRTI thermal image into FEM data structures
        ideal_system.project_solution(get_mrti_data,NULL,
                                        user->_eqnSystem->parameters);

        // read the uncertainty data
        this->ReadUncertaintyData(idMRTIfile);
        printf("group %d received mrti uncertainty time instance %d\n",
                                                   user->GroupID,idMRTIfile);

        // Get a reference to the ExplicitSystem for ideal data
        TransientFEMSystem & ideal_uncertainty_system =
           user->_eqnSystem->get_system<TransientFEMSystem>("IdealUncertaintySystem");
        // put MRTI thermal image into FEM data structures
        ideal_uncertainty_system.project_solution(get_uncertainty_data,NULL,
                                        user->_eqnSystem->parameters);

        // apply dirichlet data if any
        if(pdeSolver->m_dirichletNodes.size())
          {
            //return probe temp 
            for( unsigned int Ii = 0; 
                              Ii<pdeSolver->m_dirichletNodes.size(); Ii++) 
              {
                 Point dum;
                 ideal_system.solution->set(pdeSolver->m_dirichletNodes[Ii],
              user->pdeSolver->InitValues(0,0,dum, user->_eqnSystem->parameters)
                                            );
                 // dirichlet data should not be weighted
                 //ideal_uncertainty_system.solution->set(
                 //                       pdeSolver->dirichletNodes[Ii],
                 //                                                    sigma);
              }
            ideal_system.solution->localize(
                                          *ideal_system.current_local_solution);
            //ideal_uncertainty_system.solution->localize(
            //                *ideal_uncertainty_system.current_local_solution);
          }

        // store the solution history
        PetscInt nstephi =  idMRTIfile     * baseInfo::ISTEPS_PER_IDEAL;
        *ideal_system.vector_solution.at(nstephi)=
                                           *ideal_system.current_local_solution;
        *ideal_uncertainty_system.vector_solution.at(nstephi)=
                          *ideal_uncertainty_system.current_local_solution;
       
        // interpolate the data onto the fem data structures
        if(idMRTIfile > qoiOptimizer.IDEAL_NZERO){
           PetscInt nsteplo = (idMRTIfile-1)  * baseInfo::ISTEPS_PER_IDEAL;
           for(int jjj = 1; jjj <  baseInfo::ISTEPS_PER_IDEAL ; jjj++)
            {
               PetscInt idfemfile = nsteplo + jjj ; 
               // mrti data
               ideal_system.vector_solution.at(idfemfile)->zero();
               ideal_system.vector_solution.at(idfemfile)->add(
                       (double)(baseInfo::ISTEPS_PER_IDEAL - jjj )/
                              (double)(baseInfo::ISTEPS_PER_IDEAL) ,
                           *ideal_system.vector_solution.at( nsteplo ) );
               ideal_system.vector_solution.at(idfemfile)->add(
                       (double)(                            jjj )/
                              (double)(baseInfo::ISTEPS_PER_IDEAL) ,
                           *ideal_system.vector_solution.at( nstephi ) ); 
               // associated uncertainty
               ideal_uncertainty_system.vector_solution.at(idfemfile)->zero();
               ideal_uncertainty_system.vector_solution.at(idfemfile)->add(
                       (double)(baseInfo::ISTEPS_PER_IDEAL - jjj )/
                              (double)(baseInfo::ISTEPS_PER_IDEAL) ,
                      *ideal_uncertainty_system.vector_solution.at( nsteplo ) );
               ideal_uncertainty_system.vector_solution.at(idfemfile)->add(
                       (double)(                            jjj )/
                              (double)(baseInfo::ISTEPS_PER_IDEAL) ,
                      *ideal_uncertainty_system.vector_solution.at( nstephi ) );
            }
        }

        PetscLogEventEnd(baseInfo::logevents[3],0,0,0,0); // read disk
        MRTI_MEM.at(idMRTIfile) = 1;
     }
   }
  PetscFunctionReturn(0);
}
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
PetscErrorCode Image::GetSpacing(PetscScalar* SpacingForFortran)
{
  PetscErrorCode info;
  // used for error handling
  PetscFunctionBegin; 
  SpacingForFortran[0] = spacing[0] ; 
  SpacingForFortran[1] = spacing[1] ; 
  SpacingForFortran[2] = spacing[2] ; 

  PetscFunctionReturn(0);
}
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
void Image::GetControlData(PetscInt Ntime)
{
  PetscErrorCode info;

/*
  no point of openmp here any more b/c this is trying to fork an 
      already forked thread (for the thermal images)
      which doesn't seem to be working

  omp_set_num_threads(NDATATHRD); 
#pragma omp parallel for private(file_i,file,size,tid,ipixfil,iii,jjj) schedule(dynamic)
*/
  // filtering parameters
  GetPot controlfile(dfltINIFile);
     
  //median filter neighborhood size
  try
  {  // note if file not open will return default
     MEDIANCROPSIZE= controlfile("mrti/mediancropsize",2 );
  }
  catch(const std::exception& e) //catch bad lexical cast 
  {
     std::cout << e.what() << std::endl;
     std::cout << "setting default MEDIANCROPSIZE\n"; MEDIANCROPSIZE=2;
  }
  //deriche smoothing coefficient
  try
  {  // note if file not open will return default
     DERICHESIGMA = controlfile("mrti/derichesigma",0.5 );
  }
  catch(const std::exception& e) //catch bad lexical cast 
  {
     std::cout << e.what() << std::endl;
     std::cout << "setting default DERICHESIGMA \n"; DERICHESIGMA=0.5;
  }
  // maximum temperature difference allowed 
  // between two successive thermal images
  try
  {  // note if file not open will return default
     MAXDIFF= controlfile("mrti/maxdiff",11.0);
  }
  catch(const std::exception& e) //catch bad lexical cast 
  {
     std::cout << e.what() << std::endl;
     std::cout << "setting default MAXDIFF \n"; MAXDIFF=11.0;
  }
  //cropping border
  try
  {  // note if file not open will return default
     IXLO= controlfile("mrti/ixlo",            0            );
     IXHI= controlfile("mrti/ixhi", (PetscInt) ImagingData[0].width );
     IYLO= controlfile("mrti/iylo",            0            );
     IYHI= controlfile("mrti/iyhi", (PetscInt) ImagingData[0].height);
  }
  catch(const std::exception& e) //catch bad lexical cast 
  {
     std::cout << e.what() << std::endl;
     std::cout << "not cropping use full image \n";
     IXLO = 0; IXHI = ImagingData[0].width -1;
     IYLO = 0; IYHI = ImagingData[0].height-1;
  }

  // array bounds check
  if(IXLO >= ImagingData[0].width  || IXLO < 0 ) IXLO = 0;
  if(IXHI >= ImagingData[0].width  || IXHI < 0 ) IXHI = ImagingData[0].width -1;
  if(IYLO >= ImagingData[0].height || IYLO < 0 ) IYLO = 0;
  if(IYHI >= ImagingData[0].height || IYHI < 0 ) IYHI = ImagingData[0].height-1;

  // transpose the data buffer
  // strcasecmp != 0   ==> strings don't match (case insensitive)
  if(controlfile("mrti/transpose",false)){ 
     TRANSPOSE = PETSC_TRUE ; 
  } else { 
     TRANSPOSE = PETSC_FALSE; 
  }

  //conversion factor
  try
  {  // note if file not open will return default
     CONVERSIONFACT= controlfile("mrti/conversionfact", 0.0);
  }
  catch(const std::exception& e) //catch bad lexical cast 
  {
     std::cout << e.what() << std::endl;
     std::cout << "using initial temperature for conversion\n"; 
     CONVERSIONFACT=0.0;
  }

  // echo input data 
  //int tid = omp_get_thread_num(); /* Obtain thread number */
  std::cout << ":ReadImageData opening time instance " << Ntime
       << "\n             IXLO= " << IXLO  << " IXHI= " <<  IXHI 
       << "\n             IYLO= " << IYLO  << " IYHI= " <<  IYHI 
       << "  sigma = "            << sigma 
       << "\n             DERICHESIGMA  =" << DERICHESIGMA
       <<             "   MEDIANCROPSIZE=" << MEDIANCROPSIZE
       <<             "   MAXDIFF       =" << MAXDIFF
       << "\n             TRANSPOSE     =" << TRANSPOSE
       <<             "   APPLICATORDIFF=" << APPLICATORDIFF 
       <<             "   CONVERSIONFACT=" << CONVERSIONFACT << std::endl;
  PetscFunctionReturnVoid();
}
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
PetscErrorCode RealTimeImage::checkMRTI_data4Pow(AppSolve *user)
{
  qoiBaseClass   &qoiOptimizer = *user->qoiOptimizer; // QOI data
  PetscErrorCode info;
  // used for error handling
  PetscFunctionBegin; 

  PetscInt  IDfile = qoiOptimizer.IDEAL_NTIME ;
  //  wait for the final image of this optimization step to be available
  //  before proceeding. This is so the group will wait to recieve the power
  //  profile and can calibrate from real-time treatment power update
  //  from the user
  PetscTruth GoToNextOptStep=PETSC_FALSE ;
  while(!GoToNextOptStep)
         GoToNextOptStep=CheckImageData(IDfile);

  PetscFunctionReturn(0);
}
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
// write data buffer file size  and coord data
void ImageServer::ImageSetupFile(){

 /* this routine is called by ONLY by visualization task */

 int ierr;
 MPI_File fh; // output file
 MPI_Status state; 

 /*----------buffer to write out avs .fld file----------------------*/
 if(WRITEAVS){
    const unsigned int sizmri    = size[0]*size[1]*(size[2]+1);
    const unsigned int sizcoord  = size[0]+size[1]+(size[2]+1);
    std::vector<float> coordbuffer(sizcoord,0.0);/* explicitly give the physical
                                                  coordinates for avs field */
    /* using rectilinear feild type in avs to explicitly give the coordinates
       of the point b/c uniform field type was not correctly doing this for
       when cropping */
    for(unsigned int i=0;i<size[0] ;i++) // x-coords
       coordbuffer.at(                             i ) = origin[0]+i*spacing[0];
    for(unsigned int i=0;i<size[1] ;i++) // y-coords
       coordbuffer.at( size[0]               + i ) = origin[1]+i*spacing[1];
    /* notice that avs buffer is padded by one in the z-dimension
       this is b/c for nslice = 1 avs crop breaks, must have at least 2 */
    for(unsigned int i=0;i<size[2] + 1 ;i++) // z-coords
       coordbuffer.at( size[0]+size[1] + i )  = origin[2]+i*spacing[2];
    if(BYTESWAP) swapByteOrder(&coordbuffer[0],coordbuffer.size()); 

    // allocate space and initialize avs buffer
    avsbuffer.reserve( (2*sizmri+sizcoord)*sizeof(float));  
    avsbuffer.resize( (2*sizmri+sizcoord)*sizeof(float),0);  

    // pack coordinates
    int coordposition= 2*sizmri*sizeof(float);  
    ierr=MPI_Pack(&coordbuffer[0],coordbuffer.size(),MPI_FLOAT,&avsbuffer[0],
                               avsbuffer.size(),&coordposition,PETSC_COMM_SELF);

    // assume that once the file has been written out to the buffer
    // size then the full file has been written
    std::ofstream sizefile; //output file
    sizefile.open (VIS_SIZEFILE.c_str(), std::ios::out);
    sizefile <<avsbuffer.size()<< std::endl;
    sizefile.close();
 }
 /*----------buffer to write out cvc .rawiv file-------------------*/
 if(WRITERAWIV){
    // cvc rawiv file data
    float spans[3],rawivorigin[3],bounds[6];
    unsigned int numVerts, numCells,dims[3];
    const unsigned int sizmri    = size[0]*size[1]*size[2];
    /* .rawiv header info */
    bounds[0]=origin[0]; 
    bounds[1]=origin[1]; 
    bounds[2]=origin[2]; 
    bounds[3]=origin[0] + size[0] * spacing[0];
    bounds[4]=origin[1] + size[1] * spacing[1];
    bounds[5]=origin[2] + size[2] * spacing[2];
    dims[0]=size[0]; 
    dims[1]=size[1]; 
    dims[2]=size[2]; 
    spans[0]=(bounds[3]-bounds[0])/(dims[0]-1);
    spans[1]=(bounds[4]-bounds[1])/(dims[1]-1);
    spans[2]=(bounds[5]-bounds[2])/(dims[2]-1);
    rawivorigin[0]=RAWORIGIN[0];
    rawivorigin[1]=RAWORIGIN[1];
    rawivorigin[2]=RAWORIGIN[2];
    numVerts = (size[0]  )*(size[1]  )*(size[2]  );
    numCells = (size[0]-1)*(size[1]-1)*(size[2]-1);
    // always big endian
    swapByteOrder(bounds,(unsigned int)6); 
    swapByteOrder(rawivorigin,(unsigned int)3);
    swapByteOrder(&numVerts,(unsigned int)1); 
    swapByteOrder(&numCells,(unsigned int)1);
    swapByteOrder(dims,(unsigned int)3); 
    swapByteOrder(spans,(unsigned int)3);
    // initialize position in buffer for packing
    IPOS=0;
    // allocate space and initialize .rawiv buffer
    rawivbuffer.reserve( 68+sizmri*sizeof(float));  
    rawivbuffer.resize( 68+sizmri*sizeof(float),0);  
    // pack header
    ierr=MPI_Pack(bounds,6,MPI_FLOAT,&rawivbuffer[0],
                  rawivbuffer.size(),&IPOS,PETSC_COMM_SELF);
    ierr=MPI_Pack(&numVerts,1,MPI_INT,&rawivbuffer[0],
                  rawivbuffer.size(),&IPOS,PETSC_COMM_SELF);
    ierr=MPI_Pack(&numCells,1,MPI_INT,&rawivbuffer[0],
                  rawivbuffer.size(),&IPOS,PETSC_COMM_SELF);
    ierr=MPI_Pack(dims,3,MPI_INT,&rawivbuffer[0],
                  rawivbuffer.size(),&IPOS,PETSC_COMM_SELF);
    ierr=MPI_Pack(rawivorigin,3,MPI_FLOAT,&rawivbuffer[0],
                  rawivbuffer.size(),&IPOS,PETSC_COMM_SELF);
    ierr=MPI_Pack(spans,3,MPI_FLOAT,&rawivbuffer[0],
                  rawivbuffer.size(),&IPOS,PETSC_COMM_SELF);
 }

 return ;
}
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*  read image data from ITK assume prefiltered data*/
void Image::ReadImageData(PetscInt Ntime)
{
  this->ReadSingleImage(Ntime,FILEIN,baseImage); 
}
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*  read image data from ITK assume prefiltered data*/
void Image::ReadUncertaintyData(PetscInt Ntime)
{
  this->ReadSingleImage(Ntime,UNCERTAINTYIN,sigmaImage); 
}
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*  read image data from ITK assume prefiltered data*/
void Image::ReadSingleImage(PetscInt Ntime, std::string &fileName, 
                                            ImageType::Pointer &refImage)
{
  PetscLogEventBegin(baseInfo::logevents[3],0,0,0,0); // read disk

  // write out for verification
  typedef float PixelType;
  typedef itk::Image<    PixelType   , 3 >  OutImageType;
  typedef itk::ImageFileWriter<             OutImageType > WriterType;

  switch( libMesh::processor_id() )
    {
     case 0:  // read on root then broadcast
      {
       // Create data buffer to work with for this time instance
       CImg<AVSDATA_TYPE> DataBuffer(ImagingData[0].data,
                                     size[0],size[1],size[2],1,true);
       this->GetControlData(Ntime);

       // read in data
       char file_i[MAXLEN];

       // file Name
       sprintf(file_i,fileName.c_str(),Ntime); 

       // open using ITK library of file types
       typedef itk::ImageSeriesReader< ImageType >   LocalReaderType;
       LocalReaderType::Pointer  reader = LocalReaderType::New();

       reader->SetFileName( file_i );
       try
         {
         reader->Update();

         // pointer to image in file
         ImageType::Pointer localImage = reader->GetOutput();

         #if defined(PETSC_USE_DEBUG)
         // write out for verification
         WriterType::Pointer writer = WriterType::New();
         // set the file output name
         OStringStream file_o;
         file_o << "mrivis/orig"<<fileName.substr(fileName.rfind('/') +1) ;
         // read in data
         char fileOut[MAXLEN];

         // file Name
         sprintf(fileOut,file_o.str().c_str(),Ntime); 
         writer->SetFileName( fileOut );
         std::cout << "writing " << fileOut << std::endl;
         writer->SetInput ( reader->GetOutput() );
         writer->Update();
         #endif

	 // Software Guide : BeginLatex setup real, imaginary, base phase, and
	 // temperature map iterators The const slice iterator walks the 3D
	 // input image, and the non-const linear iterator walks the 2D output
	 // image.  The iterators are initialized to walk the same linear path
	 // through a slice.  Remember that the \emph{second} direction of the
	 // slice iterator defines the direction that linear iteration walks
	 // within a slice

         typedef itk::ImageSliceIteratorWithIndex<ImageType>  LocalIteratorType;
         LocalIteratorType  tempIt(localImage,localImage->GetRequestedRegion()),
                            baseIt(  refImage,  refImage->GetRequestedRegion());
         
         tempIt.SetFirstDirection(  0 );  tempIt.SetSecondDirection( 1 );
         baseIt.SetFirstDirection(  0 );  baseIt.SetSecondDirection( 1 );
         
         // Software Guide : EndCodeSnippet 
         // set iterators to the beginning
         tempIt.GoToBegin();
         baseIt.GoToBegin();

         PetscPrintf(PETSC_COMM_WORLD,"Extracting Data...\n");

         /* loop through data structures  and store data */ 
         tempIt.GoToBegin(); baseIt.GoToBegin();
         while( !tempIt.IsAtEnd() )
           {
           const int kkk = tempIt.GetIndex()[2];
           while ( !tempIt.IsAtEndOfSlice() )
             {
             const int jjj = tempIt.GetIndex()[1];
             while ( !tempIt.IsAtEndOfLine() )
               {
               const int iii = tempIt.GetIndex()[0];
               DataBuffer(iii,jjj,kkk) = (AVSDATA_TYPE) baseIt.Get();
               //crop region of interest
               if(IYLO <= jjj && jjj<=IYHI && IXLO <= iii && iii<=IXHI)
                   DataBuffer(iii,jjj,kkk) = (AVSDATA_TYPE) tempIt.Get() ;
               ++tempIt;
               ++baseIt;
               }
             baseIt.NextLine();
             tempIt.NextLine();
             }
           baseIt.NextSlice();
           tempIt.NextSlice();
           }

         CHKMEMA; // check for memory corruption use -malloc_debug to enable

         std::cout << " time instance "<< Ntime << " in memory" << std::endl;
         }
       catch (itk::ExceptionObject &excp)
         {
         std::cerr << "Exception thrown" << excp << std::endl;
         std::cerr << "using default image" << std::endl;
         // Software Guide : BeginCodeSnippet
         typedef itk::ImageRegionIterator< ImageType >  RegionIteratorType;
         RegionIteratorType baseIt( refImage, refImage->GetRequestedRegion() );
         baseIt.GoToBegin(); 
         while ( !baseIt.IsAtEnd())
           {
             const int iii = baseIt.GetIndex()[0];
             const int jjj = baseIt.GetIndex()[1];
             const int kkk = baseIt.GetIndex()[2];
             DataBuffer(iii,jjj,kkk) = (AVSDATA_TYPE) baseIt.Get() ;
             ++baseIt; 
           }
         // Software Guide : EndCodeSnippet
         }

       #if defined(PETSC_USE_DEBUG)
       sprintf(file_i,"DataBuffer Ntime=%d",Ntime); 
       DataBuffer.print(file_i);
       try // write out for verification
         {
         WriterType::Pointer writer = WriterType::New();
         // set the file output name
         OStringStream file_o;
         file_o << "mrivis/CImg"<<fileName.substr(fileName.rfind('/') +1) ;
         // read in data
         char fileOut[MAXLEN];

         // file Name
         sprintf(fileOut,file_o.str().c_str(),Ntime); 
         writer->SetFileName( fileOut );
         std::cout << "writing " << fileOut << std::endl;
         writer->SetInput ( importFilter->GetOutput() );
         writer->Update();
         }
       catch (itk::ExceptionObject &excp)
         {
         std::cout << "Exception thrown" << excp << std::endl;
         }
       #endif

       if(TRANSPOSE) DataBuffer.transpose(); //data is transposed
       PetscLogEventEnd(baseInfo::logevents[3],0,0,0,0); // read disk

      }// don't break... continue and broadcast
     default:
       // Broadcast locally (MPI-1 implementation)
       MPI_Bcast(ImagingData[0].data, ImagingData[0].size(),
                                      MPI_FLOAT, 0, PETSC_COMM_WORLD);
       // need to update?
       importFilter->Update();
    }
  return;
}
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*  get image data from server*/
void RealTimeImage::ReadImageData(PetscInt Ntime )
{
  /* Data Server is always rank 1 in DDDAS_COMM_WORLD 
     must do intra-communicator Isend then inter-communicator 
     broadcast... there is no MPI_Ibcast */
  MPI_Status status;
  // Broadcast locally (MPI-1 implementation)
  //if( !libMesh::processor_id() )
  //MPI_Recv(ImagingData[Ntime].data,ImagingData[Ntime].size(),MPI_FLOAT,
  //                                0,idMRTIfile,user->DataComm,&status);
  MPI_Bcast(ImagingData[Ntime].data, ImagingData[Ntime].size(),
                                          MPI_FLOAT, 0, PETSC_COMM_WORLD);
  return;
}
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*  read image data slice by slice */
void ImageServer::ReadImageData(PetscInt Ntime)
{
  if(!WRITEAVS && !WRITEMHA && !WRITERAWIV) return;

  PetscLogEventBegin(baseInfo::logevents[3],0,0,0,0); // read disk

  // Create data buffer to work with for this time instance
  CImg<AVSDATA_TYPE> DataBuffer(ImagingData[Ntime].data,
                                size[0],size[1],size[2],2,true);
  this->GetControlData(Ntime);

  // read in data
  char file_i[MAXLEN];
  std::ifstream file;

  for(unsigned int kkk=0 ; kkk < DataBuffer.depth ; kkk++){
      sprintf(file_i,FILEIN.c_str(),kkk+1,Ntime); 
      if(FILEIN.substr(FILEIN.find('.')+1) == "bin") { // raw binary data
         // if made it to this point file should be ready to be read in w/o any
         // additional error checks, open at end to get file size
         file.open(file_i, std::ios::in|std::ios::binary|std::ios::ate); 
         std::ifstream::pos_type size = file.tellg(); // get file size in bytes
         file.seekg (0, std::ios::beg); //set to beginning
         /* allocate buffer space to read in data */
         if(PIXEL_SIZE == 1){// 8 bits
            CImg<char> raw(DataBuffer.width,DataBuffer.height,1,1);
            file.read (reinterpret_cast<char *>(raw.data), size);
            #if defined(PETSC_USE_DEBUG)
            raw.print(file_i);
            #endif
            /* store data and if needed, convert intensity value to degK */ 
            cimg_forXY(raw,iii,jjj) DataBuffer(iii,jjj,kkk,1)=
                               (AVSDATA_TYPE) raw(iii,jjj) + CONVERSIONFACT;
         }else if(PIXEL_SIZE == 2){// 16 bits
            CImg<short> raw(DataBuffer.width,DataBuffer.height,1,1);
            file.read (reinterpret_cast<char *>(raw.data), size);
            #if defined(PETSC_USE_DEBUG)
            raw.print(file_i);
            #endif
            //if (!cimg::endian()) cimg::endian_swap(raw.data,width*height);
            /* store data and if needed, convert intensity value to degK */ 
            cimg_forXY(raw,iii,jjj) DataBuffer(iii,jjj,kkk,1)=
                               (AVSDATA_TYPE) raw(iii,jjj) + CONVERSIONFACT;
         }else {
            printf("unknown pixel_size = %d\n",PIXEL_SIZE); abort();
         }
         file.close();
      } else { // open using CImg library of file types
         CImg<AVSDATA_TYPE> slicebuffer(file_i);
         /* store data and if needed, convert intensity value to degK */ 
         cimg_forXY(slicebuffer,iii,jjj) DataBuffer(iii,jjj,kkk,1)=
                      (AVSDATA_TYPE) slicebuffer(iii,jjj) + CONVERSIONFACT;
      }
      for(unsigned int iii=IXLO; iii<=IXHI; iii++) // crop region of interest
         for(unsigned int jjj=IYLO; jjj<=IYHI; jjj++)
                        DataBuffer(iii,jjj,kkk,0) = DataBuffer(iii,jjj,kkk,1);
  }
  std::cout << " time instance "<< Ntime << " in memory\n";
  #if defined(PETSC_USE_DEBUG)
  sprintf(file_i,"DataBuffer Ntime=%d",Ntime); 
  DataBuffer.print(file_i);
  #endif
  //data is transposed
  if(TRANSPOSE) DataBuffer.transpose();
  PetscLogEventEnd(baseInfo::logevents[3],0,0,0,0); // read disk
  /* apply filter by creating data structure with access
     to the memory in MRTI data. apply filters to this
     data structure. when the data structure is destroyed, the 
     shared = true parameter prevents the memory from being freed. 
     The filter is implemented as a series of 2d filters instead
     of 1 3d filter b/c the 3d filter decrease the maximum temperature
     value too much*/
  PetscLogEventBegin(baseInfo::logevents[18],0,0,0,0); // filter
  for(unsigned int kkk=0 ; kkk < DataBuffer.depth ; kkk++){
     CImg<AVSDATA_TYPE> FILTER2D( DataBuffer.ptr(0,0,kkk),DataBuffer.width,
                                  DataBuffer.height      ,     1      ,1,true);
     cimg_forXY(FILTER2D,iii,jjj) 
      if( abs(FILTER2D(iii,jjj) - 
                 ImagingData[Ntime-1](iii,jjj,kkk,0)) > MAXDIFF )
                      FILTER2D(iii,jjj)=ImagingData[Ntime-1](iii,jjj,kkk,0);
     FILTER2D.blur_median(MEDIANCROPSIZE);//MEDIANCROPSIZE <= 1 --> no filter
     FILTER2D.blur(DERICHESIGMA);//DERICHESSIGMA = 0.0 --> no filter
  }
  PetscLogEventEnd(baseInfo::logevents[18],0,0,0,0); // filter
  return;
}
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*  write image data */
PetscScalar ImageServer::WriteImageData(PetscInt ntime)
{
  char file_o[MAXLEN];
  MPI_Status state; 
  int mode = MPI_MODE_CREATE | MPI_MODE_WRONLY;
  int ierr;
  MPI_File fh; // output file
  std::ofstream avsfld; //output file

  PetscLogEventBegin(baseInfo::logevents[16],0,0,0,0); // data write
  AVSDATA_TYPE pmaxfilt=0.0;

  // write out avs format
  if(WRITEAVS){
     AVSDATA_TYPE pminfilt,pminorig,pmaxorig;
     // bounds for filtered data
     CImg<AVSDATA_TYPE> Filter( &ImagingData[ntime](0,0,0,0),size[0],
                                 size[1],size[2],1,true);
     pminfilt=ImagingData[ntime](0,0,0,0), pmaxfilt=pminfilt;
     cimg_for(Filter,ptr,AVSDATA_TYPE) {
       const AVSDATA_TYPE& a=*ptr;
       if (a<pminfilt) pminfilt=a;
       if (a>pmaxfilt) pmaxfilt=a;
     }
     // bounds for original data
     CImg<AVSDATA_TYPE> Original( &ImagingData[ntime](0,0,0,1),size[0],
                                     size[1],size[2],1,true);
     pminorig=ImagingData[ntime](0,0,0,1), pmaxorig=pminorig;
     cimg_for(Original,ptr,AVSDATA_TYPE) {
       const AVSDATA_TYPE& a=*ptr;
       if (a<pminorig) pminorig=a;
       if (a>pmaxorig) pmaxorig=a;
     }
     int nwidth,nheight,ndepth;
     nwidth=size[0]; nheight=size[1]; ndepth=size[2];
     // filename
     sprintf(file_o,FILEOUT.c_str(),ntime,"fld");
     std::cout << "writing " << file_o << std::endl;

     // open the file
     avsfld.open (file_o, std::ios::out);
     avsfld << "# AVS field file " << std::endl;
     avsfld << "ndim=3        # number of dimensions in the field" << std::endl;
     avsfld << "nspace=3      # no. of physical coordinates per point" << std::endl;
     avsfld << "dim1= " << nwidth  << "     # dimension of axis 1"<< std::endl;
     avsfld << "dim2= " << nheight << "     # dimension of axis 2"<< std::endl;
    /* notice that avs buffer is padded by one in the z-dimension
       this is b/c for nslice = 1 avs crop breaks, must have at least 2 */
     avsfld << "dim3= " << ndepth+1<< "     # dimension of axis 3"<< std::endl;
     avsfld << "veclen=2        # number of components at each point" << std::endl;
     avsfld << "label=filtered-mrti        # data labels" << std::endl;
     avsfld << "label=original-mrti        # data labels" << std::endl;
     avsfld << "unit =kelvin               # data units"  << std::endl;
     avsfld << "unit =kelvin               # data units"  << std::endl;
     avsfld << "min_val = "<<pminfilt<<" "<<pminorig<< " # " << std::endl;
     avsfld << "max_val = "<<pmaxfilt<<" "<<pmaxorig<< " # " << std::endl;
     avsfld << "data=float       # native format binary" << std::endl;
     avsfld << "field=rectilinear # field type" << std::endl;
     avsfld <<'\f'<<'\f';
     avsfld.close();

     // write out pixel data and coordinate data
     avsfld.open (file_o, std::ios::out|std::ios::binary|std::ios::app);
     // pack the data
     int ipos = 0 ;
     cimg_forXYZ( ImagingData[ntime],iii,jjj,kkk){
        ierr=MPI_Pack(&ImagingData[ntime](iii,jjj,kkk,0),1,MPI_FLOAT,
                                          &avsbuffer[0],avsbuffer.size(),
                                          &ipos,PETSC_COMM_SELF);
        ierr=MPI_Pack(&ImagingData[ntime](iii,jjj,kkk,1),1,MPI_FLOAT,
                                          &avsbuffer[0],avsbuffer.size(),
                                          &ipos,PETSC_COMM_SELF);
     }
     if(BYTESWAP)swapByteOrder((float*)&avsbuffer[0],ImagingData[ntime].size());
     avsfld.write (reinterpret_cast<char *>(&avsbuffer[0]),avsbuffer.size());
     avsfld.close();
  }
  // write out mha format
  if(WRITEMHA){
     // useful typedef's
     //typedef itk::Vector<   float       , 2 >  PixelType;
     typedef float PixelType;
     typedef itk::Image<    PixelType   , 3 >  OutImageType;
     typedef itk::ImageSliceIteratorWithIndex< OutImageType > OutIterType;
     typedef itk::ImageFileWriter<             OutImageType > WriterType;
     // allocate memory for the output image
     OutImageType::RegionType region;
     region.SetSize(  size  );
     OutImageType::Pointer outputImage = OutImageType::New();
     outputImage->SetRegions( region );
     outputImage->Allocate();
     // set image dimensions
     outputImage->SetSpacing(spacing);
     outputImage->SetOrigin(origin);
     /*  set fastest iteration direction to the x direction
         and second fastest iteration direction to the y direction
                   ________
                  |->|->|->|   end 
                   --------
                  |->|->|->|        
          y        --------         
           ^      |->|->|->|        
           |       --------         
           |      |->|->|->|        
           |       --------
           | begin|->|->|->|  
                  |--------|   ----------> x                        */
     OutIterType    outputIt( outputImage, outputImage->GetRequestedRegion() );
     outputIt.SetFirstDirection(  0 );
     outputIt.SetSecondDirection( 1 );

     // Software Guide : EndCodeSnippet
     int kkk=0; outputIt.GoToBegin();
     while( !outputIt.IsAtEnd() )
       { int jjj=0;
       while ( !outputIt.IsAtEndOfSlice() )
         { int iii=0;
         while ( !outputIt.IsAtEndOfLine() )
           {
           OutImageType::PixelType   pixelValue;
           pixelValue =  ImagingData[ntime](iii,jjj,kkk,0) ; // filtered
           //pixelValue[0] =  ImagingData[ntime](iii,jjj,kkk,0) ; // filtered
           //pixelValue[1] =  ImagingData[ntime](iii,jjj,kkk,1) ; // unfiltered
           outputIt.Set( pixelValue ) ; 
           iii++;
           ++outputIt;
           }
         outputIt.NextLine(); jjj++;
         }
       outputIt.NextSlice(); kkk++;
       }

     WriterType::Pointer writer = WriterType::New();
     sprintf(file_o,FILEOUT.c_str(),ntime,"mha");
     writer->SetFileName( file_o );
     std::cout << "writing " << file_o << std::endl;
     writer->SetInput(outputImage);
     try
       {
       writer->Update();
       }
     catch ( itk::ExceptionObject &err)
       {
       std::cout << "ExceptionObject caught !" << std::endl; 
       std::cout << err << std::endl; 
       //return -1;   
       }
  }
  // write out rawiv format
  if(WRITERAWIV){
     int ipos = IPOS; 

     // pack the data
     ierr=MPI_Pack(ImagingData[ntime].data,ImagingData[ntime].size()/2,
                   MPI_FLOAT,&rawivbuffer[0],rawivbuffer.size(),&ipos,PETSC_COMM_SELF);
     // assume working on little endian machine
     swapByteOrder( (float*) &rawivbuffer[IPOS],ImagingData[ntime].size()/2);

     // open the file
     sprintf(file_o,FILEOUT.c_str(),ntime,"rawiv");
     std::cout << "writing " << file_o << std::endl;
     ierr=MPI_File_open(PETSC_COMM_SELF,file_o,mode,MPI_INFO_NULL,&fh);

     // write the file
     ierr=MPI_File_write(fh,&rawivbuffer[0],rawivbuffer.size(),MPI_BYTE,&state);

     // close the file
     ierr=MPI_File_close(&fh) ;
  }
  
  PetscLogEventEnd(baseInfo::logevents[16],0,0,0,0); // data write
  PetscScalar MaxTemp = (PetscScalar) pmaxfilt;
  return MaxTemp;
}

/* -------------------------------------------------------------------- */
void ImageServer::LaserImageControl(PetscInt &IDimgfile,
                       std::vector<qoiBaseClass> *qoiOptimizer,
                       std::vector<MPI_Comm> *DataComm,
                       std::vector<MPI_Request> *Send_Request){
  MPI_Request send_request;

  if(CheckImageData(IDimgfile)){
     if(IDimgfile > baseInfo::MRTI_nzero) ReadImageData(IDimgfile);
   
     //for(unsigned int iii=0; iii < qoiOptimizer->size() ; iii++){
     //  if(Need_to_send.at(IDimgfile).at(iii)){//group needs data
     //     // compute groups only have a buffer for the filtered data
     //     // but pixel buffer stored continuously in memory the first have
     //     // of the buffer should be the filtered data
     //     int compute_group_buffer_size = ImagingData[IDimgfile].size()/2 ;
     //     MPI_Isend( ImagingData[IDimgfile].data,compute_group_buffer_size,
     //                 MPI_FLOAT,0,IDimgfile,DataComm->at(iii),&send_request);
     //     // build a list of send request
     //     Send_Request->push_back(send_request);
     //  }
     //}
     PetscScalar MaxTemp=WriteImageData(IDimgfile);
     // when an image representing t_i is sent the power 
     // that should be output for [t_i,t_{i+1}) is immediately 
     // retrieved by the MRTI server. In order for the immediate
     // retrieval of the power, the power file for [t_i,t_{i+1}) must
     // be written before the the image at t_i is sent. Hence, after
     // reading and writing the image at t_i, the power for 
     // [t_{i+1},t_{i+2}) is written so it can already be there ready
     // to be retrieved by the MRTI server. NOTE that the power stored
     // in module constit_data in POW is of the form
     //
     //               line 1     line 2 ... line n-1       line n
     //   |----P(t1)----|---P(t2)---|----------|-----P(tn)-----|
     //  t=0           t1          t2   ...   tn-1             tn
     //
     //  if the image just read in represents time t_i then the
     //  power that should be output is for [t_i,t_{i+1}) which
     //  is stored under index IDimgfile={i+1}
     IDimgfile= IDimgfile + 1 ; 
     PetscInt IDpow = (IDimgfile+1); 
     //if(IDpow <= baseInfo::MRTI_ntime)
     //    FORTRAN_NAME(write_visualase_file)(&IDpow,&MaxTemp);
  }
}
} /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
                         end dddas namespace
         since calling c++ class member functions may not be 
         portable, i.e.

     http://www.parashift.com/c++-faq-lite/pointers-to-members.html

         below are wrappers to call the above class...
   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* Function to get MRTI data */
Number get_mrti_data (const Point& p,
                      const Parameters& parameters,
                      const std::string& sys_name,
                      const std::string& unknown_name)
{
  
  AppSolve *user = parameters.get<AppSolve*> ("AppSolve");
  TransientFEMSystem*  pdeSolver = user->pdeSolver;  // abstract pde solver

  InterpolatorType::PointType point;
  point[0] = p(0);
  point[1] = p(1);
  point[2] = p(2);

  Number temp;
  PetscScalar dum;

  if( user->MRTI->interpolator->IsInsideBuffer(point))
   {
     temp =  user->MRTI->interpolator->Evaluate( point ); 
   }
  else
   {
     temp =  pdeSolver->InitValues(0,0,p,parameters);
   }

  return temp;
}
/* Function to get MRTI uncertainty */
Number get_uncertainty_data (const Point& p,
                             const Parameters& parameters,
                             const std::string& sys_name,
                             const std::string& unknown_name)
{
  AppSolve *user = parameters.get<AppSolve*> ("AppSolve");

  InterpolatorType::PointType point;
  point[0] = p(0);
  point[1] = p(1);
  point[2] = p(2);

  Number temp;

  if( user->MRTI->interpolator->IsInsideBuffer(point))
   {
     temp =  user->MRTI->interpolator->Evaluate( point ); 
   }
  else
   {
     temp =  user->MRTI->sigma ; 
   }

  return temp;
}
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
PetscTruth checkMRTI_data4OPT(void *ctx)
{
  AppSolve *user = (AppSolve *) ctx;  
  qoiBaseClass   &qoiOptimizer = *user->qoiOptimizer; // QOI data
  dddas::Image   &MRTI    = *user->MRTI;    // image data
  PetscInt   IDfile;
  PetscTruth GoToNextOptStep=PETSC_FALSE ;
  PetscErrorCode info;

  IDfile = qoiOptimizer.IDEAL_NTIME + qoiOptimizer.NUMIDEAL ;
  // if the next time window is ready to be read in, then 
  //   kill this optimization step. NUMIDEAL is used to check for the
  //   next time window and may be used to:
  //        1) finish the optimization completely if IDfile will never exist.
  //           i.e. IDfile >  MRTI_NTIME 
  //        2) move to next optimization window if IDfile detected
  //           i.e. IDfile <= MRTI_NTIME 
  GoToNextOptStep=MRTI.CheckImageData(IDfile);

  return  GoToNextOptStep;
}
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
PetscTruth checkMRTI_data4_IC(void *ctx)
{
  AppSolve      *user = (AppSolve *) ctx;  
  qoiBaseClass  &qoiOptimizer = *user->qoiOptimizer; // QOI data
  dddas::Image  &MRTI    = *user->MRTI; // image data
  PetscTruth GoToNextOptStep=PETSC_FALSE ;
  PetscErrorCode info;

  PetscInt   IDfile = qoiOptimizer.IDEAL_NZERO + qoiOptimizer.TIMELAG ;
  // This assumes that the optimal control computations have 
  //   been given a head start through the use of TIMELAG_IC
  //   when the image IDfile has been detected the image representing
  //   the beginning of the time window has been detected
  //   finish and plot. Notice that the computations will run to completion
  //   if IDfile is never detected
  GoToNextOptStep=MRTI.CheckImageData(IDfile);

  return  GoToNextOptStep;
}
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     read MRTI data for objective function evaluations   
  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
PetscErrorCode readMRTI_data( AppSolve *user ){
  dddas::Image &MRTI = *user->MRTI; // image data
  PetscErrorCode info;

  PetscFunctionBegin; 
  info = MRTI.ReadThermalData(user);
  PetscFunctionReturn(0);
}

