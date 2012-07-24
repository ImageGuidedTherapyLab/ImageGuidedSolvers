//prevent multiple inclusions of header file
#ifndef READ_MRIDATA_H
#define READ_MRIDATA_H

#include "mri_params.h" // mri params
// ITK include files
#include "itkOrientedImage.h"
#include "itkImportImageFilter.h"
//#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "baseInfo.h"


// global typedefs
typedef itk::Image< AVSDATA_TYPE, Dimension > ImageType;
//typedef itk::NearestNeighborInterpolateImageFunction<
typedef itk::LinearInterpolateImageFunction<
                                 ImageType, double >  InterpolatorType;
typedef itk::ImportImageFilter< AVSDATA_TYPE, Dimension >   ImportFilterType;
typedef std::vector<PetscTruth> NeedData;

namespace dddas
{
// base case image with dicom header info
// used to manage file global CImg buffer
// base class is for retrospective data assumed that the data is already
// filtered
class Image{

// allow ImageServer to share all data
friend class ImageServer; 

public : 
   Image(GetPot &controlfile) ; //default constructor
   void ImageSetup(AppSolve &); //setup main data structure
   void SetImagePointers(AppSolve &); //setup pointers for use with ITK

   // check file availability
   virtual PetscTruth  CheckImageData(PetscInt)
     { libmesh_error(); return PETSC_FALSE; }

   // control code execution
   virtual PetscErrorCode checkMRTI_data4Pow(AppSolve *)
     { libmesh_error(); return 0; }

   // load thermal images as IC
   virtual PetscErrorCode loadMRTItempfield(AppSolve *)
     { libmesh_error(); return 0; }

   // load thermal images into QOI
   PetscErrorCode ReadThermalData(AppSolve *);
 
   // turn off image write if images not setup
   void NoWrite();

   // ITK data structures
   InterpolatorType::Pointer interpolator; 
   ImportFilterType::Pointer importFilter;
   // This filter requires the user to specify the size of the image to be
   // produced as output.  The \code{SetRegion()} method is used to this end.
   // The image size should exactly match the number of pixels available in the
   // locally allocated buffer. 
   ImportFilterType::SizeType  size;

   PetscErrorCode GetSpacing(PetscScalar*);
   PetscScalar    sigma;// default uncertainty in image
protected : 
   void GetControlData(PetscInt);

   // read Image functions
   virtual void ReadImageData(PetscInt);
   virtual void ReadUncertaintyData(PetscInt);
   void ReadSingleImage(PetscInt, std::string &, ImageType::Pointer &); 

   // physical space dimensions
   ImageType::SpacingType spacing;
   ImageType::PointType    origin;

   // file name
   std::string FILEIN,UNCERTAINTYIN; 

   PetscTruth TRANSPOSE;// logical variable to transpose the data buffer

   // keep track of which mrti data is in memory
   std::vector<PetscInt> MRTI_MEM; 

   PetscTruth WRITEAVS, // logical variable to determine if avs file is written
              WRITEMHA, // determines if mha file is written
              WRITERAWIV;// determines if rawiv file is written
   unsigned int PIXEL_SIZE;

   // filtering data
   int MEDIANCROPSIZE;
   float DERICHESIGMA;
   AVSDATA_TYPE MAXDIFF; // used for filtering, this is the maximum allowed 
                         // temperature difference between two 
                         // successive thermal images
   unsigned int IXLO,IXHI,IYLO,IYHI;
   /* convert temp images to degK, used for both absolute and relative 
      temperature images. 
      absolute ==> absolute degC given in the temperature images
      relative ==> delta Temp (in degC) given in the temperature images */
   AVSDATA_TYPE CONVERSIONFACT,APPLICATORDIFF;
private : 
   ImageType::Pointer baseImage,sigmaImage;
};

// derived class of dddas::Image used for real time
// assume that the recieved data is already filtered
class RealTimeImage: public Image
{
public : 
   //default constructor
   RealTimeImage(GetPot &controlfile):Image(controlfile){};

   // check file availability
   virtual PetscTruth  CheckImageData(PetscInt);

   // control code execution
   virtual PetscErrorCode checkMRTI_data4Pow(AppSolve *);

   // load thermal images as IC
   virtual PetscErrorCode loadMRTItempfield(AppSolve *);

private : 
   // read Image function
   virtual void ReadImageData(PetscInt );
};

// derived class of dddas::Image used to act as an image server
// used to manage file global CImgList buffer
class ImageServer : public Image
{
public : 
 ImageServer(GetPot &controlfile,std::vector<qoiBaseClass> &);

 void ServerSetup(AppSolve &,std::vector<qoiBaseClass> &); //setup data structures

 // write data buffer file size  and coord data
 void ImageSetupFile();
 // write Image function
 PetscScalar WriteImageData(PetscInt);
 // control timing of Image Read Laser parameter write 
 void LaserImageControl(PetscInt &,  std::vector<qoiBaseClass> *,
                        std::vector<MPI_Comm> *, std::vector<MPI_Request> *);
private : 

   // read Image function
   virtual void ReadImageData(PetscInt);
   // byte swapping
   PetscTruth BYTESWAP; // intel&amd =little endian   SUN=bid endian
   //origin used for rawiv data
   float RAWORIGIN[3];
   // store parameters to read/write in the mri/MRTI files
   std::string FILEOUT, VIS_SIZEFILE;
   std::vector<unsigned char> avsbuffer,rawivbuffer;// data that gets written to disk
   int IPOS, // stores ending position of header within the buffer
       NDATATHRD;// number of threads to use to read in mrti data slices
   // Create array that contains precomputed groups that need data
   std::vector<NeedData > Need_to_send;
};

} // end dddas namespace

PetscErrorCode readMRTI_data(     AppSolve *);
PetscTruth     checkMRTI_data4OPT(void *);
PetscTruth     checkMRTI_data4_IC(void *);
Number         get_mrti_data (       const Point&, const Parameters&,
                                     const std::string&, const std::string& );
Number         get_uncertainty_data (const Point&, const Parameters&,
                                     const std::string&, const std::string& );
#endif
