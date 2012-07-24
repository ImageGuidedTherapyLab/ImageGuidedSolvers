#ifndef __realtimeimagingImaging_h
#define __realtimeimagingImaging_h

// itk includes
#include "itkExtractImageFilter.h"
#include "itkAtan2ImageFilter.h"
#include "itkBinaryMagnitudeImageFilter.h"
#include "itkScalarToArrayCastImageFilter.h"
#include "itkImageSliceConstIteratorWithIndex.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkImportImageFilter.h"

// itk includes
#include <itksys/SystemTools.hxx>
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"

// local includes
#include "tttkUtilities.h"


/** default pixel type is double */
typedef PetscScalar                                              InputPixelType;
typedef float                                                   OutputPixelType;
const unsigned int          Dimension = 3;
typedef itk::Image< InputPixelType, Dimension >                  InputImageType;
typedef itk::ImageSeriesReader<    InputImageType >             InputReaderType;

const PetscInt nvarplot = 1;
//typedef itk::Image< itk::Vector< OutputPixelType , nvarplot > , Dimension > 
typedef itk::Image< OutputPixelType , Dimension >          VecOutputImageType;
typedef itk::Image< OutputPixelType, Dimension >              OutputImageType;
typedef itk::ImageFileWriter<    VecOutputImageType >           VecWriterType;
//typedef itk::ImageSliceIteratorWithIndex<OutputImageType > OutputIteratorType;
typedef itk::ImageFileWriter<       OutputImageType >              WriterType;

//template< class FilterType >
//InputImageType::Pointer 
template <class TPixel , class FilterType>
typename itk::Image< TPixel, Dimension >::Pointer 
GetImageData( FilterType ITKPointer, 
              const PetscInt rank,  std::vector<std::string>  &filenames )
{
  // Software Guide : BeginLatex
  // 
  // We can trigger the reading process by calling the \code{Update()}
  // method on the series reader. It is wise to put this invocation inside
  // a \code{try/catch} block since the process may eventually throw
  // exceptions.
  //
  // Software Guide : EndLatex
     
  bool DataNotRead = true;
  while(DataNotRead)
    {
     for (unsigned int i = 0; i < filenames.size(); i++)
      {
       if(access(filenames[i].c_str(),R_OK))
        {
         if(!rank) std::cout << filenames[i] << " NOT FOUND" << std::endl;
        }
       else 
        {
         //if(!rank) std::cout << filenames[i].substr( 
         //          filenames[i].rfind('/')+1, std::string::npos ) 
         if(!rank) std::cout << filenames[i]
                             << " found...." << std::endl;
        }
      }
     try
       {
       ITKPointer->Update();
       DataNotRead = false;
       }
     catch (itk::ExceptionObject &excp)
       {
       std::cout << "Exception thrown" << excp << std::endl;
       sleep(1);
       }
    }
  std::cout << std::flush ;
  return ITKPointer->GetOutput();
}

// useful typedefs
typedef itk::BinaryMagnitudeImageFilter<InputImageType,InputImageType,InputImageType>
                                                          MagnitudeFilterType;
typedef itk::ExtractImageFilter< InputImageType, InputImageType> 
                                                               ProcFilterType;
typedef itk::ImageSliceConstIteratorWithIndex< InputImageType >  
                                                            InputIteratorType;
typedef itk::ImageSliceIteratorWithIndex< InputImageType > BufferIteratorType;
typedef itk::ImageRegionIterator< InputImageType >         RegionIteratorType;

typedef itk::InterpolateImageFunction<
                             InputImageType , InputPixelType>  InterpolatorType;
typedef itk::LinearInterpolateImageFunction<
                             InputImageType , InputPixelType>  LinearInterpolatorType;
typedef itk::NearestNeighborInterpolateImageFunction<
                             InputImageType , InputPixelType>  NearestInterpolatorType;

/**
 * Input image acquistion
 */
enum ImageAcquisitionType {DICOM2D         =0,       
	                   VOLUME3D        =1,      
	                   BACKGROUNDPHASE =2,      
	                   EIGHT_ECHO_CSI  =3,  
	                   SIXTEEN_ECHO_CSI=4,
                           NO_IMAGING  };


// Abstract base class data structure for Kalman Filtering
//  Notice that ABC is not templated. ITK follows a similar design pattern 
class Imaging
{
public:
  Imaging(const GetPot &controlfile); //constructor

  void DebugOn(); // enable debugging

  // update imaging data structure container
  virtual PetscErrorCode UpdateImaging(const int) 
   { libmesh_not_implemented(); return 0;} 

  PetscErrorCode UpdateImageStats(const int); 

  // tmap mask
  void GetImageMask( );

  /** 
   * generate the expected file of a given image acquisition set for a 2d
   * Series of DICOM DATA
   */
  PetscErrorCode RealTimeGenerateFileNames( const int , const int , 
                              std::vector < std::vector < std::string > > &);

  /** generate the expected file of a single 3d image file */
  PetscErrorCode GetImageFileName( const int , std::string, std::string &);

  /** 
    * generate the expected file of a single 3d image file 
    *  possibly change file location from command line
    */
  PetscErrorCode GetImageFileNameFromCmdLine( const int , std::string, std::string &);

  /** get dicom header info */
  virtual PetscErrorCode GetHeaderData(ImageAcquisitionType ); 

  /** set header info */
  virtual PetscErrorCode SetHeaderData(int , int , int ,
				       double , double , double , 
                                       double , double , double  );

  /** project imaging to system */
  PetscErrorCode ProjectImagingToFEMMesh (System *, Vec , double );
  
  // filtering parameters
  InputImageType::SizeType zeroFilterRadius, medianFilterRadius;

  // to hold image and variance
  InputImageType::Pointer   net_Image;
  InputImageType::Pointer   var_Image;
  InputImageType::Pointer   snr_Image;
  InputImageType::Pointer   arr_Image;

  /** write an ini file containing dicom header info */
  PetscErrorCode WriteIni();
 
  /** echo local data */
  void printSelf(std::ostream& );

  PetscScalar GetTmapFactor() 
   {
    return 1.0 / ( 2.0 * libMesh::pi * m_imagfreq 
                       * m_alpha * m_echotime * 1.e-3 );
   }

  PetscScalar MagneticPotentialDenominator() 
   {
    return 1.0 / ( 2.0 * libMesh::pi * m_alpha * 1.5 );
   }

  PetscScalar GetEchoTime(){return m_echotime;}
  PetscScalar GetImagFreq(){return m_imagfreq;}
  PetscInt    get_size(PetscInt iI){return this->size[iI];}
  PetscScalar get_orgn(PetscInt iI){return this->orgn[iI];}
  PetscScalar get_space(PetscInt iI){return this->sp[ iI];}

  // return debug status
  PetscTruth Debug(){return m_Debug;}

  // magnitude file name
  std::string GetMagnitudeImageFileName(){ return MagnFile ;}

  typedef itk::ImportImageFilter< InputPixelType,Dimension >  ImportFilterType;

  // ITK data for interpolation
  InterpolatorType::Pointer m_interpolator; 

  // ITK data for importing data
  ImportFilterType::Pointer m_importFilter;
private:
  // image statistics
  InputImageType::Pointer   meanImage;

  /* counters for when to stop collecting 
     updates of magnitude and temperature statistics */
  PetscInt statCount, maxStat;

  // initial measurement data covariance [deg C]
  // @todo {is this used anywhere important... could prob delete}
  PetscScalar m_meascov ;   

protected:
  // dimensions for all images
  InputImageType::SpacingType sp;
  InputImageType::PointType orgn;
  InputImageType::RegionType::SizeType  size;
  InputImageType::RegionType::SizeType  n_roi;
  InputImageType::RegionType::IndexType index;
  /** initial temp should be redundant w/ pennes model */ 
  PetscScalar m_bodyTemp ;

  // tmap factor
  PetscScalar m_imagfreq, m_alpha , m_echotime;

  // acquisition time step 
  PetscScalar m_acquisition_dt;

  // for spacing if needed
  std::string MagnFile ;

  /** output directory */
  const char *OutputDir;
  // image location
  const char *ExamPath;
  int DirId;   

  //magnitude image intensity noise... (measured in an ROI in air)
  PetscScalar noise; 

  PetscInt necho,       ///< # of echo
           noffset,     ///< initial image offest 
           median;      ///< median filter

  PetscTruth m_Debug;  // control debuggin

  PetscScalar maxdiff; // IC

  // container for filenames
  std::vector< std::vector < std::string > > filenames;

  // to hold image mask
  InputImageType::Pointer   maskImage;

  // arrhenius damage
  PetscScalar  m_ActivationEnergy, m_FrequencyFactor,
               m_GasConstant     , m_BaseTemperature; 

  PetscInt m_nzero;       ///< # to begin imaging
};
/**
 * default derived class for directly reading volumetric 3d imaging data. 
 * the class reads in images of the instantiated pixel dimension
 */ 
template< unsigned int VPixelDimension = 1 >
class ImagingDirect : public Imaging 
{
public:
  //constructor
  ImagingDirect(const GetPot &controlfile);

  /** base class */
  typedef Imaging Parent;

  // update imaging data structure container
  virtual PetscErrorCode UpdateImaging(const int);

  // update imaging data structure container
  virtual PetscErrorCode GetHeaderData(ImageAcquisitionType ); 

  // class dependent typedefs
  typedef itk::Image< itk::Vector< PetscScalar, VPixelDimension >, Dimension > 
                                                                 ClassImageType;
  typedef itk::ImageSeriesReader< ClassImageType >              ClassReaderType;
  typedef itk::ScalarToArrayCastImageFilter< InputImageType, ClassImageType >
                                                                ClassFilterType;
  typedef itk::ImageSliceConstIteratorWithIndex< ClassImageType > 
                                                              ClassIteratorType;

protected:
  // store base Image
  typename ClassImageType::Pointer  baseImage;  

private:
  // file reader
  typename ClassReaderType::Pointer  imageReader; 

};
// constructor
#undef __FUNCT__
#define __FUNCT__ "ImagingDirect::ImagingDirect"
template< unsigned int VPixelDimension >
ImagingDirect< VPixelDimension  >::ImagingDirect( const GetPot &controlfile) 
    : Imaging( controlfile)
{

  // allocate file reader
  imageReader = ClassReaderType::New();

}
//get dicom header data
#undef __FUNCT__
#define __FUNCT__ "ImagingDirect::GetHeaderData"
template< unsigned int VPixelDimension >
PetscErrorCode ImagingDirect< VPixelDimension  >::
GetHeaderData(ImageAcquisitionType inputmethod)
{
  PetscFunctionBegin;
 
  // call base class first
  this->Parent::GetHeaderData(inputmethod);

  InputImageType::RegionType    fullRegion;
  fullRegion.SetSize(  this->size  );

  baseImage = ClassImageType::New();
  baseImage->SetRegions( fullRegion );
  baseImage->Allocate();
  baseImage->SetSpacing( this->sp  );
  baseImage->SetOrigin(  this->orgn);

  PetscFunctionReturn(0);
}
// update the imaging data structure and snr maps at this time instance
#undef __FUNCT__
#define __FUNCT__ "ImagingDirect::UpdateImaging"
template< unsigned int VPixelDimension >
PetscErrorCode ImagingDirect< VPixelDimension >::UpdateImaging(const int istep)
{
  PetscFunctionBegin;

  // get file name
  this->GetImageFileName(istep,std::string("inputImage."),filenames[0][0]);
  imageReader->SetFileName(filenames[0][0]);

  // update reader
  baseImage = GetImageData< typename itk::Vector< PetscScalar, VPixelDimension >, 
                            typename ClassReaderType::Pointer > (imageReader,
                                      libMesh::processor_id(), filenames[0]);
  PetscFunctionReturn(0);
}
/**************************************************************************** 
   derived class for complex images as the input
 ***************************************************************************/
class ImagingComplex : public ImagingDirect < 2 > 
{

public:
  //constructor
  ImagingComplex(const GetPot &controlfile)
    : ImagingDirect< 2 >( controlfile ){}

  /** base class */
  typedef ImagingDirect < 2 > Parent;

  /** get dicom header info */
  virtual PetscErrorCode GetHeaderData(ImageAcquisitionType ); 

  // update the tmap and snr maps at this time instance
  virtual PetscErrorCode UpdateImaging(const int); 

protected:
  // snr from magnitude
  void UpdateSNR(  InputImageType::Pointer &, const int );

  // get Noise data
  void ExtractNoise(const int );

  // generate phase images from real imaginary
  ClassImageType::Pointer GetPhaseImage(const int,ClassFilterType::Pointer);


};

/**************************************************************************** 
   derived class for input of tmaps and snr images 
 ***************************************************************************/
class ImagingPreProcessed : public ImagingDirect < 1 > 
{
public:
  //constructor
  ImagingPreProcessed(const GetPot &controlfile);

  // update the tmap and snr maps at this time instance
  virtual PetscErrorCode UpdateImaging(const int); 

};
/**************************************************************************** 
   csi imaging: # pixel dimensions = # of echos
 ***************************************************************************/
template< unsigned int VPixelDimension  >
class ImagingChemicalShift : public ImagingDirect < VPixelDimension > 
{
public:
  //constructor
  ImagingChemicalShift(const GetPot &controlfile)
    : ImagingDirect < VPixelDimension  > (controlfile){}

private:
};
void WriteImage(InputImageType::Pointer, std::ostringstream &,
                InputImageType::SizeType &, std::string  );
void ImageStats(int itime, const char *ImageName,InputImageType::Pointer Image);
#endif
