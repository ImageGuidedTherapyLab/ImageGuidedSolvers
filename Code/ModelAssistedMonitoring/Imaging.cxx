// tao includes
#include "tao.h"

// libmesh includes
#include "libmesh.h"
#include "dof_map.h"
#include "petsc_macro.h"
#include "mesh.h"
#include "o_string_stream.h"
#include "getpot.h"
#include "parameters.h"

// itk includes
#include "itkRegionOfInterestImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkSubtractImageFilter.h"
#include <itksys/SystemTools.hxx>

// itk includes
#include "itkMedianImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkScalarImageToHistogramGenerator.h"

// local includes
#include "petsc_fem_system.h"
#include "Imaging.h"
#include "itkVTKImageVariableNameIO.h"

const PetscScalar noiseTol = 5.e-2;
void Imaging::printSelf(std::ostream& os)
{
 PetscFunctionBegin; 
  if( !libMesh::processor_id() )
   {
    os << "Imaging::m_Debug            = "<< m_Debug            << std::endl;
    os << "Imaging::sp                 = "<< sp                 << std::endl;
    os << "Imaging::orgn               = "<< orgn               << std::endl;
    os << "Imaging::size               = "<< size               << std::endl;
    os << "Imaging::n_roi              = "<< n_roi              << std::endl;
    os << "Imaging::index              = "<< index              << std::endl;
    os << "Imaging::OutputDir          = "<< OutputDir          << std::endl;
    os << "Imaging::ExamPath           = "<< ExamPath           << std::endl;
    os << "Imaging::DirId              = "<< DirId              << std::endl;
    os << "Imaging::necho              = "<< necho              << std::endl;
    os << "Imaging::median             = "<< median             << std::endl;
    os << "Imaging::noffset            = "<< noffset            << std::endl;
    os << "Imaging::m_ActivationEnergy = "<< m_ActivationEnergy << std::endl;
    os << "Imaging::m_FrequencyFactor  = "<< m_FrequencyFactor  << std::endl;
    os << "Imaging::m_GasConstant      = "<< m_GasConstant      << std::endl;
    os << "Imaging::m_BaseTemperature  = "<< m_BaseTemperature  << std::endl;
    os << "Imaging::m_alpha  (ppm/degC)= "<< m_alpha            << std::endl;
    os << "Imaging::maxdiff    (degC)  = "<< maxdiff            << std::endl;
    os << "Imaging::m_echo time  (ms)  = "<< m_echotime         << std::endl;
    os << "Imaging::m_imaging freq(MHz)= "<< m_imagfreq         << std::endl;
    os << "Imaging::tmap_factor        = "<< this->GetTmapFactor() << std::endl;
    os << "Imaging::magintude file name= "<< this->GetMagnitudeImageFileName() << std::endl;
   }
 PetscFunctionReturnVoid(); 
}
// constructor
#undef __FUNCT__
#define __FUNCT__ "Imaging::Imaging"
Imaging::Imaging(const  GetPot &controlfile )
{
  /* Read options */

  // defaults are NO FILTERING and no debugging
  median  = controlfile("mrti/median" ,   0   );
  maxdiff = controlfile("mrti/maxdiff", 1.e16 );
  m_Debug=PETSC_FALSE;
  if( controlfile("mrti/debug", false ) ) m_Debug=PETSC_TRUE;

  //basic data
  noffset    = controlfile("mrti/noffset",0);
  zeroFilterRadius[0] = 0; // radius along x
  zeroFilterRadius[1] = 0; // radius along y
  zeroFilterRadius[2] = 0; // radius along z
 
  //directory info... 
  ExamPath  = controlfile("mrti/exampath" ,"/file/not/found");
  OutputDir = controlfile("mrti/outputdir","mrivis");
  DirId     = controlfile("mrti/dirid",0);

  //magnitude file for spacing info
  std::ostringstream dflt_file_name;
  dflt_file_name << "/Processed/s" << DirId << "/magnitude.mha" ;
  std::string magFileName( controlfile("setup/imagefile",
                                    dflt_file_name.str() ) );
  std::ostringstream FullMagnitudeFileName;
  FullMagnitudeFileName << ExamPath << "/" << magFileName ;
  MagnFile = FullMagnitudeFileName.str();

  // initial measurement data covariance [deg C]
  // TODO: is this used anywhere important... could prob delete
  m_meascov = 1.0;
  PetscOptionsGetScalar(PETSC_NULL,"-meascov",&m_meascov,PETSC_NULL);

  //  Software Guide : BeginLatex
  //
  //  The size of the neighborhood is defined along every dimension by
  //  passing a \code{SizeType} object with the corresponding values. The
  //  value on each dimension is used as the semi-size of a rectangular
  //  box. For example, in $2D$ a size of \(1,2\) will result in a $3 \times
  //  5$ neighborhood.
  //
  //  \index{itk::MedianImageFilter!Radius}
  //  \index{itk::MedianImageFilter!Neighborhood}
  //
  //  Software Guide : EndLatex 
  medianFilterRadius[0] = median; // radius along x
  medianFilterRadius[1] = median; // radius along y
  medianFilterRadius[2] = 0; // radius along z

  // stats
  maxStat = controlfile("mrti/maxstat",10);

  statCount = 1; // initialize statistics counter

  // Arrhenius parameters (default to biotex values)
  m_ActivationEnergy = controlfile("arrhenius/ActivationEnerg",3.1e98    );// 1/s
  m_FrequencyFactor  = controlfile("arrhenius/FrequencyFactor",6.28e5    );// J/mol
  m_GasConstant      = controlfile("arrhenius/GasConstant"    ,8.314472e0);// J/mol/K
  m_BaseTemperature  = controlfile("arrhenius/BaseTemperature",273.0e0   );// convert to Kelvin

  m_BaseTemperature  = controlfile("arrhenius/BaseTemperature",273.0e0   );// convert to Kelvin
  m_bodyTemp     = controlfile("initial_condition/u_init",37.0);//celcius

  m_nzero = controlfile("mrti/nzero",0);

  m_acquisition_dt = 1.0;

  // faster to multiply
  // ROI noise in air of magnitude image
  noise = 0.0; 

  //itk::Indent indentor =  itk::Indent::New();
  // initialize interpolator
  // default to linear interpolator
  if(controlfile("interpolate/nearest",false)) 
    m_interpolator = NearestInterpolatorType::New();
  else
    m_interpolator = LinearInterpolatorType::New();
  std::cout << m_interpolator;

  // setup import filters
  m_importFilter = ImportFilterType::New();      
}
//get dicom header data
#undef __FUNCT__
#define __FUNCT__ "Imaging::GetHeaderData"
PetscErrorCode Imaging::GetHeaderData(ImageAcquisitionType inputmethod)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  // We construct one instance of the series reader object. Set the DICOM
  // image IO object to be use with it, and set the list of filenames to
  // read.
  InputReaderType::Pointer       realReader =   InputReaderType::New();
  InputImageType::Pointer        realImage ;

  // header info should already be pre processes
  std::ostringstream HeaderIniFileName;
  HeaderIniFileName<< ExamPath << "/Processed/s" << DirId << "/imageHeader.ini"; 
  
  if( !itksys::SystemTools::FileExists(HeaderIniFileName.str().c_str(),true) )
   {std::cerr<< std::endl
             <<"######################### "<< HeaderIniFileName.str() 
             << " not found "<< std::endl << std::flush; abort();
   }
  GetPot HeaderIni( HeaderIniFileName.str() );
  std::cout << "######################### opened "<< HeaderIniFileName.str() 
            << std::endl << std::flush ; 

  if( !HeaderIni.have_variable("rawdata/necho") ) 
   {std::cerr<< std::endl
             <<"#########################Error getting # echo " 
             << std::endl << std::flush; abort();
   }
  necho = HeaderIni("rawdata/necho", 1 );

  if( !HeaderIni.have_variable("rawdata/echotime") ) 
   {std::cerr<< std::endl
             <<"#########################Error getting echo time " 
             << std::endl << std::flush; abort();
   }
  m_echotime = HeaderIni("rawdata/echotime", 0.0 );

  if( !HeaderIni.have_variable("rawdata/imagfreq") ) 
     {std::cerr<< std::endl
               << "#################Error getting Imaging Freq"
               << std::endl << std::flush; abort();
     }
  m_imagfreq = HeaderIni("rawdata/imagfreq", 63.64 );

  // update from command line if needed
  m_alpha = -0.0097; 
  ierr=PetscOptionsGetScalar(PETSC_NULL,"-alpha",&m_alpha,PETSC_NULL);

  // default is not to scale image
  PetscScalar scalefactor = 1.0;
  ierr=PetscOptionsGetScalar(PETSC_NULL,"-scaling",&scalefactor,PETSC_NULL);

  // read in dimension information from header 
  switch(inputmethod)
   {// 2d dicom series
    case DICOM2D: case BACKGROUNDPHASE: 
    case EIGHT_ECHO_CSI: case SIXTEEN_ECHO_CSI:
     {

      // setup space for file names
      // generate first set of images filenames
      for(int i = 0 ; i < 2 * necho ; i++)
       {
         filenames.push_back( std::vector< std::string >::vector(1,"") );
       }

      GetImageFileName(1,std::string("imageReal."),filenames[0][0]);

      // write out header info
      this->WriteIni();
     } break;
    case VOLUME3D:  // 3d Volume image
     {

      // setup space for file names   
      filenames.push_back( std::vector< std::string >::vector(1,"") );
      // generate magnitude image filename
      filenames[0][0] = this->GetMagnitudeImageFileName() ;

     } break;
    default: 
      std::cout << "unknown inputmethod " << inputmethod << std::endl;
      PetscFunctionReturn(0);
   }

  // only need to read in the first image to get the header info
  realReader->SetFileName( filenames[0][0] );
  // need to explicitly instantiate the template
  realImage = GetImageData< PetscScalar , InputReaderType::Pointer > ( 
                            realReader, libMesh::processor_id(), filenames[0]);

  // get size information
  size=realImage->GetRequestedRegion().GetSize();
  sp = scalefactor * realImage->GetSpacing();
  const InputImageType::PointType& orgn_mm = realImage->GetOrigin();
  orgn[0] = scalefactor * orgn_mm[0] ;
  orgn[1] = scalefactor * orgn_mm[1] ;
  orgn[2] = scalefactor * orgn_mm[2] ;

  // resize filenames
  filenames.clear();
  for(int i = 0 ; i < this->necho; i++ )
   {
    filenames.push_back( std::vector< std::string >::vector(size[2],"") );
    filenames.push_back( std::vector< std::string >::vector(size[2],"") );
   }

  // allocate full images buffers
  InputImageType::RegionType    fullRegion;
  fullRegion.SetSize(  this->size  );

  meanImage = InputImageType::New();
  meanImage->SetRegions( fullRegion );
  meanImage->Allocate();
  meanImage->SetSpacing( sp  );
  meanImage->SetOrigin(  orgn);
  meanImage->FillBuffer( this->m_bodyTemp ); // initial mean should be body temp

  snr_Image = InputImageType::New();
  snr_Image->SetRegions( fullRegion );
  snr_Image->Allocate();
  snr_Image->SetSpacing( sp  );
  snr_Image->SetOrigin(  orgn);
  snr_Image->FillBuffer( m_meascov ); // initial snr based temp error is meascov

  var_Image = InputImageType::New();
  var_Image->SetRegions( fullRegion );
  var_Image->Allocate();
  var_Image->SetSpacing( sp  );
  var_Image->SetOrigin(  orgn);
  var_Image->FillBuffer( m_meascov ); // initialize to meascov 

  // allocate full images 
  net_Image = InputImageType::New();
  net_Image->SetRegions( fullRegion );
  net_Image->Allocate();
  net_Image->SetSpacing( sp  );
  net_Image->SetOrigin(  orgn);
  net_Image->FillBuffer( this->m_bodyTemp );

  maskImage = InputImageType::New();
  maskImage->SetRegions( fullRegion );
  maskImage->Allocate();
  maskImage->SetSpacing( sp  );
  maskImage->SetOrigin(  orgn);

  arr_Image = InputImageType::New();
  arr_Image->SetRegions( fullRegion );
  arr_Image->Allocate();
  arr_Image->SetSpacing( sp  );
  arr_Image->SetOrigin(  orgn);
  arr_Image->FillBuffer( 0.0 );

  // output initial maps
  if( !libMesh::processor_id() ) 
   {
    std::cout << "initial measurement Covariance "<< m_meascov << std::endl;
    ImageStats( this->m_nzero ,"meanImage",meanImage);
    ImageStats( this->m_nzero ,"var_Image",var_Image);
    ImageStats( this->m_nzero ,"snr_Image",snr_Image);
    if( m_Debug )
     for( int iii = 0 ; iii <= this->m_nzero; iii++)
      {
       std::ostringstream unfiltered_filename;
       unfiltered_filename << OutputDir <<"/average.";
       OSSRealzeroright(unfiltered_filename,4,0,iii);
       unfiltered_filename << ".vtk" ;
       WriteImage(meanImage, unfiltered_filename , medianFilterRadius ,
                                             std::string("average")  ) ; 

       unfiltered_filename.str(""); // reset before reuse
       unfiltered_filename << OutputDir <<"/variance.";
       OSSRealzeroright(unfiltered_filename,4,0,iii);
       unfiltered_filename << ".vtk" ;
       WriteImage(var_Image, unfiltered_filename , medianFilterRadius,
                                             std::string("variance")  ) ; 

       // output SNR image
       std::ostringstream snrfilename;
       snrfilename << OutputDir <<"/snr.";
       OSSRealzeroright(snrfilename,4,0,iii);
       snrfilename << ".vtk";
       WriteImage(snr_Image , snrfilename , zeroFilterRadius,
                                             std::string("snr")  ) ; 
      }
   }
  // acquisition time
                            
  if(!libMesh::processor_id()) std::cout << "Spacing = "
                      << sp[0] << ", " << sp[1] << ", " << sp[2] << std::endl;
  if(!libMesh::processor_id()) std::cout << "Origin = " << orgn[0] << ", " 
                                     << orgn[1] << ", " << orgn[2] << std::endl;
  if(!libMesh::processor_id()) std::cout << 
                         "   ****Dimension Data****   " <<std::endl<<std::endl;

  PetscFunctionReturn(0);
}
// subroutine to write the image to disk
#undef __FUNCT__
#define __FUNCT__ "WriteImage"
void WriteImage(InputImageType::Pointer Image,std::ostringstream &filename,
                InputImageType::SizeType &filterRadius, std::string variableName)
{
  PetscFunctionBegin;

  // setup writer
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( filename.str() );
  itk::VTKImageVariableNameIO::Pointer ioBasePointer = itk::VTKImageVariableNameIO::New();
  // set variable name
  std::ostringstream annotatedVariableName;
  annotatedVariableName<< variableName ;
                       //<< FiniteElementInterface::profileID ;
  ioBasePointer->SetVariableName( annotatedVariableName.str()  ) ;

  writer->SetImageIO( ioBasePointer ); // use custom class for proper data name
  std::cout << "writing " << filename.str() << std::endl;

  //  Software Guide : BeginLatex
  //
  //  Using the image types, it is now possible to define the filter type
  //  and create the filter object.
  //
  //  \index{itk::MedianImageFilter!instantiation}
  //  \index{itk::MedianImageFilter!New()}
  //  \index{itk::MedianImageFilter!Pointer}
  // 
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  typedef itk::MedianImageFilter<
               InputImageType, OutputImageType >  OutputFilterType;

  OutputFilterType::Pointer filter = OutputFilterType::New();
  // Software Guide : EndCodeSnippet

  // Software Guide : BeginCodeSnippet
  
  filter->SetRadius( filterRadius );
  // Software Guide : EndCodeSnippet

  //  Software Guide : BeginLatex
  //
  //  The input to the filter can be taken from any other filter, for example
  //  a reader. The output can be passed down the pipeline to other filters,
  //  for example, a writer. An update call on any downstream filter will
  //  trigger the execution of the median filter.
  //
  //  \index{itk::MedianImageFilter!SetInput()}
  //  \index{itk::MedianImageFilter!GetOutput()}
  //
  //  Software Guide : EndLatex 


  // Software Guide : BeginCodeSnippet
  filter->SetInput( Image );
  writer->SetInput( filter->GetOutput() );
  // Software Guide : EndCodeSnippet
  try
    {
    writer->Update();
    }
  catch (itk::ExceptionObject &ex)
    {
    std::cout << ex; abort();
    }
  PetscFunctionReturnVoid();
}

// subroutine to write the image to disk
#undef __FUNCT__
#define __FUNCT__ "ImageStats"
void ImageStats(int itime, const char *ImageName,InputImageType::Pointer Image)
{
// Software Guide : BeginLatex
//
// We use now the image type in order to instantiate the type of the
// corresponding histogram generator class, and invoke its \code{New()} method
// in order to construct one.
//
// \index{itk::Statistics::Scalar\-Image\-To\-Histogram\-Generator!header}
//
// Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
  typedef itk::Statistics::ScalarImageToHistogramGenerator< 
                                 InputImageType >   HistogramGeneratorType;

  HistogramGeneratorType::Pointer histogramGenerator =
                                        HistogramGeneratorType::New();
// Software Guide : EndCodeSnippet
// Software Guide : BeginLatex
//
// The image to be passed as input to the histogram generator is taken in this
// case from the output of an image reader.
//
// Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
  histogramGenerator->SetInput( Image );
// Software Guide : EndCodeSnippet

// Software Guide : BeginLatex
//
// We define also the typical parameters that specify the characteristics of
// the histogram to be computed.
//
// Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
  histogramGenerator->SetNumberOfBins( 8 );
  histogramGenerator->SetMarginalScale( 10.0 );

  histogramGenerator->SetHistogramMin(-255.5 );
  histogramGenerator->SetHistogramMax( 255.5 );
// Software Guide : EndCodeSnippet


// Software Guide : BeginLatex
//
// Finally we trigger the computation of the histogram by invoking the
// \code{Compute()} method of the generator. Note again, that a generator is
// not a pipeline object and therefore it is up to you to make sure that the
// filters providing the input image have been updated.
//
// \index{itk::Statistics::Scalar\-Image\-To\-Histogram\-Generator!Compute()}
//
// Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
  histogramGenerator->Compute();
// Software Guide : EndCodeSnippet

// Software Guide : BeginLatex
//
// The resulting histogram can be obtained from the generator by invoking its
// \code{GetOutput()} method. It is also convenient to get the Histogram type
// from the traits of the generator type itself as shown in the code below.
//
// Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
  typedef HistogramGeneratorType::HistogramType  HistogramType;

  const HistogramType * histogram = histogramGenerator->GetOutput();
// Software Guide : EndCodeSnippet

  const unsigned int histogramSize = histogram->Size();


// Software Guide : BeginLatex
//
// In this case we simply print out the frequency values of the histogram.
// These values can be accessed by using iterators.
//
// \index{itk::Statistics::Histogram!Iterators}
//
// Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
//  HistogramType::ConstIterator itr = histogram->Begin();
//  HistogramType::ConstIterator end = histogram->End();
// 
//  unsigned int binNumber = 0;
//  while( itr != end )
//    {
//    std::cout << "bin = " << binNumber << " frequency = ";
//    std::cout << itr.GetFrequency() << std::endl;     
//    ++itr;
//    ++binNumber;
//    }
// Software Guide : EndCodeSnippet

  typedef itk::StatisticsImageFilter<InputImageType> FilterType;
  FilterType::Pointer filter = FilterType::New();

  filter->SetInput (Image);
  filter->UpdateLargestPossibleRegion();

  std::cout << ImageName << " time = " << itime 
            << " Minimum "  << filter->GetMinimum()              
            << " Maximum "  << filter->GetMaximum()              
            << " Sum "      << filter->GetSum()                  
            << " Mean "     << filter->GetMean()                 
            << " Variance " << filter->GetVariance()<< std::endl;

  std::cout << " Histogram size " << histogramSize << std::endl;

  //unsigned int bin;
  //for( bin=0; bin < histogramSize; bin++ )
  //  {
  //  std::cout << "bin = " << bin << " frequency = ";
  //  std::cout << histogram->GetFrequency( bin, 0 ) << std::endl;
  //  }
}

// setup header info
#undef __FUNCT__
#define __FUNCT__ "Imaging::SetHeaderData"
PetscErrorCode Imaging::SetHeaderData( 
                                       int nx_image, int ny_image, int nz_image,
				       double X0, double Y0, double Z0, 
                                       double DX, double DY, double DZ )
{
  PetscFunctionBegin; 

  // pixel dimensions
  size[0]  = nx_image;  // size along X
  size[1]  = ny_image;  // size along Y
  size[2]  = nz_image;  // size along Z

  ImportFilterType::IndexType start;
  start.Fill( 0 );

  ImportFilterType::RegionType region;
  region.SetIndex( start );
  region.SetSize(  size  );

  m_importFilter->SetRegion( region );

  // origin 
  orgn[0] = X0;    // X coordinate 
  orgn[1] = Y0;    // Y coordinate
  orgn[2] = Z0;    // Z coordinate

  m_importFilter->SetOrigin( orgn );

  // spacing
  sp[0] = DX;    // along X direction 
  sp[1] = DY;    // along Y direction
  sp[2] = DZ;    // along Z direction

  m_importFilter->SetSpacing( sp );

  PetscFunctionReturn(0);
}
// return an image mask from the real/imaginary images 
#undef __FUNCT__
#define __FUNCT__ "Imaging::GetImageMask"
void Imaging::GetImageMask()
{
  PetscFunctionBegin;

//  // containers for real and imaginary data
//  VecFilterType::Pointer   realImageFilter=VecFilterType::New(),// real images
//                           imagImageFilter=VecFilterType::New();// imag images
//
//  // pointers to real and imaginary vector images for CSI computation
//  VecImageType::Pointer realImages, imagImages;
//
//  // get image data
//  realImages = GetVecTimeImage(realImageFilter,2,0); //real images
//  imagImages = GetVecTimeImage(imagImageFilter,2,1); //imag images
//
//  // Software Guide : BeginLatex
//  //
//  // \index{Iterators!construction of} \index{Iterators!and image regions} The
//  // necessary images and region definitions are now in place.  All that is
//  // left to do is to create the iterators and perform the copy.  Note that
//  // image iterators are not accessed via smart pointers so they are
//  // light-weight objects that are instantiated on the stack.  Also notice how
//  // the input and output iterators are defined over the \emph{same
//  // corresponding region}.  Though the images are different sizes, they both
//  // contain the same target subregion.
//  //
//  // Software Guide : EndLatex
//  
//  // Software Guide : BeginCodeSnippet
//
//  // initialize ITK iterators
//  VecIteratorType    realIt( realImages, realImages->GetRequestedRegion() );
//  VecIteratorType    imagIt( imagImages, imagImages->GetRequestedRegion() );
//  BufferIteratorType maskIt( maskImage ,  maskImage->GetRequestedRegion() );
//  realIt.SetFirstDirection(  0 );   realIt.SetSecondDirection( 1 );
//  imagIt.SetFirstDirection(  0 );   imagIt.SetSecondDirection( 1 );
//  maskIt.SetFirstDirection(  0 );   maskIt.SetSecondDirection( 1 );
//
//  // data structure for max
//  std::vector< double > maxvalue(nslice,0.0);
//
//  double avgnorm; // scratch
//
//  // set iterators to the beginning
//  maskIt.GoToBegin();
//  realIt.GoToBegin();
//  imagIt.GoToBegin();
//  for(int kkk=0; kkk < nslice ; kkk++)
//    {
//    while ( !maskIt.IsAtEndOfSlice() )
//      {
//      while ( !maskIt.IsAtEndOfLine() )
//        {
//        avgnorm= std::sqrt(
//                 pow((0.5*(realIt.Get()[0] + realIt.Get()[1])),2) 
//                                 +
//                 pow((0.5*(imagIt.Get()[0] + imagIt.Get()[1])),2)
//                                 ) ;
//        if(maxvalue[kkk] < avgnorm) maxvalue[kkk] = avgnorm; 
//        ++maskIt;
//        ++realIt;
//        ++imagIt;
//        }
//        maskIt.NextLine();
//        realIt.NextLine();
//        imagIt.NextLine();
//      }
//    // get next slice
//    maskIt.NextSlice();
//    realIt.NextSlice();
//    imagIt.NextSlice();
//    }
//
//  // set iterators to the beginning
//  maskIt.GoToBegin();
//  realIt.GoToBegin();
//  imagIt.GoToBegin();
//  for(int kkk=0; kkk < nslice ; kkk++)
//    {
//    double snrthreshold = 0.05 * maxvalue[kkk];
//    snrthreshold = -1.0; // turn off
//    while ( !maskIt.IsAtEndOfSlice() )
//      {
//      while ( !maskIt.IsAtEndOfLine() )
//        {
//        avgnorm= std::sqrt(
//                 pow((0.5*(realIt.Get()[0] + realIt.Get()[1])),2) 
//                                 +
//                 pow((0.5*(imagIt.Get()[0] + imagIt.Get()[1])),2)
//                                 ) ;
//        if(avgnorm > snrthreshold) 
//          maskIt.Set(1.0);
//        else
//          maskIt.Set(0.0);
//        ++maskIt;
//        ++realIt;
//        ++imagIt;
//        }
//        maskIt.NextLine();
//        realIt.NextLine();
//        imagIt.NextLine();
//      }
//    // get next slice
//    maskIt.NextSlice();
//    realIt.NextSlice();
//    imagIt.NextSlice();
//    }
// 
//  // output image mask
//  std::ostringstream maskfilename;
//  maskfilename << OutputDir <<"/mask."<< rank << ".vtk";
//  WriteImage(maskImage , maskfilename , zeroFilterRadius);

  PetscFunctionReturnVoid();
}
// update SNR from magnitude image and write
#undef __FUNCT__
#define __FUNCT__ "ImagingComplex::UpdateSNR"
void ImagingComplex::UpdateSNR(InputImageType::Pointer &magnitudeImage, const int )
{
  PetscFunctionBegin;

  // Software Guide : EndCodeSnippet
  
  // Software Guide : BeginLatex
  //
  // The following code creates an output image and iterator.
  //   
  // Software Guide : EndLatex
    
  // Software Guide : BeginCodeSnippet
  magnitudeImage->SetRequestedRegionToLargestPossibleRegion();
  RegionIteratorType snrIt(     snr_Image,     snr_Image->GetRequestedRegion()),
                     magIt(magnitudeImage,magnitudeImage->GetRequestedRegion());
  // Software Guide : EndCodeSnippet
  
  
  // Software Guide : BeginLatex
  //
  // It is equivalent to use the six corresponding integer array indices instead.
  // For example, the offsets \code{(-1,-1)} and \code{(1, -1)} are 
  // equivalent to the integer indices \code{0} and \code{2}, respectively.
  //
  // The calculations are done in a \code{for} loop that moves the input and
  // output iterators synchronously across their respective images.  The
  // \code{sum} variable is used to sum the results of the finite differences.
  //
  // Software Guide : EndLatex
  
  // Software Guide : BeginCodeSnippet
  snrIt.GoToBegin(); 
  magIt.GoToBegin();
  while ( !snrIt.IsAtEnd())
    {
     if( magIt.Get() < noiseTol)
        snrIt.Set(noiseTol);
     else if( noise > 0.0 )
        snrIt.Set(magIt.Get()/noise);
     else 
        {std::cout << "noise <= 0.0 " << std::endl << std::flush; abort();}
     ++snrIt; 
     ++magIt;
    }
  // Software Guide : EndCodeSnippet
    
  PetscFunctionReturnVoid();
}
// update noise from the an ROI in air in the magnitude image
#undef __FUNCT__
#define __FUNCT__ "ImagingComplex::ExtractNoise"
void ImagingComplex::ExtractNoise(const int istep)
{
  PetscFunctionBegin;

  // FiniteElementInterface Space to hold two magnitude images
  InputReaderType::Pointer  realreaderOne = InputReaderType::New(),
                            imagreaderOne = InputReaderType::New(),
                            realreaderTwo = InputReaderType::New(),
                            imagreaderTwo = InputReaderType::New();
  MagnitudeFilterType::Pointer
                 magnitudeFilterOne=MagnitudeFilterType::New(),
                 magnitudeFilterTwo=MagnitudeFilterType::New();

  // get first image set 
  RealTimeGenerateFileNames(istep,0,filenames);
  realreaderOne->SetFileNames(filenames[0] );
  InputImageType::Pointer realImageOne = 
         GetImageData< PetscScalar , InputReaderType::Pointer > (realreaderOne,
                                                     libMesh::processor_id() , filenames[0]);
  imagreaderOne->SetFileNames(filenames[1] );
  InputImageType::Pointer imagImageOne = 
         GetImageData< PetscScalar , InputReaderType::Pointer > (imagreaderOne,
                                                    libMesh::processor_id() , filenames[1]) ;
  magnitudeFilterOne->SetInput1(   realImageOne );
  magnitudeFilterOne->SetInput2(   imagImageOne );

  // get second image set 
  RealTimeGenerateFileNames(istep+1,0,filenames);
  realreaderTwo->SetFileNames( filenames[0] );
  InputImageType::Pointer realImageTwo = 
         GetImageData< PetscScalar , InputReaderType::Pointer > (realreaderTwo,
                                                     libMesh::processor_id() , filenames[0]);
  imagreaderTwo->SetFileNames( filenames[1] );
  InputImageType::Pointer imagImageTwo = 
         GetImageData< PetscScalar , InputReaderType::Pointer > (imagreaderTwo,
                                                    libMesh::processor_id() , filenames[1]) ;
  magnitudeFilterTwo->SetInput1(   realImageTwo );
  magnitudeFilterTwo->SetInput2(   imagImageTwo );

  typedef itk::SubtractImageFilter< InputImageType , InputImageType ,
                                    InputImageType > SubtractFilterType;
  SubtractFilterType::Pointer
                 subtractFilter=SubtractFilterType::New();
  subtractFilter->SetInput1(   magnitudeFilterOne->GetOutput() );
  subtractFilter->SetInput2(   magnitudeFilterTwo->GetOutput() );
  try 
    { 
    subtractFilter->Update(); 
    } 
  catch( itk::ExceptionObject & err ) 
    { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl << std::flush; abort();
    } 

  // Software Guide : EndCodeSnippet

  // get ROI
  PetscErrorCode ierr;
  PetscInt magNx = 10            , magIx =230,
           magNy = 10            , magIy =230,
           magNz = size[2]       , magIz =0;
  ierr=PetscOptionsGetInt(PETSC_NULL,"-magNx",&magNx,PETSC_NULL);CHKERRV(ierr);
  ierr=PetscOptionsGetInt(PETSC_NULL,"-magNy",&magNy,PETSC_NULL);CHKERRV(ierr);
  ierr=PetscOptionsGetInt(PETSC_NULL,"-magNz",&magNz,PETSC_NULL);CHKERRV(ierr);
  ierr=PetscOptionsGetInt(PETSC_NULL,"-magIx",&magIx,PETSC_NULL);CHKERRV(ierr);
  ierr=PetscOptionsGetInt(PETSC_NULL,"-magIy",&magIy,PETSC_NULL);CHKERRV(ierr);
  ierr=PetscOptionsGetInt(PETSC_NULL,"-magIz",&magIz,PETSC_NULL);CHKERRV(ierr);

  // FIXME: create a copy to avoid set spacing and avoid errors on update
  InputImageType::Pointer subtractImage = InputImageType::New();
  subtractImage->SetRegions(subtractFilter->GetOutput()->GetRequestedRegion());
  subtractImage->Allocate();
  subtractImage->SetSpacing( sp  );
  subtractImage->SetOrigin(  orgn);
  subtractImage->FillBuffer(    0.0     );

  // FiniteElementInterface iterator to only iterate over roi
  InputImageType::RegionType roiRegion; 
  InputImageType::RegionType::IndexType roiStart;
  InputImageType::RegionType::SizeType  roiWidth;
  roiStart[0]=magIx; roiWidth[0]=magNx;
  roiStart[1]=magIy; roiWidth[1]=magNy;
  roiStart[2]=magIz; roiWidth[2]=magNz;
  roiRegion.SetSize(  roiWidth );
  roiRegion.SetIndex( roiStart );

  RegionIteratorType roiIt(subtractFilter->GetOutput() , roiRegion);
  RegionIteratorType tmpIt(subtractImage               , roiRegion);
  
  // loop over pixels and compute the average
  // store ROI for plotting verification of ROI position
  int count=0; double avgIntensity = 0.0;
  for(roiIt.GoToBegin(),tmpIt.GoToBegin(); !roiIt.IsAtEnd(); ++roiIt,++tmpIt)
    {/* {kkk,jjj,iii} = roiIt.GetIndex(); should return such that
            magIz <= kkk < magIz + magNz 
            magIy <= jjj < magIy + magNy 
            magIx <= iii < magIx + magNx */ 
     avgIntensity += roiIt.Get(); count++; 
     tmpIt.Set( roiIt.Get() );
    }
  avgIntensity = avgIntensity/count;

  // compute stddev 
  double stddevIntensity = 0.0;
  for(roiIt.GoToBegin(); !roiIt.IsAtEnd(); ++roiIt)
      stddevIntensity += std::pow( roiIt.Get() - avgIntensity , 2 ); 
  stddevIntensity = std::sqrt(stddevIntensity/count);
  /* Difference Method of calculating SNR 
      requires an additional factor of sqrt(2) 
     [Reeder(CH4)-Measurment of Signal to Noise Ratio pp55]
     FIXME - Derive the sqrt(2) factor
     [Price etal 1990] */
  stddevIntensity = stddevIntensity / std::sqrt(2.0);

  // Convert from Rayleigh Distribution to Gaussian (appendix B in Haacke)
  noise = stddevIntensity/0.655; 

  // error check
  if(noise < 0.0)
   {
     std::cout << "negative Noise?" << std::endl << std::flush ;  abort();
   }
  else if(noise < noiseTol)
   {
     std::cout << "zero Noise?" << std::endl << std::flush ;  noise=noiseTol;
   }

  // output ROI image and noise value
  if(!libMesh::processor_id() ) 
   {
     std::ostringstream roifilename;
     roifilename << OutputDir <<"/magnitudeROI.mha";
     //OSSRealzeroright(roifilename,4,0,istep);
     //roifilename << ".vtk";
     WriteImage( subtractImage , roifilename , zeroFilterRadius ,
                                          std::string("magnitude")  ) ; 
     std::cout << "Magnitude Noise[" << istep << "]="<< noise 
               << " noise conversion from Mean[" << istep << "]="
               << avgIntensity / 1.25 << std::endl << std::flush ; 
   }
  // Software Guide : EndCodeSnippet

  PetscFunctionReturnVoid();
}

// turn on debuggin
#undef __FUNCT__
#define __FUNCT__ "FiniteElementInterface::DebugOn"
void Imaging::DebugOn()
{
  PetscFunctionBegin;
  net_Image->DebugOn();
  //sp[0]= 1.00;
  //sp[1]= 1.00;
  //sp[2]= 1.00;
  //orgn[0] = 0.0;
  //orgn[1] = 0.0;
  //orgn[2] = 0.0;
  PetscFunctionReturnVoid();
}

// constructor
#undef __FUNCT__
#define __FUNCT__ "ImagingComplex::GetHeaderData"
PetscErrorCode ImagingComplex::GetHeaderData(ImageAcquisitionType inputmethod)
{
  //PetscErrorCode ierr;
  PetscFunctionBegin;
  // setup base class info...
  this->Parent::GetHeaderData(inputmethod); 

  // update noise from the an ROI in air in the magnitude image
  this->ExtractNoise( this->m_nzero );

  ClassFilterType::Pointer   buffphaseFilter=ClassFilterType::New();
  ClassImageType::Pointer  buffImage;  
  // read in everything before nzero to write out magnitude and rawdata files
  // for consistency. Read nzero last to have the image at nzero for the
  // reference
  for( int iii = 0 ; iii <= this->m_nzero  ; iii++)
       buffImage = GetPhaseImage( iii ,buffphaseFilter);
  // FiniteElementInterface iterators
  itk::ImageRegionConstIterator< ClassImageType > 
                      buffIt(buffImage, buffImage->GetRequestedRegion());
  itk::ImageRegionIterator< ClassImageType > 
                      baseIt(baseImage, baseImage->GetRequestedRegion());
  // set iterators to the beginning
  buffIt.GoToBegin();
  baseIt.GoToBegin();
  //compute net phase difference
  while( !buffIt.IsAtEnd() )
    {
     baseIt.Set( buffIt.Get() );
     ++buffIt;
     ++baseIt;
    }

  if(!libMesh::processor_id() ) 
   {
    ImageStats( this->m_nzero,"net_Image",net_Image);
    ImageStats( this->m_nzero,"arr_Image",arr_Image);
    for( int iii = 0 ; iii <= this->m_nzero  ; iii++)
     {
       // output unfiltered tmap
       std::ostringstream unfiltered_filename;
       unfiltered_filename << OutputDir <<"/temperature.";
       OSSRealzeroright(unfiltered_filename,4,0,iii);
       unfiltered_filename << ".vtk" ;
       WriteImage(net_Image, unfiltered_filename , medianFilterRadius,
                                          std::string("temperature")  ) ; 
       // output damage map
       std::ostringstream damage_filename;
       damage_filename << OutputDir <<"/damage.";
       OSSRealzeroright(damage_filename,4,0,iii);
       damage_filename << ".vtk" ;
       WriteImage(arr_Image, damage_filename , medianFilterRadius,
                                          std::string("damage")  ) ; 
     }
   }
  PetscFunctionReturn(0);
}
// return a phase image from the real/imaginary images
#undef __FUNCT__
#define __FUNCT__ "ImagingComplex::GetPhaseImage"
ImagingComplex::ClassImageType::Pointer ImagingComplex::GetPhaseImage(int istep,
                                       ClassFilterType::Pointer complexFilter)
{
  PetscFunctionBegin;

  // create image mask 
  // GetImageMask();

  RealTimeGenerateFileNames(istep,0,filenames);

  // We construct one instance of the series reader object. Set the DICOM
  // image IO object to be use with it, and set the list of filenames to read.
  InputReaderType::Pointer  realreader = InputReaderType::New(),
                            imagreader = InputReaderType::New();

  // get real images
  realreader->SetFileNames(filenames[0] );
  InputImageType::Pointer realImage = 
       GetImageData< PetscScalar , InputReaderType::Pointer > (realreader, 
                                                     libMesh::processor_id() , filenames[0]);
  // need to set spacing in source images for subsequent filters/writers
  realImage->SetSpacing( sp  );
  realImage->SetOrigin(  orgn);
  complexFilter->PushBackInput( realImage );
  // get imaginary images
  imagreader->SetFileNames(filenames[1] );
  InputImageType::Pointer imagImage = 
       GetImageData< PetscScalar , InputReaderType::Pointer > (imagreader, 
                                                    libMesh::processor_id() , filenames[1]) ;
  // need to set spacing in source images for subsequent filters/writers
  imagImage->SetSpacing( sp  );
  imagImage->SetOrigin(  orgn);
  complexFilter->PushBackInput( imagImage );

  // FiniteElementInterface magnitude filter to update snr
  MagnitudeFilterType::Pointer
                 magnitudeFilter=MagnitudeFilterType::New();
  magnitudeFilter->SetInput1(   realImage );
  magnitudeFilter->SetInput2(   imagImage );
  magnitudeFilter->GetOutput()->SetSpacing( sp  );
  magnitudeFilter->GetOutput()->SetOrigin(  orgn);

  // compute phase
  try
    {
    complexFilter->Update();
    magnitudeFilter->Update();
    }
  catch (itk::ExceptionObject &excp)
    {
    std::cout << "Exception thrown" << excp << std::endl; abort();
    }

  // update SNR measurement and write
  InputImageType::Pointer magnitudeImage = magnitudeFilter->GetOutput();
  this->UpdateSNR( magnitudeImage , istep);

  // write raw data, phase, and magnitude image as well...
  if( !libMesh::processor_id()  && m_Debug )
   {
     // write out a magnitude image for ITK snap
     std::ostringstream magnitude_filename;
     magnitude_filename << OutputDir <<"/magnitude.";
     OSSRealzeroright(magnitude_filename,4,0,istep);
     magnitude_filename<< ".vtk" ;
     WriteImage(magnitudeFilter->GetOutput(), 
                magnitude_filename , zeroFilterRadius,
                std::string("magnitude")  ) ; 

     typedef itk::Image< itk::Vector< PetscScalar, 3 >,Dimension >  
                                                                 RawImageType;
     typedef itk::ScalarToArrayCastImageFilter< InputImageType, RawImageType >
                                                                RawFilterType;
     //  A filter object is created with the New()
     //  method and assigned to a SmartPointer.
     RawFilterType::Pointer  rawfilter=RawFilterType::New();
     // write raw data
     rawfilter->PushBackInput( realImage );
     rawfilter->PushBackInput( imagImage );

     typedef itk::Atan2ImageFilter<InputImageType,InputImageType,InputImageType>
                                                                Atan2FilterType;
     Atan2FilterType::Pointer  atanFilter=Atan2FilterType::New();
     atanFilter->SetInput1(   imagImage );
     atanFilter->SetInput2(   realImage );

     rawfilter->PushBackInput( atanFilter->GetOutput() );

     // FiniteElementInterface writer
     typedef itk::ImageFileWriter<  RawImageType >           RawWriterType;
     RawWriterType::Pointer writer = RawWriterType::New();

     // output base phase image as a place holder
     std::ostringstream rawfilename;
     // output file
     rawfilename << OutputDir <<"/rawdata.";
     OSSRealzeroright(rawfilename,4,0,istep);
     rawfilename << ".vtk" ;

     writer->SetFileName( rawfilename.str() );
     std::cout << "writing " << rawfilename.str() << std::endl;
     // Software Guide : EndCodeSnippet
     try
       {
       rawfilter->Update();
       rawfilter->GetOutput()->SetSpacing( sp  );
       rawfilter->GetOutput()->SetOrigin(  orgn);
       writer->SetInput( rawfilter->GetOutput() );
       writer->Update();
       }
     catch (itk::ExceptionObject &ex)
       {
       std::cout << ex; abort();
       }

   }

  // Software Guide : EndCodeSnippet
  // FIXME not sure why. But need GetOutput(0) instead of GetOutput() in
  // debugger. should be the same though... 
  PetscFunctionReturn(complexFilter->GetOutput(0));
}
// update the tmap and snr maps at this time instance
#undef __FUNCT__
#define __FUNCT__ "ImagingComplex::UpdateImaging"
PetscErrorCode ImagingComplex::UpdateImaging(const int istep) 
{
  PetscFunctionBegin;

  // containers for phase data
  ClassFilterType::Pointer   currphaseFilter=ClassFilterType::New();
  ClassImageType::Pointer  currImage;  

  // get image data
  currImage = GetPhaseImage( istep   ,currphaseFilter);

  // Software Guide : BeginLatex
  // get header info
  // 
  // We can trigger the reading process by calling the \code{Update()}
  // method on the series reader. It is wise to put this invocation inside
  // a \code{try/catch} block since the process may eventually throw
  // exceptions.
  //
  // Software Guide : EndLatex
  
  //foundFlag == PETSC_TRUE ==> flag found on cmd line

  // FiniteElementInterface iterators
  itk::ImageRegionConstIterator< ClassImageType > 
                      currIt(currImage, currImage->GetRequestedRegion());
  itk::ImageRegionIterator< ClassImageType > 
                      baseIt(baseImage, baseImage->GetRequestedRegion());
  RegionIteratorType  snr_It(snr_Image, snr_Image->GetRequestedRegion()),
                      net_It(net_Image, net_Image->GetRequestedRegion()),
                      arr_It(arr_Image, arr_Image->GetRequestedRegion());

  //std::cout << "basePointer" << baseImage->GetBufferPointer() << " " 
  //          << "currPointer" << currImage->GetBufferPointer() 
  //          << std::endl;

  // Software Guide : EndCodeSnippet
  // set iterators to the beginning
  currIt.GoToBegin();
  baseIt.GoToBegin();
  net_It.GoToBegin();
  snr_It.GoToBegin();
  arr_It.GoToBegin();

/* 
  compute net phase difference. It is assumed that the water proton chemical
shift to lower frequenecies with higher temperatures (caused by rupture,
stretching, bending of hydrogen bonds).
@article{ishihara1995paf, 
   title={{A precise and fast temperature mapping using water proton chemical shift}}, 
   author={Ishihara, Y. and Calderon, A. and Watanabe, H. and Okamoto, K. and Suzuki, Y. and Kuroda, K. and Suzuki, Y.},
   journal={Magnetic Resonance in Medicine}, 
   volume={34}, 
   number={6}, 
   year={1995},
   publisher={Wiley Subscription Services, Inc., A Wiley Company Baltimore} 
}
  thus a phase change that corresponds to a decrease in the larmor frequency
implies a positive temperature change. 
                     y
                      |
                      |    _
                      |   |\  positive theta
                      |     |
              -----------------  x
                      |
                      |
                      |
       FIXME: double check w/ rjs
notice that in the usual convention of positive theta with increasing
angle from x-axis ACTUALLY corresponds to a NEGATIVE phase change.
This is because the ROTATING FRAME OF REFERENCE is in the NEGATIVE Z-DIRECTION.
A positive phase shift in the rotating frame corresponds to an increase
in frequency and thus negative theta in the usual convention.

To avoid temporal phase wrap artifacts, MUST accumulate phase as the
complex phase difference of the incremental sum of the Signal at time i, 
$S^i=\|S^i\|e^{i\theta^i}$,
to the signal at time i+1, $S^{i+1}=\|S^{i+1}\|e^{i\theta^{i+1}}$.
This is equivalent
to an R3 transformation of the frame of reference such that
 the signal at time i is along the x-axis in the new frame of reference.
The original coordinates in this frame of reference is given as
\[
i_1 = \cos \theta^i i_2 - \sin \theta^i j_2 
j_1 = \sin \theta^i i_2 + \cos \theta^i j_2 
\]
THe new signal in this reference frame is
\[
   S^{i+1} = S^{i+1}_x i_1  + S^{i+1}_y j_1  
           = S^{i+1}_x ( \cos \theta^i i_2 - \sin \theta^i j_2  )
           + S^{i+1}_y ( \sin \theta^i i_2 + \cos \theta^i j_2  )
\]
and the phase change is given by the angle with the i_2 axis. 
\[
      \delta \phi   = 
    - \delta \theta = atan( 
                           ( S^{i+1}_y \cos \theta^i - S^{i+1}_x \sin \theta^i )
                                                 /
                           ( S^{i+1}_x \cos \theta^i + S^{i+1}_y \sin \theta^i )
                          ) 
                    = atan( 
                           ( S^{i+1}_y S^i_x - S^{i+1}_x S^i_y )
                                                 /
                           ( S^{i+1}_x S^i_x + S^{i+1}_y S^i_y )
                          ) 
                    = atan( conj(S^i) * S^{i+1} ) 
                    = atan2(Im,Re) 
                    = atan2( S^{i+1}_y S^i_x - S^{i+1}_x S^i_y ,
                             S^{i+1}_x S^i_x + S^{i+1}_y S^i_y ) 
\]
to agree with convention
\[
    -atan(Im/Re)  = atan(-Im/Re)  = atan2(-Im,Re) = atan( conj(S^{i+1}) * S^i ) 
\]
 */
  while( !net_It.IsAtEnd() )
    {
     PetscScalar deltaTemp = this->GetTmapFactor() * std::atan2( 
             currIt.Get()[0]*baseIt.Get()[1] - baseIt.Get()[0]*currIt.Get()[1],
             baseIt.Get()[0]*currIt.Get()[0] + baseIt.Get()[1]*currIt.Get()[1]
                                       );
     //std::cout << baseIt.Get()[0] << " " << baseIt.Get()[1] << " " 
     //          << currIt.Get()[0] << " " << currIt.Get()[1] 
     //          << " " << deltaTemp << std::endl;
     snr_It.Set( std::abs(this->GetTmapFactor())*std::sqrt(2) / snr_It.Get() );
     if( std::abs( deltaTemp ) < maxdiff)
       {
        net_It.Set( net_It.Get() + deltaTemp );
        // update the base to contain the n-1 image for next time
        baseIt.Set( currIt.Get() );
        arr_It.Set( arr_It.Get() + m_ActivationEnergy * 
                    std::exp( -m_FrequencyFactor/
                              (m_BaseTemperature+net_It.Get()) / m_GasConstant) 
                  );
       }
     ++currIt;
     ++baseIt;
     ++net_It;
     ++snr_It;
     ++arr_It;
    }

  if(!libMesh::processor_id() ) 
   {
    // output unfiltered tmap
    std::ostringstream unfiltered_filename;
    ImageStats(istep,"net_Image",net_Image);
    unfiltered_filename << OutputDir <<"/temperature.";
    OSSRealzeroright(unfiltered_filename,4,0,istep);
    unfiltered_filename << ".vtk" ;
    WriteImage(net_Image , unfiltered_filename , medianFilterRadius,
                std::string("temperature")  ) ; 

    // output SNR image
    std::ostringstream snrfilename;
    ImageStats(istep,"snr_Image",snr_Image);
    snrfilename << OutputDir <<"/snr.";
    OSSRealzeroright(snrfilename,4,0,istep);
    snrfilename << ".vtk";
    WriteImage(snr_Image , snrfilename , zeroFilterRadius,
                std::string("snr")  ) ; 

    // output damage tmap
    std::ostringstream damage_filename;
    ImageStats(istep,"arr_Image",arr_Image);
    damage_filename << OutputDir <<"/damage.";
    OSSRealzeroright(damage_filename,4,0,istep);
    damage_filename << ".vtk" ;
    WriteImage(arr_Image , damage_filename , medianFilterRadius,
                std::string("damage")  ) ; 

   }

  PetscFunctionReturn(0);
}
// update image stats CALL AFTER UpdateSNR
#undef __FUNCT__
#define __FUNCT__ "Imaging::UpdateImageStats"
PetscErrorCode Imaging::UpdateImageStats(const int istep) 
{
  PetscFunctionBegin;


  // setup iterators
  itk::ImageRegionConstIterator< InputImageType > 
                      net_It(net_Image, net_Image->GetRequestedRegion());
  RegionIteratorType  meanIt(meanImage, meanImage->GetRequestedRegion()),
                      var_It(var_Image, var_Image->GetRequestedRegion());

  // Software Guide : EndCodeSnippet
  // set iterators to the beginning
  net_It.GoToBegin();
  meanIt.GoToBegin();
  var_It.GoToBegin();

  statCount++; // update counter

  //compute net phase difference
  while( !net_It.IsAtEnd() )
    {
     // update pixel wise image statistics
     //  covariance estimates should NOT include heating
     if( statCount  < maxStat)
       {
        meanIt.Set( ((statCount-1)*meanIt.Get()+net_It.Get()) / statCount );
        var_It.Set( ((statCount-2)*var_It.Get() + statCount/(statCount-1) * 
                         (net_It.Get() - meanIt.Get()) * 
                         (net_It.Get() - meanIt.Get()) )/(statCount - 1 ) );
       }
     ++net_It;
     ++meanIt;
     ++var_It;
    }

  // output unfiltered tmap
  if(!libMesh::processor_id() ) 
   {
    std::ostringstream unfiltered_filename;

    ImageStats(istep,"meanImage",meanImage);
    unfiltered_filename << OutputDir <<"/average.";
    OSSRealzeroright(unfiltered_filename,4,0,istep);
    unfiltered_filename << ".vtk" ;
    if( m_Debug ) WriteImage(meanImage, unfiltered_filename ,
                             medianFilterRadius, std::string("average")  ) ; 

    unfiltered_filename.str(""); // reset before reuse
    ImageStats(istep,"var_Image",var_Image);
    unfiltered_filename << OutputDir <<"/variance.";
    OSSRealzeroright(unfiltered_filename,4,0,istep);
    unfiltered_filename << ".vtk" ;
    if( m_Debug ) WriteImage(var_Image, unfiltered_filename , 
                             medianFilterRadius, std::string("variance")  ) ; 
   }

  PetscFunctionReturn(0);
}

// constructor
#undef __FUNCT__
#define __FUNCT__ "ImagingPreProcessed::ImagingPreProcessed"
ImagingPreProcessed::ImagingPreProcessed(const  GetPot &controlfile )
    : ImagingDirect< 1 >( controlfile )
{
}
// update the imaging data structure and snr maps at this time instance
#undef __FUNCT__
#define __FUNCT__ "ImagingPreProcessed::UpdateImaging"
PetscErrorCode ImagingPreProcessed::UpdateImaging(const int istep)
{
  //PetscErrorCode ierr;
  PetscFunctionBegin;

  // setup file readers
  InputReaderType::Pointer  tmapReader = InputReaderType::New(),
                             snrReader = InputReaderType::New();

  // get tmap 
  GetImageFileNameFromCmdLine(istep,std::string("tmap."),
                                    filenames[0][0]);
  tmapReader->SetFileName(filenames[0][0]);
  net_Image = GetImageData< PetscScalar, InputReaderType::Pointer > 
               (tmapReader, libMesh::processor_id() , filenames[0]);
  Point dum;
  RegionIteratorType  net_It(net_Image, net_Image->GetRequestedRegion());
  net_It.GoToBegin();
  // assume reading in temperature change only... need convert to absolute temp
  while( !net_It.IsAtEnd() )
    {
     net_It.Set( net_It.Get() + m_bodyTemp );
     ++net_It;
    }


  // get snr 
  GetImageFileName(istep,std::string("snruncert."),filenames[1][0]);
  snrReader->SetFileName(filenames[1][0]);
  snr_Image = GetImageData< PetscScalar, InputReaderType::Pointer > 
                                        (snrReader, libMesh::processor_id() , filenames[1]);

  // ensure spacing correct
  net_Image->SetSpacing( this->sp   );
  net_Image->SetOrigin(  this->orgn );
  snr_Image->SetSpacing( this->sp   );
  snr_Image->SetOrigin(  this->orgn );

  // output unfiltered tmap
  if( !libMesh::processor_id()  && m_Debug ) 
   {
    ImageStats(istep,"net_Image",net_Image);
    ImageStats(istep,"snr_Image",snr_Image);
    std::ostringstream unfiltered_filename;
    unfiltered_filename << OutputDir <<"/temperature.";
    OSSRealzeroright(unfiltered_filename,4,0,istep);
    unfiltered_filename << ".vtk" ;
    WriteImage(net_Image , unfiltered_filename , zeroFilterRadius,
                                           std::string("temperature")  ) ; 
   }

  PetscFunctionReturn(0);
}
// subroutine to generate dicom filenames
#undef __FUNCT__
#define __FUNCT__ "Imaging::RealTimeGenerateFileNames"
PetscErrorCode 
Imaging::RealTimeGenerateFileNames(const int istep, const int iecho,
                         std::vector < std::vector < std::string > > &filenames)
{
   PetscFunctionBegin;
   // filenames[0][0] ==> real images 
   // filenames[1][0] ==> imaginary images
   {
      std::ostringstream file_name;
      file_name << ExamPath << "/Processed/s" << DirId << "/imageReal."; 
      OSSRealzeroright(file_name,4,0,necho*istep+iecho+noffset);
      file_name << ".vtk";
      filenames[0][0] = file_name.str();
   }
   {
      std::ostringstream file_name;
      file_name << ExamPath << "/Processed/s" << DirId << "/imageImag."; 
      OSSRealzeroright(file_name,4,0,necho*istep+iecho+noffset);
      file_name << ".vtk";
      filenames[1][0] = file_name.str();
   }
   PetscFunctionReturn(0);
}

// subroutine to generate single 3d image filename
#undef __FUNCT__
#define __FUNCT__ "Imaging::GetImageFileName"
PetscErrorCode 
Imaging::GetImageFileName(const int istep, std::string   fileBase, 
                                         std::string  &inputfile)
{
  PetscFunctionBegin;

  std::ostringstream file_name;
  file_name << ExamPath << "/Processed/s" << DirId << "/" << fileBase;
  OSSRealzeroright(file_name,4,0,istep);
  file_name << ".vtk";
  inputfile = file_name.str();

  PetscFunctionReturn(0);
}
// possibly change image path at the command line
#undef __FUNCT__
#define __FUNCT__ "Imaging::GetImageFileNameFromCmdLine"
PetscErrorCode 
Imaging::GetImageFileNameFromCmdLine(const int istep, std::string   fileBase, 
                                     std::string  &inputfile)
{
  PetscFunctionBegin;

  PetscInt OriginalDirId = this->DirId;
  PetscOptionsGetInt(PETSC_NULL,"-originaldirid",
                     &OriginalDirId,PETSC_NULL); 

  std::ostringstream file_name;
  file_name << ExamPath << "/Processed/s" << OriginalDirId << "/" << fileBase;
  OSSRealzeroright(file_name,4,0,istep);
  file_name << ".vtk";
  inputfile = file_name.str();

  PetscFunctionReturn(0);
}
// write an Ini File containg dicom header info 
#undef __FUNCT__
#define __FUNCT__ "Imaging::WriteIni"
PetscErrorCode Imaging::WriteIni()
{
  PetscFunctionBegin;
  std::ofstream Inifile;
  std::ostringstream fileID ; // filename

  // output file
  fileID << OutputDir <<"/mrti.ini" ;
  Inifile.open(fileID.str().c_str(), std::ios::out);

  // dicom header parameters...
  Inifile <<"[mrti]" << std::endl;

  /* dimensions of MRTI data */
  Inifile <<"xpixel=" << size[0] << std::endl ;
  Inifile <<"ypixel=" << size[1] << std::endl ;
  Inifile <<"nslice=" << size[2] << std::endl ;

  /* number of echoes */
  Inifile <<"necho=" << necho << std::endl ;

  /* acquisition time */
  Inifile <<"dt=" << m_acquisition_dt << std::endl ;

  /* physical space dimensions */

  // origin
  Inifile <<"x0=" << orgn[ 0] << std::endl ;
  Inifile <<"y0=" << orgn[ 1] << std::endl ;
  Inifile <<"z0=" << orgn[ 2] << std::endl ;

  // spacing
  Inifile <<"dx=" << sp[0] << std::endl ;
  Inifile <<"dy=" << sp[1] << std::endl ;
  Inifile <<"dz=" << sp[2] << std::endl ;

  // close file 
  Inifile.close();

  PetscFunctionReturn(0);
}

/* Function to get MRTI data */
Number ProjectITKImageData (const Point& p,
                            const Parameters& parameters,
                            const std::string& ,
                            const std::string& )
{
  InterpolatorType::PointType point;
  point[0] = p(0);
  point[1] = p(1);
  point[2] = p(2);

  Imaging *imageData = parameters.get<Imaging*>("ImagingPointer");
  Number value = parameters.get<PetscScalar>("DefaultImageValue");
  if( imageData->m_interpolator->IsInsideBuffer(point) )
    {
     value =  imageData->m_interpolator->Evaluate( point ); 
    }

  return value;
}

// Project Imaging to FEM Mesh
#undef __FUNCT__
#define __FUNCT__ "Imaging::ProjectImagingToFEMMesh"
PetscErrorCode Imaging::ProjectImagingToFEMMesh (System *currentSystem , Vec VecData , double DefaultValue)
  
{
  PetscFunctionBegin; 

  //  The buffer is passed to the ImportImageFilter with the
  //  \code{SetImportPointer()}. Note that the last argument of this method
  //  specifies who will be responsible for deleting the memory block once it is
  //  no longer in use. A \code{false} value indicates that the
  //  ImportImageFilter will not try to delete the buffer when its destructor is
  //  called. A \code{true} value, on the other hand, will allow the filter to
  //  delete the memory block upon destruction of the import filter.
  //
  //  For the ImportImageFilter to appropriately delete the
  //  memory block, the memory must be allocated with the C++
  //  \code{new()} operator.  Memory allocated with other memory
  //  allocation mechanisms, such as C \code{malloc} or \code{calloc}, will not
  //  be deleted properly by the ImportImageFilter.  In
  //  other words, it is the application programmer's responsibility
  //  to ensure that ImportImageFilter is only given
  //  permission to delete the C++ \code{new} operator-allocated memory.
  //  Software Guide : EndLatex

  // Pass petsc Vec pointer to store the setup the image data structure
  PetscInt Nvecsize=0;
  VecGetSize(VecData, &Nvecsize );

  // TODO error check dimensions are the correct
  // if(Nvecsize != m_size[0] * m_size[1] * m_size[2]) 
  //   {
  //     std::cerr << "image dimensions do not match petsc Vec " << std::endl;
  //     libmesh_error();
  //   }

  // petsc will own the data buffer
  InputPixelType *localBuffer; 
  VecGetArray(VecData,&localBuffer);
  const bool importImageFilterWillOwnTheBuffer = false;
  this->m_importFilter->SetImportPointer( localBuffer, Nvecsize, 
                                    importImageFilterWillOwnTheBuffer );

  // set interpolator to project from this image
  this->m_importFilter->Update();
  this->m_interpolator->SetInputImage( 
                               this->m_importFilter->GetOutput() );

  // Project Data onto the system 

  // project_vector is called. data should be stored in 
  // currentSystem.solution & currentSystem.current_local_solution
  Parameters localparameters;
  localparameters.set<PetscScalar>("DefaultImageValue")=DefaultValue;
  localparameters.set<Imaging *>("ImagingPointer")=this;
  currentSystem->project_solution(ProjectITKImageData,NULL,
                                 localparameters);

  // All done.  
  VecRestoreArray(VecData,&localBuffer);
  PetscFunctionReturn(0);
}
