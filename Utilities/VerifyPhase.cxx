/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkStatisticsImageFilterTest.cxx,v $
  Language:  C++
  Date:      $Date: 2008-02-27 19:40:03 $
  Version:   $Revision: 1.11 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include <itksys/SystemTools.hxx>
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesWriter.h"

#include "petsc.h"
// globals 
const double                       pi = 3.1415926535897931;
static char help[] = "Imaging for Debugging\n\n";

// useful type defs
typedef itk::Image<short int,2> Image2DType;
typedef itk::Image<short int,3> ImageType;

// write a dicom and an mha image
void WriteImage( ImageType::Pointer image, std::ostringstream &imagefilename, 
                                       std::vector < std::string > &filenames )
{

   PetscFunctionBegin;

   std::cout << "writing "  << imagefilename.str() << std::endl;

   typedef itk::ImageSeriesWriter< 
                              ImageType, Image2DType >  SeriesWriterType;

   //  We construct a series writer and connect to its input the output from the
   //  reader. Then we pass the GDCM image IO object in order to be able to
   //  write the images in DICOM format.  

   SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();
   typedef itk::GDCMImageIO                        ImageIOType;
   ImageIOType::Pointer gdcmIO = ImageIOType::New();

   seriesWriter->SetInput( image );
   seriesWriter->SetImageIO( gdcmIO );
   seriesWriter->SetFileNames( filenames );

   //  The following line of code is extremely important for this process to
   //  work correctly.  The line is taking the MetaDataDictionary from the input
   //  reader and passing it to the output writer. The reason why this step is
   //  so important is that the MetaDataDictionary contains all the entries of
   //  the input DICOM header.
   //
   //  \index{itk::ImageSeriesReader!GetMetaDataDictionaryArray()}
   //  \index{itk::ImageSeriesWriter!SetMetaDataDictionaryArray()}
   //
   //seriesWriter->SetMetaDataDictionaryArray( 
   //                      image->GetMetaDataDictionaryArray() );

   // Finally we trigger the writing process by invoking the \code{Update()} method
   // in the series writer. We place this call inside a try/catch block, in case
   // any exception is thrown during the writing process.
   try
     {
     seriesWriter->Update();
     }
   catch( itk::ExceptionObject & excp )
     {
     std::cerr << "Exception thrown while writing the series " << std::endl;
     std::cerr << excp << std::endl;
     return;
     }

   // Software Guide : EndCodeSnippet
   //
   // We create a second writer object that will save the rescaled image into a
   // file. This time not in DICOM format. This is done only for the sake of
   // verifying the image against the one that will be saved in DICOM format
   // later on this example.
   //
   // Software Guide : EndLatex 
   
   // Software Guide : BeginCodeSnippet
   typedef itk::ImageFileWriter< ImageType >  WriterType;
   WriterType::Pointer writer2 = WriterType::New();
   
   writer2->SetFileName( imagefilename.str() );
   writer2->SetInput( image );

   // Software Guide : EndCodeSnippet
   
   
   // Software Guide : BeginLatex
   //
   // The writer can be executed by invoking the Update() method from inside a
   // try/catch block.
   //
   // Software Guide : EndLatex 
   try
     {
     writer2->Update();
     }
   catch (itk::ExceptionObject & e)
     {
     std::cerr << "exception in file writer " << std::endl;
     std::cerr << e << std::endl;
     }

   PetscFunctionReturnVoid();
}

// main driver routine
int main( int argc, char* argv[] )
{
  PetscErrorCode ierr;
  /* Initialize Petsc */
  ierr = PetscInitialize(&argc,&argv,(char*)0,help);CHKERRQ(ierr); 
  PetscFunctionBegin;

  // defaults
  PetscInt ntime = 10, nslice= 1, npixel=4; 

  //error check
  if( argc < 2 )
    {
     std::cerr << "Usage: " << argv[0] <<
                  " OutputDir [options] "   << std::endl;
     std::cerr << "Available options: "     << std::endl ; 
     std::cerr << "-ntime     Number of Time Points: "     << std::endl ;
     std::cerr << "-nslice    Number of slices: "          << std::endl ;
     std::cerr << "-npixel    Number of pixels: "          << std::endl ;
     return EXIT_FAILURE;
    }

  const char *OutputDir  =    argv[1] ;
  int DirId = 12345; 

  std::ostringstream outputPath ; // filename
  outputPath << OutputDir << "/knownphase/s" << DirId ;
  itksys::SystemTools::MakeDirectory( outputPath.str().c_str() );

  // get options
  ierr=PetscOptionsGetInt(PETSC_NULL,"-ntime",&ntime,PETSC_NULL);
  CHKERRQ(ierr);

  ierr=PetscOptionsGetInt(PETSC_NULL,"-nslice",&nslice,PETSC_NULL);
  CHKERRQ(ierr);

  ierr=PetscOptionsGetInt(PETSC_NULL,"-npixel",&npixel,PETSC_NULL);
  CHKERRQ(ierr);

  std::cout << "itkPhaseImages " << "ntime "<< ntime 
            << "npixel" << npixel << "nslice" << nslice << std::endl;

  ImageType::Pointer    image  = ImageType::New();
  ImageType::RegionType region;
  ImageType::SizeType   size; size[0]=npixel; size[1]=npixel; size[2]=nslice;
  ImageType::IndexType  index; index.Fill(0);

  region.SetIndex (index);
  region.SetSize (size);

  image->SetRegions( region );
  ImageType::SpacingType sp;
  sp[0] = 0.9375 ;
  sp[1] = 0.9375 ;
  sp[2] =   5.0  ;
  image->SetSpacing( sp  );
  ImageType::PointType orgn;
  orgn[0] = 0.0 ;
  orgn[1] = 0.0 ;
  orgn[2] = 0.0 ;
  image->SetOrigin(  orgn);
  image->Allocate();

  // set some params of the DicomHeader
  typedef itk::MetaDataDictionary   DictionaryType;
  DictionaryType & dictionary = image->GetMetaDataDictionary();

  //std::string pixelspacekey = "0028|0030";
  std::string slicethickkey = "0018|0050";
  //std::string imagfreqkey   = "0018|0084";
  //std::string echotimekey   = "0018|0081";

  //itk::EncapsulateMetaData<std::string>( dictionary, pixelspacekey, "0.9375" );
  // slice thickness still not showing up w/ DicomImageReadPrintTags but
  // seems to be working none the less...
  std::ostringstream slicethickvalue;
  slicethickvalue << sp[2];
  itk::EncapsulateMetaData<std::string>( dictionary, slicethickkey, slicethickvalue.str() );
  //itk::EncapsulateMetaData<std::string>( dictionary, imagfreqkey  , value );
  //itk::EncapsulateMetaData<std::string>( dictionary, echotimekey  , value );

  double dtheta = 5.0;
 
  // container for filenames
  std::vector < std::string >  filenames(nslice,"");

  for(int istep = 0; istep <= ntime; istep++)
    {
     // real image 
     double realFillValue = 1000.;
     image->FillBuffer( static_cast< ImageType::PixelType >( realFillValue ) );

     // write real image
     for( int jjj = 0 ; jjj < nslice ; jjj++)
       {
        std::ostringstream filename2d;
        filename2d << OutputDir << "/knownphase/s" << DirId << "/i" 
                  <<   DirId     + 2*nslice*istep + (2*jjj+1) 
                  << ".MRDC."<<    2*nslice*istep + (2*jjj+1)  ;
        filenames[jjj] = filename2d.str();
       }
     std::ostringstream realfile_name;
     realfile_name << OutputDir << "/knownphase/s" 
                   << DirId << "/real"<< istep << ".mha";
     WriteImage(image,realfile_name,filenames);

     // imaginary image 
     double imagFillValue = realFillValue * std::tan(istep * dtheta * pi/180.0);

     image->FillBuffer( static_cast< ImageType::PixelType >( imagFillValue ) );

     // write imaginary image
     for( int jjj = 0 ; jjj < nslice ; jjj++)
       {
        std::ostringstream filename2d;
        filename2d << OutputDir << "/knownphase/s" << DirId << "/i" 
                   <<   DirId     + 2*nslice*istep + (2*jjj+1) + 1
                   << ".MRDC."<<    2*nslice*istep + (2*jjj+1) + 1 ;
        filenames[jjj] = filename2d.str();
       }
     std::ostringstream imagfile_name;
     imagfile_name << OutputDir << "/knownphase/s" 
                   << DirId << "/imag"<< istep << ".mha";
     WriteImage(image,imagfile_name,filenames);

    }

  ierr = PetscFinalize();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

