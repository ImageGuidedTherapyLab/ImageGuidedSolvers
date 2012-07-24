/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkVTKImageIO.h,v $
  Language:  C++
  Date:      $Date: 2007-03-22 14:28:53 $
  Version:   $Revision: 1.16 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkVTKImageVariableNameIO_h
#define __itkVTKImageVariableNameIO_h

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include "itkVTKImageIO.h"


namespace itk
{

/** \class VTKImageIO
 *
 *  \brief ImageIO class for reading VTK images
 *
 * \ingroup IOFilters
 *
 */
class ITK_EXPORT VTKImageVariableNameIO : public VTKImageIO 

{
public:
  /** Standard class typedefs. */
  typedef VTKImageVariableNameIO    Self;
  typedef VTKImageIO          Superclass;
  typedef SmartPointer<Self>     Pointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(VTKImageVariableNameIO , Superclass);

  /** Writes the data to disk from the memory buffer provided. Make sure
   * that the IORegion has been set properly. */
  virtual void Write(const void* buffer);

  /** Specify the name of the output file to write. */
  itkSetStringMacro(VariableName);
  itkGetStringMacro(VariableName);
  
protected:
  VTKImageVariableNameIO(); // call superclass
  ~VTKImageVariableNameIO();
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  VTKImageVariableNameIO(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  std::string        m_VariableName;
};

} // end namespace itk

#endif // __itkVTKImageIO_h
