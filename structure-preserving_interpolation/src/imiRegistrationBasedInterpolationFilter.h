/** \file imiRegistrationBasedInterpolationFilter.h
 *
 *  <small> <!--Copyright Information: -->
 *  \b Author: Jan Ehrhardt \n
 *  \b Copyright (C) 2010, Jan Ehrhardt, Institute of Medical Informatics,
 *     University of Luebeck, Germany. All rights reserved.\n
 *     Please cite the following paper:\n
 *     Ehrhardt, J., Säring, D., & Handels, H. (2007). Structure-preserving
 *     interpolation of temporal and spatial image sequences using an optical
 *     flow-based method. Methods of information in medicine, 46(03), 300-307.\n *
 *  </small>
 ****************************************************************************/
/*
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 * contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef IMIREGISTRATIONBASEDINTERPOLATIONFILTER_H_
#define IMIREGISTRATIONBASEDINTERPOLATIONFILTER_H_

#include "itkImage.h"

#include "imiMacro.h"

namespace imi
{

template <typename PixelType, unsigned int Dimension=3>
class imiRegistrationBasedInterpolationFilter
{
public:
  imiRegistrationBasedInterpolationFilter();
  virtual ~imiRegistrationBasedInterpolationFilter();

  //-----------------------------------
  // Typedefs:
  //-----------------------------------
  typedef itk::Image<PixelType,Dimension> InputImageType;
  typedef typename InputImageType::Pointer InputImageTypePointer;
  typedef itk::Image<PixelType,Dimension> OutputImageType;
  typedef typename OutputImageType::Pointer OutputImageTypePointer;
  typedef itk::Vector< float, Dimension>    VectorPixelType;
  typedef itk::Image<VectorPixelType,Dimension> DeformationFieldType;
  typedef typename DeformationFieldType::Pointer DeformationFieldTypePointer;

  void SetInputReferenceImage(InputImageTypePointer referenceImage);
  void SetInputTargetImage(InputImageTypePointer referenceImage);
  void SetDeformationField(DeformationFieldTypePointer deformationField);
  void SetInverseDeformationField(DeformationFieldTypePointer invDeformationField);

  /** set a value between 0 and 1 to interpolate between the Images */
  void  SetInterpolationValue (double alpha);
  double GetInterpolationValue (void) const {return m_fInterpolationValue;};

  void SetInterpolationTypeToLinear() { m_iInterpolationType=1; };
  void SetInterpolationTypeToNearestNeighbour() { m_iInterpolationType=0; };
  void SetInterpolationTypeToCubicBSpline() { m_iInterpolationType=2; };
  void SetInterpolationTypeToLabelGaussian() { m_iInterpolationType=3; };

  bool InterpolationTypeIsLinear() const { return (m_iInterpolationType==1); };
  bool InterpolationTypeIsNearestNeighbour() const { return (m_iInterpolationType==0); };
  bool InterpolationTypeIsCubicBSpline() const { return (m_iInterpolationType==2); };
  bool InterpolationTypeIsLabelGaussian() const { return (m_iInterpolationType==3); };

  void SetWarperTypeToITKStandard() { m_iWarperType=0; };
  void SetWarperTypeToContinuousBorder() { m_iWarperType=1; };
  void SetWarperType(int type) { m_iWarperType=type; };

  bool WarperTypeIsITKStandard() const { return (m_iWarperType==0); };
  bool WarperTypeIsContinuousBorder() const { return (m_iWarperType==1); };
  int GetWarperType() { return m_iWarperType; };

  bool Execute();

  OutputImageTypePointer GetOutputImage() {return m_OutputImage;};

protected:
  bool m_bInputWasModified;
  double m_fInterpolationValue;
  int m_iInterpolationType;
  int m_iWarperType;
  InputImageTypePointer m_InputReferenceImage;
  InputImageTypePointer m_InputTargetImage;
  DeformationFieldTypePointer m_DeformationField;
  DeformationFieldTypePointer m_InverseDeformationField;
  OutputImageTypePointer m_OutputImage;

  bool Initialize();
};

} /* namespace imi */

#ifndef ITK_MANUAL_INSTANTIATION
#include "imiRegistrationBasedInterpolationFilter.hxx"
#endif

#endif /* IMIREGISTRATIONBASEDINTERPOLATIONFILTER_H_ */
