/** \file imiRegistrationBasedInterpolationFilter.txx
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
#ifndef __imiRegistrationBasedInterpolationFilter_txx
#define __imiRegistrationBasedInterpolationFilter_txx

#include "itkNeighborhoodIterator.h"

#include "imiRegistrationBasedInterpolationFilter.h"
#include "itkIterativeInverseDeformationFieldImageFilter.h"
#include "imiImageFunctions.h"

namespace imi
{

//
// Constructor of this filter
template<typename PixelType, unsigned int Dimension>
imiRegistrationBasedInterpolationFilter<PixelType, Dimension>::imiRegistrationBasedInterpolationFilter()
{
  m_bInputWasModified = true;
  m_fInterpolationValue = 0.5;
  m_iInterpolationType = 1;
  m_iWarperType = 1;
}

template<typename PixelType, unsigned int Dimension>
imiRegistrationBasedInterpolationFilter<PixelType, Dimension>::~imiRegistrationBasedInterpolationFilter()
{
  // TODO Auto-generated destructor stub
}

template<typename PixelType, unsigned int Dimension>
void imiRegistrationBasedInterpolationFilter<PixelType, Dimension>::SetInputReferenceImage( InputImageTypePointer referenceImage )
{
  m_bInputWasModified = true;
  m_InputReferenceImage = referenceImage;
}
template<typename PixelType, unsigned int Dimension>
void imiRegistrationBasedInterpolationFilter<PixelType, Dimension>::SetInputTargetImage( InputImageTypePointer targetImage )
{
  m_bInputWasModified = true;
  m_InputTargetImage = targetImage;
}
template<typename PixelType, unsigned int Dimension>
void imiRegistrationBasedInterpolationFilter<PixelType, Dimension>::SetDeformationField( DeformationFieldTypePointer deformationField )
{
  m_bInputWasModified = true;
  m_DeformationField = deformationField;
}

template<typename PixelType, unsigned int Dimension>
void imiRegistrationBasedInterpolationFilter<PixelType, Dimension>::SetInverseDeformationField( DeformationFieldTypePointer invDeformationField )
{
  m_bInputWasModified = true;
  m_InverseDeformationField = invDeformationField;
}

/** set a value between 0 and 1 to interpolate between the Images */
template<typename PixelType, unsigned int Dimension>
void imiRegistrationBasedInterpolationFilter<PixelType, Dimension>::SetInterpolationValue( double alpha )
{
  m_fInterpolationValue = alpha;
}

template<typename PixelType, unsigned int Dimension>
bool imiRegistrationBasedInterpolationFilter<PixelType, Dimension>::Initialize()
{
  if( m_fInterpolationValue < 0 || m_fInterpolationValue > 1.0 )
  {
    imiERROR( "imiRegistrationBasedInterpolationFilter: Interpolation value must be between 0 and 1 !" );
    return false;
  }
  if( !m_bInputWasModified )
  {
    return true;
  }

  //
  // Check, if input was given
  //
  if( m_InputReferenceImage.IsNull() || m_InputTargetImage.IsNull() )
  {
    imiERROR( "imiRegistrationBasedInterpolationFilter: Reference or Target image not set!" );
    return false;
  }
  if( m_DeformationField.IsNull() )
  {
    imiERROR( "imiRegistrationBasedInterpolationFilter: Deformation field not set!" );
    return false;
  }

  //
  // Check size and spacing of input images
  //
  if( !imiImageFunctions::CheckImageParameters<InputImageType, InputImageType>( m_InputReferenceImage, m_InputTargetImage ) )
  {
    imiERROR( "imiRegistrationBasedInterpolationFilter: Input images have not equal size and spacing!" );
    return false;
  }
  if( !imiImageFunctions::CheckImageParameters<InputImageType, DeformationFieldType>( m_InputReferenceImage, m_DeformationField ) )
  {
    imiERROR( "imiRegistrationBasedInterpolationFilter: Input images and deformation field have not equal size and spacing!" );
    return false;
  }
  if( m_InverseDeformationField.IsNotNull() )
  {
    if( !imiImageFunctions::CheckImageParameters<DeformationFieldType, DeformationFieldType>( m_DeformationField, m_InverseDeformationField ) )
    {
      imiERROR( "imiRegistrationBasedInterpolationFilter: Input images and inverse deformation field have not equal size and spacing!" );
      return false;
    }
  }
  else
  {
    //
    //  Generate inverse displacement field, if not available
    //
    imiINFO( "Invert the Deformation-Field (iterative) ..." );
    typedef itk::IterativeInverseDeformationFieldImageFilter<DeformationFieldType, DeformationFieldType> InvertFieldFilterType;

    typename InvertFieldFilterType::Pointer invertFilter = InvertFieldFilterType::New();
    invertFilter->SetNumberOfIterations( 10 );
    invertFilter->SetInput( m_DeformationField );
    invertFilter->Update();

    m_InverseDeformationField = invertFilter->GetOutput();
  }

  return true;
}

template<typename PixelType, unsigned int Dimension>
bool imiRegistrationBasedInterpolationFilter<PixelType, Dimension>::Execute()
{
  if( !Initialize() )
  {
    return false;
  }

  imiDEBUGINFO(5, "Execute registration-based interpolation (interpolation- value="<<m_fInterpolationValue<<" type ="<<m_iInterpolationType<<")");

  //
  // Use a different approach for linear/cubic and nearest-neighbour interpolation
  //
  if( InterpolationTypeIsLinear() || InterpolationTypeIsCubicBSpline() )
  {
    // Linear/cubic interpolation:
    // Warp target image with scaled deformation field towards reference image
    // (scaling factor is 1.0 - interpolationValue)
    //
    if( !imiImageFunctions::WarpImageWithFactorND( m_DeformationField, m_InputTargetImage, m_OutputImage, 1.0 - m_fInterpolationValue, m_iInterpolationType, m_iWarperType ) )
    {
      imiERROR( "imiRegistrationBasedInterpolationFilter: Warping input target image failed!" );
      return false;
    }

    OutputImageTypePointer invWarpedImage;

    //
    // Warp reference image with scaled inverse deformation field towards target image
    // (scaling factor is interpolationValue)
    //
    if( !imiImageFunctions::WarpImageWithFactorND( m_InverseDeformationField, m_InputReferenceImage, invWarpedImage, m_fInterpolationValue, m_iInterpolationType, m_iWarperType ) )
    {
      imiERROR( "imiRegistrationBasedInterpolationFilter: Warping input reference image failed!" );
      return false;
    }

    //
    // Weighted add of the corresponding pixel values of the warped images
    // interpolationValue * warped target + (1 - interpolationValue) * warped reference
    //
    typedef itk::ImageRegionIterator<OutputImageType> IteratorType;
    typedef typename OutputImageType::ValueType ValueType;

    IteratorType outputImageIt( m_OutputImage, m_OutputImage->GetLargestPossibleRegion() );
    IteratorType invWarpedImageIt( invWarpedImage, invWarpedImage->GetLargestPossibleRegion() );

    outputImageIt.GoToBegin();
    invWarpedImageIt.GoToBegin();

    const double interpValue = m_fInterpolationValue;
    const double invInterpValue = 1.0 - m_fInterpolationValue;
    unsigned long int count = 0;

    ValueType pixelValue;
    while( !outputImageIt.IsAtEnd() )
    {

      pixelValue = static_cast<ValueType>( outputImageIt.Get() * interpValue + invWarpedImageIt.Get() * invInterpValue );
      outputImageIt.Set( pixelValue );

      ++outputImageIt;
      ++invWarpedImageIt;
      count++;
    }
    imiDEBUGINFO(9, "Combined "<<count<<" Pixels (linear)." );
  }
  else if( InterpolationTypeIsNearestNeighbour() )
  {
    // NEAREST-NEIGHBOUR interpolation:
    // Warp only the nearest image with scaled (inverse) deformation field
    //
    // Set output value to the values in the warped image (no linear interpolation for label images)
    //
    if( m_fInterpolationValue > 0.5 )
    {
      // Warp target image with scaled deformation field towards reference image
      // (scaling factor is 1.0 - interpolationValue)
      if( !imiImageFunctions::WarpImageWithFactorND( m_DeformationField, m_InputTargetImage, m_OutputImage, 1.0 - m_fInterpolationValue, m_iInterpolationType, m_iWarperType ) )
      {
        imiERROR( "imiRegistrationBasedInterpolationFilter: Warping input target image failed!" );
        return false;
      }
    }
    else
    {
      // Warp reference image with scaled inverse deformation field towards target image
      // (scaling factor is interpolationValue)
      if( !imiImageFunctions::WarpImageWithFactorND( m_InverseDeformationField, m_InputReferenceImage, m_OutputImage, m_fInterpolationValue, m_iInterpolationType, m_iWarperType ) )
      {
        imiERROR( "imiRegistrationBasedInterpolationFilter: Warping input reference image failed!" );
        return false;
      }
    }
  }
  else if( InterpolationTypeIsLabelGaussian() )
  {
    // GAUSSIAN LABEL interpolation:
    //
    // Warp target image with scaled deformation field towards reference image
    // (scaling factor is 1.0 - interpolationValue)
    //
    OutputImageTypePointer warpedImage;

    if( !imiImageFunctions::WarpImageWithFactorND( m_DeformationField, m_InputTargetImage, warpedImage, 1.0 - m_fInterpolationValue, m_iInterpolationType, m_iWarperType ) )
    {
      imiERROR( "imiRegistrationBasedInterpolationFilter: Warping input target image failed!" );
      return false;
    }

    OutputImageTypePointer invWarpedImage;

    //
    // Warp reference image with scaled inverse deformation field towards target image
    // (scaling factor is interpolationValue)
    //
    if( !imiImageFunctions::WarpImageWithFactorND( m_InverseDeformationField, m_InputReferenceImage, invWarpedImage, m_fInterpolationValue, m_iInterpolationType, m_iWarperType ) )
    {
      imiERROR( "imiRegistrationBasedInterpolationFilter: Warping input reference image failed!" );
      return false;
    }

    m_OutputImage = imiImageFunctions::AllocateImage(warpedImage);
    if(m_OutputImage.IsNull())
    {
      imiERROR("imiRegistrationBasedInterpolationFilter: Can not allocate output image.");
      return false;
    }

    typedef itk::NeighborhoodIterator<OutputImageType> IteratorType;
    typedef typename OutputImageType::PixelType OutPixelType;
    typedef typename itk::NumericTraits<OutPixelType>::RealType RealType;

    // Set the radius for the iterator
    typename IteratorType::RadiusType labelRadius;
    labelRadius.Fill( 1 );

    itk::ImageRegionIterator<OutputImageType> outputImageIt( m_OutputImage, m_OutputImage->GetLargestPossibleRegion() );
    IteratorType warpedImageIt( labelRadius, warpedImage, warpedImage->GetLargestPossibleRegion() );
    IteratorType invWarpedImageIt( labelRadius, invWarpedImage, invWarpedImage->GetLargestPossibleRegion() );

    outputImageIt.GoToBegin();
    warpedImageIt.GoToBegin();
    invWarpedImageIt.GoToBegin();

    typedef std::less<RealType> ComparatorType;
    typedef std::map<OutPixelType, RealType, ComparatorType> WeightMapType;
    typedef typename std::map<OutPixelType, RealType, ComparatorType>::iterator WeightMapIteratorType;
    WeightMapType weightMap;

    unsigned long int count = 0;
    while( !outputImageIt.IsAtEnd() )
    {
      const OutPixelType warpedCenterValue = warpedImageIt.GetCenterPixel();
      const OutPixelType invwarpedCenterValue = invWarpedImageIt.GetCenterPixel();

      OutPixelType maxValue = warpedCenterValue;
      if( warpedCenterValue != invwarpedCenterValue )
      {
        RealType maxCount = 0.0;
        weightMap.clear();

        const unsigned int hoodSize = warpedImageIt.Size();
        for( unsigned int indct = 0; indct < hoodSize; indct++ )
        {
          {
            const OutPixelType warpedValue = warpedImageIt.GetPixel( indct );

            WeightMapIteratorType it = weightMap.find( warpedValue );
            RealType wtest = 0.0;
            const RealType weight = m_fInterpolationValue;
            if( it != weightMap.end() )
            {
              it->second += weight;
              wtest = it->second;
            }
            else
            {
              weightMap.insert( std::make_pair( warpedValue, weight ) );
              wtest = weight;
            }

            //Keep track of the max value
            if( wtest > maxCount )
            {
              maxCount = wtest;
              maxValue = warpedValue;
            }
          }

          {
            const OutPixelType invWarpedValue = invWarpedImageIt.GetPixel( indct );
            WeightMapIteratorType it = weightMap.find( invWarpedValue );
            RealType wtest = 0.0;
            const RealType weight = 1.0 - m_fInterpolationValue;
            if( it != weightMap.end() )
            {
              it->second += weight;
              wtest = it->second;
            }
            else
            {
              weightMap.insert( std::make_pair( invWarpedValue, weight ) );
              wtest = weight;
            }

            //Keep track of the max value
            if( wtest > maxCount )
            {
              maxCount = wtest;
              maxValue = invWarpedValue;
            }
          }
        } //for( unsigned int indct = 0; indct < hoodSize; indct++ )
      } //if( warpedValue != invwarpedValue )
      //
      // Set value with maximum count
      outputImageIt.Set( maxValue );

      ++outputImageIt;
      ++warpedImageIt;
      ++invWarpedImageIt;
      count++;
    } // while( !outputImageIt.IsAtEnd() )
    imiDEBUGINFO(9, "Combined "<<count<<" Pixels (gaussian label)." );

  } // else if( InterpolationTypeIsLabelGaussian() )
  else
  {
    imiERROR( "imiRegistrationBasedInterpolationFilter: Unknown InterpolationType!" );
    return false;
  }
  m_bInputWasModified = false;
  return true;
}

} /* namespace imi */

#endif // __imiRegistrationBasedInterpolationFilter_txx
