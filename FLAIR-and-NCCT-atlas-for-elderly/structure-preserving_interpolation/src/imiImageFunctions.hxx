/** \file imiImageFunctions.txx
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
// ITK includes:
#include "itkContinuousBorderWarpImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkLabelImageGaussianInterpolateImageFunction.h"

#include "itkMultiplyImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkCastImageFilter.h"

// Project includes.
#include "imiMacro.h"
#include "imiImageFunctions.h"

namespace imi
{
int imiImageFunctions::m_iGlobalDebugLevel = 1;

///////////////////////////////////////////////////////////////////////////////////////////
//
// CheckImageParameters
//
///////////////////////////////////////////////////////////////////////////////////////////
template<typename TImageType1, typename TImageType2>
bool imiImageFunctions::CheckImageParameters( itk::SmartPointer<TImageType1> inputImage1, itk::SmartPointer<TImageType2> inputImage2, bool checkSize, bool checkSpacing, bool checkOrigin )
{
  typename TImageType1::ConstPointer ptr1 = inputImage1.GetPointer();
  typename TImageType2::ConstPointer ptr2 = inputImage2.GetPointer();
  return CheckImageParameters(ptr1,ptr2);
}
///////////////////////////////////////////////////////////////////////////////////////////
//
// CheckImageParameters
//
///////////////////////////////////////////////////////////////////////////////////////////
template<typename TImageType1, typename TImageType2>
bool imiImageFunctions::CheckImageParameters( itk::SmartPointer<const TImageType1> inputImage1, itk::SmartPointer<const TImageType2> inputImage2, bool checkSize, bool checkSpacing, bool checkOrigin )
{
  const double TOLERANCE = 1e-4;
  bool dataIsEqual = true;

  // Check if target and reference are set
  if( inputImage1.IsNull() || inputImage2.IsNull() )
  {
    imiERROR( "One of the images is NULL!" );
    return false;
  }

  // Check if target and reference have same size, spacing and origin.
  if( checkSize )
  {
    typename TImageType1::RegionType::SizeType size1 = inputImage1->GetRequestedRegion().GetSize();
    typename TImageType2::RegionType::SizeType size2 = inputImage2->GetRequestedRegion().GetSize();
    if( size1 != size2 )
    {
      imiINFO( "Images do not have the same size!" );
      imiINFO( "        Data1: " << size1 );
      imiINFO( "        Data2: " << size2 );
      dataIsEqual = false;
    }
  }
  if( checkOrigin )
  {
    typename TImageType1::PointType origin1 = inputImage1->GetOrigin();
    typename TImageType2::PointType origin2 = inputImage2->GetOrigin();

    bool originIsEqual = true;
    for( unsigned int dim = 0; dim < TImageType1::GetImageDimension(); ++dim )
    {
      if( vnl_math_abs( origin1[dim] - origin2[dim] ) > TOLERANCE )
      {
        originIsEqual = false;
      }
    }

    if( !originIsEqual )
    {
      imiINFO( "Images do not have the same origin!" );
      imiINFO( "        Data1: " << origin1 );
      imiINFO( "        Data2: " << origin2 );
      dataIsEqual = false;
    }
  }
  if( checkSpacing )
  {
    typename TImageType1::SpacingType spacing1 = inputImage1->GetSpacing();
    typename TImageType2::SpacingType spacing2 = inputImage2->GetSpacing();

    bool spacingIsEqual = true;
    for( unsigned int dim = 0; dim < TImageType1::GetImageDimension(); ++dim )
    {
      if( vnl_math_abs( spacing1[dim] - spacing2[dim] ) > TOLERANCE )
      {
        spacingIsEqual = false;
      }
    }

    if( !spacingIsEqual )
    {
      imiINFO( "Images do not have the same spacing!" );
      imiINFO( "        Data1: " << spacing1 );
      imiINFO( "        Data2: " << spacing2 );
      dataIsEqual = false;
    }
  }
  return dataIsEqual;
}
///////////////////////////////////////////////////////////////////////////////////////////
//
// CheckImageParameters
//
///////////////////////////////////////////////////////////////////////////////////////////
template<typename TImageType1, typename TImageType2>
bool imiImageFunctions::CheckImageParametersQuiet( itk::SmartPointer<TImageType1> inputImage1, itk::SmartPointer<TImageType2> inputImage2, bool checkSize, bool checkSpacing, bool checkOrigin )
{
  const double TOLERANCE = 1e-4;
  bool dataIsEqual = true;

  // Check if target and reference are set
  if( inputImage1.IsNull() || inputImage2.IsNull() )
  {
    return false;
  }

  // Check if target and reference have same size, spacing and origin.
  if( checkSize )
  {
    typename TImageType1::RegionType::SizeType size1 = inputImage1->GetRequestedRegion().GetSize();
    typename TImageType2::RegionType::SizeType size2 = inputImage2->GetRequestedRegion().GetSize();
    if( size1 != size2 )
    {
      return false;
    }
  }
  if( checkOrigin )
  {
    typename TImageType1::PointType origin1 = inputImage1->GetOrigin();
    typename TImageType2::PointType origin2 = inputImage2->GetOrigin();

    for( unsigned int dim = 0; dim < TImageType1::GetImageDimension(); ++dim )
    {
      if( vnl_math_abs( origin1[dim] - origin2[dim] ) > TOLERANCE )
      {
        return false;
      }
    }
  }
  if( checkSpacing )
  {
    typename TImageType1::SpacingType spacing1 = inputImage1->GetSpacing();
    typename TImageType2::SpacingType spacing2 = inputImage2->GetSpacing();

    for( unsigned int dim = 0; dim < TImageType1::GetImageDimension(); ++dim )
    {
      if( vnl_math_abs( spacing1[dim] - spacing2[dim] ) > TOLERANCE )
      {
        return false;
      }
    }
  }
  return true;
}

////////////////////////////////////////////////////////////////
//
// WarpImage() full templated
//
////////////////////////////////////////////////////////////////
template<typename TImageType, class TDisplacementFieldType>
bool imiImageFunctions::WarpImageND( itk::SmartPointer<TDisplacementFieldType> displField, itk::SmartPointer<TImageType> inputImage, itk::SmartPointer<TImageType> &outputImage, int InterpolationType, int WarpFilterType )
{

  ///////////////////////////
  // Setup warper.
  imiDEBUGINFO( 7, "   Warping "<<TImageType::GetImageDimension()<<"D-image..." );

  // Generate the warper.
  typedef itk::WarpImageFilter<TImageType, TImageType, TDisplacementFieldType> WarpImageFilterType;
  typedef itk::ContinuousBorderWarpImageFilter<TImageType, TImageType, TDisplacementFieldType> ContinuousBorderWarpImageFilterType;

  typename WarpImageFilterType::Pointer warpFilter;
  switch( WarpFilterType )
  {
    case 0:
      warpFilter = WarpImageFilterType::New();
      break;
    default:
      warpFilter = ContinuousBorderWarpImageFilterType::New();
      break;
  }

  // Setup input.
  warpFilter->SetInput( inputImage );
  warpFilter->SetOutputParametersFromImage( displField );
  warpFilter->SetDeformationField( displField );
  warpFilter->SetEdgePaddingValue( 0 );

  //
  // Instantiate the interpolator
  typedef itk::InterpolateImageFunction<TImageType, double> InterpolatorType;
  typename InterpolatorType::Pointer interpolator = GetInterpolateImageFunctionByType<TImageType>( InterpolationType );
  warpFilter->SetInterpolator( interpolator );

  // Execute warper.
  try
  {
    warpFilter->Update();
  }
  catch( itk::ExceptionObject &err )
  {
    imiERROR( "Cannot warp image! Exception error: " << err );
    return false;
  }

  outputImage = warpFilter->GetOutput();

  return true;
}

////////////////////////////////////////////////////////////////
//
// WarpImageWithFactor() full templated
//
////////////////////////////////////////////////////////////////
template<class TImageType, class TDisplacementFieldType>
bool imiImageFunctions::WarpImageWithFactorND( itk::SmartPointer<TDisplacementFieldType> displField, itk::SmartPointer<TImageType> inputImage, itk::SmartPointer<TImageType> &outputImage, double factor, int InterpolationType, int WarpFiltertype )
{
  imiDEBUGINFO( 7, "   Warping image with factor " << factor << "..." );

  // Use factor only if necessary. In many cases
  // the factor is one so this will be skipped
  typename TDisplacementFieldType::Pointer scaledDisplField;
  if( vcl_fabs( factor - 1.0 ) > 1.0e-4 )
  {
    typedef itk::MultiplyImageFilter<TDisplacementFieldType, itk::Image<double,TDisplacementFieldType::ImageDimension >, TDisplacementFieldType> MultiplyFilterType;
    typename MultiplyFilterType::Pointer multiplier = MultiplyFilterType::New();
    multiplier->SetConstant( factor );
    multiplier->SetInput( displField );
    try
    {
      multiplier->Update();
    }
    catch( itk::ExceptionObject &err )
    {
      imiERROR( "Cannot scale deformation field: " << err );
      return false;
    }

    scaledDisplField = multiplier->GetOutput();
  }
  else
  {
    scaledDisplField = displField;
  }

  return WarpImageND( scaledDisplField, inputImage, outputImage, InterpolationType, WarpFiltertype );
}

////////////////////////////////////////////////////////////////
//
// ReadImage()
//
////////////////////////////////////////////////////////////////
template<typename TImageType>
bool imiImageFunctions::ReadImage(
    std::string filename,
    itk::SmartPointer<TImageType> &pImage)
{
  typedef itk::ImageFileReader<TImageType> ReaderType;
  typename ReaderType::Pointer imageReader = ReaderType::New();
  try
  {

    imageReader->SetFileName( filename );

    imageReader->Update();
  }
  catch( itk::ExceptionObject & excep )
  {
    imiERROR( "Reading failed with exception:" << excep << std::endl );
    return false;
  }
  pImage = imageReader->GetOutput();

  return true;
}
////////////////////////////////////////////////////////////////
//
// WriteImage()
//
////////////////////////////////////////////////////////////////
template<typename TImageType>
bool imiImageFunctions::WriteImage( itk::SmartPointer<TImageType> pImage, std::string filename )
{
  if( pImage.IsNull() )
  {
    imiERROR( "Can not write NULL image!" );
    return false;
  }

  try
  {
    typedef itk::ImageFileWriter<TImageType> WriterType;
    typename WriterType::Pointer imageWriter = WriterType::New();

    imageWriter->SetInput( pImage );
    imageWriter->SetFileName( filename );

    imageWriter->Update();
  }
  catch( itk::ExceptionObject & excep )
  {
    imiERROR( "Writing failed with exception:" << excep << std::endl );
    return false;
  }
  return true;
}


/** allocates an image with same TYPE, same size, spacing origin etc as the reference image */
template<typename TImageType>
typename TImageType::Pointer imiImageFunctions::AllocateImage( itk::SmartPointer<TImageType> referenceImage )
{
  typename TImageType::Pointer outputImage = TImageType::New();
  outputImage->CopyInformation( referenceImage );
  outputImage->SetRegions( referenceImage->GetLargestPossibleRegion() );
  outputImage->Allocate();

  return outputImage;
}

/** allocates an image with different TYPE, but same size, spacing origin etc as the reference image */
template<typename TInputImageType, typename TOutputImageType>
void imiImageFunctions::AllocateImage( itk::SmartPointer<TInputImageType> referenceImage, itk::SmartPointer<TOutputImageType> &outputImage)
{
  outputImage = TOutputImageType::New();
  outputImage->CopyInformation( referenceImage );
  outputImage->SetRegions( referenceImage->GetLargestPossibleRegion() );
  outputImage->Allocate();
}

/** allocates an image with different TYPE, but same size, spacing origin etc as the reference image */
template<typename TInputImageType, typename TOutputImageType>
void imiImageFunctions::AllocateVectorImage( itk::SmartPointer<TInputImageType> referenceImage, int vectorLength, itk::SmartPointer<TOutputImageType> &outputImage)
{
  outputImage = TOutputImageType::New();
  outputImage->CopyInformation( referenceImage );
  outputImage->SetNumberOfComponentsPerPixel( vectorLength );
  outputImage->SetRegions( referenceImage->GetLargestPossibleRegion() );
  outputImage->Allocate();
}

template<class TImageType>
typename itk::InterpolateImageFunction<TImageType, double>::Pointer imiImageFunctions::GetInterpolateImageFunctionByType( int InterpolationType )
{
  // Generate the interpolator.
  typedef itk::InterpolateImageFunction<TImageType, double> InterpolatorType;
  typedef itk::NearestNeighborInterpolateImageFunction<TImageType, double> NNInterpolatorType;
  typedef itk::LinearInterpolateImageFunction<TImageType, double> LInterpolatorType;
  typedef itk::BSplineInterpolateImageFunction<TImageType, double> BSInterpolatorType;
  typedef itk::LabelImageGaussianInterpolateImageFunction<TImageType, double> GaussianInterpolatorType;

  typename InterpolatorType::Pointer interpolator;
  switch( InterpolationType )
  {
    case 0:
      imiDEBUGINFO( 7, "   Using nearest neighbor interpolation." );
      interpolator = NNInterpolatorType::New();
      break;
    case 1:
      imiDEBUGINFO( 7, "   Using linear interpolation." );
      interpolator = LInterpolatorType::New();
      break;
    case 2:
      imiDEBUGINFO( 7, "   Using B-Spline interpolation." );
      interpolator = BSInterpolatorType::New();
      break;
    case 3:
      imiDEBUGINFO( 7, "   Using Gaussian label interpolation." );
      interpolator = GaussianInterpolatorType::New();
      break;
    default:
      imiWARNING( "imiImageFunctions::WarpImageND(): unknown interpolator type - using linear." );
      interpolator = LInterpolatorType::New();
      break;
  }

  return interpolator;
}

}
