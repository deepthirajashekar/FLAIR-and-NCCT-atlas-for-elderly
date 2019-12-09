/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkCastImageFilter.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkRoundAndCastImageFilter_h
#define __itkRoundAndCastImageFilter_h

#include "itkUnaryFunctorImageFilter.h"
#include "itkProgressReporter.h"

#include "vnl/vnl_math.h"

namespace itk
{

/** \class RoundAndCastImageFilter
 *
 * \brief Casts input pixels to output pixel type. The values are rounded before casting.
 *
 * This filter is templated over the input image type
 * and the output image type.
 *
 * If you need to perform a dimensionaly reduction, you may want
 * to use the ExtractImageFilter instead of the CastImageFilter.
 *
 * \ingroup IntensityImageFilters  Multithreaded
 * \sa UnaryFunctorImageFilter
 * \sa ExtractImageFilter
 */
namespace Functor {

template< class TInput, class TOutput>
class RoundAndCast
{
public:
	RoundAndCast() {};
  virtual ~RoundAndCast() {};
  bool operator!=( const RoundAndCast & ) const
    {
    return false;
    }
  bool operator==( const RoundAndCast & other ) const
    {
    return !(*this != other);
    }
  inline TOutput operator()( const TInput & A ) const
    {
    return static_cast<TOutput>( vnl_math_rnd( static_cast<double>(A) ) );
    }
};
}

template <class TInputImage, class TOutputImage>
class ITK_EXPORT RoundAndCastImageFilter :
    public
UnaryFunctorImageFilter<TInputImage,TOutputImage,
                        Functor::RoundAndCast<
  typename TInputImage::PixelType,
  typename TOutputImage::PixelType> >
{
public:
  /** Standard class typedefs. */
  typedef RoundAndCastImageFilter               Self;
  typedef UnaryFunctorImageFilter<TInputImage,TOutputImage,
                        Functor::RoundAndCast<
                          typename TInputImage::PixelType,
                          typename TOutputImage::PixelType>
                                     >  Superclass;
  typedef SmartPointer<Self>            Pointer;
  typedef SmartPointer<const Self>      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(RoundAndCastImageFilter, UnaryFunctorImageFilter);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputCovertibleToDoubleCheck,
    (Concept::Convertible<typename TInputImage::PixelType, double>));
  itkConceptMacro(DoubleConvertibleToOutputCheck,
    (Concept::Convertible<double, typename TOutputImage::PixelType>));
  itkConceptMacro(InputConvertibleToOutputCheck,
    (Concept::Convertible<typename TInputImage::PixelType,
                          typename TOutputImage::PixelType>));
  /** End concept checking */
#endif

protected:
  RoundAndCastImageFilter() {}
  virtual ~RoundAndCastImageFilter() {}

  void GenerateData()
    {
    if( this->GetInPlace() && this->CanRunInPlace() )
      {
      // nothing to do, so avoid iterating over all the pixels
      // for nothing! Allocate the output, generate a fake progress and exit
      this->AllocateOutputs();
      ProgressReporter progress(this, 0, 1);
      return;
      }
    Superclass::GenerateData();
    }



private:
  RoundAndCastImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};


} // end namespace itk

#endif
