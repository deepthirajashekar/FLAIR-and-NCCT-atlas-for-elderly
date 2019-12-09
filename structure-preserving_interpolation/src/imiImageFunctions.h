/** \file imiImageFunctions.h
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
#ifndef __imiImageFunctions_h
#define __imiImageFunctions_h

// ITK includes:
#include "itkNumericTraits.h"
#include "itkInterpolateImageFunction.h"
// Project includes:
#include "imiMacro.h"

namespace imi
{
/** \class imi::imiImageFunctions
 *  \brief Header file for a class collecting general static functions.
 *
 *  <small> <!--Copyright Information: -->
 *  \b Initial \b Author: Jan Ehrhardt \n
 *  \b Copyright (C) 2010 Institute of Medical Informatics,
 *     University of Luebeck\n
 *  </small>
 ***************************************************************************/
class imiImageFunctions
{
private:
  /** The global debug level. */
  static int m_iGlobalDebugLevel;

public:
  /** \brief Set the global debug level.
   *
   *  The debug level is static, i.e. the same for all instances of imiObject.
   *  \param level The debug level. */
  static void SetGlobalDebugLevel( int level ) { m_iGlobalDebugLevel = level; }

  /** \brief Get the global debug level.
   *  \return The debug level.  */
  static int GetGlobalDebugLevel() { return m_iGlobalDebugLevel; }

  /** \brief Check if the global debug level is above a given level.
   *
   *  \param level The level to check.
   *  \return Boolean if the global debug level is >= level. */
  static bool CheckGlobalDebugLevel(int level)
    { return m_iGlobalDebugLevel >= level; }

  template<typename TImageType>
  static bool ReadImage(
      std::string filename,
      itk::SmartPointer<TImageType> &pImage);

  template<typename TImageType>
  static bool WriteImage(
      itk::SmartPointer<TImageType> pImage,
      std::string filename );

  template<typename TImageType1, typename TImageType2>
  static bool CheckImageParameters(
      itk::SmartPointer<TImageType1> inputImage1,
      itk::SmartPointer<TImageType2> inputImage2,
      bool checkSize = true,
      bool checkSpacing = true,
      bool checkOrigin = true );

  template<typename TImageType1, typename TImageType2>
  static bool CheckImageParameters(
      itk::SmartPointer<const TImageType1> inputImage1,
      itk::SmartPointer<const TImageType2> inputImage2,
      bool checkSize = true,
      bool checkSpacing = true,
      bool checkOrigin = true );

  template<typename TImageType1, typename TImageType2>
  static bool CheckImageParametersQuiet(
      itk::SmartPointer<TImageType1> inputImage1,
      itk::SmartPointer<TImageType2> inputImage2,
      bool checkSize = true,
      bool checkSpacing = true,
      bool checkOrigin = true );

  template<typename TImageType, class TDisplacementFieldType>
  static bool WarpImageND(
      itk::SmartPointer<TDisplacementFieldType> displField,
      itk::SmartPointer<TImageType> inputImage,
      itk::SmartPointer<TImageType> &outputImage,
      int InterpolationType,
      int WarpFilterType );

  template<typename TImageType, class TDisplacementFieldType>
  static bool WarpImageWithFactorND(
      itk::SmartPointer<TDisplacementFieldType> displField,
      itk::SmartPointer<TImageType> inputImage,
      itk::SmartPointer<TImageType> &outputImage,
      double factor,
      int InterpolationType = 1,
      int WarpFiltertype = 1 );

  template<class TImageType>
  static typename itk::InterpolateImageFunction<TImageType, double>::Pointer GetInterpolateImageFunctionByType(int InterpolationType);

  /** allocates an image with same size, spacing origin etc as the reference image */
  template<typename TImageType>
  static typename TImageType::Pointer AllocateImage(itk::SmartPointer<TImageType> referenceImage);

  /** allocates an image with different TYPE, but same size, spacing origin etc as the reference image */
  template<typename TInputImageType, typename TOutputImageType>
  static void AllocateImage( itk::SmartPointer<TInputImageType> referenceImage, itk::SmartPointer<TOutputImageType> &outputImage);

  /** allocates an image with different TYPE, but same size, spacing origin etc as the reference image */
  template<typename TInputImageType, typename TOutputImageType>
  static void AllocateVectorImage( itk::SmartPointer<TInputImageType> referenceImage, int vectorLength, itk::SmartPointer<TOutputImageType> &outputImage);

};

} // namespace

#include "imiImageFunctions.hxx"

#endif // __imiImageFunctions_h
