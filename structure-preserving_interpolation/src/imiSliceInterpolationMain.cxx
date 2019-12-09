/** \file imiSliceInterpolationMain.cxx
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
// System includes:
#include <iostream>
#include <string>
// TCLAP command line parsing includes:
#include "tclap/CmdLine.h"
#include "tclap/VectorValueArg.h"
// ITK includes
#include "itkCastImageFilter.h"
#include "itkWarpImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkRoundAndCastImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkJoinSeriesImageFilter.h"
// ITK Variational Registration includes:
#include "itkVariationalRegistrationMultiResolutionFilter.h"
#include "itkVariationalSymmetricDiffeomorphicRegistrationFilter.h"
#include "itkVariationalRegistrationDemonsFunction.h"
#include "itkVariationalRegistrationDiffusionRegularizer.h"
#include "itkVariationalRegistrationStopCriterion.h"
#include "itkVariationalRegistrationLogger.h"
#include "itkExponentialDisplacementFieldImageFilter.h"
// IMI Project includes:
#include "imiRegistrationBasedInterpolationFilter.h"
#include "imiMacro.h"
#include "imiImageFunctions.h"

using namespace TCLAP;
using namespace imi;
/////////////////////////////////////////////////////////////////////////////////
//
// Define global variables for registration
//
// this is to avoid function calling with many parameters
/////////////////////////////////////////////////////////////////////////////////
bool bUseLinearInterpolation = false;
struct RegistrationParameters {
	double timestep = 1.0;
	double regulAlpha = 2.0;
	std::vector<unsigned int> numberOfIterationsPerLevel = {100, 100, 50};
	unsigned int numberOfLevels = 3;
	bool useImageSpacing = true;
	double stopCriterionSlope = 0.005;

	void SetClampedNumberOfIterationsPerLevel(unsigned int levels, std::vector<unsigned int> its)
	{
    assert(levels > 0);
    assert(its.size() > 0);
    numberOfLevels = levels;
	  numberOfIterationsPerLevel.resize(numberOfLevels);
	  for(size_t i=0; i<numberOfLevels; i++)
	    numberOfIterationsPerLevel[i] = its[std::min(i, its.size()-1)];
	}

	void Print() {
		imiINFO("Registration Parameters:");
		imiINFO("             time-step: "<<timestep);
		imiINFO("                 alpha: "<<regulAlpha);
		imiINFO("      number of levels: "<<numberOfLevels);
		imiINFO("  number of iterations: ");
		for (const auto it: numberOfIterationsPerLevel)
		  std::cout << it << ' ';
		imiINFO("     use image spacing: "<<((useImageSpacing) ? "true" : "false"));
		imiINFO("  stop criterion slope: "<<stopCriterionSlope);
	}
};
/////////////////////////////////////////////////////////////////////////////////
//
// Define Template functions to call
//
/////////////////////////////////////////////////////////////////////////////////
template<class TImage3DType, class TLabelImage3DType>
bool ProcessSliceInterpolation(typename TImage3DType::Pointer inputImage,
		typename TLabelImage3DType::Pointer labelImage,
		int numberSlicesToInterpolate,
		unsigned int imageInterpolationType,
		unsigned int labelInterpolationType,
		const RegistrationParameters regParams,
		typename TImage3DType::Pointer &outputImage,
		typename TLabelImage3DType::Pointer &outputLabelImage);

template<class TFixedImageType, class TMovingImageType,
		class TDisplacementFieldType>
bool ComputeVariationalDiffeomorphicRegistration(
		typename TFixedImageType::Pointer fixedImage,
		typename TMovingImageType::Pointer movingImage,
		const RegistrationParameters regParams,
		typename TDisplacementFieldType::Pointer &outputVelocityField);

template<class TDisplacementFieldType>
void PrintFieldStatistic(itk::SmartPointer<TDisplacementFieldType> inputField,
		const char *fieldName);

/////////////////////////////////////////////////////////////////////////////////
//
// Start Main Program
//
/////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
	std::string copyright = "Copyright (C) 2010, Jan Ehrhardt, Institute of Medical Informatics,";
	copyright += " University of Luebeck, Germany. All rights reserved.\n";
	copyright += "Please cite the following paper:\n";
	copyright += "  Ehrhardt, J., Säring, D., & Handels, H. (2007). Structure-preserving";
	copyright += " interpolation of temporal and spatial image sequences using an optical";
	copyright += " flow-based method. Methods of information in medicine, 46(03), 300-307.\n";

	imiINFO("==========================================");
	imiINFO("=====      imiSliceInterpolation     =====");
	imiINFO("==========================================");
	imiINFO(copyright);
	imiINFO("READING parameters...\n");

	std::string programDescription;
	programDescription =
			"This program implements a slice-wise registration-based interpolation of a 3D input image as described in:\n";
	programDescription +=
			"  Ehrhardt, J., Säring, D., & Handels, H. (2007). Structure-preserving interpolation of temporal and spatial image sequences using an optical flow-based method. Methods of information in medicine, 46(03), 300-307.\n";
	programDescription +=
			"Instead of the optical-flow based registration method, we use a symmetric diffeomorph variational registration with Thirion-like forces and diffusive regularization.\n";
	programDescription += "See \n";
	programDescription +=
			" Werner, R., Schmidt-Richberg, A., Handels, H., & Ehrhardt, J. (2014). Estimation of lung motion fields in 4D CT data by variational non-linear intensity-based registration: A comparison and evaluation study. Physics in Medicine & Biology, 59(15), 4247.\n";
	programDescription += "or \n";
	programDescription +=
			" Schmidt-Richberg, A., Werner, R., Handels, H., & Ehrhardt, J. (2014). A flexible variational registration framework. Insight Journal.\n";

	// Variable to hold registration parameters, initialized with default values
	RegistrationParameters regParams;

	std::vector<std::string> allowedInterpolationTypes;
	allowedInterpolationTypes.push_back("nearest");
	allowedInterpolationTypes.push_back("linear");
	allowedInterpolationTypes.push_back("cubic");
	ValuesConstraint<std::string> allowedInterpolationValues(
			allowedInterpolationTypes);

	std::vector<std::string> allowedLabelInterpolationTypes;
	allowedLabelInterpolationTypes.push_back("nearest");
	allowedLabelInterpolationTypes.push_back("linear");
	allowedLabelInterpolationTypes.push_back("cubic");
	allowedLabelInterpolationTypes.push_back("gaussian");
	ValuesConstraint<std::string> allowedLabelInterpolationValues(
			allowedLabelInterpolationTypes);

	CmdLine cmdParser(programDescription, ' ', "0.9");
	ValueArg<int> debugLevelArg("x", "debug-level", "Debug-level", false, 3,
			"[0-9]", cmdParser);
	SwitchArg linearInterpolationSwitch("", "linear",
			"Use only linear interpolation between slices (without inter-slice registration)",
			cmdParser);
	ValueArg<std::string> labelInterpolatorArg("", "label-interpolate",
			"Interpolation type used for warping the label slices. (default: 'nearest').",
			false, "nearest", &allowedLabelInterpolationValues, cmdParser);
	ValueArg<std::string> interpolatorArg("", "interpolate",
			"Interpolation type used for warping the image slices. (default: 'linear').",
			false, "linear", &allowedInterpolationValues, cmdParser);
	ValueArg<double> regStopConditionSlopeArg("g", "reg-slope",
			"Registration: Fitted line slope for stop criterion (default=" + std::to_string(regParams.stopCriterionSlope) + ").",
			false, regParams.stopCriterionSlope, "slope", cmdParser);
	ValueArg<double> regTauArg("t", "reg-tau",
			"Registration: Time step tau (default=" + std::to_string(regParams.timestep) + ").",
			false, regParams.timestep, "tau", cmdParser);
	ValueArg<double> regAlphaArg("a", "reg-alpha",
			"Registration: Regularization weight alpha (default=" + std::to_string(regParams.regulAlpha) + ").",
			false, regParams.regulAlpha, "alpha", cmdParser);
	ValueArg<unsigned int> regNumberOfLevelsArg("l", "reg-levels",
			"Registration: Number of multi-resolution levels (default=" + std::to_string(regParams.numberOfLevels) + ").",
			false, regParams.numberOfLevels, "levels", cmdParser);
	VectorValueArg<unsigned int> regNumberOfIterationsArg("i", "reg-iters",
			"Registration: Number of iterations (default=100).",
			false, 100, "iterations", cmdParser);
	ValueArg<unsigned int> numberOfSlicesArg("n", "num-slices",
			"Number of slices to interpolate between original slices. '0' will result in a nearly isotropic output image (default=0).",
			false, 0, "num-slices", cmdParser);
	ValueArg<std::string> outputLabelImageFilenameArg("Q", "output-label",
			"Output label image filename", false, "", "filename", cmdParser);
	ValueArg<std::string> outputImageFilenameArg("O", "output",
			"Output image filename", true, "", "filename", cmdParser);
	ValueArg<std::string> inputLabelImageFilenameArg("L", "label",
			"Filename of the 3D image (unsigned char) with labels to interpolate too.",
			false, "", "filename", cmdParser);
	ValueArg<std::string> inputImageFilenameArg("I", "input",
			"Filename of the 3D input image to interpolate slices.", true, "",
			"filename", cmdParser);

	// Parse the args.
	cmdParser.parse(argc, argv);

	std::string inputImageFilename = inputImageFilenameArg.getValue();
	std::string inputLabelImageFilename = inputLabelImageFilenameArg.getValue();
	std::string outputImageFilename = outputImageFilenameArg.getValue();
	std::string outputLabelImageFilename =
			outputLabelImageFilenameArg.getValue();
	unsigned int numberSlicesToInterpolate = numberOfSlicesArg.getValue();
	// setup registration parameters
	regParams.timestep = regTauArg.getValue();
	regParams.regulAlpha = regAlphaArg.getValue();
  regParams.numberOfLevels = regNumberOfLevelsArg.getValue(); // set number of levels first
  if(regNumberOfIterationsArg.isSet())
  {
    regParams.SetClampedNumberOfIterationsPerLevel(regParams.numberOfLevels, regNumberOfIterationsArg.getValue());
  }
	regParams.stopCriterionSlope = regStopConditionSlopeArg.getValue();

	unsigned int imageInterpolationType = 0;
	for (imageInterpolationType = 0;
			imageInterpolationType < allowedInterpolationTypes.size();
			imageInterpolationType++)
		if (allowedInterpolationTypes[imageInterpolationType]
				== interpolatorArg.getValue())
			break;
	unsigned int labelInterpolationType = 0;
	for (labelInterpolationType = 0;
			labelInterpolationType < allowedLabelInterpolationTypes.size();
			labelInterpolationType++)
		if (allowedLabelInterpolationTypes[labelInterpolationType]
				== labelInterpolatorArg.getValue())
			break;

	bUseLinearInterpolation = linearInterpolationSwitch.isSet();
	imiImageFunctions::SetGlobalDebugLevel(debugLevelArg.getValue());

	if(imiImageFunctions::CheckGlobalDebugLevel(3))
	{
		imiINFO("                   Input image: "<<inputImageFilename);
		imiINFO("                  Output image: "<<outputImageFilename);
		if(inputLabelImageFilenameArg.isSet())
			imiINFO("             Input label image: "<<inputLabelImageFilename);
		if(outputLabelImageFilenameArg.isSet())
			imiINFO("            Output label image: "<<outputLabelImageFilename);
		imiINFO("number of interpolation slices: "<<numberSlicesToInterpolate);
		imiINFO("            image interpolator: "<<allowedInterpolationTypes[imageInterpolationType]);
		if(inputLabelImageFilenameArg.isSet())
			imiINFO("      label image interpolator: "<<allowedLabelInterpolationTypes[labelInterpolationType]);
		imiINFO("      use linear interpolation: "<<((bUseLinearInterpolation)? "true" : "false"));
		regParams.Print();
	}
	//////////////////////////////////////////////
	//
	// Type Definitions
	//
	//////////////////////////////////////////////
	typedef itk::Image<short, 3> InputImageType;
	typedef itk::Image<unsigned char, 3> LabelImageType;

	//////////////////////////////////////////////
	//
	// Load input images
	//
	//////////////////////////////////////////////
	imiINFO("Loading input image...");
	InputImageType::Pointer inputImage = nullptr;

	if (!imiImageFunctions::ReadImage(inputImageFilename, inputImage)) {
		imiERROR("Cannot load input image!");
		return EXIT_FAILURE;
	}

	// load label image
	LabelImageType::Pointer labelImage = nullptr;
	if (inputLabelImageFilenameArg.isSet()) {
		if (!outputLabelImageFilenameArg.isSet()) {
			imiERROR("No output label image filename given!");
			return EXIT_FAILURE;
		}

		imiINFO("Loading label image...");
		if (!imiImageFunctions::ReadImage(inputLabelImageFilename,
				labelImage)) {
			imiERROR("Cannot load label image!");
			return EXIT_FAILURE;
		}
	}

	/////////////////////////////////////////////////////////////////////////////////
	// Compute number of slices to interpolate, if not given
	//
	if (numberSlicesToInterpolate == 0) {
		imiINFO("Compute numbers of slices to interpolate, generate (nearly) isotropic image spacing.");

		InputImageType::SpacingType inputImageSpacing =
				inputImage->GetSpacing();
		numberSlicesToInterpolate = round(
				inputImageSpacing[2] / inputImageSpacing[0]);
	}
	//////////////////////////////////////////////
	//
	// Call Interpolation dependent on image type
	//
	//////////////////////////////////////////////
	InputImageType::Pointer outputImage;
	LabelImageType::Pointer outputLabelImage;
	bool success = ProcessSliceInterpolation<InputImageType, LabelImageType>(
			inputImage, labelImage, numberSlicesToInterpolate, imageInterpolationType, labelInterpolationType, regParams, outputImage,
			outputLabelImage);
	if (!success) {
		imiERROR("Failed to interpolate image!");
		return EXIT_FAILURE;
	}

	imiINFO("Save interpolated image.");
	if (!imiImageFunctions::WriteImage(outputImage,
			outputImageFilenameArg.getValue())) {
		imiERROR("Cannot save image!");
		return EXIT_FAILURE;
	}

	if (outputLabelImageFilenameArg.isSet() && outputLabelImage.IsNotNull()) {
		imiINFO("Save interpolated Label image.");
		if (!imiImageFunctions::WriteImage(outputLabelImage,
				outputLabelImageFilenameArg.getValue())) {
			imiERROR("Cannot save label image!");
			return EXIT_FAILURE;
		}
	}

	imiINFO("Finished.");
	imiINFO("==========================================\n");
	return EXIT_SUCCESS;
}

template<class TImage3DType, class TLabelImage3DType>
bool ProcessSliceInterpolation(typename TImage3DType::Pointer inputImage,
		typename TLabelImage3DType::Pointer labelImage,
		int numberSlicesToInterpolate,
		unsigned int imageInterpolationType,
		unsigned int labelInterpolationType,
		const RegistrationParameters regParams,
		typename TImage3DType::Pointer &outputImage,
		typename TLabelImage3DType::Pointer &outputLabelImage) {
	if (numberSlicesToInterpolate <= 0)
		return false;

	imiINFO("Start interpolating "<<numberSlicesToInterpolate<<" slices between neighbouring slices...");
	//
	// Some Typedefs
	//
	// image-type for a 2D slice (we use float internal to account for rounding errors)
	typedef itk::Image<float, 2> SliceImageType;
	// image-type for the rejoint 3D image (we use float internal to account for rounding errors)
	typedef itk::Image<float, 3> VolumeImageType;
	// image-type for 2D displacement fields between two slices
	typedef itk::Image<itk::Vector<float, 2>, 2> DisplacementField2DType;
	typedef DisplacementField2DType::Pointer DisplacementField2DTypePointer;

	/////////////////////////////////////////////////////////////////////////////////
	// Get size and spacing of input image
	//
	typename TImage3DType::SizeType inputImageSize =
			inputImage->GetLargestPossibleRegion().GetSize();
	typename TImage3DType::IndexType inputImageStart =
			inputImage->GetLargestPossibleRegion().GetIndex();
	typename TImage3DType::SpacingType inputImageSpacing =
			inputImage->GetSpacing();
	typename TImage3DType::PointType inputImageOrigin = inputImage->GetOrigin();

	const unsigned int startSlice = inputImageStart[2];
	const unsigned int endSlice = startSlice + inputImageSize[2] - 1;

	/////////////////////////////////////////////////////////////////////////////////
	// Compute interpolationStepSize and z-Spacing from number of slices to interpolate
	//
	const double interpolationValueStepSize = 1.0
			/ static_cast<double>(numberSlicesToInterpolate);
	const double outputZSpacing = inputImageSpacing[2]
			/ static_cast<double>(numberSlicesToInterpolate);

	typename TImage3DType::SizeType currentSliceSize = inputImageSize;
	currentSliceSize[2] = 0; // set z-dimension to zero to extract only one slice
	typename TImage3DType::IndexType currentSliceStart = inputImageStart;
	//set extract region for the first slice
	typename TImage3DType::RegionType currentSliceRegion;

	//Define JoinSeriesFilter to build a 3D image stack for storing all single slices

	typedef itk::JoinSeriesImageFilter<SliceImageType, VolumeImageType> JoinSeriesImageFilter;
	JoinSeriesImageFilter::Pointer pJoinSeriesImageFilter =
			JoinSeriesImageFilter::New();
	pJoinSeriesImageFilter->SetOrigin(inputImageOrigin[2]);
	pJoinSeriesImageFilter->SetSpacing(outputZSpacing);

	//Define JoinSeriesFilter to build a 3D image stack for storing all single label slices
	JoinSeriesImageFilter::Pointer pJoinSeriesLabelImageFilter =
			JoinSeriesImageFilter::New();
	pJoinSeriesLabelImageFilter->SetOrigin(inputImageOrigin[2]);
	pJoinSeriesLabelImageFilter->SetSpacing(outputZSpacing);

	//
	// Go through all slices of the 3D Image (z-direction) and extract
	// two adjacent slices
	//
	unsigned int sliceCounter = 0;
	unsigned int labelSliceCounter = 0;
	for (unsigned int currentSlice = startSlice; currentSlice < endSlice;
			currentSlice++) {
		/////////////////////////////////////////////////////////////////////////////////
		// Extract two adjacent slices
		//
		imiDEBUGINFO(3,
				"  extract slice "<<currentSlice<<" and "<<currentSlice+1<<" (of "<<inputImageSize[2]<<" slices)...");

		//
		// Setup the correct region for the current slice
		currentSliceStart[2] = currentSlice;
		currentSliceRegion.SetSize(currentSliceSize);
		currentSliceRegion.SetIndex(currentSliceStart);

		//
		// Use itkExtractImageFilter to extract current slice
		typedef itk::ExtractImageFilter<TImage3DType, SliceImageType> ExtractSliceFilterType;
		typedef itk::ExtractImageFilter<TLabelImage3DType, SliceImageType> ExtractLabelSliceFilterType;

		typename ExtractSliceFilterType::Pointer pExtractCurrentSliceFilter =
				ExtractSliceFilterType::New();
		pExtractCurrentSliceFilter->SetInput(inputImage);
		pExtractCurrentSliceFilter->SetExtractionRegion(currentSliceRegion);
		pExtractCurrentSliceFilter->Update();

		SliceImageType::Pointer currentSliceImage =
				pExtractCurrentSliceFilter->GetOutput();

		//
		// Extract the current slice of the label image, if given
		SliceImageType::Pointer currentLabelSlice = nullptr;
		if (labelImage.IsNotNull()) {
			typename ExtractLabelSliceFilterType::Pointer pExtractCurrentLabelSliceFilter =
					ExtractLabelSliceFilterType::New();
			pExtractCurrentLabelSliceFilter->SetInput(labelImage);
			pExtractCurrentLabelSliceFilter->SetExtractionRegion(
					currentSliceRegion);
			pExtractCurrentLabelSliceFilter->Update();

			currentLabelSlice = pExtractCurrentLabelSliceFilter->GetOutput();
		}

		//
		// Setup the correct region for the next slice
		currentSliceStart[2] = currentSlice + 1;
		currentSliceRegion.SetSize(currentSliceSize);
		currentSliceRegion.SetIndex(currentSliceStart);

		//
		// Use itkExtractImageFilter to extract next slice
		typename ExtractSliceFilterType::Pointer pExtractNextSliceFilter =
				ExtractSliceFilterType::New();
		pExtractNextSliceFilter->SetInput(inputImage);
		pExtractNextSliceFilter->SetExtractionRegion(currentSliceRegion);
		pExtractNextSliceFilter->Update();

		SliceImageType::Pointer nextSliceImage =
				pExtractNextSliceFilter->GetOutput();

		//
		// Extract the next slice of the label image, if given
		SliceImageType::Pointer nextLabelSlice = NULL;
		if (labelImage.IsNotNull()) {
			typename ExtractLabelSliceFilterType::Pointer pExtractNextLabelSliceFilter =
					ExtractLabelSliceFilterType::New();
			pExtractNextLabelSliceFilter->SetInput(labelImage);
			pExtractNextLabelSliceFilter->SetExtractionRegion(
					currentSliceRegion);
			pExtractNextLabelSliceFilter->Update();

			nextLabelSlice = pExtractNextLabelSliceFilter->GetOutput();
		}

		/////////////////////////////////////////////////////////////////////////////////
		// Compute a registration between the two slice images
		//
		imiDEBUGINFO(5, "  register these two neighbouring slices ...");

		DisplacementField2DTypePointer velocityField;
		//typedef itk::Image<itk::Vector<float, 3>, 3> VeloFieldType;
		//VeloFieldType::Pointer velocityField;

		//
		// Perform the registration, get velocity field in return
		if (!ComputeVariationalDiffeomorphicRegistration<SliceImageType, SliceImageType, DisplacementField2DType>(
				currentSliceImage, nextSliceImage, regParams, velocityField)) {
			imiERROR(
					"Diffeomorphic registration for slice "<<currentSlice<<" failed!");
			return false;
		}

		/////////////////////////////////////////////////////////////////////////////////
		// Get forward and backward displacement field from the velocity field
		//
		imiDEBUGINFO(5, "  compute forward and inverse displacement field ...");

		typedef itk::ExponentialDisplacementFieldImageFilter<
				DisplacementField2DType, DisplacementField2DType> FieldExponentiatorType;
		FieldExponentiatorType::Pointer exponentiator =
				FieldExponentiatorType::New();

		exponentiator->SetInput(velocityField);
		exponentiator->AutomaticNumberOfIterationsOff();
		exponentiator->SetMaximumNumberOfIterations(4);
		exponentiator->Update();

		DisplacementField2DTypePointer outputForwardDisplacementField =
				exponentiator->GetOutput();

		FieldExponentiatorType::Pointer inverseExponentiator =
				FieldExponentiatorType::New();

		inverseExponentiator->SetInput(velocityField);
		inverseExponentiator->AutomaticNumberOfIterationsOff();
		inverseExponentiator->SetMaximumNumberOfIterations(4);
		inverseExponentiator->ComputeInverseOn();
		inverseExponentiator->Update();

		DisplacementField2DTypePointer outputInverseDisplacementField =
				inverseExponentiator->GetOutput();

		if (imiImageFunctions::CheckGlobalDebugLevel(7)) {
			PrintFieldStatistic(outputForwardDisplacementField,
					"forward-Displacement");
			PrintFieldStatistic(outputInverseDisplacementField,
					"inverse-Displacement");
		}
		/////////////////////////////////////////////////////////////////////////////////
		// Perform the interpolation between the slices
		//
		imiDEBUGINFO(5, "  perform interpolation between the slices...");

		typedef imiRegistrationBasedInterpolationFilter<
				SliceImageType::PixelType, 2> InterpolationFilterType;

		InterpolationFilterType* interpolationFilter =
				new InterpolationFilterType();
		interpolationFilter->SetInputReferenceImage(currentSliceImage);
		interpolationFilter->SetInputTargetImage(nextSliceImage);
		interpolationFilter->SetDeformationField(
				outputForwardDisplacementField);
		interpolationFilter->SetInverseDeformationField(
				outputInverseDisplacementField);
		switch (imageInterpolationType) {
		case 0:
			interpolationFilter->SetInterpolationTypeToNearestNeighbour();
			break;
		case 1:
			interpolationFilter->SetInterpolationTypeToLinear();
			break;
		case 2:
			interpolationFilter->SetInterpolationTypeToCubicBSpline();
			break;
		default:
			interpolationFilter->SetInterpolationTypeToLinear();
			break;
		}

		//
		// generate interpolation filter for label images, if given
		InterpolationFilterType* interpolationLabelFilter = nullptr;
		if (labelImage.IsNotNull()) {
			interpolationLabelFilter = new InterpolationFilterType();
			interpolationLabelFilter->SetInputReferenceImage(currentLabelSlice);
			interpolationLabelFilter->SetInputTargetImage(nextLabelSlice);
			interpolationLabelFilter->SetDeformationField(
					outputForwardDisplacementField);
			interpolationLabelFilter->SetInverseDeformationField(
					outputInverseDisplacementField);
			switch (labelInterpolationType) {
			case 0:
				interpolationLabelFilter->SetInterpolationTypeToNearestNeighbour();
				break;
			case 1:
				interpolationLabelFilter->SetInterpolationTypeToLinear();
				break;
			case 2:
				interpolationLabelFilter->SetInterpolationTypeToCubicBSpline();
				break;
			case 3:
				interpolationLabelFilter->SetInterpolationTypeToLabelGaussian();
				break;
			default:
				interpolationLabelFilter->SetInterpolationTypeToNearestNeighbour();
				break;
			}
		}


		////////////////////////////////////////////////////////
		// Start interpolating slices in between
		//
		for (unsigned int interpolationStep = 0;
				interpolationStep < numberSlicesToInterpolate;
				interpolationStep++) {
			// compute the interpolation value between 0 and 1 (is always LESS than one)
			const double interpolationValue =
					static_cast<double>(interpolationStep)
							* interpolationValueStepSize;

			interpolationFilter->SetInterpolationValue(interpolationValue);

			imiDEBUGINFO(7,
					"    Execute interpolation (value="<<interpolationValue<<")");

			// Execute the interpolator for this slice
			if (!interpolationFilter->Execute()) {
				imiERROR("Execution of interpolation filter failed!");
				return false;
			}

			SliceImageType::Pointer interpolatedImageslice =
					interpolationFilter->GetOutputImage();

			//
			// Add the interpolated slice to the image stack
			//
			imiDEBUGINFO(7,
					"    Insert interpolated slice into the output image...");

			pJoinSeriesImageFilter->SetInput(sliceCounter,
					interpolatedImageslice);
			pJoinSeriesImageFilter->Update();
			sliceCounter++;

			if (labelImage.IsNotNull()) {
				imiDEBUGINFO(7,
						"    Execute interpolation of label image (value="<<interpolationValue<<")");

				interpolationLabelFilter->SetInterpolationValue(
						interpolationValue);

				// Execute the interpolator for this slice
				if (!interpolationLabelFilter->Execute()) {
					imiERROR("Execution of interpolation label filter failed!");
					return false;
				}

				SliceImageType::Pointer interpolatedLabelSlice =
						interpolationLabelFilter->GetOutputImage();

				//
				// Add the interpolated slice to the image stack
				//
				imiDEBUGINFO(7,
						"    Insert interpolated label slice into the output label image...");

				pJoinSeriesLabelImageFilter->SetInput(labelSliceCounter,
						interpolatedLabelSlice);
				pJoinSeriesLabelImageFilter->Update();
				labelSliceCounter++;
			}
		}

		/////////////////////////////////////////////////////////////////////////////////
		// We have to add the last slice of the original input
		//
		if (currentSlice == endSlice - 1) {
			imiDEBUGINFO(7, "  Add last slice to the output image...");

			pJoinSeriesImageFilter->SetInput(sliceCounter, nextSliceImage);
			pJoinSeriesImageFilter->Update();
			sliceCounter++;

			if (labelImage.IsNotNull()) {
				pJoinSeriesLabelImageFilter->SetInput(labelSliceCounter,
						nextLabelSlice);
				pJoinSeriesLabelImageFilter->Update();
				labelSliceCounter++;
			}
		}

		if (interpolationFilter) {
			delete interpolationFilter;
			interpolationFilter = nullptr;
		}
		if (interpolationLabelFilter) {
			delete interpolationLabelFilter;
			interpolationLabelFilter = nullptr;
		}
	}

	/////////////////////////////////////////////////////////////////////////////////
	// Write the interpolated output image to file
	//
	typename TImage3DType::Pointer interpolatedImageVolume = NULL;

	imiDEBUGINFO(5, "  Converting interpolated image back ...");
	//
	// Cast image back to original image type (loat to short..)
	//
	typedef itk::RoundAndCastImageFilter<VolumeImageType, TImage3DType> RoundImageFilterType;
	typename RoundImageFilterType::Pointer roundImageFilter =
			RoundImageFilterType::New();
	roundImageFilter->SetInput(pJoinSeriesImageFilter->GetOutput());

	try {
		roundImageFilter->Update();
		outputImage = roundImageFilter->GetOutput();
	} catch (itk::ExceptionObject & excep) {
		imiERROR(
				" Rounded casting failed with exception:" << excep << std::endl);
		return false;
	}

	if (labelImage.IsNotNull()) {
		imiDEBUGINFO(5, "  Converting label image back ...");

		//
		// Round shift-scaled image back to original label type
		//
		typedef itk::RoundAndCastImageFilter<VolumeImageType, TLabelImage3DType> RoundImageFilterType;
		typename RoundImageFilterType::Pointer roundFilter =
				RoundImageFilterType::New();
		roundFilter->SetInput(pJoinSeriesLabelImageFilter->GetOutput());

		try {
			roundFilter->Update();
			outputLabelImage = roundFilter->GetOutput();
		} catch (itk::ExceptionObject & excep) {
			imiERROR(
					" Rounded casting failed with exception:" << excep << std::endl);
			return false;
		}
	}

	return true;

}

/////////////////////////////////////////////////////////////////////////////////
//
// ComputeVariationalDiffeomorphicRegistration
//
/////////////////////////////////////////////////////////////////////////////////
template<class TFixedImageType, class TMovingImageType,
		class TDisplacementFieldType>
bool ComputeVariationalDiffeomorphicRegistration(
		typename TFixedImageType::Pointer fixedImage,
		typename TMovingImageType::Pointer movingImage,
		const RegistrationParameters regParams,
		typename TDisplacementFieldType::Pointer &outputVelocityField) {
	//////////////////////////////////////////////
	//
	// If linear interpolation is wanted, create zero-vectro field
	//
	//////////////////////////////////////////////
	if (bUseLinearInterpolation) {
		imiDEBUGINFO(5,
				"  Linear Interpolation between slices corresponds to zero displacement field...");

		outputVelocityField = TDisplacementFieldType::New();
		outputVelocityField->CopyInformation(fixedImage);
		outputVelocityField->SetRegions(fixedImage->GetLargestPossibleRegion());
		outputVelocityField->Allocate();

		typename TDisplacementFieldType::PixelType veloVector;
		veloVector.Fill(0);

		outputVelocityField->FillBuffer(veloVector);

		return true;
	}
	//
	// Setup registration function
	//
	typedef itk::VariationalRegistrationDemonsFunction<TFixedImageType,
			TMovingImageType, TDisplacementFieldType> DemonsFunctionType;
	typename DemonsFunctionType::Pointer demonsFunction =
			DemonsFunctionType::New();
	demonsFunction->SetGradientTypeToSymmetric();
	demonsFunction->SetTimeStep(regParams.timestep);
	//
	// Setup regularizer
	//
	typedef itk::VariationalRegistrationDiffusionRegularizer<
			TDisplacementFieldType> DiffusionRegularizerType;
	typename DiffusionRegularizerType::Pointer diffRegularizer =
			DiffusionRegularizerType::New();
	diffRegularizer->SetAlpha(regParams.regulAlpha);
	diffRegularizer->InPlaceOff();
	diffRegularizer->SetUseImageSpacing(regParams.useImageSpacing);
	//
	// Setup registration filter
	//
	typedef itk::VariationalSymmetricDiffeomorphicRegistrationFilter<
			TFixedImageType, TMovingImageType, TDisplacementFieldType> SymmetricDiffeomorphicRegistrationFilterType;
	typename SymmetricDiffeomorphicRegistrationFilterType::Pointer regFilter;
	regFilter = SymmetricDiffeomorphicRegistrationFilterType::New();
	regFilter->SetRegularizer(diffRegularizer);
	regFilter->SetDifferenceFunction(demonsFunction);
	//
	// Setup multi-resolution filter
	//
	unsigned int its[regParams.numberOfLevels];
	for (int level = 0; level < regParams.numberOfLevels; ++level) {
		its[level] = regParams.numberOfIterationsPerLevel[level];
	}

	typedef itk::VariationalRegistrationMultiResolutionFilter<TFixedImageType,
			TMovingImageType, TDisplacementFieldType> MRRegistrationFilterType;
	typename MRRegistrationFilterType::Pointer mrRegFilter =
			MRRegistrationFilterType::New();
	mrRegFilter->SetRegistrationFilter(regFilter);
	mrRegFilter->SetMovingImage(movingImage);
	mrRegFilter->SetFixedImage(fixedImage);
	//mrRegFilter->SetMaskImage( maskImage );
	mrRegFilter->SetNumberOfLevels(regParams.numberOfLevels);
	mrRegFilter->SetNumberOfIterations(its);
	//mrRegFilter->SetInitialDisplacementField( initialDisplacementField );
	//
	// Setup stop criterion
	//
	typedef itk::VariationalRegistrationStopCriterion<
			SymmetricDiffeomorphicRegistrationFilterType,
			MRRegistrationFilterType> StopCriterionType;
	typename StopCriterionType::Pointer stopCriterion =
			StopCriterionType::New();
	stopCriterion->SetRegressionLineSlopeThreshold(
			regParams.stopCriterionSlope);
	stopCriterion->PerformLineFittingMaxDistanceCheckOn();
	stopCriterion->SetMultiResolutionPolicyToSimpleGraduated();

	regFilter->AddObserver(itk::IterationEvent(), stopCriterion);
	mrRegFilter->AddObserver(itk::IterationEvent(), stopCriterion);
	mrRegFilter->AddObserver(itk::InitializeEvent(), stopCriterion);
	//
	// Setup logger
	//
	typedef itk::VariationalRegistrationLogger<
			SymmetricDiffeomorphicRegistrationFilterType,
			MRRegistrationFilterType> LoggerType;
	typename LoggerType::Pointer logger = LoggerType::New();

	if (imiImageFunctions::CheckGlobalDebugLevel(7)) {
		regFilter->AddObserver(itk::IterationEvent(), logger);
	}
	if (imiImageFunctions::CheckGlobalDebugLevel(3)) {
		mrRegFilter->AddObserver(itk::IterationEvent(), logger);
	}

	//
	// Execute registration
	//
	imiDEBUGINFO(5, "  Starting variational registration...");
	try {
		mrRegFilter->Update();
	} catch (itk::ExceptionObject& err) {
		imiERROR("Registration of slices failed! Error:" << err);
		return false;
	}

	imiDEBUGINFO(5, "  Registration execution finished.");

	outputVelocityField = mrRegFilter->GetOutput();

	return true;
}

template<class TDisplacementFieldType>
void PrintFieldStatistic(itk::SmartPointer<TDisplacementFieldType> inputField,
		const char *fieldName) {
	typedef itk::ImageRegionConstIterator<TDisplacementFieldType> ConstFieldIteratorType;

	ConstFieldIteratorType inputIt(inputField,
			inputField->GetLargestPossibleRegion());

	typename TDisplacementFieldType::PixelType pixelValue;
	double min = 100000;
	double max = -100000;
	double mean = 0;
	double meansquare = 0;
	unsigned long numzeropixels = 0;

	unsigned long numpixels = 0;

	double magnitude = 0.0;
	double magnitude2 = 0.0;

	for (inputIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt) {
		pixelValue = inputIt.Get();

		magnitude2 = 0.0;
		for (unsigned int i = 0; i < TDisplacementFieldType::ImageDimension;
				i++) {
			magnitude2 += pixelValue[i] * pixelValue[i];
		}

		magnitude = sqrt(magnitude2);
		meansquare += magnitude2;
		mean += magnitude;

		if (magnitude < min)
			min = magnitude;
		if (magnitude > max)
			max = magnitude;

		numpixels++;
		if (magnitude <= 1e-10)
			numzeropixels++;
	}
	mean /= numpixels;
	meansquare /= numpixels;

	imiINFO(
			"    " <<fieldName<< ": " <<min<<", "<<max<<", "<<mean<<", "<<meansquare<<", "<<numzeropixels);
}
