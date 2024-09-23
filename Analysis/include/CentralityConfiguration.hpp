// $HEADER$
//------------------------------------------------------------------------------------------------
//                           CentralityConfiguration flags
//------------------------------------------------------------------------------------------------
// CentralityConfiguration
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Flags and constant variables that determine centrality configuration
 **/
//------------------------------------------------------------------------------------------------

#ifndef CENTRALITY_CONFIGURATION
#define CENTRALITY_CONFIGURATION

#include <string>
#include <array>

#include "GlobalConfiguration.h"

#ifdef RUN7AUAU200

#ifdef CENTRALITY_CLASSES_MB4

struct CentralityContainer
{
	static constexpr int size = 5;
   
	std::array<int, size> minIndex = {0, 2, 4, 6, 0};
	std::array<int, size> maxIndex = {1, 3, 5, 9, 9};

	std::array<std::string, size> minName = {"0", "20", "40", "60", "0"};
	std::array<std::string, size> maxName = {"20", "40", "60", "93", "93"};

	std::array<std::string, size> name = {"0-20%", "20-40%", "40-60%", "60-93%", "MB"};
	std::array<std::string, size> nameWithoutPercent = {"0-20", "20-40", "40-60", "60-93", "MB"};
	std::array<std::string, size> nameLatex = {"0-20\\%", "20-40\\%", "40-60\\%", "60-93\\%", "MB"};

	std::array<Color_t, size> color = {kAzure+5, kTeal+5, kSpring+5, kOrange+5, kBlue-2};
	std::array<Style_t, size> markerStyle = {53, 54, 55, 59, 47};

	std::array<double, size> nColls = {784, 300.8, 94.2, 14.8, 257.8};
	std::array<double, size> nCollsUncertainty = {77.0, 31.4, 13.0, 4.0, 25.4};

	const int centralIndex = 0;
	const int peripheralIndex = 3;
};

#endif /* RUN7AUAU_MB4 */

#ifdef CENTRALITY_CLASSES_MB5

struct CentralityContainer
{
	static constexpr int size = 6;
   
	std::array<int, size> minIndex = {0, 1, 2, 4, 6, 0};
	std::array<int, size> maxIndex = {0, 1, 3, 5, 9, 9};

	std::array<std::string, size> cmin_name = {"0", "10", "20", "40", "60", "0"};
	std::array<std::string, size> cmax_name = {"10", "20", "40", "60", "93", "93"};

	std::array<std::string, size> cname_nop = {"0-10", "10-20", "20-40", "40-60", "60-93", "MB"};
	std::array<std::string, size> cname = {"0-10%", "10-20%", "20-40%", "40-60%", "60-93%", "MB"};
	std::array<std::string, size> cname_latex = {"0-10\\%", "10-20\\%", "20-40\\%", "40-60\\%", "60-93\\%", "MB"};

	std::array<Color_t, size> color = {kAzure, kGray, kGreen, kOrange, kRed, kViolet};
	std::array<Style_t, size> marker_style = {71, 74, 72, 73, 77, 75};

	std::array<double, size> ncolls = {900., 784, 300.8, 94.2, 14.8, 257.8};
	std::array<double, size> ncolls_uncertainty = {80., 77.0, 31.4, 13.0, 4.0, 25.4};

	const int centralIndex = 0;
	const int peripheralIndex = 4;
};

#endif /* RUN7AUAU_MB5 */

#endif /* RUN7AUAU200 */

#ifdef RUN14HEAU200

#ifdef CENTRALITY_CLASSES_MB4

struct CentralityContainer
{
	std::array<int, 5> cmin = {0, 0, 2, 4, 6};
	std::array<int, 5> cmax = {9, 1, 3, 5, 9};

	std::array<std::string, 5> cmin_name = {"00", "0", "20", "40", "60"};
	std::array<std::string, 5> cmax_name = {"88", "20", "40", "60", "88"};

	std::array<Color_t, 5> color = {kAzure, kGreen, kOrange, kRed, kViolet};
	std::array<Style_t, 5> marker_style = {53, 54, 55, 59};
};

#endif /* CENTRALITY_CLASSES_MB4 */

#endif /* RUN14HEAU200 */

#endif /* CENTRALITY_CONFIGURATION */
