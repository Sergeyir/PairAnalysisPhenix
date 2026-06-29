/** 
 *  @file   DeadMapSys.hpp 
 *  @brief  Contains declarations of functions and variables that are used for estimation of systematic uncertainties of heatmaps (real vs simulated) 
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysis).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef DEAD_MAP_SYS_HPP
#define DEAD_MAP_SYS_HPP

#include <thread>

#include "TError.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TStyle.h"
#include "TLegend.h"

#include "StrTools.hpp"
#include "IOTools.hpp"
#include "MathTools.hpp"
#include "Table.hpp"

#include "TCanvasTools.hpp"

#include "PBar.hpp"

#include "InputYAMLReader.hpp"
#include "DeadMapCutter.hpp"

/*! @namespace DeadMapSys
 * @brief Contains all functions and variables for DeadMapSyscpp
 */
namespace DeadMapSys
{
   /// Contents of input .yaml file for run configuration
   InputYAMLReader inputYAMLMain;
   /// Name of run (e.g. Run14HeAu200 or Run7AuAu200)
   std::string runName;
   /// Output directory for pictures for deadmaps
   std::string outputDirDM;
   /// Output directory for pictures for systematics
   std::string outputDirSys;
   /// Output directory for parameters
   std::string outputDirParameters;
   /// Real data input file
   TFile *inputRealDataFile;
   /// Simulated data input file
   TFile *inputSimDataFile;
   /// Table in which systematics will be written
   CppTools::Table table{4};
   /// Cutter for bad/dead areas of heatmaps
   DeadMapCutter dmCutter;
   /*! @brief Returns the uncertainty that is originated from the difference between 2 2D 
    * histograms measured by taking numberOfProjectionDivisions regions for normalization. 
    * Also draws the real data and simulated heatmaps, and its projection for comparison, 
    * and prints the amount of data lost when the distributions are cut 
    * by the exclusion of bad/dead areas
    *
    * @param[in] realHeatmap heatmap from the real data
    * @param[in] simHeatmap heatmap from the PHENIX simulation
    * @param[in] realCutHeatmap heatmap from the real data with bad/dead areas cut
    * @param[in] simCutHeatmap heatmap from the PHENIX simulation with bad/dead areas cut
    * @param[in] detectorName name of the detector heatmap of which was passed 
    * @param[in] numberOfHeatmapDivisions number of divisions of heatmap for systematic ucnertainty evaluation 
    * @param[in] detectorName name of the detector heatmap of which was passed 
    * @param[in] title title that will be assigned to the histograms
    * @param[in] xTitle title of X axis that will be assigned to the histograms
    * @param[in] yTitle title of Y axis that will be assigned to the histograms
    * @param[in] rebinX rebin of X projection along X axis
    * @param[in] rebinY rebin of Y projection along Y axis
    * @param[in] drawYProj shows whether Y projection will be taken and drawn on the heatmap comparison canvas
    */
   double GetUncertainty(TH2F *&realHeatmap, TH2F *&simHeatmap,
                         TH2F *&realCutHeatmap, TH2F *&simCutHeatmap,
                         const int numberOfHeatmapDivisions,
                         const std::string& detectorName, const std::string& title,
                         const std::string& xTitle, const std::string& yTitle,
                         const int rebinX = 1, const int rebinY = 1, const bool drawYProj = true);

   /*! @brief Returns the uncertainty that is originated from the difference between 2 1D histograms
    *
    *  @param[in] realCutDistr real histogram projection with cut bad/dead areas
    *  @param[in] simCutDistr histogram projection from the simulation with cut bad/dead areas
    */
   void CheckHists(const TH2F *histReal, const TH2F *histSim, const std::string& name);
   /*! @brief Returns normalized ratio when the ratio should be between 0 and 1 and variables that are being ratioed can give the result bigger than 1 
    *
    *  @param[in] ratio ratio of 2 variables
    */
   double GetNormRatio(const double ratio);
};

#endif /* DEAD_MAP_SYS_HPP */
