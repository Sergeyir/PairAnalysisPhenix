/** 
 *  @file   CheckRuns.hpp
 *  @brief  Contains declarations of functions and variables that are used for the determination of good runs 
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysisPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef CHECK_RUNS_HPP
#define CHECK_RUNS_HPP

#include <regex>

#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLine.h"

#include "InputYAMLReader.hpp"

#include "IOTools.hpp"
#include "StrTools.hpp"
#include "Box.hpp"

#include "PBar.hpp"

#include "TCanvasTools.hpp"

#include "DeadMapCutter.hpp"

/* @namespace CheckRuns
 *
 * @brief Contains all functions, variables, and containers needed for CheckRuns.cpp
 *
 * This namespace is eployed so that documentation will not become a pile of variables, types, and functions from many different files that are intended to be compiled and used as executables. With this namespace finding the needed information for the given executable is easier since everything belongs to the current namespace
 */
namespace CheckRuns
{
   /// file reader for all required parameters for the current run
   InputYAMLReader inputYAMLMain;
   /// run name (e.g. Run14HeAu200)
   std::string runName;
   /// directory containing the run files
   std::string inputDir;
   /// directory for writing pictures and .root files
   std::string outputDir;
   /// simulated file name
   std::string simFileName;
   /// Contains all runs 
   std::vector<int> runs;
   /// Contains good runs
   std::vector<int> goodRuns;
   /// Contains bad runs
   std::vector<int> badRuns;
   /// cutter for deadmaps
   DeadMapCutter dmCutter;
   /// Threshold for the absolute multiplicity/event and charge+/charge- deviation from the average
   double multThreshold = 1.3;
   /// const fit parameter deviation threshold from unity for linear fit of the deviation heatmap from the simulated distributions
   double constParDeviationThreshold = 1.1;
   /// Chi2/NDF threshold for linear fit of the deviation heatmap from the simulated distributions
   double chi2NDFThreshold = 3.;
   /// Checks all run files from goodRunNames by multiplicities and finds the simulated file.
   /// goodRunNames and badRunNames will be updated accordingly after the check
   void CheckRunsByMultiplicity();
   /// Checks all run files from goodRunNames by DC projection on board.
   /// goodRunNames and badRunNames will be updated accordingly after the check
   void CheckRunsByDCBoard();
}

#endif /* CHECK_RUNS_HPP */
