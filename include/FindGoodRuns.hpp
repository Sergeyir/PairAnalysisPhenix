/** 
 *  @file   FindGoodRuns.hpp
 *  @brief  Contains declarations of functions and variables that are used for the determination of good runs 
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysisPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef FIND_GOOD_RUNS_HPP
#define FIND_GOOD_RUNS_HPP

#include <regex>

#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

#include "InputYAMLReader.hpp"

#include "IOTools.hpp"
#include "StrTools.hpp"
#include "Box.hpp"

#include "PBar.hpp"

#include "TCanvasTools.hpp"

#include "DeadMapCutter.hpp"

/* @namespace FindGoodRuns
 *
 * @brief Contains all functions, variables, and containers needed for FindGoodRuns.cpp
 *
 * This namespace is eployed so that documentation will not become a pile of variables, types, and functions from many different files that are intended to be compiled and used as executables. With this namespace finding the needed information for the given executable is easier since everything belongs to the current namespace
 */
namespace FindGoodRuns
{
   /// run name (e.g. Run14HeAu200)
   std::string runName;
   /// directory containing the run files
   std::string inputDir;
   /// Reference run number (e.g. the run with the highest mult/event ratio)
   int referenceRun;
   /// Contains all runs 
   std::vector<int> runs;
   /// Contains good runs
   std::vector<int> goodRuns;
   /// Contains bad runs
   std::vector<int> badRuns;
   /// Threshold for the absolute multiplicity/event and charge+/charge- deviation from the average
   double multThreshold = 2.;
   /// Chi2/NDF threshold for linear fit of the deviation heatmap from the reference run distributions 
   double chi2NDFThreshold = 3.;
   /// Checks all run files from goodRunNames by multiplicities and finds the reference run file.
   /// goodRunNames and badRunNames will be updated accordingly after the check
   void CheckRunsByMultiplicityAndFindReferenceRun();
}

#endif /* FIND_GOOD_RUNS_HPP */
