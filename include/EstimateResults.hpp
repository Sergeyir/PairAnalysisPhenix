/** 
 *  @file   EstimateResults.hpp 
 *  @brief  Contains declarations of functions and variables that are used for estimation of invariant pT spectra and nuclear modification factors R_{AB} and R_{CP}
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysis).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef ESTIMATE_RESULTS_HPP
#define ESTIMATE_RESULTS_HPP

#include <thread>

#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"

#include "StrTools.hpp"
#include "IOTools.hpp"
#include "MathTools.hpp"

#include "TCanvasTools.hpp"

#include "PBar.hpp"

#include "InputYAMLReader.hpp"
#include "FitFunc.hpp"

/*! @namespace EstimateResults
 * @brief Contains all functions and variables for EstimateSResultscpp
 */
namespace EstimateResults
{
   /// Contents of input .yaml file for run configuration
   InputYAMLReader inputYAMLMain;
   /// Contents of input .yaml file for the information about resonance
   InputYAMLReader inputYAMLResonance;
   /// Name of run (e.g. Run14HeAu200 or Run7AuAu200)
   std::string runName;
   /// Taxi number
   int taxiNumber;
   /// Name of the input file
   std::string inputFileName;
   /// Input file containing raw yields
   TFile *inputFile;
   /// Name of the reconstruction efficiency input file
   std::string inputRecEffFileName;
   /// Input file containing reconstruction efficiencies
   TFile *inputRecEffFile;
   /// Output file for writing invariant pT spectra
   TFile *outputFile;
   /// number of pT bins
   unsigned int pTNBins;
   /// pT bins ranges [GeV/c]
   std::vector<double> pTBinRanges;
   /// shows whether nuclear modification factors will be estimated (false for p+p, true for other collision systems)
   bool estimateFactors;
   /// TText object template for quick text insertions
   TText text;
   /// TLatex object template for quick text insertions
   TLatex texText;
   /// Number of consequent fits of spectra for better bin shift correction
   const unsigned int fitNTries = 5;
};

#endif /* ESTIMATE_RESULTS_HPP */
