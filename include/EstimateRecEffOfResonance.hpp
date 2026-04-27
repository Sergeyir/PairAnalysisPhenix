/** 
 *  @file   EstimateRecEffOfResonance.hpp 
 *  @brief  Contains declarations of functions and variables that are used for estimation of resonance reconstruction efficiency with the use of the data from MC
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysis).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef ESTIMATE_REC_EFF_OF_RESONANCE_HPP
#define ESTIMATE_REC_EFF_OF_RESONANCE_HPP

#include <thread>
#include <regex>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLine.h"
#include "TLegend.h"
#include "TMath.h"

#include "StrTools.hpp"
#include "IOTools.hpp"
#include "MathTools.hpp"

#include "TCanvasTools.hpp"

#include "PBar.hpp"

#include "InputYAMLReader.hpp"
#include "FitFunc.hpp"

/*! @namespace EstimateRecEffOfResonance
 * @brief Contains all functions and variables for EstimateRecEffOfResonance.cpp
 */
namespace EstimateRecEffOfResonance
{
   /*! Performs approximations of invariant mass distributions for all pT ranges for the given method
    *
    * @param[in] methodName name of the method that was used to extract pairs of charged tracks
    */
   void PerformMInvFitsForMethod(const std::string& methodName);
   /*! Performs approximations of invariant mass distributions for the given histogram
    *
    * @param[in] pTBin pT bin index
    * @param[in] methodName name of the method
    * @param[in] file file from which needed distributions will be read
    * @param[in] distrRecEffVsPT histogram containing information about reconstruction efficiency vs pT; for the current pTBin the information will be updated after fit is performed
    * @param[in] distrMeansVsPT histogram containing information about means vs pT; for the current pTBin the information will be updated after fit is performed
    * @param[in] distrGammasVsPT histogram containing information about gammas vs pT; for the current pTBin the information will be updated after fit is performed
    * @param[in] outputFileNameWithoutExt file name without extention in which pictures will be written (.png and .pdf). If empty string is specified (as is by default) no pictures will be saved.
    */
   void PerformMInvFit(const unsigned int pTBin, const std::string& methodName, TFile* file,
                       TH1D& distrRecEffVsPT, TH1D& distrMeansVsPT, TH1D& distrGammasVsPT,
                       const std::string& outputFileNameWithoutExt = "");
   /// Sets parameters for a function needed for estimating width of 
   /// gaus for convolution of Gaus and Breit-Wigner
   void SetGaussianBroadeningFunction();
   /* Extracts the yield by integrating the distribution and subtracting the background in the specified range
    *
    * @param[in] distrInvM invariant mass distribution from which the yield will be calculated
    * @param[in] funcBG function that approximates the background
    * @param[in] xMin minimum M_{inv} value of an extraction range [GeV/c^2]
    * @param[in] xMax maximum M_{inv} value of an extraction range [GeV/c^2]
    * @param[in] err yield statistical uncertainty
    */
   double GetYield(TH1D *distr, const TF1& funcBG, 
                   const double xMin, const double xMax, double &err);
   /// Contents of input .yaml file for run configuration
   InputYAMLReader inputYAMLMain;
   /// Contents of input .yaml file for the information about resonance
   InputYAMLReader inputYAMLResonance;
   /// Name of run (e.g. Run14HeAu200 or Run7AuAu200)
   std::string runName;
   /// Name of the input file
   std::string inputFileName;
   /// Input file
   TFile *inputFile;
   /// Names of alternative pTScale simulation input files for systematic uncertainty evaluation
   std::vector<std::string> altPTScaleSimInputFileNames;
   /// Alernative pTScale simulation input files for systematic uncertainty evaluation
   std::vector<TFile *> altPTScaleSimInputFiles;
   /// Output file (for writing means, gammas, efficiency reconstruction, etc. vs pT)
   TFile *outputFile;
   /// name of the resonance
   std::string resonanceName;
   /// mass of the resonance [GeV/c^2]
   double massResonance;
   /// gamma of the resonance [GeV/c^2]
   double gammaResonance;
   /// yield extraction range in +-(Gamma + sigma)*sigmalizedYieldExtractionRange from mean
   double sigmalizedYieldExtractionRange;
   /// number of pT bins
   unsigned int pTNBins;
   /// pT bins ranges [GeV/c]
   std::vector<double> pTBinRanges;
   /// function for estimating width of gaus for convolution of Gaus and Breit-Wigner
   TF1 *gaussianBroadeningEstimatorFunc;
   /// TText object template for quick text insertions
   TText text;
   /// TLatex object template for quick text insertions
   TLatex texText;
   /// Progress bar that shows progress in terminal
   ProgressBar pBar("FANCY", "", PBarColor::BOLD_CYAN);
   /// Number of calls in an iteration. Needed by pBar
   unsigned long numberOfCalls = 0;
   /// Overal number of iterations. Needed by pBar
   unsigned long numberOfIterations = 0;
   /// Number of consequent fits of dphi and dz distributions for better approximation results
   /// each consequent fit decreases the limits around value from previous fit for every parameter
   /// which makes bettter gradual gradient descent of approximation parameters since ROOT built in
   /// approximation algorithm has only limited resource to perform the gradient descent
   const unsigned int fitNTries = 0;
};

#endif /* ESTIMATE_REC_EFF_OF_RESONANCE_HPP */
