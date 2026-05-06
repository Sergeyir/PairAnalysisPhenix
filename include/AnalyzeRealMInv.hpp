/** 
 *  @file   AnalyzeRealMInv.hpp 
 *  @brief  Contains declarations of functions and variables that are used for analyzis of experimental invariant mass distributions
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysis).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef ANALYZE_REAL_M_INV_HPP
#define ANALYZE_REAL_M_INV_HPP

#include <thread>
#include <filesystem>

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
#include "TGraph.h"

#include "StrTools.hpp"
#include "IOTools.hpp"
#include "MathTools.hpp"

#include "TCanvasTools.hpp"

#include "PBar.hpp"

#include "InputYAMLReader.hpp"
#include "FitFunc.hpp"
#include "Constants.hpp"

#include "MInv.hpp"

/*! @namespace AnalyzeRealMInv
 * @brief Contains all functions and variables for AnalyzeRealMInv.cpp
 */
namespace AnalyzeRealMInv
{
   /*! Performs approximations of invariant mass distributions for all pT ranges and all centrality classes for the given method
    *
    * @param[in] methodName name of the method that was used to extract pairs of charged tracsk
    */
   void PerformMInvFits(const YAML::Node& method);
   /*! Performs approximations of invariant mass distributions for all pT ranges for the given method and centrality. This function is implemented and used in GUI/MInvFit.cpp
    *
    * @param[in] methodName name of the method that was used to extract pairs of charged tracsk
    * @param[in] centralityBin centrality class bin that will be processed
    */
   void PerformMInvFits(const YAML::Node& method, const unsigned int centralityBin);
   /*! Returns TFile pointer if file exits and handles warning and info outputs
    *
    * @param[in] inputFileName .root file containing fixed BG fits. Whether fits are in this file this function does not check
    * @param[in] fitTypeName name of the fit type. This variable is needed to make the warning and info output clear
    * @param[in] printFreeFitWarning prints warning that free BG fits will be used instead of fixed ones. This parameter doesn't change the behaviour of the program to perform free or fixed fit and is only needed for clarity when running the executable. The behaviour is needed to be implemented by the user.
    *
    * @param[out] TFile pointer. Returns nullptr if file doesn't exist
    */
   TFile *SetFixedBGFile(const std::string& inputFileName, const std::string& fitTypeName, 
                         const bool printFreeFitWarning = true);
   /*! Set the parameters of BG fit
    *
    * @param[in] inputFile file from which approximation parameters will be written
    * @param[in] fitBG background fit to which BG approximation parameters will be applied
    * @param[in] identifier string containing all important information about the fit type and pT range. This value is not required to be absolutely correct and is only used to print warning if something goes wrong to help pinpoint the problem if there is one.
    *
    * @param[out] readStatus value that shows whether assigning fit parameters values from inputFile was succesfull (true if everything was fine, otherwise false)
    */
   bool SetBGFit(TFile *&inputFile, TF1 *&fitBG, const std::string& identifier);
   /// Sets parameters for a function needed for estimating width of 
   /// gaus for convolution of Gaus and Breit-Wigner
   void SetGaussianBroadeningFunction();
   /* Extracts the yield by integrating the distribution and subtracting the background in the specified range
    *
    * @param[in] distrInvM invariant mass distribution from which the yield will be calculated
    * @param[in] funcBG function that approximates the background
    * @param[in] xMin minimum M_{inv} value of an extraction range [GeV/c^2]
    * @param[in] xMax maximum M_{inv} value of an extraction range [GeV/c^2]
    */
   double GetYield(TH1D *distr, TF1 *funcBG, const double xMin, const double xMax);
   /// Contents of input .yaml file for run configuration
   InputYAMLReader inputYAMLMain;
   /// Contents of input .yaml file for the information about resonance
   InputYAMLReader inputYAMLResonance;
   /// Name of run (e.g. Run14HeAu200 or Run7AuAu200)
   std::string runName;
   /// Number of a taxi
   int taxiNumber;
   /// Name of the input file (for getting invariant mass distributions from the real data)
   std::string inputFileName;
   /// Input file (for getting invariant mass distributions from the real data)
   TFile *inputFile;
   /// Directory to which files containing all parameters and yields will be written
   std::string parametersOutputDir;
   /// File in which all parameters and yields will be written
   TFile *parametersOutputFile;
   /// name of the resonance
   std::string resonanceName;
   /// mass of the resonance [GeV/c^2]
   double massResonance;
   /// gamma of the resonance [GeV/c^2]
   double gammaResonance;
   /// id of 1st decay product
   int daughter1Id;
   /// id of 2nd decay product
   int daughter2Id;
   /// M_{inv} range minimum value [GeV/c^2]
   double minMInv;
   /// M_{inv} range maximum value [GeV/c^2]
   double maxMInv;
   /// number of pT bins
   unsigned int pTNBins;
   /// Rebin value for M_{inv} i.e. x axis 
   int rebinX = 1;
   /// yield extraction range in +-(Gamma + sigma)*sigmalizedYieldExtractionRange from mean
   double sigmalizedYieldExtractionRange;
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
   /// Number of consequent fits of invariant mass distributions for better approximation results
   /// each consequent fit decreases the limits around value from previous fit for every parameter
   /// which makes bettter gradual gradient descent of approximation parameters since ROOT built in
   /// approximation algorithm has only limited resource to perform the gradient descent
   const unsigned int fitNTries = 1;
};

#endif /* ANALYZE_REAL_M_INV_HPP */
