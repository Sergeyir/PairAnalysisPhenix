/** 
 *  @file   EstimateResonanceEff.hpp 
 *  @brief  Contains declarations of functions and variables that are used for estimation of resonance reconstruction efficiency with the use of the data from MC
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysis).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef GAUSSIAN_BROADENING_HPP
#define GAUSSIAN_BROADENING_HPP

#include <thread>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TGraphErrors.h"
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

/*! @namespace EstimateResonanceEff
 * @brief Contains all functions and variables for EstimateResonanceEff.cpp
 */
namespace EstimateResonanceEff
{
   /*! Performs approximations of invariant mass distributions for all pT ranges for the given method
    *
    * @param[in] methodName name of the method that was used to extract pairs of charged tracsk
    */
   void PerformMInvFitsForMethod(const std::string& methodName);
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
   /// File in which widths of gausses will be written
   std::ofstream parametersOutput;
   /// Name of the input file
   std::string inputFileName;
   /// Input file
   TFile *inputFile;
   /// unscaled pT distribution of original generated particles
   TH1F *distrOrigUnscaledPT;
   /// pT distribution of original generated particles
   TH1F *distrOrigPT;
   /// name of the resonance
   std::string resonanceName;
   /// mass of the resonance [GeV/c^2]
   double massResonance;
   /// gamma of the resonance [GeV/c^2]
   double gammaResonance;
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
   const unsigned int fitNTries = 3;
};

int main(int argc, char **argv);

#endif /* ESTIMATE_GAUSSIAN_BROADENING_HPP */
