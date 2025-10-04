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
#include "Constants.hpp"

/*! @namespace AnalyzeRealMInv
 * @brief Contains all functions and variables for AnalyzeRealMInv.cpp
 */
namespace AnalyzeRealMInv
{
   /*! Performs approximations of invariant mass distributions for all pT ranges for the given method
    *
    * @param[in] methodName name of the method that was used to extract pairs of charged tracsk
    */
   void PerformMInvFitsForMethod(const YAML::Node& method);
   /*! Merges invariant mass distributions with subtracted background for all centrality (c in CabanaBoy), z_{vtx} (z in CabanaBoy) and r_{vtx} (r in CabanaBoy)
    *
    * @param[in] methodName name of the method that was used to extract pairs of charged tracsk
    * @param[in] centralityBin centrality bin which will be used for invariant mass histogram merging
    * @param[in] pTBin pT bin which will be used for invariant mass histogram merging
    * @param[in] distrMInvMergedFG histogram to pass that will be filled with contents of all scaled foreground histograms. The value that is passed for this histogram must be nullptr.
    * @param[in] distrMInvMergedBG histogram to pass that will be filled with contents of all scaled background histograms. The value that is passed for this histogram must be nullptr.
    * @param[in] distrMInvMergedFGLR histogram to pass that will be filled with contents of all scaled foreground histograms with low resolution (for background scaling). The value that is passed for this histogram must be nullptr.
    * @param[in] distrMInvMergedBGLR histogram to pass that will be filled with contents of all scaled background histograms with low resolution (for background scaling). The value that is passed for this histogram must be nullptr.
    * @param[out] merged invariant mass distribution with background extracted
    */
   TH1D *MergeMInv(const std::string& methodName, const YAML::Node& centralityBin,
                   const int pTBin, TH1D*& distrMInvMergedFG, TH1D*& distrMInvMergedBG,
                   TH1D*& distrMInvMergedFGLR, TH1D*& distrMInvMergedBGLR);
   /*! Subtracts background for the specified histogram 
    *
    * @param[in] distrMInvFG foreground M_{inv} distribution from which background will be extracted
    * @param[in] distrMInvFG background M_{inv} distribution which will be extracted from foreground; in the process scaling will be applied
    * @param[out] invariant mass distribution with background subtracted
    */
   TH1D *SubtractBG(TH1D*& distrMInvFG, TH1D*& distrMInvBG, 
                    TH1D*& distrMInvFGLR, TH1D*& distrMInvBGLR);
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
   /// Number of a taxi
   int taxiNumber;
   /// Name of the input file
   std::string inputFileName;
   /// Input file
   TFile *inputFile;
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
   /// shows whether a resonance has an antiparticel (for which idDaugther1 and idDaughter2 are swapped)
   bool hasAntiparticle;
   /// M_{inv} range minimum value [GeV/c^2]
   double minMInv;
   /// M_{inv} range maximum value [GeV/c^2]
   double maxMInv;
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

#endif /* ANALYZE_REAL_M_INV_HPP */
