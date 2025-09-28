/** 
 *  @file   EstimateGaussianBroadening.hpp 
 *  @brief  Contains declarations of functions and variables that are used for estimation gaussian broadening parameter sigma used as a parameter in gaus which is convoluted with breit-wigner
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysis).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef ESTIMATE_GAUSSIAN_BROADENING_HPP
#define ESTIMATE_GAUSSIAN_BROADENING_HPP

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

/*! @namespace EstimateGaussianBroadening
 * @brief Contains all functions and variables for EstimateGaussianBroadening.cpp
 */
namespace EstimateGaussianBroadening
{
   /*! Performs approximation for the invariant mass distribution for the given pT range
    *
    * @param[in] pTMin minimum bin of a pT range
    * @param[in] pTMax maximum bin of a pT range
    */
   void PerformMInvFit(const int pTBinMin, const int pTBinMax);
   /// Contents of input .yaml file for run configuration
   InputYAMLReader inputYAMLMain;
   /// Contents of input .yaml file for the information about resonance
   InputYAMLReader inputYAMLResonance;
   /// Name of run (e.g. Run14HeAu200 or Run7AuAu200)
   std::string runName;
   /// Output directory
   std::string outputDir;
   /// File in which widths of gausses will be written
   std::ofstream parametersOutput;
   /// Input file
   TFile *inputFile;
   /// name of the resonance
   std::string resonanceName;
   /// mass of the resonance [GeV/c^2]
   double massResonance;
   /// Graph containing widths of gausses of resonance signals [GeV/c^2]
   TGraphErrors grSigmas;
   /// histogram with counts vs invariant mass vs pT distribution
   TH2F *distr2DMInv;
   /// TText object template for quick text insertions
   TText text;
   /// TLatex object template for quick text insertions
   TLatex texText;
   /// Number of consequent fits of dphi and dz distributions for better approximation results
   /// each consequent fit decreases the limits around value from previous fit for every parameter
   /// which makes bettter gradual gradient descent of approximation parameters since ROOT built in
   /// approximation algorithm has only limited resource to perform the gradient descent
   const unsigned int fitNTries = 3;
   /// formula that is uded in TF1 constructor for approximation of sigmas
   const std::string fitSigmasFormula = "pol1";
};

int main(int argc, char **argv);

#endif /* ESTIMATE_ESTIMATE_GAUSSIAN_BROADENING_HPP */
