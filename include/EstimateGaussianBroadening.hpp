/** 
 *  @file   EstimateGaussianBroadening.hpp 
 *  @brief  Contains declarations of functions and variables that are used for estimation gaussian broadening parameter sigma that used as a parameter in gaus which is convoluted with breit-wigner
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
#include "TH3.h"
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
    * @param[in] pTMin minimum value of a pT range [GeV/c]
    * @param[in] pTMax maximum value of a pT range [GeV/c]
    */
   void PerformInvMassFit(const double pTmin, const double pTMax);
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
   /// Graph containing widths of gausses
   TGraphErrors grWidths;
   /// Number of consequent fits of dphi and dz distributions for better approximation results
   /// each consequent fit decreases the limits around value from previous fit for every parameter
   /// which makes bettter gradual gradient descent of approximation parameters since ROOT built in
   /// approximation algorithm has only limited resource to perform the gradient descent
   const unsigned int fitNTries = 5;
};

int main(int argc, char **argv);

#endif /* ESTIMATE_GAUSSIAN_BROADENING_HPP */
