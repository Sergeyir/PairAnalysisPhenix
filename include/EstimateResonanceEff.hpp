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

/*! @namespace EstimateResonanceEff
 * @brief Contains all functions and variables for EstimateResonanceEff.cpp
 */
namespace EstimateResonanceEff
{
   /*! Performs approximations of invariant mass distributions for all pT ranges for the given method
    *
    * @param[in] methodName name of the method that was used to extract pairs of charged tracsk
    */
   void PerformInvMassFits(const std::string& methodName);
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
   std::string nameResonance;
   /// mass of the resonance [GeV/c^2]
   double massResonance;
   /// Graph containing widths of Breit-Wigner approximations of resonance signals [GeV/c^2]
   TGraphErrors grGammas;
   /// Progress bar that shows progress in terminal
   ProgressBar pBar("FANCY", "", PBarColor::BOLD_CYAN);
   /// Number of consequent fits of dphi and dz distributions for better approximation results
   /// each consequent fit decreases the limits around value from previous fit for every parameter
   /// which makes bettter gradual gradient descent of approximation parameters since ROOT built in
   /// approximation algorithm has only limited resource to perform the gradient descent
   const unsigned int fitNTries = 5;
};

int main(int argc, char **argv);

#endif /* ESTIMATE_GAUSSIAN_BROADENING_HPP */
