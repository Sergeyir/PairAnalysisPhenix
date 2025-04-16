/** 
 *  @file   SigmalizedResiduals.hpp 
 *  @brief  Contains declarations of functions and variables that are used for estimation of values for calibration of sigmalized residuals sdphi and sdz from dphi and dz values from the PHENIX simulation
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysis).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef SIGMALIZED_RESIDUALS_HPP
#define SIGMALIZED_RESIDUALS_HPP

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

#include "IOTools.hpp"
#include "MathTools.hpp"

#include "TCanvasTools.hpp"

#include "InputYAMLReader.hpp"

/*! @namespace SigmalizedResiduals
 * @brief Contains all functions and containers for SigmalizedResiduals.cpp and CheckSigmalizedResiduals.cpp
 */
namespace SigmalizedResiduals
{
/*! @brief Performs sdz and sdphi calibrations from dphi and dz values for the given detector from the PHENIX simulation
 * @param[in] detectorName name of the detector for which the calibrations will be performed
 * @param[in] variableName name of the variable the calibration will be performed for 
 * @param[in] charge the charge of the tracks the calibration will be performed for
 */
   void PerformCalibrationsForDetector(const std::string& detectorName, 
                                       const std::string& variableName, const int charge);
   /// Contents of input .yaml file for run configuration
   InputYAMLReader inputYAMLMain;
   /// Name of run (e.g. Run14HeAu200 or Run7AuAu200)
   std::string runName;
   // Charges of particles to be analyzed independently
   const std::array<int, 2> particleCharges{1, -1};
   /// Names of variables to be calibrated
   std::array<std::string, 2> variableName{"dphi", "dz"};
   /// Names of variables to be calibrated in LaTex format
   std::array<std::string, 2> variableNameTex{"d#varphi", "dz_{DC}"};
   /// Output directory
   std::string outputDir;
   /// Input file
   std::unique_ptr<TFile> inputFile;
   /// Minimum pT of the whole pT range
   double pTMin;
   /// Maximum pT of the whole pT range
   double pTMax;
   /// Number of consequent fits of dphi and dz distributions for better approximation results
   /// each consequent fit decreases the limits around value from previous fit for every parameter
   /// which makes bettter gradual gradient descent of approximation parameters since ROOT built in
   /// approximation algorithm has only limited resource to perform the gradient descent
   /// This value will be read and updated from .json calibration input file
   const unsigned int fitNTries = 5;
};
int main(int argc, char **argv);

#endif /* SIGMALIZED_RESIDUALS_HPP */
