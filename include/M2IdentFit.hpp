/** 
 *  @file   M2IdentFit.hpp 
 *  @brief  Contains declarations of functions and variables that are used for approximation of m2 distributions for the estimation of timing parameters that are used for identification via m2
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/CalPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef M2_IDENT_FIT_HPP
#define M2_IDENT_FIT_HPP

#include <cmath>
#include <string>
#include <filesystem>

#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

#include "ErrorHandler.hpp"
#include "IOTools.hpp"
#include "StrTools.hpp"
#include "Box.hpp"

#include "TCanvasTools.hpp"

#include "InputYAMLReader.hpp"


int main(int argc, char **argv);

/* @namespace M2IdentFit
 *
 * @brief Contains all functions, variables, and containers needed for M2IdentFit
 *
 * This namespace is eployed so that documentation will not become a pile of variables, types, and functions from many different files that are intended to be compiled and used as executables. With this namespace finding the needed information for the given executable is easier since everything belongs to the current namespace
 */
namespace M2IdentFit
{
   /* @brief Performs all fits for charged hadrons m2 distribution 
    * for the given detector and for the given centrality class
    *
    * @param[in] detector yaml node for the given detector
    * @param[in] centralityMin lowest bound of centrality range
    * @param[in] centralityMax highest bound of centrality range
    */
   void PerformFitsForDetector(const YAML::Node& detector, 
                               const double centralityMin,
                               const double centralityMax);
   /* @brief Performs m2 fits for charged hadrons for the given histogram
    *
    * @param[in] pTMin minimum pT for the current bin [GeV/c]
    * @param[in] pTMax maximum pT for the current bin [GeV/c]
    * @param[in] massProj projection of m2 histogram distribution taken from the real data
    * @param[in] fitPar approximation data container in which the data for the current pT bin will be written to
    */
   void PerformSingleM2Fit(const double pTMin, const double pTMax, TH1F *massProj, 
                           FitParameters& fitPar, const std::string& funcBG, const double distance, 
                           const double sigmaAlpha, const double sigmaMS, const double sigmaT);
   /* @brief Calculates the yield of the particle from the m2 distribution. Since the distribution
    * is discrete the yield will may be extracted in the range that is narrower than specified to
    * avoid subtraction of the wider range. The difference in range is corrected by the ratio of
    * error functions of the extracted range to the specified range. The latter is performed under
    * the assumption that all signals have purely gaussian shape which is not true and needs to
    * be accounted by the evaluation of systematic uncertainty. Yields must be extracted with
    * this function in multiple different ranges for this estimation.
    *
    * @param[in] hist the histogram the yield from which will be calculated
    * @param[in] mean mean of the gaussian distribution for the given particle 
    * @param[in] sigma sigma of the gaussian distribution for the given particle 
    * @param[in] sigmalizedYieldExtractionRange sigmalized yield extraction range 
    * (i.e. if this variable is set to 2 then yield will be calculated from the integral
    * in the range from mean - 2*sigma to mean + 2*sigma)
    * @param[in] fitGaus1 gaussian approximation of a particle different than specified
    * @param[in] fitGaus2 gaussian approximation of a particle different than specified
    * and different than the particle which was used for fitGaus1
    * @param[in] fitBG approximation of background
    * @param[in] vetoLow lowest bound on m2 identification range that comes from other particle
    * @param[in] vetoHigh highest bound on m2 identification range that comes from other particle
    * @param[in] err error of the yield
    * @param[out] yield of the particle in the given range
    */
   double GetYield(TH1F *hist, const double mean, const double sigma,
                   const double sigmalizedYieldExtractionRange,
                   TF1 *fitGaus1, TF1 *fitGaus2, TF1 *fitBG,
                   const double vetoLow, const double vetoHigh);
   /* @struct FitParameters
    *
    * @brief Contains fit parameters for different pT
    */
   struct FitParameters
   {
      /// @brief Default deleted constructor
      FitParameters() = delete;
      /* @brief Constructor for defining the fit parameters
       * @param[in] isChargePositive shows whether the charge is positive
       */
      FitParameters(const std::string& particleName);

      /// means vs pT
      TGraph meansVsPT;
      /// sigmas vs pT
      TGraph sigmasVsPT;
      /// means vs pT fit for
      std::unique_ptr<TF1> meansVsPTFit;
      /// sigmas vs pT fit for
      std::unique_ptr<TF1> sigmasVsPTFit;
      /// file in which raw yields
      ofstream rawYieldsOutputFile;
   };
   /// file with real data
   TFile *inputDataFile;
   /// file reader for all required parameters for the m2 identification
   InputYAMLReader inputYAMLM2Id;
   /// file reader for all required parameters for the current run
   InputYAMLReader inputYAMLMain;
   /// name of a run (i.e. Run14HeAu200)
   std::string runName;
   /// name of a collision system (i.e. HeAu200)
   std::string collisionSystemName;
   /// directory in which all output files will be written
   std::string outputDir;
   /// directory in which all approximation parameters will be written
   std::string parametersDir;
   /// directory in which all yields will be written
   std::string rawYieldsDir;
   /// number of sequential fits with regressive parameter limiter 
   /// for the improvement of approximation
   const unsigned int nFitTries = 5.;
}

#endif /* M2_IDENT_FIT_HPP */
