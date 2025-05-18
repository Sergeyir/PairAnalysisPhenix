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
/* @struct M2IdentFit
 *
 * @brief Contains fit parameters for different pT
 */
   struct FitParameters
   {
      TGraph meansVsPT;
      TGraph sigmasVsPT;

      TF1 *meansVsPTFit;
      TF1 *sigmasVsPTFit;
   }
   /// @brief Calculates the yield of the particle from the m2 distribution
   double GetYield(TH1F *hist, const double mean, const double sigma, 
                   const double vetoLow, const double vetoHigh);

   void PerformFitsForDetector(const YAML::Node& detector);
   void PerformSingleM2Fit(const double pT, TH1F *massProj, FitParameters& fitPar);
}

#endif /* M2_IDENT_FIT_HPP */
