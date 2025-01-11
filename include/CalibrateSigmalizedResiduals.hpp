// $HEADER$
//------------------------------------------------------------------------------------------------
//                      CalibrateSigmailzedResiduals functions declarations
//------------------------------------------------------------------------------------------------
// CalibrateSigmalizedResiduals
//
// ** Free and open code for anyone to use **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic code for calibrating sigmalized residuals (sdphi, sdz)
 **/
//------------------------------------------------------------------------------------------------

#ifndef CALIBRATE_SIGMALIZED_RESIDUALS_HPP
#define CALIBRATE_SIGMALIZED_RESIDUALS_HPP

#include <memory>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"

#include "IOTools.hpp"
#include "MathTools.hpp"

#include "TCanvasPrinter.hpp"

#include "PBar.hpp"

#include "InputReader.hpp"

struct
{
   const std::array<double, 10> zDCMin{-75, -60, -45, -30, -15, 0, 15, 30, 45, 60};
   const std::array<double, 10> zDCMax{-60, -45, -30, -15, 0, 15, 30, 45, 60, 75};

   const std::array<Color_t, 10> markerColor{kP10Cyan, kP10Ash, kP10Green, kP10Orange, kP10Brown, 
                                             kP10Violet, kP10Gray, kP10Red, kP10Yellow, kP10Blue};
   const std::array<Style_t, 10> markerStyle{24, 25, 27, 28, 26, 32, 24, 25, 27, 28};
   
   // Since it is difficult to approximate all distributions and get a reliable result
   // means and sigmas are sequentially approximated with more intricate functions and 
   // in wider pT range.
   const std::vector<std::string> meansFitFunc =
      {"[0] - [1]*exp(pol1(2)) + [4]*exp(pol1(5))", 
       "[0] + [1]/x + [2]/x^2 + [3]/x^3 + [4]/x^4 + pol4(5)"};
   const std::vector<std::string> sigmasFitFunc = 
      {"[0] - [1]*exp(pol1(2)) + [4]*exp(pol1(5))", 
       "pol2(0) + [3]*exp(pol2(4))"};
   
   /*
   const std::vector<double> pTMin{0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 
                                   1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 
                                   2.3, 2.6, 2.8, 3.0, 3.5, 4.0, 6.0, 8.0};
   const std::vector<double> pTMax{0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 
                                   1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.3, 
                                   2.6, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0};
                                   */
   const std::vector<double> pTMin{0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 
                                   1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 
                                   2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 3.0,
                                   3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0,
                                   5.2, 5.6, 6.0, 6.4, 6.8, 7.2, 7.6, 8.0,
                                   8.4, 9.0, 9.5};
   const std::vector<double> pTMax{0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 
                                   1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 
                                   2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 3.0,
                                   3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0,
                                   5.2, 5.6, 6.0, 6.4, 6.8, 7.2, 7.6, 8.0,
                                   8.4, 9.0, 9.5, 10.0};
   // number of pT bins whose ranges are listed above
   // this is needed for the canvas division
   const int pTXNBins = 8;
   const int pTYNBins = 6;

   // ranges of pT for sequential approximation of sigmas and means
   const std::vector<double> pTMinFit = {pTMin[2], pTMin.front(), pTMin.front()};
   const std::vector<double> pTMaxFit = {2.5, 4.0, pTMax.back()};
   
   const double minIntegralValue = 1e2; // minimum number of entries for 
                                        // the histogram to be approximated
                                        // if the requirement for this value is not met
                                        // warning will be printed
   double centralityMin;
   double centralityMax;
   int centralityNBins;
   
   const unsigned short fitNTries = 5;

   TLatex texText;
} Par;

int main(int argc, char **argv);
void PerformFits(TH3F *hist, TGraphErrors& grMeans, TGraphErrors& grSigmas, 
                 const std::string& outputFileNameNoExt, const std::string& dValName,
                 const std::string& detectorName, const std::string& zDCRangeName, 
                 const double centralityMin, const double centralityMax,
                 const std::string& chargeName);

#endif /* CALIBRATE_SIGMALIZED_RESIDUALS_HPP */
