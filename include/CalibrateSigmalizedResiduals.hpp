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
#include "TH3.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TGraphErrors.h"

#include "IOTools.hpp"
#include "MathTools.hpp"

#include "TCanvasPrinter.hpp"

#include "PBar.hpp"

#include "InputReader.hpp"

struct
{
   std::unique_ptr<TFile> inputFile;
   const std::vector<int> zDCMin{-75, -60, -45, -30, -15, 0, 15, 30, 45, 60};
   const std::vector<int> zDCMax{-60, -45, -30, -15, 0, 15, 30, 45, 60, 75};
   
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
                                   2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0,
                                   3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0,
                                   5.2, 5.4, 5.6, 5.8, 6.0, 6.4, 6.8, 7.2, 7.6, 8.0,
                                   8.4, 8.8, 9.2, 9.6};
   const std::vector<double> pTMax{0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 
                                   1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 
                                   2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0,
                                   3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0,
                                   5.2, 5.4, 5.6, 5.8, 6.0, 6.4, 6.8, 7.2, 7.6, 8.0,
                                   8.4, 8.8, 9.2, 9.6, 10.0};
   
   const double minIntegralValue = 1e2; // minimum number of entries for 
                                        // the histogram to be approximated
   double centralityMin;
   double centralityMax;
   int centralityNBins;
} Par;

int main(int argc, char **argv);
void PerformFits(const std::string& runName, const std::string& detectorName, 
                 const std::string& variableName, const int zDCMin, const int zDCMax,
                 const bool isPositive);

#endif /* CALIBRATE_SIGMALIZED_RESIDUALS_HPP */
