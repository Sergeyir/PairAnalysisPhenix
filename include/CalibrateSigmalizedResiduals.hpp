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

#include <array>

#include "json/json.h"
#include "json/value.h"

#include "TFile.h"
#include "TH1.h"
#include "TH3.h"
#include "TF1.h"
#include "TROOT.h"
#include "TGraphErrors.h"

#include "ErrorHandler.hpp"
#include "IOTools.hpp"

#include "MathTools.hpp"

#include "TCanvasPrinter.hpp"

int main(int argc, char **argv);
void PerformFits(const std::string& runName, const std::string& detectorName, 
                 const std::string& variableName, const int zDCMin, const int zDCMax,
                 const bool isPositive);

#endif /* CALIBRATE_SIGMALIZED_RESIDUALS_HPP */
