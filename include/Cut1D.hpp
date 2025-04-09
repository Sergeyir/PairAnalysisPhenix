/** 
 *  @file   Cut1D.hpp 
 *  @brief  Contains declarations of functions to cut bad/dead areas from 1D distributions (slat, strip)
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/CalPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef CUT_1D_HPP
#define CUT_1D_HPP

#include <string>
#include <vector>
#include <fstream>

#include "TROOT.h"
#include "TError.h"
#include "TFile.h"
#include "TH1.h"
#include "TColor.h"
#include "TStyle.h"
#include "TCanvas.h"

#include "ErrorHandler.hpp"
#include "IOTools.hpp"

#include "TCanvasTools.hpp"

namespace Cut1D
{
   void CutDistribution(TH1D *distr, const std::string& name, 
                        const double minValueFromAverage);
}
int main();

#endif /* CUT_1D_HPP */
