/** 
 *  @file   GoodRunsSelector.hpp
 *  @brief  Contains declarations of functions and variables that are used for the determination of good runs 
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysisPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef GOOD_RUNS_SELECTOR_HPP
#define GOOD_RUNS_SELECTOR_HPP

#include <signal.h>
#include <thread>

#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"

#include "InputYAMLReader.hpp"

#include "IOTools.hpp"
#include "StrTools.hpp"
#include "Box.hpp"

#include "ThrObj.hpp"

#include "SingleTrackFunc.hpp"
#include "SimTreeReader.hpp"
#include "DeadMapCutter.hpp"
#include "SimCalibrator.hpp"

#include "PBar.hpp"

#include "ROOT/TTreeProcessorMT.hxx"

int main(int argc, char **argv);

/* @namespace GoodRunsSelector
 *
 * @brief Contains all functions, variables, and containers needed for GoodRunsSelector
 *
 * This namespace is eployed so that documentation will not become a pile of variables, types, and functions from many different files that are intended to be compiled and used as executables. With this namespace finding the needed information for the given executable is easier since everything belongs to the current namespace
 */
namespace GoodRunsSelector
{
}

#endif /* GOOD_RUNS_SELECTOR_HPP */
