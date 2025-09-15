/** 
 *  @file   EstimateSingleTrackEff.hpp 
 *  @brief  Contains declarations of functions and variables that are used for estimation of registering/identification of charged tracks of pions, kaons, and protons in MC
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/CalPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef ESTIMATE_SINGLE_TRACK_EFF_HPP
#define ESTIMATE_SINGLE_TRACK_EFF_HPP

#include <cmath>
#include <string>
#include <filesystem>

#include "TROOT.h"
#include "TError.h"
#include "TFile.h"
#include "TStyle.h"
#include "TH1.h"
#include "TLine.h"
#include "TAttLine.h"
#include "TColor.h"
#include "TText.h"
#include "TLegend.h"

#include "ErrorHandler.hpp"
#include "IOTools.hpp"
#include "StrTools.hpp"
#include "MathTools.hpp"

#include "TCanvasTools.hpp"

#include "InputYAMLReader.hpp"

int main(int argc, char **argv);

/* @namespace EstimateSingleTrackEff
 *
 * @brief Contains all functions, variables, and containers needed for EstimateSingleTrackEff
 *
 * This namespace is eployed so that documentation will not become a pile of variables, types, and functions from many different files that are intended to be compiled and used as executables. With this namespace finding the needed information for the given executable is easier since everything belongs to the current namespace
 */
namespace EstimateSingleTrackEff
{
   /* @brief Performs estimation of registering/identification efficiency for the given detector
    *
    * @param[in] distrName name of the detector
    * @param[in] isIdentification if true distributions of identified tracks will be used; otherwise distributions of registered tracks
    */
   void EstimateEffForSingleDetector(const std::string& detectorName, 
                                     const bool isIdentication = false);
   /// file reader for all required parameters for single track MC info
   InputYAMLReader inputYAMLSim;
   /// file reader for all required parameters for the current run
   InputYAMLReader inputYAMLMain;
   /// name of a run (i.e. Run14HeAu200)
   std::string runName;
   /// pT range lower bound
   double pTMin;
   /// pT range upper bound
   double pTMax;
   /// file with processed MC \pi^{+} data
   TFile *inputDataFilePiPlus;
   /// file with processed MC K^{+} data
   TFile *inputDataFileKPlus;
   /// file with processed MC p data
   TFile *inputDataFileP;
   /// file with processed MC \pi^{-} data
   TFile *inputDataFilePiMinus;
   /// file with processed MC K^{-} data
   TFile *inputDataFileKMinus;
   /// file with processed MC \bar{p} data
   TFile *inputDataFilePBar;
   /// directory in which all output files will be written
   std::string outputDir;
   /// original pT distibution of \pi^{+}
   TH1F *distrOrigPTPiPlus;
   /// original pT distibution of K^{+}
   TH1F *distrOrigPTKPlus;
   /// original pT distibution of p
   TH1F *distrOrigPTP;
   /// original pT distibution of \pi^{-}
   TH1F *distrOrigPTPiMinus;
   /// original pT distibution of K^{-}
   TH1F *distrOrigPTKMinus;
   /// original pT distibution of \bar{p}
   TH1F *distrOrigPTPBar;
}

#endif /* ESTIMATE_SINGLE_TRACK_EFF_HPP */
