// $HEADER$
//------------------------------------------------------------------------------------------------
//                        Par struct delcalarion and realisation
//------------------------------------------------------------------------------------------------
// ParAnalyzeHeatMaps
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic container for storing AnalyzeHeatMaps macro parameters
 **/
//------------------------------------------------------------------------------------------------

#ifndef PAR_ANALYZE_HEAT_MAPS_HPP
#define PAR_ANALYZE_HEAT_MAPS_HPP

#define DEAD_AREAS_CUTS_RUN7AUAU

#include <thread>

#include "TFile.h"
#include "TH1.h"

#include "Particles.hpp"

#include "../../Analysis/include/DeadAreasCuts.hpp"

struct
{
   const std::string system = "AuAu200";
   const std::string run = "Run7";
   const bool doUseWeightFunc = true;
   bool doReweightAlpha = true;
   
   const std::string simDataDir = "../data/Sim/";
   const std::string realDataDir = "../data/Analysis/";
   const std::string outputDir = "../data/PostSim/";
   
   std::vector<std::string> partQueue = {"pion", "kaon", "apion", "akaon", "proton", "aproton"};
   std::vector<std::string> magfQueue = {"+-", "-+"};
   std::vector<std::string> auxNameQueue = {"_lpt", "_hpt"};
   
   const double ptMin = 0.3;
   const double ptMax = 8.;
   
   const int ptNBins = static_cast<int>((ptMax - ptMin)*10.);
   
   const int nthreads = std::thread::hardware_concurrency();
   
   const std::string runName = run + system;
   
   std::unique_ptr<TFile> realDataInputFile;
   std::unique_ptr<TFile> simInputFile;
   std::unique_ptr<TFile> alphaReweightInputFile;
   
   TH1F *alphaReweightDCe0;
   TH1F *alphaReweightDCe1;
   TH1F *alphaReweightDCw0;
   TH1F *alphaReweightDCw1;
} Par;

#endif /*PAR_ANALYZE_HEAT_MAPS_HPP*/
