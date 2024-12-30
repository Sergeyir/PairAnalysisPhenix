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

#include <thread>

#include "TFile.h"
#include "TH1.h"

using namespace Run14HeAu200Cuts;

struct
{
   const std::string run = "Run14";
   const std::string system = "HeAu200";
   const bool doUseWeightFunc = true;
   bool doReweightAlpha = true;
   
   const std::string simDataDir = "data/Sim/";
   const std::string realDataDir = "data/Real/";
   const std::string outputDir = "data/PostSim/";
   
   std::vector<std::string> partQueue = {"pion", "apion"};//, "kaon", "apion", "akaon", "proton", "aproton"};
   std::vector<std::string> magfQueue = {""};
   //std::vector<std::string> magfQueue = {"+-", "-+"};
   std::vector<std::string> auxNameQueue = {"_lpt", "_hpt"};
   
   const double pTMin = 0.3;
   const double pTMax = 8.;

   //adc cut + efficiency correction of TOFw form AN814
   const double correctionTOFw = 0.799;
   
   const int pTNBins = static_cast<int>((pTMax - pTMin)*10.);
   
   const int nThreads = std::thread::hardware_concurrency();
   
   const std::string runName = run + system;
   
   TH1F *alphaReweightDCe0;
   TH1F *alphaReweightDCe1;
   TH1F *alphaReweightDCw0;
   TH1F *alphaReweightDCw1;
} Par;

#endif /*PAR_ANALYZE_HEAT_MAPS_HPP*/
