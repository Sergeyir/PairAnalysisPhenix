// $HEADER$
//------------------------------------------------------------------------------------------------
//                        Par struct declaration and realisation
//------------------------------------------------------------------------------------------------
// ParAnalyzeSingleTrack
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic container for storing AnalyzeSingleTrack macro parameters
 **/
//------------------------------------------------------------------------------------------------

#ifndef PAR_ANALYZE_SINGLE_TRACK_HPP
#define PAR_ANALYZE_SINGLE_TRACK_HPP

#define RUN7AUAU

#include <thread>

#include "TFile.h"

#include "Particles.hpp"

#include "DeadAreasCuts.hpp"

struct
{   
   const std::string system = "AuAu200";
   const std::string run = "Run7";
   bool doUseWeightFunc = true;
   
   const std::string simDataDir = "../data/Sim/";
   const std::string realDataDir = "../data/Analysis/";
   const std::string outputDir = "../data/PostSim/";
   
   std::vector<std::string> magfQueue = {"+-", "-+"};
   std::vector<std::string> partQueue = {"pion", "apion", "kaon", "akaon"};
   std::vector<std::string> auxNameQueue = {"_lpt", "_hpt"};
   
   const double ptMin = 0.3;
   const double ptMax = 8.;
   
   const int nthreads = std::thread::hardware_concurrency();

   const int ptNBins = static_cast<int>((ptMax - ptMin)*10.);
   
   const std::string runName = run + system;
} Par;

#endif /*PAR_ANALYZE_SINGLE_TRACK_HPP*/
