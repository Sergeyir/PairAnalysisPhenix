// $HEADER$
//------------------------------------------------------------------------------------------------
//                         Par struct declaration and realisation
//------------------------------------------------------------------------------------------------
// ParAnalyzeEmbedding - analyze embedding parameters
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic container for storing AnalyzeEmbedding macro parameters
 **/
//------------------------------------------------------------------------------------------------

#ifndef PAR_ANALYZE_EMBEDDING_HPP
#define PAR_ANALYZE_EMBEDDING_HPP

#include <thread>

#include "RunConfiguration.hpp"
#include "DeadAreasCuts.hpp"

struct
{   
   const std::string run = NUMBERED_RUN_NAME;
   const std::string system = COLLISION_SYSTEM_NAME;

   const std::string dataDir = "../data/Sim/";
   const std::string outputDir = "../data/PostSim/";

   std::vector<std::string> magfQueue = {"+-", "-+"};
   std::vector<std::string> partQueue = {"pion", "kaon"};
   std::vector<std::string> centrQueue = {"00-20", "20-40", "40-60", "60-93", "00-93"};
   std::vector<std::string> partChargeQueue = {"", "a"};

   const double ptMin = 0.3;
   const double ptMax = 8.;

   const int nthreads = std::thread::hardware_concurrency();

   const unsigned int centrNBins = centrQueue.size();

   const std::string runName = run + system;
} Par;

#endif /*PAR_ANALYZE_EMBEDDING_HPP*/
