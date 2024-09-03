// $HEADER$
//------------------------------------------------------------------------------------------------
//                                  ParEmb declaration
//------------------------------------------------------------------------------------------------
// ParEmb - parameters of embedding
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic global container for storing parameters for embedding analysis
 **/
//------------------------------------------------------------------------------------------------

#ifndef PAR_EMB_HPP
#define PAR_EMB_HPP

#define DEAD_AREAS_CUTS_RUN7AUAU

#include <thread>

#include "../../Analysis/include/DeadAreasCuts.hpp"

struct
{   
   const std::string system = "AuAu200";
   const std::string run = "Run7";

   const std::string dataDir = "../data/Sim/";
   const std::string outputDir = "../data/postSim/";

   std::vector<std::string> magfQueue = {"+-", "-+"};
   std::vector<std::string> partQueue = {"pion", "kaon"};
   std::vector<std::string> centrQueue = {"00-20", "20-40", "40-60", "60-93", "00-93"};
   std::vector<std::string> partChargeQueue = {"", "a"};

   const double ptmin = 0.3;
   const double ptmax = 8.;

   const int nthreads = std::thread::hardware_concurrency();

   const unsigned int centrNBins = centrQueue.size();

   const std::string runName = run + system;
} Par;

#endif
