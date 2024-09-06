// $HEADER$
//------------------------------------------------------------------------------------------------
//                        Par struct declaration and realisation
//------------------------------------------------------------------------------------------------
// ParAnalyzeResonance
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic container for storing AnalyzeResonance macro parameters
 **/
//------------------------------------------------------------------------------------------------

#ifndef PAR_ANALYZE_RESONANCE_HPP
#define PAR_ANALYZE_RESONANCE_HPP

#define RUN7UAU //for deadmaps and centrality
#define RUN7UAU_MB5 //for centrality

#include <thread>

#include "Particles.hpp"
#include "M2Par.hpp"

#include "../../Analysis/include/DeadAreasCuts.hpp"
#include "../../Analysis/include/CentralityTypes.hpp"

struct
{
   const std::string system = "AuAu200";
   const std::string run = "Run7";
   
   KStar892 origParticle;

   std::vector<std::string> magfQueue = {"+-", "-+"};
   std::vector<std::string> auxNameQueue = {"lpt", "hpt"};

   std::vector<std::string> daughter1Queue = {"pion", "kaon"};
   std::vector<std::string> daughter2Queue = {"akaon", "apion"};

   bool doUseWeightFunc = true;
   
   std::vector<double> ptDeviationQueue = {1., 0.995, 1., 1.005};
   
   const double ptMin = 0.3;
   const double ptMax = 8.;

   const double ptMinPair = 0.9;
   const double ptMaxPair = 8.5;

   const int invMNBins = 300;

   //adc cut + efficiency correction of TOFw form AN814
   const double correctionTOFw = 0.799;

   std::array<std::string, 13> detectors = 
      {"dc_pc1", "pc2", "pc3", "tofe", "tofw", 
       "emcale0", "emcale1", "emcale2", "emcale3", 
       "emcalw0", "emcalw1", "emcalw2", "emcalw3"};
   
   const std::string runName = run + system ;

   const int nThreads = std::thread::hardware_concurrency();
   
   const int ptNBins = static_cast<int>((pair_ptmax - pair_ptmin)*10.);
} Par;

#endif /* PAR_ANALYZE_RESONANCE_HPP */
