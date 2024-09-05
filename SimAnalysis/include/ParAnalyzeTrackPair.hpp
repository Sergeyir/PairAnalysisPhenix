#pragma once

#define RUN7_AUAU

#include <thread>

#include "Particles.h"
#include "Run7AuAuM2Par.h"

#include "../../Analysis/include/DeadAreasCuts.hpp"
#include "../../Analysis/include/CentralityTypes.hpp"

struct
{
   const std::string system = "AuAu200";
   const std::string run = "Run7";
   
   KStar892 orig_part;

   std::vector<std::string> magf_queue = {"+-", "-+"};
   std::vector<std::string> aux_name_queue = {"lpt", "hpt"};

   std::vector<std::string> daughter1_queue = {"pion", "kaon"};
   std::vector<std::string> daughter2_queue = {"akaon", "apion"};

   AuAu200CTypeMB4 CType;

   bool do_use_weight_func = true;
   
   std::vector<double> pt_deviation_queue = {1., 0.995, 1., 1.005};
   
   const double ptmin = 0.3;
   const double ptmax = 8.;

   const double pair_ptmin = 0.9;
   const double pair_ptmax = 8.5;

   const int invm_nbins = 300;

   //adc cut + efficiency correction of TOFw form AN814
   const double tofw_correction = 0.799;

   std::array<std::string, 13> detectors = {"dc_pc1", "pc2", "pc3", "tofe", "tofw", 
      "emcale0", "emcale1", "emcale2", "emcale3", "emcalw0", "emcalw1", "emcalw2", "emcalw3"};
   
   const std::string run_name = run + system ;

   const int nthreads = std::thread::hardware_concurrency();
   
   const int pt_nbins = static_cast<int>((pair_ptmax - pair_ptmin)*10.);
} Par;
