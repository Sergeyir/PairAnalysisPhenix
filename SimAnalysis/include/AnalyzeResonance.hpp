// $HEADER$
//------------------------------------------------------------------------------------------------
//                     AnalyzeTrackPair function declaration
//------------------------------------------------------------------------------------------------
// AnalyzeTrackPair
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic macro used for evaluation of registration of pair of tracks for different methods
 * from simulation output of event-like TTrees to processed histograms 
 * for further track pair registering correction evaluation
 **/
//------------------------------------------------------------------------------------------------

#ifndef ANALYZE_RESONANCE_HPP
#define ANALYZE_RESONANCE_HPP

#include "STrackFun.hpp"
#include "PTrackFun.hpp"
#include "ParAnalyzeResonance.hpp"
#include "../lib/Ident.h"

#include "../lib/EffTreeReader.h"

#include "PBar.hpp"
#include "../lib/Tools.h"
#include "../lib/InputTools.h"
#include "../lib/Box.h"

#include "ROOT/TTreeProcessorMT.hxx"

struct ThrContainer
{
   ThrObj<TH1F> origPtDistr = ThrObj<TH1F>
      ("orig","orig", Par.pt_nbins, Par.pair_ptmin, Par.pair_ptmax);

   ThrObj<TH2F> orig_pt_vs_pt = ThrObj<TH2F>
      ("orig_pt_vs_pt", "orig vs pt", 
       Par.pt_nbins, Par.ptmin, Par.ptmax, 
       Par.pt_nbins, Par.pair_ptmin, Par.pair_ptmax);
   
   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> InvM_noPID;
   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> InvM_noPID_acc_decreased;
   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> InvM_noPID_acc_increased;

   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> InvM_1PID;
   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> InvM_1PID_acc_decreased;
   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> InvM_1PID_acc_increased;

   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> InvM_2PID;
   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> InvM_2PID_acc_decreased;
   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> InvM_2PID_acc_increased;
   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> InvM_2PID_m2_eff_decreased;
   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> InvM_2PID_m2_eff_increased;

   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> InvM_TOF2PID;
   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> InvM_TOF2PID_acc_decreased;
   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> InvM_TOF2PID_acc_increased;

   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> InvM_EMC2PID;
   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> InvM_EMC2PID_acc_decreased;
   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> InvM_EMC2PID_acc_increased;
   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> InvM_EMC2PID_m2_eff_decreased;
   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> InvM_EMC2PID_m2_eff_increased;
};

struct ParticleContainer
{
   double mass;
   double iter;
   int orig_id;
   int geant_id;
   
   std::unique_ptr<double> dc_pc1_emb;
   std::unique_ptr<double> pc2_emb;
   std::unique_ptr<double> pc3_emb;
   std::unique_ptr<double> tofe_emb;
   std::unique_ptr<double> tofw_emb;
   
   std::array<std::unique_ptr<double>, 4> emcale_emb;
   std::array<std::unique_ptr<double>, 4> emcalw_emb;

   std::unique_ptr<double> emcale_m2_eff;
   std::unique_ptr<double> emcalw_m2_eff;
   std::unique_ptr<double> emcale_m2_eff_sys;
   std::unique_ptr<double> emcalw_m2_eff_sys;
   
   int size;
   double mom[3];
   std::array<int, 50> index;
   
   std::array<int, 50> pc2_id;
   std::array<int, 50> pc3_id;
   std::array<int, 50> tof_id;
   std::array<int, 50> emc_id;
   
   std::array<std::array<double, 50>, Par.CType.size> pc2_weight;
   std::array<std::array<double, 50>, Par.CType.size> pc3_weight;
   std::array<std::array<double, 50>, Par.CType.size> tof_weight;
   std::array<std::array<double, 50>, Par.CType.size> emc_weight;
   
   std::array<std::array<double, 50>, Par.CType.size> tof_id_weight;
   std::array<std::array<double, 50>, Par.CType.size> emc_id_weight;
   
   void ResetTrack(const int i)
   {
      index[i] = 0.;
      
      pc2_id[i] = PartId.junk;
      pc3_id[i] = PartId.junk;
      tof_id[i] = PartId.junk;
      emc_id[i] = PartId.junk;
      
      for (int j = 0; j < Par.CType.size; j++)
      {
         pc2_weight[j][i] = 0.;
         pc3_weight[j][i] = 0.;
         tof_weight[j][i] = 0.;
         emc_weight[j][i] = 0.;
         
         tof_id_weight[j][i] = 0.;
         emc_id_weight[j][i] = 0.;
      }
   }
};

void AnalyzeConfiguration(ThrContainer *thrContainer, const std::string& daughter1, 
                          const std::string& daughter2, const std::string& magf, 
                          const std::string& auxName, const double ptDeviation, 
                          const int procNum);
void AnalyzeResonance();
int main();

#endif /* ANALYZE_RESONANCE_HPP */
