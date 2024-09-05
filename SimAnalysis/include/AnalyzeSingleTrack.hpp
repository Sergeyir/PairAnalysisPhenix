// $HEADER$
//------------------------------------------------------------------------------------------------
//                              AnalyzeSingleTrack function declaration
//------------------------------------------------------------------------------------------------
// STrackStudy - single track study
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic macro used for evaluation of registration and/or identification of single tracks
 * from simulation output of event-like TTrees to processed histograms
 * for further efficiency evaluation
 **/
//------------------------------------------------------------------------------------------------

#ifndef ANALYZE_SINGLE_TRACK_HPP
#define ANALYZE_SINGLE_TRACK_HPP

#include "ROOT/TTreeProcessorMT.hxx"

#include "Box.hpp"
#include "IOTools.hpp"

#include "ThrObj.hpp"

#include "PBar.hpp"

#include "STrackFun.hpp"
#include "EffTreeReader.hpp"
#include "ParAnalyzeSingleTrack.hpp"

struct ThrContainer
{
   ThrObj<TH2F> nPartDistr = ThrObj<TH2F>
      ("npart", "npart", 50, 0, 50, 100, Par.ptMin, Par.ptMax);

   ThrObj<TH1F> origPtDistr = ThrObj<TH1F>
      ("orig","orig", Par.ptNBins, Par.ptMin, Par.ptMax);

   ThrObj<TH2F> m2TOFeDistr = ThrObj<TH2F>
      ("m2_tofe", "m2", Par.ptNBins, Par.ptMin, Par.ptMax, 2000, -4., 4.);
   ThrObj<TH2F> m2TOFwDistr = ThrObj<TH2F>
      ("m2_tofw", "m2", Par.ptNBins, Par.ptMin, Par.ptMax, 2000, -4., 4.);

   ThrObj<TH2F> timeTOFeDistr = ThrObj<TH2F>
      ("t_tofe", "t", Par.ptNBins, Par.ptMin, Par.ptMax, 2000, -40, 40.);
   ThrObj<TH2F> timeTOFwDistr = ThrObj<TH2F>
      ("t_tofw", "t", Par.ptNBins, Par.ptMin, Par.ptMax, 2000, -40, 40.);

   ThrObj<TH2F> origPtVsRecPtDistr = ThrObj<TH2F>
      ("orig_pt_vs_pt", "pt", Par.ptNBins, Par.ptMin, Par.ptMax, 
       Par.ptNBins, Par.ptMin, Par.ptMax);

   ThrObj<TH1F> regTOFeDistr = ThrObj<TH1F>
      ("reg_tofe", "tofe", Par.ptNBins, Par.ptMin, Par.ptMax);
   ThrObj<TH1F> regTOFwDistr = ThrObj<TH1F>
      ("reg_tofw", "tofw", Par.ptNBins, Par.ptMin, Par.ptMax);
   
   std::array<ThrObj<TH1F>, 4> regEMCaleDistr = {
      ThrObj<TH1F>("reg_emcale0", "emc", Par.ptNBins, Par.ptMin, Par.ptMax),
      ThrObj<TH1F>("reg_emcale1", "emc", Par.ptNBins, Par.ptMin, Par.ptMax),
      ThrObj<TH1F>("reg_emcale2", "emc", Par.ptNBins, Par.ptMin, Par.ptMax),
      ThrObj<TH1F>("reg_emcale3", "emc", Par.ptNBins, Par.ptMin, Par.ptMax)};

   std::array<ThrObj<TH1F>, 4> regEMCalwDistr = {
      ThrObj<TH1F>("reg_emcalw0", "emc", Par.ptNBins, Par.ptMin, Par.ptMax),
      ThrObj<TH1F>("reg_emcalw1", "emc", Par.ptNBins, Par.ptMin, Par.ptMax),
      ThrObj<TH1F>("reg_emcalw2", "emc", Par.ptNBins, Par.ptMin, Par.ptMax),
      ThrObj<TH1F>("reg_emcalw3", "emc", Par.ptNBins, Par.ptMin, Par.ptMax)};
};

//auxName can be used to specify different statistics of the same dataset e.g. low pt or high pt
void AnalyzeParticle(ThrContainer *thrContainer, const std::string& particle, 
                     const std::string& magf, const std::string& auxName, const int procNum); 
void AnalyzeSingleTrack();
int main();

#endif /*ANALYZE_SINGLE_TRACK_HPP*/
