// $HEADER$
//------------------------------------------------------------------------------------------------
//                          AnalyzeHeatMaps function delcalarion
//------------------------------------------------------------------------------------------------
// AnalyzeHeatMaps
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic macro for evaluating heat maps distributions for different detectors
 * from simulation output of event-like TTrees to processed histograms
 **/
//------------------------------------------------------------------------------------------------

#ifndef ANALYZE_HEAT_MAPS_HPP
#define ANALYZE_HEAT_MAPS_HPP

#include <signal.h>

#include "ParAnalyzeHeatMaps.hpp"

#include "IOTools.hpp"
#include "StrTools.hpp"
#include "Box.hpp"

#include "PBar.hpp"

#include "ThrObj.hpp"
#include "STrackFun.hpp"
#include "EffTreeReader.hpp"

#include "ROOT/TTreeProcessorMT.hxx"

//container for storing heat maps ThrObj histograms for different detectors
struct ThrContainer
{
   ThrObj<TH1F> orig = ThrObj<TH1F>
      ("orig","orig", Par.ptNBins, Par.ptMin, Par.ptMax);

   //original generated pT before simulation vs reconstructed pT in the simulation
   ThrObj<TH2F> origPTVsRecPT = ThrObj<TH2F>
      ("orig_pt_vs_pt","orig_pt vs pt", 
      100, 0., 10., Par.ptNBins, Par.ptMin, Par.ptMax);
      
   ThrObj<TH2F> dce0 = ThrObj<TH2F>
      ("dceast0", "map", 400, -1.5, 80.5, 200, -0.3, 0.3);
   ThrObj<TH2F> dcw0 = ThrObj<TH2F>
      ("dcwest0", "map", 400, -1.5, 80.5, 200, -0.3, 0.3);
   ThrObj<TH2F> dce1 = ThrObj<TH2F>
      ("dceast1", "map", 400, -1.5, 80.5, 200, -0.3, 0.3);
   ThrObj<TH2F> dcw1 = ThrObj<TH2F>
      ("dcwest1", "map", 400, -1.5, 80.5, 200, -0.3, 0.3);

   //unscaled alpha
   ThrObj<TH2F> dce0Unscaled = ThrObj<TH2F>
      ("unscaled_dceast0", "map", 400, -1.5, 80.5, 200, -0.3, 0.3);
   ThrObj<TH2F> dcw0Unscaled = ThrObj<TH2F>
      ("unscaled_dcwest0", "map", 400, -1.5, 80.5, 200, -0.3, 0.3);
   ThrObj<TH2F> dce1Unscaled = ThrObj<TH2F>
      ("unscaled_dceast1", "map", 400, -1.5, 80.5, 200, -0.3, 0.3);
   ThrObj<TH2F> dcw1Unscaled = ThrObj<TH2F>
      ("unscaled_dcwest1", "map", 400, -1.5, 80.5, 200, -0.3, 0.3);

   ThrObj<TH2F> zVsPhiPc1e = ThrObj<TH2F>
      ("pc1e_z_vs_phi", "z vs phi", 180, -90., 90., 330, 2.1, 3.75);
   ThrObj<TH2F> zVsPhiPc1w = ThrObj<TH2F>
      ("pc1w_z_vs_phi", "z vs phi", 180, -90., 90., 330, -0.6, 1.05);

   ThrObj<TH2F> zVsPhiPc2 = ThrObj<TH2F>
      ("pc2_z_vs_phi", "z vs phi", 160, -160., 160., 330, -0.6, 1.05);

   ThrObj<TH2F> zVsPhiPc3e = ThrObj<TH2F>
      ("pc3e_z_vs_phi", "z vs phi", 190, -190., 190., 320, 2.15, 3.75);
   ThrObj<TH2F> zVsPhiPc3w = ThrObj<TH2F>
      ("pc3w_z_vs_phi", "z vs phi", 190, -190., 190., 330, -0.6, 1.05);

   //positive tracks in EMCale
   std::array<ThrObj<TH2F>, 4> emcalePos = 
   {
      ThrObj<TH2F>("emcale0_pos", "y vs z", 100, -300., -90., 210, -30., 210.),
      ThrObj<TH2F>("emcale1_pos", "y vs z", 100, -110., 110., 210, -30., 210.),
      ThrObj<TH2F>("emcale2_pos", "y vs z", 100, 90., 300., 210, -30., 210.),
      ThrObj<TH2F>("emcale3_pos", "y vs z", 100, 280., 440., 210, -30., 210.)
   };

   std::array<ThrObj<TH2F>, 4> emcaleNeg = 
   {
      ThrObj<TH2F>("emcale0_neg", "y vs z", 100, -300., -90., 210, -210., 30.),
      ThrObj<TH2F>("emcale1_neg", "y vs z", 100, -110., 110., 210, -210., 30.),
      ThrObj<TH2F>("emcale2_neg", "y vs z", 100, 90., 300., 210, -210., 30.),
      ThrObj<TH2F>("emcale3_neg", "y vs z", 100, 280., 440., 210, -210., 30.)
   };
   
   std::array<ThrObj<TH2F>, 4> emcalwPos = 
   {
      ThrObj<TH2F>("emcalw0_pos", "y vs z", 100, -300., -90., 210, -30., 210.),
      ThrObj<TH2F>("emcalw1_pos", "y vs z", 100, -110., 110., 210, -30., 210.),
      ThrObj<TH2F>("emcalw2_pos", "y vs z", 100, 90., 300., 210, -30., 210.),
      ThrObj<TH2F>("emcalw3_pos", "y vs z", 100, 280., 440., 210, -30., 210.)
   };

   std::array<ThrObj<TH2F>, 4> emcalwNeg = 
   {
      ThrObj<TH2F>("emcalw0_neg", "y vs z", 100, -300., -90., 210, -210., 30.),
      ThrObj<TH2F>("emcalw1_neg", "y vs z", 100, -110., 110., 210, -210., 30.),
      ThrObj<TH2F>("emcalw2_neg", "y vs z", 100, 90., 300., 210, -210., 30.),
      ThrObj<TH2F>("emcalw3_neg", "y vs z", 100, 280., 440., 210, -210., 30.)
   };

   ThrObj<TH1F> stripTofw = ThrObj<TH1F>
      ("strip_tofw", "strip", 550, 0., 550.);
   ThrObj<TH1F> slatTofe = ThrObj<TH1F>
      ("slat", "slat", 1000, 0., 1000.);

   ThrObj<TH2F> tofe0 = ThrObj<TH2F>
      ("tofe0", "tofe0", 200, -300., 100., 200, -40., 210.);
   ThrObj<TH2F> tofe1 = ThrObj<TH2F>
      ("tofe1", "tofe1", 200, -300., 100., 200, -210., 40.);

   ThrObj<TH2F> tofw0 = ThrObj<TH2F>
      ("tofw0", "tofw0", 200, 15, 68, 200, -0.3, 0.3);
   ThrObj<TH2F> tofw1 = ThrObj<TH2F>
      ("tofw1", "tofw1", 200, 15, 68, 200, -0.3, 0.3);
};

//auxName can be used to specify different statistics of the same dataset e.g. low pt or high pt
void Analyze(ThrContainer *thrContainer, const std::string& part, const std::string& magf, 
             const std::string& auxName, const int procNum); 
void HeatMapper(); //for ROOT CINT call
int main(); //main calls HeatMapper - the same CINT would called; used in compiled binary

#endif /*ANALYZE_HEAT_MAPS_HPP*/
