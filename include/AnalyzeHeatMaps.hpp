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
#include <thread>

#include "TFile.h"

#include "json/json.h"
#include "json/value.h"

#include "IOTools.hpp"
#include "StrTools.hpp"
#include "Box.hpp"

#include "PBar.hpp"

#include "ThrObj.hpp"
#include "STrackFun.hpp"
#include "EffTreeReader.hpp"

#include "DataMethodsSelector.hpp"

#include "ROOT/TTreeProcessorMT.hxx"

struct Parameters
{
   const std::string outputDir = "data/PostSim/";
   
   std::string collisionSystemName;
   std::string runName;
   
   bool reweightForSpectra;
   bool reweightForAlpha;
   
   std::vector<std::string> partQueue;
   std::vector<std::string> magfQueue;
   std::vector<std::string> pTRangeQueue;

   double pTMin, pTMax;

   double correctionTOFw;

   int numberOfThreads;

   DataMethodsSelector dms;

   TH1F *alphaReweightDCe0;
   TH1F *alphaReweightDCe1;
   TH1F *alphaReweightDCw0;
   TH1F *alphaReweightDCw1;
   
   void Init(const std::string inputFileName, const int nThr);
};

//container for storing ThrObj objects that store TThreadedObject objects
struct ThrContainer
{
   ThrObj<TH1F> distrOrigPT = 
      ThrObj<TH1F>("orig pT", "p_{T}", 100., 0., 10.);

   //original generated pT vs reconstructed pT in the simulation
   ThrObj<TH2F> distrOrigPTVsRecPT = 
      ThrObj<TH2F>("orig pT vs rec pT","p_{T}^{orig} vs p_{T}^{rec}", 
                   100, 0., 10., 100, 0., 10.);
      
   ThrObj<TH2F> heatmapDCe0 = 
      ThrObj<TH2F>("Heatmap: DCe, zDC>=0", "board vs alpha", 400, 0., 80., 195, -0.39, 0.39);
   ThrObj<TH2F> heatmapDCe1 = 
      ThrObj<TH2F>("Heatmap: DCe, zDC<0", "board vs alpha", 400, 0., 80., 195, -0.39, 0.39);
   ThrObj<TH2F> heatmapDCw0 = 
      ThrObj<TH2F>("Heatmap: DCw, zDC>=0", "board vs alpha", 400, 0., 80., 195, -0.39, 0.39);
   ThrObj<TH2F> heatmapDCw1 = 
      ThrObj<TH2F>("Heatmap: DCw, zDC<0", "board vs alpha", 400, 0., 80., 195, -0.39, 0.39);

   // DC heatmaps without alpha scaling: needed for alpha scaling
   ThrObj<TH2F> heatmapUnscaledDCe0 = 
      ThrObj<TH2F>("Unscaled heatmap: DCe, zDC>=0", "board vs alpha", 
                   400, 0., 80., 195, -0.39, 0.39);
   ThrObj<TH2F> heatmapUnscaledDCe1 = 
      ThrObj<TH2F>("Unscaled heatmap: DCe, zDC<0", "board vs alpha", 
                   400, 0., 80., 195, -0.39, 0.39);
   ThrObj<TH2F> heatmapUnscaledDCw0 = 
      ThrObj<TH2F>("Unscaled heatmap: DCw, zDC>=0", "board vs alpha", 
                   400, 0., 80., 195, -0.39, 0.39);
   ThrObj<TH2F> heatmapUnscaledDCw1 = 
      ThrObj<TH2F>("Unscaled heatmap: DCw, zDC<0", "board vs alpha", 
                   400, 0., 80., 195, -0.39, 0.39);

   ThrObj<TH2F> heatmapPC1e = 
      ThrObj<TH2F>("Heatmap: PC1e", "pc1z vs pc1phi", 360, -90., 90.,170, 2.05, 3.75);
   ThrObj<TH2F> heatmapPC1w = 
      ThrObj<TH2F>("Heatmap: PC1w", "pc1z vs pc1phi", 360, -90., 90., 160, -0.6, 1.05);

   ThrObj<TH2F> heatmapPC2 = ThrObj<TH2F>
      ("Heatmap: PC2", "pc3z vs pc3phi", 320, -160., 160., 330, -0.6, 1.05);

   ThrObj<TH2F> heatmapPC3e = ThrObj<TH2F>
      ("Heatmap: PC3e", "pc3z vs pc3phi", 380, -190., 190., 170, 2.1, 3.8);
   ThrObj<TH2F> heatmapPC3w = ThrObj<TH2F>
      ("Heatmap: PC3w", "pc3z vs pc3phi", 380, -190., 190., 170, -0.65, 1.05);

   ThrObj<TH2F> heatmapTOFe = 
      ThrObj<TH2F>("Heatmap: TOFe", "ptofy vs ptofz", 185, -280., 90., 195, -195., 195.);

   ThrObj<TH2F> heatmapTOFw0 = 
      ThrObj<TH2F>("Heatmap: TOFw, ptofy<100", "ptofy vs ptofz", 90, -55., 35., 195, -195., 195.);
   ThrObj<TH2F> heatmapTOFw1 = 
      ThrObj<TH2F>("Heatmap: TOFw, ptofy>100", "ptofy vs ptofz", 85, 185., 270., 195, -195., 195.);

   std::array<ThrObj<TH2F>, 4> heatmapEMCale = 
   {
      ThrObj<TH2F>("Heatmap: EMCale0", "pemcy vs pemcz", 200, -300., -90., 205, -205., 205.),
      ThrObj<TH2F>("Heatmap: EMCale1", "pemcy vs pemcz", 200, -110., 110., 205, -205., 205.),
      ThrObj<TH2F>("Heatmap: EMCale2", "pemcy vs pemcz", 200, 90., 300., 205, -205., 205.),
      ThrObj<TH2F>("Heatmap: EMCale3", "pemcy vs pemcz", 200, 280., 440., 205, -205., 205.)
   };

   std::array<ThrObj<TH2F>, 4> heatmapEMCalw = 
   {
      ThrObj<TH2F>("Heatmap: EMCalw0", "pemcy vs pemcz", 200, -300., -90., 205, -205., 205.),
      ThrObj<TH2F>("Heatmap: EMCalw1", "pemcy vs pemcz", 200, -110., 110., 205, -205., 205.),
      ThrObj<TH2F>("Heatmap: EMCalw2", "pemcy vs pemcz", 200, 90., 300., 205, -205., 205.),
      ThrObj<TH2F>("Heatmap: EMCalw3", "pemcy vs pemcz", 200, 280., 440., 205, -205., 205.)
   };


   ThrObj<TH1F> distrStripTOFw = ThrObj<TH1F>
      ("strip: TOFw", "striptofw", 512, 0., 512.);
   ThrObj<TH1F> distrSlatTOFe = ThrObj<TH1F>
      ("slat: TOFe", "slat", 960, 0., 960.);
};

TH2F *GetDCHeatmap(TFile *file, const std::string& histName);
void CheckHistsAxis(TH2F *hist1, TH2F *hist2);

void AnalyzeConfiguration(ThrContainer *thrContainer, const std::string& part, 
                          const std::string& magf, const std::string& auxName, const int procNum); 
int main(int argc, char **argv);

#endif /*ANALYZE_HEAT_MAPS_HPP*/
