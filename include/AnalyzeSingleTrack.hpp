/** 
 *  @file   AnalyzeSingleTrack.hpp
 *  @brief  Contains declarations of functions and variables that are used for analysis of a single track from a trees acquired from the PHENIX simulation
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysisPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef ANALYZE_SINGLE_TRACK_HPP
#define ANALYZE_SINGLE_TRACK_HPP

#include <signal.h>
#include <thread>

#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"

#include "InputYAMLReader.hpp"

#include "IOTools.hpp"
#include "StrTools.hpp"
#include "Box.hpp"

#include "ThrObj.hpp"

#include "SingleTrackFunc.hpp"
#include "SimTreeReader.hpp"
#include "DeadMapCutter.hpp"

#include "PBar.hpp"

#include "ROOT/TTreeProcessorMT.hxx"

int main(int argc, char **argv);

/* @namespace AnalyzeSingleTrack
 *
 * @brief Contains all functions, variables, and containers needed for AnalyzeSingleTrack 
 *
 * This namespace is eployed so that documentation will not become a pile of variables, types, and functions from many different files that are intended to be compiled and used as executables. With this namespace finding the needed information for the given executable is easier since everything belongs to the current namespace
 */
namespace AnalyzeSingleTrack
{
   /// output directory
   std::string outputDir;
   /// collision system name
   std::string collisionSystemName;
   /// run name
   std::string runName;
   /// minimum pT of a charged track
   double pTMin;
   /// maximum pT of a charged track
   double pTMax;
   /// shows whether the  particles will be reweighted to the corresponding spectra
   bool reweightForSpectra;
   /// shows whether heatmaps will be reweighted so that the discrepancy between real data and simulation is canceled
   bool reweightHeatmapsForAlpha;
   /// shows whether the sigmalized residuals dphi and dz for all detectors are calibrated
   bool areSigmalizedResidualsCalibrated;
   /// correction for TOFw due to ADC and efficiency correction
   double correctionTOFw;
   /// file reader with all required parameters for the simulation processing
   InputYAMLReader inputYAMLSim;
   /// file reader with all required parameters for the current run
   InputYAMLReader inputYAMLMain;
   /// number of threads
   int numberOfThreads;
   /// histogram with alpha scaling for DCe, zDC>=0 that is used when reweightForAlpha is true
   TH1F *alphaReweightDCe0;
   /// histogram with alpha scaling for DCe, zDC<0 that is used when reweightForAlpha is true
   TH1F *alphaReweightDCe1;
   /// histogram with alpha scaling for DCw, zDC>=0 that is used when reweightForAlpha is true
   TH1F *alphaReweightDCw0;
   /// histogram with alpha scaling for DCw, zDC<0 that is used when reweightForAlpha is true
   TH1F *alphaReweightDCw1;
   /// number of events across all trees
   unsigned long numberOfEvents = 0;
   /// parameter for monitoring the progress
   unsigned long numberOfCalls = 0;
   /// cutter for deadmaps
   DeadMapCutter dmCutter;

   /* @brief Get the DC heatmap from the specified file
    * @param[in] file file from which the histogram will be read
    * @param[in] histName name of the histogram which will be read
    * @param[out] hist read histogram
    */
   TH2F *GetHistogramFromFile(TFile& file, const std::string& histName);
   /* @brief Compares axis of 2 histograms and if the axis are not the same prints error
    *
    * If the number of bins and axis ranges are equal for both histograms 
    * this function does nothing; else it prints an error and the program will be closed.
    * This action is required for the integrity of alpha rescale of heatmaps
    *
    * @param[in] hist1 1st histogram to be compared with the 2nd one
    * @param[in] hist2 2nd histogram to be compared with the 1st one
    * @param[out] hist read histogram
    */
   void CheckHistsAxis(const TH2F *hist1, const TH2F *hist2);

   /* @struct ThrContainerCopy
    * @brief Container for storing local ThrContainer copies (at least 1 for each thread) 
    * since copying speeds up the calculation because local copies 
    * do not require synchronization from each other
    *
    * Contains the same histograms ThrContainer but as shared_ptr<T>
    */
   struct ThrContainerCopy
   {
      /// distribution of original generated pT
      std::shared_ptr<TH1F> distrOrigPT;
      // distribution of original generated pT vs reconstructed pT in the simulation
      std::shared_ptr<TH2F> distrOrigPTVsRecPT;
      /// unscaled by alpha heatmap of DCe, zDC>=0
      std::shared_ptr<TH2F> heatmapUnscaledDCe0;
      /// unscaled by alpha heatmap of DCe, zDC<0
      std::shared_ptr<TH2F> heatmapUnscaledDCe1;
      /// unscaled by alpha heatmap of DCw, zDC>=0
      std::shared_ptr<TH2F> heatmapUnscaledDCw0;
      /// unscaled by alpha heatmap of DCw, zDC<0
      std::shared_ptr<TH2F> heatmapUnscaledDCw1;
      /// heatmap of DCe, zDC>=0
      std::shared_ptr<TH2F> heatmapDCe0;
      /// heatmap of DCe, zDC<0
      std::shared_ptr<TH2F> heatmapDCe1;
      /// heatmap of DCw, zDC>=0
      std::shared_ptr<TH2F> heatmapDCw0;
      /// heatmap of DCw, zDC<0
      std::shared_ptr<TH2F> heatmapDCw1;
      /// heatmap of PC1e
      std::shared_ptr<TH2F> heatmapPC1e;
      /// heatmap of PC1w
      std::shared_ptr<TH2F> heatmapPC1w;
      /// heatmap of PC2
      std::shared_ptr<TH2F> heatmapPC2;
      /// heatmap of PC3e
      std::shared_ptr<TH2F> heatmapPC3e;
      /// heatmap of PC3w
      std::shared_ptr<TH2F> heatmapPC3w;
      /// heatmap of TOFe
      std::shared_ptr<TH2F> heatmapTOFe;
      /// heatmap of TOFw, ptofy<100
      std::shared_ptr<TH2F> heatmapTOFw0;
      /// heatmap of TOFw, ptofy>100
      std::shared_ptr<TH2F> heatmapTOFw1;
      /// heatmaps of EMCale(0-3)
      std::array<std::shared_ptr<TH2F>, 4> heatmapEMCale;
      /// heatmaps of EMCalw(0-3)
      std::array<std::shared_ptr<TH2F>, 4> heatmapEMCalw;
      /// ecore in EMCale(0-3) vs pT distributions
      std::array<std::shared_ptr<TH2F>, 4> distrECoreVsPTEMCale;
      /// ecore in EMCalw(0-3) vs pT distributions
      std::array<std::shared_ptr<TH2F>, 4> distrECoreVsPTEMCalw;
      /// strip distribution in TOFw
      std::shared_ptr<TH1F> distrStripTOFw;
      /// slat distribution in TOFe
      std::shared_ptr<TH1F> distrSlatTOFe;
      /// eloss vs beta distribution in TOFe
      std::shared_ptr<TH2F> distrELossTOFe;
      /// pc2dphi vs pT distribution for positive tracks
      std::shared_ptr<TH2F> distrDPhiVsPTPC2Pos;
      /// pc2dz vs pT distribution for positive tracks
      std::shared_ptr<TH2F> distrDZVsPTPC2Pos;
      /// pc2dphi vs pT distribution for negative tracks
      std::shared_ptr<TH2F> distrDPhiVsPTPC2Neg;
      /// pc2dz vs pT distribution for negative tracks
      std::shared_ptr<TH2F> distrDZVsPTPC2Neg;
      /// pc3dphi vs pT distribution for positive tracks
      std::shared_ptr<TH2F> distrDPhiVsPTPC3Pos;
      /// pc3dz vs pT distribution for positive tracks
      std::shared_ptr<TH2F> distrDZVsPTPC3Pos;
      /// pc3dphi vs pT distribution for negative tracks
      std::shared_ptr<TH2F> distrDPhiVsPTPC3Neg;
      /// pc3dz vs pT distribution for negative tracks
      std::shared_ptr<TH2F> distrDZVsPTPC3Neg;
      /// tofdphi vs pT distribution for positive tracks
      std::shared_ptr<TH2F> distrDPhiVsPTTOFePos;
      /// tofdz vs pT distribution for positive tracks
      std::shared_ptr<TH2F> distrDZVsPTTOFePos;
      /// tofdphi vs pT distribution for negative tracks
      std::shared_ptr<TH2F> distrDPhiVsPTTOFeNeg;
      /// tofdz vs pT distribution for negative tracks
      std::shared_ptr<TH2F> distrDZVsPTTOFeNeg;
      /// tofwdphi vs pT distribution for positive tracks
      std::shared_ptr<TH2F> distrDPhiVsPTTOFwPos;
      /// tofwdz vs pT distribution for positive tracks
      std::shared_ptr<TH2F> distrDZVsPTTOFwPos;
      /// tofwdphi vs pT distribution for negative tracks
      std::shared_ptr<TH2F> distrDPhiVsPTTOFwNeg;
      /// tofwdz vs pT distribution for negative tracks
      std::shared_ptr<TH2F> distrDZVsPTTOFwNeg;
      /// emcdphi vs pT distributions for (0-3) sectors in east arm for positive tracks
      std::array<std::shared_ptr<TH2F>, 4> distrDPhiVsPTEMCalePos;
      /// emcdz vs pT distributions for (0-3) sectors in east arm for positive tracks
      std::array<std::shared_ptr<TH2F>, 4> distrDZVsPTEMCalePos;
      /// emcdphi vs pT distributions for (0-3) sectors in east arm for negative tracks
      std::array<std::shared_ptr<TH2F>, 4> distrDPhiVsPTEMCaleNeg;
      /// emcdz vs pT distributions for (0-3) sectors in east arm for negative tracks
      std::array<std::shared_ptr<TH2F>, 4> distrDZVsPTEMCaleNeg;
      /// emcdphi vs pT distributions for (0-3) sectors in west arm for positive tracks
      std::array<std::shared_ptr<TH2F>, 4> distrDPhiVsPTEMCalwPos;
      /// emcdz vs pT distributions for (0-3) sectors in west arm for positive tracks
      std::array<std::shared_ptr<TH2F>, 4> distrDZVsPTEMCalwPos;
      /// emcdphi vs pT distributions for (0-3) sectors in west arm for negative tracks
      std::array<std::shared_ptr<TH2F>, 4> distrDPhiVsPTEMCalwNeg;
      /// emcdz vs pT distributions for (0-3) sectors in west arm for negative tracks
      std::array<std::shared_ptr<TH2F>, 4> distrDZVsPTEMCalwNeg;
   };
   /* @struct ThrContainer
    * @brief Container for storing ROOTTools::ThrObj variables 
    * (histograms for multithreading with TTreeProcessorMT)
    */
   struct ThrContainer
   {
      /// @brief returns ThrContainerCopy with std::shared_ptr copies of ROOTTools::ThrObj histograms
      ThrContainerCopy GetCopy();
      /// distribution of original generated pT
      ROOTTools::ThrObj<TH1F> distrOrigPT = 
         ROOTTools::ThrObj<TH1F>("orig pT", "p_{T}", 100., 0., 10.);
      // distribution of original generated pT vs reconstructed pT in the simulation
      ROOTTools::ThrObj<TH2F> distrOrigPTVsRecPT = 
         ROOTTools::ThrObj<TH2F>("orig pT vs rec pT","p_{T}^{orig} vs p_{T}^{rec}", 
                                 100, 0., 10., 100, 0., 10.);
      /// unscaled by alpha heatmap of DCe, zDC>=0
      ROOTTools::ThrObj<TH2F> heatmapUnscaledDCe0 = 
         ROOTTools::ThrObj<TH2F>("Unscaled heatmap: DCe, zDC>=0", "board vs alpha", 
                                 400, 0., 80., 195, -0.39, 0.39);
      /// unscaled by alpha heatmap of DCe, zDC<0
      ROOTTools::ThrObj<TH2F> heatmapUnscaledDCe1 = 
         ROOTTools::ThrObj<TH2F>("Unscaled heatmap: DCe, zDC<0", "board vs alpha", 
                                 400, 0., 80., 195, -0.39, 0.39);
      /// unscaled by alpha heatmap of DCw, zDC>=0
      ROOTTools::ThrObj<TH2F> heatmapUnscaledDCw0 = 
         ROOTTools::ThrObj<TH2F>("Unscaled heatmap: DCw, zDC>=0", "board vs alpha", 
                                 400, 0., 80., 195, -0.39, 0.39);
      /// unscaled by alpha heatmap of DCw, zDC<0
      ROOTTools::ThrObj<TH2F> heatmapUnscaledDCw1 = 
         ROOTTools::ThrObj<TH2F>("Unscaled heatmap: DCw, zDC<0", "board vs alpha", 
                                 400, 0., 80., 195, -0.39, 0.39);
      /// heatmap of DCe, zDC>=0
      ROOTTools::ThrObj<TH2F> heatmapDCe0 = 
         ROOTTools::ThrObj<TH2F>("Heatmap: DCe,zDC>=0", "board vs alpha", 
                                 400, 0., 80., 195, -0.39, 0.39);
      /// heatmap of DCe, zDC<0
      ROOTTools::ThrObj<TH2F> heatmapDCe1 = 
         ROOTTools::ThrObj<TH2F>("Heatmap: DCe,zDC<0", "board vs alpha", 
                                 400, 0., 80., 195, -0.39, 0.39);
      /// heatmap of DCw, zDC>=0
      ROOTTools::ThrObj<TH2F> heatmapDCw0 = 
         ROOTTools::ThrObj<TH2F>("Heatmap: DCw,zDC>=0", "board vs alpha", 
                                 400, 0., 80., 195, -0.39, 0.39);
      /// heatmap of DCw, zDC<0
      ROOTTools::ThrObj<TH2F> heatmapDCw1 = 
         ROOTTools::ThrObj<TH2F>("Heatmap: DCw,zDC<0", "board vs alpha", 
                                 400, 0., 80., 195, -0.39, 0.39);
      /// heatmap of PC1e
      ROOTTools::ThrObj<TH2F> heatmapPC1e = 
         ROOTTools::ThrObj<TH2F>("Heatmap: PC1e", "pc1z vs pc1phi", 
                                 380, -95., 95., 170, 2.05, 3.75);
      /// heatmap of PC1w
      ROOTTools::ThrObj<TH2F> heatmapPC1w = 
         ROOTTools::ThrObj<TH2F>("Heatmap: PC1w", "pc1z vs pc1phi", 
                                 380, -95., 95., 165, -0.6, 1.05);
      /// heatmap of PC2
      ROOTTools::ThrObj<TH2F> heatmapPC2 = ROOTTools::ThrObj<TH2F>
         ("Heatmap: PC2", "pc3z vs pc3phi", 330, -165., 165., 165, -0.6, 1.05);
      /// heatmap of PC3e
      ROOTTools::ThrObj<TH2F> heatmapPC3e = ROOTTools::ThrObj<TH2F>
         ("Heatmap: PC3e", "pc3z vs pc3phi", 390, -195., 195., 170, 2.1, 3.8);
      /// heatmap of PC3w
      ROOTTools::ThrObj<TH2F> heatmapPC3w = ROOTTools::ThrObj<TH2F>
         ("Heatmap: PC3w", "pc3z vs pc3phi", 390, -195., 195., 170, -0.65, 1.05);
      /// heatmap of TOFe
      ROOTTools::ThrObj<TH2F> heatmapTOFe = 
         ROOTTools::ThrObj<TH2F>("Heatmap: TOFe", "ptofy vs ptofz", 
                                 185, -280., 90., 200, -200., 200.);
      /// heatmap of TOFw, ptofy<100
      ROOTTools::ThrObj<TH2F> heatmapTOFw0 = 
         ROOTTools::ThrObj<TH2F>("Heatmap: TOFw, ptofy<100", "ptofy vs ptofz", 
                                 90, -55., 35., 195, -195., 195.);
      /// heatmap of TOFw, ptofy>100
      ROOTTools::ThrObj<TH2F> heatmapTOFw1 = 
         ROOTTools::ThrObj<TH2F>("Heatmap: TOFw, ptofy>100", "ptofy vs ptofz", 
                                 85, 185., 270., 195, -195., 195.);
      /// heatmaps of EMCale(0-3)
      std::array<ROOTTools::ThrObj<TH2F>, 4> heatmapEMCale = 
      {
         ROOTTools::ThrObj<TH2F>("Heatmap: EMCale0", "ytower vs ztower", 48, 0., 48., 97, 0., 97.),
         ROOTTools::ThrObj<TH2F>("Heatmap: EMCale1", "ytower vs ztower", 48, 0., 48., 97, 0., 97.),
         ROOTTools::ThrObj<TH2F>("Heatmap: EMCale2", "ytower vs ztower", 36, 0., 36, 72, 0., 72.),
         ROOTTools::ThrObj<TH2F>("Heatmap: EMCale3", "ytower vs ztower", 36, 0., 36, 72, 0., 72.)
      };
      /// heatmaps of EMCalw(0-3)
      std::array<ROOTTools::ThrObj<TH2F>, 4> heatmapEMCalw = 
      {
         ROOTTools::ThrObj<TH2F>("Heatmap: EMCalw0", "pemcy vs pemcz", 36, 0., 36, 72, 0., 72.),
         ROOTTools::ThrObj<TH2F>("Heatmap: EMCalw1", "pemcy vs pemcz", 36, 0., 36, 72, 0., 72.),
         ROOTTools::ThrObj<TH2F>("Heatmap: EMCalw2", "pemcy vs pemcz", 36, 0., 36, 72, 0., 72.),
         ROOTTools::ThrObj<TH2F>("Heatmap: EMCalw3", "pemcy vs pemcz", 36, 0., 36, 72, 0., 72.)
      };
      /// ecore in EMCale(0-3) vs pT distributions
      std::array<ROOTTools::ThrObj<TH2F>, 4> distrECoreVsPTEMCale = 
      {
         ROOTTools::ThrObj<TH2F>("ecore vs pT: EMCale0", "E_{core} vs p_{T}", 
                                 200, 0., 10., 200, 0., 5.),
         ROOTTools::ThrObj<TH2F>("ecore vs pT: EMCale1", "E_{core} vs p_{T}", 
                                 200, 0., 10., 200, 0., 5.),
         ROOTTools::ThrObj<TH2F>("ecore vs pT: EMCale2", "E_{core} vs p_{T}", 
                                 200, 0., 10., 200, 0., 5.),
         ROOTTools::ThrObj<TH2F>("ecore vs pT: EMCale3", "E_{core} vs p_{T}", 
                                 200, 0., 10., 200, 0., 5.)
      };
      /// ecore in EMCalw(0-3) vs pT distributions
      std::array<ROOTTools::ThrObj<TH2F>, 4> distrECoreVsPTEMCalw = 
      {
         ROOTTools::ThrObj<TH2F>("ecore vs pT: EMCalw0", "E_{core} vs p_{T}", 
                                 200, 0., 10., 200, 0., 5.),
         ROOTTools::ThrObj<TH2F>("ecore vs pT: EMCalw1", "E_{core} vs p_{T}", 
                                 200, 0., 10., 200, 0., 5.),
         ROOTTools::ThrObj<TH2F>("ecore vs pT: EMCalw2", "E_{core} vs p_{T}", 
                                 200, 0., 10., 200, 0., 5.),
         ROOTTools::ThrObj<TH2F>("ecore vs pT: EMCalw3", "E_{core} vs p_{T}", 
                                 200, 0., 10., 200, 0., 5.)
      };
      /// strip distribution in TOFw
      ROOTTools::ThrObj<TH1F> distrStripTOFw = ROOTTools::ThrObj<TH1F>
         ("striptofw: TOFw", "strip number", 512, 0., 512.);
      /// slat distribution in TOFe
      ROOTTools::ThrObj<TH1F> distrSlatTOFe = ROOTTools::ThrObj<TH1F>
         ("slat: TOFe", "slat number", 960, 0., 960.);
      /// eloss vs beta distribution in TOFe
      ROOTTools::ThrObj<TH2F> distrELossTOFe = ROOTTools::ThrObj<TH2F>
         ("eloss: TOFe", "#beta vs E_{TOFe}", 100, 0., 1., 100, 0., 0.03);
      /// pc2dphi vs pT distribution for positive tracks
      ROOTTools::ThrObj<TH2F> distrDPhiVsPTPC2Pos = ROOTTools::ThrObj<TH2F>
         ("dphi vs pT: PC2, charge>0", "d#varphi_{PC2} vs p_{T}", 200, -0.15, 0.15, 100, 0., 10.);
      /// pc2dz vs pT distribution for positive tracks
      ROOTTools::ThrObj<TH2F> distrDZVsPTPC2Pos = ROOTTools::ThrObj<TH2F>
         ("2dz vs pT: PC2, charge>0", "dz_{PC2} vs p_{T}", 200, -60., 60., 100, 0., 10.);
      /// pc2dphi vs pT distribution for negative tracks
      ROOTTools::ThrObj<TH2F> distrDPhiVsPTPC2Neg = ROOTTools::ThrObj<TH2F>
         ("dphi vs pT: PC2, charge<0", "d#varphi_{PC2} vs p_{T}", 200, -0.15, 0.15, 100, 0., 10.);
      /// pc2dz vs pT distribution for negative tracks
      ROOTTools::ThrObj<TH2F> distrDZVsPTPC2Neg = ROOTTools::ThrObj<TH2F>
         ("2dz vs pT: PC2, charge<0", "dz_{PC2} vs p_{T}", 200, -60., 60., 100, 0., 10.);
      /// pc3dphi vs pT distribution for positive tracks
      ROOTTools::ThrObj<TH2F> distrDPhiVsPTPC3Pos = ROOTTools::ThrObj<TH2F>
         ("dphi vs pT: PC3, charge>0", "d#varphi_{PC3} vs p_{T}", 200, -0.15, 0.15, 100, 0., 10.);
      /// pc3dz vs pT distribution for positive tracks
      ROOTTools::ThrObj<TH2F> distrDZVsPTPC3Pos = ROOTTools::ThrObj<TH2F>
         ("dz vs pT: PC3, charge>0", "dz_{PC3} vs p_{T}", 200, -60., 60., 100, 0., 10.);
      /// pc3dphi vs pT distribution for negative tracks
      ROOTTools::ThrObj<TH2F> distrDPhiVsPTPC3Neg = ROOTTools::ThrObj<TH2F>
         ("dphi vs pT: PC3, charge<0", "d#varphi_{PC3} vs p_{T}", 200, -0.15, 0.15, 100, 0., 10.);
      /// pc3dz vs pT distribution for negative tracks
      ROOTTools::ThrObj<TH2F> distrDZVsPTPC3Neg = ROOTTools::ThrObj<TH2F>
         ("dz vs pT: PC3, charge<0", "dz_{PC3} vs p_{T}", 200, -60., 60., 100, 0., 10.);
      /// tofdphi vs pT distribution for positive tracks
      ROOTTools::ThrObj<TH2F> distrDPhiVsPTTOFePos = ROOTTools::ThrObj<TH2F>
         ("dphi vs pT: TOFe, charge>0", "d#varphi_{TOFe} vs p_{T}", 200, -0.15, 0.15, 100, 0., 10.);
      /// tofdz vs pT distribution for positive tracks
      ROOTTools::ThrObj<TH2F> distrDZVsPTTOFePos = ROOTTools::ThrObj<TH2F>
         ("dz vs pT: TOFe, charge>0", "dz_{TOFe} vs p_{T}", 200, -60., 60., 100, 0., 10.);
      /// tofdphi vs pT distribution for negative tracks
      ROOTTools::ThrObj<TH2F> distrDPhiVsPTTOFeNeg = ROOTTools::ThrObj<TH2F>
         ("dphi vs pT: TOFe, charge<0", "d#varphi_{TOFe} vs p_{T}", 200, -0.15, 0.15, 100, 0., 10.);
      /// tofdz vs pT distribution for negative tracks
      ROOTTools::ThrObj<TH2F> distrDZVsPTTOFeNeg = ROOTTools::ThrObj<TH2F>
         ("dz vs pT: TOFe, charge<0", "dz_{TOFe} vs p_{T}", 200, -60., 60., 100, 0., 10.);
      /// tofwdphi vs pT distribution for positive tracks
      ROOTTools::ThrObj<TH2F> distrDPhiVsPTTOFwPos = ROOTTools::ThrObj<TH2F>
         ("dphi vs pT: TOFw, charge>0", "d#varphi_{TOFw} vs p_{T}", 200, -0.15, 0.15, 100, 0., 10.);
      /// tofwdz vs pT distribution for positive tracks
      ROOTTools::ThrObj<TH2F> distrDZVsPTTOFwPos = ROOTTools::ThrObj<TH2F>
         ("dz vs pT: TOFw, charge>0", "dz_{TOFw} vs p_{T}", 200, -60., 60., 100, 0., 10.);
      /// tofwdphi vs pT distribution for negative tracks
      ROOTTools::ThrObj<TH2F> distrDPhiVsPTTOFwNeg = ROOTTools::ThrObj<TH2F>
         ("dphi vs pT: TOFw, charge<0", "d#varphi_{TOFw} vs p_{T}", 200, -0.15, 0.15, 100, 0., 10.);
      /// tofwdz vs pT distribution for negative tracks
      ROOTTools::ThrObj<TH2F> distrDZVsPTTOFwNeg = ROOTTools::ThrObj<TH2F>
         ("dz vs pT: TOFw, charge<0", "dz_{TOFw} vs p_{T}", 200, -60., 60., 100, 0., 10.);
      /// emcdphi vs pT distributions for (0-3) sectors in east arm for positive tracks
      std::array<ROOTTools::ThrObj<TH2F>, 4> distrDPhiVsPTEMCalePos = 
      {
         ROOTTools::ThrObj<TH2F>("dphi vs pT: EMCale0, charge>0", "d#varphi_{EMCale0} vs p_{T}", 
                                 200, -0.15, 0.15, 100, 0., 10.),
         ROOTTools::ThrObj<TH2F>("dphi vs pT: EMCale1, charge>0", "d#varphi_{EMCale1} vs p_{T}", 
                                200, -0.15, 0.15, 100, 0., 10.),
         ROOTTools::ThrObj<TH2F>("dphi vs pT: EMCale2, charge>0", "d#varphi_{EMCale2} vs p_{T}", 
                                 200, -0.15, 0.15, 100, 0., 10.),
         ROOTTools::ThrObj<TH2F>("dphi vs pT: EMCale3, charge>0", "d#varphi_{EMCale3} vs p_{T}", 
                                 200, -0.15, 0.15, 100, 0., 10.)
      };
      /// emcdz vs pT distributions for (0-3) sectors in east arm for positive tracks
      std::array<ROOTTools::ThrObj<TH2F>, 4> distrDZVsPTEMCalePos = 
      {
         ROOTTools::ThrObj<TH2F>("dz vs pT: EMCale0, charge>0", "dz_{EMCale0} vs p_{T}", 
                                 200, -60., 60., 100, 0., 10.),
         ROOTTools::ThrObj<TH2F>("dz vs pT: EMCale1, charge>0", "dz_{EMCale1} vs p_{T}", 
                                 200, -60., 60., 100, 0., 10.),
         ROOTTools::ThrObj<TH2F>("dz vs pT: EMCale2, charge>0", "dz_{EMCale2} vs p_{T}", 
                                 200, -60., 60., 100, 0., 10.),
         ROOTTools::ThrObj<TH2F>("dz vs pT: EMCale3, charge>0", "dz_{EMCale3} vs p_{T}", 
                                 200, -60., 60., 100, 0., 10.)
      };
      /// emcdphi vs pT distributions for (0-3) sectors in east arm for negative tracks
      std::array<ROOTTools::ThrObj<TH2F>, 4> distrDPhiVsPTEMCaleNeg = 
      {
         ROOTTools::ThrObj<TH2F>("dphi vs pT: EMCale0, charge<0", "d#varphi_{EMCale0} vs p_{T}", 
                                 200, -0.15, 0.15, 100, 0., 10.),
         ROOTTools::ThrObj<TH2F>("dphi vs pT: EMCale1, charge<0", "d#varphi_{EMCale1} vs p_{T}", 
                                 200, -0.15, 0.15, 100, 0., 10.),
         ROOTTools::ThrObj<TH2F>("dphi vs pT: EMCale2, charge<0", "d#varphi_{EMCale2} vs p_{T}", 
                                 200, -0.15, 0.15, 100, 0., 10.),
         ROOTTools::ThrObj<TH2F>("dphi vs pT: EMCale3, charge<0", "d#varphi_{EMCale3} vs p_{T}", 
                                 200, -0.15, 0.15, 100, 0., 10.)
      };
      /// emcdz vs pT distributions for (0-3) sectors in east arm for negative tracks
      std::array<ROOTTools::ThrObj<TH2F>, 4> distrDZVsPTEMCaleNeg = 
      {
         ROOTTools::ThrObj<TH2F>("dz vs pT: EMCale0, charge<0", "dz_{EMCale0} vs p_{T}", 
                                 200, -60., 60., 100, 0., 10.),
         ROOTTools::ThrObj<TH2F>("dz vs pT: EMCale1, charge<0", "dz_{EMCale1} vs p_{T}", 
                                 200, -60., 60., 100, 0., 10.),
         ROOTTools::ThrObj<TH2F>("dz vs pT: EMCale2, charge<0", "dz_{EMCale2} vs p_{T}", 
                                 200, -60., 60., 100, 0., 10.),
         ROOTTools::ThrObj<TH2F>("dz vs pT: EMCale3, charge<0", "dz_{EMCale3} vs p_{T}", 
                                 200, -60., 60., 100, 0., 10.)
      };
      /// emcdphi vs pT distributions for (0-3) sectors in west arm for positive tracks
      std::array<ROOTTools::ThrObj<TH2F>, 4> distrDPhiVsPTEMCalwPos = 
      {
         ROOTTools::ThrObj<TH2F>("dphi vs pT: EMCalw0, charge>0", "d#varphi_{EMCale0} vs p_{T}", 
                                 200, -0.15, 0.15, 100, 0., 10.),
         ROOTTools::ThrObj<TH2F>("dphi vs pT: EMCalw1, charge>0", "d#varphi_{EMCale1} vs p_{T}", 
                                 200, -0.15, 0.15, 100, 0., 10.),
         ROOTTools::ThrObj<TH2F>("dphi vs pT: EMCalw2, charge>0", "d#varphi_{EMCale2} vs p_{T}", 
                                 200, -0.15, 0.15, 100, 0., 10.),
         ROOTTools::ThrObj<TH2F>("dphi vs pT: EMCalw3, charge>0", "d#varphi_{EMCale3} vs p_{T}", 
                                 200, -0.15, 0.15, 100, 0., 10.)
      };
      /// emcdz vs pT distributions for (0-3) sectors in west arm for positive tracks
      std::array<ROOTTools::ThrObj<TH2F>, 4> distrDZVsPTEMCalwPos = 
      {
         ROOTTools::ThrObj<TH2F>("dz vs pT: EMCalw0, charge>0", "dz_{EMCalw0} vs p_{T}", 
                                 200, -60., 60., 100, 0., 10.),
         ROOTTools::ThrObj<TH2F>("dz vs pT: EMCalw1, charge>0", "dz_{EMCalw1} vs p_{T}", 
                                 200, -60., 60., 100, 0., 10.),
         ROOTTools::ThrObj<TH2F>("dz vs pT: EMCalw2, charge>0", "dz_{EMCalw2} vs p_{T}", 
                                 200, -60., 60., 100, 0., 10.),
         ROOTTools::ThrObj<TH2F>("dz vs pT: EMCalw3, charge>0", "dz_{EMCalw3} vs p_{T}", 
                                 200, -60., 60., 100, 0., 10.)
      };
      /// emcdphi vs pT distributions for (0-3) sectors in west arm for negative tracks
      std::array<ROOTTools::ThrObj<TH2F>, 4> distrDPhiVsPTEMCalwNeg = 
      {
         ROOTTools::ThrObj<TH2F>("dphi vs pT: EMCalw0, charge<0", "d#varphi_{EMCale0} vs p_{T}", 
                                 200, -0.15, 0.15, 100, 0., 10.),
         ROOTTools::ThrObj<TH2F>("dphi vs pT: EMCalw1, charge<0", "d#varphi_{EMCale1} vs p_{T}", 
                                 200, -0.15, 0.15, 100, 0., 10.),
         ROOTTools::ThrObj<TH2F>("dphi vs pT: EMCalw2, charge<0", "d#varphi_{EMCale2} vs p_{T}", 
                                 200, -0.15, 0.15, 100, 0., 10.),
         ROOTTools::ThrObj<TH2F>("dphi vs pT: EMCalw3, charge<0", "d#varphi_{EMCale3} vs p_{T}", 
                                 200, -0.15, 0.15, 100, 0., 10.)
      };
      /// emcdz vs pT distributions for (0-3) sectors in west arm for negative tracks
      std::array<ROOTTools::ThrObj<TH2F>, 4> distrDZVsPTEMCalwNeg = 
      {
         ROOTTools::ThrObj<TH2F>("dz vs pT: EMCalw0, charge<0", "dz_{EMCalw0} vs p_{T}", 
                                 200, -60., 60., 100, 0., 10.),
         ROOTTools::ThrObj<TH2F>("dz vs pT: EMCalw1, charge<0", "dz_{EMCalw1} vs p_{T}", 
                                 200, -60., 60., 100, 0., 10.),
         ROOTTools::ThrObj<TH2F>("dz vs pT: EMCalw2, charge<0", "dz_{EMCalw2} vs p_{T}", 
                                 200, -60., 60., 100, 0., 10.),
         ROOTTools::ThrObj<TH2F>("dz vs pT: EMCalw3, charge<0", "dz_{EMCalw3} vs p_{T}", 
                                 200, -60., 60., 100, 0., 10.)
      };
   };
   /* @brief Processes the single configuration (for the given particle, 
    * magnetic field, and pT range) from one file
    *
    * @param[in] thrContainer the current ThrContainer in which the data will be written to 
    * @param[in] particleName name of the particle to be analyzed 
    * @param[in] magneticFieldName name of the magnetic field to be analyzed 
    * @param[in] pTRangeName name of the pT range to be analyzed 
    */
   void AnalyzeConfiguration(ThrContainer& thrContainer, 
                             const std::string& particleName, const std::string& magneticFieldName, 
                             const std::string& pTRangeName); 
}

#endif /* ANALYZE_SINGLE_TRACK_HPP */
