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
#include "SimCalibrator.hpp"

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
   /// file reader for all required parameters for the simulation processing
   InputYAMLReader inputYAMLSim;
   /// file reader for all required parameters for the current run
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
   /// calibrator for simulated data
   SimCalibrator simCalibrator;

   /* @brief Get the DC heatmap from the specified file
    * @param[in] file file from which the histogram will be read
    * @param[in] histName name of the histogram which will be read
    * @param[out] hist read histogram
    */
   TH2F *GetHistogramFromFile(TFile &file, const std::string& histName);
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
      /// heatmap of DCe X1, zDC>=0
      std::shared_ptr<TH2F> heatmapDCe0X1;
      /// heatmap of DCe X1, zDC<0
      std::shared_ptr<TH2F> heatmapDCe1X1;
      /// heatmap of DCw X1, zDC>=0
      std::shared_ptr<TH2F> heatmapDCw0X1;
      /// heatmap of DCw X1, zDC<0
      std::shared_ptr<TH2F> heatmapDCw1X1;
      /// heatmap of DCe X2, zDC>=0
      std::shared_ptr<TH2F> heatmapDCe0X2;
      /// heatmap of DCe X2, zDC<0
      std::shared_ptr<TH2F> heatmapDCe1X2;
      /// heatmap of DCw X2, zDC>=0
      std::shared_ptr<TH2F> heatmapDCw0X2;
      /// heatmap of DCw X2, zDC<0
      std::shared_ptr<TH2F> heatmapDCw1X2;
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
      /// heatmap of TOFw
      std::shared_ptr<TH2F> heatmapTOFw;
      /// heatmaps of EMCale(0-3)
      std::array<std::shared_ptr<TH2F>, 4> heatmapEMCale;
      /// heatmaps of EMCalw(0-3)
      std::array<std::shared_ptr<TH2F>, 4> heatmapEMCalw;
      /// heatmaps of EMCale(0-3) hits
      std::array<std::shared_ptr<TH2F>, 4> heatmapEMCaleHit;
      /// heatmaps of EMCalw(0-3) hits
      std::array<std::shared_ptr<TH2F>, 4> heatmapEMCalwHit;
      /// prob in EMCale(0-3) vs pT distributions
      std::array<std::shared_ptr<TH2F>, 4> distrProbVsPTEMCale;
      /// prob in EMCalw(0-3) vs pT distributions
      std::array<std::shared_ptr<TH2F>, 4> distrProbVsPTEMCalw;
      /// ecore in EMCale(0-3) vs pT distributions
      std::array<std::shared_ptr<TH2F>, 4> distrECoreVsPTEMCale;
      /// ecore in EMCalw(0-3) vs pT distributions
      std::array<std::shared_ptr<TH2F>, 4> distrECoreVsPTEMCalw;
      /// ecore in EMCale(0-3) vs pT distributions for original particles only
      std::array<std::shared_ptr<TH2F>, 4> distrECoreVsPTEMCaleOrig;
      /// ecore in EMCalw(0-3) vs pT distributions for original particles only
      std::array<std::shared_ptr<TH2F>, 4> distrECoreVsPTEMCalwOrig;
      /// eloss vs beta distribution in TOFe
      std::shared_ptr<TH2F> distrBetaVsETOFe;
      /// pc2dphi vs pT distribution for positive tracks
      std::shared_ptr<TH2F> distrDPhiVsPTPC2Pos;
      /// pc2dz vs pT distribution for positive tracks
      std::shared_ptr<TH2F> distrDZVsPTPC2Pos;
      /// pc2dphi vs pT distribution for negative tracks
      std::shared_ptr<TH2F> distrDPhiVsPTPC2Neg;
      /// pc2dz vs pT distribution for negative tracks
      std::shared_ptr<TH2F> distrDZVsPTPC2Neg;
      /// pc3dphi vs pT distribution for positive tracks for east arm
      std::shared_ptr<TH2F> distrDPhiVsPTPC3ePos;
      /// pc3dz vs pT distribution for positive tracks for east arm
      std::shared_ptr<TH2F> distrDZVsPTPC3ePos;
      /// pc3dphi vs pT distribution for negative tracks for east arm
      std::shared_ptr<TH2F> distrDPhiVsPTPC3eNeg;
      /// pc3dz vs pT distribution for negative tracks for east arm
      std::shared_ptr<TH2F> distrDZVsPTPC3eNeg;
      /// pc3dphi vs pT distribution for positive tracks for west arm
      std::shared_ptr<TH2F> distrDPhiVsPTPC3wPos;
      /// pc3dz vs pT distribution for positive tracks for west arm
      std::shared_ptr<TH2F> distrDZVsPTPC3wPos;
      /// pc3dphi vs pT distribution for negative tracks for west arm
      std::shared_ptr<TH2F> distrDPhiVsPTPC3wNeg;
      /// pc3dz vs pT distribution for negative tracks for west arm
      std::shared_ptr<TH2F> distrDZVsPTPC3wNeg;
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
      /// pc2dphi vs pT distribution for positive tracks
      std::shared_ptr<TH2F> distrSDPhiVsPTPC2Pos;
      /// pc2dz vs pT distribution for positive tracks
      std::shared_ptr<TH2F> distrSDZVsPTPC2Pos;
      /// pc2dphi vs pT distribution for negative tracks
      std::shared_ptr<TH2F> distrSDPhiVsPTPC2Neg;
      /// pc2dz vs pT distribution for negative tracks
      std::shared_ptr<TH2F> distrSDZVsPTPC2Neg;
      /// pc3dphi vs pT distribution for positive tracks for east arm
      std::shared_ptr<TH2F> distrSDPhiVsPTPC3ePos;
      /// pc3dz vs pT distribution for positive tracks for east arm
      std::shared_ptr<TH2F> distrSDZVsPTPC3ePos;
      /// pc3dphi vs pT distribution for negative tracks for east arm
      std::shared_ptr<TH2F> distrSDPhiVsPTPC3eNeg;
      /// pc3dz vs pT distribution for negative tracks for east arm
      std::shared_ptr<TH2F> distrSDZVsPTPC3eNeg;
      /// pc3dphi vs pT distribution for positive tracks for west arm
      std::shared_ptr<TH2F> distrSDPhiVsPTPC3wPos;
      /// pc3dz vs pT distribution for positive tracks for west arm
      std::shared_ptr<TH2F> distrSDZVsPTPC3wPos;
      /// pc3dphi vs pT distribution for negative tracks for west arm
      std::shared_ptr<TH2F> distrSDPhiVsPTPC3wNeg;
      /// pc3dz vs pT distribution for negative tracks for west arm
      std::shared_ptr<TH2F> distrSDZVsPTPC3wNeg;
      /// tofdphi vs pT distribution for positive tracks
      std::shared_ptr<TH2F> distrSDPhiVsPTTOFePos;
      /// tofdz vs pT distribution for positive tracks
      std::shared_ptr<TH2F> distrSDZVsPTTOFePos;
      /// tofdphi vs pT distribution for negative tracks
      std::shared_ptr<TH2F> distrSDPhiVsPTTOFeNeg;
      /// tofdz vs pT distribution for negative tracks
      std::shared_ptr<TH2F> distrSDZVsPTTOFeNeg;
      /// tofwdphi vs pT distribution for positive tracks
      std::shared_ptr<TH2F> distrSDPhiVsPTTOFwPos;
      /// tofwdz vs pT distribution for positive tracks
      std::shared_ptr<TH2F> distrSDZVsPTTOFwPos;
      /// tofwdphi vs pT distribution for negative tracks
      std::shared_ptr<TH2F> distrSDPhiVsPTTOFwNeg;
      /// tofwdz vs pT distribution for negative tracks
      std::shared_ptr<TH2F> distrSDZVsPTTOFwNeg;
      /// emcdphi vs pT distributions for (0-3) sectors in east arm for positive tracks
      std::array<std::shared_ptr<TH2F>, 4> distrSDPhiVsPTEMCalePos;
      /// emcdz vs pT distributions for (0-3) sectors in east arm for positive tracks
      std::array<std::shared_ptr<TH2F>, 4> distrSDZVsPTEMCalePos;
      /// emcdphi vs pT distributions for (0-3) sectors in east arm for negative tracks
      std::array<std::shared_ptr<TH2F>, 4> distrSDPhiVsPTEMCaleNeg;
      /// emcdz vs pT distributions for (0-3) sectors in east arm for negative tracks
      std::array<std::shared_ptr<TH2F>, 4> distrSDZVsPTEMCaleNeg;
      /// emcdphi vs pT distributions for (0-3) sectors in west arm for positive tracks
      std::array<std::shared_ptr<TH2F>, 4> distrSDPhiVsPTEMCalwPos;
      /// emcdz vs pT distributions for (0-3) sectors in west arm for positive tracks
      std::array<std::shared_ptr<TH2F>, 4> distrSDZVsPTEMCalwPos;
      /// emcdphi vs pT distributions for (0-3) sectors in west arm for negative tracks
      std::array<std::shared_ptr<TH2F>, 4> distrSDPhiVsPTEMCalwNeg;
      /// emcdz vs pT distributions for (0-3) sectors in west arm for negative tracks
      std::array<std::shared_ptr<TH2F>, 4> distrSDZVsPTEMCalwNeg;
      /// TOFe t-t_{exp}^{pi} distribution
      std::shared_ptr<TH2F> distrTTOFe;
      /// TOFw t-t_{exp}^{pi} distribution
      std::shared_ptr<TH2F> distrTTOFw;
      /// EMCale(0-3) t-t_{exp}^{pi} distributions
      std::array<std::shared_ptr<TH2F>, 4> distrTEMCale;
      /// EMCalw(0-3) t-t_{exp}^{pi} distributions
      std::array<std::shared_ptr<TH2F>, 4> distrTEMCalw;
      /// TOFe m2 distribution for positive tracks
      std::shared_ptr<TH2F> distrM2TOFePosCharge;
      /// TOFe m2 distribution for negative tracks
      std::shared_ptr<TH2F> distrM2TOFeNegCharge;
      /// TOFw m2 distribution for positive tracks
      std::shared_ptr<TH2F> distrM2TOFwPosCharge;
      /// TOFw m2 distribution for negative tracks
      std::shared_ptr<TH2F> distrM2TOFwNegCharge;
      /// EMCale(0-3) m2 distributions for positive tracks
      std::array<std::shared_ptr<TH2F>, 4> distrM2EMCalePosCharge;
      /// EMCale(0-3) m2 distributions for negative tracks
      std::array<std::shared_ptr<TH2F>, 4> distrM2EMCaleNegCharge;
      /// EMCalw(0-3) m2 distributions for positive tracks
      std::array<std::shared_ptr<TH2F>, 4> distrM2EMCalwPosCharge;
      /// EMCalw(0-3) m2 distributions for negative tracks
      std::array<std::shared_ptr<TH2F>, 4> distrM2EMCalwNegCharge;
   };
   /* @struct ThrContainer
    * @brief Container for storing ROOT::TThreadedObject variables 
    * (histograms for multithreading with TTreeProcessorMT)
    */
   struct ThrContainer
   {
      /// @brief returns ThrContainerCopy with std::shared_ptr 
      /// copies of ROOT::TThreadedObject histograms
      ThrContainerCopy GetCopy();
      /// @brief writes the merged histograms across all threads into the file with a specified name
      void Write(const std::string& outputFileName);
      /// distribution of original generated pT
      std::unique_ptr<ROOT::TThreadedObject<TH1F>> distrOrigPT = 
         std::make_unique<ROOT::TThreadedObject<TH1F>>("orig pT", "p_{T}", 100., 0., 10.);
      // distribution of original generated pT vs reconstructed pT in the simulation
      ROOT::TThreadedObject<TH2F> 
         distrOrigPTVsRecPT{"orig pT vs rec pT", "p_{T}^{orig} vs p_{T}^{rec}", 
                            100, 0., 10., 100, 0., 10.};
      /// unscaled by alpha heatmap of DCe, zDC>=0
      ROOT::TThreadedObject<TH2F> 
         heatmapUnscaledDCe0{"Unscaled heatmap: DCe, zDC>=0", "board vs alpha", 
                             400, 0., 80., 195, -0.39, 0.39};
      /// unscaled by alpha heatmap of DCe, zDC<0
      ROOT::TThreadedObject<TH2F> 
         heatmapUnscaledDCe1{"Unscaled heatmap: DCe, zDC<0", "board vs alpha", 
                             400, 0., 80., 195, -0.39, 0.39};
      /// unscaled by alpha heatmap of DCw, zDC>=0
      ROOT::TThreadedObject<TH2F> 
         heatmapUnscaledDCw0{"Unscaled heatmap: DCw, zDC>=0", "board vs alpha", 
                             400, 0., 80., 195, -0.39, 0.39};
      /// unscaled by alpha heatmap of DCw, zDC<0
      ROOT::TThreadedObject<TH2F> 
         heatmapUnscaledDCw1{"Unscaled heatmap: DCw, zDC<0", "board vs alpha", 
                             400, 0., 80., 195, -0.39, 0.39};
      /// heatmap of DCe, zDC>=0
      ROOT::TThreadedObject<TH2F> heatmapDCe0{"_Heatmap: DCe, zDC>=0", "board vs alpha", 
                                              400, 0., 80., 195, -0.39, 0.39};
      /// heatmap of DCe, zDC<0
      ROOT::TThreadedObject<TH2F> heatmapDCe1{"_Heatmap: DCe, zDC<0", "board vs alpha", 
                                              400, 0., 80., 195, -0.39, 0.39};
      /// heatmap of DCw, zDC>=0
      ROOT::TThreadedObject<TH2F> heatmapDCw0{"_Heatmap: DCw, zDC>=0", "board vs alpha", 
                                              400, 0., 80., 195, -0.39, 0.39};
      /// heatmap of DCw, zDC<0
      ROOT::TThreadedObject<TH2F> heatmapDCw1{"_Heatmap: DCw, zDC<0", "board vs alpha", 
                                              400, 0., 80., 195, -0.39, 0.39};
      /// heatmap of DCeX1, zDC>=0
      ROOT::TThreadedObject<TH2F> heatmapDCe0X1{"_Heatmap: DCeX1, zDC>=0", "board vs alpha", 
                                                400, 0., 80., 195, -0.39, 0.39};
      /// heatmap of DCeX1, zDC<0
      ROOT::TThreadedObject<TH2F> heatmapDCe1X1{"Heatmap: DCeX1, zDC<0", "board vs alpha", 
                                                400, 0., 80., 195, -0.39, 0.39};
      /// heatmap of DCwX1, zDC>=0
      ROOT::TThreadedObject<TH2F> heatmapDCw0X1{"Heatmap: DCwX1, zDC>=0", "board vs alpha", 
                                                400, 0., 80., 195, -0.39, 0.39};
      /// heatmap of DCwX1, zDC<0
      ROOT::TThreadedObject<TH2F> heatmapDCw1X1{"Heatmap: DCwX1, zDC<0", "board vs alpha", 
                                                400, 0., 80., 195, -0.39, 0.39};
      /// heatmap of DCeX2, zDC>=0
      ROOT::TThreadedObject<TH2F> heatmapDCe0X2{"Heatmap: DCeX2, zDC>=0", "board vs alpha", 
                                                400, 0., 80., 195, -0.39, 0.39};
      /// heatmap of DCeX2, zDC<0
      ROOT::TThreadedObject<TH2F> heatmapDCe1X2{"Heatmap: DCeX2, zDC<0", "board vs alpha", 
                                                400, 0., 80., 195, -0.39, 0.39};
      /// heatmap of DCwX2, zDC>=0
      ROOT::TThreadedObject<TH2F> heatmapDCw0X2{"Heatmap: DCwX2, zDC>=0", "board vs alpha", 
                                                400, 0., 80., 195, -0.39, 0.39};
      /// heatmap of DCwX2, zDC<0
      ROOT::TThreadedObject<TH2F> heatmapDCw1X2{"Heatmap: DCwX2, zDC<0", "board vs alpha", 
                                                400, 0., 80., 195, -0.39, 0.39};
      /// heatmap of PC1e
      ROOT::TThreadedObject<TH2F> heatmapPC1e{"Heatmap: PC1e", "pc1z vs pc1phi", 
                                              380, -95., 95., 170, 2.05, 3.75};
      /// heatmap of PC1w
      ROOT::TThreadedObject<TH2F> heatmapPC1w{"Heatmap: PC1w", "pc1z vs pc1phi", 
                                              380, -95., 95., 165, -0.6, 1.05};
      /// heatmap of PC2
      ROOT::TThreadedObject<TH2F> 
         heatmapPC2{"Heatmap: PC2", "pc3z vs pc3phi", 330, -165., 165., 165, -0.6, 1.05};
      /// heatmap of PC3e
      ROOT::TThreadedObject<TH2F> 
         heatmapPC3e{"Heatmap: PC3e", "pc3z vs pc3phi", 390, -195., 195., 170, 2.1, 3.8};
      /// heatmap of PC3w
      ROOT::TThreadedObject<TH2F> 
         heatmapPC3w{"Heatmap: PC3w", "pc3z vs pc3phi", 390, -195., 195., 170, -0.65, 1.05};
      /// heatmap of TOFe
      ROOT::TThreadedObject<TH2F> heatmapTOFe{"Heatmap: TOFe", "chamber vs slat", 
                                              10, 0., 10., 96, 0., 96.};
      /// heatmap of TOFw, ptofy<100
      ROOT::TThreadedObject<TH2F> heatmapTOFw{"Heatmap: TOFw", "chamber vs strip", 
                                              8, 0., 8., 64, 0., 64.};
      /// heatmaps of EMCale(0-3)
      std::array<ROOT::TThreadedObject<TH2F>, 4> heatmapEMCale
      {
         ROOT::TThreadedObject<TH2F>("Heatmap: EMCale0", "ytower vs ztower", 
                                     48, 0., 48., 97, 0., 97.),
         ROOT::TThreadedObject<TH2F>("Heatmap: EMCale1", "ytower vs ztower", 
                                     48, 0., 48., 97, 0., 97.),
         ROOT::TThreadedObject<TH2F>("Heatmap: EMCale2", "ytower vs ztower", 
                                     36, 0., 36, 72, 0., 72.),
         ROOT::TThreadedObject<TH2F>("Heatmap: EMCale3", "ytower vs ztower", 
                                     36, 0., 36, 72, 0., 72.)
      };
      /// heatmaps of EMCalw(0-3)
      std::array<ROOT::TThreadedObject<TH2F>, 4> heatmapEMCalw
      {
         ROOT::TThreadedObject<TH2F>("Heatmap: EMCalw0", "pemcy vs pemcz", 36, 0., 36, 72, 0., 72.),
         ROOT::TThreadedObject<TH2F>("Heatmap: EMCalw1", "pemcy vs pemcz", 36, 0., 36, 72, 0., 72.),
         ROOT::TThreadedObject<TH2F>("Heatmap: EMCalw2", "pemcy vs pemcz", 36, 0., 36, 72, 0., 72.),
         ROOT::TThreadedObject<TH2F>("Heatmap: EMCalw3", "pemcy vs pemcz", 36, 0., 36, 72, 0., 72.)
      };
      /// heatmaps of EMCale(0-3) hits
      std::array<ROOT::TThreadedObject<TH2F>, 4> heatmapEMCaleHit
      {
         ROOT::TThreadedObject<TH2F>("Heatmap: EMCale0 hit", "ytower vs ztower", 
                                     48, 0., 48., 97, 0., 97.),
         ROOT::TThreadedObject<TH2F>("Heatmap: EMCale1 hit", "ytower vs ztower", 
                                     48, 0., 48., 97, 0., 97.),
         ROOT::TThreadedObject<TH2F>("Heatmap: EMCale2 hit", "ytower vs ztower", 
                                     36, 0., 36, 72, 0., 72.),
         ROOT::TThreadedObject<TH2F>("Heatmap: EMCale3 hit", "ytower vs ztower", 
                                     36, 0., 36, 72, 0., 72.)
      };
      /// heatmaps of EMCalw(0-3) hits
      std::array<ROOT::TThreadedObject<TH2F>, 4> heatmapEMCalwHit
      {
         ROOT::TThreadedObject<TH2F>("Heatmap: EMCalw0 hit", "pemcy vs pemcz", 
                                     36, 0., 36, 72, 0., 72.),
         ROOT::TThreadedObject<TH2F>("Heatmap: EMCalw1 hit", "pemcy vs pemcz", 
                                     36, 0., 36, 72, 0., 72.),
         ROOT::TThreadedObject<TH2F>("Heatmap: EMCalw2 hit", "pemcy vs pemcz", 
                                     36, 0., 36, 72, 0., 72.),
         ROOT::TThreadedObject<TH2F>("Heatmap: EMCalw3 hit", "pemcy vs pemcz", 
                                     36, 0., 36, 72, 0., 72.)
      };
      /// prob in EMCale(0-3) vs pT distributions
      std::array<ROOT::TThreadedObject<TH2F>, 4> distrProbVsPTEMCale
      {
         ROOT::TThreadedObject<TH2F>("prob vs pT, EMCale0", "prob vs p_{T}", 
                                     200, 0., 10., 200, 0., 1.),
         ROOT::TThreadedObject<TH2F>("prob vs pT, EMCale1", "prob vs p_{T}", 
                                     200, 0., 10., 200, 0., 1.),
         ROOT::TThreadedObject<TH2F>("prob vs pT, EMCale2", "prob vs p_{T}", 
                                     200, 0., 10., 200, 0., 1.),
         ROOT::TThreadedObject<TH2F>("prob vs pT, EMCale3", "prob vs p_{T}", 
                                     200, 0., 10., 200, 0., 1.)
      };
      /// prob in EMCalw(0-3) vs pT distributions
      std::array<ROOT::TThreadedObject<TH2F>, 4> distrProbVsPTEMCalw
      {
         ROOT::TThreadedObject<TH2F>("prob vs pT, EMCalw0", "prob vs p_{T}", 
                                     200, 0., 10., 200, 0., 1.),
         ROOT::TThreadedObject<TH2F>("prob vs pT, EMCalw1", "prob vs p_{T}", 
                                     200, 0., 10., 200, 0., 1.),
         ROOT::TThreadedObject<TH2F>("prob vs pT, EMCalw2", "prob vs p_{T}", 
                                     200, 0., 10., 200, 0., 1.),
         ROOT::TThreadedObject<TH2F>("prob vs pT, EMCalw3", "prob vs p_{T}", 
                                     200, 0., 10., 200, 0., 1.)
      };
      /// ecore in EMCale(0-3) vs pT distributions
      std::array<ROOT::TThreadedObject<TH2F>, 4> distrECoreVsPTEMCale
      {
         ROOT::TThreadedObject<TH2F>("ecore vs pT, EMCale0", "E_{core} vs p_{T}", 
                                     200, 0., 10., 200, 0., 2.),
         ROOT::TThreadedObject<TH2F>("ecore vs pT, EMCale1", "E_{core} vs p_{T}", 
                                     200, 0., 10., 200, 0., 2.),
         ROOT::TThreadedObject<TH2F>("ecore vs pT, EMCale2", "E_{core} vs p_{T}", 
                                     200, 0., 10., 200, 0., 2.),
         ROOT::TThreadedObject<TH2F>("ecore vs pT, EMCale3", "E_{core} vs p_{T}", 
                                     200, 0., 10., 200, 0., 2.)
      };
      /// ecore in EMCalw(0-3) vs pT distributions
      std::array<ROOT::TThreadedObject<TH2F>, 4> distrECoreVsPTEMCalw
      {
         ROOT::TThreadedObject<TH2F>("ecore vs pT, EMCalw0", "E_{core} vs p_{T}", 
                                     200, 0., 10., 200, 0., 2.),
         ROOT::TThreadedObject<TH2F>("ecore vs pT, EMCalw1", "E_{core} vs p_{T}", 
                                     200, 0., 10., 200, 0., 2.),
         ROOT::TThreadedObject<TH2F>("ecore vs pT, EMCalw2", "E_{core} vs p_{T}", 
                                     200, 0., 10., 200, 0., 2.),
         ROOT::TThreadedObject<TH2F>("ecore vs pT, EMCalw3", "E_{core} vs p_{T}", 
                                     200, 0., 10., 200, 0., 2.)
      };
      /// ecore in EMCale(0-3) vs pT distributions for original particles only
      std::array<ROOT::TThreadedObject<TH2F>, 4> distrECoreVsPTEMCaleOrig
      {
         ROOT::TThreadedObject<TH2F>("ecore vs pT, EMCale0, orig only", "E_{core} vs p_{T}", 
                                     200, 0., 10., 200, 0., 2.),
         ROOT::TThreadedObject<TH2F>("ecore vs pT, EMCale1, orig only", "E_{core} vs p_{T}", 
                                     200, 0., 10., 200, 0., 2.),
         ROOT::TThreadedObject<TH2F>("ecore vs pT, EMCale2, orig only", "E_{core} vs p_{T}", 
                                     200, 0., 10., 200, 0., 2.),
         ROOT::TThreadedObject<TH2F>("ecore vs pT, EMCale3, orig only", "E_{core} vs p_{T}", 
                                     200, 0., 10., 200, 0., 2.)
      };
      /// ecore in EMCalw(0-3) vs pT distributions
      std::array<ROOT::TThreadedObject<TH2F>, 4> distrECoreVsPTEMCalwOrig
      {
         ROOT::TThreadedObject<TH2F>("ecore vs pT, EMCalw0, orig only", "E_{core} vs p_{T}", 
                                     200, 0., 10., 200, 0., 2.),
         ROOT::TThreadedObject<TH2F>("ecore vs pT, EMCalw1, orig only", "E_{core} vs p_{T}", 
                                     200, 0., 10., 200, 0., 2.),
         ROOT::TThreadedObject<TH2F>("ecore vs pT, EMCalw2, orig only", "E_{core} vs p_{T}", 
                                     200, 0., 10., 200, 0., 2.),
         ROOT::TThreadedObject<TH2F>("ecore vs pT, EMCalw3, orig only", "E_{core} vs p_{T}", 
                                     200, 0., 10., 200, 0., 2.)
      };
      /// eloss vs beta distribution in TOFe
      ROOT::TThreadedObject<TH2F> distrBetaVsETOFe{"beta vs E, TOFe", "#beta vs E_{TOFe}", 
                                                   100, 0., 1., 100, 0., 0.03};
      /// pc2dphi vs pT distribution for positive tracks
      ROOT::TThreadedObject<TH2F> distrDPhiVsPTPC2Pos
         {"dphi vs pT, PC2, charge>0", "d#varphi_{PC2} vs p_{T}", 200, -0.1, 0.1, 100, 0., 10.};
      /// pc2dz vs pT distribution for positive tracks
      ROOT::TThreadedObject<TH2F> distrDZVsPTPC2Pos
         {"dz vs pT, PC2, charge>0", "dz_{PC2} vs p_{T}", 200, -50., 50., 100, 0., 10.};
      /// pc2dphi vs pT distribution for negative tracks
      ROOT::TThreadedObject<TH2F> distrDPhiVsPTPC2Neg
         {"dphi vs pT, PC2, charge<0", "d#varphi_{PC2} vs p_{T}", 200, -0.1, 0.1, 100, 0., 10.};
      /// pc2dz vs pT distribution for negative tracks
      ROOT::TThreadedObject<TH2F> distrDZVsPTPC2Neg
         {"dz vs pT, PC2, charge<0", "dz_{PC2} vs p_{T}", 200, -50., 50., 100, 0., 10.};
      /// pc3dphi vs pT distribution for positive tracks for east arm
      ROOT::TThreadedObject<TH2F> distrDPhiVsPTPC3ePos
         {"dphi vs pT, PC3e, charge>0", "d#varphi_{PC3e} vs p_{T}", 200, -0.1, 0.1, 100, 0., 10.};
      /// pc3dz vs pT distribution for positive tracks for east arm
      ROOT::TThreadedObject<TH2F> distrDZVsPTPC3ePos
         {"dz vs pT, PC3e, charge>0", "dz_{PC3e} vs p_{T}", 200, -50., 50., 100, 0., 10.};
      /// pc3dphi vs pT distribution for negative tracks for east arm
      ROOT::TThreadedObject<TH2F> distrDPhiVsPTPC3eNeg
         {"dphi vs pT, PC3e, charge<0", "d#varphi_{PC3e} vs p_{T}", 200, -0.1, 0.1, 100, 0., 10.};
      /// pc3dz vs pT distribution for negative tracks for east arm
      ROOT::TThreadedObject<TH2F> distrDZVsPTPC3eNeg
         {"dz vs pT, PC3e, charge<0", "dz_{PC3e} vs p_{T}", 200, -50., 50., 100, 0., 10.};
      /// pc3dphi vs pT distribution for positive tracks for west arm
      ROOT::TThreadedObject<TH2F> distrDPhiVsPTPC3wPos
         {"dphi vs pT, PC3w, charge>0", "d#varphi_{PC3w} vs p_{T}", 200, -0.1, 0.1, 100, 0., 10.};
      /// pc3dz vs pT distribution for positive tracks for west arm
      ROOT::TThreadedObject<TH2F> distrDZVsPTPC3wPos
         {"dz vs pT, PC3w, charge>0", "dz_{PC3w} vs p_{T}", 200, -50., 50., 100, 0., 10.};
      /// pc3dphi vs pT distribution for negative tracks for west arm
      ROOT::TThreadedObject<TH2F> distrDPhiVsPTPC3wNeg
         {"dphi vs pT, PC3w, charge<0", "d#varphi_{PC3w} vs p_{T}", 200, -0.1, 0.1, 100, 0., 10.};
      /// pc3dz vs pT distribution for negative tracks for west arm
      ROOT::TThreadedObject<TH2F> distrDZVsPTPC3wNeg
         {"dz vs pT, PC3w, charge<0", "dz_{PC3w} vs p_{T}", 200, -50., 50., 100, 0., 10.};
      /// tofdphi vs pT distribution for positive tracks
      ROOT::TThreadedObject<TH2F> distrDPhiVsPTTOFePos
         {"dphi vs pT, TOFe, charge>0", "d#varphi_{TOFe} vs p_{T}", 200, -0.1, 0.1, 100, 0., 10.};
      /// tofdz vs pT distribution for positive tracks
      ROOT::TThreadedObject<TH2F> distrDZVsPTTOFePos
         {"dz vs pT, TOFe, charge>0", "dz_{TOFe} vs p_{T}", 200, -50., 50., 100, 0., 10.};
      /// tofdphi vs pT distribution for negative tracks
      ROOT::TThreadedObject<TH2F> distrDPhiVsPTTOFeNeg
         {"dphi vs pT, TOFe, charge<0", "d#varphi_{TOFe} vs p_{T}", 200, -0.1, 0.1, 100, 0., 10.};
      /// tofdz vs pT distribution for negative tracks
      ROOT::TThreadedObject<TH2F> distrDZVsPTTOFeNeg
         {"dz vs pT, TOFe, charge<0", "dz_{TOFe} vs p_{T}", 200, -50., 50., 100, 0., 10.};
      /// tofwdphi vs pT distribution for positive tracks
      ROOT::TThreadedObject<TH2F> distrDPhiVsPTTOFwPos
         {"dphi vs pT, TOFw, charge>0", "d#varphi_{TOFw} vs p_{T}", 200, -0.1, 0.1, 100, 0., 10.};
      /// tofwdz vs pT distribution for positive tracks
      ROOT::TThreadedObject<TH2F> distrDZVsPTTOFwPos
         {"dz vs pT, TOFw, charge>0", "dz_{TOFw} vs p_{T}", 200, -50., 50., 100, 0., 10.};
      /// tofwdphi vs pT distribution for negative tracks
      ROOT::TThreadedObject<TH2F> distrDPhiVsPTTOFwNeg
         {"dphi vs pT, TOFw, charge<0", "d#varphi_{TOFw} vs p_{T}", 200, -0.1, 0.1, 100, 0., 10.};
      /// tofwdz vs pT distribution for negative tracks
      ROOT::TThreadedObject<TH2F> distrDZVsPTTOFwNeg
         {"dz vs pT, TOFw, charge<0", "dz_{TOFw} vs p_{T}", 200, -50., 50., 100, 0., 10.};
      /// emcdphi vs pT distributions for (0-3) sectors in east arm for positive tracks
      std::array<ROOT::TThreadedObject<TH2F>, 4> distrDPhiVsPTEMCalePos
      {
         ROOT::TThreadedObject<TH2F>("dphi vs pT, EMCale0, charge>0", "d#varphi_{EMCale0} vs p_{T}", 
                                     200, -0.1, 0.1, 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("dphi vs pT, EMCale1, charge>0", "d#varphi_{EMCale1} vs p_{T}", 
                                     200, -0.1, 0.1, 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("dphi vs pT, EMCale2, charge>0", "d#varphi_{EMCale2} vs p_{T}", 
                                     200, -0.1, 0.1, 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("dphi vs pT, EMCale3, charge>0", "d#varphi_{EMCale3} vs p_{T}", 
                                     200, -0.1, 0.1, 100, 0., 10.)
      };
      /// emcdz vs pT distributions for (0-3) sectors in east arm for positive tracks
      std::array<ROOT::TThreadedObject<TH2F>, 4> distrDZVsPTEMCalePos
      {
         ROOT::TThreadedObject<TH2F>("dz vs pT, EMCale0, charge>0", "dz_{EMCale0} vs p_{T}", 
                                     200, -50., 50., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("dz vs pT, EMCale1, charge>0", "dz_{EMCale1} vs p_{T}", 
                                     200, -50., 50., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("dz vs pT, EMCale2, charge>0", "dz_{EMCale2} vs p_{T}", 
                                     200, -50., 50., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("dz vs pT, EMCale3, charge>0", "dz_{EMCale3} vs p_{T}", 
                                     200, -50., 50., 100, 0., 10.)
      };
      /// emcdphi vs pT distributions for (0-3) sectors in east arm for negative tracks
      std::array<ROOT::TThreadedObject<TH2F>, 4> distrDPhiVsPTEMCaleNeg
      {
         ROOT::TThreadedObject<TH2F>("dphi vs pT, EMCale0, charge<0", "d#varphi_{EMCale0} vs p_{T}", 
                                     200, -0.1, 0.1, 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("dphi vs pT, EMCale1, charge<0", "d#varphi_{EMCale1} vs p_{T}", 
                                     200, -0.1, 0.1, 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("dphi vs pT, EMCale2, charge<0", "d#varphi_{EMCale2} vs p_{T}", 
                                     200, -0.1, 0.1, 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("dphi vs pT, EMCale3, charge<0", "d#varphi_{EMCale3} vs p_{T}", 
                                     200, -0.1, 0.1, 100, 0., 10.)
      };
      /// emcdz vs pT distributions for (0-3) sectors in east arm for negative tracks
      std::array<ROOT::TThreadedObject<TH2F>, 4> distrDZVsPTEMCaleNeg
      {
         ROOT::TThreadedObject<TH2F>("dz vs pT, EMCale0, charge<0", "dz_{EMCale0} vs p_{T}", 
                                     200, -50., 50., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("dz vs pT, EMCale1, charge<0", "dz_{EMCale1} vs p_{T}", 
                                     200, -50., 50., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("dz vs pT, EMCale2, charge<0", "dz_{EMCale2} vs p_{T}", 
                                     200, -50., 50., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("dz vs pT, EMCale3, charge<0", "dz_{EMCale3} vs p_{T}", 
                                     200, -50., 50., 100, 0., 10.)
      };
      /// emcdphi vs pT distributions for (0-3) sectors in west arm for positive tracks
      std::array<ROOT::TThreadedObject<TH2F>, 4> distrDPhiVsPTEMCalwPos
      {
         ROOT::TThreadedObject<TH2F>("dphi vs pT, EMCalw0, charge>0", "d#varphi_{EMCale0} vs p_{T}", 
                                     200, -0.1, 0.1, 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("dphi vs pT, EMCalw1, charge>0", "d#varphi_{EMCale1} vs p_{T}", 
                                     200, -0.1, 0.1, 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("dphi vs pT, EMCalw2, charge>0", "d#varphi_{EMCale2} vs p_{T}", 
                                     200, -0.1, 0.1, 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("dphi vs pT, EMCalw3, charge>0", "d#varphi_{EMCale3} vs p_{T}", 
                                     200, -0.1, 0.1, 100, 0., 10.)
      };
      /// emcdz vs pT distributions for (0-3) sectors in west arm for positive tracks
      std::array<ROOT::TThreadedObject<TH2F>, 4> distrDZVsPTEMCalwPos
      {
         ROOT::TThreadedObject<TH2F>("dz vs pT, EMCalw0, charge>0", "dz_{EMCalw0} vs p_{T}", 
                                     200, -50., 50., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("dz vs pT, EMCalw1, charge>0", "dz_{EMCalw1} vs p_{T}", 
                                     200, -50., 50., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("dz vs pT, EMCalw2, charge>0", "dz_{EMCalw2} vs p_{T}", 
                                     200, -50., 50., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("dz vs pT, EMCalw3, charge>0", "dz_{EMCalw3} vs p_{T}", 
                                     200, -50., 50., 100, 0., 10.)
      };
      /// emcdphi vs pT distributions for (0-3) sectors in west arm for negative tracks
      std::array<ROOT::TThreadedObject<TH2F>, 4> distrDPhiVsPTEMCalwNeg
      {
         ROOT::TThreadedObject<TH2F>("dphi vs pT, EMCalw0, charge<0", "d#varphi_{EMCale0} vs p_{T}", 
                                     200, -0.1, 0.1, 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("dphi vs pT, EMCalw1, charge<0", "d#varphi_{EMCale1} vs p_{T}", 
                                     200, -0.1, 0.1, 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("dphi vs pT, EMCalw2, charge<0", "d#varphi_{EMCale2} vs p_{T}", 
                                     200, -0.1, 0.1, 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("dphi vs pT, EMCalw3, charge<0", "d#varphi_{EMCale3} vs p_{T}", 
                                     200, -0.1, 0.1, 100, 0., 10.)
      };
      /// emcdz vs pT distributions for (0-3) sectors in west arm for negative tracks
      std::array<ROOT::TThreadedObject<TH2F>, 4> distrDZVsPTEMCalwNeg
      {
         ROOT::TThreadedObject<TH2F>("dz vs pT, EMCalw0, charge<0", "dz_{EMCalw0} vs p_{T}", 
                                     200, -50., 50., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("dz vs pT, EMCalw1, charge<0", "dz_{EMCalw1} vs p_{T}", 
                                     200, -50., 50., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("dz vs pT, EMCalw2, charge<0", "dz_{EMCalw2} vs p_{T}", 
                                     200, -50., 50., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("dz vs pT, EMCalw3, charge<0", "dz_{EMCalw3} vs p_{T}", 
                                     200, -50., 50., 100, 0., 10.)
      };
      /// pc2sdphi vs pT distribution for positive tracks
      ROOT::TThreadedObject<TH2F> distrSDPhiVsPTPC2Pos
         {"sdphi vs pT, PC2, charge>0", "sd#varphi_{PC2} vs p_{T}", 200, -5., 5., 100, 0., 10.};
      /// pc2sdz vs pT distribution for positive tracks
      ROOT::TThreadedObject<TH2F> distrSDZVsPTPC2Pos
         {"sdz vs pT, PC2, charge>0", "sdz_{PC2} vs p_{T}", 200, -5., 5., 100, 0., 10.};
      /// pc2sdphi vs pT distribution for negative tracks
      ROOT::TThreadedObject<TH2F> distrSDPhiVsPTPC2Neg
         {"sdphi vs pT, PC2, charge<0", "sd#varphi_{PC2} vs p_{T}", 200, -5., 5., 100, 0., 10.};
      /// pc2sdz vs pT distribution for negative tracks
      ROOT::TThreadedObject<TH2F> distrSDZVsPTPC2Neg
         {"sdz vs pT, PC2, charge<0", "sdz_{PC2} vs p_{T}", 200, -5., 5., 100, 0., 10.};
      /// pc3sdphi vs pT distribution for positive tracks for east arm
      ROOT::TThreadedObject<TH2F> distrSDPhiVsPTPC3ePos
         {"sdphi vs pT, PC3e, charge>0", "sd#varphi_{PC3e} vs p_{T}", 200, -5., 5., 100, 0., 10.};
      /// pc3sdz vs pT distribution for positive tracks for east arm
      ROOT::TThreadedObject<TH2F> distrSDZVsPTPC3ePos
         {"sdz vs pT, PC3e, charge>0", "sdz_{PC3e} vs p_{T}", 200, -5., 5., 100, 0., 10.};
      /// pc3sdphi vs pT distribution for negative tracks for east arm
      ROOT::TThreadedObject<TH2F> distrSDPhiVsPTPC3eNeg
         {"sdphi vs pT, PC3e, charge<0", "sd#varphi_{PC3e} vs p_{T}", 200, -5., 5., 100, 0., 10.};
      /// pc3sdz vs pT distribution for negative tracks for east arm
      ROOT::TThreadedObject<TH2F> distrSDZVsPTPC3eNeg
         {"sdz vs pT, PC3e, charge<0", "sdz_{PC3e} vs p_{T}", 200, -5., 5., 100, 0., 10.};
      /// pc3sdphi vs pT distribution for positive tracks for west arm
      ROOT::TThreadedObject<TH2F> distrSDPhiVsPTPC3wPos
         {"sdphi vs pT, PC3w, charge>0", "sd#varphi_{PC3w} vs p_{T}", 200, -5., 5., 100, 0., 10.};
      /// pc3sdz vs pT distribution for positive tracks for west arm
      ROOT::TThreadedObject<TH2F> distrSDZVsPTPC3wPos
         {"sdz vs pT, PC3w, charge>0", "sdz_{PC3w} vs p_{T}", 200, -5., 5., 100, 0., 10.};
      /// pc3sdphi vs pT distribution for negative tracks for west arm
      ROOT::TThreadedObject<TH2F> distrSDPhiVsPTPC3wNeg
         {"sdphi vs pT, PC3w, charge<0", "sd#varphi_{PC3w} vs p_{T}", 200, -5., 5., 100, 0., 10.};
      /// pc3sdz vs pT distribution for negative tracks for west arm
      ROOT::TThreadedObject<TH2F> distrSDZVsPTPC3wNeg
         {"sdz vs pT, PC3w, charge<0", "sdz_{PC3w} vs p_{T}", 200, -5., 5., 100, 0., 10.};
      /// tofsdphi vs pT distribution for positive tracks
      ROOT::TThreadedObject<TH2F> distrSDPhiVsPTTOFePos
         {"sdphi vs pT, TOFe, charge>0", "sd#varphi_{TOFe} vs p_{T}", 200, -5., 5., 100, 0., 10.};
      /// tofsdz vs pT distribution for positive tracks
      ROOT::TThreadedObject<TH2F> distrSDZVsPTTOFePos
         {"sdz vs pT, TOFe, charge>0", "sdz_{TOFe} vs p_{T}", 200, -5., 5., 100, 0., 10.};
      /// tofsdphi vs pT distribution for negative tracks
      ROOT::TThreadedObject<TH2F> distrSDPhiVsPTTOFeNeg
         {"sdphi vs pT, TOFe, charge<0", "sd#varphi_{TOFe} vs p_{T}", 200, -5., 5., 100, 0., 10.};
      /// tofsdz vs pT distribution for negative tracks
      ROOT::TThreadedObject<TH2F> distrSDZVsPTTOFeNeg
         {"sdz vs pT, TOFe, charge<0", "sdz_{TOFe} vs p_{T}", 200, -5., 5., 100, 0., 10.};
      /// tofwsdphi vs pT distribution for positive tracks
      ROOT::TThreadedObject<TH2F> distrSDPhiVsPTTOFwPos
         {"sdphi vs pT, TOFw, charge>0", "sd#varphi_{TOFw} vs p_{T}", 200, -5., 5., 100, 0., 10.};
      /// tofwsdz vs pT distribution for positive tracks
      ROOT::TThreadedObject<TH2F> distrSDZVsPTTOFwPos
         {"sdz vs pT, TOFw, charge>0", "sdz_{TOFw} vs p_{T}", 200, -5., 5., 100, 0., 10.};
      /// tofwsdphi vs pT distribution for negative tracks
      ROOT::TThreadedObject<TH2F> distrSDPhiVsPTTOFwNeg
         {"sdphi vs pT, TOFw, charge<0", "sd#varphi_{TOFw} vs p_{T}", 200, -5., 5., 100, 0., 10.};
      /// tofwsdz vs pT distribution for negative tracks
      ROOT::TThreadedObject<TH2F> distrSDZVsPTTOFwNeg
         {"sdz vs pT, TOFw, charge<0", "sdz_{TOFw} vs p_{T}", 200, -5., 5., 100, 0., 10.};
      /// emcsdphi vs pT distributions for (0-3) sectors in east arm for positive tracks
      std::array<ROOT::TThreadedObject<TH2F>, 4> distrSDPhiVsPTEMCalePos
      {
         ROOT::TThreadedObject<TH2F>("sdphi vs pT, EMCale0, charge>0", 
                                     "sd#varphi_{EMCale0} vs p_{T}", 200, -5., 5., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("sdphi vs pT, EMCale1, charge>0", 
                                     "sd#varphi_{EMCale1} vs p_{T}", 200, -5., 5., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("sdphi vs pT, EMCale2, charge>0", 
                                     "sd#varphi_{EMCale2} vs p_{T}", 200, -5., 5., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("sdphi vs pT, EMCale3, charge>0", 
                                     "sd#varphi_{EMCale3} vs p_{T}", 200, -5., 5., 100, 0., 10.)
      };
      /// emcsdz vs pT distributions for (0-3) sectors in east arm for positive tracks
      std::array<ROOT::TThreadedObject<TH2F>, 4> distrSDZVsPTEMCalePos
      {
         ROOT::TThreadedObject<TH2F>("sdz vs pT, EMCale0, charge>0", "sdz_{EMCale0} vs p_{T}", 
                                     200, -5., 5., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("sdz vs pT, EMCale1, charge>0", "sdz_{EMCale1} vs p_{T}", 
                                     200, -5., 5., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("sdz vs pT, EMCale2, charge>0", "sdz_{EMCale2} vs p_{T}", 
                                     200, -5., 5., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("sdz vs pT, EMCale3, charge>0", "sdz_{EMCale3} vs p_{T}", 
                                     200, -5., 5., 100, 0., 10.)
      };
      /// emcsdphi vs pT distributions for (0-3) sectors in east arm for negative tracks
      std::array<ROOT::TThreadedObject<TH2F>, 4> distrSDPhiVsPTEMCaleNeg
      {
         ROOT::TThreadedObject<TH2F>("sdphi vs pT, EMCale0, charge<0", 
                                     "sd#varphi_{EMCale0} vs p_{T}", 200, -5., 5., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("sdphi vs pT, EMCale1, charge<0", 
                                     "sd#varphi_{EMCale1} vs p_{T}", 200, -5., 5., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("sdphi vs pT, EMCale2, charge<0", 
                                     "sd#varphi_{EMCale2} vs p_{T}", 200, -5., 5., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("sdphi vs pT, EMCale3, charge<0", 
                                     "sd#varphi_{EMCale3} vs p_{T}", 200, -5., 5., 100, 0., 10.)
      };
      /// emcsdz vs pT distributions for (0-3) sectors in east arm for negative tracks
      std::array<ROOT::TThreadedObject<TH2F>, 4> distrSDZVsPTEMCaleNeg
      {
         ROOT::TThreadedObject<TH2F>("sdz vs pT, EMCale0, charge<0", "sdz_{EMCale0} vs p_{T}", 
                                     200, -5., 5., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("sdz vs pT, EMCale1, charge<0", "sdz_{EMCale1} vs p_{T}", 
                                     200, -5., 5., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("sdz vs pT, EMCale2, charge<0", "sdz_{EMCale2} vs p_{T}", 
                                     200, -5., 5., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("sdz vs pT, EMCale3, charge<0", "sdz_{EMCale3} vs p_{T}", 
                                     200, -5., 5., 100, 0., 10.)
      };
      /// emcsdphi vs pT distributions for (0-3) sectors in west arm for positive tracks
      std::array<ROOT::TThreadedObject<TH2F>, 4> distrSDPhiVsPTEMCalwPos
      {
         ROOT::TThreadedObject<TH2F>("sdphi vs pT, EMCalw0, charge>0", 
                                     "sd#varphi_{EMCale0} vs p_{T}", 200, -5., 5., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("sdphi vs pT, EMCalw1, charge>0", 
                                     "sd#varphi_{EMCale1} vs p_{T}", 200, -5., 5., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("sdphi vs pT, EMCalw2, charge>0", 
                                     "sd#varphi_{EMCale2} vs p_{T}", 200, -5., 5., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("sdphi vs pT, EMCalw3, charge>0", 
                                     "sd#varphi_{EMCale3} vs p_{T}", 200, -5., 5., 100, 0., 10.)
      };
      /// emcsdz vs pT distributions for (0-3) sectors in west arm for positive tracks
      std::array<ROOT::TThreadedObject<TH2F>, 4> distrSDZVsPTEMCalwPos
      {
         ROOT::TThreadedObject<TH2F>("sdz vs pT, EMCalw0, charge>0", "sdz_{EMCalw0} vs p_{T}", 
                                     200, -5., 5., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("sdz vs pT, EMCalw1, charge>0", "sdz_{EMCalw1} vs p_{T}", 
                                     200, -5., 5., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("sdz vs pT, EMCalw2, charge>0", "sdz_{EMCalw2} vs p_{T}", 
                                     200, -5., 5., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("sdz vs pT, EMCalw3, charge>0", "sdz_{EMCalw3} vs p_{T}", 
                                     200, -5., 5., 100, 0., 10.)
      };
      /// emcsdphi vs pT distributions for (0-3) sectors in west arm for negative tracks
      std::array<ROOT::TThreadedObject<TH2F>, 4> distrSDPhiVsPTEMCalwNeg
      {
         ROOT::TThreadedObject<TH2F>("sdphi vs pT, EMCalw0, charge<0", 
                                     "sd#varphi_{EMCale0} vs p_{T}", 200, -5., 5., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("sdphi vs pT, EMCalw1, charge<0", 
                                     "sd#varphi_{EMCale1} vs p_{T}", 200, -5., 5., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("sdphi vs pT, EMCalw2, charge<0", 
                                     "sd#varphi_{EMCale2} vs p_{T}", 200, -5., 5., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("sdphi vs pT, EMCalw3, charge<0", 
                                     "sd#varphi_{EMCale3} vs p_{T}", 200, -5., 5., 100, 0., 10.)
      };
      /// emcsdz vs pT distributions for (0-3) sectors in west arm for negative tracks
      std::array<ROOT::TThreadedObject<TH2F>, 4> distrSDZVsPTEMCalwNeg
      {
         ROOT::TThreadedObject<TH2F>("sdz vs pT, EMCalw0, charge<0", "sdz_{EMCalw0} vs p_{T}", 
                                     200, -5., 5., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("sdz vs pT, EMCalw1, charge<0", "sdz_{EMCalw1} vs p_{T}", 
                                     200, -5., 5., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("sdz vs pT, EMCalw2, charge<0", "sdz_{EMCalw2} vs p_{T}", 
                                     200, -5., 5., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("sdz vs pT, EMCalw3, charge<0", "sdz_{EMCalw3} vs p_{T}", 
                                     200, -5., 5., 100, 0., 10.)
      };
      /// TOFe t-t_{exp}^{pi} distribution
      ROOT::TThreadedObject<TH2F> distrTTOFe{"t - t_exp^pi, TOFe", "t - t_{exp}^{#pi^{+}} vs p_{T}",
                                             1000., -20., 20., 100, 0., 10.};
      /// TOFw t-t_{exp}^{pi} distribution
      ROOT::TThreadedObject<TH2F> distrTTOFw{"t - t_exp^pi, TOFw", "t - t_{exp}^{#pi^{+}} vs p_{T}",
                                             1000., -20., 20., 100, 0., 10.};
      /// EMCale(0-3) t-t_{exp}^{pi} distributions
      std::array<ROOT::TThreadedObject<TH2F>, 4> distrTEMCale
      {
         ROOT::TThreadedObject<TH2F>("t - t_exp^pi, EMCale0", "t - t_{exp}^{#pi^{+}} vs p_{T}", 
                                     1000., -40., 40., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("t - t_exp^pi, EMCale1", "t - t_{exp}^{#pi^{+}} vs p_{T}", 
                                     1000., -40., 40., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("t - t_exp^pi, EMCale2", "t - t_{exp}^{#pi^{+}} vs p_{T}", 
                                     1000., -40., 40., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("t - t_exp^pi, EMCale3", "t - t_{exp}^{#pi^{+}} vs p_{T}", 
                                     1000., -40., 40., 100, 0., 10.)
      };
      /// EMCalw(0-3) t-t_{exp}^{pi} distributions
      std::array<ROOT::TThreadedObject<TH2F>, 4> distrTEMCalw
      {
         ROOT::TThreadedObject<TH2F>("t - t_exp^pi, EMCalw0", "t - t_{exp}^{#pi^{+}} vs p_{T}", 
                                     1000., -40., 40., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("t - t_exp^pi, EMCalw1", "t - t_{exp}^{#pi^{+}} vs p_{T}", 
                                     1000., -40., 40., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("t - t_exp^pi, EMCalw2", "t - t_{exp}^{#pi^{+}} vs p_{T}", 
                                     1000., -40., 40., 100, 0., 10.),
         ROOT::TThreadedObject<TH2F>("t - t_exp^pi, EMCalw3", "t - t_{exp}^{#pi^{+}} vs p_{T}", 
                                     1000., -40., 40., 100, 0., 10.)
      };
      /// TOFe m2 distribution for positive tracks
      ROOT::TThreadedObject<TH2F> distrM2TOFePosCharge{"m2, TOFe, charge>0", "p_{T} vs m^{2}",
                                                       100, 0., 10., 500, -0.2, 1.5};
      /// TOFe m2 distribution for negative tracks
      ROOT::TThreadedObject<TH2F> distrM2TOFeNegCharge{"m2, TOFe, charge<0", "p_{T} vs m^{2}",
                                                       100, 0., 10., 500, -0.2, 1.5};
      /// TOFw m2 distribution for positive tracks
      ROOT::TThreadedObject<TH2F> distrM2TOFwPosCharge{"m2, TOFw, charge>0", "p_{T} vs m^{2}",
                                                       100, 0., 10., 500, -0.2, 1.5};
      /// TOFw m2 distribution for negative tracks
      ROOT::TThreadedObject<TH2F> distrM2TOFwNegCharge{"m2, TOFw, charge<0", "p_{T} vs m^{2}",
                                                       100, 0., 10., 500, -0.2, 1.5};
      /// EMCale(0-3) m2 distributions for positive tracks
      std::array<ROOT::TThreadedObject<TH2F>, 4> distrM2EMCalePosCharge
      {
         ROOT::TThreadedObject<TH2F>("m2, EMCale0, charge>0", "p_{T} vs m^{2}", 
                                     100, 0., 10., 250, -0.5, 2.0),
         ROOT::TThreadedObject<TH2F>("m2, EMCale1, charge>0", "p_{T} vs m^{2}", 
                                     100, 0., 10., 250, -0.5, 2.0),
         ROOT::TThreadedObject<TH2F>("m2, EMCale2, charge>0", "p_{T} vs m^{2}", 
                                     100, 0., 10., 250, -0.5, 2.0),
         ROOT::TThreadedObject<TH2F>("m2, EMCale3, charge>0", "p_{T} vs m^{2}", 
                                     100, 0., 10., 250, -0.5, 2.0)
      };
      /// EMCale(0-3) m2 distributions for negative tracks
      std::array<ROOT::TThreadedObject<TH2F>, 4> distrM2EMCaleNegCharge
      {
         ROOT::TThreadedObject<TH2F>("m2, EMCale0, charge<0", "p_{T} vs m^{2}", 
                                     100, 0., 10., 250, -0.5, 2.0),
         ROOT::TThreadedObject<TH2F>("m2, EMCale1, charge<0", "p_{T} vs m^{2}", 
                                     100, 0., 10., 250, -0.5, 2.0),
         ROOT::TThreadedObject<TH2F>("m2, EMCale2, charge<0", "p_{T} vs m^{2}", 
                                     100, 0., 10., 250, -0.5, 2.0),
         ROOT::TThreadedObject<TH2F>("m2, EMCale3, charge<0", "p_{T} vs m^{2}", 
                                     100, 0., 10., 250, -0.5, 2.0)
      };
      /// EMCalw(0-3) m2 distributions for positive tracks
      std::array<ROOT::TThreadedObject<TH2F>, 4> distrM2EMCalwPosCharge
      {
         ROOT::TThreadedObject<TH2F>("m2, EMCalw0, charge>0", "p_{T} vs m^{2}", 
                                     100, 0., 10., 250, -0.5, 2.0),
         ROOT::TThreadedObject<TH2F>("m2, EMCalw1, charge>0", "p_{T} vs m^{2}", 
                                     100, 0., 10., 250, -0.5, 2.0),
         ROOT::TThreadedObject<TH2F>("m2, EMCalw2, charge>0", "p_{T} vs m^{2}", 
                                     100, 0., 10., 250, -0.5, 2.0),
         ROOT::TThreadedObject<TH2F>("m2, EMCalw3, charge>0", "p_{T} vs m^{2}", 
                                     100, 0., 10., 250, -0.5, 2.0)
      };
      /// EMCalw(0-3) m2 distributions for negative tracks
      std::array<ROOT::TThreadedObject<TH2F>, 4> distrM2EMCalwNegCharge
      {
         ROOT::TThreadedObject<TH2F>("m2, EMCalw0, charge<0", "p_{T} vs m^{2}", 
                                     100, 0., 10., 250, -0.5, 2.0),
         ROOT::TThreadedObject<TH2F>("m2, EMCalw1, charge<0", "p_{T} vs m^{2}", 
                                     100, 0., 10., 250, -0.5, 2.0),
         ROOT::TThreadedObject<TH2F>("m2, EMCalw2, charge<0", "p_{T} vs m^{2}", 
                                     100, 0., 10., 250, -0.5, 2.0),
         ROOT::TThreadedObject<TH2F>("m2, EMCalw3, charge<0", "p_{T} vs m^{2}", 
                                     100, 0., 10., 250, -0.5, 2.0)
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
   void AnalyzeConfiguration(ThrContainer& thrContainer, const std::string& particleName,
                             const short particleGeantId, const std::string& magneticFieldName, 
                             const std::string& pTRangeName); 
}

#endif /* ANALYZE_SINGLE_TRACK_HPP */
