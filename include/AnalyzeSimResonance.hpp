/** 
 *  @file   AnalyzeSimResonance.hpp
 *  @brief  Contains declarations of functions and variables that are used for analysis of a simulated resonance obtained from a trees acquired from the PHENIX simulation
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysisPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef ANALYZE_SIM_RESONANCE_HPP
#define ANALYZE_SIM_RESONANCE_HPP

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

#include "Constants.hpp"
#include "SingleTrackFunc.hpp"
#include "PairTrackFunc.hpp"
#include "SimTreeReader.hpp"
#include "DeadMapCutter.hpp"
#include "SimSigmalizedResiduals.hpp"
#include "SimM2Identificator.hpp"

#include "PBar.hpp"

#include "ROOT/TTreeProcessorMT.hxx"

/* @namespace AnalyzeSimResonance
 *
 * @brief Contains all functions, variables, and containers needed for AnalyzeSimResonance 
 *
 * This namespace is eployed so that documentation will not become a pile of variables, types, and functions from many different files that are intended to be compiled and used as executables. With this namespace finding the needed information for the given executable is easier since everything belongs to the current namespace
 */
namespace AnalyzeSimResonance
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
   /// file reader for all required parameters for the resonance and for its simulation processing
   InputYAMLReader inputYAMLResonance;
   /// file reader for all required parameters for the current run
   InputYAMLReader inputYAMLMain;
   /// file reader for all required parameters for the simulation processing 
   /// of widthless resonance decay products
   InputYAMLReader inputYAMLSimSingleTrack;
   /// pT scale (original pT as well as all single tracks) for systematics evaluation
   double pTScale = 1.;
   /// acceptance variation by the value of acceptance systematic uncertainty
   /// can have any value, but only 3 interpretations by the program:
   /// 1. 0 (by default) - no variation
   /// 2. negative - acceptance is decreased by its systematic uncertainty
   /// 2. positive - acceptance is increased by its systematic uncertainty
   int acceptanceVar = 0;
   /// cuts variation by the value of detectors resolution systematic uncertainty, by the width of sigmalized variables (sdphi, sdz), and by other fiducialy chosen variations
   /// can have any value, but only 3 interpretations by the program:
   /// 1. 0 (by default) - no variation
   /// 2. negative - cuts are tightened
   /// 2. positive - cuts are loosened
   int cutsVar = 0;
   /// sigmalized variation to the cuts
   const double cutsSigmOffset = 0.;
   /// number of threads
   int numberOfThreads;
   /// number of events across all trees
   unsigned long numberOfEvents = 0;
   /// parameter for monitoring the progress
   unsigned long numberOfCalls = 0;
   /// cutter for deadmaps
   DeadMapCutter dmCutter;
   /// calibrator for simulated sigmalized residuals of charged tracks sdphi and sdz
   SimSigmalizedResiduals simSigmRes;
   /// identificator for m2 hadron identification procedure in MC
   SimM2Identificator simM2Id;
   /// shows whether PC2 is used
   bool usePC2;
   /// shows whether PC3 is used
   bool usePC3;
   /// shows whether TOFe is used
   bool useTOFe;
   /// shows whether TOFw is used
   bool useTOFw;
   /// shows whether EMCal is used
   bool useEMCal;
   /// shows whether identification with EMCal is used
   bool useEMCalId;

   /* @struct ThrContainerCopy
    * @brief Container for storing local ThrContainer copies (at least 1 for each thread) 
    * since copying speeds up the calculation because local copies 
    * do not require synchronization from each other
    *
    * Contains the same histograms ThrContainer but as shared_ptr<T>
    */
   struct ThrContainerCopy
   {
      /// distribution of original unscaled generated pT
      std::shared_ptr<TH1D> distrOrigUnscaledPT;
      /// distribution of original generated pT
      std::shared_ptr<TH1D> distrOrigPT;
      /// distribution of original generated pT vs reconstructed pT of daughters
      std::shared_ptr<TH2F> distrOrigPTVsRecDaughtersPT;
      /// distribution of original generated pT vs reconstructed pT in the simulation
      std::shared_ptr<TH2F> distrOrigPTVsRecPT;
      /// distribution of original resonance pT vs pT of all particles in events (decay products and secondary particles)
      std::shared_ptr<TH2F> distrOrigPTVsDecayRecPT;
      /// invariant mass distribution with the cut opposite to one arm cut
      std::shared_ptr<TH2F> distrMInvOneArmAntiCut;
      /// TOFeNoPID invariant mass distribution with ghost anti cut
      std::shared_ptr<TH2F> distrMInvTOFeGhostNoPID;
      /// TOFwNoPID invariant mass distribution with ghost anti cut
      std::shared_ptr<TH2F> distrMInvTOFwGhostNoPID;
      /// EMCalNoPID invariant mass distribution with ghost anti cut
      std::shared_ptr<TH2F> distrMInvEMCalGhostNoPID;
      /// DCPC1NoPID invariant mass distribution
      std::shared_ptr<TH2F> distrMInvDCPC1NoPID;
      /// NoPID invariant mass distribution
      std::shared_ptr<TH2F> distrMInvNoPID;
      /// PC2NoPID invariant mass distribution
      std::shared_ptr<TH2F> distrMInvPC2NoPID;
      /// PC3NoPID invariant mass distribution
      std::shared_ptr<TH2F> distrMInvPC3NoPID;
      /// TOFeNoPID invariant mass distribution
      std::shared_ptr<TH2F> distrMInvTOFeNoPID;
      /// TOFwNoPID invariant mass distribution
      std::shared_ptr<TH2F> distrMInvTOFwNoPID;
      /// EMCalNoPID invariant mass distribution
      std::shared_ptr<TH2F> distrMInvEMCalNoPID;
      /// DCPC11PID invariant mass distribution
      std::shared_ptr<TH2F> distrMInvDCPC11PID;
      /// 1TOFDCPC11PID invariant mass distribution
      std::shared_ptr<TH2F> distrMInv1TOFDCPC11PID;
      /// 1EMCalDCPC11PID invariant mass distribution
      std::shared_ptr<TH2F> distrMInv1EMCalDCPC11PID;
      /// 1PID invariant mass distribution
      std::shared_ptr<TH2F> distrMInv1PID;
      /// 1TOF1PID invariant mass distribution
      std::shared_ptr<TH2F> distrMInv1TOF1PID;
      /// 1EMCal1PID invariant mass distribution
      std::shared_ptr<TH2F> distrMInv1EMCal1PID;
      /// 2PID invariant mass distribution
      std::shared_ptr<TH2F> distrMInv2PID;
      /// TOF2PID invariant mass distribution
      std::shared_ptr<TH2F> distrMInvTOF2PID;
      /// TOFe2PID invariant mass distribution
      std::shared_ptr<TH2F> distrMInvTOFe2PID;
      /// TOFw2PID invariant mass distribution
      std::shared_ptr<TH2F> distrMInvTOFw2PID;
      /// EMCal2PID invariant mass distribution
      std::shared_ptr<TH2F> distrMInvEMCal2PID;
      /// abs(p1 - p2)/(p1 + p2) vs pT vs Minv for pairs of tracks
      std::shared_ptr<TH3F> distrPAsymVsPT;
      /// phi1 - phi2 vs pT vs Minv for pairs of tracks
      std::shared_ptr<TH3F> distrDPhiVsPT;
      /// alpha1 - alpha2 vs pT vs Minv for pairs of tracks
      std::shared_ptr<TH3F> distrDAlphaVsPT;
      /// zed1 - zed2 vs pT vs Minv for pairs of tracks
      std::shared_ptr<TH3F> distrDZedVsPT;
      /// pc2phi1 - pc2phi2 vs pc2z1 - pc2z2 vs pT within 
      /// 2*Gamma + 10 MeV of the center of the signal for PC2
      std::shared_ptr<TH2F> distrDPC2PhiDPC2ZVsPT;
      /// pc2phi1 - pc2phi2 vs pc3z1 - pc3z2 vs pT in the same arm within
      /// 2*Gamma + 10 MeV of the center of the signal for PC3
      std::shared_ptr<TH2F> distrDPC3PhiDPC3ZVsPT;
      /// chamber1 - chamber2 vs slat1 - slat2 vs pT within 
      /// 2*Gamma + 10 MeV of the center of the signal for TOFe
      std::shared_ptr<TH3F> distrDChamberDSlatVsPT;
      /// chamber1 - chamber2 vs strip1 - strip2 vs pT within 
      /// 2*Gamma + 10 MeV of the center of the signal for TOFw
      std::shared_ptr<TH3F> distrDChamberDStripVsPT;
      /// ytower1 - ytower2 vs ztower1 - ztower2 vs pT within 
      /// 2*Gamma + 10 MeV of the center of the signal for the same sector of EMCal
      std::shared_ptr<TH3F> distrDYTowerDZTowerVsPT;
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
      /// distribution of original generated unscaled pT
      std::unique_ptr<ROOT::TThreadedObject<TH1D>> distrOrigUnscaledPT = 
         std::make_unique<ROOT::TThreadedObject<TH1D>>("orig unscaled pT", "p_{T}", 100., 0., 10.);
      /// distribution of original generated pT
      std::unique_ptr<ROOT::TThreadedObject<TH1D>> distrOrigPT = 
         std::make_unique<ROOT::TThreadedObject<TH1D>>("orig pT", "p_{T}", 100., 0., 10.);
      // distribution of original generated pT vs reconstructed pT of daughters in the simulation
      ROOT::TThreadedObject<TH2F> 
         distrOrigPTVsRecDaughtersPT{"orig pT vs rec daughters pT", "p_{T}^{orig} vs p_{T}^{rec}", 
                            100, 0., 10., 100, 0., 10.};
      // distribution of original generated pT vs reconstructed pT in the simulation
      ROOT::TThreadedObject<TH2F> 
         distrOrigPTVsRecPT{"orig pT vs rec pT", "p_{T}^{orig} vs p_{T}^{rec}", 
                            100, 0., 10., 100, 0., 10.};
      /// distribution of original resonance pT vs pT of all particles in events (decay products and secondary particles)
      ROOT::TThreadedObject<TH2F> 
         distrOrigPTVsDecayRecPT{"orig pT vs decay rec pT", "p_{T}^{orig} vs p_{T}^{rec}", 
                                 100, 0., 10., 100, 0., 10.};
      /// invariant mass distribution with the cut opposite to one arm cut
      ROOT::TThreadedObject<TH2F> distrMInvOneArmAntiCut{"M_inv: DCPC1NoPID, one arm anti cut", 
                                                         "M_{inv} vs p_{T}", 
                                                         200, 0., 20., 1000, 0., 5.};
      /// TOFeNoPID invariant mass distribution with TOFe ghost anticut
      ROOT::TThreadedObject<TH2F> distrMInvTOFeGhostNoPID{"M_inv: TOFeNoPID, ghost anticut", 
                                                          "M_{inv} vs p_{T}", 
                                                          200, 0., 20., 1000, 0., 5.};
      /// TOFwNoPID invariant mass distribution with TOFw ghost anticut
      ROOT::TThreadedObject<TH2F> distrMInvTOFwGhostNoPID{"M_inv: TOFwNoPID, ghost anticut", 
                                                          "M_{inv} vs p_{T}", 
                                                          200, 0., 20., 1000, 0., 5.};
      /// EMCalNoPID invariant mass distribution with EMCal ghost anticut
      ROOT::TThreadedObject<TH2F> distrMInvEMCalGhostNoPID{"M_inv: EMCalNoPID, ghost anticut", 
                                                           "M_{inv} vs p_{T}", 
                                                           200, 0., 20., 1000, 0., 5.};
      /// DCPC1NoPID invariant mass distribution
      ROOT::TThreadedObject<TH2F> distrMInvDCPC1NoPID{"M_inv: DCPC1NoPID", "M_{inv} vs p_{T}", 
                                                      200, 0., 20., 1000, 0., 5.};
      /// NoPID invariant mass distribution
      ROOT::TThreadedObject<TH2F> distrMInvNoPID{"M_inv: NoPID", "M_{inv} vs p_{T}", 
                                                 200, 0., 20., 1000, 0., 5.};
      /// PC2NoPID invariant mass distribution
      ROOT::TThreadedObject<TH2F> distrMInvPC2NoPID{"M_inv: PC2NoPID", "M_{inv} vs p_{T}", 
                                                    200, 0., 20., 1000, 0., 5.};
      /// PC3NoPID invariant mass distribution
      ROOT::TThreadedObject<TH2F> distrMInvPC3NoPID{"M_inv: PC3NoPID", "M_{inv} vs p_{T}", 
                                                    200, 0., 20., 1000, 0., 5.};
      /// TOFeNoPID invariant mass distribution
      ROOT::TThreadedObject<TH2F> distrMInvTOFeNoPID{"M_inv: TOFeNoPID", "M_{inv} vs p_{T}", 
                                                     200, 0., 20., 1000, 0., 5.};
      /// TOFwNoPID invariant mass distribution
      ROOT::TThreadedObject<TH2F> distrMInvTOFwNoPID{"M_inv: TOFwNoPID", "M_{inv} vs p_{T}", 
                                                     200, 0., 20., 1000, 0., 5.};
      /// EMCalNoPID invariant mass distribution
      ROOT::TThreadedObject<TH2F> distrMInvEMCalNoPID{"M_inv: EMCalNoPID", "M_{inv} vs p_{T}", 
                                                      200, 0., 20., 1000, 0., 5.};
      /// DCPC11PID invariant mass distribution
      ROOT::TThreadedObject<TH2F> distrMInvDCPC11PID{"M_inv: DCPC11PID", "M_{inv} vs p_{T}", 
                                                     200, 0., 20., 1000, 0., 5.};
      /// 1TOFDCPC11PID invariant mass distribution
      ROOT::TThreadedObject<TH2F> distrMInv1TOFDCPC11PID{"M_inv: 1TOFDCPC11PID", "M_{inv} vs p_{T}", 
                                                         200, 0., 20., 1000, 0., 5.};
      /// 1EMCalDCPC11PID invariant mass distribution
      ROOT::TThreadedObject<TH2F> distrMInv1EMCalDCPC11PID{"M_inv: 1EMCalDCPC11PID", 
                                                           "M_{inv} vs p_{T}", 
                                                           200, 0., 20., 1000, 0., 5.};
      /// 1PID invariant mass distribution
      ROOT::TThreadedObject<TH2F> distrMInv1PID{"M_inv: 1PID", "M_{inv} vs p_{T}", 
                                                200, 0., 20., 1000, 0., 5.};
      /// 1TOF1PID invariant mass distribution
      ROOT::TThreadedObject<TH2F> distrMInv1TOF1PID{"M_inv: 1TOF1PID", "M_{inv} vs p_{T}", 
                                                    200, 0., 20., 1000, 0., 5.};
      /// 1EMCal1PID invariant mass distribution
      ROOT::TThreadedObject<TH2F> distrMInv1EMCal1PID{"M_inv: 1EMCal1PID", "M_{inv} vs p_{T}", 
                                                      200, 0., 20., 1000, 0., 5.};
      /// 2PID invariant mass distribution
      ROOT::TThreadedObject<TH2F> distrMInv2PID{"M_inv: 2PID", "M_{inv} vs p_{T}", 
                                                200, 0., 20., 1000, 0., 5.};
      /// TOF2PID invariant mass distribution
      ROOT::TThreadedObject<TH2F> distrMInvTOF2PID{"M_inv: TOF2PID", "M_{inv} vs p_{T}", 
                                                   200, 0., 20., 1000, 0., 5.};
      /// TOFe2PID invariant mass distribution
      ROOT::TThreadedObject<TH2F> distrMInvTOFe2PID{"M_inv: TOFe2PID", "M_{inv} vs p_{T}", 
                                                    200, 0., 20., 1000, 0., 5.};
      /// TOFw2PID invariant mass distribution
      ROOT::TThreadedObject<TH2F> distrMInvTOFw2PID{"M_inv: TOFw2PID", "M_{inv} vs p_{T}", 
                                                    200, 0., 20., 1000, 0., 5.};
      /// EMCal2PID invariant mass distribution
      ROOT::TThreadedObject<TH2F> distrMInvEMCal2PID{"M_inv: EMCal2PID", "M_{inv} vs p_{T}", 
                                                     200, 0., 20., 1000, 0., 5.};
      /// abs(p1 - p2)/(p1 + p2) vs pT vs Minv for pairs of tracks
      ROOT::TThreadedObject<TH3F> 
         distrPAsymVsPT{"p asym", "(p_{pos} - p_{neg})/(p_{pos} + p_{neg}) vs p_{T} vs M_{inv}", 
                        100, 0., 20., 200, -1., 1., 100, 0., 5.};
      /// phi1 - phi2 vs pT vs Minv for pairs of tracks
      ROOT::TThreadedObject<TH3F> 
         distrDPhiVsPT{"delta phi", "#varphi_{pos} - #varphi_{neg} vs p_{T} vs M_{inv}", 
                       100, 0., 20., 200, -M_PI, M_PI, 100, 0., 5.};
      /// alpha1 - alpha2 vs pT vs Minv for pairs of tracks
      ROOT::TThreadedObject<TH3F> 
         distrDAlphaVsPT{"delta alpha", "#alpha_{pos} - #alpha_{neg} vs p_{T} vs M_{inv}", 
                         100, 0., 20., 200, -1., 0., 100, 0., 5.};
      /// zed1 - zed2 vs pT vs Minv for pairs of tracks
      ROOT::TThreadedObject<TH3F> 
         distrDZedVsPT{"delta zed", "zed_{pos} - zed_{neg} vs p_{T} vs M_{inv}", 
                       100, 0., 20., 200, -150., 150., 100, 0., 5.};
      /// pc2phi1 - pc2phi2 vs pc2z1 - pc2z vs pT within 
      /// 2*Gamma + 10 MeV of the center of the signal for the same sector of EMCal
      ROOT::TThreadedObject<TH2F> 
         distrDPC2PhiDPC2ZVsPT{"delta pc2phi vs delta pc2z", 
                                "z_{1}^{PC2} - z_{2}^{PC2} vs "\
                                "#varphi_{1}^{PC2} - #varphi_{2}^{PC2}",
                                200, -10., 10., 200, -0.1, 0.1};
      /// pc3phi1 - pc3phi2 vs pc3z1 - pc3z vs pT in the same arm within 
      /// 2*Gamma + 10 MeV of the center of the signal
      ROOT::TThreadedObject<TH2F> 
         distrDPC3PhiDPC3ZVsPT{"delta pc3phi vs delta pc3z", 
                                "z_{1}^{PC3} - z_{2}^{PC3} vs "\
                                "#varphi_{1}^{PC3} - #varphi_{2}^{PC3}",
                                200, -10., 10., 200., -0.1, 0.1};
      /// ytower1 - ytower2 vs ztower1 - ztower2 vs pT for the pair in the same sector within 
      /// 2*Gamma + 10 MeV of the center of the signal for the same sector of EMCal
      ROOT::TThreadedObject<TH3F> 
         distrDYTowerDZTowerVsPT{"delta ytower vs delta ztower", 
                                 "ytower_{1} - ytower_{2} vs ztower_{1} - ztower_{2} vs p_{T}",
                                 96, -48., 48., 194, -97., 97, 10, 0., 10.};
      /// chamber1 - chamber2 vs slat1 - slat2 vs pT within 
      /// 2*Gamma + 10 MeV of the center of the signal for TOFe
      ROOT::TThreadedObject<TH3F> 
         distrDChamberDSlatVsPT{"delta chamber vs delta slat", 
                                "chamber_{1} - chamber_{2} vs slat_{1} - slat_{2} vs p_{T}",
                                20, -10., 10., 192, -96., 96, 10, 0., 10.};
      /// chamber1 - chamber2 vs strip1 - strip2 vs pT within 
      /// 2*Gamma + 10 MeV of the center of the signal for TOFw
      ROOT::TThreadedObject<TH3F> 
         distrDChamberDStripVsPT{"delta chamber vs delta strip", 
                                "chamber_{1} - chamber_{2} vs strip_{1} - strip_{2} vs p_{T}",
                                20, -10., 10., 128, -64., 64, 10, 0., 10.};
   };
   /* @brief Processes the single configuration (for the given particle, 
    * magnetic field, and pT range) from one file
    *
    * @param[in] thrContainer the current ThrContainer in which the data will be written to 
    * @param[in] particleName name of the particle to be analyzed 
    * @param[in] pTRangeName name of the pT range to be analyzed 
    */
   void AnalyzeConfiguration(ThrContainer& thrContainer, const std::string& particleName,
                             const int daughter1Id, const int daugther2Id,
                             const std::string& magneticFieldName, const std::string& pTRangeName);
   /* @struct AcceptanceVar
    * @brief Container for storing acceptance variations for different detectors
    */
   struct AcceptanceVar
   {
      /// @brief sets acceptance variation values from the specified file
      void Set(const std::string& fileName);
      /// value of acceptance variation for DCe0
      double DCe0 = 0.;
      /// value of acceptance variation for DCe1
      double DCe1 = 0.;
      /// value of acceptance variation for DCw0
      double DCw0 = 0.;
      /// value of acceptance variation for DCw1
      double DCw1 = 0.;
      /// values of acceptance variation for PC1e and PC1w
      std::array<double, 2> PC1 = {0., 0.};
      /// value of acceptance variation for PC2
      double PC2 = 0.;
      /// values of acceptance variation for PC3e and PC3w
      std::array<double, 2> PC3 = {0., 0.};
      /// value of acceptance variation for TOFe
      double TOFe = 0.;
      /// value of acceptance variation for TOFw
      double TOFw = 0.;
      /// values of acceptance variation for EMCale(0-3) and EMCalw(0-3)
      std::array<std::array<double, 4>, 2> EMCal{{{{0., 0., 0., 0.}}, {{0., 0., 0., 0.}}}};
   };
   /* @struct CutsVar
    * @brief Container for storing cuts variations for different detectors
    */
   struct CutsVar
   {
      /// @brief sets acceptance variation values from the specified file
      void Set();
      /// value of bbcz cut variation
      double bbcz = 0.;
      /// value of zed cut variation
      double zed = 0.;
      /// value of sdphi variation for matching cut (for all detectors)
      double sdphi = 0.;
      /// value of sdz variation for matching cut (for all detectors)
      double sdz = 0.;
      /// value of sigmalized m^2 variation for m^2_{id} identification region cut (for all detectors)
      double m2Id = 0.;
      /// value of sigmalized m^2 variation for m^2_{id} identification veto cut (for all detectors)
      double m2Veto = 0.;
      /// value of ecore variation for E_{core} energy cut in EMCal
      double ecore = 0.;
      /// value of etof variation for E_{loss} energy cut in TOFe
      double etof = 0.;
   };
}

#endif /* ANALYZE_SIM_RESONANCE_HPP */
