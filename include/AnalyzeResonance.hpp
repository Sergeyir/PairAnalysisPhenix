/** 
 *  @file   AnalyzeResonance.hpp
 *  @brief  Contains declarations of functions and variables that are used for analysis of a simulated resonance obtained from a trees acquired from the PHENIX simulation
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysisPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef ANALYZE_RESONANCE_HPP
#define ANALYZE_RESONANCE_HPP

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
#include "SimCalibrator.hpp"

#include "PBar.hpp"

#include "ROOT/TTreeProcessorMT.hxx"

int main(int argc, char **argv);

/* @namespace AnalyzeResonance
 *
 * @brief Contains all functions, variables, and containers needed for AnalyzeResonance 
 *
 * This namespace is eployed so that documentation will not become a pile of variables, types, and functions from many different files that are intended to be compiled and used as executables. With this namespace finding the needed information for the given executable is easier since everything belongs to the current namespace
 */
namespace AnalyzeResonance
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
   /// file reader for all required parameters for the simulation processing of widthless resonance
   InputYAMLReader inputYAMLSim;
   /// file reader for all required parameters for the current run
   InputYAMLReader inputYAMLMain;
   /// file reader for all required parameters for the simulation processing 
   /// of widthless resonance decay products
   InputYAMLReader inputYAMLSimSingleTrack;
   /// number of threads
   int numberOfThreads;
   /// number of events across all trees
   unsigned long numberOfEvents = 0;
   /// parameter for monitoring the progress
   unsigned long numberOfCalls = 0;
   /// cutter for deadmaps
   DeadMapCutter dmCutter;
   /// calibrator for simulated data
   SimCalibrator simCalibrator;

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
      std::shared_ptr<TH1D> distrOrigPT;
      /// distribution of original generated pT vs reconstructed pT in the simulation
      std::shared_ptr<TH2F> distrOrigPTVsRecPT;
      /// invariant mass distribution for all pairs without any cuts on pairs
      std::shared_ptr<TH2F> distrMInv;
      /// NoPID invariant mass distribution
      std::shared_ptr<TH2F> distrMInvNoPID;
      /// NoPID invariant mass distribution with the cut opposite to one arm cut
      std::shared_ptr<TH2F> distrMInvNoPIDOneArmAntiCut;
      /// NoPID invariant mass distribution with the cut opposite to ghost cut
      std::shared_ptr<TH2F> distrMInvNoPIDGhostAntiCut;
      /// NoPID invariant mass distribution with sailor cut applied
      std::shared_ptr<TH2F> distrMInvNoPIDSailorCut;
      /// NoPID invariant mass distribution with cowboy cut applied
      std::shared_ptr<TH2F> distrMInvNoPIDCowboyCut;
      /// abs(E1 - E2)/(E1 + E2) for a pair of tracks within 2*Gamma + 10 MeV 
      /// of the center of the signal
      std::shared_ptr<TH3F> distrEAsymVsPT;
      /// abs(p1 - p2)/(p1 + p2) for a pair of tracks within 2*Gamma + 10 MeV 
      /// of the center of the signal
      std::shared_ptr<TH3F> distrPAsymVsPT;
      /// E1 - E2 for a pair of tracks within 2*Gamma + 10 MeV of the center of the signal
      std::shared_ptr<TH3F> distrDEVsPT;
      /// p1 - p2 for a pair of tracks within 2*Gamma + 10 MeV of the center of the signal
      std::shared_ptr<TH3F> distrDPVsPT;
      /// phi1 - phi2 for a pair of tracks within 2*Gamma + 10 MeV of the center of the signal
      std::shared_ptr<TH3F> distrDPhiVsPT;
      /// alpha1 - alpha2 for a pair of tracks within 2*Gamma + 10 MeV of the center of the signal
      std::shared_ptr<TH3F> distrDAlphaVsPT;
      /// zed1 - zed2 for a pair of tracks within 2*Gamma + 10 MeV of the center of the signal
      std::shared_ptr<TH3F> distrDZedVsPT;
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
      std::unique_ptr<ROOT::TThreadedObject<TH1D>> distrOrigPT = 
         std::make_unique<ROOT::TThreadedObject<TH1D>>("orig pT", "p_{T}", 100., 0., 10.);
      // distribution of original generated pT vs reconstructed pT in the simulation
      ROOT::TThreadedObject<TH2F> 
         distrOrigPTVsRecPT{"orig pT vs rec pT", "p_{T}^{orig} vs p_{T}^{rec}", 
                            100, 0., 10., 100, 0., 10.};
      /// invariant mass distribution for all pair without any cuts on pairs
      ROOT::TThreadedObject<TH2F> distrMInv{"M_inv: all pairs", "M_{inv} vs p_{T}", 
                                            200, 0., 20., 1000, 0., 5.};
      /// NoPID invariant mass distribution
      ROOT::TThreadedObject<TH2F> distrMInvNoPID{"M_inv: NoPID", "M_{inv} vs p_{T}", 
                                                 200, 0., 20., 1000, 0., 5.};
      /// NoPID invariant mass distribution with the cut opposite to one arm cut
      ROOT::TThreadedObject<TH2F> distrMInvNoPIDOneArmAntiCut{"M_inv: NoPID, one arm anti cut", 
                                                              "M_{inv} vs p_{T}", 
                                                              200, 0., 20., 1000, 0., 5.};
      /// NoPID invariant mass distribution with the cut opposite to ghost cut
      ROOT::TThreadedObject<TH2F> distrMInvNoPIDGhostAntiCut{"M_inv: NoPID, ghost anti cut", 
                                                             "M_{inv} vs p_{T}", 
                                                             200, 0., 20., 1000, 0., 5.};
      /// NoPID invariant mass distribution with sailor cut applied
      ROOT::TThreadedObject<TH2F> distrMInvNoPIDSailorCut{"M_inv: NoPID, sailor cut", 
                                                          "M_{inv} vs p_{T}", 
                                                          200, 0., 20., 1000, 0., 5.};
      /// NoPID invariant mass distribution with cowboy cut applied
      ROOT::TThreadedObject<TH2F> distrMInvNoPIDCowboyCut{"M_inv: NoPID, cowboy cut", 
                                                          "M_{inv} vs p_{T}", 
                                                          200, 0., 20., 1000, 0., 5.};
      /// abs(E1 - E2)/(E1 + E2) for a pair of tracks within 2*Gamma + 10 MeV 
      /// of the center of the signal
      ROOT::TThreadedObject<TH3F> 
         distrEAsymVsPT{"E asym", "(E_{pos} - E_{neg})/(E_{pos} + E_{neg}) vs p_{T}", 
                        100, 0., 20., 200, -1., 1., 100, 0., 5.};
      /// abs(p1 - p2)/(p1 + p2) for a pair of tracks within 2*Gamma + 10 MeV 
      /// of the center of the signal
      ROOT::TThreadedObject<TH3F> 
         distrPAsymVsPT{"p asym", "(p_{pos} - p_{neg})/(p_{pos} + p_{neg}) vs p_{T}", 
                        100, 0., 20., 200, -1., 1., 100, 0., 5.};
      /// E1 - E2 for a pair of tracks within 2*Gamma + 10 MeV of the center of the signal
      ROOT::TThreadedObject<TH3F> distrDEVsPT{"delta E", "E_{pos} - E_{neg} vs p_{T}", 
                                              100, 0., 20., 200, -20., 20., 100, 0., 5.};
      /// p1 - p2 for a pair of tracks within 2*Gamma + 10 MeV of the center of the signal
      ROOT::TThreadedObject<TH3F> distrDPVsPT{"delta p", "p_{pos} - p_{neg} vs p_{T}", 
                                              100, 0., 20., 200, -20., 20., 100, 0., 5.};
      /// phi1 - phi2 for a pair of tracks within 2*Gamma + 10 MeV of the center of the signal
      ROOT::TThreadedObject<TH3F> 
         distrDPhiVsPT{"delta phi", "#varphi_{pos} - #varphi_{neg} vs p_{T}", 
                       100, 0., 20., 200, -M_PI, M_PI, 100, 0., 5.};
      /// alpha1 - alpha2 for a pair of tracks within 2*Gamma + 10 MeV of the center of the signal
      ROOT::TThreadedObject<TH3F> 
         distrDAlphaVsPT{"delta alpha", "#alpha_{pos} - #alpha_{neg} vs p_{T}", 
                         100, 0., 20., 200, -1., 0., 100, 0., 5.};
      /// zed1 - zed2 for a pair of tracks within 2*Gamma + 10 MeV of the center of the signal
      ROOT::TThreadedObject<TH3F> distrDZedVsPT{"delta zed", "zed_{pos} - zed_{neg} vs p_{T}", 
                                                100, 0., 20., 200, -150., 150., 100, 0., 5.};
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
}

#endif /* ANALYZE_RESONANCE_HPP */
