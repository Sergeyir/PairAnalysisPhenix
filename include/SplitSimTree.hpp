/** 
 *  @file   SplitSimTree.hpp
 *  @brief  Contains declarations of functions and variables that are used for splitting one TFile with simulated tree into TFile containing TTree with smaller pT range. This is useful when more statistics is needed on the given pT range while the simulated Tree pT range is wider.
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysisPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef SPLIT_SIM_TREE_HPP
#define SPLIT_SIM_TREE_HPP

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

int main(int argc, char **argv);

/* @namespace SplitSimTree
 *
 * @brief Contains all functions, variables, and containers needed for SplitSimTree 
 *
 * This namespace is eployed so that documentation will not become a pile of variables, types, and functions from many different files that are intended to be compiled and used as executables. With this namespace finding the needed information for the given executable is easier since everything belongs to the current namespace
 */
namespace SplitSimTree
{
   /// output directory
   std::string outputDir;
   /// collision system name
   std::string collisionSystemName;
   /// run name
   std::string runName;
   /// minimum pT cut for smaller TTree
   double pTMin;
   /// maximum pT cut for smaller TTree
   double pTMax;
   /// number of threads
   int numberOfThreads;
   /// number of events across all trees
   unsigned long numberOfEvents = 0;
   /// parameter for monitoring the progress
   unsigned long numberOfCalls = 0;

   /* @struct ThrContainerCopy
    * @brief Container for storing local ThrContainer copies (at least 1 for each thread) 
    * since copying speeds up the calculation because local copies 
    * do not require synchronization from each other
    *
    * Names of <TH1D> variables are identical to those of PHCentralTrack and those, 
    * that can be obtained with SimTreeReader
    *
    * Contains the same histograms ThrContainer but as shared_ptr<T>
    */
   struct ThrContainerCopy
   {
      /// smaller pT region tree
      std::shared_ptr<TTree> smallTree;
      std::shared_ptr<TH1D> nchDistr;
      std::shared_ptr<TH1D> bbczDistr;
      std::shared_ptr<TH1D> mom_origDistr;
      std::shared_ptr<TH1D> dcarmDistr;
      std::shared_ptr<TH1D> phiDistr;
      std::shared_ptr<TH1D> alphaDistr;
      std::shared_ptr<TH1D> zedDistr;
      std::shared_ptr<TH1D> momDistr;
      std::shared_ptr<TH1D> the0Distr;
      std::shared_ptr<TH1D> phi0Distr;
      std::shared_ptr<TH1D> nx1hitsDistr;
      std::shared_ptr<TH1D> nx2hitsDistr;
      std::shared_ptr<TH1D> qualDistr;
      std::shared_ptr<TH1D> chargeDistr;
      std::shared_ptr<TH1D> parent_idDistr;
      std::shared_ptr<TH1D> primary_idDistr;
      std::shared_ptr<TH1D> particle_idDistr;
      std::shared_ptr<TH1D> ttofDistr;
      std::shared_ptr<TH1D> ttofwDistr;
      std::shared_ptr<TH1D> temcDistr;
      std::shared_ptr<TH1D> pltofDistr;
      std::shared_ptr<TH1D> pltofwDistr;
      std::shared_ptr<TH1D> plemcDistr;
      std::shared_ptr<TH1D> ptofxDistr;
      std::shared_ptr<TH1D> ptofyDistr;
      std::shared_ptr<TH1D> ptofzDistr;
      std::shared_ptr<TH1D> ptofwxDistr;
      std::shared_ptr<TH1D> ptofwyDistr;
      std::shared_ptr<TH1D> ptofwzDistr;
      std::shared_ptr<TH1D> pemcxDistr;
      std::shared_ptr<TH1D> pemcyDistr;
      std::shared_ptr<TH1D> pemczDistr;
      std::shared_ptr<TH1D> ppc1xDistr;
      std::shared_ptr<TH1D> ppc1yDistr;
      std::shared_ptr<TH1D> ppc1zDistr;
      std::shared_ptr<TH1D> ppc2xDistr;
      std::shared_ptr<TH1D> ppc2yDistr;
      std::shared_ptr<TH1D> ppc2zDistr;
      std::shared_ptr<TH1D> ppc3xDistr;
      std::shared_ptr<TH1D> ppc3yDistr;
      std::shared_ptr<TH1D> ppc3zDistr;
      std::shared_ptr<TH1D> ptecxDistr;
      std::shared_ptr<TH1D> ptecyDistr;
      std::shared_ptr<TH1D> pteczDistr;
      std::shared_ptr<TH1D> tofdzDistr;
      std::shared_ptr<TH1D> tofdphiDistr;
      std::shared_ptr<TH1D> tofwdzDistr;
      std::shared_ptr<TH1D> tofwdphiDistr;
      std::shared_ptr<TH1D> emcdzDistr;
      std::shared_ptr<TH1D> emcdphiDistr;
      std::shared_ptr<TH1D> pc2dzDistr;
      std::shared_ptr<TH1D> pc2dphiDistr;
      std::shared_ptr<TH1D> pc3dzDistr;
      std::shared_ptr<TH1D> pc3dphiDistr;
      std::shared_ptr<TH1D> striptofwDistr;
      std::shared_ptr<TH1D> slatDistr;
      std::shared_ptr<TH1D> etofDistr;
      std::shared_ptr<TH1D> ecoreDistr;
      std::shared_ptr<TH1D> emceDistr;
      std::shared_ptr<TH1D> ecentDistr;
      std::shared_ptr<TH1D> e9Distr;
      std::shared_ptr<TH1D> emcchi2Distr;
      std::shared_ptr<TH1D> twrhitDistr;
      std::shared_ptr<TH1D> emcdispyDistr;
      std::shared_ptr<TH1D> emcdispzDistr;
      std::shared_ptr<TH1D> probDistr;
      std::shared_ptr<TH1D> sectDistr;
      std::shared_ptr<TH1D> ysectDistr;
      std::shared_ptr<TH1D> zsectDistr;
      std::shared_ptr<TH1D> n0Distr;
      std::shared_ptr<TH1D> npe0Distr;
      std::shared_ptr<TH1D> n1Distr;
      std::shared_ptr<TH1D> npe1Distr;
      std::shared_ptr<TH1D> n2Distr;
      std::shared_ptr<TH1D> npe2Distr;
      std::shared_ptr<TH1D> n3Distr;
      std::shared_ptr<TH1D> npe3Distr;
      std::shared_ptr<TH1D> center_phiDistr;
      std::shared_ptr<TH1D> center_zDistr;
      std::shared_ptr<TH1D> cross_phiDistr;
      std::shared_ptr<TH1D> cross_zDistr;
      std::shared_ptr<TH1D> dispDistr;
      std::shared_ptr<TH1D> chi2Distr;
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
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> nchDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("nch", "nch", 20, 0., 20.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> bbczDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("bbcz", "bbcz", 100, -50., 50.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> mom_origDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("mom_orig", "mom_orig", 200, -20., 20.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> dcarmDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("dcarm", "dcarm", 2, 0., 2.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> phiDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("phi", "phi", 64, -1., 4.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> alphaDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("alpha", "alpha", 100, -1.5, 1.5);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> zedDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("zed", "zed", 100, -100., 100.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> momDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("mom", "mom", 200, 0., 20.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> the0Distr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("the0", "the0", 100, 1., 2.2);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> phi0Distr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("phi0", "phi0", 100, -2., 5.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> nx1hitsDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("nx1hits", "nx1hits", 20, 0., 20.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> nx2hitsDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("nx2hits", "nx2hits", 20, 0., 20.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> qualDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("qual", "qual", 100, 0., 100.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> chargeDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("charge", "charge", 6, -3., 3.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> parent_idDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("parent_id", "parent_id", 20, 0., 20.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> primary_idDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("primary_id", "primary_id", 20, 0., 20.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> particle_idDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("particle_id", "particle_id", 20, 0., 20.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> ttofDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("ttof", "ttof", 200, 0., 50.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> ttofwDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("ttofw", "ttofw", 200, 0., 50.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> temcDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("temc", "temc", 200, 0., 50.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> pltofDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("pltof", "pltof", 100, 480., 650.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> pltofwDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("pltofw", "pltofw", 100, 460., 600.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> plemcDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("plemc", "plemc", 100, 500., 660.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> ptofxDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("ptofx", "ptofx", 100, -520., -420.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> ptofyDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("ptofy", "ptofy", 100, -300., 100.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> ptofzDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("ptofz", "ptofz", 100, -210., 210.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> ptofwxDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("ptofwx", "ptofwx", 100, 390., 500.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> ptofwyDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("ptofwy", "ptofwy", 100, -100., 300.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> ptofwzDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("ptofwz", "ptofwz", 100, -200., 200.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> pemcxDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("pemcx", "pemcx", 100, -800., 800.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> pemcyDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("pemcy", "pemcy", 100, -320., 460.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> pemczDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("pemcz", "pemcz", 100, -220., 220.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> ppc1xDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("ppc1x", "ppc1x", 100, -300., 300.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> ppc1yDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("ppc1y", "ppc1y", 100, -200., 250.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> ppc1zDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("ppc1z", "ppc1z", 100, -100., 100.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> ppc2xDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("ppc2x", "ppc2x", 100, -450., 450.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> ppc2yDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("ppc2y", "ppc2y", 100, -300., 400.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> ppc2zDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("ppc2z", "ppc2z", 100, -200., 200.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> ppc3xDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("ppc3x", "ppc3x", 100, -550., 550.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> ppc3yDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("ppc3y", "ppc3y", 100, -300., 300.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> ppc3zDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("ppc3z", "ppc3z", 100, -200., 200.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> ptecxDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("ptecx", "ptecx", 100, -550., 550.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> ptecyDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("ptecy", "ptecy", 100, -300., 300.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> pteczDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("ptecz", "ptecz", 100, -200., 200.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> tofdzDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("tofdz", "tofdz", 100, -50., 50.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> tofdphiDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("tofdphi", "tofdphi", 100, -0.2, 0.2);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> tofwdzDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("tofwdz", "tofwdz", 100, -50., 50.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> tofwdphiDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("tofwdphi", "tofwdphi", 100, -0.2, 0.2);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> emcdzDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("emcdz", "emcdz", 100, -60., 60.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> emcdphiDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("emcdphi", "emcdphi", 100, -0.2, 0.2);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> pc2dzDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("pc2dz", "pc2dz", 100 , -50., 50.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> pc2dphiDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("pc2dphi", "pc2dphi", 100, -0.2, 0.2);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> pc3dzDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("pc3dz", "pc3dz", 100, -50., 50.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> pc3dphiDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("pc3dphi", "pc3dphi", 100, - 0.2, 0.2);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> striptofwDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("striptofw", "striptofw", 200, 0., 550.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> slatDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("slat", "slat", 200, 0., 1000.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> etofDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("etof", "etof", 100, 0., 0.03);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> ecoreDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("ecore", "ecore", 100, 0., 5.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> emceDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("emce", "emce", 100, 0., 5.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> ecentDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("ecent", "ecent", 100, 0., 5.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> e9Distr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("e9", "e9", 100, 0., 5.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> emcchi2Distr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("emcchi2", "emcchi2", 100, 0., 200.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> twrhitDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("twrhit", "twrhit", 20, 0., 20.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> emcdispyDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("emcdispy", "emcdispy", 100, 0., 20.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> emcdispzDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("emcdispz", "emcdispz", 100, 0., 20.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> probDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("prob", "prob", 100, 0., 1.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> sectDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("sect", "sect", 4, 0., 4.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> ysectDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("ysect", "ysect", 48, 0., 48.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> zsectDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("zsect", "zsect", 97, 0., 97.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> n0Distr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("n0", "n0", 20, 0., 20.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> npe0Distr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("npe0", "npe0", 20, 0., 20.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> n1Distr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("n1", "n1", 20, 0., 20.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> npe1Distr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("npe1", "npe1", 20, 0., 20.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> n2Distr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("n2", "n2", 20, 0., 20.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> npe2Distr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("npe2", "npe2", 20, 0., 20.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> n3Distr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("n3", "n3", 20, 0., 20.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> npe3Distr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("npe3", "npe3", 20, 0., 20.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> center_phiDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("center_phi", "center_phi", 
                                                        100, -2.*M_PI, 2.*M_PI);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> center_zDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("center_z", "center_z", 100, -500., 500.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> cross_phiDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("cross_phi", "cross_phi", 
                                                        100, -2.*M_PI, 2.*M_PI);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> cross_zDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("cross_z", "cross_z", 100, -500., 500.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> dispDistr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("disp", "disp", 100, 0., 10.);
       std::unique_ptr<ROOT::TThreadedObject<TH1D>> chi2Distr = 
          std::make_unique<ROOT::TThreadedObject<TH1D>>("chi2", "chi2", 100, 0., 200.);
}

#endif /* SPLIT_SIM_TREE_HPP */
