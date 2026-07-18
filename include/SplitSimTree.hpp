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
#include "TTreeReader.h"

#include "IOTools.hpp"
#include "StrTools.hpp"

#include "SimTreeReader.hpp"

#include "PBar.hpp"

/* @namespace SplitSimTree
 *
 * @brief Contains all functions, variables, and containers needed for SplitSimTree 
 *
 * This namespace is eployed so that documentation will not become a pile of variables, types, and functions from many different files that are intended to be compiled and used as executables. With this namespace finding the needed information for the given executable is easier since everything belongs to the current namespace
 */
namespace SplitSimTree
{
   /// minimum pT cut for smaller TTree
   double pTMin;
   /// maximum pT cut for smaller TTree
   double pTMax;
   /// number of events across all trees
   unsigned long numberOfEvents = 0;
   /// parameter for monitoring the progress
   unsigned long numberOfCalls = 0;

   // the following code was copied from /gpfs/mnt/gpfs02/phenix/plhf/plhf1/dlario/sim/tree_maker/DataAnalyzer.C to avoid any possible mismatch

   TTree tree("Tree", "Tree");

   TH1D distrOrigPT = TH1D("orig_pt", "pt", 200, 0., 20);
   TH2D distrOrigPTVsPT = TH2D("distrOrigPTVsPT", "pt", 200, 0., 20, 200, 0., 20);;
   
   // namespace that stores all histograms added with structures below;
   // this namespace allows straighforward access to the histograms and 
   // provides an easy way to write them
   namespace HistHolder
   {
      std::vector<TH1D *> vec;
      void Write();
      void Clear();
      unsigned int CurrentPos();
   };

   // below structs hold variable(s) and help fill the histogarms for this/these variable(s)
   /// structure for storing int variable
   struct SimVarInt
   {
      SimVarInt(const std::string& name, const double xMin, 
                const double xMax, const int nbins = 100);
      void SetValue(const int newVal);
      int val;
      unsigned int histHolderIndex;
   }; 
   /// structure for storing short variable
   struct SimVarShort
   {	
      SimVarShort(const std::string& name, const double xMin, 
                  const double xMax, const int nbins = 100);
      void SetValue(const short newVal);
      short val;
      unsigned int histHolderIndex;
   }; 
   /// structure for storing float variable
   struct SimVarF
   {	
      SimVarF(const std::string& name, const double xMin, 
              const double xMax, const int nbins = 100);
      void SetValue(const float newVal);
      float val;
      unsigned int histHolderIndex;
   }; 
   /// structure for storing int array
   struct SimVarIntArr
   {	
      SimVarIntArr(const std::string& name, const double xMin, const double xMax, 
                   const std::string& size, const int nbins = 100);
      void SetValue(const int newVal, const unsigned int index);
      int val[50];
      unsigned int histHolderIndex;
   }; 
   /// structure for storing short array
   struct SimVarShortArr
   {	
      SimVarShortArr(const std::string& name, const double xMin, const double xMax, 
                     const std::string& size, const int nbins = 100);
      void SetValue(const short newVal, const unsigned int index);
      short val[50];
      unsigned int histHolderIndex;
   }; 
   /// structure for storing float array
   struct SimVarFArr
   {	
      SimVarFArr(const std::string& name, const double xMin, const double xMax, 
                 const std::string& size, const int nbins = 100);
      void SetValue(const float newVal, const unsigned int index);
      float val[50];
      unsigned int histHolderIndex;
   }; 

   namespace Container
   {
      // name of variables and branches are written similarly to the names in PHCentralTrack
      //event variables
      SimVarInt nch("nch", 0., 20., 20);
      SimVarF bbcz("bbcz", -50, 50, 100);
      SimVarFArr mom_orig("mom_orig", -20, 20, "3", 200);
      
      SimVarShortArr dcarm("dcarm", 0., 2., "nch", 2);
      SimVarFArr phi("phi", -1., 4., "nch", 64);
      SimVarFArr alpha("alpha", -1.5, 1.5, "nch", 100);
      SimVarFArr zed("zed", -100., 100., "nch", 100);
      SimVarFArr mom("mom", 0., 20., "nch", 200);
      SimVarFArr the0("the0", 1., 2.2, "nch", 100);
      SimVarFArr phi0("phi0", -2., 5., "nch", 100);

      SimVarShortArr nx1hits("nx1hits", 0., 20., "nch", 20);
      SimVarShortArr nx2hits("nx2hits", 0., 20., "nch", 20);
      
      SimVarShortArr qual("qual", 0., 100., "nch", 100);
      SimVarShortArr charge("charge", -3., 3., "nch", 6);
      
      SimVarShortArr parent_id("parent_id", 0., 20, "nch", 20);
      SimVarShortArr primary_id("primary_id", 0., 20, "nch", 20);
      SimVarShortArr particle_id("particle_id", 0., 20, "nch", 20);

      SimVarFArr ttof("ttof", 0., 50., "nch", 200);
      SimVarFArr ttofw("ttofw", 0., 50., "nch", 200);
      SimVarFArr temc("temc", 0., 50., "nch", 200);

      SimVarFArr pltof("pltof", 480., 650., "nch", 100);
      SimVarFArr pltofw("pltofw", 460., 600., "nch", 100);
      SimVarFArr plemc("plemc", 500., 660., "nch", 100);

      SimVarFArr ptofx("ptofx", -520., -420., "nch", 100);
      SimVarFArr ptofy("ptofy", -300., 100., "nch", 100);
      SimVarFArr ptofz("ptofz", -210., 210., "nch", 100);

      SimVarFArr ptofwx("ptofwx", 390., 500., "nch", 100);
      SimVarFArr ptofwy("ptofwy", -100., 300., "nch", 100);
      SimVarFArr ptofwz("ptofwz", -200., 200., "nch", 100);

      SimVarFArr pemcx("pemcx", -800., 800., "nch", 100);
      SimVarFArr pemcy("pemcy", -320., 460., "nch", 100);
      SimVarFArr pemcz("pemcz", -220., 220., "nch", 100);

      SimVarFArr ppc1x("ppc1x", -300., 300., "nch", 100);
      SimVarFArr ppc1y("ppc1y", -200., 250., "nch", 100);
      SimVarFArr ppc1z("ppc1z", -100., 100., "nch", 100);

      SimVarFArr ppc2x("ppc2x", -450., 450., "nch", 100);
      SimVarFArr ppc2y("ppc2y", -300., 400., "nch", 100);
      SimVarFArr ppc2z("ppc2z", -200., 200., "nch", 100);
      
      SimVarFArr ppc3x("ppc3x", -550., 550., "nch", 100);
      SimVarFArr ppc3y("ppc3y", -300., 300., "nch", 100);
      SimVarFArr ppc3z("ppc3z", -200., 200., "nch", 100);

      SimVarFArr ptecx("ptecx", -550., 550., "nch", 100);
      SimVarFArr ptecy("ptecy", -300., 300., "nch", 100);
      SimVarFArr ptecz("ptecz", -200., 200., "nch", 100);

      SimVarFArr tofdz("tofdz", -50., 50., "nch", 100);
      SimVarFArr tofdphi("tofdphi", -0.2, 0.2, "nch", 100);

      SimVarFArr tofwdz("tofwdz", -50., 50., "nch", 100);
      SimVarFArr tofwdphi("tofwdphi", -0.2, 0.2, "nch", 100);
      
      SimVarFArr emcdz("emcdz", -60., 60., "nch", 100);
      SimVarFArr emcdphi("emcdphi", -0.2, 0.2, "nch", 100);
      
      SimVarFArr pc2dz("pc2dz", -50., 50., "nch", 100);
      SimVarFArr pc2dphi("pc2dphi", -0.2, 0.2, "nch", 100);

      SimVarFArr pc3dz("pc3dz", -50., 50., "nch", 100);
      SimVarFArr pc3dphi("pc3dphi", -0.2, 0.2, "nch", 100);

      SimVarShortArr striptofw("striptofw", 0., 550., "nch", 200);
      SimVarShortArr slat("slat", 0., 1000., "nch", 200);
      
      SimVarFArr etof("etof", 0., 0.03, "nch", 100);
      SimVarFArr ecore("ecore", 0., 5., "nch", 100);
      SimVarFArr emce("emce", 0., 5., "nch", 100);
      SimVarFArr ecent("ecent", 0., 5., "nch", 100);
      SimVarFArr e9("e9", 0., 5., "nch", 100);
      SimVarFArr emcchi2("emcchi2", 0., 20., "nch", 100);
      SimVarShortArr twrhit("twrhit", 0., 20., "nch", 20);
      SimVarFArr emcdispy("emcdispy", 0., 20., "nch", 100);
      SimVarFArr emcdispz("emcdispz", 0., 20., "nch", 100);
      SimVarFArr prob("prob", 0., 1., "nch", 100);

      SimVarShortArr sect("sect", 0., 4., "nch", 4);
      SimVarShortArr ysect("ysect", 0, 48., "nch", 48);
      SimVarShortArr zsect("zsect", 0, 97., "nch", 97);
      
      SimVarShortArr n0("n0", 0., 20., "nch", 20);
      SimVarShortArr npe0("npe0", 0., 20., "nch", 20);
      SimVarShortArr n1("n1", 0., 20., "nch", 20);
      SimVarShortArr npe1("npe1", 0., 20., "nch", 20);
      SimVarShortArr n2("n2", 0., 20., "nch", 20);
      SimVarShortArr npe2("npe2", 0., 20., "nch", 20);
      SimVarShortArr n3("n3", 0., 20., "nch", 20);
      SimVarShortArr npe3("npe3", 0., 20., "nch", 20);
      SimVarFArr center_phi("center_phi", -2.*M_PI, 2.*M_PI, "nch", 100);
      SimVarFArr center_z("center_z", -500., 500., "nch", 100);
      SimVarFArr cross_phi("cross_phi", -2.*M_PI, 2*M_PI, "nch", 100);
      SimVarFArr cross_z("cross_z", -500., 500., "nch", 100);
      SimVarFArr disp("disp", 0., 10., "nch", 100);
      SimVarFArr chi2("chi2", 0., 200., "nch", 100);
   };
}

#endif /* SPLIT_SIM_TREE_HPP */
