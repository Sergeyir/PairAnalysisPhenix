// $HEADER$
//------------------------------------------------------------------------------------------------
//                        AnalyzeResonance function declaration
//------------------------------------------------------------------------------------------------
// AnalyzeResonance
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic macro used for evaluation of registration of resonance for different methods
 * from simulation output of event-like TTrees to processed histograms 
 * for further track pair from resonance registering correction evaluation
 **/
//------------------------------------------------------------------------------------------------

#ifndef ANALYZE_RESONANCE_HPP
#define ANALYZE_RESONANCE_HPP

#include "PBar.hpp"
#include "IOTools.hpp"
#include "Box.hpp"

#include "ThrObj.hpp"

#include "STrackFun.hpp"
#include "PTrackFun.hpp"
#include "ParAnalyzeResonance.hpp"
#include "Ident.hpp"
#include "EffTreeReader.hpp"

#include "ROOT/TTreeProcessorMT.hxx"

struct ThrContainer
{
   ThrObj<TH1F> distrOrigPT = ThrObj<TH1F>(
      "orig","orig", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair);

   ThrObj<TH2F> distrOrigPTVsRecPT = ThrObj<TH2F>(
      "orig_pt_vs_pt", "orig vs pt", 
      Par.ptNBins, Par.ptMin, Par.ptMax, 
      Par.ptNBins, Par.ptMinPair, Par.ptMaxPair);
   
   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> distrInvMNoPID;
   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> distrInvMNoPIDDecreasedAcceptance;
   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> distrInvMNoPIDIncreasedAcceptance;

   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> distrInvM1PID;
   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> distrInvM1PIDDecreasedAcceptance;
   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> distrInvM1PIDIncreasedAcceptance;

   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> distrInvM2PID;
   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> distrInvM2PIDDecreasedAcceptance;
   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> distrInvM2PIDIncreasedAcceptance;
   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> distrInvM2PIDDecreasedM2Eff;
   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> distrInvM2PIDIncreasedM2Eff;

   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> distrInvMTOF2PID;
   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> distrInvMTOF2PIDDecreasedAcceptance;
   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> distrInvMTOF2PIDIncreasedAcceptance;

   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> distrInvMEMC2PID;
   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> distrInvMEMC2PIDDecreasedAcceptance;
   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> distrInvMEMC2PIDIncreasedAcceptance;
   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> distrInvMEMC2PIDDecreasedM2Eff;
   std::array<std::unique_ptr<ThrObj<TH2F>>, Par.CType.size> distrInvMEMC2PIDIncreasedM2Eff;
};

struct ParticleContainer
{
   int size;
   
   double mass;
   double iter;
   int origId;
   int geantId;
   
   std::unique_ptr<double> embDCPC1;
   std::unique_ptr<double> embPC2;
   std::unique_ptr<double> embPC3;
   std::unique_ptr<double> embTOFe;
   std::unique_ptr<double> embTOFw;
   
   std::array<std::unique_ptr<double>, 4> embEMCale;
   std::array<std::unique_ptr<double>, 4> embEMCalw;

   std::unique_ptr<double> m2EffEMCale;
   std::unique_ptr<double> m2EffEMCalw;
   std::unique_ptr<double> m2EffSysEMCale;
   std::unique_ptr<double> m2EffSysEMCalw;
   
   double mom[3];
   std::array<int, 50> index;
   
   std::array<int, 50> idPC2;
   std::array<int, 50> idPC3;
   std::array<int, 50> idTOF;
   std::array<int, 50> idEMC;
   
   std::array<std::array<double, 50>, Par.CType.size> weightPC2;
   std::array<std::array<double, 50>, Par.CType.size> weightPC3;
   std::array<std::array<double, 50>, Par.CType.size> weightTOF;
   std::array<std::array<double, 50>, Par.CType.size> weightEMC;
   
   std::array<std::array<double, 50>, Par.CType.size> weightIdTOF;
   std::array<std::array<double, 50>, Par.CType.size> weightIdEMC;
   
   void ResetTrack(const int i)
   {
      index[i] = 0.;
      
      idPC2[i] = PartId.junk;
      idPC3[i] = PartId.junk;
      idTOF[i] = PartId.junk;
      idEMC[i] = PartId.junk;
      
      for (int j = 0; j < Par.CType.size; j++)
      {
         weightPC2[j][i] = 0.;
         weightPC3[j][i] = 0.;
         weightTOF[j][i] = 0.;
         weightEMC[j][i] = 0.;
         
         weightIdTOF[j][i] = 0.;
         weightIdEMC[j][i] = 0.;
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
