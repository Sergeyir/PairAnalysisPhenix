/** 
 *  @file   AnalyzeSingleTrack.hpp
 *  @brief  Contains declarations of functions and variables that are used for analysis of a single track from a trees acquired from the PHENIX simulation
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysisPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef ANALYZE_SINGLE_TRACK_CPP
#define ANALYZE_SINGLE_TRACK_CPP

#include "../include/AnalyzeSingleTrack.hpp"

// this namespace is only used so that documentation will not become a mess
// so there is no need to enforce the contents inside of it 
// being accessed only via the scope resolution operator in this file
using namespace AnalyzeSingleTrack;

void AnalyzeSingleTrack::AnalyzeConfiguration(ThrContainer &thrContainer, 
                                              const std::string& particleName, 
                                              const std::string& magneticFieldName, 
                                              const std::string &pTRangeName)
{ 

   std::string simInputFileName = 
      "data/SimTrees/" + runName + "/SingleTrack/" + 
      particleName + "_" + pTRangeName + magneticFieldName + ".root";
   std::string realDataFileName = 
      "data/Real/" + runName + "/SingleTrack/sum" + magneticFieldName + ".root";

   TFile simInputFile = TFile(simInputFileName.c_str());
   TFile realDataFile = TFile(realDataFileName.c_str());

   TH1F *origPTHist = static_cast<TH1F *>(simInputFile.Get("orig_pt"));

   // weight function for spectra
   std::unique_ptr<TF1> weightFunc;
  
   // normalization of the number of particles to the number of events
   // this normalization is needed to seamlessly merge 2 files with 
   // flat pT distribution with different ranges
   double eventNormWeight = 1.;
   if (reweightForSpectra)
   {
      // threshold is needed since there can be a little noise in the histogram
      const double origPTThreshold = 
         origPTHist->Integral()/static_cast<double>(origPTHist->GetXaxis()->GetNbins())/2.;
 
      // pT bounds of original pT distribution for spectra scale
      const double lowPTBound = 
         origPTHist->GetXaxis()->GetBinLowEdge(origPTHist->FindFirstBinAbove(origPTThreshold));
      const double upPTBound = 
         origPTHist->GetXaxis()->GetBinUpEdge(origPTHist->FindLastBinAbove(origPTThreshold));

      TH1F *centrHist = static_cast<TH1F *>(realDataFile.Get("centrality"));
      if (!centrHist) 
      {
         CppTools::PrintError("Histogram \"centrality\" does not exist in file" + 
                              static_cast<std::string>(realDataFile.GetName()));
      }

      eventNormWeight = origPTHist->Integral(origPTHist->GetXaxis()->FindBin(pTMin), 
                                            origPTHist->GetXaxis()->FindBin(pTMax))/
                        centrHist->Integral(1, centrHist->GetXaxis()->GetNbins());

      eventNormWeight *= (pTMax - pTMin)/(upPTBound-lowPTBound);

      InputYAMLReader inputYAMLSpectraFit("data/Parameters/SpectraFit/" + collisionSystemName + 
                                       "/" + particleName + ".yaml");
      inputYAMLSpectraFit.CheckStatus("spectra_fit");

      weightFunc = std::make_unique<TF1>
         ("weightFunc", inputYAMLSpectraFit["fit_function"].as<std::string>().c_str());

      for (unsigned int i = 0; i < inputYAMLSpectraFit["fit_parameters"].size(); i++)
      {
         weightFunc->SetParameter(i, inputYAMLSpectraFit["fit_parameters"][i].as<double>());
      }
   }

   ROOT::TTreeProcessorMT tp(simInputFileName.c_str());
  
   auto ProcessMP = [&](TTreeReader &reader)
   { 
      ThrContainerCopy histContainer = thrContainer.GetCopy();

      SimTreeReader STR(reader);

      while (reader.Next())
      { 
         numberOfCalls++;
         const double origPT = sqrt(pow(STR.mom_orig(0), 2) + pow(STR.mom_orig(1), 2));

         double eventWeight;
 
         if (reweightForSpectra) 
         {
            eventWeight = weightFunc->Eval(origPT)/eventNormWeight;
         }
         else eventWeight = 1.;
 
         histContainer.distrOrigPT->Fill(origPT, eventWeight);
 
         const double bbcz = STR.bbcz();
         if (fabs(bbcz) > 30.) continue;

         for(int i = 0; i < STR.nch(); i++)
         {
            const double the0 = STR.the0(i);
            const double pT = (STR.mom(i))*sin(the0);

            if (pT < pTMin) continue;
            if (IsQualityCut(STR.qual(i))) continue;

            const int charge = STR.charge(i);
            if (charge != -1 && charge != 1) continue;

            const int dcarm = STR.dcarm(i);

            const double zed = STR.zed(i);
            if (fabs(zed) > 75. && fabs(zed) < 3.) continue;

            if (!(fabs(the0) < 100. &&
               ((bbcz > 0. && ((bbcz - 250.*tan(the0 - 3.1416/2.)) > 2. ||
               (bbcz - 200.*tan(the0 - 3.1416/2)) < -2.)) ||
               (bbcz < 0. && ((bbcz - 250.*tan(the0 - 3.1416/2.))< -2. ||
               (bbcz - 200.*tan(the0 - 3.1416/2)) > 2.))))) continue;
 
            //end of basic cuts

            const double alpha = STR.alpha(i);
            const double phi = STR.phi(i);
            double board;

            if (phi > M_PI/2.) board = ((3.72402 - phi + 0.008047*cos(phi + 0.87851))/0.01963496);
            else board = ((0.573231 + phi - 0.0046 * cos(phi + 0.05721))/0.01963496);

            double alphaReweight = 1.;

            if (dcarm == 0)
            {
               if (zed >= 0) 
               {
                  histContainer.heatmapUnscaledDCe0->Fill(board, alpha, eventWeight);
                  if (reweightHeatmapsForAlpha) alphaReweight = 
                     alphaReweightDCe0->GetBinContent(alphaReweightDCe0->FindBin(alpha));
                  histContainer.heatmapDCe0->Fill(board, alpha, eventWeight*alphaReweight);
               }
               else 
               {
                  histContainer.heatmapUnscaledDCe1->Fill(board, alpha, eventWeight);
                  if (reweightHeatmapsForAlpha) alphaReweight = 
                     alphaReweightDCe1->GetBinContent(alphaReweightDCe1->FindBin(alpha));
                  histContainer.heatmapDCe1->Fill(board, alpha, eventWeight*alphaReweight);
               }
            }
            else
            {
               if (zed >= 0) 
               {
                  histContainer.heatmapUnscaledDCw0->Fill(board, alpha, eventWeight);
                  if (reweightHeatmapsForAlpha) alphaReweight = 
                     alphaReweightDCw0->GetBinContent(alphaReweightDCw0->FindBin(alpha));
                  histContainer.heatmapDCw0->Fill(board, alpha, eventWeight*alphaReweight);
               }
               else 
               {
                  histContainer.heatmapUnscaledDCw1->Fill(board, alpha, eventWeight);
                  if (reweightHeatmapsForAlpha) alphaReweight = 
                     alphaReweightDCw1->GetBinContent(alphaReweightDCw1->FindBin(alpha));
                  histContainer.heatmapDCw1->Fill(board, alpha, eventWeight*alphaReweight);
               }
            }

            bool isDCPC1TrackCut = dmCutter.IsDeadDC(dcarm, zed, board, alpha);

            const double ppc1phi = atan2(STR.ppc1y(i), STR.ppc1x(i));

            if (dcarm == 1) 
            {
               histContainer.heatmapPC1w->Fill(STR.ppc1z(i), ppc1phi, eventWeight*alphaReweight);
            }
            else
            {
               if (ppc1phi < 0) 
               {
                  histContainer.heatmapPC1e->Fill(STR.ppc1z(i), ppc1phi + 2.*M_PI, 
                                                  eventWeight*alphaReweight);
                  if (dmCutter.IsDeadPC1(dcarm, STR.ppc1z(i), ppc1phi + 2.*M_PI))
                  {
                     isDCPC1TrackCut = false;
                  }
               }
               else 
               {
                  histContainer.heatmapPC1e->Fill(STR.ppc1z(i), ppc1phi, 
                                                  eventWeight*alphaReweight);
                  if (dmCutter.IsDeadPC1(dcarm, STR.ppc1z(i), ppc1phi + 2.*M_PI))
                  {
                     isDCPC1TrackCut = false;
                  }
               }
            }

            if (!isDCPC1TrackCut) histContainer.distrOrigPTVsRecPT->Fill(origPT, pT, eventWeight);

            if (IsHit(STR.pc2dphi(i)))
            {
               if (charge == 1) 
               {
                  histContainer.distrDPhiVsPTPC2Pos->Fill(STR.pc2dphi(i), pT, eventWeight);
                  histContainer.distrDZVsPTPC2Pos->Fill(STR.pc2dz(i), pT, eventWeight);

                  histContainer.distrSDPhiVsPTPC2Pos->
                     Fill(simCalibrator.PC2SDPhi(STR.pc2dphi(i), pT, charge), 
                          pT, eventWeight);
                  histContainer.distrSDZVsPTPC2Pos->
                     Fill(simCalibrator.PC2SDZ(STR.pc2dz(i), pT, charge), 
                          pT, eventWeight);
               }
               else
               {
                  histContainer.distrDPhiVsPTPC2Neg->Fill(STR.pc2dphi(i), pT, eventWeight);
                  histContainer.distrDZVsPTPC2Neg->Fill(STR.pc2dz(i), pT, eventWeight);

                  histContainer.distrSDPhiVsPTPC2Neg->
                     Fill(simCalibrator.PC2SDPhi(STR.pc2dphi(i), pT, charge), 
                          pT, eventWeight);
                  histContainer.distrSDZVsPTPC2Neg->
                     Fill(simCalibrator.PC2SDZ(STR.pc2dz(i), pT, charge), 
                          pT, eventWeight);
               }

               if (IsMatch(simCalibrator.PC2SDPhi(STR.pc2dphi(i), pT, charge), 
                           simCalibrator.PC2SDZ(STR.pc2dz(i), pT, charge)))
               {
                  const double pc2z = STR.ppc2z(i) - STR.pc2dz(i);
                  const double pc2phi = atan2(STR.ppc2y(i), STR.ppc2x(i)) - STR.pc2dphi(i);
 
                  histContainer.heatmapPC2->Fill(pc2z, pc2phi, eventWeight*alphaReweight);
               }
            }

            if (IsHit(STR.pc3dphi(i)))
            {
               if (dcarm == 0)
               {
                  if (charge == 1) 
                  {
                     histContainer.distrDPhiVsPTPC3ePos->Fill(STR.pc3dphi(i), pT, eventWeight);
                     histContainer.distrDZVsPTPC3ePos->Fill(STR.pc3dz(i), pT, eventWeight);

                     histContainer.distrSDPhiVsPTPC3ePos->
                        Fill(simCalibrator.PC3SDPhi(STR.pc3dphi(i), pT, charge, dcarm), 
                             pT, eventWeight);
                     histContainer.distrSDZVsPTPC3ePos->
                        Fill(simCalibrator.PC3SDZ(STR.pc3dz(i), pT, charge, dcarm), 
                             pT, eventWeight);
                  }
                  else
                  {
                     histContainer.distrDPhiVsPTPC3eNeg->Fill(STR.pc3dphi(i), pT, eventWeight);
                     histContainer.distrDZVsPTPC3eNeg->Fill(STR.pc3dz(i), pT, eventWeight);

                     histContainer.distrSDPhiVsPTPC3eNeg->
                        Fill(simCalibrator.PC3SDPhi(STR.pc3dphi(i), pT, charge, dcarm), 
                             pT, eventWeight);
                     histContainer.distrSDZVsPTPC3eNeg->
                        Fill(simCalibrator.PC3SDZ(STR.pc3dz(i), pT, charge, dcarm), 
                             pT, eventWeight);
                  }
               }
               else
               {
                  if (charge == 1) 
                  {
                     histContainer.distrDPhiVsPTPC3wPos->Fill(STR.pc3dphi(i), pT, eventWeight);
                     histContainer.distrDZVsPTPC3wPos->Fill(STR.pc3dz(i), pT, eventWeight);

                     histContainer.distrSDPhiVsPTPC3wPos->
                        Fill(simCalibrator.PC3SDPhi(STR.pc3dphi(i), pT, charge, dcarm), 
                             pT, eventWeight);
                     histContainer.distrSDZVsPTPC3wPos->
                        Fill(simCalibrator.PC3SDZ(STR.pc3dz(i), pT, charge, dcarm), 
                             pT, eventWeight);
                  }
                  else
                  {
                     histContainer.distrDPhiVsPTPC3wNeg->Fill(STR.pc3dphi(i), pT, eventWeight);
                     histContainer.distrDZVsPTPC3wNeg->Fill(STR.pc3dz(i), pT, eventWeight);

                     histContainer.distrSDPhiVsPTPC3wNeg->
                        Fill(simCalibrator.PC3SDPhi(STR.pc3dphi(i), pT, charge, dcarm), 
                             pT, eventWeight);
                     histContainer.distrSDZVsPTPC3wNeg->
                        Fill(simCalibrator.PC3SDZ(STR.pc3dz(i), pT, charge, dcarm), 
                             pT, eventWeight);
                  }
               }
               if (IsMatch(simCalibrator.PC3SDPhi(STR.pc3dphi(i), pT, charge, dcarm), 
                           simCalibrator.PC3SDZ(STR.pc3dz(i), pT, charge, dcarm)))
               {
                  const double pc3z = STR.ppc3z(i) - STR.pc3dz(i);
                  double pc3phi = atan2(STR.ppc3y(i), STR.ppc3x(i) - STR.pc3dphi(i));
 
                  if (dcarm == 0) 
                  {
                     if (pc3phi < 0) 
                     {
                        histContainer.heatmapPC3e->Fill(pc3z, pc3phi + 2.*M_PI, 
                                                        eventWeight*alphaReweight);
                     }
                     else 
                     {
                        histContainer.heatmapPC3e->Fill(pc3z, pc3phi, eventWeight*alphaReweight);
                     }
                  }
                  else histContainer.heatmapPC3w->Fill(pc3z, pc3phi, eventWeight*alphaReweight);
               }
            }

            if (IsHit(STR.tofdz(i)))
            {
               if (charge == 1) 
               {
                  histContainer.distrDPhiVsPTTOFePos->Fill(STR.tofdphi(i), pT, eventWeight);
                  histContainer.distrDZVsPTTOFePos->Fill(STR.tofdz(i), pT, eventWeight);

                  histContainer.distrSDPhiVsPTTOFePos->
                     Fill(simCalibrator.TOFeSDPhi(STR.tofdphi(i), pT, charge),
                          pT, eventWeight);
                  histContainer.distrSDZVsPTTOFePos->
                     Fill(simCalibrator.TOFeSDZ(STR.tofdz(i), pT, charge),
                          pT, eventWeight);
               }
               else
               {
                  histContainer.distrDPhiVsPTTOFeNeg->Fill(STR.tofdphi(i), pT, eventWeight);
                  histContainer.distrDZVsPTTOFeNeg->Fill(STR.tofdz(i), pT, eventWeight);

                  histContainer.distrSDPhiVsPTTOFeNeg->
                     Fill(simCalibrator.TOFeSDPhi(STR.tofdphi(i), pT, charge), pT, eventWeight);
                  histContainer.distrSDZVsPTTOFeNeg->
                     Fill(simCalibrator.TOFeSDZ(STR.tofdz(i), pT, charge), pT, eventWeight);
               }

               const double beta = STR.pltof(i)/STR.ttof(i)/29.97;
               //const double eloss = 0.0014*pow(beta, -1.66);

               histContainer.distrSlatTOFe->Fill(STR.slat(i), eventWeight*alphaReweight);
               histContainer.distrELossTOFe->Fill(beta, STR.etof(i));
 
               if (IsMatch(simCalibrator.TOFeSDPhi(STR.tofdphi(i), pT, charge), 
                           simCalibrator.TOFeSDZ(STR.tofdz(i), pT, charge)))
               {
                  histContainer.heatmapTOFe->
                     Fill(STR.ptofy(i), STR.ptofz(i), STR.etof(i)*eventWeight*alphaReweight);
               }
            }
            else if (IsHit(STR.tofwdz(i)))
            {
               if (charge == 1) 
               {
                  histContainer.distrDPhiVsPTTOFwPos->Fill(STR.tofwdphi(i), pT, eventWeight);
                  histContainer.distrDZVsPTTOFwPos->Fill(STR.tofwdz(i), pT, eventWeight);

                  histContainer.distrSDPhiVsPTTOFwPos->
                     Fill(simCalibrator.TOFwSDPhi(STR.tofwdphi(i), pT, charge), 
                          pT, eventWeight);
                  histContainer.distrSDZVsPTTOFwPos->
                     Fill(simCalibrator.TOFwSDZ(STR.tofwdz(i), pT, charge), 
                          pT, eventWeight);
               }
               else
               {
                  histContainer.distrDPhiVsPTTOFwNeg->Fill(STR.tofwdphi(i), pT, eventWeight);
                  histContainer.distrDZVsPTTOFwNeg->Fill(STR.tofwdz(i), pT, eventWeight);

                  histContainer.distrSDPhiVsPTTOFwNeg->
                     Fill(simCalibrator.TOFwSDPhi(STR.tofwdphi(i), pT, charge),
                          pT, eventWeight);
                  histContainer.distrSDZVsPTTOFwNeg->
                     Fill(simCalibrator.TOFwSDZ(STR.tofwdz(i), pT, charge),
                          pT, eventWeight);
               }
               if (IsMatch(simCalibrator.TOFwSDPhi(STR.tofwdphi(i), pT, charge), 
                           simCalibrator.TOFwSDZ(STR.tofwdz(i), pT, charge)))
               {
                  histContainer.distrStripTOFw->Fill(STR.striptofw(i), eventWeight*alphaReweight);
                  if (STR.ptofwy(i) < 100.) 
                  {
                     histContainer.heatmapTOFw0->Fill(STR.ptofwy(i), STR.ptofwz(i), 
                                                      eventWeight*alphaReweight);
                  }
                  else 
                  {
                     histContainer.heatmapTOFw1->Fill(STR.ptofwy(i), STR.ptofwz(i), 
                                                      eventWeight*alphaReweight);
                  }
               }
            }

            if (IsHit(STR.emcdz(i)))
            {
               if (dcarm == 0)
               {
                  if (charge == 1) 
                  {
                     histContainer.distrDPhiVsPTEMCalePos[STR.sect(i)]->
                        Fill(STR.emcdphi(i), pT, eventWeight);
                     histContainer.distrDZVsPTEMCalePos[STR.sect(i)]->
                        Fill(STR.emcdz(i), pT, eventWeight);

                     histContainer.distrSDPhiVsPTEMCalePos[STR.sect(i)]->
                        Fill(simCalibrator.EMCalSDPhi(STR.emcdphi(i), pT, charge, 
                                                      dcarm, STR.sect(i)),
                             pT, eventWeight);
                     histContainer.distrSDZVsPTEMCalePos[STR.sect(i)]->
                        Fill(simCalibrator.EMCalSDZ(STR.emcdz(i), pT, charge, 
                                                    dcarm, STR.sect(i)),
                             pT, eventWeight);
                  }
                  else
                  {
                     histContainer.distrDPhiVsPTEMCaleNeg[STR.sect(i)]->
                        Fill(STR.emcdphi(i), pT, eventWeight);
                     histContainer.distrDZVsPTEMCaleNeg[STR.sect(i)]->
                        Fill(STR.emcdz(i), pT, eventWeight);

                     histContainer.distrSDPhiVsPTEMCaleNeg[STR.sect(i)]->
                        Fill(simCalibrator.EMCalSDPhi(STR.emcdphi(i), pT, charge, 
                                                      dcarm, STR.sect(i)),
                             pT, eventWeight);
                     histContainer.distrSDZVsPTEMCaleNeg[STR.sect(i)]->
                        Fill(simCalibrator.EMCalSDZ(STR.emcdz(i), pT, charge, 
                                                    dcarm, STR.sect(i)),
                             pT, eventWeight);
                  }
               }
               else
               {
                  if (charge == 1) 
                  {
                     histContainer.distrDPhiVsPTEMCalwPos[STR.sect(i)]->
                        Fill(STR.emcdphi(i), pT, eventWeight);
                     histContainer.distrDZVsPTEMCalwPos[STR.sect(i)]->
                        Fill(STR.emcdz(i), pT, eventWeight);

                     histContainer.distrSDPhiVsPTEMCalwPos[STR.sect(i)]->
                        Fill(simCalibrator.EMCalSDPhi(STR.emcdphi(i), pT, charge, 
                                                      dcarm, STR.sect(i)),
                             pT, eventWeight);
                     histContainer.distrSDZVsPTEMCalwPos[STR.sect(i)]->
                        Fill(simCalibrator.EMCalSDZ(STR.emcdz(i), pT, charge, 
                                                    dcarm, STR.sect(i)),
                             pT, eventWeight);
                  }
                  else
                  {
                     histContainer.distrDPhiVsPTEMCalwNeg[STR.sect(i)]->
                        Fill(STR.emcdphi(i), pT, eventWeight);
                     histContainer.distrDZVsPTEMCalwNeg[STR.sect(i)]->
                        Fill(STR.emcdz(i), pT, eventWeight);

                     histContainer.distrSDPhiVsPTEMCalwNeg[STR.sect(i)]->
                        Fill(simCalibrator.EMCalSDPhi(STR.emcdphi(i), pT, charge, 
                                                      dcarm, STR.sect(i)),
                             pT, eventWeight);
                     histContainer.distrSDZVsPTEMCalwNeg[STR.sect(i)]->
                        Fill(simCalibrator.EMCalSDZ(STR.emcdz(i), pT, charge, 
                                                    dcarm, STR.sect(i)),
                             pT, eventWeight);
                  }
               }

               if (IsMatch(simCalibrator.EMCalSDPhi(STR.emcdphi(i), pT, charge, 
                                                    dcarm, STR.sect(i)), 
                           simCalibrator.EMCalSDZ(STR.emcdz(i), pT, charge, 
                                                  dcarm, STR.sect(i))))
               { 

                  if (dcarm == 0) 
                  {
                     histContainer.distrECoreVsPTEMCale[STR.sect(i)]->Fill(pT, STR.ecore(i));
                     histContainer.heatmapEMCale[STR.sect(i)]->
                        Fill(STR.ysect(i), STR.zsect(i), STR.ecore(i)*eventWeight*alphaReweight);
                  }
                  else 
                  {
                     histContainer.distrECoreVsPTEMCalw[STR.sect(i)]->Fill(pT, STR.ecore(i));
                     histContainer.heatmapEMCalw[STR.sect(i)]->
                        Fill(STR.ysect(i), STR.zsect(i), STR.ecore(i)*eventWeight*alphaReweight);
                  }
               }
            }
         }
      }
   };

   tp.Process(ProcessMP);
}

int main(int argc, char **argv)
{
   if (argc < 2 || argc > 3) 
   {
      std::string errMsg = "Expected 1-2 parameters while " + std::to_string(argc) + " ";
      errMsg += "parameter(s) were provided \n Usage: bin/AnalyzeSingleTrack ";
      errMsg += "inputYAMLName numberOfThreads=std::thread::hardware_concurrency()";
      CppTools::PrintError(errMsg);
   }
 
   CppTools::CheckInputFile(argv[1]);
 
   if (argc == 2) numberOfThreads = std::thread::hardware_concurrency();
   else numberOfThreads = std::stoi(argv[2]);

   ROOT::EnableImplicitMT(numberOfThreads);

   inputYAMLSim.OpenFile(argv[1], "single_track_sim");
   inputYAMLSim.CheckStatus("single_track_sim");

   runName = inputYAMLSim["run_name"].as<std::string>();

   inputYAMLMain.OpenFile("input/" + runName + "/main.yaml");
   inputYAMLMain.CheckStatus("main");
 
   collisionSystemName = inputYAMLMain["collision_system_name"].as<std::string>();

   outputDir = "data/PostSim/" + runName + "/SingleTrack/";
   system(("mkdir -p " + outputDir).c_str());

   pTMin = inputYAMLSim["pt_min"].as<double>();
   pTMax = inputYAMLSim["pt_max"].as<double>();

   reweightForSpectra = inputYAMLSim["reweight_for_spectra"].as<bool>();
   reweightHeatmapsForAlpha = inputYAMLSim["reweight_heatmaps_for_alpha"].as<bool>();

   dmCutter.Initialize(runName, inputYAMLMain["detectors_configuration"].as<std::string>());
   simCalibrator.Initialize(runName, inputYAMLMain["detectors_configuration"].as<std::string>());

   if (reweightForSpectra)
   {
      for (const auto& particle : inputYAMLSim["particles"])
      {
         CppTools::CheckInputFile("data/Parameters/SpectraFit/" + collisionSystemName + 
                                  "/" + particle["name"].as<std::string>() + ".yaml");
      }
   }
 
   for (const auto& magneticField : inputYAMLMain["magnetic_field_configurations"])
   {
      CppTools::CheckInputFile("data/Real/" + runName + "/SingleTrack/sum" + 
                               magneticField["name"].as<std::string>() + ".root");

      for (const auto& particle : inputYAMLSim["particles"])
      {
         for (const auto& pTRange : inputYAMLSim["pt_ranges"])
         {
            const std::string simInputFileName = 
               "data/SimTrees/" + runName + "/SingleTrack/" + 
               particle["name"].as<std::string>() + "_" + pTRange["name"].as<std::string>() + 
               magneticField["name"].as<std::string>() + ".root";

            CppTools::CheckInputFile(simInputFileName);

            const unsigned long currentConfigurationNumberOfEvents = static_cast<unsigned long>
               ((static_cast<TTree *>(TFile::Open(simInputFileName.c_str())->
                                      Get("Tree"))->GetEntries()));
            if (currentConfigurationNumberOfEvents <= 0)
            {
               CppTools::PrintError("Number of events is equal or less than 0 in file " + 
                                    simInputFileName);
            }
            numberOfEvents += currentConfigurationNumberOfEvents;
         }
      }
   }

   if (reweightHeatmapsForAlpha)
   {
      // real data file in which real data DC heatmaps are stored
      std::string realDataInputFileName = "data/Real/" + runName + "/SingleTrack/sum.root";
      // file in which simulated unscaled DC heatmaps are stored
      std::string postSimInputFileName = "data/PostSim/" + runName + "/SingleTrack/all.root";

      CppTools::CheckInputFile(realDataInputFileName);

      if (CppTools::CheckInputFile(postSimInputFileName, false)) 
      {
         TFile realDataInputFile(realDataInputFileName.c_str());
         TFile postSimInputFile(postSimInputFileName.c_str());

         TH2F *realDataDCe0 = GetHistogramFromFile(realDataInputFile, "Heatmap: DCe, zDC>=0");
         TH2F *realDataDCe1 = GetHistogramFromFile(realDataInputFile, "Heatmap: DCe, zDC<0");
         TH2F *realDataDCw0 = GetHistogramFromFile(realDataInputFile, "Heatmap: DCw, zDC>=0");
         TH2F *realDataDCw1 = GetHistogramFromFile(realDataInputFile, "Heatmap: DCw, zDC<0");

         TH2F *simUnscaledDCe0 = 
            GetHistogramFromFile(postSimInputFile, "Unscaled heatmap: DCe, zDC>=0");
         TH2F *simUnscaledDCe1 = 
            GetHistogramFromFile(postSimInputFile, "Unscaled heatmap: DCe, zDC<0");
         TH2F *simUnscaledDCw0 = 
            GetHistogramFromFile(postSimInputFile, "Unscaled heatmap: DCw, zDC>=0");
         TH2F *simUnscaledDCw1 = 
            GetHistogramFromFile(postSimInputFile, "Unscaled heatmap: DCw, zDC<0");
 
         CheckHistsAxis(realDataDCe0, simUnscaledDCe0);
         CheckHistsAxis(realDataDCe1, simUnscaledDCe1);
         CheckHistsAxis(realDataDCw0, simUnscaledDCw0);
         CheckHistsAxis(realDataDCw1, simUnscaledDCw1);
 
         for (int i = 1; i <= realDataDCe0->GetXaxis()->GetNbins(); i++)
         {
            const double board = realDataDCe0->GetXaxis()->GetBinCenter(i);
            for (int j = 1; j <= realDataDCe0->GetYaxis()->GetNbins(); j++)
            {
               const double alpha = realDataDCe0->GetYaxis()->GetBinCenter(j);
               if (dmCutter.IsDeadDC(1, 1., board, alpha))
               {
                  realDataDCe0->SetBinContent(i, j, 0.);
                  simUnscaledDCe0->SetBinContent(i, j, 0.);
               }
            }
         }
         for (int i = 1; i <= realDataDCe1->GetXaxis()->GetNbins(); i++)
         {
            const double board = realDataDCe1->GetXaxis()->GetBinCenter(i);
            for (int j = 1; j < realDataDCe1->GetYaxis()->GetNbins(); j++)
            {
               const double alpha = realDataDCe1->GetYaxis()->GetBinCenter(j);
               if (dmCutter.IsDeadDC(1, -1., board, alpha))
               {
                  realDataDCe1->SetBinContent(i, j, 0.);
                  simUnscaledDCe1->SetBinContent(i, j, 0.);
               }
            }
         }
         for (int i = 1; i <= realDataDCw0->GetXaxis()->GetNbins(); i++)
         {
            const double board = realDataDCe1->GetXaxis()->GetBinCenter(i);
            for (int j = 1; j < realDataDCw0->GetYaxis()->GetNbins(); j++)
            {
               const double alpha = realDataDCe1->GetYaxis()->GetBinCenter(j);
               if (dmCutter.IsDeadDC(0, 1., board, alpha))
               {
                  realDataDCw0->SetBinContent(i, j, 0.);
                  simUnscaledDCw0->SetBinContent(i, j, 0.);
               }
            }
         }
         for (int i = 1; i <= realDataDCw1->GetXaxis()->GetNbins(); i++)
         {
            const double board = realDataDCw1->GetXaxis()->GetBinCenter(i);
            for (int j = 1; j < realDataDCw1->GetYaxis()->GetNbins(); j++)
            {
               const double alpha = realDataDCw1->GetYaxis()->GetBinCenter(j);
               if (dmCutter.IsDeadDC(0, -1., board, alpha))
               {
                  realDataDCw1->SetBinContent(i, j, 0.);
                  simUnscaledDCw1->SetBinContent(i, j, 0.);
               }
            }
         }

         alphaReweightDCe0 = 
            static_cast<TH1F *>(realDataDCe0->ProjectionY("Alpha reweight: DCe, zDC>=0",
                                1, realDataDCe0->GetXaxis()->GetNbins())->Clone());
         alphaReweightDCe1 = 
            static_cast<TH1F *>(realDataDCe1->ProjectionY("Alpha reweight: DCe, zDC>=0",
                                1, realDataDCe1->GetXaxis()->GetNbins())->Clone());
         alphaReweightDCw0 = 
            static_cast<TH1F *>(realDataDCw0->ProjectionY("Alpha reweight: DCe, zDC>=0",
                                1, realDataDCw0->GetXaxis()->GetNbins())->Clone());
         alphaReweightDCw1 = 
            static_cast<TH1F *>(realDataDCw1->ProjectionY("Alpha reweight: DCe, zDC>=0",
                                1, realDataDCw1->GetXaxis()->GetNbins())->Clone());

         alphaReweightDCe0->Scale(simUnscaledDCe0->Integral()/alphaReweightDCe0->Integral());
         alphaReweightDCe1->Scale(simUnscaledDCe1->Integral()/alphaReweightDCe1->Integral());
         alphaReweightDCw0->Scale(simUnscaledDCw0->Integral()/alphaReweightDCw0->Integral());
         alphaReweightDCw1->Scale(simUnscaledDCw1->Integral()/alphaReweightDCw1->Integral());

         alphaReweightDCe0->Divide(simUnscaledDCe0->ProjectionY("DCe0 proj",
                                   1, simUnscaledDCe0->GetXaxis()->GetNbins()));
         alphaReweightDCe1->Divide(simUnscaledDCe1->ProjectionY("DCe1 proj",
                                   1, simUnscaledDCe1->GetXaxis()->GetNbins()));
         alphaReweightDCw0->Divide(simUnscaledDCw0->ProjectionY("DCw0 proj",
                                   1, simUnscaledDCw0->GetXaxis()->GetNbins()));
         alphaReweightDCw1->Divide(simUnscaledDCw1->ProjectionY("DCw1 proj",
                                   1, simUnscaledDCw1->GetXaxis()->GetNbins()));

         // capping alpha reweight since experiments extends to higher pT than simulation
         // capped values do not affect the simulation since they are statistically insufficient
         // capping is only needed to exclude very big weights in some points in the DC map
         // that are not used for better visibility
         // also setting 0 bins into 1 since some of the areas might have been cut 
         // when they are assumed to be bad/dead
         for (int i = 0; i < alphaReweightDCe0->GetXaxis()->GetNbins(); i++)
         {
            if (alphaReweightDCe0->GetBinContent(i) > 100.)
            {
               alphaReweightDCe0->SetBinContent(i, 100.);
            }
            if (alphaReweightDCe0->GetBinContent(i) < 1e-15)
            {
               alphaReweightDCe0->SetBinContent(i, 1.);
            }
         }
         for (int i = 0; i < alphaReweightDCe1->GetXaxis()->GetNbins(); i++)
         {
            if (alphaReweightDCe1->GetBinContent(i) > 100.)
            {
               alphaReweightDCe1->SetBinContent(i, 100.);
            }
            if (alphaReweightDCe1->GetBinContent(i) < 1e-15)
            {
               alphaReweightDCe1->SetBinContent(i, 1.);
            }
         }
         for (int i = 0; i < alphaReweightDCw0->GetXaxis()->GetNbins(); i++)
         {
            if (alphaReweightDCw0->GetBinContent(i) > 100.)
            {
               alphaReweightDCw0->SetBinContent(i, 100.);
            }
            if (alphaReweightDCw0->GetBinContent(i) < 1e-15)
            {
               alphaReweightDCw0->SetBinContent(i, 1.);
            }
         }
         for (int i = 0; i < alphaReweightDCw1->GetXaxis()->GetNbins(); i++)
         {
            if (alphaReweightDCw1->GetBinContent(i) > 100.)
            {
               alphaReweightDCw1->SetBinContent(i, 100.);
            }
            if (alphaReweightDCw1->GetBinContent(i) < 1e-15)
            {
               alphaReweightDCw1->SetBinContent(i, 1.);
            }
         }

         std::string alphaReweightOutputFileName = 
            "data/PostSim/" + runName + "/SingleTrack/alpha_reweight.root";

         TFile alphaReweightOutputFile = TFile(alphaReweightOutputFileName.c_str(), "RECREATE");

         alphaReweightOutputFile.cd();
 
         alphaReweightDCe0->Write();
         alphaReweightDCe1->Write();
         alphaReweightDCw0->Write();
         alphaReweightDCw1->Write();

         alphaReweightDCe0->SetDirectory(0);
         alphaReweightDCe1->SetDirectory(0);
         alphaReweightDCw0->SetDirectory(0);
         alphaReweightDCw1->SetDirectory(0);
 
         alphaReweightOutputFile.Close();
         CppTools::PrintInfo("File " + alphaReweightOutputFileName + " was written");
      }
      else 
      {
         CppTools::PrintInfo("Alpha reweight is now disabled. The missing file will be created"\
                             "after the current process is finished. Try again after the"\
                             "missing file is written.");
         reweightHeatmapsForAlpha = false;
      }
   }

   CppTools::PrintInfo("Clearing output directory: " + outputDir);
   system(("find " + outputDir + " ! -name 'all.root' ! -name 'alpha_reweight.root' -type f -exec rm -f {} +").c_str());

   std::vector<std::string> particleList;
   for (const auto& particle : inputYAMLSim["particles"])
   {
      particleList.emplace_back(particle["name"].as<std::string>());
   }

   std::vector<std::string> magneticFieldsList;
   for (const auto& magneticField : inputYAMLMain["magnetic_field_configurations"])
   {
      magneticFieldsList.emplace_back(magneticField["name"].as<std::string>());
   }

   std::vector<std::string> pTRangesList;
   for (const auto& pTRange : inputYAMLSim["pt_ranges"])
   {
      pTRangesList.emplace_back(pTRange["name"].as<std::string>());
   }

   CppTools::Box box{"Parameters"};
 
   box.AddEntry("Run name", runName);
   box.AddEntry("Particles", particleList);
   if (magneticFieldsList.size() == 1 && magneticFieldsList[0] == "")
   {
      box.AddEntry("Magnetic field", "run default");
   }
   else box.AddEntry("Magnetic fields", magneticFieldsList);
   box.AddEntry("pT ranges", pTRangesList);
   box.AddEntry("Charged track minimum pT, GeV", pTMin);
   box.AddEntry("Charged track maximum pT, GeV", pTMax);
   box.AddEntry("Reweight for pT spectra", reweightForSpectra);
   box.AddEntry("Reweight for alpha", reweightHeatmapsForAlpha);
   box.AddEntry("Number of threads", numberOfThreads);
   box.AddEntry("Number of events to be analyzed, 1e6", 
                static_cast<double>(numberOfEvents)/1e6, 3);
   box.Print();

   bool isProcessFinished = false;

   auto pBarCall = [&]()
   {
      ProgressBar pBar{"BLOCK"};
      while (!isProcessFinished)
      {
         pBar.Print(static_cast<double>(numberOfCalls)/
                    static_cast<double>(numberOfEvents));
         std::this_thread::sleep_for(std::chrono::milliseconds(100));
      }
      pBar.Print(1.);
   };
 
   std::thread pBarThread(pBarCall);
 
   for (const auto& particle : inputYAMLSim["particles"])
   {
      for (const auto& magneticField : inputYAMLMain["magnetic_field_configurations"])
      {
         ThrContainer thrContainer;
         for (const auto& pTRange : inputYAMLSim["pt_ranges"])
         {
            AnalyzeConfiguration(thrContainer, 
                                 particle["name"].as<std::string>(), 
                                 magneticField["name"].as<std::string>(), 
                                 pTRange["name"].as<std::string>());
         }
         // writing the result
         std::string outputFileName = "data/PostSim/" + runName + "/SingleTrack/" + 
                                      particle["name"].as<std::string>();
         if (magneticField["name"].as<std::string>() != "") 
         {
            outputFileName += "magf" + magneticField["name"].as<std::string>();
         }
         outputFileName += ".root";
         ROOTTools::ThrObjHolder::Write(outputFileName);
      }
   }

   isProcessFinished = true;
   pBarThread.join();

   CppTools::PrintInfo("Merging output files into one");

   std::string haddCommand = "hadd -f " + outputDir + "all.root ";
   for (const auto& particle : inputYAMLSim["particles"])
   {
      haddCommand += outputDir + particle["name"].as<std::string>() + ".root ";
   }
   system(haddCommand.c_str());

   CppTools::PrintInfo("If you see an error below it's a ROOT thing (ROOT handles weirdly its "\
                       "objects when they are already deleted with unique_ptr); no worries it's "\
                       "just an error of destructor, all calculations are finished and all files "\
                       "are already written");

   return 0;
}

TH2F *AnalyzeSingleTrack::GetHistogramFromFile(TFile& file, const std::string& histName)
{
   /// why the f doesn't ROOT allow to read the thing from read only TFile when its const
   TH2F *hist = static_cast<TH2F *>(file.Get(histName.c_str()));
   if (!hist) CppTools::PrintError("Histogram " + histName + " was not found in file " + 
                                   static_cast<std::string>(file.GetName()));
   // ROOT thing: making histograms not be automatically deleted when the file is closed
   hist->SetDirectory(0);
   return hist;
}

void AnalyzeSingleTrack::CheckHistsAxis(const TH2F *hist1, const TH2F *hist2)
{
   if (hist1->GetXaxis()->GetNbins() != hist2->GetXaxis()->GetNbins())
   {
      CppTools::PrintError("Histograms \"" + static_cast<std::string>(hist1->GetName()) + 
                           "\" and \"" + static_cast<std::string>(hist2->GetName()) + 
                           "\" have different number of bins on X axis");
   }
   if (hist1->GetYaxis()->GetNbins() != hist2->GetYaxis()->GetNbins())
   {
      CppTools::PrintError("Histograms \"" + static_cast<std::string>(hist1->GetName()) + 
                           "\" and \"" + static_cast<std::string>(hist2->GetName()) + 
                           "\" have different number of bins on Y axis");
   }
   if (fabs(hist1->GetXaxis()->GetBinLowEdge(1) - 
            hist2->GetXaxis()->GetBinLowEdge(1)) > 1e-15 ||
       fabs(hist1->GetXaxis()->GetBinLowEdge(hist1->GetXaxis()->GetNbins()) - 
            hist2->GetXaxis()->GetBinLowEdge(hist2->GetXaxis()->GetNbins())) > 1e-15)
   {
      CppTools::PrintError("Histograms \"" + static_cast<std::string>(hist1->GetName()) + 
                           "\" and \"" + static_cast<std::string>(hist2->GetName()) + 
                           "\" have different ranges on X axis");
   }
   if (fabs(hist1->GetYaxis()->GetBinLowEdge(1) - 
            hist2->GetYaxis()->GetBinLowEdge(1)) > 1e-15 ||
       fabs(hist1->GetYaxis()->GetBinLowEdge(hist1->GetYaxis()->GetNbins()) - 
            hist2->GetYaxis()->GetBinLowEdge(hist2->GetYaxis()->GetNbins())) > 1e-15)
   {
      CppTools::PrintError("Histograms \"" + static_cast<std::string>(hist1->GetName()) + 
                           "\" and \"" + static_cast<std::string>(hist2->GetName()) + 
                           "\" have different ranges on Y axis");
   }
}

ThrContainerCopy AnalyzeSingleTrack::ThrContainer::GetCopy()
{
   ThrContainerCopy copy;

   copy.distrOrigPT = distrOrigPT.Get();
   copy.distrOrigPTVsRecPT = distrOrigPTVsRecPT.Get();
   copy.heatmapUnscaledDCe0 = heatmapUnscaledDCe0.Get();
   copy.heatmapUnscaledDCe1 = heatmapUnscaledDCe1.Get();
   copy.heatmapUnscaledDCw0 = heatmapUnscaledDCw0.Get();
   copy.heatmapUnscaledDCw1 = heatmapUnscaledDCw1.Get();
   copy.heatmapDCe0 = heatmapDCe0.Get();
   copy.heatmapDCe1 = heatmapDCe1.Get();
   copy.heatmapDCw0 = heatmapDCw0.Get();
   copy.heatmapDCw1 = heatmapDCw1.Get();
   copy.heatmapPC1e = heatmapPC1e.Get();
   copy.heatmapPC1w = heatmapPC1w.Get();
   copy.heatmapPC2 = heatmapPC2.Get();
   copy.heatmapPC3e = heatmapPC3e.Get();
   copy.heatmapPC3w = heatmapPC3w.Get();
   copy.heatmapTOFe = heatmapTOFe.Get();
   copy.heatmapTOFw0 = heatmapTOFw0.Get();
   copy.heatmapTOFw1 = heatmapTOFw1.Get();
   copy.distrStripTOFw = distrStripTOFw.Get();
   copy.distrSlatTOFe = distrSlatTOFe.Get();
   copy.distrELossTOFe = distrELossTOFe.Get();
   copy.distrDPhiVsPTPC2Pos = distrDPhiVsPTPC2Pos.Get();
   copy.distrDZVsPTPC2Pos = distrDZVsPTPC2Pos.Get();
   copy.distrDPhiVsPTPC2Neg = distrDPhiVsPTPC2Neg.Get();
   copy.distrDZVsPTPC2Neg = distrDZVsPTPC2Neg.Get();
   copy.distrDPhiVsPTPC3ePos = distrDPhiVsPTPC3ePos.Get();
   copy.distrDZVsPTPC3ePos = distrDZVsPTPC3ePos.Get();
   copy.distrDPhiVsPTPC3eNeg = distrDPhiVsPTPC3eNeg.Get();
   copy.distrDZVsPTPC3eNeg = distrDZVsPTPC3eNeg.Get();
   copy.distrDPhiVsPTPC3wPos = distrDPhiVsPTPC3wPos.Get();
   copy.distrDZVsPTPC3wPos = distrDZVsPTPC3wPos.Get();
   copy.distrDPhiVsPTPC3wNeg = distrDPhiVsPTPC3wNeg.Get();
   copy.distrDZVsPTPC3wNeg = distrDZVsPTPC3wNeg.Get();
   copy.distrDPhiVsPTTOFePos = distrDPhiVsPTTOFePos.Get();
   copy.distrDZVsPTTOFePos = distrDZVsPTTOFePos.Get();
   copy.distrDPhiVsPTTOFeNeg = distrDPhiVsPTTOFeNeg.Get();
   copy.distrDZVsPTTOFeNeg = distrDZVsPTTOFeNeg.Get();
   copy.distrDPhiVsPTTOFwPos = distrDPhiVsPTTOFwPos.Get();
   copy.distrDZVsPTTOFwPos = distrDZVsPTTOFwPos.Get();
   copy.distrDPhiVsPTTOFwNeg = distrDPhiVsPTTOFwNeg.Get();
   copy.distrDZVsPTTOFwNeg = distrDZVsPTTOFwNeg.Get();
   copy.distrSDPhiVsPTPC2Pos = distrSDPhiVsPTPC2Pos.Get();
   copy.distrSDZVsPTPC2Pos = distrSDZVsPTPC2Pos.Get();
   copy.distrSDPhiVsPTPC2Neg = distrSDPhiVsPTPC2Neg.Get();
   copy.distrSDZVsPTPC2Neg = distrSDZVsPTPC2Neg.Get();
   copy.distrSDPhiVsPTPC3ePos = distrSDPhiVsPTPC3ePos.Get();
   copy.distrSDZVsPTPC3ePos = distrSDZVsPTPC3ePos.Get();
   copy.distrSDPhiVsPTPC3eNeg = distrSDPhiVsPTPC3eNeg.Get();
   copy.distrSDZVsPTPC3eNeg = distrSDZVsPTPC3eNeg.Get();
   copy.distrSDPhiVsPTPC3wPos = distrSDPhiVsPTPC3wPos.Get();
   copy.distrSDZVsPTPC3wPos = distrSDZVsPTPC3wPos.Get();
   copy.distrSDPhiVsPTPC3wNeg = distrSDPhiVsPTPC3wNeg.Get();
   copy.distrSDZVsPTPC3wNeg = distrSDZVsPTPC3wNeg.Get();
   copy.distrSDPhiVsPTTOFePos = distrSDPhiVsPTTOFePos.Get();
   copy.distrSDZVsPTTOFePos = distrSDZVsPTTOFePos.Get();
   copy.distrSDPhiVsPTTOFeNeg = distrSDPhiVsPTTOFeNeg.Get();
   copy.distrSDZVsPTTOFeNeg = distrSDZVsPTTOFeNeg.Get();
   copy.distrSDPhiVsPTTOFwPos = distrSDPhiVsPTTOFwPos.Get();
   copy.distrSDZVsPTTOFwPos = distrSDZVsPTTOFwPos.Get();
   copy.distrSDPhiVsPTTOFwNeg = distrSDPhiVsPTTOFwNeg.Get();
   copy.distrSDZVsPTTOFwNeg = distrSDZVsPTTOFwNeg.Get();

   // iterating over EMCal sectors
   for (int i = 0; i < 4; i++)
   {
      copy.heatmapEMCale[i] = heatmapEMCale[i].Get();
      copy.heatmapEMCalw[i] = heatmapEMCalw[i].Get();
      copy.distrECoreVsPTEMCale[i] = distrECoreVsPTEMCale[i].Get();
      copy.distrECoreVsPTEMCalw[i] = distrECoreVsPTEMCalw[i].Get();
      copy.distrDPhiVsPTEMCalePos[i] = distrDPhiVsPTEMCalePos[i].Get();
      copy.distrDZVsPTEMCalePos[i] = distrDZVsPTEMCalePos[i].Get();
      copy.distrDPhiVsPTEMCaleNeg[i] = distrDPhiVsPTEMCaleNeg[i].Get();
      copy.distrDZVsPTEMCaleNeg[i] = distrDZVsPTEMCaleNeg[i].Get();
      copy.distrDPhiVsPTEMCalwPos[i] = distrDPhiVsPTEMCalwPos[i].Get();
      copy.distrDZVsPTEMCalwPos[i] = distrDZVsPTEMCalwPos[i].Get();
      copy.distrDPhiVsPTEMCalwNeg[i] = distrDPhiVsPTEMCalwNeg[i].Get();
      copy.distrDZVsPTEMCalwNeg[i] = distrDZVsPTEMCalwNeg[i].Get();
      copy.distrSDPhiVsPTEMCalePos[i] = distrSDPhiVsPTEMCalePos[i].Get();
      copy.distrSDZVsPTEMCalePos[i] = distrSDZVsPTEMCalePos[i].Get();
      copy.distrSDPhiVsPTEMCaleNeg[i] = distrSDPhiVsPTEMCaleNeg[i].Get();
      copy.distrSDZVsPTEMCaleNeg[i] = distrSDZVsPTEMCaleNeg[i].Get();
      copy.distrSDPhiVsPTEMCalwPos[i] = distrSDPhiVsPTEMCalwPos[i].Get();
      copy.distrSDZVsPTEMCalwPos[i] = distrSDZVsPTEMCalwPos[i].Get();
      copy.distrSDPhiVsPTEMCalwNeg[i] = distrSDPhiVsPTEMCalwNeg[i].Get();
      copy.distrSDZVsPTEMCalwNeg[i] = distrSDZVsPTEMCalwNeg[i].Get();
   }
   return copy;
}

#endif /* ANALYZE_SINGLE_TRACK_CPP */
