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
                                              const short particleGeantId, 
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

      SimTreeReader simCNT(reader);

      while (reader.Next())
      { 
         numberOfCalls++;
         const double origPT = sqrt(pow(simCNT.mom_orig(0), 2) + pow(simCNT.mom_orig(1), 2));

         double eventWeight;
 
         if (reweightForSpectra) 
         {
            eventWeight = weightFunc->Eval(origPT)/eventNormWeight;
         }
         else eventWeight = 1.;
 
         histContainer.distrOrigPT->Fill(origPT, eventWeight);
 
         const double bbcz = simCNT.bbcz();
         if (fabs(bbcz) > 30.) continue;

         for(int i = 0; i < simCNT.nch(); i++)
         {
            const double the0 = simCNT.the0(i);
            const double pT = (simCNT.mom(i))*sin(the0);

            if (pT < pTMin || pT > pTMax) continue;
            if (IsQualityCut(simCNT.qual(i))) continue;

            const int charge = simCNT.charge(i);
            if (charge != -1 && charge != 1) continue;

            const int dcarm = simCNT.dcarm(i);

            const double zed = simCNT.zed(i);
            if (fabs(zed) > 75. && fabs(zed) < 3.) continue;

            if (!(fabs(the0) < 100. &&
               ((bbcz > 0. && ((bbcz - 250.*tan(the0 - 3.1416/2.)) > 2. ||
               (bbcz - 200.*tan(the0 - 3.1416/2)) < -2.)) ||
               (bbcz < 0. && ((bbcz - 250.*tan(the0 - 3.1416/2.))< -2. ||
               (bbcz - 200.*tan(the0 - 3.1416/2)) > 2.))))) continue;
 
            //end of basic cuts

            const double alpha = simCNT.alpha(i);
            const double phi = simCNT.phi(i);
            double board;

            if (phi > M_PI/2.) board = ((3.72402 - phi + 0.008047*cos(phi + 0.87851))/0.01963496);
            else board = ((0.573231 + phi - 0.0046 * cos(phi + 0.05721))/0.01963496);

            double alphaReweight = 1.;

            if (dcarm == 0) // DCe
            {
               if (zed >= 0) 
               {
                  histContainer.heatmapUnscaledDCe0->Fill(board, alpha, eventWeight);
                  if (reweightHeatmapsForAlpha) alphaReweight = 
                     alphaReweightDCe0->GetBinContent(alphaReweightDCe0->FindBin(alpha));
                  histContainer.heatmapDCe0->Fill(board, alpha, eventWeight*alphaReweight);
                  histContainer.heatmapDCe0X1->Fill(board, alpha, eventWeight*alphaReweight*
                                                    static_cast<double>(simCNT.nx1hits(i)));
                  histContainer.heatmapDCe0X2->Fill(board, alpha, eventWeight*alphaReweight*
                                                    static_cast<double>(simCNT.nx2hits(i)));
               }
               else 
               {
                  histContainer.heatmapUnscaledDCe1->Fill(board, alpha, eventWeight);
                  if (reweightHeatmapsForAlpha) alphaReweight = 
                     alphaReweightDCe1->GetBinContent(alphaReweightDCe1->FindBin(alpha));
                  histContainer.heatmapDCe1->Fill(board, alpha, eventWeight*alphaReweight);
                  histContainer.heatmapDCe1X1->Fill(board, alpha, eventWeight*alphaReweight*
                                                    static_cast<double>(simCNT.nx1hits(i)));
                  histContainer.heatmapDCe1X2->Fill(board, alpha, eventWeight*alphaReweight*
                                                    static_cast<double>(simCNT.nx2hits(i)));
               }
            } // DCw
            else
            {
               if (zed >= 0) 
               {
                  histContainer.heatmapUnscaledDCw0->Fill(board, alpha, eventWeight);
                  if (reweightHeatmapsForAlpha) alphaReweight = 
                     alphaReweightDCw0->GetBinContent(alphaReweightDCw0->FindBin(alpha));
                  histContainer.heatmapDCw0->Fill(board, alpha, eventWeight*alphaReweight);
                  histContainer.heatmapDCw0X1->Fill(board, alpha, static_cast<double>
                                                    (simCNT.nx1hits(i))*eventWeight*alphaReweight);
                  histContainer.heatmapDCw0X2->Fill(board, alpha, eventWeight*alphaReweight*
                                                    static_cast<double>(simCNT.nx2hits(i)));
               }
               else 
               {
                  histContainer.heatmapUnscaledDCw1->Fill(board, alpha, eventWeight);
                  if (reweightHeatmapsForAlpha) alphaReweight = 
                     alphaReweightDCw1->GetBinContent(alphaReweightDCw1->FindBin(alpha));
                  histContainer.heatmapDCw1->Fill(board, alpha, eventWeight*alphaReweight);
                  histContainer.heatmapDCw1X1->Fill(board, alpha, eventWeight*alphaReweight*
                                                    static_cast<double>(simCNT.nx1hits(i)));
                  histContainer.heatmapDCw1X2->Fill(board, alpha, eventWeight*alphaReweight*
                                                    static_cast<double>(simCNT.nx2hits(i)));
               }
            }

            const bool isParticleOrig = (simCNT.particle_id(i) == particleGeantId && 
                                         simCNT.primary_id(i) == -999);

            if (dmCutter.IsDeadDC(dcarm, zed, board, alpha)) continue;

            double ppc1phi = atan2(simCNT.ppc1y(i), simCNT.ppc1x(i));

            if (dcarm == 1) // PC1w
            {
               histContainer.heatmapPC1w->Fill(simCNT.ppc1z(i), ppc1phi, 
                                               eventWeight*alphaReweight);
            }
            else // PC1e
            {
               if (ppc1phi < 0) ppc1phi += 2.*M_PI;
               histContainer.heatmapPC1e->Fill(simCNT.ppc1z(i), ppc1phi, 
                                               eventWeight*alphaReweight);
            }
            if (dmCutter.IsDeadPC1(dcarm, simCNT.ppc1z(i), ppc1phi)) continue;

            histContainer.distrOrigPTVsRecPT->Fill(origPT, pT, eventWeight);

            if (IsHit(simCNT.pc2dphi(i)))
            {
               const double sdphi = simCalibrator.PC2SDPhi(simCNT.pc2dphi(i), pT, charge);
               const double sdz = simCalibrator.PC2SDZ(simCNT.pc2dz(i), pT, charge);

               if (charge == 1) 
               {
                  histContainer.distrDPhiVsPTPC2Pos->Fill(simCNT.pc2dphi(i), pT, eventWeight);
                  histContainer.distrDZVsPTPC2Pos->Fill(simCNT.pc2dz(i), pT, eventWeight);

                  histContainer.distrSDPhiVsPTPC2Pos->Fill(sdphi, pT, eventWeight);
                  histContainer.distrSDZVsPTPC2Pos->Fill(sdz, pT, eventWeight);
               }
               else
               {
                  histContainer.distrDPhiVsPTPC2Neg->Fill(simCNT.pc2dphi(i), pT, eventWeight);
                  histContainer.distrDZVsPTPC2Neg->Fill(simCNT.pc2dz(i), pT, eventWeight);

                  histContainer.distrSDPhiVsPTPC2Neg->Fill(sdphi, pT, eventWeight);
                  histContainer.distrSDZVsPTPC2Neg->Fill(sdz, pT, eventWeight);
               }

               if (IsMatch(pT, sdphi, sdz, 0.25))
               {
                  const double pc2phi = atan2(simCNT.ppc2y(i), simCNT.ppc2x(i));
                  histContainer.heatmapPC2->Fill(simCNT.ppc2z(i), pc2phi, 
                                                 eventWeight*alphaReweight);
               }
            }

            bool isMatchAndGoodPC3 = false;

            if (IsHit(simCNT.pc3dphi(i)))
            {
               const double sdphi = simCalibrator.PC3SDPhi(simCNT.pc3dphi(i), pT, charge, dcarm);
               const double sdz = simCalibrator.PC3SDZ(simCNT.pc3dz(i), pT, charge, dcarm);

               if (dcarm == 0) // PC3e
               {
                  if (charge == 1) 
                  {
                     histContainer.distrDPhiVsPTPC3ePos->Fill(simCNT.pc3dphi(i), pT, eventWeight);
                     histContainer.distrDZVsPTPC3ePos->Fill(simCNT.pc3dz(i), pT, eventWeight);

                     histContainer.distrSDPhiVsPTPC3ePos->Fill(sdphi, pT, eventWeight);
                     histContainer.distrSDZVsPTPC3ePos->Fill(sdz, pT, eventWeight);
                  }
                  else
                  {
                     histContainer.distrDPhiVsPTPC3eNeg->Fill(simCNT.pc3dphi(i), pT, eventWeight);
                     histContainer.distrDZVsPTPC3eNeg->Fill(simCNT.pc3dz(i), pT, eventWeight);

                     histContainer.distrSDPhiVsPTPC3eNeg->Fill(sdphi, pT, eventWeight);
                     histContainer.distrSDZVsPTPC3eNeg->Fill(sdz, pT, eventWeight);
                  }
               }
               else // PC3w
               {
                  if (charge == 1) 
                  {
                     histContainer.distrDPhiVsPTPC3wPos->Fill(simCNT.pc3dphi(i), pT, eventWeight);
                     histContainer.distrDZVsPTPC3wPos->Fill(simCNT.pc3dz(i), pT, eventWeight);

                     histContainer.distrSDPhiVsPTPC3wPos->Fill(sdphi, pT, eventWeight);
                     histContainer.distrSDZVsPTPC3wPos->Fill(sdz, pT, eventWeight);
                  }
                  else
                  {
                     histContainer.distrDPhiVsPTPC3wNeg->Fill(simCNT.pc3dphi(i), pT, eventWeight);
                     histContainer.distrDZVsPTPC3wNeg->Fill(simCNT.pc3dz(i), pT, eventWeight);

                     histContainer.distrSDPhiVsPTPC3wNeg->Fill(sdphi, pT, eventWeight);
                     histContainer.distrSDZVsPTPC3wNeg->Fill(sdz, pT, eventWeight);
                  }
               }

               if (IsMatch(pT, sdphi, sdz))
               {
                  double pc3phi = atan2(simCNT.ppc3y(i), simCNT.ppc3x(i));

                  if (dcarm == 0 && pc3phi < 0) pc3phi += 2.*M_PI;

                  if (IsMatch(pT, sdphi, sdz, 0.25))
                  {
                     if (dcarm == 0) // PC3e
                     {
                        histContainer.heatmapPC3e->Fill(simCNT.ppc3z(i), pc3phi, 
                                                        eventWeight*alphaReweight);
                     }
                     else // PC3w
                     {
                        histContainer.heatmapPC3w->Fill(simCNT.ppc3z(i), pc3phi, 
                                                        eventWeight*alphaReweight);
                     }
                  }
                  if (!dmCutter.IsDeadPC3(dcarm, simCNT.ppc3z(i), pc3phi)) isMatchAndGoodPC3 = true;
               }
            }

            bool isMatchAndGoodEMCal = false;

            if (IsHit(simCNT.emcdz(i)))
            {
               const double sdphi = 
                  simCalibrator.EMCalSDPhi(simCNT.emcdphi(i), pT, charge, dcarm, simCNT.sect(i));
               const double sdz = 
                  simCalibrator.EMCalSDZ(simCNT.emcdz(i), pT, charge, dcarm, simCNT.sect(i));

               if (dcarm == 0) // EMCale
               {
                  histContainer.distrECoreVsPTEMCale[simCNT.sect(i)]->Fill(pT, simCNT.ecore(i));

                  if (isParticleOrig)
                  {
                     histContainer.distrECoreVsPTEMCaleOrig[simCNT.sect(i)]->
                        Fill(pT, simCNT.ecore(i));
                  }
                  if (charge == 1) 
                  {
                     histContainer.distrProbVsPTEMCale[simCNT.sect(i)]->Fill(pT, simCNT.prob(i));

                     histContainer.distrDPhiVsPTEMCalePos[simCNT.sect(i)]->
                        Fill(simCNT.emcdphi(i), pT, eventWeight);
                     histContainer.distrDZVsPTEMCalePos[simCNT.sect(i)]->
                        Fill(simCNT.emcdz(i), pT, eventWeight);

                     histContainer.distrSDPhiVsPTEMCalePos[simCNT.sect(i)]->
                        Fill(sdphi, pT, eventWeight);
                     histContainer.distrSDZVsPTEMCalePos[simCNT.sect(i)]->
                        Fill(sdz, pT, eventWeight);
                  }
                  else
                  {
                     histContainer.distrDPhiVsPTEMCaleNeg[simCNT.sect(i)]->
                        Fill(simCNT.emcdphi(i), pT, eventWeight);
                     histContainer.distrDZVsPTEMCaleNeg[simCNT.sect(i)]->
                        Fill(simCNT.emcdz(i), pT, eventWeight);

                     histContainer.distrSDPhiVsPTEMCaleNeg[simCNT.sect(i)]->
                        Fill(sdphi, pT, eventWeight);
                     histContainer.distrSDZVsPTEMCaleNeg[simCNT.sect(i)]->
                        Fill(sdz, pT, eventWeight);
                  }
               }
               else // EMCalw
               {
                  histContainer.distrProbVsPTEMCalw[simCNT.sect(i)]->Fill(pT, simCNT.prob(i));
                  histContainer.distrECoreVsPTEMCalw[simCNT.sect(i)]->Fill(pT, simCNT.ecore(i));

                  if (isParticleOrig)
                  {
                     histContainer.distrECoreVsPTEMCalwOrig[simCNT.sect(i)]->
                        Fill(pT, simCNT.ecore(i));
                  }

                  if (charge == 1) 
                  {
                     histContainer.distrDPhiVsPTEMCalwPos[simCNT.sect(i)]->
                        Fill(simCNT.emcdphi(i), pT, eventWeight);
                     histContainer.distrDZVsPTEMCalwPos[simCNT.sect(i)]->
                        Fill(simCNT.emcdz(i), pT, eventWeight);

                     histContainer.distrSDPhiVsPTEMCalwPos[simCNT.sect(i)]->
                        Fill(sdphi, pT, eventWeight);
                     histContainer.distrSDZVsPTEMCalwPos[simCNT.sect(i)]->
                        Fill(sdz, pT, eventWeight);
                  }
                  else
                  {
                     histContainer.distrDPhiVsPTEMCalwNeg[simCNT.sect(i)]->
                        Fill(simCNT.emcdphi(i), pT, eventWeight);
                     histContainer.distrDZVsPTEMCalwNeg[simCNT.sect(i)]->
                        Fill(simCNT.emcdz(i), pT, eventWeight);

                     histContainer.distrSDPhiVsPTEMCalwNeg[simCNT.sect(i)]->
                        Fill(sdphi, pT, eventWeight);
                     histContainer.distrSDZVsPTEMCalwNeg[simCNT.sect(i)]->
                        Fill(sdz, pT, eventWeight);
                  }
               }

               if (IsMatch(pT, sdphi, sdz))
               {
                  bool isCutByECore;
                  if (dcarm == 0 && simCNT.sect(i) < 2) isCutByECore = (simCNT.ecore(i) < 0.35);
                  else isCutByECore = (simCNT.ecore(i) < 0.25); // PbSc

                  if (!isCutByECore && IsMatch(pT, sdphi, sdz, 0.25))
                  { 
                     if (dcarm == 0) // EMCale
                     {
                        histContainer.heatmapEMCale[simCNT.sect(i)]->
                           Fill(static_cast<double>(simCNT.ysect(i)) + 0.5, 
                                static_cast<double>(simCNT.zsect(i)) + 0.5, 
                                simCNT.ecore(i)*eventWeight*alphaReweight);
                        histContainer.heatmapEMCaleHit[simCNT.sect(i)]->
                           Fill(static_cast<double>(simCNT.ysect(i)) + 0.5, 
                                static_cast<double>(simCNT.zsect(i)) + 0.5, 
                                eventWeight*alphaReweight);
                     }
                     else // EMCalw
                     {
                        histContainer.heatmapEMCalw[simCNT.sect(i)]->
                           Fill(static_cast<double>(simCNT.ysect(i)) + 0.5, 
                                static_cast<double>(simCNT.zsect(i)) + 0.5, 
                                simCNT.ecore(i)*eventWeight*alphaReweight);
                        histContainer.heatmapEMCalwHit[simCNT.sect(i)]->
                           Fill(static_cast<double>(simCNT.ysect(i)) + 0.5, 
                                static_cast<double>(simCNT.zsect(i)) + 0.5, 
                                eventWeight*alphaReweight);
                     }
                  }

                  if (!dmCutter.IsDeadEMCal(dcarm, simCNT.sect(i), 
                                            simCNT.ysect(i), simCNT.zsect(i)))
                  {
                     isMatchAndGoodEMCal = true;

                     if (!isCutByECore)
                     {
                        const double tExpPi = sqrt(simCNT.plemc(i)*simCNT.plemc(i)/
                                                   (SPEED_OF_LIGHT*SPEED_OF_LIGHT)*
                                                   (MASS_PION*MASS_PION/
                                                    simCNT.mom(i)*simCNT.mom(i) + 1.));

                        const double m2 = simCNT.mom(i)*simCNT.mom(i)*
                                          (simCNT.temc(i)*simCNT.temc(i)*
                                           SPEED_OF_LIGHT*SPEED_OF_LIGHT/
                                           (simCNT.plemc(i)*simCNT.plemc(i) - 1.));

                        if (dcarm == 0)
                        {
                           histContainer.distrTEMCale[simCNT.sect(i)]->
                              Fill(simCNT.temc(i) - tExpPi, pT, eventWeight);
                           if (charge == 1)
                           {
                              histContainer.distrM2EMCalePosCharge[simCNT.sect(i)]->
                                 Fill(pT, m2, eventWeight);
                           }
                           else
                           {
                              histContainer.distrM2EMCaleNegCharge[simCNT.sect(i)]->
                                 Fill(pT, m2, eventWeight);
                           }
                        }
                        else
                        {
                           histContainer.distrTEMCalw[simCNT.sect(i)]->
                              Fill(simCNT.temc(i) - tExpPi, pT, eventWeight);
                           if (charge == 1)
                           {
                              histContainer.distrM2EMCalwPosCharge[simCNT.sect(i)]->
                                 Fill(pT, m2, eventWeight);
                           }
                           else
                           {
                              histContainer.distrM2EMCalwNegCharge[simCNT.sect(i)]->
                                 Fill(pT, m2, eventWeight);
                           }
                        }
                     }
                  }
               }
            }

            if (IsHit(simCNT.tofdz(i)))
            {
               if (charge == 1) 
               {
                  histContainer.distrDPhiVsPTTOFePos->Fill(simCNT.tofdphi(i), pT, eventWeight);
                  histContainer.distrDZVsPTTOFePos->Fill(simCNT.tofdz(i), pT, eventWeight);

                  histContainer.distrSDPhiVsPTTOFePos->
                     Fill(simCalibrator.TOFeSDPhi(simCNT.tofdphi(i), pT, charge),
                          pT, eventWeight);
                  histContainer.distrSDZVsPTTOFePos->
                     Fill(simCalibrator.TOFeSDZ(simCNT.tofdz(i), pT, charge),
                          pT, eventWeight);
               }
               else
               {
                  histContainer.distrDPhiVsPTTOFeNeg->Fill(simCNT.tofdphi(i), pT, eventWeight);
                  histContainer.distrDZVsPTTOFeNeg->Fill(simCNT.tofdz(i), pT, eventWeight);

                  histContainer.distrSDPhiVsPTTOFeNeg->
                     Fill(simCalibrator.TOFeSDPhi(simCNT.tofdphi(i), pT, charge), pT, eventWeight);
                  histContainer.distrSDZVsPTTOFeNeg->
                     Fill(simCalibrator.TOFeSDZ(simCNT.tofdz(i), pT, charge), pT, eventWeight);
               }

               const double sdphi = simCalibrator.TOFeSDPhi(simCNT.tofdphi(i), pT, charge);
               const double sdz = simCalibrator.TOFeSDZ(simCNT.tofdz(i), pT, charge);

               const double beta = simCNT.pltof(i)/simCNT.ttof(i)/29.9792;
               const double eloss = 0.0005*pow(beta, -2.5);
               histContainer.distrBetaVsETOFe->Fill(beta, simCNT.etof(i));

               if (simCNT.etof(i) > eloss && IsMatch(pT, sdphi, sdz))
               {
                  // slats are organized in 10 lines of 96 we define as chambers
                  const int chamber = simCNT.slat(i)/96;
                  // slat number for the current chamber
                  const int slat = simCNT.slat(i) % 96;

                  if (IsMatch(pT, sdphi, sdz, 0.25))
                  {
                     histContainer.heatmapTOFe->Fill(static_cast<double>(chamber) + 0.5, 
                                                     static_cast<double>(slat) + 0.5, 
                                                     simCNT.etof(i)*eventWeight*alphaReweight);
                  }

                  if (isMatchAndGoodPC3 && isMatchAndGoodEMCal &&
                      !dmCutter.IsDeadTOFe(chamber, slat))
                  {
                     const double tExpPi = sqrt(simCNT.pltof(i)*simCNT.pltof(i)/
                                                (SPEED_OF_LIGHT*SPEED_OF_LIGHT)*
                                                (MASS_PION*MASS_PION/
                                                 simCNT.mom(i)*simCNT.mom(i) + 1.));

                     const double m2 = simCNT.mom(i)*simCNT.mom(i)*
                                       (simCNT.ttof(i)*simCNT.ttof(i)*
                                        SPEED_OF_LIGHT*SPEED_OF_LIGHT/
                                        (simCNT.pltof(i)*simCNT.pltof(i) - 1.));

                     histContainer.distrTTOFe->Fill(simCNT.ttof(i) - tExpPi, pT, eventWeight);

                     if (charge == 1) histContainer.distrM2TOFePosCharge->Fill(pT, m2, eventWeight);
                     else histContainer.distrM2TOFeNegCharge->Fill(pT, m2, eventWeight);
                  }
               }
            }
            else if (IsHit(simCNT.tofwdz(i)))
            {
               if (charge == 1) 
               {
                  histContainer.distrDPhiVsPTTOFwPos->Fill(simCNT.tofwdphi(i), pT, eventWeight);
                  histContainer.distrDZVsPTTOFwPos->Fill(simCNT.tofwdz(i), pT, eventWeight);

                  histContainer.distrSDPhiVsPTTOFwPos->
                     Fill(simCalibrator.TOFwSDPhi(simCNT.tofwdphi(i), pT, charge), 
                          pT, eventWeight);
                  histContainer.distrSDZVsPTTOFwPos->
                     Fill(simCalibrator.TOFwSDZ(simCNT.tofwdz(i), pT, charge), 
                          pT, eventWeight);
               }
               else
               {
                  histContainer.distrDPhiVsPTTOFwNeg->Fill(simCNT.tofwdphi(i), pT, eventWeight);
                  histContainer.distrDZVsPTTOFwNeg->Fill(simCNT.tofwdz(i), pT, eventWeight);

                  histContainer.distrSDPhiVsPTTOFwNeg->
                     Fill(simCalibrator.TOFwSDPhi(simCNT.tofwdphi(i), pT, charge),
                          pT, eventWeight);
                  histContainer.distrSDZVsPTTOFwNeg->
                     Fill(simCalibrator.TOFwSDZ(simCNT.tofwdz(i), pT, charge),
                          pT, eventWeight);
               }

               const double sdphi = simCalibrator.TOFwSDPhi(simCNT.tofwdphi(i), pT, charge);
               const double sdz = simCalibrator.TOFwSDZ(simCNT.tofwdz(i), pT, charge);

               // strips are organized in 8 lines of 64 we define as chambers
               const int chamber = simCNT.striptofw(i)/64;
               // strip number for the current chamber
               const int strip = simCNT.striptofw(i) % 64;

               if (IsMatch(pT, sdphi, sdz))
               {
                  if (IsMatch(pT, sdphi, sdz, 0.25))
                  {
                     histContainer.heatmapTOFw->Fill(static_cast<double>(chamber) + 0.5, 
                                                     static_cast<double>(strip) + 0.5, 
                                                     eventWeight*alphaReweight);
                  }

                  if (isMatchAndGoodPC3 && isMatchAndGoodEMCal &&
                      !dmCutter.IsDeadTOFw(chamber, strip))
                  {
                     const double tExpPi = sqrt(simCNT.pltofw(i)*simCNT.pltofw(i)/
                                                (SPEED_OF_LIGHT*SPEED_OF_LIGHT)*
                                                (MASS_PION*MASS_PION/
                                                 simCNT.mom(i)*simCNT.mom(i) + 1.));

                     const double m2 = simCNT.mom(i)*simCNT.mom(i)*
                                       (simCNT.ttofw(i)*simCNT.ttofw(i)*
                                        SPEED_OF_LIGHT*SPEED_OF_LIGHT/
                                        (simCNT.pltofw(i)*simCNT.pltofw(i) - 1.));

                     histContainer.distrTTOFw->Fill(simCNT.ttofw(i) - tExpPi, pT, eventWeight);

                     if (charge == 1) histContainer.distrM2TOFwPosCharge->Fill(pT, m2, eventWeight);
                     else histContainer.distrM2TOFwNegCharge->Fill(pT, m2, eventWeight);
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

         TH2F *realDataDCe0 = GetHistogramFromFile(realDataInputFile, "_Heatmap: DCe, zDC>=0");
         TH2F *realDataDCe1 = GetHistogramFromFile(realDataInputFile, "_Heatmap: DCe, zDC<0");
         TH2F *realDataDCw0 = GetHistogramFromFile(realDataInputFile, "_Heatmap: DCw, zDC>=0");
         TH2F *realDataDCw1 = GetHistogramFromFile(realDataInputFile, "_Heatmap: DCw, zDC<0");

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
         CppTools::PrintInfo("Alpha reweight is now disabled. The missing file will be created "\
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
      pBar.Finish();
      CppTools::PrintInfo("AnalyzerSingleTrack has finished processing simulated data");
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
                                 particle["geant_id"].as<short>(), 
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
         thrContainer.Write(outputFileName);
      }
   }

   isProcessFinished = true;
   pBarThread.join();

   CppTools::PrintInfo("Merging output files into one");

   std::string haddCommand = "hadd -f6 " + outputDir + "all.root ";
   for (const auto& particle : inputYAMLSim["particles"])
   {
      haddCommand += outputDir + particle["name"].as<std::string>() + ".root ";
   }
   system(haddCommand.c_str());

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

   copy.distrOrigPT = distrOrigPT->Get();
   copy.distrOrigPTVsRecPT = distrOrigPTVsRecPT.Get();
   copy.heatmapUnscaledDCe0 = heatmapUnscaledDCe0.Get();
   copy.heatmapUnscaledDCe1 = heatmapUnscaledDCe1.Get();
   copy.heatmapUnscaledDCw0 = heatmapUnscaledDCw0.Get();
   copy.heatmapUnscaledDCw1 = heatmapUnscaledDCw1.Get();
   copy.heatmapDCe0 = heatmapDCe0.Get();
   copy.heatmapDCe1 = heatmapDCe1.Get();
   copy.heatmapDCw0 = heatmapDCw0.Get();
   copy.heatmapDCw1 = heatmapDCw1.Get();
   copy.heatmapDCe0X1 = heatmapDCe0X1.Get();
   copy.heatmapDCe1X1 = heatmapDCe1X1.Get();
   copy.heatmapDCw0X1 = heatmapDCw0X1.Get();
   copy.heatmapDCw1X1 = heatmapDCw1X1.Get();
   copy.heatmapDCe0X2 = heatmapDCe0X2.Get();
   copy.heatmapDCe1X2 = heatmapDCe1X2.Get();
   copy.heatmapDCw0X2 = heatmapDCw0X2.Get();
   copy.heatmapDCw1X2 = heatmapDCw1X2.Get();
   copy.heatmapPC1e = heatmapPC1e.Get();
   copy.heatmapPC1w = heatmapPC1w.Get();
   copy.heatmapPC2 = heatmapPC2.Get();
   copy.heatmapPC3e = heatmapPC3e.Get();
   copy.heatmapPC3w = heatmapPC3w.Get();
   copy.heatmapTOFe = heatmapTOFe.Get();
   copy.heatmapTOFw = heatmapTOFw.Get();
   copy.distrBetaVsETOFe = distrBetaVsETOFe.Get();
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
   copy.distrTTOFe = distrTTOFe.Get();
   copy.distrTTOFw = distrTTOFw.Get();
   copy.distrM2TOFePosCharge = distrM2TOFePosCharge.Get();
   copy.distrM2TOFeNegCharge = distrM2TOFeNegCharge.Get();
   copy.distrM2TOFwPosCharge = distrM2TOFwPosCharge.Get();
   copy.distrM2TOFwNegCharge = distrM2TOFwNegCharge.Get();

   // iterating over EMCal sectors
   for (int i = 0; i < 4; i++)
   {
      copy.heatmapEMCale[i] = heatmapEMCale[i].Get();
      copy.heatmapEMCalw[i] = heatmapEMCalw[i].Get();
      copy.heatmapEMCaleHit[i] = heatmapEMCaleHit[i].Get();
      copy.heatmapEMCalwHit[i] = heatmapEMCalwHit[i].Get();
      copy.distrECoreVsPTEMCale[i] = distrECoreVsPTEMCale[i].Get();
      copy.distrECoreVsPTEMCalw[i] = distrECoreVsPTEMCalw[i].Get();
      copy.distrECoreVsPTEMCaleOrig[i] = distrECoreVsPTEMCaleOrig[i].Get();
      copy.distrECoreVsPTEMCalwOrig[i] = distrECoreVsPTEMCalwOrig[i].Get();
      copy.distrProbVsPTEMCale[i] = distrProbVsPTEMCale[i].Get();
      copy.distrProbVsPTEMCalw[i] = distrProbVsPTEMCalw[i].Get();
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
      copy.distrTEMCale[i] = distrTEMCale[i].Get();
      copy.distrTEMCalw[i] = distrTEMCalw[i].Get();
      copy.distrM2EMCalePosCharge[i] = distrM2EMCalePosCharge[i].Get();
      copy.distrM2EMCaleNegCharge[i] = distrM2EMCaleNegCharge[i].Get();
      copy.distrM2EMCalwPosCharge[i] = distrM2EMCalwPosCharge[i].Get();
      copy.distrM2EMCalwNegCharge[i] = distrM2EMCalwNegCharge[i].Get();
   }
   return copy;
}

void AnalyzeSingleTrack::ThrContainer::Write(const std::string& outputFileName)
{
   TFile outputFile(outputFileName.c_str(), "RECREATE");
   outputFile.SetCompressionLevel(6);
   outputFile.cd();

   static_cast<std::shared_ptr<TH1F>>(distrOrigPT->Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrOrigPTVsRecPT.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(heatmapUnscaledDCe0.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(heatmapUnscaledDCe1.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(heatmapUnscaledDCw0.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(heatmapUnscaledDCw1.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(heatmapDCe0.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(heatmapDCe1.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(heatmapDCw0.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(heatmapDCw1.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(heatmapDCe0X1.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(heatmapDCe1X1.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(heatmapDCw0X1.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(heatmapDCw1X1.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(heatmapDCe0X2.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(heatmapDCe1X2.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(heatmapDCw0X2.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(heatmapDCw1X2.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(heatmapPC1e.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(heatmapPC1w.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(heatmapPC2.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(heatmapPC3e.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(heatmapPC3w.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(heatmapTOFe.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(heatmapTOFw.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrBetaVsETOFe.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrDZVsPTPC2Pos.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrTTOFe.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrTTOFw.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrM2TOFePosCharge.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrM2TOFeNegCharge.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrM2TOFwPosCharge.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrM2TOFwNegCharge.Merge())->Write();

   // iterating over EMCal sectors
   for (int i = 0; i < 4; i++)
   {
      static_cast<std::shared_ptr<TH2F>>(heatmapEMCale[i].Merge())->Write();
      static_cast<std::shared_ptr<TH2F>>(heatmapEMCalw[i].Merge())->Write();
      static_cast<std::shared_ptr<TH2F>>(heatmapEMCaleHit[i].Merge())->Write();
      static_cast<std::shared_ptr<TH2F>>(heatmapEMCalwHit[i].Merge())->Write();
      static_cast<std::shared_ptr<TH2F>>(distrECoreVsPTEMCale[i].Merge())->Write();
      static_cast<std::shared_ptr<TH2F>>(distrECoreVsPTEMCalw[i].Merge())->Write();
      static_cast<std::shared_ptr<TH2F>>(distrECoreVsPTEMCaleOrig[i].Merge())->Write();
      static_cast<std::shared_ptr<TH2F>>(distrECoreVsPTEMCalwOrig[i].Merge())->Write();
      static_cast<std::shared_ptr<TH2F>>(distrProbVsPTEMCale[i].Merge())->Write();
      static_cast<std::shared_ptr<TH2F>>(distrProbVsPTEMCalw[i].Merge())->Write();
      static_cast<std::shared_ptr<TH2F>>(distrTEMCale[i].Merge())->Write();
      static_cast<std::shared_ptr<TH2F>>(distrTEMCalw[i].Merge())->Write();
      static_cast<std::shared_ptr<TH2F>>(distrM2EMCalePosCharge[i].Merge())->Write();
      static_cast<std::shared_ptr<TH2F>>(distrM2EMCaleNegCharge[i].Merge())->Write();
      static_cast<std::shared_ptr<TH2F>>(distrM2EMCalwPosCharge[i].Merge())->Write();
      static_cast<std::shared_ptr<TH2F>>(distrM2EMCalwNegCharge[i].Merge())->Write();
   }

   outputFile.mkdir("dphi");
   outputFile.cd("dphi");

   static_cast<std::shared_ptr<TH2F>>(distrDPhiVsPTPC2Pos.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrDPhiVsPTPC2Neg.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrDPhiVsPTPC3ePos.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrDPhiVsPTPC3eNeg.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrDPhiVsPTPC3wPos.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrDPhiVsPTPC3wNeg.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrDPhiVsPTTOFePos.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrDPhiVsPTTOFeNeg.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrDPhiVsPTTOFwPos.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrDPhiVsPTTOFwNeg.Merge())->Write();

   for (int i = 0; i < 4; i++)
   {
      static_cast<std::shared_ptr<TH2F>>(distrDPhiVsPTEMCalePos[i].Merge())->Write();
      static_cast<std::shared_ptr<TH2F>>(distrDPhiVsPTEMCaleNeg[i].Merge())->Write();
      static_cast<std::shared_ptr<TH2F>>(distrDPhiVsPTEMCalwPos[i].Merge())->Write();
      static_cast<std::shared_ptr<TH2F>>(distrDPhiVsPTEMCalwNeg[i].Merge())->Write();
   }

   outputFile.mkdir("dz");
   outputFile.cd("dz");

   static_cast<std::shared_ptr<TH2F>>(distrDZVsPTPC2Neg.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrDZVsPTPC3ePos.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrDZVsPTPC3eNeg.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrDZVsPTPC3wPos.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrDZVsPTPC3wNeg.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrDZVsPTTOFePos.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrDZVsPTTOFeNeg.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrDZVsPTTOFwPos.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrDZVsPTTOFwNeg.Merge())->Write();

   for (int i = 0; i < 4; i++)
   {
      static_cast<std::shared_ptr<TH2F>>(distrDZVsPTEMCalePos[i].Merge())->Write();
      static_cast<std::shared_ptr<TH2F>>(distrDZVsPTEMCaleNeg[i].Merge())->Write();
      static_cast<std::shared_ptr<TH2F>>(distrDZVsPTEMCalwPos[i].Merge())->Write();
      static_cast<std::shared_ptr<TH2F>>(distrDZVsPTEMCalwNeg[i].Merge())->Write();
   }

   outputFile.mkdir("sdphi");
   outputFile.cd("sdphi");

   static_cast<std::shared_ptr<TH2F>>(distrSDPhiVsPTPC2Pos.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrSDPhiVsPTPC2Neg.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrSDPhiVsPTPC3ePos.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrSDPhiVsPTPC3eNeg.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrSDPhiVsPTPC3wPos.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrSDPhiVsPTPC3wNeg.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrSDPhiVsPTTOFePos.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrSDPhiVsPTTOFeNeg.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrSDPhiVsPTTOFwPos.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrSDPhiVsPTTOFwNeg.Merge())->Write();

   for (int i = 0; i < 4; i++)
   {
      static_cast<std::shared_ptr<TH2F>>(distrSDPhiVsPTEMCalePos[i].Merge())->Write();
      static_cast<std::shared_ptr<TH2F>>(distrSDPhiVsPTEMCaleNeg[i].Merge())->Write();
      static_cast<std::shared_ptr<TH2F>>(distrSDPhiVsPTEMCalwPos[i].Merge())->Write();
      static_cast<std::shared_ptr<TH2F>>(distrSDPhiVsPTEMCalwNeg[i].Merge())->Write();
   }

   outputFile.mkdir("sdz");
   outputFile.cd("sdz");

   static_cast<std::shared_ptr<TH2F>>(distrSDZVsPTPC2Pos.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrSDZVsPTPC2Neg.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrSDZVsPTPC3ePos.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrSDZVsPTPC3eNeg.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrSDZVsPTPC3wPos.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrSDZVsPTPC3wNeg.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrSDZVsPTTOFePos.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrSDZVsPTTOFeNeg.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrSDZVsPTTOFwPos.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrSDZVsPTTOFwNeg.Merge())->Write();

   for (int i = 0; i < 4; i++)
   {
      static_cast<std::shared_ptr<TH2F>>(distrSDZVsPTEMCalePos[i].Merge())->Write();
      static_cast<std::shared_ptr<TH2F>>(distrSDZVsPTEMCaleNeg[i].Merge())->Write();
      static_cast<std::shared_ptr<TH2F>>(distrSDZVsPTEMCalwPos[i].Merge())->Write();
      static_cast<std::shared_ptr<TH2F>>(distrSDZVsPTEMCalwNeg[i].Merge())->Write();
   }

   outputFile.Close();
}

#endif /* ANALYZE_SINGLE_TRACK_CPP */
