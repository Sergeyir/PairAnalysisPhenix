/** 
 *  @file   AnalyzeResonance.hpp
 *  @brief  Contains declarations of functions and variables that are used for analysis of a resonance from a trees acquired from the PHENIX simulation
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysisPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef ANALYZE_RESONANCE_CPP
#define ANALYZE_RESONANCE_CPP

#include "../include/AnalyzeResonance.hpp"

// this namespace is only used so that documentation does not become a mess
// so there is no need to enforce the contents inside of it 
// being accessed only via the scope resolution operator in this file
using namespace AnalyzeResonance;

void AnalyzeResonance::AnalyzeConfiguration(ThrContainer &thrContainer, 
                                                     const std::string& particleName, 
                                                     const int daughter1Id,
                                                     const int daughter2Id,
                                                     const std::string& magneticFieldName, 
                                                     const std::string &pTRangeName)
{ 
   std::string simInputFileName = "data/SimTrees/" + runName + "/Resonance/" + 
                                  particleName + "_" + ParticleMap::name[daughter1Id] + 
                                  ParticleMap::name[daughter2Id] + "_" + 
                                  pTRangeName + magneticFieldName + ".root";

   TFile simInputFile = TFile(simInputFileName.c_str());

   TH1D *origPTHist = static_cast<TH1D *>(simInputFile.Get("orig_pt"));

   // weight function for spectra
   std::unique_ptr<TF1> weightFunc;
  
   // normalization of the number of particles to the number of events
   // this normalization is needed to seamlessly merge 2 files with 
   // flat pT distribution with different ranges
   double eventNormWeight = 1.;

   std::string realDataFileName = "data/Real/" + runName + "/SingleTrack/sum" + 
                                  magneticFieldName + ".root";
   TFile realDataFile = TFile(realDataFileName.c_str());

   // threshold is needed since there can be a little noise in the histogram
   const double origPTThreshold = 
      origPTHist->Integral()/static_cast<double>(origPTHist->GetXaxis()->GetNbins())/2.;

   // pT bounds of original pT distribution for spectra scale
   const double lowPTBound = 
      origPTHist->GetXaxis()->GetBinLowEdge(origPTHist->FindFirstBinAbove(origPTThreshold));
   const double upPTBound = 
      origPTHist->GetXaxis()->GetBinUpEdge(origPTHist->FindLastBinAbove(origPTThreshold));

   TH1D *centrHist = static_cast<TH1D *>(realDataFile.Get("centrality"));
   if (!centrHist) 
   {
      CppTools::PrintError("Histogram \"centrality\" does not exist in file" + 
                           static_cast<std::string>(realDataFile.GetName()));
   }

   const double resonancePTMin = inputYAMLResonance["pt_bins"][0]["min"].as<double>();
   const double resonancePTMax = 
      inputYAMLResonance["pt_bins"][inputYAMLResonance["pt_bins"].size() - 1]["max"].as<double>();

   eventNormWeight = origPTHist->Integral(1, origPTHist->GetXaxis()->GetNbins())/
                     centrHist->Integral(1, centrHist->GetXaxis()->GetNbins());

   eventNormWeight *= (resonancePTMax - resonancePTMin)/(upPTBound - lowPTBound);

   if (reweightForSpectra)
   {
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
   else
   {
      weightFunc = std::make_unique<TF1>("weightFunc", "expo");
      weightFunc->SetParameters(0., -1.);
   }

   const double resonanceMass = inputYAMLResonance["mass"].as<double>();
   const double resonanceGamma = inputYAMLResonance["gamma"].as<double>();

   const double daughter1Mass = ParticleMap::mass[daughter1Id];
   const double daughter2Mass = ParticleMap::mass[daughter2Id];

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
 
         eventWeight = weightFunc->Eval(origPT)/eventNormWeight;
 
         histContainer.distrOrigPT->Fill(origPT, eventWeight);
 
         const double bbcz = simCNT.bbcz();
         if (fabs(bbcz) > 30.) continue;

         std::vector<ChargedTrack> positiveTracks;
         std::vector<ChargedTrack> negativeTracks;

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

            if (dmCutter.IsDeadDC(dcarm, zed, board, alpha)) continue;

            double ppc1phi = atan2(simCNT.ppc1y(i), simCNT.ppc1x(i));
            if (dcarm == 0 && ppc1phi < 0) ppc1phi += 2.*M_PI;
            if (dmCutter.IsDeadPC1(dcarm, simCNT.ppc1z(i), ppc1phi)) continue;

            switch (charge)
            {
               case 1:
                  positiveTracks.emplace_back(daughter1Mass, simCNT, i);
                  break;
               case -1:
                  negativeTracks.emplace_back(daughter2Mass, simCNT, i);
                  break;
            }
         }

         for (const auto& posTrack : positiveTracks)
         {
            for (const auto& negTrack : negativeTracks)
            {
               // invariant mass [GeV/c^2]
               const double mInv = GetPairMass(posTrack, negTrack);
               // pT of a pair [GeV/c]
               const double pT = GetPairPT(posTrack, negTrack);

               // check that shows whether invariant mass is within 2 gamma from mean of the signal
               // 10 is a rough estimation for gaussian widening 
               // due to finite momentum resolution of a detector system 
               const bool isWithin2Gamma = (mInv > resonanceMass - resonanceGamma*2. - 10. && 
                                            mInv < resonanceMass + resonanceGamma*2. + 10.);

               thrContainer.distrMInv->Fill(pT, mInv, eventWeight);

               if (IsGhostCut(posTrack, negTrack))
               {
                  thrContainer.distrMInvNoPIDGhostAntiCut->Fill(pT, mInv, eventWeight);
                  continue;
               }

               if (IsOneArmCut(posTrack, negTrack)) 
               {
                  thrContainer.distrMInvNoPIDOneArmAntiCut->Fill(pT, mInv, eventWeight);
                  continue;
               }

               if (isWithin2Gamma)
               {
                  histContainer.distrPAsymVsPT->Fill(pT, (posTrack.p - negTrack.p)/
                                                     (posTrack.p + negTrack.p), mInv, eventWeight);
                  histContainer.distrDPhiVsPT->Fill(pT, posTrack.phi - negTrack.phi, 
                                                    mInv, eventWeight);
                  histContainer.distrDAlphaVsPT->Fill(pT, posTrack.alpha - negTrack.alpha, 
                                                      mInv, eventWeight);
                  histContainer.distrDZedVsPT->Fill(pT, posTrack.zed - negTrack.zed, 
                                                    mInv, eventWeight);
               }

               if (isWithin2Gamma)
               {
                  if (simCNT.sect(posTrack.index) == simCNT.sect(negTrack.index) &&
                      simCNT.ysect(posTrack.index) == simCNT.ysect(negTrack.index) &&
                      simCNT.zsect(posTrack.index) == simCNT.zsect(negTrack.index))
                  {
                     if (IsHit(simCNT.pc2dz(posTrack.index)) && IsHit(simCNT.pc2dz(negTrack.index)))
                     {
                        thrContainer.distrDPC2PhiDPC2ZVsPT->
                           Fill(simCNT.ppc2z(posTrack.index) - simCNT.ppc2z(negTrack.index),
                                atan2(simCNT.ppc2y(posTrack.index), simCNT.ppc2x(posTrack.index)) - 
                                atan2(simCNT.ppc2y(negTrack.index), simCNT.ppc2x(negTrack.index)),
                                eventWeight);
                     }

                     if (IsHit(simCNT.pc3dz(posTrack.index)) && IsHit(simCNT.pc3dz(negTrack.index)))
                     {
                        double pc3phi1 = atan2(simCNT.ppc3y(posTrack.index), 
                                               simCNT.ppc3x(posTrack.index));
                        double pc3phi2 = atan2(simCNT.ppc3y(negTrack.index), 
                                               simCNT.ppc3x(negTrack.index));

                        if (simCNT.dcarm(posTrack.index) == 0 && pc3phi1 < 0) pc3phi1 += 2.*M_PI;
                        if (simCNT.dcarm(negTrack.index) == 0 && pc3phi2 < 0) pc3phi2 += 2.*M_PI;

                        thrContainer.distrDPC3PhiDPC3ZVsPT->
                           Fill(simCNT.ppc3z(posTrack.index) - simCNT.ppc3z(negTrack.index), 
                                pc3phi1 - pc3phi2, eventWeight);
                     }

                     if (IsHit(simCNT.tofdz(posTrack.index)) && 
                         IsHit(simCNT.tofdz(negTrack.index)))
                     {
                        thrContainer.distrDChamberDSlatVsPT->
                           Fill(static_cast<double>(simCNT.slat(posTrack.index)/96 - 
                                                    simCNT.slat(negTrack.index)/96) + 0.5,
                                static_cast<double>((simCNT.slat(posTrack.index) % 96) - 
                                                    (simCNT.slat(negTrack.index) % 96)) + 0.5,
                                pT, eventWeight);
                     }
                     else if (IsHit(simCNT.tofwdz(posTrack.index)) && 
                              IsHit(simCNT.tofwdz(negTrack.index)))
                     {
                        thrContainer.distrDChamberDStripVsPT->
                           Fill(static_cast<double>(simCNT.striptofw(posTrack.index)/96 - 
                                                    simCNT.striptofw(negTrack.index)/96) + 0.5,
                                static_cast<double>((simCNT.striptofw(posTrack.index) % 96) - 
                                                    (simCNT.striptofw(negTrack.index) % 96)) + 0.5,
                                pT, eventWeight);
                     }
                  }

                  if (IsHit(simCNT.emcdz(posTrack.index)) &&
                      simCNT.sect(posTrack.index) == simCNT.sect(negTrack.index))
                  {
                     thrContainer.distrDYTowerDZTowerVsPT->
                        Fill(static_cast<double>(simCNT.ysect(posTrack.index) - 
                                                 simCNT.ysect(negTrack.index)) + 0.5, 
                             static_cast<double>(simCNT.zsect(posTrack.index) - 
                                                 simCNT.zsect(negTrack.index)) + 0.5, 
                             pT, eventWeight);
                  }
               }

               if (!IsNoPID(posTrack, negTrack)) continue;

               thrContainer.distrMInvNoPID->Fill(pT, mInv, eventWeight);

               if (mInv > resonanceMass - resonanceGamma*2. - 10. && 
                   mInv < resonanceMass + resonanceGamma*2. + 10.)
               {
                  histContainer.distrOrigPTVsRecPT->Fill(origPT, pT, eventWeight);
               }

               if (IsSailorCut(posTrack, negTrack))
               {
                  thrContainer.distrMInvNoPIDSailorCut->Fill(pT, mInv, eventWeight);
               }
               else
               {
                  thrContainer.distrMInvNoPIDCowboyCut->Fill(pT, mInv, eventWeight);
               }

               //if (!Is2PID(posTrack, negTrack, daughter1Id, daughter2Id)) continue;
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
      errMsg += "parameter(s) were provided \n Usage: bin/AnalyzeResonance ";
      errMsg += "inputYAMLName numberOfThreads=std::thread::hardware_concurrency()";
      CppTools::PrintError(errMsg);
   }
 
   CppTools::CheckInputFile(argv[1]);
 
   if (argc == 2) numberOfThreads = std::thread::hardware_concurrency();
   else numberOfThreads = std::stoi(argv[2]);

   ROOT::EnableImplicitMT(numberOfThreads);

   inputYAMLResonance.OpenFile(argv[1]);
   inputYAMLResonance.CheckStatus("resonance");

   runName = inputYAMLResonance["run_name"].as<std::string>();

   inputYAMLMain.OpenFile("input/" + runName + "/main.yaml");
   inputYAMLMain.CheckStatus("main");

   inputYAMLSimSingleTrack.OpenFile("input/" + runName + "/single_track_sim.yaml");
   inputYAMLSimSingleTrack.CheckStatus("single_track_sim");
 
   collisionSystemName = inputYAMLMain["collision_system_name"].as<std::string>();

   outputDir = "data/PostSim/" + runName + "/Resonance/";
   system(("mkdir -p " + outputDir).c_str());

   pTMin = inputYAMLSimSingleTrack["pt_min"].as<double>();
   pTMax = inputYAMLSimSingleTrack["pt_max"].as<double>();

   dmCutter.Initialize(runName, inputYAMLMain["detectors_configuration"].as<std::string>());

   if (CppTools::FileExists("data/Parameters/SpectraFit/" + collisionSystemName + 
                            "/" + inputYAMLResonance["name"].as<std::string>() + ".yaml"))
   {
      CppTools::PrintInfo("Fit parameters for spectra were found");
      reweightForSpectra = true;
   }
   else 
   {
      CppTools::PrintInfo("Fit parameters for spectra were not found; setting reweight to e^{-1*pT}");
      reweightForSpectra = false;
   }
 
   for (const auto& magneticField : inputYAMLMain["magnetic_field_configurations"])
   {
      CppTools::CheckInputFile("data/Real/" + runName + "/SingleTrack/sum" + 
                               magneticField["name"].as<std::string>() + ".root");
      for (const auto& pTRange : inputYAMLResonance["pt_ranges"])
      {
         std::string simInputFileName = 
            "data/SimTrees/" + runName + "/Resonance/" + 
            inputYAMLResonance["name"].as<std::string>() + "_" + 
            ParticleMap::name[inputYAMLResonance["daughter1_id"].as<int>()] + 
            ParticleMap::name[inputYAMLResonance["daughter2_id"].as<int>()] + 
            "_" + pTRange["name"].as<std::string>() + 
            magneticField["name"].as<std::string>() + ".root";

         CppTools::CheckInputFile(simInputFileName);

         unsigned long currentConfigurationNumberOfEvents = static_cast<unsigned long>
            ((static_cast<TTree *>(TFile::Open(simInputFileName.c_str())->
                                   Get("Tree"))->GetEntries()));
         if (currentConfigurationNumberOfEvents <= 0)
         {
            CppTools::PrintError("Number of events is equal or less than 0 in file " + 
                                 simInputFileName);
         }
         numberOfEvents += currentConfigurationNumberOfEvents;

         if (inputYAMLResonance["has_antiparticle"].as<bool>())
         {
            simInputFileName = 
               "data/SimTrees/" + runName + "/Resonance/" + 
               inputYAMLResonance["name"].as<std::string>() + "_" + 
               ParticleMap::name[-1*inputYAMLResonance["daughter2_id"].as<int>()] + 
               ParticleMap::name[-1*inputYAMLResonance["daughter1_id"].as<int>()] + 
               "_" + pTRange["name"].as<std::string>() + 
               magneticField["name"].as<std::string>() + ".root";

            CppTools::CheckInputFile(simInputFileName);

            currentConfigurationNumberOfEvents = static_cast<unsigned long>
               ((static_cast<TTree *>(TFile::Open(simInputFileName.c_str())->
                                      Get("Tree"))->GetEntries()));
            if (currentConfigurationNumberOfEvents <= 0)
            {
               CppTools::PrintError("Number of events is equal or less than 0 in file " + 
                                    simInputFileName);
            }
         }
         numberOfEvents += currentConfigurationNumberOfEvents;
      }
   }

   std::vector<std::string> magneticFieldsList;
   for (const auto& magneticField : inputYAMLMain["magnetic_field_configurations"])
   {
      magneticFieldsList.emplace_back(magneticField["name"].as<std::string>());
   }

   std::vector<std::string> pTRangesList;
   for (const auto& pTRange : inputYAMLResonance["pt_ranges"])
   {
      pTRangesList.emplace_back(pTRange["name"].as<std::string>());
   }

   CppTools::Box box{"Parameters"};
 
   box.AddEntry("Run name", runName);
   box.AddEntry("Particle", inputYAMLResonance["name"].as<std::string>());
   if (magneticFieldsList.size() == 1 && magneticFieldsList[0] == "")
   {
      box.AddEntry("Magnetic field", "run default");
   }
   else box.AddEntry("Magnetic fields", magneticFieldsList);
   box.AddEntry("pT ranges", pTRangesList);
   box.AddEntry("Charged track minimum pT [GeV/c]", pTMin);
   box.AddEntry("Charged track maximum pT [GeV/c]", pTMax);
   box.AddEntry("Reweight for pT spectra", reweightForSpectra);
   box.AddEntry("Number of threads", numberOfThreads);
   box.AddEntry("Number of events to be analyzed, 1e6", static_cast<double>(numberOfEvents)/1e6, 3);
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
      CppTools::PrintInfo("AnalyzerResonance has finished processing simulated data");
   };
 
   std::thread pBarThread(pBarCall);
 
   ThrContainer thrContainer;

   for (const auto& magneticField : inputYAMLMain["magnetic_field_configurations"])
   {
      for (const auto& pTRange : inputYAMLResonance["pt_ranges"])
      {
         AnalyzeConfiguration(thrContainer, inputYAMLResonance["name"].as<std::string>(), 
                              inputYAMLResonance["daughter1_id"].as<int>(),
                              inputYAMLResonance["daughter2_id"].as<int>(),
                              magneticField["name"].as<std::string>(), 
                              pTRange["name"].as<std::string>());
         if (inputYAMLResonance["has_antiparticle"].as<bool>())
         {
            AnalyzeConfiguration(thrContainer, inputYAMLResonance["name"].as<std::string>(), 
                                 -1*inputYAMLResonance["daughter2_id"].as<int>(),
                                 -1*inputYAMLResonance["daughter1_id"].as<int>(),
                                 magneticField["name"].as<std::string>(), 
                                 pTRange["name"].as<std::string>());
         }
      }
   }


   isProcessFinished = true;
   pBarThread.join();

   // writing the result
   std::string outputFileName = "data/PostSim/" + runName + "/Resonance/" + 
                                inputYAMLResonance["name"].as<std::string>() + ".root";
   thrContainer.Write(outputFileName);

   return 0;
}

ThrContainerCopy AnalyzeResonance::ThrContainer::GetCopy()
{
   ThrContainerCopy copy;

   copy.distrOrigPT = distrOrigPT->Get();
   copy.distrOrigPTVsRecPT = distrOrigPTVsRecPT.Get();
   copy.distrMInv = distrMInv.Get();
   copy.distrMInvNoPID = distrMInvNoPID.Get();
   copy.distrMInvNoPIDOneArmAntiCut = distrMInvNoPIDOneArmAntiCut.Get();
   copy.distrMInvNoPIDGhostAntiCut = distrMInvNoPIDGhostAntiCut.Get();
   copy.distrMInvNoPIDSailorCut = distrMInvNoPIDSailorCut.Get();
   copy.distrMInvNoPIDCowboyCut = distrMInvNoPIDCowboyCut.Get();
   copy.distrPAsymVsPT = distrPAsymVsPT.Get();
   copy.distrDPhiVsPT = distrDPhiVsPT.Get();
   copy.distrDAlphaVsPT = distrDAlphaVsPT.Get();
   copy.distrDZedVsPT = distrDZedVsPT.Get();
   copy.distrDPC2PhiDPC2ZVsPT = distrDPC2PhiDPC2ZVsPT.Get();
   copy.distrDPC3PhiDPC3ZVsPT = distrDPC3PhiDPC3ZVsPT.Get();
   copy.distrDChamberDSlatVsPT = distrDChamberDSlatVsPT.Get();
   copy.distrDChamberDStripVsPT = distrDChamberDStripVsPT.Get();
   copy.distrDYTowerDZTowerVsPT = distrDYTowerDZTowerVsPT.Get();

   return copy;
}

void AnalyzeResonance::ThrContainer::Write(const std::string& outputFileName)
{
   TFile outputFile(outputFileName.c_str(), "RECREATE");
   outputFile.SetCompressionLevel(6);
   outputFile.cd();

   static_cast<std::shared_ptr<TH1D>>(distrOrigPT->Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrOrigPTVsRecPT.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrMInv.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrMInvNoPID.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrMInvNoPIDOneArmAntiCut.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrMInvNoPIDGhostAntiCut.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrMInvNoPIDSailorCut.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrMInvNoPIDCowboyCut.Merge())->Write();
   static_cast<std::shared_ptr<TH3F>>(distrPAsymVsPT.Merge())->Write();
   static_cast<std::shared_ptr<TH3F>>(distrDPhiVsPT.Merge())->Write();
   static_cast<std::shared_ptr<TH3F>>(distrDAlphaVsPT.Merge())->Write();
   static_cast<std::shared_ptr<TH3F>>(distrDZedVsPT.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrDPC2PhiDPC2ZVsPT.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrDPC3PhiDPC3ZVsPT.Merge())->Write();
   static_cast<std::shared_ptr<TH3F>>(distrDChamberDSlatVsPT.Merge())->Write();
   static_cast<std::shared_ptr<TH3F>>(distrDChamberDStripVsPT.Merge())->Write();
   static_cast<std::shared_ptr<TH3F>>(distrDYTowerDZTowerVsPT.Merge())->Write();

   outputFile.Close();
}

#endif /* ANALYZE_RESONANCE_CPP */
