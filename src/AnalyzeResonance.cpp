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
      weightFunc = std::make_unique<TF1>("weightFunc", "exp(-x)");
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

         for(int i = 0; i < simCNT.nch(); i++) // loop over particles in one event
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

            int idPC2 = PART_ID::JUNK;
            int idPC3 = PART_ID::JUNK;
            int idEMCal = PART_ID::JUNK;
            int idTOFe = PART_ID::JUNK;
            int idTOFw = PART_ID::JUNK;

            double weightPC2 = 0.;
            double weightPC3 = 0.;
            double weightEMCal = 0.;
            double weightTOFe = 0.;
            double weightTOFw = 0.;

            double weightIdEMCal = 0.;
            double weightIdTOFe = 0.;
            double weightIdTOFw = 0.;

            if (IsHit(simCNT.pc2dphi(i)))
            {
               const double sdphi = simSigmRes.PC2SDPhi(simCNT.pc2dphi(i), pT, charge);
               const double sdz = simSigmRes.PC2SDZ(simCNT.pc2dz(i), pT, charge);
               const double pc2phi = atan2(simCNT.ppc2y(i), simCNT.ppc2x(i));

               if (IsMatch(pT, sdphi, sdz) && !dmCutter.IsDeadPC2(simCNT.ppc2z(i), pc2phi))
               {
                  idPC2 = PART_ID::NONE;
                  weightPC2 = 1.;
               }
            }

            if (IsHit(simCNT.pc3dphi(i)))
            {
               const double sdphi = simSigmRes.PC3SDPhi(simCNT.pc3dphi(i), pT, charge, dcarm);
               const double sdz = simSigmRes.PC3SDZ(simCNT.pc3dz(i), pT, charge, dcarm);

               double pc3phi = atan2(simCNT.ppc3y(i), simCNT.ppc3x(i));
               if (dcarm == 0 && pc3phi < 0) pc3phi += 2.*M_PI;

               if (IsMatch(pT, sdphi, sdz) && !dmCutter.IsDeadPC3(dcarm, simCNT.ppc2z(i), pc3phi))
               {
                  idPC3 = PART_ID::NONE;
                  weightPC3 = 1.;
               }
            }

            if (IsHit(simCNT.emcdz(i)))
            {
               const double sdphi = 
                  simSigmRes.EMCalSDPhi(simCNT.emcdphi(i), pT, charge, dcarm, simCNT.sect(i));
               const double sdz = 
                  simSigmRes.EMCalSDZ(simCNT.emcdz(i), pT, charge, dcarm, simCNT.sect(i));

               bool isCutByECore;
               if (dcarm == 0 && simCNT.sect(i) < 2) isCutByECore = (simCNT.ecore(i) < 0.35);
               else isCutByECore = (simCNT.ecore(i) < 0.25); // PbSc

               if (IsMatch(pT, sdphi, sdz) && !isCutByECore && 
                   !dmCutter.IsDeadEMCal(dcarm, simCNT.sect(i), simCNT.ysect(i), simCNT.zsect(i)))
               {
                  weightEMCal = 1.;
                  if (!(dcarm == 0 && simCNT.sect(i) < 2) &&
                      !dmCutter.IsDeadTimingEMCal(dcarm, simCNT.sect(i), 
                                                  simCNT.ysect(i), simCNT.zsect(i)))
                  {
                     switch (charge)
                     {
                        case 1:
                           idEMCal = daughter1Id;
                           break;
                        case -1:
                           idEMCal = daughter12d;
                           break;
                     }
                     weightIdEMCal = simM2Id.GetEMCalIdProb(simCNT.dcarm(i), simCNT.sect(i), 
                                                            pT, idEMCal, 2., 2.)*weightEMCal;
                     if (weightIdEMCal <= 0.)
                     {
                        idEMCal = PART_ID::NONE;
                        weightIdEMCal = 0.;
                     }
                  }
                  else idEMCal = PART_ID::NONE;
               }
            }

            if (IsHit(simCNT.tofdz(i)))
            {
               const double sdphi = simSigmRes.TOFeSDPhi(simCNT.tofdphi(i), pT, charge);
               const double sdz = simSigmRes.TOFeSDZ(simCNT.tofdz(i), pT, charge);

               const double beta = simCNT.pltof(i)/simCNT.ttof(i)/29.9792;
               const double eloss = 0.0005*pow(beta, -2.5);

               // slats are organized in 10 lines of 96 we define as chambers
               const int chamber = simCNT.slat(i)/96;
               // slat number for the current chamber
               const int slat = simCNT.slat(i) % 96;

               if (simCNT.etof(i) > eloss && IsMatch(pT, sdphi, sdz) && 
                   !dmCutter.IsDeadTOFe(chamber, slat))
               {
                  weightTOFe = 1.;
                  if (!dmCutter.IsDeadTimingTOFe(chamber, slat))
                  {
                     switch (charge)
                     {
                        case 1:
                           idTOFe = daughter1Id;
                           break;
                        case -1:
                           idTOFe = daughter2Id;
                           break;
                     }
                     weightIdTOFe = simM2Id.GetTOFeIdProb(idTOFe, pT, 2., 2.)*weightTOFe;

                     if (weightIdTOFe <= 0.)
                     {
                        idTOFe = PART_ID::NONE;
                        weightIdTOFe = 0.;
                     }
                  }
                  else idTOFe = PART_ID::NONE;
               }
            }
            else if (IsHit(simCNT.tofwdz(i)))
            {
               const double sdphi = simSigmRes.TOFwSDPhi(simCNT.tofwdphi(i), pT, charge);
               const double sdz = simSigmRes.TOFwSDZ(simCNT.tofwdz(i), pT, charge);

               // strips are organized in 8 lines of 64 we define as chambers
               const int chamber = simCNT.striptofw(i)/64;
               // strip number for the current chamber
               const int strip = simCNT.striptofw(i) % 64;

               if (IsMatch(pT, sdphi, sdz) && !dmCutter.IsDeadTOFw(chamber, strip))
               {
                  weightTOFw = 1.;
                  if (!dmCutter.IsDeadTimingTOFw(chamber, strip))
                  {
                     switch (charge)
                     {
                        case 1:
                           idTOFw = daughter1Id;
                           break;
                        case -1:
                           idTOFw = daughter2Id;
                           break;
                     }
                     weightIdTOFw = simM2Id.GetTOFwIdProb(idTOFw, pT, 2., 2.)*weightTOFw;

                     if (weightIdTOFw <= 0.)
                     {
                        idTOFw = PART_ID::NONE;
                        weightIdTOFw = 0.;
                     }
                  }
                  else idTOFw = PART_ID::NONE;
               }
            }

            if (idPC2 == PART_ID::JUNK && idPC3 == PART_ID::JUNK && idEMCal == PART_ID::JUNK && 
                idTOFe == PART_ID::JUNK && idTOFw == PART_ID::JUNK) continue;

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

         // looping over pairs of tracks
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

               if (IsGhostCut(posTrack, negTrack)) continue;

               if (IsOneArmCut(posTrack, negTrack)) 
               {
                  thrContainer.distrMInvOneArmAntiCut->Fill(pT, mInv, eventWeight);
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
                  if (posTrack.idPC2 != PART_ID::JUNK && negTrack.idPC2 != PART_ID::JUNK)
                  {
                     thrContainer.distrDPC2PhiDPC2ZVsPT->Fill(posTrack.pc2z - negTrack.pc2z, 
                                                              posTrack.pc2phi - negTrack.pc2phi,
                                                              eventWeight);
                  }

                  if (posTrack.idPC3 != PART_ID::JUNK && negTrack.idPC3 != PART_ID::JUNK)
                  {
                     thrContainer.distrDPC3PhiDPC3ZVsPT->Fill(posTrack.pc3z - negTrack.pc3z, 
                                                              posTrack.pc3phi - negTrack.pc3phi, 
                                                              eventWeight);
                  }

                  if (posTrack.idTOFe != PART_ID::JUNK && negTrack.idTOFe != PART_ID::JUNK)
                  {
                     thrContainer.distrDChamberDSlatVsPT->
                        Fill(static_cast<double>(posTrack.slat/96 - negTrack.slat/96) + 0.5,
                             static_cast<double>((posTrack.slat % 96) - 
                                                 (negTrack.slat % 96)) + 0.5,
                             pT, eventWeight);
                  }
                  else if (posTrack.idTOFw != PART_ID::JUNK && negTrack.idTOFw != PART_ID::JUNK)
                  {
                     thrContainer.distrDChamberDStripVsPT->
                        Fill(static_cast<double>(posTrack.strip/96 - negTrack.strip/96) + 0.5,
                             static_cast<double>((posTrack.strip % 96) - 
                                                 (negTrack.strip % 96)) + 0.5,
                             pT, eventWeight);
                  }

                  if (posTrack.idEMCal != PART_ID::JUNK && negTrack.idEMCal != PART_ID::JUNK &&
                      posTrack.sector == negTrack.sector)
                  {
                     thrContainer.distrDYTowerDZTowerVsPT->
                        Fill(static_cast<double>(posTrack.yTower - negTrack.yTower) + 0.5, 
                             static_cast<double>(posTrack.zTower - negTrack.zTower) + 0.5, 
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

               if (Is1PID(posTrack, negTrack, daughter1Id, daughter2Id))
               {
                  thrContainer.distrMInv1PID->Fill(pT, mInv, eventWeight);
               }

               if (!Is2PID(posTrack, negTrack, daughter1Id, daughter2Id)) continue;

               thrContainer.distrMInv2PID->Fill(pT, mInv, eventWeight);

               if (IsTOF2PID(posTrack, negTrack, daughter1Id, daughter2Id))
               {
                  thrContainer.distrMInvTOF2PID->Fill(pT, mInv, eventWeight);
               }
               if (IsEMCal2PID(posTrack, negTrack, daughter1Id, daughter2Id))
               {
                  thrContainer.distrMInvEMCal2PID->Fill(pT, mInv, eventWeight);
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
   simSigmRes.Initialize(runName, inputYAMLMain["detectors_configuration"].as<std::string>());
   simM2Id.Initialize(runName, true);

   if (CppTools::FileExists("data/Parameters/SpectraFit/" + collisionSystemName + 
                            "/" + inputYAMLResonance["name"].as<std::string>() + ".yaml"))
   {
      CppTools::PrintInfo("Fit parameters for spectra were found");
      reweightForSpectra = true;
   }
   else 
   {
      CppTools::PrintInfo("Fit parameters for spectra were not found;" \
                          "setting reweight to e^{-p_{T}}");
      reweightForSpectra = false;
   }
 
   for (const auto& magneticField : inputYAMLMain["magnetic_field_configurations"])
   {
      CppTools::CheckInputFile("data/Real/" + runName + "/SingleTrack/sum" + 
                               magneticField["name"].as<std::string>() + ".root");
      for (const auto& pTRange : inputYAMLResonance["sim_pt_ranges"])
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
   for (const auto& pTRange : inputYAMLResonance["sim_pt_ranges"])
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
      for (const auto& pTRange : inputYAMLResonance["sim_pt_ranges"])
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
   copy.distrMInvOneArmAntiCut = distrMInvOneArmAntiCut.Get();
   copy.distrMInvNoPID = distrMInvNoPID.Get();
   copy.distrMInv1PID = distrMInv1PID.Get();
   copy.distrMInv2PID = distrMInv2PID.Get();
   copy.distrMInvTOF2PID = distrMInvTOF2PID.Get();
   copy.distrMInvEMCal2PID = distrMInvEMCal2PID.Get();
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
   static_cast<std::shared_ptr<TH2F>>(distrMInvOneArmAntiCut.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrMInvNoPID.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrMInv1PID.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrMInv2PID.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrMInvTOF2PID.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrMInvEMCal2PID.Merge())->Write();
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
