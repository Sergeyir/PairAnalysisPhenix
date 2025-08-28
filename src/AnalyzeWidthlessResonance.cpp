/** 
 *  @file   AnalyzeWidthlessResonance.hpp
 *  @brief  Contains declarations of functions and variables that are used for analysis of a widthless resonance from a trees acquired from the PHENIX simulation. Widthless resonances are used for determining the determination of a mass resolution of a detector system for a given resonance.
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysisPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef ANALYZE_WIDTHLESS_RESONANCE_CPP
#define ANALYZE_WIDTHLESS_RESONANCE_CPP

#include "../include/AnalyzeWidthlessResonance.hpp"

// this namespace is only used so that documentation does not become a mess
// so there is no need to enforce the contents inside of it 
// being accessed only via the scope resolution operator in this file
using namespace AnalyzeWidthlessResonance;

void AnalyzeWidthlessResonance::AnalyzeConfiguration(ThrContainer &thrContainer, 
                                                     const std::string& particleName, 
                                                     const int daughter1Id,
                                                     const int daughter2Id,
                                                     const std::string& magneticFieldName, 
                                                     const std::string &pTRangeName)
{ 
   std::string simInputFileName = "data/SimTrees/" + runName + "/WidthlessResonance/" + 
                                  particleName + "_" + ParticleMap::name[daughter1Id] + 
                                  ParticleMap::name[daughter2Id] + "_" + 
                                  pTRangeName + magneticFieldName + ".root";

   TFile simInputFile = TFile(simInputFileName.c_str());

   TH1F *origPTHist = static_cast<TH1F *>(simInputFile.Get("orig_pt"));

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

   TH1F *centrHist = static_cast<TH1F *>(realDataFile.Get("centrality"));
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

            histContainer.distrOrigPTVsRecPT->Fill(origPT, pT, eventWeight);

            int idPC2 = PART_ID::JUNK;
            int idPC3 = PART_ID::JUNK;
            int idEMCal = PART_ID::JUNK;
            int idTOFe = PART_ID::JUNK;
            int idTOFw = PART_ID::JUNK;

            if (IsHit(simCNT.pc2dphi(i)))
            {
               const double sdphi = simCalibrator.PC2SDPhi(simCNT.pc2dphi(i), pT, charge);
               const double sdz = simCalibrator.PC2SDZ(simCNT.pc2dz(i), pT, charge);
               const double pc2phi = atan2(simCNT.ppc2y(i), simCNT.ppc2x(i));

               if (IsMatch(pT, sdphi, sdz) && !dmCutter.IsDeadPC2(simCNT.ppc2z(i), pc2phi))
               {
                  idPC2 = PART_ID::NONE;
               }
            }

            if (IsHit(simCNT.pc3dphi(i)))
            {
               const double sdphi = simCalibrator.PC3SDPhi(simCNT.pc3dphi(i), pT, charge, dcarm);
               const double sdz = simCalibrator.PC3SDZ(simCNT.pc3dz(i), pT, charge, dcarm);

               double pc3phi = atan2(simCNT.ppc3y(i), simCNT.ppc3x(i));
               if (dcarm == 0 && pc3phi < 0) pc3phi += 2.*M_PI;

               if (IsMatch(pT, sdphi, sdz) && !dmCutter.IsDeadPC3(dcarm, simCNT.ppc2z(i), pc3phi))
               {
                  idPC3 = PART_ID::NONE;
               }
            }

            if (IsHit(simCNT.emcdz(i)))
            {
               const double sdphi = 
                  simCalibrator.EMCalSDPhi(simCNT.emcdphi(i), pT, charge, dcarm, simCNT.sect(i));
               const double sdz = 
                  simCalibrator.EMCalSDZ(simCNT.emcdz(i), pT, charge, dcarm, simCNT.sect(i));

               bool isCutByECore;
               if (dcarm == 0 && simCNT.sect(i) < 2) isCutByECore = (simCNT.ecore(i) < 0.35);
               else isCutByECore = (simCNT.ecore(i) < 0.25); // PbSc

               if (IsMatch(pT, sdphi, sdz) && !isCutByECore && 
                   !dmCutter.IsDeadEMCal(dcarm, simCNT.sect(i), simCNT.ysect(i), simCNT.zsect(i)))
               {
                  idEMCal = PART_ID::NONE;
               }
            }

            if (IsHit(simCNT.tofdz(i)))
            {
               const double sdphi = simCalibrator.TOFeSDPhi(simCNT.tofdphi(i), pT, charge);
               const double sdz = simCalibrator.TOFeSDZ(simCNT.tofdz(i), pT, charge);

               const double beta = simCNT.pltof(i)/simCNT.ttof(i)/29.9792;
               const double eloss = 0.0005*pow(beta, -2.5);

               // slats are organized in 10 lines of 96 we define as chambers
               const int chamber = simCNT.slat(i)/96;
               // slat number for the current chamber
               const int slat = simCNT.slat(i) % 96;

               if (simCNT.etof(i) > eloss && IsMatch(pT, sdphi, sdz) && 
                   !dmCutter.IsDeadTOFe(chamber, slat))
               {
                  idTOFe = PART_ID::NONE;
               }
            }
            else if (IsHit(simCNT.tofwdz(i)))
            {
               const double sdphi = simCalibrator.TOFwSDPhi(simCNT.tofwdphi(i), pT, charge);
               const double sdz = simCalibrator.TOFwSDZ(simCNT.tofwdz(i), pT, charge);

               // strips are organized in 8 lines of 64 we define as chambers
               const int chamber = simCNT.striptofw(i)/64;
               // strip number for the current chamber
               const int strip = simCNT.striptofw(i) % 64;

               if (IsMatch(pT, sdphi, sdz) && !dmCutter.IsDeadTOFw(chamber, strip))
               {
                  idTOFw = PART_ID::NONE;
               }
            }

            if (idPC2 == PART_ID::JUNK && idPC3 == PART_ID::JUNK && 
                idEMCal == PART_ID::JUNK && 
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

         for (const auto& posTrack : positiveTracks)
         {
            for (const auto& negTrack : negativeTracks)
            {
               if (IsOneArmCut(posTrack, negTrack) || IsGhostCut(posTrack, negTrack) ||
                   !IsNoPID(posTrack, negTrack)) continue;

               const double mInv = GetPairMass(posTrack, negTrack);
               const double pT = GetPairPT(posTrack, negTrack);

               thrContainer.distrMInvNoPID->Fill(pT, mInv, eventWeight);
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
      errMsg += "parameter(s) were provided \n Usage: bin/AnalyzeWidthlessResonance ";
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

   outputDir = "data/PostSim/" + runName + "/WidthlessResonance/";
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
      CppTools::PrintInfo("Fit parameters for spectra were not found; setting reweight to e^{-1*x}");
      reweightForSpectra = false;
   }
 
   for (const auto& magneticField : inputYAMLMain["magnetic_field_configurations"])
   {
      CppTools::CheckInputFile("data/Real/" + runName + "/SingleTrack/sum" + 
                               magneticField["name"].as<std::string>() + ".root");
      for (const auto& pTRange : inputYAMLResonance["pt_ranges"])
      {
         std::string simInputFileName = 
            "data/SimTrees/" + runName + "/WidthlessResonance/" + 
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
               "data/SimTrees/" + runName + "/WidthlessResonance/" + 
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
      CppTools::PrintInfo("AnalyzerWidthlessResonance has finished processing simulated data");
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
   std::string outputFileName = "data/PostSim/" + runName + "/WidthlessResonance/" + 
                                inputYAMLResonance["name"].as<std::string>() + ".root";
   thrContainer.Write(outputFileName);

   return 0;
}

ThrContainerCopy AnalyzeWidthlessResonance::ThrContainer::GetCopy()
{
   ThrContainerCopy copy;

   copy.distrOrigPT = distrOrigPT->Get();
   copy.distrOrigPTVsRecPT = distrOrigPTVsRecPT.Get();
   copy.distrMInvNoPID = distrMInvNoPID.Get();

   return copy;
}

void AnalyzeWidthlessResonance::ThrContainer::Write(const std::string& outputFileName)
{
   TFile outputFile(outputFileName.c_str(), "RECREATE");
   outputFile.SetCompressionLevel(6);
   outputFile.cd();

   static_cast<std::shared_ptr<TH1F>>(distrOrigPT->Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrOrigPTVsRecPT.Merge())->Write();
   static_cast<std::shared_ptr<TH2F>>(distrMInvNoPID.Merge())->Write();

   outputFile.Close();
}

#endif /* ANALYZE_WIDTHLESS_RESONANCE_CPP */
