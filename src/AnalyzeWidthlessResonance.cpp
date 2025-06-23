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

// this namespace is only used so that documentation will not become a mess
// so there is no need to enforce the contents inside of it 
// being accessed only via the scope resolution operator in this file
using namespace AnalyzeWidthlessResonance;

void AnalyzeWidthlessResonance::AnalyzeConfiguration(ThrContainer &thrContainer, 
                                                     const std::string& particleName, 
                                                     const short particleGeantId, 
                                                     const std::string& magneticFieldName, 
                                                     const std::string &pTRangeName)
{ 
   std::string simInputFileName = "data/SimTrees/" + runName + "/WidthlessResonance/" + 
                                  particleName + "_" + pTRangeName + magneticFieldName + ".root";

   TFile simInputFile = TFile(simInputFileName.c_str());

   TH1F *origPTHist = static_cast<TH1F *>(simInputFile.Get("orig_pt"));

   // weight function for spectra
   std::unique_ptr<TF1> weightFunc;
  
   // normalization of the number of particles to the number of events
   // this normalization is needed to seamlessly merge 2 files with 
   // flat pT distribution with different ranges
   double eventNormWeight = 1.;
   if (reweightForSpectra)
   {
      std::string realDataFileName = "data/Real/" + runName + "/WidthlessResonance/sum" + 
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

            switch (charge)
            {
               case 1:
                  positiveTracks.emplace_back(simCNT.mom[i]*sin(theta0)*cos(simCNT.phi0(i)), 
                                              simCNT.mom[i]*sin(theta0)*sin(simCNT.phi0(i)), 
                                              simCNT.mom[i]*cos(theta0), phi, alpha, zed);
                  break;
               case -1:
                  negativeTracks.emplace_back(simCNT.mom[i]*sin(theta0)*cos(simCNT.phi0(i)), 
                                              simCNT.mom[i]*sin(theta0)*sin(simCNT.phi0(i)), 
                                              simCNT.mom[i]*cos(theta0), phi, alpha, zed);
                  break;
            }
         }

         for (const auto& posTrack : positiveTracks)
         {
            for (const auto& negTrack : negativeTracks)
            {
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

   inputYAMLSim.OpenFile(argv[1]);
   inputYAMLSim.CheckStatus("resonance");

   runName = inputYAMLSim["run_name"].as<std::string>();

   inputYAMLMain.OpenFile("input/" + runName + "/main.yaml");
   inputYAMLMain.CheckStatus("main");

   inputYAMLSimSingleTrack.OpenFile("input/" + runName + "/single_track.yaml");
   inputYAMLSimSingleTrack.CheckStatus("single_track");
 
   collisionSystemName = inputYAMLMain["collision_system_name"].as<std::string>();

   outputDir = "data/PostSim/" + runName + "/WidthlessResonance/";
   system(("mkdir -p " + outputDir).c_str());

   pTMin = inputYAMLSimSingleTrack["pt_min"].as<double>();
   pTMax = inputYAMLSimSingleTrack["pt_max"].as<double>();

   reweightForSpectra = inputYAMLSim["reweight_for_spectra"].as<bool>();

   dmCutter.Initialize(runName, inputYAMLMain["detectors_configuration"].as<std::string>());

   if (reweightForSpectra)
   {
      CppTools::CheckInputFile("data/Parameters/SpectraFit/" + collisionSystemName + 
                               "/" + inpytYAMLSim["name"].as<std::string>() + ".yaml");
   }
 
   for (const auto& magneticField : inputYAMLMain["magnetic_field_configurations"])
   {
      for (const auto& pTRange : inputYAMLSim["pt_ranges"])
      {
         const std::string simInputFileName = 
            "data/SimTrees/" + runName + "/WidthlessResonance/" + 
            inputYAMLSim["name"].as<std::string>() + "_" + pTRange["name"].as<std::string>() + 
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
   box.AddEntry("Particle", inputYAMLSim["name"]);
   if (magneticFieldsList.size() == 1 && magneticFieldsList[0] == "")
   {
      box.AddEntry("Magnetic field", "run default");
   }
   else box.AddEntry("Magnetic fields", magneticFieldsList);
   box.AddEntry("pT ranges", pTRangesList);
   box.AddEntry("Charged track minimum pT, GeV", pTMin);
   box.AddEntry("Charged track maximum pT, GeV", pTMax);
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
      for (const auto& pTRange : inputYAMLSim["pt_ranges"])
      {
         AnalyzeConfiguration(thrContainer, inputYAMLSim["name"].as<std::string>(), 
                              magneticField["name"].as<std::string>(), 
                              pTRange["name"].as<std::string>());
      }
   }


   isProcessFinished = true;
   pBarThread.join();

   // writing the result
   std::string outputFileName = "data/PostSim/" + runName + "/WidthlessResonance/" + 
                                inputYAMLSim["name"].as<std::string>() + ".root";
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
