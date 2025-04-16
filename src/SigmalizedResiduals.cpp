/** 
 *  @file   SigmalizedResiduals.cpp 
 *  @brief  Contains realisation of functions that are used for estimation of values for calibration of sigmalized residuals dphi and dz from the PHENIX simulation
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysisPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef SIGMALIZED_RESIDUALS_CPP
#define SIGMALIZED_RESIDUALS_CPP

#include "../include/SigmalizedResiduals.hpp"

int main(int argc, char **argv)
{
   using namespace SigmalizedResiduals;

   if (argc < 2 || argc > 3) 
   {
      std::string errMsg = "Expected 1-2 parameters while " + std::to_string(argc - 1) + 
                           " parameter(s) were provided \n";
      errMsg += "Usage: bin/SigmalizedResiduals inputFile numberOfThreads=" + 
                std::to_string(std::thread::hardware_concurrency()) + "*\n";
      errMsg += "*: default argument is the number of threads on the current machine \n";
      CppTools::PrintError(errMsg);
   }

   // initializing this program parameters
   inputYAMLCal.OpenFile(argv[1], "single_track_sim");
   inputYAMLCal.CheckStatus("single_track_sim");

   runName = inputYAMLCal["run_name"].as<std::string>();

   // opening input file with parameters of a run
   inputYAMLMain.OpenFile("input/" + runName + "/main.yaml");
   inputYAMLMain.CheckStatus("main");

   unsigned int numberOfThreads;
   if (argc > 2) numberOfThreads = std::stoi(argv[2]);
   else numberOfThreads = std::thread::hardware_concurrency();
   if (numberOfThreads == 0) CppTools::PrintError("Number of threads must be bigger than 0");

   // initializing ROOT parameters
   ROOT::EnableThreadSafety();
   gErrorIgnoreLevel = kWarning;
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);

   TDirectory::AddDirectory(kFALSE);

   ROOT::EnableImplicitMT(numberOfThreads);

   outputDir = "output/SigmalizedResiduals/" + runName + "/";
   system(("mkdir -p " + outputDir).c_str());

   inputFile = std::unique_ptr<TFile>(TFile::Open(("data/PostSim/" + runName + 
                                                   "SingleTrack/all.root").c_str(), "READ"));

   PerformCalibrationForDetector("PC2");
   PerformCalibrationForDetector("PC3e");
   PerformCalibrationForDetector("PC3w");
   PerformCalibrationForDetector("PC3e");
   PerformCalibrationForDetector("TOFe");
   PerformCalibrationForDetector("TOFw");

   for (int i = 0; i < 4; i++)
   {
      PerformCalibrationForDetector("EMCale" + std::to_string(i));
      PerformCalibrationForDetector("EMCalw" + std::to_string(i));
   }
 
   return 0;
}

void SigmalizedResiduals::PerformCalibrationForDetector(const std::string& detectorName, 
                                                        const std::string& variableName,
                                                        const int charge)
{
   const std::string chargeName = ((charge > 0) ? "charge>0" : "charge<0");
   const std::string chargeNameShort = ((charge > 0) ? "pos" : "neg");

   TH2F *distrDValVsPT = 
      static_cast<TH2F *>(inputFile->Get((variableName + " vs pT: " + detectorName + 
                                          ", " + chargeName).c_str()));

}

#endif /* SIGMALIZED_RESIDUALS_CPP */
