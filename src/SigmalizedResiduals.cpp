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
   inputYAMLMain.OpenFile(argv[1], "single_track_sim");
   inputYAMLMain.CheckStatus("single_track_sim");

   runName = inputYAMLMain["run_name"].as<std::string>();

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

   const std::string detectorsConfiguration = 
      inputYAMLMain["detectors_configuration"].as<std::string>();

   if (detectorsConfiguration[2] == '1') 
   {
      PerformCalibrationForDetector("PC2");
   }
   if (detectorsConfiguration[3] == '1') 
   {
      PerformCalibrationForDetector("PC3e");
      PerformCalibrationForDetector("PC3w");
   }
   if (detectorsConfiguration[4] == '1') 
   {
      PerformCalibrationForDetector("TOFe");
   }
   if (detectorsConfiguration[5] == '1') 
   {
      PerformCalibrationForDetector("TOFw");
   }
   if (detectorsConfiguration[6] == '1') 
   {
      for (int i = 0; i < 4; i++)
      {
         PerformCalibrationForDetector("EMCale" + std::to_string(i));
         PerformCalibrationForDetector("EMCalw" + std::to_string(i));
      }
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

   TGraphErrors meansDVal;
   TGraphErrors sigmasDVal;

   const double xMin = 
      distrDValVsPT->GetXaxis()->GetBinLowEdge(1);
   const double maxX = 
      distrDValVsPT->GetXaxis()->GetBinUpEdge(distrDValVsPT->GetXaxis()->GetNbins());

   for (int i = 1; i < distrDValVsPT->GetYaxis()->FindBin(3.); i++) // pT < 3 GeV/c
   {
      if (distrDValVsPT->Integral(1, distrDValVsPT->GetXaxis()->GetNbins(), i, i) < 1e-15) continue;

      TH1F *distrDValProj = 
         distrDValVsPT->ProjectionX((distrDValVsPT->GetName() + std::to_string(i)).c_str(), i, i);

      TF1 fitDistrDVal((detectorName + variableName + chargeName + std::to_string(i)).c_str(), 
                       "gaus(0) + gaus(3)");

      const double maxBinVal = distrDValProj->GetBinContent(distrDValProj->GetMaximumBin());

      fitDistrDVal.SetParameters(maxBinVal, 0., distrDValProj->GetXaxis()->GetBinWidth(1)*10.);

      fitDistrDVal.SetParLimits(0, maxBinVal/2., maxBinVal);
      fitDistrDVal.SetParLimits(1, xMin/5., xMax/5.);
      fitDistrDVal.SetParLimits(2, distrDValProj->GetXaxis()->GetBinWidth(1), xMax/2.);
      fitDistrDVal.SetParLimits(4, xMax/3., xMax*3.);

      fitDistrDVal.SetRange(xMin, xMax);

      for (int i = 1; i <= fitNTries; i++)
      {
         distrDValProj->Fit(&fitDistrDVal, "RQMBN");

         fitDistrDVal.SetParLimits(0, fitdistrDVal.Getparameter(0)/
                                   (1. + 2./static_cast<double>(i*i*i)),
                                   fitdistrDVal.Getparameter(0)*
                                   (1. + 2./static_cast<double>(i*i*i)));
         fitDistrDVal.SetParLimits(1, fitdistrDVal.Getparameter(1)*
                                   (1. - 6./static_cast<double>(i*i*i)),
                                   fitdistrDVal.Getparameter(1)*
                                   (1. + 4./static_cast<double>(i*i*i)));
         fitDistrDVal.SetParLimits(2, fitdistrDVal.Getparameter(2)/
                                   (1. + 5./static_cast<double>(i*i*i)),
                                   fitdistrDVal.Getparameter(2)*
                                   (1. + 5./static_cast<double>(i*i*i)));
         fitDistrDVal.SetParLimits(3, fitdistrDVal.Getparameter(3)/
                                   (1. + 2./static_cast<double>(i*i*i)),
                                   fitdistrDVal.Getparameter(3)*
                                   (1. + 2./static_cast<double>(i*i*i)));
         fitDistrDVal.SetParLimits(4, fitdistrDVal.Getparameter(4)*
                                   (1. - 6./static_cast<double>(i*i*i)),
                                   fitdistrDVal.Getparameter(4)*
                                   (1. + 4./static_cast<double>(i*i*i)));
         fitDistrDVal.SetParLimits(5, fitdistrDVal.Getparameter(5)/
                                   (1. + 5./static_cast<double>(i*i*i)),
                                   fitdistrDVal.Getparameter(5)*
                                   (1. + 5./static_cast<double>(i*i*i)));

         fitDistrDVal.SetRange(fitDistrDVal.GetParameter(1) - 5.*fitDistrDVal.GetParameter(2), 
                               fitDistrDVal.GetParameter(1) + 5.*fitDistrDVal.GetParameter(2));
      }

      distrDValProj->Fit(&fitDistrDVal, "RQMB");
      
      meansDVal.AddPoint(distrDValVsPT->GetXaxis()->GetBinCenter(i), fitDValProj.GetParameter(1));
      sigmasDVal.AddPoint(distrDValVsPT->GetXaxis()->GetBinCenter(i), fitDValProj.GetParameter(2));
   }
}

#endif /* SIGMALIZED_RESIDUALS_CPP */
