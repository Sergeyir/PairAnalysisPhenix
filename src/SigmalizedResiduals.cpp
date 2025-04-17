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

// this namespace is only used so that documentation will not become a mess
// so there is no need to enforce the contents inside of it 
// being accessed only via the scope resolution operator in this file
using namespace SigmalizedResiduals;

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
   inputYAMLMain.OpenFile(argv[1], "main");
   inputYAMLMain.CheckStatus("main");

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

   const std::string inputFileName = "data/PostSim/" + runName + "/SingleTrack/all.root";
   CppTools::CheckInputFile(inputFileName);

   inputFile = std::unique_ptr<TFile>(TFile::Open(inputFileName.c_str(), "READ"));

   const std::string detectorsConfiguration = 
      inputYAMLMain["detectors_configuration"].as<std::string>();

   if (detectorsConfiguration[2] == '1') numberOfIterations += 1;
   if (detectorsConfiguration[3] == '1') numberOfIterations += 2;
   if (detectorsConfiguration[4] == '1') numberOfIterations += 1;
   if (detectorsConfiguration[5] == '1') numberOfIterations += 1;
   if (detectorsConfiguration[6] == '1') numberOfIterations += 8;

   numberOfIterations *= 4;

   for (const std::string& variableName : variableNames)
   {
      for (const int charge : particleCharges)
      {
         if (detectorsConfiguration[2] == '1') 
         {
            PerformCalibrationsForDetector("PC2", variableName, charge);
         }
         if (detectorsConfiguration[3] == '1') 
         {
            PerformCalibrationsForDetector("PC3e", variableName, charge);
            PerformCalibrationsForDetector("PC3w", variableName, charge);
         }
         if (detectorsConfiguration[4] == '1') 
         {
            PerformCalibrationsForDetector("TOFe", variableName, charge);
         }
         if (detectorsConfiguration[5] == '1') 
         {
            PerformCalibrationsForDetector("TOFw", variableName, charge);
         }
         if (detectorsConfiguration[6] == '1') 
         {
            for (int i = 0; i < 4; i++)
            {
               PerformCalibrationsForDetector("EMCale" + std::to_string(i), variableName, charge);
               PerformCalibrationsForDetector("EMCalw" + std::to_string(i), variableName, charge);
            }
         }
      }
   }

   pBar.Print(1.);
 
   return 0;
}

void SigmalizedResiduals::PerformCalibrationsForDetector(const std::string& detectorName, 
                                                         const std::string& variableName,
                                                         const int charge)
{
   pBar.Print(static_cast<double>(numberOfCalls)/static_cast<double>(numberOfIterations));

   const std::string chargeName = ((charge > 0) ? "charge>0" : "charge<0");
   const std::string chargeNameShort = ((charge > 0) ? "pos" : "neg");

   TH2F *distrDValVsPT = 
      static_cast<TH2F *>(inputFile->Get((variableName + " vs pT: " + detectorName + 
                                          ", " + chargeName).c_str()));

   TGraphErrors grMeansDVal;
   TGraphErrors grSigmasDVal;

   const double xMin = 
      distrDValVsPT->GetXaxis()->GetBinLowEdge(1);
   const double xMax = 
      distrDValVsPT->GetXaxis()->GetBinUpEdge(distrDValVsPT->GetXaxis()->GetNbins());

   for (int i = 1; i < distrDValVsPT->GetYaxis()->FindBin(3.); i++) // pT < 3 GeV/c
   {
      if (distrDValVsPT->Integral(1, distrDValVsPT->GetXaxis()->GetNbins(), i, i) < 1e-15) continue;

      TH1D *distrDValProj = 
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

      for (unsigned int i = 1; i <= fitNTries; i++)
      {
         distrDValProj->Fit(&fitDistrDVal, "RQMBN");

         fitDistrDVal.SetParLimits(0, fitDistrDVal.GetParameter(0)/
                                   (1. + 2./static_cast<double>(i*i*i)),
                                   fitDistrDVal.GetParameter(0)*
                                   (1. + 2./static_cast<double>(i*i*i)));
         fitDistrDVal.SetParLimits(1, fitDistrDVal.GetParameter(1)*
                                   (1. - 6./static_cast<double>(i*i*i)),
                                   fitDistrDVal.GetParameter(1)*
                                   (1. + 4./static_cast<double>(i*i*i)));
         fitDistrDVal.SetParLimits(2, fitDistrDVal.GetParameter(2)/
                                   (1. + 5./static_cast<double>(i*i*i)),
                                   fitDistrDVal.GetParameter(2)*
                                   (1. + 5./static_cast<double>(i*i*i)));
         fitDistrDVal.SetParLimits(3, fitDistrDVal.GetParameter(3)/
                                   (1. + 2./static_cast<double>(i*i*i)),
                                   fitDistrDVal.GetParameter(3)*
                                   (1. + 2./static_cast<double>(i*i*i)));
         fitDistrDVal.SetParLimits(4, fitDistrDVal.GetParameter(4)*
                                   (1. - 6./static_cast<double>(i*i*i)),
                                   fitDistrDVal.GetParameter(4)*
                                   (1. + 4./static_cast<double>(i*i*i)));
         fitDistrDVal.SetParLimits(5, fitDistrDVal.GetParameter(5)/
                                   (1. + 5./static_cast<double>(i*i*i)),
                                   fitDistrDVal.GetParameter(5)*
                                   (1. + 5./static_cast<double>(i*i*i)));

         fitDistrDVal.SetRange(fitDistrDVal.GetParameter(1) - 5.*fitDistrDVal.GetParameter(2), 
                               fitDistrDVal.GetParameter(1) + 5.*fitDistrDVal.GetParameter(2));
      }

      distrDValProj->Fit(&fitDistrDVal, "RQMB");
      
      grMeansDVal.AddPoint(distrDValVsPT->GetYaxis()->GetBinCenter(i), 
                           fitDistrDVal.GetParameter(1));
      grSigmasDVal.AddPoint(distrDValVsPT->GetYaxis()->GetBinCenter(i), 
                            fitDistrDVal.GetParameter(2));
   }

   if (grMeansDVal.GetN() == 0) 
   {
      CppTools::PrintError("Number of points in a graph for " + detectorName + 
                           ", " + variableName + ", " + chargeName + " equals 0");
   }

   TCanvas meansVsPTCanv("means vs pT canvas", "", 800, 800);

   TH1F fitParVsPTFrame("means and sigmas frame", "", 10, 0., 10.);

   fitParVsPTFrame.SetMinimum(CppTools::Minimum(TMath::MinElement(grMeansDVal.GetN(), 
                                                                  grMeansDVal.GetY())));
   if (fitParVsPTFrame.GetMinimum() < 0.) 
   {
      fitParVsPTFrame.SetMinimum(fitParVsPTFrame.GetMinimum()*1.2);
   }
   else
   {
      fitParVsPTFrame.SetMinimum(fitParVsPTFrame.GetMinimum()/1.2);
   }

   fitParVsPTFrame.SetMaximum(CppTools::Maximum(TMath::MaxElement(grMeansDVal.GetN(), 
                                                                  grMeansDVal.GetY())));
   if (fitParVsPTFrame.GetMaximum() < 0.) 
   {
      fitParVsPTFrame.SetMaximum(fitParVsPTFrame.GetMaximum()/1.2);
   }
   else
   {
      fitParVsPTFrame.SetMaximum(fitParVsPTFrame.GetMaximum()*1.2);
   }

   fitParVsPTFrame.Draw("AXIS");
   fitParVsPTFrame.Draw("SAME AXIS X+ Y+");

   grMeansDVal.SetMarkerStyle(20);
   grMeansDVal.SetMarkerColor(kBlack);
   grMeansDVal.SetMarkerSize(1.);

   grMeansDVal.Clone()->Draw("P");

   ROOTTools::PrintCanvas(&meansVsPTCanv, outputDir + "means_" + variableName + "_" + 
                          detectorName + "_" + chargeNameShort);

   TCanvas sigmasVsPTCanv("means vs pT canvas", "", 800, 800);

   fitParVsPTFrame.SetMinimum(CppTools::Minimum(TMath::MinElement(grSigmasDVal.GetN(), 
                                                                  grSigmasDVal.GetY())));
   if (fitParVsPTFrame.GetMinimum() < 0.) 
   {
      fitParVsPTFrame.SetMinimum(fitParVsPTFrame.GetMinimum()*1.5);
   }
   else
   {
      fitParVsPTFrame.SetMinimum(fitParVsPTFrame.GetMinimum()/1.5);
   }

   fitParVsPTFrame.SetMaximum(CppTools::Maximum(TMath::MaxElement(grSigmasDVal.GetN(), 
                                                                  grSigmasDVal.GetY())));
   if (fitParVsPTFrame.GetMaximum() < 0.) 
   {
      fitParVsPTFrame.SetMaximum(fitParVsPTFrame.GetMaximum()*1.5);
   }
   else
   {
      fitParVsPTFrame.SetMaximum(fitParVsPTFrame.GetMaximum()/1.5);
   }

   grSigmasDVal.SetMarkerStyle(20);
   grSigmasDVal.SetMarkerColor(kBlack);
   grSigmasDVal.SetMarkerSize(1);

   fitParVsPTFrame.Draw("AXIS");
   fitParVsPTFrame.Draw("SAME AXIS X+ Y+");

   grSigmasDVal.Clone()->Draw("P");

   ROOTTools::PrintCanvas(&sigmasVsPTCanv, outputDir + "sigmas_" +variableName + "_" + 
                          detectorName + "_" + chargeNameShort);

   numberOfCalls++;
}

#endif /* SIGMALIZED_RESIDUALS_CPP */
