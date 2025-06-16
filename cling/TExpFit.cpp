#pragma once

#include "IOTools.hpp"

// not finished yet

void TExpFit()
{
   CppTools::PrintInfo("Available run names with real data:");
   system("ls data/Real");
   CppTools::Print("Enter the run name");
   std::string runName;
   std::cout << ">> ";
   std::cin >> runName;

   const std::string inputFileName = "data/Real/" + runName + "/SingleTrack/sum.root";
   CppTools::CheckInputFile(inputFileName);
   TFile inputFile(inputFileName.c_str());

   CppTools::PrintInfo("Availble t - t_exp distributions:");
   inputFile.ls("t - t_exp^pi,*");

   CppTools::Print("Enter the detector name");
   std::string detectorName;
   std::cout << ">> ";
   std::cin >> detectorName;

   TH2F *distrT = static_cast<TH2F *>(inputFile.Get(("t - t_exp^pi, " + detectorName).c_str()));

   if (!distrT) 
   {
      CppTools::PrintError("No histogram named \"t - t_exp^pi, " + 
                           detectorName + "\" in file: " + inputFileName);
   }

   TGraph means;

   for (int i = 1; i <= distrT->GetXaxis()->GetNbins(); i++)
   {
      TH1D *distrTProj distrT->ProjectionX((distrT->GetName + std::to_string(i)).c_str(), i, i);

      TF1 gausFit("gaus fit", "gaus");
   }
}
