#pragma once

#include "IOTools.hpp"

void TOFwADC()
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

   TH2F *distrTOFwADC = static_cast<TH2F *>(inputFile.Get("ADC: TOFw"));
   if (!distrTOFwADC) 
   {
      CppTools::PrintError("No histogram named \"ADC: TOFw\" in file: " + inputFileName);
   }

   const double fullIntegral = 
      distrTOFwADC->Integral(1, distrTOFwADC->GetXaxis()->GetNbins(),
                             1, distrTOFwADC->GetYaxis()->GetNbins());

   const double cutADCIntegral = 
      distrTOFwADC->Integral(1, distrTOFwADC->GetXaxis()->GetNbins(),
                             distrTOFwADC->GetYaxis()->FindBin(61.), 
                             distrTOFwADC->GetYaxis()->FindBin(599));

   CppTools::PrintInfo("TOFw ADC correction for simulation is:");
   CppTools::Print(cutADCIntegral/fullIntegral);
}
