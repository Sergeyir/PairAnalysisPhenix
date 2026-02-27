#pragma once

#include "IOTools.hpp"

void ConstantTimingOffset()
{
   gROOT->SetBatch(kTRUE);
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);

   CppTools::PrintInfo("Specify the input file:");
   std::string inputFileName;
   std::cout << ">> ";
   std::cin >> inputFileName;

   CppTools::CheckInputFile(inputFileName);
   TFile inputFile(inputFileName.c_str());

   std::array<std::string, 10> detectorNames = {"EMCale0", "EMCale1", "EMCale2", "EMCale3",
                                                "EMCalw0", "EMCalw1", "EMCalw2", "EMCalw3",
                                                "TOFe", "TOFw"};

   for (const std::string& detectorName : detectorNames)
   {
      const std::string distrTimeVsPTName = "t - t_exp^pi, " + detectorName;
      TH2F *distrTimeVsPT = static_cast<TH2F *>(inputFile.Get(distrTimeVsPTName.c_str()));

      if (!distrTimeVsPT) 
      {
         CppTools::PrintError("No histogram named \"" + distrTimeVsPTName + 
                              "\" in file: " + inputFileName);
      }

      TH1D *distrTime = distrTimeVsPT->ProjectionX();

      TF1 fit("gaus", "gaus");

      const double maxBinVal = distrTime->GetBinContent(distrTime->GetMaximumBin());
      const double binWidth = distrTime->GetXaxis()->GetBinWidth(1);
      const double maxBinX = distrTime->GetXaxis()->GetBinUpEdge(distrTime->GetMaximumBin());

      fit.SetRange(maxBinX - binWidth*5., maxBinX + binWidth*5.);
      fit.SetParameters(maxBinVal, maxBinX, 0.2);
      fit.SetParLimits(1, maxBinX - binWidth*2., maxBinX + binWidth*2.);

      distrTime->Fit(&fit, "RQM");

      CppTools::Print(detectorName, -fit.GetParameter(1));
   }
}
