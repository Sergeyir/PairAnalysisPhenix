#pragma once

#include "IOTools.hpp"

void RunByRunTimingOffset()
{
   gROOT->SetBatch(kTRUE);
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);

   /*
   CppTools::PrintInfo("Specify the input directory:");
   std::string inputDir;
   std::cout << ">> ";
   std::cin >> inputDir;
   */
   const std::string inputDir = "data/Real/Run15pp200/SingleTrack";


   std::array<std::string, 1> detectorNames = {/*"EMCale0", "EMCale1", "EMCale2", "EMCale3",
                                                "EMCalw0", "EMCalw1", "EMCalw2", "EMCalw3",
                                                "TOFe", */"TOFw"};

   TH1D distrT("t", "", 1000, 0., 1000.);

   int i = 1;
   for (const auto &file : std::filesystem::directory_iterator(inputDir))
   {
      const std::string fileName = static_cast<std::string>(file.path());
      std::cout << fileName << " ";

      CppTools::CheckInputFile(fileName);
      TFile inputFile(fileName.c_str());

      for (const std::string& detectorName : detectorNames)
      {
         const std::string distrTimeVsPTName = "t - t_exp^pi, " + detectorName;
         TH2F *distrTimeVsPT = static_cast<TH2F *>(inputFile.Get(distrTimeVsPTName.c_str()));

         if (!distrTimeVsPT) 
         {
            CppTools::PrintError("No histogram named \"" + distrTimeVsPTName + 
                                 "\" in file: " + fileName);
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

         std::cout << -fit.GetParameter(1) << " ";
         distrT.SetBinContent(i, -fit.GetParameter(1));
      }
      std::cout << std::endl;
      i++;
   }

   TCanvas canv("canv", "", 800, 800);
   distrT.Draw();
   canv.SaveAs("tmp/t.png");
}
