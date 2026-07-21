#pragma once

#include "IOTools.hpp"
#include "InputYAMLReader.hpp"

#include "TCanvasTools.hpp"

void FracCorr()
{
   gROOT->SetBatch(kTRUE);
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);

   CppTools::PrintInfo("List of runs in data/Real directory");
   void(system("ls data/Real/"));

   CppTools::Print("Choose the run from the above and type it in");
   std::string runName;
   std::cout << ">> ";
   std::cin >> runName;

   const std::string inputFileName = "data/Real/" + runName + "/SingleTrack/sum.root";
   CppTools::CheckInputFile(inputFileName);

   
   InputYAMLReader inputYAMLMain("input/" + runName + "/main.yaml", "main");

   TFile inputFile(inputFileName.c_str());

   std::string distrName; 
   if (inputYAMLMain["is_pp"].as<bool>()) distrName = "multiplicity vs frac";
   else distrName = "frac vs centrality";

   TH2F *distr = static_cast<TH2F *>(inputFile.Get(distrName.c_str()));

   if (!distr)
   {
      CppTools::PrintError("No histogram named " + distrName + 
                           " in file " + inputFileName);
   }

   distr->SetTitle("");

   const std::string fitFunc = "[0] - [1]*exp([2]*x*x)";

   CppTools::Print("Setting fit function \"" + fitFunc + "\"");
   TF1 fit("frac fit", fitFunc.c_str());

   const std::string parFileName = "data/Parameters/FracCut/" + runName + ".txt";

   CppTools::CheckInputFile(parFileName);
   std::ifstream parFile(parFileName);

   double par;
   int i = 0;
   while (parFile >> par)
   {
      CppTools::Print(par);
      fit.SetParameter(i, par);
      i++;
   }

   const std::string outputDir = "output/FracCorr/" + runName;
   std::filesystem::create_directories(outputDir);

   for (int i = 1; i <= distr->GetYaxis()->GetNbins(); i++)
   {
      TH1D *distrProj = distr->ProjectionX(std::to_string(i).c_str(), i, i);

      TF1 fitBG("bg frac fit", "gaus(0)");

      const double fracThreshold = fit.Eval(distr->GetYaxis()->GetBinCenter(i));
      //if (fracThreshold < 0.3) continue;

      fitBG.SetRange(0.3, 0.9);

      distrProj->Fit(&fitBG, "RQM");

      double integral = 0.;
      double integralBG = 0.;

      for (int j = distrProj->GetXaxis()->FindBin(fracThreshold); 
           j <= distrProj->GetXaxis()->GetNbins(); j++)
      {
         integral += distrProj->GetBinContent(j);
         integralBG += fitBG.Eval(distrProj->GetYaxis()->GetBinCenter(j));
      }

      CppTools::Print(i, fracThreshold, integralBG/integral);

      fitBG.SetRange(0.3, 1.);

      TCanvas canv("canv", "canv", 600, 600);

      gPad->SetLogy();

      distrProj->Draw();
      fitBG.Draw("SAME");
      canv.SaveAs((outputDir + "/" + std::to_string(i) + ".png").c_str());
   }

   exit(1);
}
