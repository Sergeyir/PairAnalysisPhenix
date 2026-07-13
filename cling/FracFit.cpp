#pragma once

#include "IOTools.hpp"
#include "InputYAMLReader.hpp"

#include "TCanvasTools.hpp"

void FracFit()
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

   const std::string fitFunc = "[0] - [1]*exp([2]*x*x)";

   CppTools::Print("Setting fit function \"" + fitFunc + "\"");
   TF1 fit("frac fit", fitFunc.c_str());

   fit.SetRange(distr->GetXaxis()->GetBinLowEdge(1), 
                distr->GetXaxis()->GetBinUpEdge(distr->GetXaxis()->GetNbins()));

   bool firstIteration = true;

   while (true)
   {
      if (firstIteration)
      {
         distr->Fit(&fit, "QMN");
         std::filesystem::create_directories("output/FracFit/");
      }
      else
      {
         CppTools::Print("Current parameters are:", fit.GetParameter(0), fit.GetParameter(1), fit.GetParameter(2));
         CppTools::Print("Type in fit parameters (divide them by spaces). To exit type any character(s)");
         double par1, par2, par3;

         if (!(std::cin >> par1 >> par2 >> par3))
         {
            CppTools::PrintInfo("Exiting the program");
            exit(0);
         }

         fit.SetParameters(par1, par2, par3);
      }

      TCanvas canv("c", "", 800, 800);
      gPad->SetLogz(true);

      distr->Draw("COLZ");
      fit.Draw("SAME");

      ROOTTools::PrintCanvas(&canv, "output/FracFit/" + runName);

      if (firstIteration)
      {
         void (system(("xdg-open output/FracFit/" + runName + ".png").c_str()));
      }
      firstIteration = false;
   }

   exit(1);
}
