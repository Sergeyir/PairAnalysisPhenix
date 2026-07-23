#pragma once

#include "IOTools.hpp"

#include "TCanvasTools.hpp"

void ElossFit()
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
   TFile inputFile(inputFileName.c_str());

   const std::string distrName = "beta vs E, TOFe";
   TH2F *distr = static_cast<TH2F *>(inputFile.Get(distrName.c_str()));

   if (!distr)
   {
      CppTools::PrintError("No histogram named " + distrName + 
                           " in file " + inputFileName);
   }

   const std::string fitFunc = "[0]*pow(x, [1])";

   CppTools::Print("Setting fit function \"" + fitFunc + "\"");
   TF1 fit("E_{loss} fit", fitFunc.c_str());

   bool firstIteration = true;

   while (true)
   {
      if (firstIteration)
      {
         if (std::filesystem::exists("data/Parameters/ElossCut/" + runName + ".txt"))
         {
            std::ifstream parFile("data/Parameters/ElossCut/" + runName + ".txt");

            double par;
            int i = 0;
            while (parFile >> par)
            {
               fit.SetParameter(i, par);
               i++;
            }
         }
         else
         {
            distr->Fit(&fit, "QMN");
            std::filesystem::create_directories("output/ElossFit/");
         }
      }
      else
      {
         CppTools::Print("Current parameters are:", fit.GetParameter(0), fit.GetParameter(1));
         CppTools::Print("Type in fit parameters (divide them by spaces). "\
                         "To exit type any character(s)");
         double par1, par2;

         if (!(std::cin >> par1 >> par2))
         {
            CppTools::PrintInfo("Exiting the program");
            exit(0);
         }

         fit.SetParameters(par1, par2);
      }

      TCanvas canv("c", "", 800, 800);
      gPad->SetLogz(true);

      gPad->SetRightMargin(0.11); gPad->SetTopMargin(0.02); 
      gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.08);

      ROOTTools::DrawFrame(distr, "", "#it{#beta}", "#it{E}_{loss}", 0.9, 1.9, 0.04, 0.04, 
                           true, true, "COLZ");

      fit.Draw("SAME");

      ROOTTools::PrintCanvas(&canv, "output/ElossFit/" + runName);

      if (firstIteration)
      {
         void (system(("xdg-open output/ElossFit/" + runName + ".png").c_str()));
      }
      firstIteration = false;
   }

   exit(1);
}
