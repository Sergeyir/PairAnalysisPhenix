#pragma once

#include "IOTools.hpp"
#include "MathTools.hpp"

#include "TCanvasTools.hpp"

void FracCorr()
{
   gROOT->SetBatch(kTRUE);
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);

   gErrorIgnoreLevel = kWarning;

   CppTools::PrintInfo("List of runs in data/Real directory");
   void(system("ls data/Real/"));

   CppTools::Print("Choose the run from the above and type it in");
   std::string runName;
   std::cout << ">> ";
   std::cin >> runName;

   const std::string inputFileName = "data/Real/" + runName + "/SingleTrack/sum.root";
   CppTools::CheckInputFile(inputFileName);

   TFile inputFile(inputFileName.c_str());

   std::string distrName; 
   distrName = "centrality vs frac";

   TH2F *distr = static_cast<TH2F *>(inputFile.Get(distrName.c_str()));

   if (!distr)
   {
      CppTools::PrintError("No histogram named " + distrName + 
                           " in file " + inputFileName);
   }

   distr->SetTitle("");

   const std::string fitFunc = 
      "([0] - exp((x/[1])^4) - x/[2])*(([0] - exp((x/[1])^4) - x/[2]) > 0.5) +"\
      "0.5*(([0] - exp((x/[1])^4) - x/[2]) <= 0.5)";


   TF1 fit("frac fit", fitFunc.c_str());

   const std::string parFileName = "data/Parameters/FracCut/" + runName + ".txt";

   CppTools::CheckInputFile(parFileName);
   std::ifstream parFile(parFileName);

   double par;
   int i = 0;
   while (parFile >> par)
   {
      fit.SetParameter(i, par);
      i++;
   }

   const std::string outputDir = "output/Frac/" + runName;
   std::filesystem::create_directories(outputDir);

   TGraph grContamination;
   TGraphErrors grFullContamination;

   double maxContamination = 0.;

   for (int i = 1; i <= distr->GetXaxis()->FindBin(39.); i++)
   {
      TH1D *distrProj = distr->ProjectionY(std::to_string(i).c_str(), i, i);

      const double fullIntegral = distrProj->Integral(1, distrProj->GetXaxis()->GetNbins());
      if (fullIntegral < 1.) continue;

      TF1 fitBG("bg frac fit", "gaus(0)");

      const double fracThreshold = fit.Eval(distr->GetXaxis()->GetBinCenter(i));
      //if (fracThreshold < 0.3) continue;

      fitBG.SetRange(fracThreshold/2., fracThreshold - 0.1);

      fitBG.SetParameters(100000., 1., 2.);
      fitBG.SetParLimits(0., 1., distrProj->GetMaximum());
      fitBG.SetParLimits(1, 0., 1.);
      fitBG.SetParLimits(2, 0., 1.);

      distrProj->Fit(&fitBG, "RQMBN");

      // integral of a signal within frac > fracThreshold
      double integralSignal = 0.;
      // background integral within frac > fracThreshold
      double integralBG = 0.;

      for (int j = distrProj->GetXaxis()->FindBin(fracThreshold); 
           j <= distrProj->GetXaxis()->GetNbins(); j++)
      {
         integralSignal += distrProj->GetBinContent(j);
         integralBG += fitBG.Eval(distrProj->GetXaxis()->GetBinCenter(j));
      }

      // background integral within the whole frac range
      double fullIntegralBG = integralBG;

      for (int j = 1; 
           j < distrProj->GetXaxis()->FindBin(fracThreshold); j++)
      {
         fullIntegralBG += distrProj->GetBinContent(j);
      }

      // contamination rate within frac > fracThreshold
      const double contamination = integralBG/integralSignal;
      // 1st estimate of a contamination rate within the whole frac distribution
      const double fullContamination1 = fullIntegralBG/fullIntegral;
      // 2nd estimate of a contamination rate within the whole frac distribution
      const double fullContamination2 = (fullIntegralBG - integralBG)/fullIntegral;

      maxContamination = 
         CppTools::Maximum(maxContamination, contamination, fullContamination1, 
                           fullContamination2);

      grContamination.AddPoint(distr->GetXaxis()->GetBinCenter(i), contamination);

      grFullContamination.AddPoint(distr->GetXaxis()->GetBinCenter(i), 
                                   (fullContamination1 + fullContamination2)/2.);
      grFullContamination.
         SetPointError(grFullContamination.GetN() - 1, distr->GetXaxis()->GetBinWidth(i)/2., 
                       fabs(fullContamination1 - fullContamination2)/2.);

      CppTools::Print(i, fracThreshold, contamination, fullContamination1, fullContamination2);

      //fitBG.SetRange(0., 1.);
      fitBG.SetRange(fracThreshold/2., 1.);

      distrProj->SetLineWidth(2);
      distrProj->SetLineColor(kBlack);

      fitBG.SetLineWidth(2);
      fitBG.SetLineColorAlpha(kRed, 0.8);

      TCanvas canv("canv", "canv", 800, 800);

      gPad->SetRightMargin(0.01); gPad->SetTopMargin(0.02); 
      gPad->SetLeftMargin(0.095); gPad->SetBottomMargin(0.08);

      gPad->SetLogy();

      ROOTTools::DrawFrame(distrProj, "", "frac", "Counts", 0.9, 1.2, 0.04, 0.04);
      //ROOTTools::DrawFrame(distrProj, );
      fitBG.Draw("SAME");

      ROOTTools::DrawLine(fracThreshold, 0., fracThreshold, distrProj->GetMaximum(), 
                          kGray + 1, 0.8, 2, 2);

      ROOTTools::PrintCanvas(&canv, outputDir + "/" + std::to_string(i));
   }

   grContamination.SetMarkerStyle(20);
   grContamination.SetMarkerSize(1.4);
   grContamination.SetMarkerColor(kRed + 1);

   grFullContamination.SetLineColor(kBlack);
   grFullContamination.SetMarkerStyle(0);
   grFullContamination.SetMarkerColorAlpha(0, 0.);
   grFullContamination.SetFillStyle(1001);
   grFullContamination.SetFillColorAlpha(kBlack, 0.2);

   TCanvas canv("canv", "canv", 800, 800);

   canv.Clear();

   gPad->SetRightMargin(0.01); gPad->SetTopMargin(0.02); 
   gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.08);

   ROOTTools::DrawFrame(0., 0.01, 39.9, maxContamination*1.1, "", "centrailty", 
                        "Contamination rate", 0.9, 1.9, 0.04, 0.04);


   grFullContamination.Draw("5");
   grContamination.Draw("P");

   ROOTTools::PrintCanvas(&canv, outputDir + "/contamination_rate");

   exit(1);
}
