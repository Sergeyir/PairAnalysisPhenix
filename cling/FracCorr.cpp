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

   double centralityMax = 100.;
   if (runName == "Run14HeAu200") centralityMax = 40.;
   else if (runName == "Run15pAl200") centralityMax = 30.;
   else if (runName == "Run15pAu200") centralityMax = 40.;
   else CppTools::PrintError("System non known for maximum centrality parameter");

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

   TGraphErrors grContamination;
   TGraphErrors grFullContamination;

   double maxContamination = 0.;

   for (int i = 1; i <= distr->GetXaxis()->GetNbins(); i++)
   {
      TH1D *distrProj = distr->ProjectionY(std::to_string(i).c_str(), i, i);

      const double fullIntegral = distrProj->Integral(1, distrProj->GetXaxis()->GetNbins());
      if (fullIntegral < 1.) continue;

      TF1 fitBG1("bg frac fit 1", "gaus");
      TF1 fitBG2("bg frac fit 2", "expo(0) + pol1(2)");

      const double fracThreshold = fit.Eval(distr->GetXaxis()->GetBinCenter(i));
      //if (fracThreshold < 0.3) continue;

      fitBG1.SetRange(fracThreshold/2., fracThreshold - 0.1);
      fitBG2.SetRange(fracThreshold - 0.2, fracThreshold - 0.1);

      fitBG1.SetParameters(100000., 1., 2.);
      fitBG1.SetParLimits(0., 1., distrProj->GetMaximum());
      fitBG1.SetParLimits(1, 0., 1.);
      fitBG1.SetParLimits(2, 0., 1.);

      distrProj->Fit(&fitBG1, "RQMBN");
      distrProj->Fit(&fitBG2, "RQMBN");

      // integral of a signal within frac > fracThreshold
      double integralSignal = 0.;
      // background integrals within frac > fracThreshold obtanied from approximations
      double integralBG1 = 0.;
      double integralBG2 = 0.;

      for (int j = distrProj->GetXaxis()->FindBin(fracThreshold); 
           j <= distrProj->GetXaxis()->GetNbins(); j++)
      {
         integralSignal += distrProj->GetBinContent(j);
         integralBG1 += fitBG1.Eval(distrProj->GetXaxis()->GetBinCenter(j));
         integralBG2 += fitBG2.Eval(distrProj->GetXaxis()->GetBinCenter(j));
      }

      // background integral within the whole frac range
      double fullIntegralBG = integralBG1;

      for (int j = 1; 
           j < distrProj->GetXaxis()->FindBin(fracThreshold); j++)
      {
         fullIntegralBG += distrProj->GetBinContent(j);
      }

      // 1st estimate of contamination rate within frac > fracThreshold
      const double contamination1 = integralBG1/integralSignal;
      // 2nd estimate of contamination rate within frac > fracThreshold
      const double contamination2 = integralBG2/integralSignal;
      // 1st estimate of a contamination rate within the whole frac distribution
      const double fullContamination1 = fullIntegralBG/fullIntegral;
      // 2nd estimate of a contamination rate within the whole frac distribution
      const double fullContamination2 = (fullIntegralBG - integralBG1)/fullIntegral;

      if (i <= distr->GetXaxis()->FindBin(centralityMax - 1e-6)) 
      {
         maxContamination = 
            CppTools::Maximum(maxContamination, contamination1, contamination2, fullContamination1, 
                              fullContamination2);

         grContamination.AddPoint(distr->GetXaxis()->GetBinCenter(i), 
                                  (contamination1 + contamination2)/2.);

         grContamination.
            SetPointError(grContamination.GetN() - 1, distr->GetXaxis()->GetBinWidth(i)/2., 
                          fabs(contamination1 - contamination2)/2.);

         grFullContamination.AddPoint(distr->GetXaxis()->GetBinCenter(i), 
                                      (fullContamination1 + fullContamination2)/2.);
         grFullContamination.
            SetPointError(grFullContamination.GetN() - 1, distr->GetXaxis()->GetBinWidth(i)/2., 
                          fabs(fullContamination1 - fullContamination2)/2.);

         CppTools::Print(i, fracThreshold, contamination1, contamination2, 
                         fullContamination1, fullContamination2, integralSignal/fullIntegral);
      }

      //fitBG1.SetRange(0., 1.);
      fitBG1.SetRange(fracThreshold/2., 1.);
      fitBG2.SetRange(fracThreshold - 0.2, 1.);

      distrProj->SetLineWidth(2);
      distrProj->SetLineColor(kBlack);

      fitBG1.SetLineWidth(2);
      fitBG1.SetLineColorAlpha(kRed, 0.8);

      fitBG2.SetLineWidth(2);
      fitBG2.SetLineStyle(2);
      fitBG2.SetLineColorAlpha(kAzure, 0.8);

      TCanvas canv("canv", "canv", 800, 800);

      gPad->SetRightMargin(0.01); gPad->SetTopMargin(0.02); 
      gPad->SetLeftMargin(0.095); gPad->SetBottomMargin(0.08);

      gPad->SetLogy();

      ROOTTools::DrawFrame(distrProj, "", "frac", "Counts", 0.9, 1.2, 0.04, 0.04);
      //ROOTTools::DrawFrame(distrProj, );
      fitBG1.Draw("SAME");
      fitBG2.Draw("SAME");

      ROOTTools::DrawLine(fracThreshold, 0., fracThreshold, distrProj->GetMaximum(), 
                          kGray + 1, 0.8, 2, 2);

      ROOTTools::PrintCanvas(&canv, outputDir + "/" + std::to_string(i));
   }

   grContamination.SetLineColor(kRed);
   grContamination.SetMarkerStyle(0);
   grContamination.SetMarkerColorAlpha(0, 0.);
   grContamination.SetFillStyle(1001);
   grContamination.SetFillColorAlpha(kRed - 3, 0.3);

   grFullContamination.SetLineColor(kBlack);
   grFullContamination.SetMarkerStyle(0);
   grFullContamination.SetMarkerColorAlpha(0, 0.);
   grFullContamination.SetFillStyle(1001);
   grFullContamination.SetFillColorAlpha(kBlack, 0.2);

   TF1 fitContamination("contamination rate fit", "pol0");
   
   fitContamination.SetRange(0.01, centralityMax - 0.1);
   fitContamination.SetLineStyle(2);
   fitContamination.SetLineWidth(4);
   fitContamination.SetLineColor(kRed + 1);

   grContamination.Fit(&fitContamination, "RQMN");

   CppTools::Print(fitContamination.GetParameter(0));

   TCanvas canv("canv", "canv", 800, 800);

   canv.Clear();

   gPad->SetRightMargin(0.01); gPad->SetTopMargin(0.02); 
   gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.08);

   ROOTTools::DrawFrame(0.01, 0., centralityMax - 0.1, maxContamination*1.1, "", "centrality", 
                        "Contamination rate", 0.9, 1.9, 0.04, 0.04);

   fitContamination.Draw("SAME");

   grFullContamination.Draw("5");
   grContamination.Draw("5");

   ROOTTools::PrintCanvas(&canv, outputDir + "/contamination_rate");

   exit(1);
}
