#pragma once

#include "IOTools.hpp"
#include "StrTools.hpp"
#include "MathTools.hpp"

#include "HistTools.hpp"
#include "TCanvasTools.hpp"

#include "PBar.hpp"

void SlewingCorrectionTOFw()
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

   const std::string distrName = "t - t_exp^pi vs ADC, TOFw";
   TH2F *distr = static_cast<TH2F *>(inputFile.Get(distrName.c_str()));

   if (!distr)
   {
      CppTools::PrintError("No histogram named " + distrName + " in file " + inputFileName);
   }

   ROOTTools::SwapAxis(distr);

   const std::string outputDir = "output/SlewingTOFw/" + runName;
   std::filesystem::create_directories(outputDir);

   TGraph grMeans;
   TGraph grSigmas;

   double minY = 1e31;
   double maxY = -1e31;

   ProgressBar pBar("FANCY");

   const int xNBins = distr->GetXaxis()->GetNbins();
   for (int i = 1; i <= xNBins; i++)
   {
      pBar.Print(static_cast<double>(i - 1)/static_cast<double>(xNBins));

      TH1D *distrProj = distr->ProjectionY(std::to_string(i).c_str(), i, i);

      TF1 fit("fit", "gaus(0) + gaus(3)");
      TF1 fitFG("fit FG", "gaus(0)");
      TF1 fitBG("fit BG", "gaus(0)");

      fit.SetLineColor(kRed - 3);
      fitFG.SetLineColor(kAzure - 3);
      fitBG.SetLineColor(kGreen - 3);

      fitFG.SetLineStyle(2);
      fitBG.SetLineStyle(2);

      const int maxBin = distrProj->GetMaximumBin();

      fit.SetParameters(distrProj->GetBinContent(maxBin)*0.99,
                        distrProj->GetBinCenter(maxBin),
                        0.15,
                        distrProj->GetBinContent(maxBin)/20.,
                        distrProj->GetBinCenter(maxBin),
                        1.);

      fit.SetParLimits(0, distrProj->GetBinContent(maxBin)/2.,
                       distrProj->GetBinContent(maxBin));
      fit.SetParLimits(1, distrProj->GetBinCenter(maxBin - 30),
                       distrProj->GetBinCenter(maxBin + 30));
      fit.SetParLimits(2, 0.05, 0.3);
      fit.SetParLimits(3, 0., distrProj->GetBinContent(maxBin)/5.);
      fit.SetParLimits(4, -5., 5.);
      fit.SetParLimits(5, 0.3, 2.);

      fit.SetRange(distrProj->GetBinCenter(maxBin) - 2., distrProj->GetBinCenter(maxBin) + 2.);

      distrProj->Fit(&fit, "RQMNB");


      fit.SetRange(fit.GetParameter(1) - fit.GetParameter(2)*5., 
                   fit.GetParameter(1) + fit.GetParameter(2)*5.);

      distrProj->Fit(&fit, "RQMNB");

      fit.SetRange(fit.GetParameter(1) - fit.GetParameter(2)*5., 
                   fit.GetParameter(1) + fit.GetParameter(2)*5.);
      fitFG.SetRange(fit.GetParameter(1) - fit.GetParameter(2)*5., 
                     fit.GetParameter(1) + fit.GetParameter(2)*5.);
      fitBG.SetRange(fit.GetParameter(1) - fit.GetParameter(2)*5., 
                     fit.GetParameter(1) + fit.GetParameter(2)*5.);

      distrProj->Fit(&fit, "RQMNB");

      for (int k = 0; k < fitFG.GetNpar(); k++)
      {
         fitFG.SetParameter(k, fit.GetParameter(k));
      }
      for (int k = 0; k < fitBG.GetNpar(); k++)
      {
         fitBG.SetParameter(k, fit.GetParameter(k + 3));
      }

      grMeans.AddPoint(distr->GetXaxis()->GetBinCenter(i), fit.GetParameter(1));
      grSigmas.AddPoint(distr->GetXaxis()->GetBinCenter(i), fit.GetParameter(2));

      minY = CppTools::Minimum(minY, fit.GetParameter(1));
      maxY = CppTools::Maximum(maxY, fit.GetParameter(1));

      TCanvas canv("fit canv", "", 500, 500);

      gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.09); 
      gPad->SetLeftMargin(0.17); gPad->SetBottomMargin(0.112);

      ROOTTools::DrawFrame(distrProj, CppTools::DtoStr(distr->GetXaxis()->GetBinLowEdge(i), 0) + 
                           "<q_{TOFw}<" + CppTools::DtoStr(distr->GetXaxis()->GetBinUpEdge(i), 0), 
                           "t - t_{exp}^{#pi^{#pm}} [ns]", "Counts", 1., 1.9);

      fit.Draw("SAME");
      fitFG.Draw("SAME");
      fitBG.Draw("SAME");

      ROOTTools::PrintCanvas(&canv, outputDir + "/" + std::to_string(i));
   }

   pBar.Finish();

   if (minY < 0) minY *= 1.1;
   else minY /= 1.1;

   if (maxY > 0) maxY *= 1.1;
   else maxY /= 1.1;

   distr->GetYaxis()->SetRange(distr->GetYaxis()->FindBin(minY - 0.5),
                               distr->GetYaxis()->FindBin(maxY + 0.5));

   grMeans.SetMarkerStyle(20);
   grSigmas.SetMarkerStyle(20);

   grMeans.SetMarkerSize(1.);
   grSigmas.SetMarkerSize(1.);

   grMeans.SetMarkerColorAlpha(kRed, 0.9);
   grSigmas.SetMarkerColorAlpha(kRed, 0.9);

   TF1 fitMeans("means vs adc fit", "[0] + ([1]/x^0.4)");

   fitMeans.SetRange(60., 600.);
   fitMeans.SetLineStyle(2);
   fitMeans.SetLineWidth(3);
   fitMeans.SetLineColor(kGray + 3);

   grMeans.Fit(&fitMeans, "RQMN");

   TCanvas canv("par canv", "", 800, 800);

   gPad->SetRightMargin(0.01); gPad->SetTopMargin(0.02); 
   gPad->SetLeftMargin(0.17); gPad->SetBottomMargin(0.112);

   ROOTTools::DrawFrame(60, minY, 600, maxY, "", "q_{TOFw}",
                        "#mu_{t - t_{exp}^{#pi^{#pm}}} [ns]", 1., 1.6);

   fitMeans.Draw("SAME");
   grMeans.Draw("P");

   ROOTTools::PrintCanvas(&canv, outputDir + "/means");

   canv.Clear();

   gPad->SetRightMargin(0.01); gPad->SetTopMargin(0.02); 
   gPad->SetLeftMargin(0.17); gPad->SetBottomMargin(0.112);

   ROOTTools::DrawFrame(60, 0., 600, 0.3, "", "q_{TOFw}", 
                        "#sigma_{t - t_{exp}^{#pi^{#pm}}} [ns]", 1., 1.6);

   grSigmas.Draw("P");

   ROOTTools::PrintCanvas(&canv, outputDir + "/sigmas");

   canv.Clear();

   gPad->SetRightMargin(0.11); gPad->SetTopMargin(0.02); 
   gPad->SetLeftMargin(0.09); gPad->SetBottomMargin(0.09);

   gPad->SetLogz();

   ROOTTools::DrawFrame(distr, "", "q_{TOFw}", "t - t_{exp}^{#pi^{#pm}} [ns]", 
                        1.0, 1.0, 0.04, 0.04, true, true, "COLZ");

   fitMeans.SetLineColor(kRed - 3);
   fitMeans.SetLineWidth(4);
   fitMeans.Draw("SAME");

   ROOTTools::PrintCanvas(&canv, outputDir + "/distr_with_fit");

   CppTools::Print("Fit parameters:", fitMeans.GetParameter(0), fitMeans.GetParameter(1));

   exit(1);
}
