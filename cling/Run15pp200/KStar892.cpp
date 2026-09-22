#pragma once

#include "IOTools.hpp"

#include "TCanvasTools.hpp"
#include "MathTools.hpp"

#include "PainterHelper.hpp"
#include "InputYAMLReader.hpp"

const std::string runName = "Run15pp200";
const std::string resonanceName = "KStar892";
const int taxiNumber = 20488;

TFile *resultsInputFile;

void KStar892()
{
   gROOT->SetBatch(kTRUE);
   gStyle->SetOptStat(kFALSE);

   gErrorIgnoreLevel = kWarning;
   gStyle->SetOptStat(0);

   InputYAMLReader inputYAMLResonance("input/" + runName + "/" + resonanceName + ".yaml");
   inputYAMLResonance.CheckStatus("resonance");

   InputYAMLReader inputYAMLMain("input/" + runName + "/main.yaml");
   inputYAMLMain.CheckStatus("main");

   const std::string resultsInputFileName = "data/Results/" + runName + "/" + 
                                            std::to_string(taxiNumber) + "_" + 
                                            resonanceName + ".root";

   CppTools::CheckInputFile(resultsInputFileName);

   resultsInputFile = TFile::Open(resultsInputFileName.c_str());

   const double resonanceMass = inputYAMLResonance["mass"].as<double>();

   const std::string outputDir = "output/Results/" + runName + "/" + std::to_string(taxiNumber);
   std::filesystem::create_directories(outputDir);

   TLegend legend(0.7, 0.8, 0.95, 0.95);
   legend.SetLineColorAlpha(0, 0.);
   legend.SetFillColorAlpha(0, 0.);

   // Spectra (old and new)
   {
      PainterHelper spectra(&legend);
      spectra.SetMarkerSize(1.4);
      spectra.SetLineWidth(2);
      spectra.SetDefaultSysWidth(0.1);

      const std::string oldSpectraFileName = "data/Spectra/pp200/KStar892_averaged.root";
      CppTools::CheckInputFile(oldSpectraFileName);

      TFile *oldSpectraFile = TFile::Open(oldSpectraFileName.c_str());

      TH1D *distrSpectraVsPTStatErr = static_cast<TH1D *>
         (resultsInputFile->Get("MB/spectra vs pT with stat errors"));
      TH1D *distrSpectraVsPTSysErr = static_cast<TH1D *>
         (resultsInputFile->Get("MB/spectra vs pT with sys errors"));

      TH1D *distrOldSpectraVsPTStatErr = static_cast<TH1D *>
         (oldSpectraFile->Get("spectra vs pT with stat err"));
      TH1D *distrOldSpectraVsPTSysErr = static_cast<TH1D *>
         (oldSpectraFile->Get("spectra vs pT with sys err"));

      distrOldSpectraVsPTStatErr->Scale(1./42.2);
      distrOldSpectraVsPTSysErr->Scale(1./42.2);

      const double xMin = CppTools::Minimum(distrOldSpectraVsPTStatErr->GetBinLowEdge(1),
                                            distrSpectraVsPTStatErr->GetBinLowEdge(1));
      const double xMax = 
         CppTools::Maximum(distrOldSpectraVsPTStatErr->GetXaxis()->
                           GetBinUpEdge(distrOldSpectraVsPTStatErr->GetXaxis()->GetNbins()),
                           distrSpectraVsPTStatErr->GetXaxis()->
                           GetBinUpEdge(distrSpectraVsPTStatErr->GetXaxis()->GetNbins()));

      double yMin = CppTools::Minimum(distrOldSpectraVsPTStatErr->GetMinimum(), 
                                      distrSpectraVsPTStatErr->GetMinimum());
      double yMax = CppTools::Maximum(distrOldSpectraVsPTStatErr->GetMaximum(), 
                                      distrSpectraVsPTStatErr->GetMaximum());

      TF1 *tsallisFit = static_cast<TF1 *>(resultsInputFile->Get("MB/tsallis fit"));

      TCanvas canv("canv", "canv", 800, 800);

      canv.SetFillStyle(0);
      canv.SetFrameFillColor(0);
      canv.SetFrameFillStyle(0);
      canv.SetFrameBorderMode(0);

      gPad->SetRightMargin(0.002); gPad->SetTopMargin(0.002); 
      gPad->SetLeftMargin(0.204); gPad->SetBottomMargin(0.105);

      gPad->SetLogy();

      ROOTTools::DrawFrame(xMin - 0.1, yMin/5., xMax + 0.1, yMax*5., "", "#it{p}_{T} [GeV/#it{c}]",
                           "#frac{1}{2#pi#it{p}_{T}} #frac{#it{d}^{2}"\
                           "#it{N}}{#it{dp}_{T}#it{dy}} [GeV/#it{c}]^{-2}", 0.91, 1.85);

      tsallisFit->Draw("SAME");

      spectra.DrawHistogram(distrOldSpectraVsPTStatErr, distrOldSpectraVsPTSysErr, kBlack,
                            0.9, 21, "old");
      spectra.DrawHistogram(distrSpectraVsPTStatErr, distrSpectraVsPTSysErr, kRed - 3,
                            0.9, 20, "new");

      legend.Draw();

      ROOTTools::PrintCanvas(&canv, outputDir + "/" + resonanceName + "_spectra");

      canv.Clear();
      legend.Clear();

      gPad->SetLogy(false);
      gPad->SetRightMargin(0.002); gPad->SetTopMargin(0.002); 
      gPad->SetLeftMargin(0.115); gPad->SetBottomMargin(0.112);

      distrOldSpectraVsPTStatErr->Divide(tsallisFit);
      distrOldSpectraVsPTSysErr->Divide(tsallisFit);
      distrSpectraVsPTStatErr->Divide(tsallisFit);
      distrSpectraVsPTSysErr->Divide(tsallisFit);

      yMin = CppTools::Minimum(distrOldSpectraVsPTStatErr->GetMinimum(), 
                               distrSpectraVsPTStatErr->GetMinimum());
      yMax = CppTools::Maximum(distrOldSpectraVsPTStatErr->GetMaximum(), 
                               distrSpectraVsPTStatErr->GetMaximum());

      ROOTTools::DrawFrame(xMin - 0.1, yMin/1.3, xMax + 0.1, yMax*1.3, 
                           "", "#it{p}_{T} [GeV/#it{c}]", "Data/Fit", 1., 1.2);

      if (yMin/1.1 < 1. && yMax*1.1 > 1.)
      {
         TLine line(xMin - 0.1, 1., xMax + 0.1, 1.);
         line.SetLineColorAlpha(kBlack, 0.5);
         line.SetLineStyle(2);
         line.SetLineWidth(4);
         line.Clone()->Draw();
      }

      spectra.DrawHistogram(distrOldSpectraVsPTStatErr, distrOldSpectraVsPTSysErr, kBlack,
                            0.9, 21, "old");
      spectra.DrawHistogram(distrSpectraVsPTStatErr, distrSpectraVsPTSysErr, kRed - 3,
                            0.9, 20, "new");

      ROOTTools::PrintCanvas(&canv, outputDir + "/" + resonanceName + 
                              + "_spectra_ratio_comp");
   }
}
