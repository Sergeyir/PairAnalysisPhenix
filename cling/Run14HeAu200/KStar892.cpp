#pragma once

#include "IOTools.hpp"

#include "TCanvasTools.hpp"
#include "MathTools.hpp"

#include "PainterHelper.hpp"
#include "InputYAMLReader.hpp"

const std::string runName = "Run14HeAu200";
const std::string resonanceName = "KStar892";
const int taxiNumber = 20292;

const double xMin = 0.8;
const double xMax = 8.6;
const double yMin = 0.;
const double yMax = 1.999;

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

   const std::string outputDir = "output/Results/" + runName + "/" + std::to_string(taxiNumber);
   std::filesystem::create_directories(outputDir);

   TLegend rabLegend(0.4, 0.8, 0.95, 0.95);
   rabLegend.SetLineColorAlpha(0, 0.);
   rabLegend.SetFillColorAlpha(0, 0.);

   // RCP for K*
   {
      const unsigned int centralIndex = 1;
      const unsigned int peripheralIndex = inputYAMLResonance["centrality_bins"].size() - 1;

      const YAML::Node central = inputYAMLResonance["centrality_bins"][centralIndex];
      const YAML::Node peripheral = inputYAMLResonance["centrality_bins"][peripheralIndex];

      const std::string centralName = central["name"].as<std::string>();
      const std::string centralNameTex = central["name_tex"].as<std::string>();

      const std::string peripheralName = peripheral["name"].as<std::string>();

      TH1D *distrCentralSpectraVsPTStatErr = static_cast<TH1D *>
         (resultsInputFile->Get((centralName + "/spectra vs pT with stat errors").c_str()));
      TH1D *distrCentralSpectraVsPTSysErr = static_cast<TH1D *>
         (resultsInputFile->Get((centralName + "/spectra vs pT with sys errors").c_str()));
      TH1D *distrPeripheralSpectraVsPTStatErr = static_cast<TH1D *>
         (resultsInputFile->Get((peripheralName + "/spectra vs pT with stat errors").c_str()));
      TH1D *distrPeripheralSpectraVsPTSysErr = static_cast<TH1D *>
         (resultsInputFile->Get((peripheralName + "/spectra vs pT with sys errors").c_str()));

      if (!distrCentralSpectraVsPTStatErr) 
      {
         CppTools::PrintError("Spectra distribution for " + centralName + 
                              " with statistical errors was not "\
                              "found in file " + resultsInputFileName);
      }
      if (!distrCentralSpectraVsPTSysErr) 
      {
         CppTools::PrintError("Spectra distribution for " + centralName + 
                              " with systematic errors was not "\
                              "found in file " + resultsInputFileName);
      }
      if (!distrPeripheralSpectraVsPTStatErr) 
      {
         CppTools::PrintError("Spectra distribution for " + peripheralName + 
                              " with statistical errors was not "\
                              "found in file " + resultsInputFileName);
      }
      if (!distrPeripheralSpectraVsPTSysErr) 
      {
         CppTools::PrintError("Spectra distribution for " + peripheralName + 
                              " with systematic errors was not "\
                              "found in file " + resultsInputFileName);
      }

      const double scalingUncertainty = 
         CppTools::UncertaintyProp(central["N_coll_uncertainty"].as<double>()/
                                   central["N_coll"].as<double>(),
                                   central["bias_factor_uncertainty"].as<double>()/
                                   central["bias_factor"].as<double>(),
                                   peripheral["N_coll_uncertainty"].as<double>()/
                                   peripheral["N_coll"].as<double>(),
                                   peripheral["bias_factor_uncertainty"].as<double>()/
                                   peripheral["bias_factor"].as<double>());

      TCanvas canv("canv", "canv", 800, 800);

      canv.SetFillStyle(4000);
      canv.SetFrameFillColor(0);
      canv.SetFrameFillStyle(0);
      canv.SetFrameBorderMode(0);

      gPad->SetRightMargin(0.002); gPad->SetTopMargin(0.002); 
      gPad->SetLeftMargin(0.1); gPad->SetBottomMargin(0.112);

      ROOTTools::DrawFrame(xMin, yMin, 6.5, yMax, "", "p_{T} [GeV/c]", "R_{CP}", 1., 0.95);

      TLine line(xMin, 1., xMax, 1.);
      line.SetLineColor(kGray + 1);
      line.SetLineWidth(3);
      line.SetLineStyle(2);

      line.Draw();

      PainterHelper rcp(&rabLegend);
      rcp.SetMarkerSize(1.5);
      rcp.SetLineWidth(2);
      rcp.SetSysWidth(0.08);

      distrCentralSpectraVsPTStatErr->Divide(distrPeripheralSpectraVsPTStatErr);
      distrCentralSpectraVsPTSysErr->Divide(distrPeripheralSpectraVsPTSysErr);

      const double rcpScaling = central["bias_factor"].as<double>()/
                                peripheral["bias_factor"].as<double>()*
                                peripheral["N_coll"].as<double>()/
                                central["N_coll"].as<double>();

      distrCentralSpectraVsPTStatErr->Scale(rcpScaling);
      distrCentralSpectraVsPTSysErr->Scale(rcpScaling);

      rcp.DrawHistogram(distrCentralSpectraVsPTStatErr, distrCentralSpectraVsPTSysErr, kP6Red, 
                        0.9, 72, "(K^{*0}(892) + #bar{K}^{*0}(892))/2");

      rcp.DrawLegend();

      rcp.DrawTypeCUncertainty(scalingUncertainty, 6.4, 1., kBlack, 0.3);

      TLatex tlText;

      tlText.SetTextFont(52);
      tlText.SetTextSize(0.05);

      tlText.DrawLatexNDC(0.8, 0.15, "#cbar#eta#cbar < 0.5");
      tlText.DrawLatexNDC(0.12, 0.15, (inputYAMLMain["collision_system_name_tex"].as<std::string>() 
                                       + "  " + centralName + "/" + peripheralName).c_str());

      ROOTTools::PrintCanvas(&canv, outputDir + "/" + resonanceName + "_RCP");

      rabLegend.Clear();
   }

   // RAB for all particles
   for (const auto& centrality : inputYAMLResonance["centrality_bins"])
   {
      const std::string centralityName = centrality["name"].as<std::string>();
      const std::string centralityNameTex = centrality["name_tex"].as<std::string>();

      const double scalingUncertainty = 
         CppTools::UncertaintyProp(centrality["N_coll_uncertainty"].as<double>()/
                                   centrality["N_coll"].as<double>(),
                                   centrality["bias_factor_uncertainty"].as<double>()/
                                   centrality["bias_factor"].as<double>());

      TH1D *distrRABVsPTStatErr = static_cast<TH1D *>
         (resultsInputFile->Get((centralityName + "/RAB vs pT with stat errors").c_str()));
      TH1D *distrRABVsPTSysErr = static_cast<TH1D *>
         (resultsInputFile->Get((centralityName + "/RAB vs pT with sys errors").c_str()));

      if (!distrRABVsPTStatErr) 
      {
         CppTools::PrintError("RAB distribution with statistical errors was not "\
                              "found in file " + resultsInputFileName);
      }
      if (!distrRABVsPTSysErr) 
      {
         CppTools::PrintError("RAB distribution with systematic errors was not "\
                              "found in file " + resultsInputFileName);
      }

      TCanvas canv("canv", "canv", 800, 800);

      canv.SetFillStyle(4000);
      canv.SetFrameFillColor(0);
      canv.SetFrameFillStyle(0);
      canv.SetFrameBorderMode(0);

      gPad->SetRightMargin(0.002); gPad->SetTopMargin(0.002); 
      gPad->SetLeftMargin(0.1); gPad->SetBottomMargin(0.112);

      ROOTTools::DrawFrame(xMin, yMin, xMax, yMax, "", "p_{T} [GeV/c]", "R_{AB}", 1., 0.95);

      TLine line(xMin, 1., xMax, 1.);
      line.SetLineColor(kGray + 1);
      line.SetLineWidth(3);
      line.SetLineStyle(2);

      line.Draw();

      PainterHelper rab(&rabLegend);
      rab.SetMarkerSize(1.5);
      rab.SetLineWidth(2);
      rab.SetSysWidth(0.08);

      rab.DrawGraphFromTXTFile("data/RAB/HeAu200/KStar892_" + centralityName + "_Vlad.txt", 
                               kBlack, 0.9, 75, "K_{Vlad}^{*0}");

      //rab.DrawGraphFromTXTFile("data/RAB/HeAu200/ppbar" + centralityName + "PHENIX.txt", 
      //                         kBlack, 0.4, 75, "(p+#bar{p})/2, PRC109, 054910");

      rab.DrawGraphFromYAMLFile("data/RAB/HeAu200/phi1020PHENIX.yaml", centralityName, 
                                kP6Blue, 0.9, 74, "#varphi(1020), PRC106, 014982");

      rab.DrawGraphFromYAMLFile("data/RAB/HeAu200/pi0PHENIX.yaml", centralityName, 
                                kP6Violet, 0.9, 77, "#pi^{0}, PRC105 064902");

      rab.DrawHistogram(distrRABVsPTStatErr, distrRABVsPTSysErr, kP6Red, 
                        0.9, 72, "(K^{*0}(892) + #bar{K}^{*0}(892))/2");

      rab.DrawLegend();

      rab.DrawTypeCUncertainty(scalingUncertainty, 8.5, 1., kBlack, 0.3);

      TLatex tlText;

      tlText.SetTextFont(52);
      tlText.SetTextSize(0.05);

      tlText.DrawLatexNDC(0.8, 0.15, "#cbar#eta#cbar < 0.5");
      tlText.DrawLatexNDC(0.12, 0.15, (centralityNameTex + "  " + centralityNameTex).c_str());
 
      ROOTTools::PrintCanvas(&canv, outputDir + "/" + resonanceName + "_RAB_comp_" + centralityName);

      rabLegend.Clear();
   }
}
