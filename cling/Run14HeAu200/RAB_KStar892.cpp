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

void RAB_KStar892()
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

   TFile *resultsInputFile = TFile::Open(resultsInputFileName.c_str());

   const std::string outputDir = "output/Results/" + runName + "/" + std::to_string(taxiNumber);
   std::filesystem::create_directories(outputDir);

   for (const auto& centrality : inputYAMLResonance["centrality_bins"])
   {
      const std::string centralityName = centrality["name"].as<std::string>();
      const std::string centralityNameTex = centrality["name_tex"].as<std::string>();

      const double scalingUncertainty = 
         CppTools::UncertaintyProp(centrality["N_coll_uncertainty"].as<double>()/
                                   centrality["N_coll"].as<double>(),
                                   centrality["bias_factor_uncertainty"].as<double>()/
                                   centrality["bias_factor"].as<double>());

      TCanvas canv("canv", "canv", 800, 800);

      canv.SetFillStyle(4000);
      canv.SetFrameFillColor(0);
      canv.SetFrameFillStyle(0);
      canv.SetFrameBorderMode(0);

      gPad->SetRightMargin(0.002); gPad->SetTopMargin(0.002); 
      gPad->SetLeftMargin(0.1); gPad->SetBottomMargin(0.112);

      ROOTTools::DrawFrame(xMin, yMin, xMax, yMax, "", "p_{T} [GeV/c]", "R_{AB}", 1., 0.95);

      TLegend legend(0.4, 0.8, 0.95, 0.95);
      legend.SetLineColorAlpha(0, 0.);
      legend.SetFillColorAlpha(0, 0.);

      TLine line(xMin, 1., xMax, 1.);
      line.SetLineColor(kGray + 1);
      line.SetLineWidth(3);
      line.SetLineStyle(2);

      line.Draw();

      PainterHelper rab(&legend);
      rab.SetMarkerSize(1.6);
      rab.SetLineWidth(2);
      rab.SetSysWidth(0.08);

      TH1D *distrRABVsPT = static_cast<TH1D *>
         (resultsInputFile->Get((centralityName + "/RAB vs pT with stat errors").c_str()));

      if (!distrRABVsPT) CppTools::PrintError("RAB distribution was not found in file " + 
                                              resultsInputFileName);

      rab.DrawGraphFromTXTFile("data/RAB/HeAu200/KStar892_" + centralityName + "_Vlad.txt", 
                               kGray + 1, 0.9, 75, "K_{Vlad}^{*0}(892)");

      //rab.DrawGraphFromTXTFile("data/RAB/HeAu200/ppbar" + centralityName + "PHENIX.txt", 
      //                         kBlack, 0.4, 75, "(p+#bar{p})/2 PRC109, 054910");

      rab.DrawGraphFromYAMLFile("data/RAB/HeAu200/phi1020PHENIX.yaml", centralityName, 
                                kAzure - 3, 0.9, 74, "#varphi(1020) PRC106, 014982");

      rab.DrawGraphFromYAMLFile("data/RAB/HeAu200/pi0PHENIX.yaml", centralityName, 
                                kGreen - 3, 0.9, 77, "#pi^{0} PRC105 064902");

      rab.DrawHistogram(distrRABVsPT, distrRABVsPT, kRed - 3, 
                        0.9, 72, "K^{*0}(892)");

      rab.DrawLegend();

      rab.DrawTypeCUncertainty(scalingUncertainty, 8.5, 1., kBlack, 0.3);

      TLatex tlText;

      tlText.SetTextFont(52);
      tlText.SetTextSize(0.05);

      tlText.DrawLatexNDC(0.8, 0.15, "#cbar#eta#cbar < 0.5");
      tlText.DrawLatexNDC(0.12, 0.15, (inputYAMLMain["collision_system_name_tex"].as<std::string>() + "  " + centralityNameTex).c_str());
 
      ROOTTools::PrintCanvas(&canv, outputDir + "/RAB_comp_" + centralityName);
   }
}
