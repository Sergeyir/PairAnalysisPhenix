#pragma once

#include "IOTools.hpp"

#include "TCanvasTools.hpp"
#include "MathTools.hpp"

#include "PainterHelper.hpp"
#include "InputYAMLReader.hpp"

const std::string runName = "Run14HeAu200";
const std::string resonanceName = "KStar892";
const int taxiNumber = 20484;

const double pTMinRAB = 0.8;
const double pTMaxRAB = 8.6;
const double pTMinRCP = 0.4;
const double pTMaxRCP = 8.6;
const double rMin = 0.;
const double rMax = 1.999;

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

   // Spectra (old and new)
   {
      TLegend spectraLegend(0.65, 0.6, 0.95, 0.95);
      spectraLegend.SetLineColorAlpha(0, 0.);
      spectraLegend.SetFillColorAlpha(0, 0.);

      PainterHelper spectra(&spectraLegend);
      spectra.SetMarkerSize(1.4);
      spectra.SetLineWidth(2);
      spectra.SetDefaultSysWidth(0.1);

      std::vector<TH1D *> histsSpectraVsPTStatErr;
      std::vector<TH1D *> histsSpectraVsPTSysErr;
      std::vector<TGraphErrors *> graphsVladSpectraVsPTStatErr;
      std::vector<TGraphErrors *> graphsVladSpectraVsPTSysErr;

      std::vector<TF1 *> tsallisFits;

      double xMin = 1e31;
      double xMax = -1e31;
      double yMin = 1e31;
      double yMax = -1e31;

      TCanvas canv("canv", "canv", 800, 800);

      canv.SetFillStyle(0);
      canv.SetFrameFillColor(0);
      canv.SetFrameFillStyle(0);
      canv.SetFrameBorderMode(0);

      gPad->SetRightMargin(0.002); gPad->SetTopMargin(0.002); 
      gPad->SetLeftMargin(0.204); gPad->SetBottomMargin(0.105);

      gPad->SetLogy();

      for (unsigned int i = 0; i < inputYAMLResonance["centrality_bins"].size(); i++)
      {
         const std::string centralityName = 
            inputYAMLResonance["centrality_bins"][i]["name"].as<std::string>();

         const std::string inputVladSpectraFileName = "data/Spectra/pp200/KStar892.root";
         CppTools::CheckInputFile(inputVladSpectraFileName);

         TFile *inputVladSpectraFile = TFile::Open(inputVladSpectraFileName.c_str());

         histsSpectraVsPTStatErr.push_back(static_cast<TH1D *>(resultsInputFile->
            Get((centralityName + "/spectra vs pT with stat errors").c_str())->Clone()));
         histsSpectraVsPTSysErr.push_back(static_cast<TH1D *>(resultsInputFile->
            Get((centralityName + "/spectra vs pT with sys errors").c_str())->Clone()));

         graphsVladSpectraVsPTSysErr.push_back(nullptr);

         graphsVladSpectraVsPTStatErr.push_back(spectra.
               GetGraphFromTXTFile("data/Spectra/HeAu200/KStar892_" + centralityName + "_Vlad.txt", 
                                   graphsVladSpectraVsPTSysErr.back(), false, true, 0.08));

         tsallisFits.push_back(static_cast<TF1 *>(resultsInputFile->
                               Get((centralityName + "/tsallis fit").c_str())->Clone()));

         const double mult = 10./pow(10., i);

         tsallisFits.back()->SetParameter(0, tsallisFits.back()->GetParameter(0)*mult);
         histsSpectraVsPTStatErr.back()->Scale(mult);
         histsSpectraVsPTSysErr.back()->Scale(mult);

         for (int j = 0; j < graphsVladSpectraVsPTStatErr.back()->GetN(); j++)
         {
            graphsVladSpectraVsPTStatErr.back()->
               SetPointY(j, graphsVladSpectraVsPTStatErr.back()->GetPointY(j)*mult);
            graphsVladSpectraVsPTStatErr.back()->
               SetPointError(j, 0., graphsVladSpectraVsPTStatErr.back()->GetErrorY(j)*mult);
            graphsVladSpectraVsPTSysErr.back()->
               SetPointY(j, graphsVladSpectraVsPTSysErr.back()->GetPointY(j)*mult);
            graphsVladSpectraVsPTSysErr.back()->
               SetPointError(j, graphsVladSpectraVsPTSysErr.back()->GetErrorX(j),
                             graphsVladSpectraVsPTSysErr.back()->GetErrorY(j)*mult);
         }

         xMin = CppTools::Minimum(xMin, histsSpectraVsPTStatErr.back()->GetBinLowEdge(1));
         xMax = CppTools::Maximum(xMax, histsSpectraVsPTStatErr.back()->GetXaxis()->
                                  GetBinUpEdge(histsSpectraVsPTStatErr.back()->
                                  GetXaxis()->GetNbins()));

         yMin = CppTools::Minimum(yMin, histsSpectraVsPTStatErr.back()->GetMinimum());
         yMax = CppTools::Maximum(yMax, histsSpectraVsPTStatErr.back()->GetMaximum());

         tsallisFits[i]->SetLineWidth(3);
         tsallisFits[i]->SetLineStyle(2);

         tsallisFits[i]->SetLineColor(kGray + 1);
      }

      ROOTTools::DrawFrame(xMin - 0.1, yMin/5., xMax + 0.1, yMax*5., "", "#it{p}_{T} [GeV/#it{c}]",
                           "#frac{1}{2#pi#it{p}_{T}} #frac{#it{d}^{2}"\
                           "#it{N}}{#it{dp}_{T}#it{dy}} [GeV/#it{c}]^{-2}", 0.91, 1.85);

      for (int i = 0; i < static_cast<int>(inputYAMLResonance["centrality_bins"].size()); i++)
      {
         const double mult = 10./pow(10., i);

         const std::string centralityNameTex = 
            inputYAMLResonance["centrality_bins"][i]["name_tex"].as<std::string>();
         const std::string multName = 
            (i == 1) ? "" : 
            "#times10^{" + std::to_string(1 - i) + "}";

         const int color = TColor::
            GetColor(inputYAMLResonance["centrality_bins"][i]["color"].as<std::string>().c_str());
         const int markerStyle = 
            inputYAMLResonance["centrality_bins"][i]["marker_style"].as<int>();

         tsallisFits[i]->Draw("SAME");

         spectra.DrawHistogram(histsSpectraVsPTStatErr[i], histsSpectraVsPTSysErr[i], color,
                               0.9, markerStyle, centralityNameTex + multName, 0.05, false, true);
      }

      spectraLegend.AddEntry(tsallisFits.back(), "Scaled Tsallis fit", "L");

      spectraLegend.Draw();

      ROOTTools::PrintCanvas(&canv, outputDir + "/" + resonanceName + "_spectra");

      canv.Clear();

      TLegend ratioLegend(0.4, 0.8, 0.95, 0.95);
      ratioLegend.SetLineColorAlpha(0, 0.);
      ratioLegend.SetFillColorAlpha(0, 0.);

      PainterHelper ratio(&ratioLegend);
      ratio.SetMarkerSize(1.4);
      ratio.SetLineWidth(2);
      ratio.SetDefaultSysWidth(0.1);

      for (int i = 0; i < static_cast<int>(inputYAMLResonance["centrality_bins"].size()); i++)
      {
         const std::string centralityName = 
            inputYAMLResonance["centrality_bins"][i]["name"].as<std::string>();

         yMin = 1e31;
         yMax = 1e31;

         histsSpectraVsPTStatErr[i]->Divide(tsallisFits[i]);
         histsSpectraVsPTSysErr[i]->Divide(tsallisFits[i]);

         for (int j = 0; j < graphsVladSpectraVsPTStatErr[i]->GetN(); j++)
         {
            const double div = tsallisFits[i]->Eval(graphsVladSpectraVsPTStatErr[i]->GetPointX(j));

            graphsVladSpectraVsPTStatErr[i]->
               SetPointY(j, graphsVladSpectraVsPTStatErr[i]->GetPointY(j)/div);
            graphsVladSpectraVsPTStatErr[i]->
               SetPointError(j, 0., graphsVladSpectraVsPTStatErr[i]->GetErrorY(j)/div);
            graphsVladSpectraVsPTSysErr[i]->
               SetPointY(j, graphsVladSpectraVsPTSysErr[i]->GetPointY(j)/div);
            graphsVladSpectraVsPTSysErr[i]->
               SetPointError(j, graphsVladSpectraVsPTSysErr[i]->GetErrorX(j),
                             graphsVladSpectraVsPTSysErr[i]->GetErrorY(j)/div);

            yMin = CppTools::Minimum(yMin, graphsVladSpectraVsPTStatErr[i]->GetPointY(j));
            yMax = CppTools::Maximum(yMin, graphsVladSpectraVsPTStatErr[i]->GetPointY(j));
         }

         yMin = CppTools::Minimum(yMin, histsSpectraVsPTStatErr[i]->GetMinimum());
         yMax = CppTools::Maximum(yMax, histsSpectraVsPTStatErr[i]->GetMaximum());

         gPad->SetLogy(false);
         gPad->SetRightMargin(0.002); gPad->SetTopMargin(0.002); 
         gPad->SetLeftMargin(0.115); gPad->SetBottomMargin(0.112);

         ROOTTools::DrawFrame(xMin - 0.1, yMin/1.5, xMax + 0.1, yMax*1.5, 
                              "", "#it{p}_{T} [GeV/#it{c}]", "Data/Fit", 1., 1.2);

         if (yMin/1.1 < 1. && yMax*1.1 > 1.)
         {
            TLine line(xMin - 0.1, 1., xMax + 0.1, 1.);
            line.SetLineColorAlpha(kBlack, 0.5);
            line.SetLineStyle(2);
            line.SetLineWidth(4);
            line.Clone()->Draw();
         }

         ratio.DrawGraph(graphsVladSpectraVsPTStatErr[i], graphsVladSpectraVsPTSysErr[i], 
                           kBlack, 0.9, 75, "#it{K}_{Vlad}^{*0}(892)");
         ratio.DrawHistogram(histsSpectraVsPTStatErr[i], histsSpectraVsPTSysErr[i], kRed - 3,
                               0.9, 72, "(#it{K}^{*0}(892) + #bar{#it{K}}^{*0}(892))/2");

         ratioLegend.Draw();

         ROOTTools::PrintCanvas(&canv, outputDir + "/" + resonanceName + 
                                 + "_spectra_ratio_comp_" + centralityName);
         ratioLegend.Clear();
      }
   }

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
         (resultsInputFile->Get((centralName + "/spectra vs pT with stat errors").c_str())->Clone());
      TH1D *distrCentralSpectraVsPTSysErr = static_cast<TH1D *>
         (resultsInputFile->Get((centralName + "/spectra vs pT with sys errors").c_str())->Clone());
      TH1D *distrPeripheralSpectraVsPTStatErr = static_cast<TH1D *>
         (resultsInputFile->Get((peripheralName + "/spectra vs pT with stat errors").c_str())->Clone());
      TH1D *distrPeripheralSpectraVsPTSysErr = static_cast<TH1D *>
         (resultsInputFile->Get((peripheralName + "/spectra vs pT with sys errors").c_str())->Clone());

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

      const double rcpScaling = central["bias_factor"].as<double>()/
                                peripheral["bias_factor"].as<double>()*
                                peripheral["N_coll"].as<double>()/
                                central["N_coll"].as<double>();

      const double scalingUncertainty = 
         CppTools::UncertaintyProp(central["N_coll_uncertainty"].as<double>()/
                                   central["N_coll"].as<double>(),
                                   central["bias_factor_uncertainty"].as<double>()/
                                   central["bias_factor"].as<double>(),
                                   peripheral["N_coll_uncertainty"].as<double>()/
                                   peripheral["N_coll"].as<double>(),
                                   peripheral["bias_factor_uncertainty"].as<double>()/
                                   peripheral["bias_factor"].as<double>());

      distrCentralSpectraVsPTStatErr->Divide(distrPeripheralSpectraVsPTStatErr);
      distrCentralSpectraVsPTSysErr->Divide(distrPeripheralSpectraVsPTSysErr);

      distrCentralSpectraVsPTStatErr->Scale(rcpScaling);
      distrCentralSpectraVsPTSysErr->Scale(rcpScaling);

      // Vlad results
      TGraphErrors grRCPVladStatErr;
      TGraphErrors grRCPVladSysErr;

      CppTools::CheckInputFile("data/Spectra/HeAu200/KStar892_0-20_Vlad.txt");
      CppTools::CheckInputFile("data/Spectra/HeAu200/KStar892_60-88_Vlad.txt");

      std::ifstream inputFileCVlad("data/Spectra/HeAu200/KStar892_0-20_Vlad.txt");
      std::ifstream inputFilePVlad("data/Spectra/HeAu200/KStar892_60-88_Vlad.txt");

      double tmp[8];
      while (inputFileCVlad >> tmp[0] >> tmp[1] >> tmp[2] >> tmp[3] && 
             inputFilePVlad >> tmp[4] >> tmp[5] >> tmp[6] >> tmp[7])
      {
         if (tmp[0] != tmp[4]) CppTools::PrintError("Something wrong with pT bins of Vlad results");

         const double valueRCP = tmp[1]/tmp[5]*rcpScaling;
         const double statErrRCP = CppTools::UncertaintyProp(tmp[2]/tmp[1], tmp[6]/tmp[5])*valueRCP;
         const double sysErrRCP = CppTools::UncertaintyProp(tmp[3]/tmp[1], tmp[7]/tmp[5])*valueRCP;

         grRCPVladStatErr.AddPoint(tmp[0], valueRCP);
         grRCPVladSysErr.AddPoint(tmp[0], valueRCP);

         grRCPVladStatErr.SetPointError(grRCPVladStatErr.GetN() - 1, 0., statErrRCP);
         grRCPVladSysErr.SetPointError(grRCPVladSysErr.GetN() - 1, 0.08, sysErrRCP);
      }

      TCanvas canv("canv", "canv", 800, 800);

      canv.SetFillStyle(4000);
      canv.SetFrameFillColor(0);
      canv.SetFrameFillStyle(0);
      canv.SetFrameBorderMode(0);

      gPad->SetRightMargin(0.002); gPad->SetTopMargin(0.002); 
      gPad->SetLeftMargin(0.1); gPad->SetBottomMargin(0.112);

      ROOTTools::DrawFrame(pTMinRCP, rMin, pTMaxRCP, rMax, "", 
                           "#it{p}_{T} [GeV/#it{c}]", "#it{R}_{CP}", 1., 0.95);

      TLine line(pTMinRCP, 1., pTMaxRCP, 1.);
      line.SetLineColor(kGray + 1);
      line.SetLineWidth(3);
      line.SetLineStyle(2);

      line.Draw();

      PainterHelper rcp(&rabLegend);
      rcp.SetMarkerSize(1.5);
      rcp.SetLineWidth(2);
      rcp.SetDefaultSysWidth(0.1);

      rcp.DrawGraph(&grRCPVladStatErr, &grRCPVladSysErr, kBlack,
                    0.9, 75, "#it{K}_{Vlad}^{*0}(892)");

      rcp.DrawHistogram(distrCentralSpectraVsPTStatErr, distrCentralSpectraVsPTSysErr, kP6Red,
                        0.9, 72, "(#it{K}^{*0}(892) + #bar{#it{K}}^{*0}(892))/2");

      rcp.DrawLegend();

      rcp.DrawTypeCUncertainty(scalingUncertainty, pTMinRCP - 0.1, 1., kBlack, 0.3);

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

      ROOTTools::DrawFrame(pTMinRAB, rMin, pTMaxRAB, rMax, 
                           "", "#it{p}_{T} [GeV/#it{c}]", "#it{R}_{AB}", 1., 0.95);

      TLine line(pTMinRAB, 1., pTMaxRAB, 1.);
      line.SetLineColor(kGray + 1);
      line.SetLineWidth(3);
      line.SetLineStyle(2);

      line.Draw();

      PainterHelper rab(&rabLegend);
      rab.SetMarkerSize(1.5);
      rab.SetLineWidth(2);
      rab.SetDefaultSysWidth(0.08);

      rab.DrawGraphFromTXTFile("data/RAB/HeAu200/KStar892_" + centralityName + "_Vlad.txt", 
                               kBlack, 0.9, 75, "#it{K}_{Vlad}^{*0}", false, true, 0.06);

      //rab.DrawGraphFromTXTFile("data/RAB/HeAu200/ppbar" + centralityName + "PHENIX.txt", 
      //                         kBlack, 0.4, 75, "(p+#bar{p})/2, PRC109, 054910");

      rab.DrawGraphFromYAMLFile("data/RAB/HeAu200/phi1020PHENIX.yaml", centralityName, 
                                kP6Blue, 0.9, 74, "#it{#varphi}(1020), PRC106, 014982");

      rab.DrawGraphFromYAMLFile("data/RAB/HeAu200/pi0PHENIX.yaml", centralityName, 
                                kP6Violet, 0.9, 77, "#it{#pi}^{0}, PRC105 064902");

      rab.DrawHistogram(distrRABVsPTStatErr, distrRABVsPTSysErr, kP6Red, 
                        0.9, 72, "(#it{K}^{*0}(892) + #bar{#it{K}}^{*0}(892))/2");

      rab.DrawLegend();

      rab.DrawTypeCUncertainty(scalingUncertainty, 8.5, 1., kBlack, 0.3);

      TLatex tlText;

      tlText.SetTextFont(52);
      tlText.SetTextSize(0.05);

      tlText.DrawLatexNDC(0.8, 0.15, "#cbar#eta#cbar < 0.5");
      tlText.DrawLatexNDC(0.12, 0.15, (inputYAMLMain["collision_system_name_tex"].as<std::string>() 
                                       + "  " + centralityNameTex).c_str());
 
      ROOTTools::PrintCanvas(&canv, outputDir + "/" + resonanceName + "_RAB_comp_" + centralityName);

      rabLegend.Clear();
   }

   // RAB for new pp200 spectra
   const std::string spectraPPNewFileName = "data/Results/Run15pp200/20488_KStar892.root";

   if (std::filesystem::exists(spectraPPNewFileName))
   {
      TFile *resultsInputFilePP = TFile::Open(spectraPPNewFileName.c_str());

      TH1D *distrSpectraVsPTStatErrPP = static_cast<TH1D *>
         (resultsInputFilePP->Get("MB/spectra vs pT with stat errors"));
      TH1D *distrSpectraVsPTSysErrPP = static_cast<TH1D *>
         (resultsInputFilePP->Get("MB/spectra vs pT with sys errors"));

      TCanvas canv("canv", "canv", 800, 800);

      canv.SetFillStyle(4000);
      canv.SetFrameFillColor(0);
      canv.SetFrameFillStyle(0);
      canv.SetFrameBorderMode(0);

      gPad->SetRightMargin(0.002); gPad->SetTopMargin(0.002); 
      gPad->SetLeftMargin(0.1); gPad->SetBottomMargin(0.112);

      ROOTTools::DrawFrame(0.49, rMin, pTMaxRAB, rMax, 
                           "", "#it{p}_{T} [GeV/#it{c}]", "#it{R}_{AB}", 1., 0.95);

      TLine line(0.49, 1., pTMaxRAB, 1.);
      line.SetLineColor(kGray + 1);
      line.SetLineWidth(3);
      line.SetLineStyle(2);

      line.Draw();

      PainterHelper rab(&rabLegend);
      rab.SetMarkerSize(1.5);
      rab.SetLineWidth(2);
      rab.SetDefaultSysWidth(0.08);

      TLatex tlText;

      tlText.SetTextFont(52);
      tlText.SetTextSize(0.05);

      tlText.DrawLatexNDC(0.8, 0.15, "#cbar#eta#cbar < 0.5");
      tlText.DrawLatexNDC(0.12, 0.15, inputYAMLMain["collision_system_name_tex"].as<std::string>().c_str());

      for (const auto& centrality : inputYAMLResonance["centrality_bins"])
      {
         const std::string centralityName = centrality["name"].as<std::string>();
         const std::string centralityNameTex = centrality["name_tex"].as<std::string>();

         const double scaleRAB = centrality["bias_factor"].as<double>()/
                                 centrality["N_coll"].as<double>();

         /*
         const double scalingUncertainty = 
            CppTools::UncertaintyProp(centrality["N_coll_uncertainty"].as<double>()/
                                      centrality["N_coll"].as<double>(),
                                      centrality["bias_factor_uncertainty"].as<double>()/
                                      centrality["bias_factor"].as<double>());
                                      */

         TH1D *distrSpectraVsPTStatErr = static_cast<TH1D *>
            (resultsInputFile->Get((centralityName + "/spectra vs pT with stat errors").c_str()));
         TH1D *distrSpectraVsPTSysErr = static_cast<TH1D *>
            (resultsInputFile->Get((centralityName + "/spectra vs pT with sys errors").c_str()));

         distrSpectraVsPTStatErr->Divide(distrSpectraVsPTStatErrPP);
         distrSpectraVsPTSysErr->Divide(distrSpectraVsPTSysErrPP);

         distrSpectraVsPTStatErr->Scale(scaleRAB);
         distrSpectraVsPTSysErr->Scale(scaleRAB);

         const Color_t color = TColor::GetColor(centrality["color"].as<std::string>().c_str()); 

         rab.DrawHistogram(static_cast<TH1D *>(distrSpectraVsPTStatErr->Clone()), 
                           static_cast<TH1D *>(distrSpectraVsPTSysErr->Clone()), color,
                           0.8, centrality["marker_style"].as<int>(), centralityNameTex.c_str());
      }

      rab.DrawLegend();
      //rab.DrawTypeCUncertainty(scalingUncertainty, 8.5, 1., kBlack, 0.3);

      ROOTTools::PrintCanvas(&canv, outputDir + "/" + resonanceName + "_RAB_all_new_pp");

      rabLegend.Clear();
   }
   else CppTools::PrintWarning("File " + spectraPPNewFileName + " does not exists");
}
