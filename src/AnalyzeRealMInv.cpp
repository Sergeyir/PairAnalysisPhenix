/** 
 *  @file   AnalyzeRealMInv.cpp 
 *  @brief  Contains realisations of functions that are used for analyzis of expeirmental invariant mass distributions 
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysis).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef ANALYZE_REAL_M_INV_CPP
#define ANALYZE_REAL_M_INV_CPP

#include "AnalyzeRealMInv.hpp"

using namespace AnalyzeRealMInv;

int main(int argc, char **argv)
{
   if (argc < 3 || argc > 5) 
   {
      CppTools::PrintError("Expected 2-3 parameters while " + std::to_string(argc - 1) + " "\
                           "parameter(s) were provided \n Usage: bin/AnalyzeRealMInv "\
                           "inputYAMLName taxiNumber methodName=all "\
                           "numberOfThreads=std::thread::hardware_concurrency()");
   }
 
   CppTools::CheckInputFile(argv[1]);

   taxiNumber = std::stoi(argv[2]);
 
   if (argc == 4) ROOT::EnableImplicitMT(std::thread::hardware_concurrency());
   else ROOT::EnableImplicitMT(std::stoi(argv[4]));

   std::string methodToAnalyze;
   if (argc >= 3) methodToAnalyze = argv[3];
   else methodToAnalyze = "all";

   gStyle->SetOptStat(0);
   gErrorIgnoreLevel = kWarning;

   TH1::AddDirectory(kFALSE);
   TH2::AddDirectory(kFALSE);

   inputYAMLResonance.OpenFile(argv[1]);
   inputYAMLResonance.CheckStatus("resonance");

   runName = inputYAMLResonance["run_name"].as<std::string>();

   inputYAMLMain.OpenFile("input/" + runName + "/main.yaml");
   inputYAMLMain.CheckStatus("main");

   resonanceName = inputYAMLResonance["name"].as<std::string>();
   massResonance = inputYAMLResonance["mass"].as<double>();
   gammaResonance = inputYAMLResonance["gamma"].as<double>();

   daughter1Id = inputYAMLResonance["daughter1_id"].as<int>();
   daughter2Id = inputYAMLResonance["daughter2_id"].as<int>();

   hasAntiparticle = inputYAMLResonance["has_antiparticle"].as<bool>();

   minMInv = inputYAMLResonance["m_inv_range_min"].as<double>();
   maxMInv = inputYAMLResonance["m_inv_range_max"].as<double>();

   SetGaussianBroadeningFunction();

   inputFileName = "data/Real/" + runName + "/Resonance/" + std::to_string(taxiNumber) + ".root";

   CppTools::CheckInputFile(inputFileName);
   inputFile = TFile::Open(inputFileName.c_str(), "READ");

   text.SetTextFont(43);
   text.SetTextSize(45);
   texText.SetTextFont(43);
   texText.SetTextSize(45);

   pTNBins = inputYAMLResonance["pt_bins"].size();

   for (unsigned int i = 0; i < pTNBins; i++)
   {
      pTBinRanges.push_back(inputYAMLResonance["pt_bins"][i]["min"].as<double>());
   }
   pTBinRanges.push_back(inputYAMLResonance["pt_bins"][pTNBins - 1]["max"].as<double>());

   gStyle->SetPalette(kRainbow);

   int methodIndex = -1;
   if (methodToAnalyze == "all")
   {
      for (const YAML::Node& method : inputYAMLResonance["pair_selection_methods"])
      {
         numberOfIterations += (method["pt_bin_max"].as<int>() - 
                                method["pt_bin_min"].as<int>() + 1)*
                               inputYAMLResonance["centrality_bins"].size();
      }
   }
   else
   {
      for (int i = 0; i < static_cast<int>(inputYAMLResonance["pair_selection_methods"].size()); i++)
      {
         if (inputYAMLResonance["pair_selection_methods"][i]["name"].as<std::string>() ==
             methodToAnalyze) methodIndex = i;
      }
      if (methodIndex == -1)
      {
         CppTools::PrintError("No method named " + methodToAnalyze + " in input file");
      }

      numberOfIterations = 
         (inputYAMLResonance["pair_selection_methods"][methodIndex]["pt_bin_max"].as<int>() - 
          inputYAMLResonance["pair_selection_methods"][methodIndex]["pt_bin_min"].as<int>() + 1)*
         inputYAMLResonance["centrality_bins"].size();
   }

   const std::string yieldOutputDir = "data/Parameters/ResonanceEff/" + runName;
   void(system(("mkdir -p " + yieldOutputDir).c_str()));

   // file in which all important data will be written
   TFile yieldOutput((yieldOutputDir + "/" + resonanceName + ".root").c_str(), "RECREATE");

   // performing fits for specified pair selection methods
   if (methodToAnalyze == "all")
   {
      for (const YAML::Node& method : inputYAMLResonance["pair_selection_methods"])
      {
         PerformMInvFitsForMethod(method);
      }
   }
   else PerformMInvFitsForMethod(inputYAMLResonance["pair_selection_methods"][methodIndex]);

   pBar.Finish();

   CppTools::PrintInfo("AnalyzeRealMInv executable has finished running succesfully");
}

void AnalyzeRealMInv::PerformMInvFitsForMethod(const YAML::Node& method)
{
   const std::string methodName = method["name"].as<std::string>();

   const std::string outputDir = "output/MInv/" + runName + "/" + methodName;
   void(system(("mkdir -p " + outputDir).c_str()));

   /*
   TH1D distrMeansVsPT((methodName + ": means vs pT").c_str(), "", pTNBins, &pTBinRanges[0]);
   TH1D distrGammasVsPT((methodName + ": gammas vs pT").c_str(), "", pTNBins, &pTBinRanges[0]);
   TH1D distrRawYieldVsPT((methodName + ": raw yield vs pT").c_str(), 
                          "", pTNBins, &pTBinRanges[0]);
   TH1D distrRawSpectraVsPT((methodName + ": raw spectra vs pT").c_str(), 
                            "", pTNBins, &pTBinRanges[0]);
                        */

   text.SetTextAngle(270.);

   for (const YAML::Node &centralityBin : inputYAMLResonance["centrality_bins"])
   {
      for (int i = method["pt_bin_min"].as<int>(); 
           i <= method["pt_bin_max"].as<int>(); i++)
      {
         pBar.Print(static_cast<double>(numberOfCalls)/static_cast<double>(numberOfIterations));

         TH1D *distrMInvFG = nullptr;
         TH1D *distrMInvBG = nullptr;
         TH1D *distrMInvFGLR = nullptr;
         TH1D *distrMInvBGLR = nullptr;

         TH1D *distrMInv = MergeMInv(methodName, centralityBin, i, distrMInvFG, distrMInvBG, 
                                     distrMInvFGLR, distrMInvBGLR);

         if (!distrMInv)
         {
            pBar.Clear();
            CppTools::PrintError("Resulting M_{inv} histogram could not be constructed for " + 
                                 centralityBin["name"].as<std::string>() + " " +
                                 CppTools::DtoStr(pTBinRanges[i], 2) + "<pT<" + 
                                 CppTools::DtoStr(pTBinRanges[i + 1], 2));
         }
         else if (distrMInv->GetEntries() == 0)
         {
            pBar.Clear();
            CppTools::PrintWarning("Resulting histogram is empty in " + 
                                   centralityBin["name"].as<std::string>() + " " +
                                   CppTools::DtoStr(pTBinRanges[i], 2) + "<pT<" + 
                                   CppTools::DtoStr(pTBinRanges[i + 1], 2));
            pBar.RePrint();
         }

         if (!distrMInvFG)
         {
            pBar.Clear();
            CppTools::PrintError("Resulting M_{inv} foreground histogram "\
                                 "could not be constructed for " + 
                                 centralityBin["name"].as<std::string>() + " " +
                                 CppTools::DtoStr(pTBinRanges[i], 2) + "<pT<" + 
                                 CppTools::DtoStr(pTBinRanges[i + 1], 2));
         }
         else if (distrMInvFG->GetEntries() == 0)
         {
            pBar.Clear();
            CppTools::PrintWarning("Resulting foreground histogram is empty in " + 
                                   centralityBin["name"].as<std::string>() + " " +
                                   CppTools::DtoStr(pTBinRanges[i], 2) + "<pT<" + 
                                   CppTools::DtoStr(pTBinRanges[i + 1], 2));
            pBar.RePrint();
         }

         if (!distrMInvBG)
         {
            pBar.Clear();
            CppTools::PrintWarning("Resulting M_{inv} foreground histogram "\
                                   "could not be constructed for " + 
                                   centralityBin["name"].as<std::string>() + " " +
                                   CppTools::DtoStr(pTBinRanges[i], 2) + "<pT<" + 
                                   CppTools::DtoStr(pTBinRanges[i + 1], 2));
            pBar.RePrint();
         }
         else if (distrMInvBG->GetEntries() == 0)
         {
            pBar.Clear();
            CppTools::PrintWarning("Resulting background histogram is empty in " + 
                                   centralityBin["name"].as<std::string>() + 
                                   CppTools::DtoStr(pTBinRanges[i], 2) + "<pT<" + 
                                   CppTools::DtoStr(pTBinRanges[i + 1], 2));
            pBar.RePrint();
         }

         /*
         // sigma of a gaus that is convoluted with Breit-Wigner
         const double gaussianBroadeningSigma = 
            gaussianBroadeningEstimatorFunc->Eval((pTBinRanges[i] + pTBinRanges[i + 1])/2.);

         // fit for resonance+bg approximation
         TF1 fit("resonance + bg fit", &FitFunc::RBWConvGausBGGaus, 
                 massResonance - gammaResonance*3., massResonance + gammaResonance*3., 7);
         // fit for resonance approximation
         TF1 fitResonance("resonance fit", &FitFunc::RBWConvGaus, massResonance - gammaResonance*3., 
                          massResonance + gammaResonance*3., 4);
         // fit for bg approximation
         TF1 fitBG("bg fit", &FitFunc::Gaus, massResonance - gammaResonance*3., 
                   massResonance + gammaResonance*3., 3);

         const double maxBinVal = distrMInv->GetBinContent(distrMInv->GetMaximumBin());

         fit.SetParameters(maxBinVal, massResonance, gammaResonance, gaussianBroadeningSigma,
                           maxBinVal/20., massResonance, gammaResonance*4.);

         fit.SetParLimits(0, maxBinVal/1.2, maxBinVal);
         fit.SetParLimits(1, massResonance/1.1, massResonance*1.1);
         fit.SetParLimits(2, gammaResonance/1.02, gammaResonance*1.05);
         //fit.FixParameter(2, gammaResonance);
         fit.FixParameter(3, gaussianBroadeningSigma);
         fit.SetParLimits(4, 0., maxBinVal/3.);
         fit.SetParLimits(5, 0., massResonance*10.);
         fit.SetParLimits(6, gammaResonance*2., gammaResonance*20.);

         distrMInv->Fit(&fit, "RQMNBLC");

         fit.SetParLimits(0, maxBinVal/3., maxBinVal);

         for (unsigned int j = 1; j <= fitNTries; j++)
         {
            fit.SetParLimits(0, fit.GetParameter(0)/(1. + 0.05/static_cast<double>(j*j)), 
                             fit.GetParameter(0)*(1. + 0.05/static_cast<double>(j*j)));
            fit.SetParLimits(1, fit.GetParameter(1)/(1. + 0.05/static_cast<double>(j*j)), 
                             fit.GetParameter(1)*(1. + 0.05/static_cast<double>(j*j)));
            fit.SetParLimits(2, fit.GetParameter(2)/(1. + 0.1/static_cast<double>(j*j)),
                             fit.GetParameter(2)*(1. + 0.2/static_cast<double>(j*j)));
            fit.SetParLimits(5, fit.GetParameter(5)/(1. + 2./static_cast<double>(j*j)),
                             fit.GetParameter(5)*(1. + 2./static_cast<double>(j*j)));
            fit.SetParLimits(6, fit.GetParameter(6)/(1. + 1./static_cast<double>(j*j)),
                             fit.GetParameter(6)*(1. + 1./static_cast<double>(j*j)));

            fit.SetRange(fit.GetParameter(1) - (fit.GetParameter(2) + gaussianBroadeningSigma)*3., 
                         fit.GetParameter(1) + (fit.GetParameter(2) + gaussianBroadeningSigma)*3.);

            distrMInv->Fit(&fit, "RQMNBLC");
            if (j == fitNTries) distrMInv->Fit(&fit, "RQMNBLE");
         }

         distrMeansVsPT.SetBinContent(i + 1, fit.GetParameter(1));
         distrMeansVsPT.SetBinError(i + 1, fit.GetParError(1));

         distrGammasVsPT.SetBinContent(i + 1, fit.GetParameter(2));
         distrGammasVsPT.SetBinError(i + 1, fit.GetParError(2));

         for (int j = 0; j < fitResonance.GetNpar(); j++)
         {
            fitResonance.SetParameter(j, fit.GetParameter(j));
         }

         for (int j = 0; j < fitBG.GetNpar(); j++)
         {
            fitBG.SetParameter(j, fit.GetParameter(j + fitResonance.GetNpar()));
         }

         double recYieldErr;
         const double recYield = 
            GetYield(distrMInv, fitBG, fit.GetParameter(1) - fit.GetParameter(2)*2., 
                     fit.GetParameter(1) + fit.GetParameter(2)*2., recYieldErr);

         distrRawYieldVsPT.SetBinContent(i + 1, recYield/(pTBinRanges[i + 1] - pTBinRanges[i]));
         distrRawYieldVsPT.SetBinError(i + 1, recYieldErr/recYield/
                                       (pTBinRanges[i + 1] - pTBinRanges[i]));

         const double rawYieldNorm = (2.*M_PI*(pTBinRanges[i] + pTBinRanges[i + 1])/2.*
                                      (pTBinRanges[i + 1] - pTBinRanges[i]));
         distrRawSpectraVsPT.
            SetBinContent(i + 1, recYield/rawYieldNorm);
         distrRawSpectraVsPT.SetBinError(i + 1, recYieldErr/recYield/rawYieldNorm);

         fitResonance.SetRange(fit.GetParameter(1) - 
                               (fit.GetParameter(2) + gaussianBroadeningSigma)*3., 
                               fit.GetParameter(1) + 
                               (fit.GetParameter(2) + gaussianBroadeningSigma)*3.);
         fitBG.SetRange(fit.GetParameter(1) - (fit.GetParameter(2) + gaussianBroadeningSigma)*3., 
                        fit.GetParameter(1) + (fit.GetParameter(2) + gaussianBroadeningSigma)*3.);

         fit.SetLineWidth(4);
         fitResonance.SetLineWidth(4);
         fitBG.SetLineWidth(4);

         fit.SetLineColorAlpha(kRed - 3, 0.8);
         fitResonance.SetLineColorAlpha(kAzure - 3, 0.8);
         fitBG.SetLineColorAlpha(kGreen - 3, 0.8);

         fitResonance.SetLineStyle(2);
         fitBG.SetLineStyle(7);

         */

         distrMInv->GetXaxis()->SetRange(distrMInv->GetXaxis()->FindBin(minMInv), 
                                         distrMInv->GetXaxis()->FindBin(maxMInv));

         distrMInv->SetLineWidth(2);
         distrMInvFG->SetLineWidth(2);
         distrMInvBG->SetLineWidth(2);

         distrMInv->SetLineColorAlpha(kBlack, 0.8);
         distrMInv->SetMarkerColorAlpha(kBlack, 0.8);

         TCanvas canvMInv("canv MInv", "", 800, 800);

         gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05); 
         gPad->SetLeftMargin(0.173); gPad->SetBottomMargin(0.112);

         ROOTTools::DrawFrame(distrMInv, "", "#it{M}_{inv} [GeV/c^{2}]", "Counts", 1., 1.9);

         text.DrawTextNDC(0.9, 0.95, (methodName).c_str());
         texText.DrawLatexNDC(0.2, 0.88, (CppTools::DtoStr(pTBinRanges[i], 1) + " < #it{p}_{T} < " + 
                              CppTools::DtoStr(pTBinRanges[i + 1], 1)).c_str());
         /*
         texText.DrawLatexNDC(0.2, 0.83, 
                              ("#it{#chi}^{2}/NDF = " + 
                               CppTools::DtoStr(fit.GetChisquare()/fit.GetNDF(), 2)).c_str());
                               */
         /*

         texText.DrawLatexNDC(0.6, 0.9, ("#it{#mu}=" + CppTools::DtoStr(fit.GetParameter(1)*1000., 1) + 
                                         " [MeV/c^{2}]").c_str());
         texText.DrawLatexNDC(0.6, 0.85, ("#it{#Gamma}=" + 
                                          CppTools::DtoStr(fit.GetParameter(2)*1000., 1) + 
                                          " [MeV/c^{2}]").c_str());
                                          */

         /*
         fitBG.Draw("SAME");
         fitResonance.Draw("SAME");
         fit.Draw("SAME");
         */

         ROOTTools::PrintCanvas(&canvMInv, outputDir + "/" + resonanceName + "_" + 
                                centralityBin["name"].as<std::string>() + "_" +
                                CppTools::DtoStr(pTBinRanges[i], 1) + "-" + 
                                CppTools::DtoStr(pTBinRanges[i + 1], 1));

         distrMInv->SetFillStyle(3001);
         distrMInvFG->SetFillStyle(3001);
         distrMInvBG->SetFillStyle(3001);

         TCanvas canvMInvSummary("canv MInv summary", "", 800, 800);

         canvMInvSummary.Divide(2, 2);
         
         canvMInvSummary.cd(1);

         gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05); 
         gPad->SetLeftMargin(0.173); gPad->SetBottomMargin(0.112);

         ROOTTools::DrawFrame(static_cast<TH1D *>(distrMInv->Clone()), 
                              "", "#it{M}_{inv} [GeV/c^{2}]", "Counts", 1., 1.9);

         canvMInvSummary.cd(3);

         gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05); 
         gPad->SetLeftMargin(0.173); gPad->SetBottomMargin(0.112);

         ROOTTools::DrawFrame(static_cast<TH1D *>(distrMInv->Clone()), 
                              "", "#it{M}_{inv} [GeV/c^{2}]", "Counts", 1., 1.9);

         canvMInvSummary.cd(4);

         gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05); 
         gPad->SetLeftMargin(0.173); gPad->SetBottomMargin(0.112);

         if (distrMInvBG && distrMInvFG)
         {
            TH1D *ratioFGBG = static_cast<TH1D *>(distrMInvFG->Clone());
            ratioFGBG->Divide(distrMInvBG);

            ratioFGBG->GetXaxis()->SetRange(ratioFGBG->GetXaxis()->FindBin(minMInv), 
                                            ratioFGBG->GetXaxis()->FindBin(maxMInv));

            ratioFGBG->SetLineWidth(2);
            ratioFGBG->SetLineColor(kBlack);
            ratioFGBG->SetMarkerColor(kBlack);

            ROOTTools::DrawFrame(ratioFGBG, "", "#it{M}_{inv} [GeV/c^{2}]", "FG/BG", 1., 1.9);
         }
         else
         {
            text.DrawTextNDC(0.75, 0.9, "No data");
         }

         text.DrawTextNDC(0.85, 0.9, (methodName).c_str());

         canvMInvSummary.cd(2);

         gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05); 
         gPad->SetLeftMargin(0.173); gPad->SetBottomMargin(0.112);

         distrMInv->GetXaxis()->SetRange(1, distrMInv->GetXaxis()->GetNbins());

         if (distrMInvFG)
         {
            ROOTTools::DrawFrame(distrMInvFG, "", "#it{M}_{inv} [GeV/c^{2}]", "Counts", 1., 1.9,
                                 0.05, 0.05, true, false);
         }
         else
         {
            ROOTTools::DrawFrame(distrMInv, "", "#it{M}_{inv} [GeV/c^{2}]", "Counts", 1., 1.9,
                                 0.05, 0.05, true, false);
         }

         distrMInv->Sumw2(false);
         distrMInvBG->Sumw2(false);
         distrMInvFG->Sumw2(false);

         if (distrMInvFG) 
         {
            distrMInvFG->SetLineColorAlpha(kAzure + 2, 0.8);
            distrMInvFG->Draw("SAME PFC");
         }
         else text.DrawTextNDC(0.85, 0.9, "No data on foreground");

         if (distrMInvBG) 
         {
            distrMInvBG->SetLineColorAlpha(kGreen + 2, 0.8);
            distrMInvBG->Draw("SAME PFC");
         }
         else text.DrawTextNDC(0.75, 0.9, "No data on background");

         distrMInv->SetLineColorAlpha(kRed + 2, 0.8);
         distrMInv->Draw("SAME PFC");

         texText.DrawLatexNDC(0.3, 0.85, (CppTools::DtoStr(pTBinRanges[i], 1) + " < #it{p}_{T} < " + 
                              CppTools::DtoStr(pTBinRanges[i + 1], 1)).c_str());

         ROOTTools::PrintCanvas(&canvMInvSummary, outputDir + "/Summary_" + resonanceName + "_" + 
                                centralityBin["name"].as<std::string>() + "_" +
                                CppTools::DtoStr(pTBinRanges[i], 1) + "-" + 
                                CppTools::DtoStr(pTBinRanges[i + 1], 1), true, false);

         TCanvas canvMInvFGVsBG("canv MInv FG vs BG", "", 800, 800);

         numberOfCalls++;
      }
   }

   text.SetTextAngle(0.);

   /*
   distrMeansVsPT.SetLineColor(kRed - 2);
   distrMeansVsPT.SetMarkerColor(kRed - 2);
   distrMeansVsPT.SetLineWidth(4);

   distrGammasVsPT.SetLineColor(kRed - 2);
   distrGammasVsPT.SetMarkerColor(kRed - 2);
   distrGammasVsPT.SetLineWidth(4);

   distrRawYieldVsPT.SetLineColor(kRed - 2);
   distrRawYieldVsPT.SetMarkerColor(kRed - 2);
   distrRawYieldVsPT.SetLineWidth(4);

   distrRawSpectraVsPT.SetLineColor(kRed - 2);
   distrRawSpectraVsPT.SetMarkerColor(kRed - 2);
   distrRawSpectraVsPT.SetLineWidth(4);

   TLine massResonancePDG(pTBinRanges[0], massResonance, pTBinRanges[pTNBins], massResonance);
   massResonancePDG.SetLineColorAlpha(kBlack, 0.5);
   massResonancePDG.SetLineStyle(2);
   massResonancePDG.SetLineWidth(4);

   TLine gammaResonancePDG(pTBinRanges[0], gammaResonance, pTBinRanges[pTNBins], gammaResonance);
   gammaResonancePDG.SetLineColorAlpha(kBlack, 0.5);
   gammaResonancePDG.SetLineStyle(2);
   gammaResonancePDG.SetLineWidth(4);

   distrMeansVsPT.SetMaximum(massResonance*1.025);
   distrMeansVsPT.SetMinimum(massResonance*0.975);

   distrGammasVsPT.SetMaximum(gammaResonance*1.5);
   distrGammasVsPT.SetMinimum(gammaResonance/2.);

   TCanvas canvMeansVsPT("canv means vs pT", "", 800, 800);

   gPad->SetRightMargin(0.03);
   gPad->SetTopMargin(0.02);
   gPad->SetLeftMargin(0.172);
   gPad->SetBottomMargin(0.112);

   ROOTTools::DrawFrame(&distrMeansVsPT, "", "#it{p}_{T} [GeV/#it{c}]", 
                        "#it{#mu} [GeV/#it{c}^{2}]", 1., 1.82);

   text.DrawTextNDC(0.38, 0.9, (methodName).c_str());
   massResonancePDG.Draw();

   ROOTTools::PrintCanvas(&canvMeansVsPT, outputDir + "/" + resonanceName + "_means");

   TCanvas canvGammasVsPT("canv gammas vs pT", "", 800, 800);

   gPad->SetRightMargin(0.03);
   gPad->SetTopMargin(0.02);
   gPad->SetLeftMargin(0.172);
   gPad->SetBottomMargin(0.112);

   ROOTTools::DrawFrame(&distrGammasVsPT, "", "#it{p}_{T} [GeV/#it{c}]", 
                        "#it{#Gamma} [GeV/#it{c}^{2}]", 1., 1.82);

   text.DrawTextNDC(0.38, 0.9, (methodName).c_str());
   gammaResonancePDG.Draw();

   ROOTTools::PrintCanvas(&canvGammasVsPT, outputDir + "/" + resonanceName + "_gammas");

   TCanvas canvRawYieldVsPT("canv raw yield vs pT", "", 800, 800);

   gPad->SetLogy();

   gPad->SetRightMargin(0.03);
   gPad->SetTopMargin(0.02);
   gPad->SetLeftMargin(0.141);
   gPad->SetBottomMargin(0.112);

   ROOTTools::DrawFrame(&distrRawYieldVsPT, "", "#it{p}_{T} [GeV/#it{c}]", 
                        "#it{dY}_{raw}/#it{dp}_{T} [(GeV/#it{c})^{-1}]", 1., 1.35);

   text.DrawTextNDC(0.38, 0.9, (methodName).c_str());

   ROOTTools::PrintCanvas(&canvRawYieldVsPT, outputDir + "/" + resonanceName + "_raw_yield");

   TCanvas canvRawSpectraVsPT("canv raw spectra vs pT", "", 800, 800);

   gPad->SetLogy();

   gPad->SetRightMargin(0.03);
   gPad->SetTopMargin(0.02);
   gPad->SetLeftMargin(0.141);
   gPad->SetBottomMargin(0.112);

   ROOTTools::DrawFrame(&distrRawSpectraVsPT, "", "#it{p}_{T} [GeV/#it{c}]", 
                        "1/(2 #pi #it{p}_{T}) #it{d}^{2}#it{Y}_{raw}/#it{dp}_{T}#it{dy} "\
                        "[(GeV/#it{c})^{-2})]", 1., 1.35);

   text.DrawTextNDC(0.38, 0.9, (methodName).c_str());

   ROOTTools::PrintCanvas(&canvRawSpectraVsPT, outputDir + "/" + resonanceName + "_raw_spectra");
   */
}

TH1D *AnalyzeRealMInv::MergeMInv(const std::string& methodName, const YAML::Node& centralityBin,
                                 const int pTBin, TH1D*& distrMInvMergedFG, TH1D*& distrMInvMergedBG,
                                 TH1D*& distrMInvMergedFGLR, TH1D*& ditrMInvmMergedBGLR)
{
   if (distrMInvMergedFG || distrMInvMergedFG || distrMInvMergedFGLR || distrMInvMergedBGLR)
   {
      CppTools::PrintError("TH1D *AnalyzeRealMInv::MergeMInv(const std::string& methodName, "
                           "const YAML::Node& centralityBin, TH1D*& distrMInvMergedFG, "
                           "TH1D*& distrMInvMergedBG, TH1D*& distrMInvMergedFGLR, "
                           "TH1D*& distrMInvMergedBGLR)) : \n" 
                           (!distrMInvMergedFG ? "" : "distrMInvMergedFG ")
                           (!distrMInvMergedBG ? "" : "distrMInvMergedBG ")
                           (!distrMInvMergedFGLR ? "" : "distrMInvMergedFGLR ")
                           (!distrMInvMergedBGLR ? "" : "distrMInvMergedBGLR ")
                           "value(s) must be nullptr");
   }

   // decay channel name
   const std::string channelName = ParticleMap::nameShort[daughter1Id] +
                                   ParticleMap::nameShort[daughter2Id];

   TH1D *distrMInvMerged = nullptr;
   // iterating over CabanaBoy centrality bins
   for (int c = centralityBin["cb_c_bin_min"].as<int>(); 
        c <= centralityBin["cb_c_bin_max"].as<int>(); c++)
   {
      const std::string cName = (c > 9) ? std::to_string(c) : "0" + std::to_string(c);

      // iterating over CabanaBoy z_{vtx} bins
      for (int z = 0; z < inputYAMLResonance["cb_z_bins"].as<int>(); z++)
      {
         const std::string zName = (z > 9) ? std::to_string(z) : "0" + std::to_string(z);

         // iterating over CabanaBoy r_{vtx} bins
         for (int r = 0; r < inputYAMLResonance["cb_r_bins"].as<int>(); r++)
         {
            const std::string rName = (r > 9) ? std::to_string(r) : "0" + std::to_string(r);

            const std::string distrMInvVsPTFGName = 
               "c" + cName + "_z" + zName + "_r" + rName + "/" + 
                methodName + ": " + channelName + "_FG12";

            TH2F *distrMInvVsPTFG = 
               static_cast<TH2F *>(inputFile->Get(distrMInvVsPTFGName.c_str()));

            if (!distrMInvVsPTFG)
            {
               pBar.Clear();
               CppTools::PrintError("Histogram named " + distrMInvVsPTFGName + 
                                    " does not exist in file " + inputFileName);
            }

            const std::string distrMInvVsPTBGName = 
               "c" + cName + "_z" + zName + "_r" + rName + "/" + 
                methodName + ": " + channelName + "_BG12";

            TH2F *distrMInvVsPTBG = 
               static_cast<TH2F *>(inputFile->Get(distrMInvVsPTBGName.c_str()));

            if (!distrMInvVsPTBG)
            {
               pBar.Clear();
               CppTools::PrintError("Histogram named " + distrMInvVsPTBGName + 
                                    " does not exist in file " + inputFileName);
            }

            for (int i = distrMInvVsPTFG->GetXaxis()->FindBin(pTBinRanges[pTBin] + 1e-6); 
                 i <= distrMInvVsPTFG->GetXaxis()->FindBin(pTBinRanges[pTBin + 1] - 1e-6); i++)
            {
               TH1D *distrMInvFG = 
                  distrMInvVsPTFG->ProjectionY((distrMInvVsPTFG->GetName() + 
                                                std::to_string(c) + std::to_string(z) + 
                                                std::to_string(r) + std::to_string(i)).c_str(), 
                                               i, i);
               TH1D *distrMInvBG = 
                  distrMInvVsPTBG->ProjectionY((distrMInvVsPTBG->GetName() + 
                                                std::to_string(c) + std::to_string(z) + 
                                                std::to_string(r) + std::to_string(i)).c_str(), 
                                               i, i);

               if (distrMInvFG->GetEntries() < 1e-3) continue;

               distrMInvFG->Sumw2();
               distrMInvBG->Sumw2();

               if (!distrMInvMerged) 
               {
                  if (distrMInvBG->GetEntries() < 1e-3) 
                  {
                     distrMInvMerged = static_cast<TH1D *>(distrMInvFG->Clone());
                  }
                  else distrMInvMerged = SubtractBG(distrMInvFG, distrMInvBG);

                  distrMInvMergedFG = static_cast<TH1D *>(distrMInvFG->Clone());
                  distrMInvMergedBG = static_cast<TH1D *>(distrMInvBG->Clone());
               }
               else 
               {
                  if (distrMInvBG->GetEntries() < 1e-3) distrMInvMerged->Add(distrMInvFG);
                  else distrMInvMerged->Add(SubtractBG(distrMInvFG, distrMInvBG));

                  distrMInvMergedFG->Add(distrMInvFG);
                  distrMInvMergedBG->Add(distrMInvBG);
               }
            }
         }
      }
   }

   return distrMInvMerged;
}

TH1D *AnalyzeRealMInv::SubtractBG(TH1D*& distrMInvFG, TH1D*& distrMInvBG)
{
   const double integralMInvFG = distrMInvFG->Integral(1, distrMInvFG->GetXaxis()->GetNbins());
   const double integralMInvBG = distrMInvBG->Integral(1, distrMInvBG->GetXaxis()->GetNbins());

   double partIntegralMInvFG = 0.;
   double partIntegralMInvBG = 0.;

   for (int i = distrMInvFG->GetXaxis()->GetNbins(); i >= 1; i--)
   {
      if (partIntegralMInvFG > integralMInvFG*0.1 && 
          partIntegralMInvBG > integralMInvBG*0.1) break;

      partIntegralMInvFG += distrMInvFG->GetBinContent(i);
      partIntegralMInvBG += distrMInvBG->GetBinContent(i);
   }

   double scaleFactorBG = partIntegralMInvFG/partIntegralMInvBG;

   partIntegralMInvFG = 0.;
   partIntegralMInvBG = 0.;

   // if background was overestimated, some bins will be negative after BG subtraction
   // this needs to be resolved by rescaling (if needed)
   for (int i = 100; i <= distrMInvFG->GetXaxis()->GetNbins(); i++)
   {
      if (distrMInvBG->GetBinContent(i) < 1e-3) continue;
      if (distrMInvFG->GetBinContent(i) < 1e-3) continue;

      if (partIntegralMInvFG > integralMInvFG*0.1 && 
          partIntegralMInvBG > integralMInvBG*0.1) break;

      partIntegralMInvFG += distrMInvFG->GetBinContent(i);
      partIntegralMInvBG += distrMInvBG->GetBinContent(i);

      if (distrMInvBG->GetBinContent(i)*scaleFactorBG > distrMInvFG->GetBinContent(i))
      {
         CppTools::Print(distrMInvFG->GetBinContent(i)/distrMInvBG->GetBinContent(i));
         scaleFactorBG *= distrMInvFG->GetBinContent(i)/distrMInvBG->GetBinContent(i);
      }
   }

   distrMInvBG->Scale(scaleFactorBG);

   TH1D *distrMInvSubtr = static_cast<TH1D *>(distrMInvFG->Clone());
   distrMInvSubtr->Add(distrMInvBG, -1.);

   return distrMInvSubtr;
}

void AnalyzeRealMInv::SetGaussianBroadeningFunction()
{
   const std::string inputFileName = "data/Parameters/GaussianBroadening/" + 
                                     runName + "/" + resonanceName + ".root";

   if (!CppTools::FileExists(inputFileName))
   {
      CppTools::PrintError(inputFileName + " does not exists. "\
                           "Run executable bin/EstimateGassianBroadening first");
   }

   gaussianBroadeningEstimatorFunc = 
      static_cast<TF1 *>(TFile::Open(inputFileName.c_str())->Get("gaussian broadening sigma fit"));
}

double AnalyzeRealMInv::GetYield(TH1D *distrMInv, const TF1& funcBG, 
                                 const double xMin, const double xMax, double &err)
{
   // integral over the signal
   double integral = 0.;
   // integral over the background
   double integralBG = 0.;

   for (int i = distrMInv->GetXaxis()->FindBin(xMin); 
        i <= distrMInv->GetXaxis()->FindBin(xMax); i++)
   {
      integral += distrMInv->GetBinContent(i);

      // integrating background over the single bin for better estimation
      for (double m = distrMInv->GetXaxis()->GetBinLowEdge(i); 
           m <= distrMInv->GetXaxis()->GetBinUpEdge(i); 
           m += distrMInv->GetXaxis()->GetBinWidth(i)/100.)
      {
         integralBG += funcBG.Eval(m);
      }
   }
   // due to the spectra scaling statistical uncertainty is not tied to the integral 
   // but rather to the number of entries of the signal
   err = sqrt(distrMInv->GetEntries()*integral/
              distrMInv->Integral(1, distrMInv->GetXaxis()->GetNbins()));
   // normalizing background integral by the number of integration steps
   integral -= integralBG/101.;

   return integral;
}

#endif /* ANALYZE_REAL_M_INV_CPP */
