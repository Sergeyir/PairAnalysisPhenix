/** 
 *  @file   MInvFit.cpp 
 *  @brief  Contains implementation for usage of ROOTTools::GUIFit for tweaking bad background approximations on PHENIX processed data
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysis).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#pragma once

#include "ErrorHandler.hpp"
#include "IOTools.hpp"
#include "StrTools.hpp"

#include "GUIFit.hpp"

#include "AnalyzeRealMInv.hpp"

using namespace AnalyzeRealMInv;

void MInvFit()
{
	gStyle->SetOptStat(0);
   gErrorIgnoreLevel = kWarning;

   TH1::AddDirectory(kFALSE);
   TH2::AddDirectory(kFALSE);

   ROOT::EnableImplicitMT();
   gROOT->SetBatch(true);

   /*
   CppTools::PrintInfo("List of directories in data/Real directory that contain "\
                       "non empty \"Resonance\" directory");
   system("find data/Real/ -name Resonance -type d -not -empty");

   CppTools::Print("Choose the run (part of a directory, example: Run14HeAu200) "\
                   "from the above directory list and type it in");
   std::cout << ">> ";
   std::cin >> runName;
   */
   runName = "Run14HeAu200"; // temporary for testing

   /*
   CppTools::PrintInfo("List of .root files in " + runName + " taxi directory");
   system(("find data/Real/" + runName + "/Resonance/ -type f -name *.root").c_str());

   CppTools::Print("Choose taxi number (part of .root file, example: 20025) "\
                   "from the above directory list and type it in");
   int taxiNumber;
   std::cout << ">> ";
   std::cin >> taxiNumber;
   */
   int taxiNumber = 20025; // temporary for testing

   inputFileName = "data/Real/" + runName + "/Resonance/" + std::to_string(taxiNumber) + ".root";
   CppTools::CheckInputFile(inputFileName);

   /* 
   CppTools::Print("Choose the particle");
   std::string particleName;
   std::cout << ">> ";
   std::cin >> particleName;
   */
   const std::string particleName = "KStar892";

   inputYAMLResonance.OpenFile("input/" + runName + "/" + particleName + ".yaml");
   inputYAMLResonance.CheckStatus("resonance");

   InputYAMLReader inputYAMLMain("input/" + runName + "/main.yaml");
   inputYAMLMain.CheckStatus("main");

   int methodBinIndex;
   while (true) // ininite loop until exit or valid input is specified
   {
      CppTools::PrintInfo("List of pair selection method bins for " + 
                          particleName + " in " + runName);

      for (int i = 0; i < inputYAMLResonance["pair_selection_methods"].size(); i++)
      {
         CppTools::Print(i, inputYAMLResonance["pair_selection_methods"][i]
                                              ["name"].as<std::string>());
      }

      CppTools::Print("Choose the pair selection method bin index from the list above "\
                      "(typing in any text will exit the program)");
      std::cout << ">> ";
      if (!(std::cin >> methodBinIndex))
      {
         CppTools::PrintInfo("Exiting the program");
         exit(1);
      }

      if (methodBinIndex < inputYAMLResonance["pair_selection_methods"].size() && 
          methodBinIndex >= 0) break;
      else CppTools::PrintWarning("Chosen method bin is out of range");
   }

   int centralityBinIndex;
   while (true) // ininite loop until exit or valid input is specified
   {
      CppTools::PrintInfo("List of centrality bins");

      for (int i = 0; i < inputYAMLResonance["centrality_bins"].size(); i++)
      {
         CppTools::Print(i, inputYAMLResonance["centrality_bins"][i]["name"].as<std::string>());
      }

      CppTools::Print("Choose the centrality bin index from the list above "\
                      "(typing in any text will exit the program)");
      std::cout << ">> ";
      if (!(std::cin >> centralityBinIndex))
      {
         CppTools::PrintInfo("Exiting the program");
         exit(1);
      }

      if (centralityBinIndex < inputYAMLResonance["centrality_bins"].size() &&
          centralityBinIndex >= 0) break;
      else CppTools::PrintWarning("Chosen centrality bin is out of range");
   }

   resonanceName = inputYAMLResonance["name"].as<std::string>();
   massResonance = inputYAMLResonance["mass"].as<double>();
   gammaResonance = inputYAMLResonance["gamma"].as<double>();

   daughter1Id = inputYAMLResonance["daughter1_id"].as<int>();
   daughter2Id = inputYAMLResonance["daughter2_id"].as<int>();

   minMInv = inputYAMLResonance["m_inv_range_min"].as<double>();
   maxMInv = inputYAMLResonance["m_inv_range_max"].as<double>();

   SetGaussianBroadeningFunction();

   inputFile = TFile::Open(inputFileName.c_str(), "READ");

   text.SetTextFont(43);
   text.SetTextSize(45);

   texText.SetTextFont(43);
   texText.SetTextSize(45);

   text.SetTextAngle(270.);

   pTNBins = inputYAMLResonance["pt_bins"].size();

   for (unsigned int i = 0; i < pTNBins; i++)
   {
      pTBinRanges.push_back(inputYAMLResonance["pt_bins"][i]["min"].as<double>());
   }
   pTBinRanges.push_back(inputYAMLResonance["pt_bins"][pTNBins - 1]["max"].as<double>());

   numberOfIterations = pTNBins;

   const std::string parametersOutputDir = "data/Parameters/ResonanceBGFit/" + runName;
   std::filesystem::create_directories(parametersOutputDir);

   PerformMInvFits(inputYAMLResonance["pair_selection_methods"][methodBinIndex], 
                   inputYAMLResonance["centrality_bins"][centralityBinIndex]);

   gROOT->SetBatch(false);

   /*
	TCanvas *canv = new TCanvas("", "", 900, 900);
   
   GUIDistrCutter2D::AddHistogram(realHist);
   GUIDistrCutter2D::AddHistogram(static_cast<TH2D *>(simHist->Clone("sim")));
   system(("mkdir -p data/Parameters/Deadmaps/" + runName).c_str());

   while (detectorName.find(" ") < detectorName.size())
   {
      const unsigned int spacePos = detectorName.find(" ");
      detectorName.erase(spacePos, 1);
   }

   const std::string outputCutsFileName = "data/Parameters/Deadmaps/" + 
                                          runName + "/" + detectorName + ".txt";

   if (CppTools::FileExists(outputCutsFileName))
   {
      GUIDistrCutter2D::ReadCutAreas(outputCutsFileName);
   }
   GUIDistrCutter2D::SetOutputFile(outputCutsFileName);

	gPad->AddExec("exec", "GUIDistrCutter2D::Exec()");
   */
}

void AnalyzeRealMInv::PerformMInvFits(const YAML::Node& method, const YAML::Node& centrality)
{
   const std::string methodName = method["name"].as<std::string>();

   const std::string outputDir = "output/MInv/" + runName + "/" + 
                                 std::to_string(taxiNumber) + "/" + methodName;
   std::filesystem::create_directories(outputDir);

   unsigned int pTBinFitMin = method["pt_bin_min"].as<int>();
   unsigned int pTBinFitMax = method["pt_bin_max"].as<int>();

   const std::string centralityName = centrality["name"].as<std::string>();

   pBar.SetText("Preparing M_{inv}");

   for (unsigned int i = 0; i < pTNBins; i++)
   {
      pBar.Print(static_cast<double>(numberOfCalls)/static_cast<double>(numberOfIterations));

      // only bins for which approximation is done needed
      if (i < pTBinFitMin && i > pTBinFitMax) continue;

      TH1D *distrMInvFG = nullptr;
      TH1D *distrMInvBG = nullptr;
      TH1D *distrMInvFGLR = nullptr; // low resolution (for BG normalization)
      TH1D *distrMInvBGLR = nullptr; // low resolution (for BG normalization)

      std::string decayMode = ParticleMap::nameShort[daughter1Id] +
                              ParticleMap::nameShort[daughter2Id];

      double numberOfEvents = 0.;

      TH1D *distrMInv = 
         MInv::Merge(inputFile, methodName, decayMode, 
                     centrality["cb_c_min"].as<int>(), centrality["cb_c_max"].as<int>(),
                     0, inputYAMLResonance["cb_z_bins"].as<int>() - 1, 
                     0, inputYAMLResonance["cb_r_bins"].as<int>() - 1,
                     pTBinRanges[i], pTBinRanges[i + 1],
                     distrMInvFG, distrMInvBG, distrMInvFGLR, distrMInvBGLR, 
                     numberOfEvents);

      if (inputYAMLResonance["has_antiparticle"].as<bool>() && 
          !inputYAMLResonance["separate_antiparticle"].as<bool>())
      {
         decayMode = ParticleMap::nameShort[daughter2Id] +
                     ParticleMap::nameShort[daughter1Id];
         distrMInv->Add(MInv::Merge(inputFile, methodName, decayMode, 
                                    centrality["cb_c_min"].as<int>(), 
                                    centrality["cb_c_max"].as<int>(),
                                    0, inputYAMLResonance["cb_z_bins"].as<int>() - 1, 
                                    0, inputYAMLResonance["cb_r_bins"].as<int>() - 1,
                                    pTBinRanges[i], pTBinRanges[i + 1],
                                    distrMInvFG, distrMInvBG, distrMInvFGLR, distrMInvBGLR,
                                    numberOfEvents));
      }

      if (!distrMInv)
      {
         pBar.Clear();
         CppTools::PrintError("Resulting M_{inv} histogram could not be constructed for " + 
                              CppTools::DtoStr(pTBinRanges[i], 2) + "<pT<" + 
                              CppTools::DtoStr(pTBinRanges[i + 1], 2));
      }
      else if (distrMInv->Integral(1, distrMInv->GetXaxis()->GetNbins()) < 1e-7)
      {
         pBar.Clear();
         CppTools::PrintWarning("Resulting histogram is empty in " + 
                                centrality["name"].as<std::string>() + " " +
                                CppTools::DtoStr(pTBinRanges[i], 2) + "<pT<" + 
                                CppTools::DtoStr(pTBinRanges[i + 1], 2));
         pBar.RePrint();
         continue;
      }

      if (!distrMInvFG)
      {
         pBar.Clear();
         CppTools::PrintError("Resulting M_{inv} foreground histogram "\
                              "could not be constructed for " + 
                              centrality["name"].as<std::string>() + " " +
                              CppTools::DtoStr(pTBinRanges[i], 2) + "<pT<" + 
                              CppTools::DtoStr(pTBinRanges[i + 1], 2));
      }

      if (!distrMInvBG)
      {
         pBar.Clear();
         CppTools::PrintWarning("Resulting M_{inv} foreground histogram "\
                                "could not be constructed for " + 
                                centrality["name"].as<std::string>() + " " +
                                CppTools::DtoStr(pTBinRanges[i], 2) + "<pT<" + 
                                CppTools::DtoStr(pTBinRanges[i + 1], 2));
         pBar.RePrint();
      }
      else if (distrMInvBG->Integral(1, distrMInvBG->GetXaxis()->GetNbins()) < 1e-7)
      {
         pBar.Clear();
         CppTools::PrintWarning("Resulting background histogram is empty in " + 
                                centrality["name"].as<std::string>() + 
                                CppTools::DtoStr(pTBinRanges[i], 2) + "<pT<" + 
                                CppTools::DtoStr(pTBinRanges[i + 1], 2));
         pBar.RePrint();
      }

      if (rebinX != 1)
      {
         distrMInv->Rebin(rebinX);
         distrMInvFG->Rebin(rebinX);
         distrMInvBG->Rebin(rebinX);
      }

      // sigma of a gaus that is convoluted with Breit-Wigner
      const double gaussianBroadeningSigma = 
         gaussianBroadeningEstimatorFunc->Eval((pTBinRanges[i] + pTBinRanges[i + 1])/2.);

      // fit for resonance+bg approximation
      TF1 fit("resonance + bg fit", &FitFunc::RBWConvGausBGPol3, 
              massResonance - gammaResonance*3., massResonance + gammaResonance*3., 8);
      // fit for resonance approximation
      TF1 fitResonance("resonance fit", &FitFunc::RBWConvGaus, massResonance - gammaResonance*3., 
                       massResonance + gammaResonance*3., 4);
      // fit for bg approximation
      TF1 fitBG("bg fit", &FitFunc::Pol3, massResonance - gammaResonance*3., 
                massResonance + gammaResonance*3., 4);

      const double maxBinVal = distrMInv->GetBinContent(distrMInv->GetMaximumBin());
      const double minBinVal = distrMInv->GetBinContent(distrMInv->GetMinimumBin());

      fit.SetParameters(maxBinVal, massResonance, gammaResonance, gaussianBroadeningSigma);

      fit.SetParLimits(0, 1., maxBinVal - minBinVal);
      fit.SetParLimits(1, massResonance/1.05, massResonance*1.05);
      fit.SetParLimits(2, gammaResonance/1.02, gammaResonance*1.05);
      fit.FixParameter(3, gaussianBroadeningSigma);

      distrMInv->Fit(&fit, "RQMNBLC");

      for (unsigned int j = 1; j <= fitNTries; j++)
      {
         fit.SetParLimits(1, fit.GetParameter(1)/(1. + 0.05/static_cast<double>(j*j)), 
                          fit.GetParameter(1)*(1. + 0.05/static_cast<double>(j*j)));
         fit.SetParLimits(2, fit.GetParameter(2)/(1. + 0.05/static_cast<double>(j*j)),
                          fit.GetParameter(2)*(1. + 0.1/static_cast<double>(j*j)));

         fit.SetRange(fit.GetParameter(1) - (fit.GetParameter(2) + 
                                             gaussianBroadeningSigma)*3., 
                      fit.GetParameter(1) + (fit.GetParameter(2) + 
                                             gaussianBroadeningSigma)*3.);

         distrMInv->Fit(&fit, "RQMNBLC");
         //if (j == fitNTries) distrMInv->Fit(&fit, "RQMNBLE");
      }

      for (int j = 0; j < fitResonance.GetNpar(); j++)
      {
         fitResonance.SetParameter(j, fit.GetParameter(j));
      }

      for (int j = 0; j < fitBG.GetNpar(); j++)
      {
         fitBG.SetParameter(j, fit.GetParameter(j + fitResonance.GetNpar()));
      }

      fitResonance.SetRange(fit.GetParameter(1) - 
                            (fit.GetParameter(2) + gaussianBroadeningSigma)*3., 
                            fit.GetParameter(1) + 
                            (fit.GetParameter(2) + gaussianBroadeningSigma)*3.);

      fitBG.SetRange(fit.GetParameter(1) - (fit.GetParameter(2) + gaussianBroadeningSigma)*3., 
                     fit.GetParameter(1) + (fit.GetParameter(2) + gaussianBroadeningSigma)*3.);

      fit.SetLineWidth(4);
      fitBG.SetLineWidth(4);

      fit.SetLineColorAlpha(kRed - 3, 0.8);
      fitBG.SetLineColorAlpha(kGray + 2, 0.8);

      fitBG.SetLineStyle(7);

      distrMInv->SetMaximum(distrMInv->GetMaximum()*1.2);

      distrMInv->GetXaxis()->SetRange(distrMInv->GetXaxis()->FindBin(minMInv + 1e-7), 
                                      distrMInv->GetXaxis()->FindBin(maxMInv - 1e-7));
      distrMInvFG->GetXaxis()->SetRange(distrMInvFG->GetXaxis()->FindBin(minMInv + 1e-7), 
                                        distrMInvFG->GetXaxis()->FindBin(maxMInv - 1e-7));
      distrMInvBG->GetXaxis()->SetRange(distrMInvBG->GetXaxis()->FindBin(minMInv + 1e-7), 
                                        distrMInvBG->GetXaxis()->FindBin(maxMInv - 1e-7));

      distrMInv->SetLineWidth(2);
      distrMInvFG->SetLineWidth(3);
      distrMInvBG->SetLineWidth(2);

      distrMInv->SetLineColor(kGray + 3);
      distrMInv->SetMarkerColor(kGray + 3);

      distrMInv->SetMarkerStyle(20);
      distrMInv->SetMarkerSize(0.7);
      distrMInv->SetMarkerColor(kGray + 3);

      { /* canvas with invariant mass distribution with subtracted background only */
         TCanvas canvMInv("canv MInv", "", 800, 800);

         gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05); 
         gPad->SetLeftMargin(0.173); gPad->SetBottomMargin(0.112);

         ROOTTools::DrawFrame(distrMInv, "", "#it{M}_{inv} [GeV/#it{c}^{2}]", "Counts", 1., 1.9);

         text.DrawTextNDC(0.9, 0.93, (methodName).c_str());
         texText.DrawLatexNDC(0.2, 0.88, (CppTools::DtoStr(pTBinRanges[i], 1) + 
                              " < #it{p}_{T} < " + 
                              CppTools::DtoStr(pTBinRanges[i + 1], 1)).c_str());
         text.DrawTextNDC(0.85, 0.93, centrality["name_tex"].as<std::string>().c_str());
         /*
         texText.DrawLatexNDC(0.2, 0.81, 
                              ("#it{#chi}^{2}/NDF=" + 
                               CppTools::DtoStr(fit.GetChisquare()/fit.GetNDF(), 1)).c_str());
                               */
         fitBG.Draw("SAME");
         fit.Draw("SAME");

         ROOTTools::PrintCanvas(&canvMInv, outputDir + "/" + resonanceName + "_" + 
                                centrality["name"].as<std::string>() + "_" +
                                CppTools::DtoStr(pTBinRanges[i], 1) + "-" + 
                                CppTools::DtoStr(pTBinRanges[i + 1], 1));
      } /* canvas with invariant mass distribution with subtracted background only */

      numberOfCalls++;
   }
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
