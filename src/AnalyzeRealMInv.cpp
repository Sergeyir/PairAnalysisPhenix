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
   if (argc < 3 || argc > 6) 
   {
      CppTools::PrintError("Expected 2-3 parameters while " + std::to_string(argc - 1) + " "\
                           "parameter(s) were provided \n Usage: bin/AnalyzeRealMInv "\
                           "inputYAMLName taxiNumber methodName=all rebinX=1 "\
                           "numberOfThreads=std::thread::hardware_concurrency()");
   }
 
   CppTools::CheckInputFile(argv[1]);

   taxiNumber = std::stoi(argv[2]);
 
   if (argc == 6) ROOT::EnableImplicitMT(std::stoi(argv[5]));
   else ROOT::EnableImplicitMT(std::thread::hardware_concurrency());

   if (argc == 5) 
   {
      rebinX = std::stoi(argv[4]);
      if (rebinX <= 0) CppTools::PrintError("Rebin value cannot be 0 or negative");
   }

   std::string methodToAnalyze;
   if (argc >= 4) methodToAnalyze = argv[3];
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

   minMInv = inputYAMLResonance["m_inv_range_min"].as<double>();
   maxMInv = inputYAMLResonance["m_inv_range_max"].as<double>();

   SetGaussianBroadeningFunction();

   inputFileName = "data/Real/" + runName + "/Resonance/" + std::to_string(taxiNumber) + ".root";

   CppTools::CheckInputFile(inputFileName);
   inputFile = TFile::Open(inputFileName.c_str(), "READ");

   const std::string inputFileRecEffName = "data/Parameters/ResonanceEff/" + runName + 
                                           "/" + resonanceName + ".root";

   if (CppTools::FileExists(inputFileRecEffName))
   {
      inputFileRecEff = TFile::Open(inputFileRecEffName.c_str(), "READ");
   }
   else
   {
      CppTools::PrintWarning("File " + inputFileRecEffName + " does not exists; spectra will "\
                             "not be normalized by the reconstruction efficiency");
   }

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

   gStyle->SetPalette(kRainbow);

   int methodIndex = -1;
   if (methodToAnalyze == "all")
   {
      for (const YAML::Node& method : inputYAMLResonance["pair_selection_methods"])
      {
         numberOfIterations += pTNBins*inputYAMLResonance["centrality_bins"].size();
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

      numberOfIterations = pTNBins*inputYAMLResonance["centrality_bins"].size();
   }

   const std::string parametersOutputDir = "data/RawYields/Resonance/" + runName;
   std::filesystem::create_directories(parametersOutputDir);

   parametersOutputFile = TFile::Open((parametersOutputDir + "/" + std::to_string(taxiNumber) + 
                                       "_" + resonanceName + ".root").c_str(), "RECREATE");

   // performing fits for specified pair selection methods
   if (methodToAnalyze == "all")
   {
      for (const YAML::Node& method : inputYAMLResonance["pair_selection_methods"])
      {
         PerformMInvFits(method);
      }
   }
   else 
   {
      PerformMInvFits(inputYAMLResonance["pair_selection_methods"][methodIndex]);
   }

   parametersOutputFile->Close();

   pBar.Finish();

   CppTools::PrintInfo("AnalyzeRealMInv executable has finished running succesfully");
}

void AnalyzeRealMInv::PerformMInvFits(const YAML::Node& method)
{
   const std::string methodName = method["name"].as<std::string>();

   const std::string outputDir = "output/MInv/" + runName + "/" + 
                                 std::to_string(taxiNumber) + "/" + methodName;
   std::filesystem::create_directories(outputDir);

   unsigned int pTBinFitMin = method["pt_bin_min"].as<int>();
   unsigned int pTBinFitMax = method["pt_bin_max"].as<int>();

   parametersOutputFile->mkdir(methodName.c_str());
   parametersOutputFile->cd(methodName.c_str());

   TH1D* distrRecEffVsPT = nullptr;

   if (inputFileRecEff)
   {
      const std::string distrRecEffVsPTName = methodName + "/reconstruction efficiency vs pT";
      distrRecEffVsPT = 
         static_cast<TH1D *>(inputFileRecEff->Get(distrRecEffVsPTName.c_str()));

      if (!distrRecEffVsPT)
      {
         CppTools::PrintWarning("Histogram named \"" +  distrRecEffVsPTName + 
                                "\" does not exist in file" + 
                                static_cast<std::string>(inputFileRecEff->GetName()) + 
                                "; correction on reconstruction efficiency will not be applied");
      }
   }

   for (const YAML::Node &centrality : inputYAMLResonance["centrality_bins"])
   {
      const std::string centralityName = centrality["name"].as<std::string>();

      parametersOutputFile->mkdir((methodName + "/" + centralityName).c_str());
      parametersOutputFile->cd((methodName + "/" + centralityName).c_str());

      TH1D distrMeansVsPT("means vs pT", "", pTNBins, &pTBinRanges[0]);
      TH1D distrGammasVsPT("gammas vs pT", "", pTNBins, &pTBinRanges[0]);
      TH1D distrRawYieldVsPT("raw yield vs pT", "", pTNBins, &pTBinRanges[0]);
      TH1D distrSpectraVsPT("spectra vs pT", "", pTNBins, &pTBinRanges[0]);

      for (unsigned int i = 0; i < pTNBins; i++)
      {
         pBar.Print(static_cast<double>(numberOfCalls)/static_cast<double>(numberOfIterations));

         TH1D *distrMInvFG = nullptr;
         TH1D *distrMInvBG = nullptr;
         TH1D *distrMInvFGLR = nullptr; // low resolution (for BG normalization)
         TH1D *distrMInvBGLR = nullptr; // low resolution (for BG normalization)

         // shows whether to perform approximations for the current bin
         // if false the histograms will not be approximated but will be printed anyways
         const bool performFit = (i >= pTBinFitMin && i <= pTBinFitMax);

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
                                 centrality["name"].as<std::string>() + " " +
                                 CppTools::DtoStr(pTBinRanges[i], 2) + "<pT<" + 
                                 CppTools::DtoStr(pTBinRanges[i + 1], 2));
         }
         else if (performFit && distrMInv->Integral(1, distrMInv->GetXaxis()->GetNbins()) < 1e-7)
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
         else if (performFit && distrMInvBG->Integral(1, distrMInvBG->GetXaxis()->GetNbins()) < 1e-7)
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

         if (performFit)
         {
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

            double rawYieldErr;
            const double rawYield = 
               GetYield(distrMInv, fitBG, fit.GetParameter(1) - fit.GetParameter(2)*2., 
                        fit.GetParameter(1) + fit.GetParameter(2)*2., rawYieldErr);

            distrRawYieldVsPT.SetBinContent(i + 1, rawYield/(pTBinRanges[i + 1] - pTBinRanges[i]));
            distrRawYieldVsPT.SetBinError(i + 1, rawYieldErr/rawYield/
                                          (pTBinRanges[i + 1] - pTBinRanges[i]));

            const double rawYieldNorm = 2.*M_PI*(pTBinRanges[i] + pTBinRanges[i + 1])/2.*
                                        (pTBinRanges[i + 1] - pTBinRanges[i])*numberOfEvents;

            distrSpectraVsPT.SetBinContent(i + 1, rawYield/rawYieldNorm);
            distrSpectraVsPT.SetBinError(i + 1, rawYieldErr/rawYield/rawYieldNorm);

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
         }

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
         distrMInvFGLR->SetLineWidth(3);
         distrMInvBGLR->SetLineWidth(2);

         distrMInv->SetLineColor(kGray + 3);
         distrMInv->SetMarkerColor(kGray + 3);

         distrMInv->SetMarkerStyle(20);
         distrMInv->SetMarkerSize(0.7);
         distrMInv->SetMarkerColor(kGray + 3);

         int lastNonZeroBinLR = 1;
         for (int i = distrMInvFGLR->GetXaxis()->GetNbins(); i >= 1; i--)
         {
            if (distrMInvFGLR->GetBinContent(i) > 1e-7 || distrMInvBGLR->GetBinContent(i) > 1e-7) 
            {
               lastNonZeroBinLR = i;
               break;
            }
         }

         distrMInvFGLR->GetXaxis()->SetRange(1, lastNonZeroBinLR + 2);
         distrMInvBGLR->GetXaxis()->SetRange(1, lastNonZeroBinLR + 2);

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
            if (performFit)
            {
               /*
               texText.DrawLatexNDC(0.2, 0.81, 
                                    ("#it{#chi}^{2}/NDF=" + 
                                     CppTools::DtoStr(fit.GetChisquare()/fit.GetNDF(), 1)).c_str());
                                     */
               fitBG.Draw("SAME");
               fit.Draw("SAME");
            }

            ROOTTools::PrintCanvas(&canvMInv, outputDir + "/" + resonanceName + "_" + 
                                   centrality["name"].as<std::string>() + "_" +
                                   CppTools::DtoStr(pTBinRanges[i], 1) + "-" + 
                                   CppTools::DtoStr(pTBinRanges[i + 1], 1));
         } /* canvas with invariant mass distribution with subtracted background only */


         { /* summary canvas: MInv, FGBG, FG/BG */
            TCanvas canvMInvSummary("canv MInv summary", "", 1920, 1080);

            distrMInv->SetFillStyle(3003);
            distrMInvFG->SetFillStyle(3005);
            distrMInvBG->SetFillStyle(3004);
            distrMInvFGLR->SetFillStyle(3005);
            distrMInvBGLR->SetFillStyle(3004);

            canvMInvSummary.Divide(3, 2);
            
            canvMInvSummary.cd(1);

            gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05); 
            gPad->SetLeftMargin(0.173); gPad->SetBottomMargin(0.112);

            ROOTTools::DrawFrame(static_cast<TH1D *>(distrMInv->Clone()), 
                                 "", "#it{M}_{inv} [GeV/#it{c}^{2}]", "Counts", 1., 1.7);

            texText.DrawLatexNDC(0.2, 0.85, (CppTools::DtoStr(pTBinRanges[i], 1) + 
                                 " < #it{p}_{T} < " + 
                                 CppTools::DtoStr(pTBinRanges[i + 1], 1)).c_str());
            text.DrawTextNDC(0.8, 0.92, centrality["name_tex"].as<std::string>().c_str());
            text.DrawTextNDC(0.88, 0.92, (methodName).c_str());

            if (performFit)
            {
               /*
               texText.DrawLatexNDC(0.2, 0.74, 
                                    ("#it{#chi}^{2}/NDF = " + 
                                     CppTools::DtoStr(fit.GetChisquare()/fit.GetNDF(), 1)).c_str());
                                     */
               /*

               texText.DrawLatexNDC(0.6, 0.9, ("#it{#mu}=" + CppTools::DtoStr(fit.GetParameter(1)*1000., 1) + 
                                               " [MeV/c^{2}]").c_str());
               texText.DrawLatexNDC(0.6, 0.85, ("#it{#Gamma}=" + 
                                                CppTools::DtoStr(fit.GetParameter(2)*1000., 1) + 
                                                " [MeV/c^{2}]").c_str());
                                                */

               fitBG.Draw("SAME");
               fit.Draw("SAME");
            }

            canvMInvSummary.cd(4);

            gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05); 
            gPad->SetLeftMargin(0.173); gPad->SetBottomMargin(0.112);

            ROOTTools::DrawFrame(static_cast<TH1D *>(distrMInv->Clone()), 
                                 "", "#it{M}_{inv} [GeV/#it{c}^{2}]", "Counts", 1., 1.7);

            canvMInvSummary.cd(5);

            gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05); 
            gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.112);

            if (distrMInvBG->GetEntries() > 1e-3 && distrMInvFG->GetEntries() > 1e-3)
            {
               TH1D *ratioFGBG = static_cast<TH1D *>(distrMInvFG->Clone());
               ratioFGBG->Divide(distrMInvBG);

               ratioFGBG->SetLineWidth(2);
               ratioFGBG->SetLineColor(kBlack);
               ratioFGBG->SetMarkerColor(kBlack);

               ROOTTools::DrawFrame(ratioFGBG, "", "#it{M}_{inv} [GeV/#it{c}^{2}]", 
                                    "FG/BG", 1., 1.7);
            }
            else
            {
               text.DrawTextNDC(0.88, 0.92, "No data");
            }

            canvMInvSummary.cd(6);

            gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05); 
            gPad->SetLeftMargin(0.148); gPad->SetBottomMargin(0.112);

            if (distrMInvBGLR->GetEntries() > 1e-3 && distrMInvFGLR->GetEntries() > 1e-3)
            {
               TH1D *ratioFGBGLR = static_cast<TH1D *>(distrMInvFGLR->Clone());
               ratioFGBGLR->Divide(distrMInvBGLR);

               ratioFGBGLR->SetLineWidth(2);
               ratioFGBGLR->SetLineColor(kBlack);
               ratioFGBGLR->SetMarkerColor(kBlack);

               ROOTTools::DrawFrame(ratioFGBGLR, "", "#it{M}_{inv} [GeV/#it{c}^{2}]", 
                                    "FG/BG", 1., 1.7);
            }
            else
            {
               text.DrawTextNDC(0.88, 0.92, "No data");
            }

            distrMInvFG->SetMaximum(distrMInvFG->GetMaximum()*1.15);
            distrMInvFGLR->SetMaximum(distrMInvFGLR->GetMaximum()*3.);

            canvMInvSummary.cd(2);

            gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05); 
            gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.112);

            distrMInv->Sumw2(false);
            distrMInvBG->Sumw2(false);
            distrMInvFG->Sumw2(false);

            if (distrMInvFG->GetEntries() > distrMInv->GetEntries())
            {
               distrMInvFG->SetMinimum(0.);
               ROOTTools::DrawFrame(distrMInvFG, "", "#it{M}_{inv} [GeV/#it{c}^{2}]", 
                                    "Counts", 1., 1.7, 0.05, 0.05, true, false);
            }
            else
            {
               ROOTTools::DrawFrame(distrMInv, "", "#it{M}_{inv} [GeV/#it{c}^{2}]", 
                                    "Counts", 1., 1.7, 0.05, 0.05, true, false);
            }

            if (distrMInvFG->GetEntries() > 1e-3) 
            {
               distrMInvFG->SetLineColorAlpha(kAzure + 2, 0.8);
               distrMInvFG->Draw("SAME PFC");
            }
            else text.DrawTextNDC(0.88, 0.9, "No data on foreground");

            if (distrMInvBG->GetEntries() > 1e-3) 
            {
               distrMInvBG->SetLineColorAlpha(kGreen + 2, 0.8);
               distrMInvBG->Draw("SAME PFC");
            }
            else text.DrawTextNDC(0.78, 0.9, "No data on background");

            distrMInv->SetLineColorAlpha(kRed + 2, 0.8);
            distrMInv->Draw("SAME PFC");

            canvMInvSummary.cd(3);

            gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05); 
            gPad->SetLeftMargin(0.148); gPad->SetBottomMargin(0.112);

            gPad->SetLogy();

            distrMInvBGLR->Sumw2(false);

            if (distrMInvFGLR->GetEntries() > 1e-3)
            {
               ROOTTools::DrawFrame(distrMInvFGLR, "", "#it{M}_{inv} [GeV/#it{c}^{2}]", 
                                    "Counts", 1., 1.7, 0.05, 0.05, true, false);
               distrMInvFGLR->SetLineColorAlpha(kAzure + 2, 0.8);
               distrMInvFGLR->Draw("SAME PFC");

               if (distrMInvBGLR->GetEntries() > 1e-3) 
               {
                  distrMInvBGLR->SetLineColorAlpha(kRed + 2, 0.8);
                  distrMInvBGLR->Draw("SAME PFC");
               }
               else text.DrawTextNDC(0.78, 0.9, "No data on background");
            }
            else text.DrawTextNDC(0.88, 0.9, "No data on foreground");

            ROOTTools::PrintCanvas(&canvMInvSummary, outputDir + "/Summary_" + resonanceName + "_" + 
                                   centrality["name"].as<std::string>() + "_" +
                                   CppTools::DtoStr(pTBinRanges[i], 1) + "-" + 
                                   CppTools::DtoStr(pTBinRanges[i + 1], 1), true, false);
         } /* summary canvas: MInv, FGBG, FG/BG */

         { /* FG, BG, and signal on the same canvas */
            TCanvas canvMInvFGBG("canv MInv FGBG", "", 800, 800);

            gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05); 
            gPad->SetLeftMargin(0.173); gPad->SetBottomMargin(0.112);

            distrMInv->Sumw2(false);
            distrMInvBG->Sumw2(false);
            distrMInvFG->Sumw2(false);

            if (distrMInvFG->GetEntries() > distrMInv->GetEntries())
            {
               distrMInvFG->SetMinimum(0.);
               ROOTTools::DrawFrame(distrMInvFG, "", "#it{M}_{inv} [GeV/#it{c}^{2}]", 
                                    "Counts", 1., 1.7, 0.05, 0.05, true, false);
            }
            else
            {
               ROOTTools::DrawFrame(distrMInv, "", "#it{M}_{inv} [GeV/#it{c}^{2}]", 
                                    "Counts", 1., 1.7, 0.05, 0.05, true, false);
            }

            if (distrMInvFG->GetEntries() > 1e-3) 
            {
               distrMInvFG->SetLineColorAlpha(kAzure + 2, 0.8);
               distrMInvFG->Draw("SAME PFC");
            }
            else text.DrawTextNDC(0.88, 0.9, "No data on foreground");

            if (distrMInvBG->GetEntries() > 1e-3) 
            {
               distrMInvBG->SetLineColorAlpha(kGreen + 2, 0.8);
               distrMInvBG->Draw("SAME PFC");
            }
            else text.DrawTextNDC(0.78, 0.9, "No data on background");

            distrMInv->SetLineColorAlpha(kRed + 2, 0.8);
            distrMInv->Draw("SAME PFC");

            TLegend legend(0.3, 0.86, 0.98, 0.94);
            legend.SetLineColorAlpha(0, 0.);
            legend.SetFillColorAlpha(0, 0.);
            legend.SetNColumns(3);

            legend.AddEntry(distrMInvFG, "FG", "PFC");
            legend.AddEntry(distrMInvBG, "BG", "PFC");
            legend.AddEntry(distrMInv, "FG-BG", "PFC");

            legend.Draw();

            ROOTTools::PrintCanvas(&canvMInvFGBG, outputDir + "/FGBG_" + resonanceName + "_" + 
                                   centrality["name"].as<std::string>() + "_" +
                                   CppTools::DtoStr(pTBinRanges[i], 1) + "-" + 
                                   CppTools::DtoStr(pTBinRanges[i + 1], 1), false);
         } /* FG, BG, and signal on the same canvas */

         numberOfCalls++;
      }

      if (distrRecEffVsPT)
      {
         distrSpectraVsPT.Divide(distrRecEffVsPT);
      }

      distrMeansVsPT.SetLineColor(kRed - 2);
      distrMeansVsPT.SetMarkerColor(kRed - 2);
      distrMeansVsPT.SetLineWidth(4);

      distrGammasVsPT.SetLineColor(kRed - 2);
      distrGammasVsPT.SetMarkerColor(kRed - 2);
      distrGammasVsPT.SetLineWidth(4);

      distrRawYieldVsPT.SetLineColor(kRed - 2);
      distrRawYieldVsPT.SetMarkerColor(kRed - 2);
      distrRawYieldVsPT.SetLineWidth(4);

      distrSpectraVsPT.SetLineColor(kRed - 2);
      distrSpectraVsPT.SetMarkerColor(kRed - 2);
      distrSpectraVsPT.SetLineWidth(4);

      TLine massResonancePDG(pTBinRanges[0], massResonance, pTBinRanges[pTNBins], massResonance);
      massResonancePDG.SetLineColorAlpha(kBlack, 0.5);
      massResonancePDG.SetLineStyle(2);
      massResonancePDG.SetLineWidth(4);

      TLine gammaResonancePDG(pTBinRanges[0], gammaResonance, pTBinRanges[pTNBins], gammaResonance);
      gammaResonancePDG.SetLineColorAlpha(kBlack, 0.5);
      gammaResonancePDG.SetLineStyle(2);
      gammaResonancePDG.SetLineWidth(4);

      distrMeansVsPT.SetMaximum(massResonance*1.05);
      distrMeansVsPT.SetMinimum(massResonance*0.95);

      distrGammasVsPT.SetMaximum(gammaResonance*1.5);
      distrGammasVsPT.SetMinimum(gammaResonance/2.);

      TCanvas canvMeansVsPT("canv means vs pT", "", 800, 800);

      gPad->SetRightMargin(0.03);
      gPad->SetTopMargin(0.02);
      gPad->SetLeftMargin(0.172);
      gPad->SetBottomMargin(0.112);

      ROOTTools::DrawFrame(&distrMeansVsPT, "", "#it{p}_{T} [GeV/#it{c}]", 
                           "#it{#mu} [GeV/#it{c}^{2}]", 1., 1.82);

      text.DrawTextNDC(0.9, 0.95, (methodName).c_str());
      massResonancePDG.Draw();

      ROOTTools::PrintCanvas(&canvMeansVsPT, outputDir + "/" + resonanceName + 
                             "_means_" + centralityName);

      TCanvas canvGammasVsPT("canv gammas vs pT", "", 800, 800);

      gPad->SetRightMargin(0.03);
      gPad->SetTopMargin(0.02);
      gPad->SetLeftMargin(0.172);
      gPad->SetBottomMargin(0.112);

      ROOTTools::DrawFrame(&distrGammasVsPT, "", "#it{p}_{T} [GeV/#it{c}]", 
                           "#it{#Gamma} [GeV/#it{c}^{2}]", 1., 1.82);

      text.DrawTextNDC(0.9, 0.95, (methodName).c_str());
      gammaResonancePDG.Draw();

      ROOTTools::PrintCanvas(&canvGammasVsPT, outputDir + "/" + resonanceName + 
                             "_gammas_" + centralityName);

      TCanvas canvRawYieldVsPT("canv raw yield vs pT", "", 800, 800);

      gPad->SetLogy();

      gPad->SetRightMargin(0.03);
      gPad->SetTopMargin(0.02);
      gPad->SetLeftMargin(0.141);
      gPad->SetBottomMargin(0.112);

      ROOTTools::DrawFrame(&distrRawYieldVsPT, "", "#it{p}_{T} [GeV/#it{c}]", 
                           "#it{dY}_{raw}/#it{dp}_{T} [(GeV/#it{c})^{-1}]", 1., 1.35);

      text.DrawTextNDC(0.9, 0.95, (methodName).c_str());

      ROOTTools::PrintCanvas(&canvRawYieldVsPT, outputDir + "/" + resonanceName + 
                             "_raw_yield_" + centralityName);

      TCanvas canvSpectraVsPT("canv spectra vs pT", "", 800, 800);

      gPad->SetLogy();

      gPad->SetRightMargin(0.03);
      gPad->SetTopMargin(0.02);
      gPad->SetLeftMargin(0.141);
      gPad->SetBottomMargin(0.112);

      if (distrRecEffVsPT)
      {
         ROOTTools::DrawFrame(&distrSpectraVsPT, "", "#it{p}_{T} [GeV/#it{c}]", 
                              "1/(2 #pi #it{p}_{T}) #it{d}^{2}#it{N}/#it{dp}_{T}#it{dy} "\
                              "[(GeV/#it{c})^{-2})]", 1., 1.35);
      }
      else
      {
         ROOTTools::DrawFrame(&distrSpectraVsPT, "", "#it{p}_{T} [GeV/#it{c}]", 
                              "1/(2 #pi #it{p}_{T}) #it{d}^{2}#it{Y}_{raw}/#it{dp}_{T}#it{dy} "\
                              "[(GeV/#it{c})^{-2})]", 1., 1.35);
      }

      text.DrawTextNDC(0.9, 0.95, (methodName).c_str());

      ROOTTools::PrintCanvas(&canvSpectraVsPT, outputDir + "/" + resonanceName + 
                             "_spectra_" + centralityName);

      distrMeansVsPT.Write();
      distrGammasVsPT.Write();
      distrRawYieldVsPT.Write();
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
