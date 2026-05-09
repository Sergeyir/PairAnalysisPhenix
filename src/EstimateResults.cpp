/** 
 *  @file   EstimateResults.hpp 
 *  @brief  Contains realisations of functions and variables that are used for estimation of invariant pT spectra and nuclear modification factors R_{AB} and R_{CP}
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysis).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef ESTIMATE_RESULTS_CPP
#define ESTIMATE_RESULTS_CPP

#include "EstimateResults.hpp"

using namespace EstimateResults;

int main(int argc, char **argv)
{
   if (argc < 3) 
   {
      CppTools::PrintError("Expected 2 parameters while " + std::to_string(argc - 1) + " "\
                           "parameter(s) were provided \n Usage: bin/EstimateResults "\
                           "inputYAMLName taxiNumber");
   }


   gStyle->SetPalette(kSouthWest);

   CppTools::CheckInputFile(argv[1]);

   taxiNumber = std::stoi(argv[2]);

   gStyle->SetOptStat(0);
   gErrorIgnoreLevel = kWarning;

   gROOT->SetBatch(kTRUE);
   gStyle->SetOptStat(kFALSE);

   TH1::AddDirectory(kFALSE);

   ROOT::EnableImplicitMT(std::thread::hardware_concurrency());

   inputYAMLResonance.OpenFile(argv[1]);
   inputYAMLResonance.CheckStatus("resonance");

   runName = inputYAMLResonance["run_name"].as<std::string>();

   inputYAMLMain.OpenFile("input/" + runName + "/main.yaml");
   inputYAMLMain.CheckStatus("main");

   const std::string resonanceName = inputYAMLResonance["name"].as<std::string>();
   const double resonanceMass = inputYAMLResonance["mass"].as<double>();

   inputRecEffFileName = "data/Parameters/RecEffResonance/" + 
                         runName + "/" + resonanceName + ".root";
   CppTools::CheckInputFile(inputRecEffFileName);
   inputRecEffFile = TFile::Open(inputRecEffFileName.c_str());

   pTNBins = inputYAMLResonance["pt_bins"].size();

   for (unsigned int i = 0; i < pTNBins; i++)
   {
      pTBinRanges.push_back(inputYAMLResonance["pt_bins"][i]["min"].as<double>());
   }
   pTBinRanges.push_back(inputYAMLResonance["pt_bins"][pTNBins - 1]["max"].as<double>());

   const std::string resultsOutputDir = "data/Results/" + runName;
   std::filesystem::create_directories(resultsOutputDir);

   const std::string outputDir = "output/Results/" + runName + "/" + std::to_string(taxiNumber);
   std::filesystem::create_directories(outputDir);

   double spectraNorm = inputYAMLResonance["branching_ratio"].as<double>();

   if (inputYAMLResonance["has_antiparticle"] && 
       inputYAMLResonance["separate_antiparticle"]) spectraNorm *= 2.;

   const std::string resultsOutputFileName = resultsOutputDir + "/" + std::to_string(taxiNumber) + 
                                             "_" + resonanceName + ".root";

   TFile *resultsOutputFile = TFile::Open(resultsOutputFileName.c_str(), "RECREATE");

   estimateFactors = !(inputYAMLMain["is_pp"].as<bool>());

   TH1D *distrSpectraPPVsPTStatErr = nullptr;
   TH1D *distrSpectraPPVsPTSysErr = nullptr;
   const std::string spectraPPFileName = "data/Spectra/pp200/" + resonanceName + ".root";

   if (estimateFactors)
   {
      if (!std::filesystem::exists(spectraPPFileName))
      {
         CppTools::PrintWarning("File " + spectraPPFileName + " does not exits; nuclear "\
                                "modification factors will not be evaluated");
         estimateFactors = false;
      }
      else
      {
         TFile *spectraPPFile = TFile::Open(spectraPPFileName.c_str());
         distrSpectraPPVsPTStatErr = 
            static_cast<TH1D *>(spectraPPFile->Get("spectra vs pT with stat err"));
         distrSpectraPPVsPTSysErr = 
            static_cast<TH1D *>(spectraPPFile->Get("spectra vs pT with sys err"));
         if (!distrSpectraPPVsPTStatErr)
         {
            CppTools::PrintError("No p+p spectra with statistical uncertainty was "\
                                 "found in file " + spectraPPFileName);
         }
         if (!distrSpectraPPVsPTSysErr)
         {
            CppTools::PrintError("No p+p spectra with systematic uncertainty was "\
                                 "found in file " + spectraPPFileName);
         }
      }
   }

   for (const YAML::Node& centralityBin : inputYAMLResonance["centrality_bins"])
   {
      const std::string centralityName = centralityBin["name"].as<std::string>();
      const std::string centralityNameTex = centralityBin["name_tex"].as<std::string>();

      std::vector<TH1D *> spectrasVsPTStatErr;
      std::vector<TH1D *> spectrasVsPTSysErr;

      double yMin = 1e31;
      double yMax = 1e-31;

      TLegend legend(0.6, 0.7, 0.95, 0.95);
      legend.SetLineColorAlpha(0, 0.);
      legend.SetFillColorAlpha(0, 0.);

      std::vector<int> methodColors;

      for (const YAML::Node& method : inputYAMLResonance["pair_selection_methods"])
      {
         const std::string methodName = method["name"].as<std::string>();
         methodColors.emplace_back(TColor::GetColor(method["color"].as<std::string>().c_str()));

         inputFileName = "data/RawYields/" + runName + "/Resonance/" + 
                         std::to_string(taxiNumber) + "_" + resonanceName + 
                         "_" + methodName + ".root";
         CppTools::CheckInputFile(inputFileName);
         inputFile = TFile::Open(inputFileName.c_str());

         TH1D *rawYieldVsPTStatErr = 
            static_cast<TH1D *>(inputFile->Get((centralityName + 
                                                "/raw yield vs pT with stat errors").c_str()));
         TH1D *rawYieldVsPTSysErr = 
            static_cast<TH1D *>(inputFile->Get((centralityName + 
                                                "/raw yield vs pT with sys errors").c_str()));

         if (!rawYieldVsPTStatErr) 
         {
            CppTools::PrintError("No raw yield distribution with statistical errors was "\
                                 "found in file " + inputFileName + " for " + methodName + 
                                 " in centrality " + centralityName);
         }
         if (!rawYieldVsPTSysErr) 
         {
            CppTools::PrintError("No raw yield distribution with systematic errors "\
                                 "was found in file " + inputFileName + " for " + methodName + 
                                 " in centrality " + centralityName);
         }

         TH1D *recEffVsPTStatErr = static_cast<TH1D *>
            (inputRecEffFile->Get((methodName + "/reconstruction efficiency "\
                                   "vs pT with stat errors").c_str()));
         TH1D *recEffVsPTSysErr = static_cast<TH1D *>
            (inputRecEffFile->Get((methodName + "/reconstruction efficiency "\
                                   "vs pT with stat errors").c_str()));

         if (!recEffVsPTStatErr)
         {
            CppTools::PrintError("No reconstruction efficiency with statistical errors was "\
                                 "found in file " + inputRecEffFileName + " for " + methodName);
         }
         if (!recEffVsPTSysErr)
         {
            CppTools::PrintError("No reconstruction efficiency with systematic errors was "\
                                 "found in file " + inputRecEffFileName + " for " + methodName);
         }

         rawYieldVsPTStatErr->Divide(recEffVsPTStatErr);
         rawYieldVsPTSysErr->Divide(recEffVsPTSysErr);

         rawYieldVsPTStatErr->Scale(1./spectraNorm);
         rawYieldVsPTSysErr->Scale(1./spectraNorm);

         spectrasVsPTStatErr.emplace_back(rawYieldVsPTStatErr);
         spectrasVsPTSysErr.emplace_back(rawYieldVsPTSysErr);

         for (unsigned int i = 1; i <= pTNBins; i++)
         {
            if (rawYieldVsPTStatErr->GetBinContent(i) < 1e-31) continue;

            yMin = CppTools::Minimum(yMin, rawYieldVsPTStatErr->GetBinContent(i));
            yMax = CppTools::Maximum(yMax, rawYieldVsPTStatErr->GetBinContent(i));
         }
      }

      resultsOutputFile->mkdir(centralityName.c_str());
      resultsOutputFile->cd(centralityName.c_str());

      TH1D distrResultingSpectraVsPTStatErr("spectra vs pT with stat errors", "", 
                                            pTNBins, &pTBinRanges[0]);
      TH1D distrResultingSpectraVsPTSysErr("spectra vs pT with sys errors", "", 
                                           pTNBins, &pTBinRanges[0]);

      distrResultingSpectraVsPTStatErr.SetLineColor(kBlack);
      distrResultingSpectraVsPTStatErr.SetLineWidth(4);

      distrResultingSpectraVsPTSysErr.SetFillColorAlpha(kBlack, 0.3);
      distrResultingSpectraVsPTSysErr.SetFillStyle(1001);

      // graph containing all points across different methods for better Tsallis fit
      TGraphErrors graphSpectraVsPTForTsallisFit;

      for (unsigned int i = 1; i <= pTNBins; i++)
      {
         double minErr = 1e31;
         double minStatErr = 1e31;
         double minSysErr = 1e31;

         double valueOfMinErrHist = -1.;
         for (unsigned int j = 0; j < spectrasVsPTStatErr.size(); j++)
         {
            if (spectrasVsPTStatErr[j]->GetBinContent(i) < 1e-31) continue;

            const double value = spectrasVsPTStatErr[j]->GetBinContent(i);
            const double statErr = spectrasVsPTStatErr[j]->GetBinError(i);
            const double sysErr = spectrasVsPTSysErr[j]->GetBinError(i);
            const double fullErr = sqrt(statErr*statErr + sysErr*sysErr);

            if (fullErr < minErr)
            {
               valueOfMinErrHist = value;
               minStatErr = statErr;
               minSysErr = sysErr;
               minErr = fullErr;
            }

            graphSpectraVsPTForTsallisFit.
               AddPoint((pTBinRanges[i - 1] + pTBinRanges[i])/2., value);
            graphSpectraVsPTForTsallisFit.
               SetPointError(graphSpectraVsPTForTsallisFit.GetN() - 1, 0., fullErr);
         }

         if (valueOfMinErrHist < 0.) continue;

         distrResultingSpectraVsPTStatErr.SetBinContent(i, valueOfMinErrHist);
         distrResultingSpectraVsPTSysErr.SetBinContent(i, valueOfMinErrHist);

         distrResultingSpectraVsPTStatErr.SetBinError(i, minStatErr);
         distrResultingSpectraVsPTSysErr.SetBinError(i, minSysErr);
      }

      double xMin = pTBinRanges.front();
      double xMax = pTBinRanges.back();

      for (unsigned int i = 1; i <= pTNBins; i++)
      {
         if (distrResultingSpectraVsPTStatErr.GetBinContent(i) < 1e-31) xMin = pTBinRanges[i - 1];
         else break;
      }

      for (unsigned int i = pTNBins; i >= 1; i--)
      {
         if (distrResultingSpectraVsPTStatErr.GetBinContent(i) < 1e-31) xMax = pTBinRanges[i - 1];
         else break;
      }

      TF1 tsallisFit("tsallis", "0.5/pi*[0]*([1] - 1.)*([1] - 2.)/([2] + [3]*([1] - 1.))/"\
                                "([2] + [3])*([2] + sqrt(x^2 + [3]^2)/([2] + [3]))^(-[1])");
      tsallisFit.SetParameters(1., 2.5, 10.);
      tsallisFit.SetParLimits(1, 2., 30.);
      tsallisFit.FixParameter(3, resonanceMass);

      tsallisFit.SetRange(xMin + 0.1, xMax - 0.1);

      tsallisFit.SetLineStyle(2);
      tsallisFit.SetLineWidth(4);
      tsallisFit.SetLineColor(kRed - 3);

      for (unsigned int i = 0; i < fitNTries; i++)
      {
         graphSpectraVsPTForTsallisFit.Fit(&tsallisFit, "RQMBN EX0");

         // clearing previous points so that corrected ones can be written
         for (int j = graphSpectraVsPTForTsallisFit.GetN() - 1; j >= 0; j--)
         {
            graphSpectraVsPTForTsallisFit.RemovePoint(j);
         }

         for (unsigned int j = 0; j < pTNBins; j++)
         {
            if (distrResultingSpectraVsPTStatErr.GetBinContent(j + 1) < 1e-31) continue;

            const double norm = tsallisFit.Integral(pTBinRanges[j], pTBinRanges[j + 1])/
                                tsallisFit.Eval((pTBinRanges[j] + pTBinRanges[j + 1])/2.)/
                                (pTBinRanges[j + 1] - pTBinRanges[j]);

            distrResultingSpectraVsPTStatErr.
               SetBinContent(j + 1, distrResultingSpectraVsPTStatErr.GetBinContent(j + 1)/norm);
            distrResultingSpectraVsPTSysErr.
               SetBinContent(j + 1, distrResultingSpectraVsPTSysErr.GetBinContent(j + 1)/norm);

            distrResultingSpectraVsPTStatErr.
               SetBinError(j + 1, distrResultingSpectraVsPTStatErr.GetBinError(j + 1)/norm);
            distrResultingSpectraVsPTSysErr.
               SetBinError(j + 1, distrResultingSpectraVsPTSysErr.GetBinError(j + 1)/norm);

            for (unsigned int k = 0; k < spectrasVsPTStatErr.size(); k++)
            {
               if (spectrasVsPTStatErr[k]->GetBinContent(j + 1) < 1e-31) continue;

               spectrasVsPTStatErr[k]->
                  SetBinContent(j + 1, spectrasVsPTStatErr[k]->GetBinContent(j + 1)/norm);
               spectrasVsPTSysErr[k]->
                  SetBinContent(j + 1, spectrasVsPTSysErr[k]->GetBinContent(j + 1)/norm);

               spectrasVsPTStatErr[k]->
                  SetBinError(j + 1, spectrasVsPTStatErr[k]->GetBinError(j + 1)/norm);
               spectrasVsPTSysErr[k]->
                  SetBinError(j + 1, spectrasVsPTSysErr[k]->GetBinError(j + 1)/norm);

               graphSpectraVsPTForTsallisFit.
                  AddPoint((pTBinRanges[j] + pTBinRanges[j + 1])/2., 
                           spectrasVsPTStatErr[k]->GetBinContent(j + 1));

               const double fullErr = sqrt(pow(spectrasVsPTSysErr[k]->GetBinError(j + 1), 2) +
                                           pow(spectrasVsPTStatErr[k]->GetBinError(j + 1), 2));
               graphSpectraVsPTForTsallisFit.
                  SetPointError(graphSpectraVsPTForTsallisFit.GetN() - 1, 0., fullErr);
            }
         }
      }

      TCanvas canvSpectra("resulting spectra canv", "", 800, 800);

      gPad->SetRightMargin(0.002); gPad->SetTopMargin(0.002); 
      gPad->SetLeftMargin(0.152); gPad->SetBottomMargin(0.112);

      gPad->SetLogy();

      ROOTTools::DrawFrame(xMin - 0.1, yMin/5., xMax + 0.1, yMax*5., "", "p_{T} [GeV/c]", 
                           "1/(2#pip_{T}) d^{2} N/dp_{T}/dy [(GeV/c)^{-2}]");

      canvSpectra.SetFillStyle(4000);
      canvSpectra.SetFrameFillColor(0);
      canvSpectra.SetFrameFillStyle(0);
      canvSpectra.SetFrameBorderMode(0);

      tsallisFit.Draw("SAME");
      distrResultingSpectraVsPTStatErr.Draw("SAME");
      distrResultingSpectraVsPTSysErr.Draw("SAME E2");

      distrResultingSpectraVsPTStatErr.Clone()->Write();
      distrResultingSpectraVsPTSysErr.Clone()->Write();

      ROOTTools::PrintCanvas(&canvSpectra, outputDir + "/" + resonanceName + 
                             "_spectra_" + centralityName);

      std::vector<TH1D *> spectraRatiosVsPTStatErr;
      std::vector<TH1D *> spectraRatiosVsPTSysErr;

      double ratioMin = 1e31;
      double ratioMax = 1e-31;

      for (unsigned int i = 0; i < spectrasVsPTStatErr.size(); i++)
      {
         spectraRatiosVsPTStatErr.
            emplace_back(static_cast<TH1D *>(spectrasVsPTStatErr[i]->Clone()));
         spectraRatiosVsPTSysErr.
            emplace_back(static_cast<TH1D *>(spectrasVsPTSysErr[i]->Clone()));

         spectraRatiosVsPTStatErr[i]->SetLineWidth(4);
         spectraRatiosVsPTSysErr[i]->SetLineWidth(4);

         spectraRatiosVsPTStatErr[i]->Divide(&tsallisFit);
         spectraRatiosVsPTSysErr[i]->Divide(&tsallisFit);

         for (unsigned int j = 1; j < pTNBins; j++)
         {
            if (spectraRatiosVsPTStatErr[i]->GetBinContent(j) < 1e-15) continue;

            ratioMin = CppTools::Minimum(ratioMin, spectraRatiosVsPTStatErr[i]->GetBinContent(j));
            ratioMax = CppTools::Maximum(ratioMax, spectraRatiosVsPTStatErr[i]->GetBinContent(j));
         }
      }

      TCanvas canvAllSpectra("all spectra canv", "", 800, 1000);

      canvAllSpectra.SetFillStyle(4000);
      canvAllSpectra.SetFrameFillColor(0);
      canvAllSpectra.SetFrameFillStyle(0);
      canvAllSpectra.SetFrameBorderMode(0);

      canvAllSpectra.Divide(1, 2, 0., 0.);

      canvAllSpectra.cd(1);

      gPad->SetLogy();

      gPad->SetPad(0., 0.3, 1., 1.);
      gPad->SetRightMargin(0.002); gPad->SetTopMargin(0.002); 
      gPad->SetLeftMargin(0.152); gPad->SetBottomMargin(0.);

      ROOTTools::DrawFrame(xMin - 0.1, yMin/5., xMax + 0.1, yMax*5., "", "", 
                           "1/(2#pip_{T}) d^{2} N/dp_{T}/dy [(GeV/c)^{-2}]");

      tsallisFit.Draw("SAME");

      for (unsigned int i = 0; i < spectrasVsPTStatErr.size(); i++)
      {
         spectrasVsPTStatErr[i]->SetLineColor(methodColors[i]);
         spectrasVsPTStatErr[i]->Draw("SAME");

         spectrasVsPTSysErr[i]->SetFillColorAlpha(methodColors[i], 0.3);
         spectrasVsPTSysErr[i]->SetFillStyle(1001);
         spectrasVsPTSysErr[i]->Draw("SAME E2");

         const std::string methodName = 
            inputYAMLResonance["pair_selection_methods"][i]["name"].as<std::string>();

         legend.AddEntry(spectrasVsPTStatErr[i], methodName.c_str(), "L");
      }

      legend.Draw();

      canvAllSpectra.cd(2);

      gPad->SetPad(0., 0., 1., 0.3);
      gPad->SetRightMargin(0.002); gPad->SetTopMargin(0.); 
      gPad->SetLeftMargin(0.152); gPad->SetBottomMargin(0.25);

      ROOTTools::DrawFrame(xMin - 0.1, ratioMin/1.1, xMax + 0.1, ratioMax*1.1, 
                           "", "p_{T} [GeV/c]", "Data/Fit", 1., 0.7, 0.11, 0.11);

      if (ratioMin/1.1 < 1. && ratioMax*1.1 > 1.)
      {
         TLine line(xMin - 0.1, 1., xMax + 0.1, 1.);
         line.SetLineColorAlpha(kBlack, 0.5);
         line.SetLineStyle(2);
         line.SetLineWidth(4);
         line.Clone()->Draw();
      }

      for (unsigned int i = 0; i < spectraRatiosVsPTStatErr.size(); i++)
      {
         spectraRatiosVsPTStatErr[i]->SetLineColor(methodColors[i]);
         spectraRatiosVsPTStatErr[i]->Draw("SAME");

         spectraRatiosVsPTSysErr[i]->SetFillColorAlpha(methodColors[i], 0.3);
         spectraRatiosVsPTSysErr[i]->SetFillStyle(1001);
         spectraRatiosVsPTSysErr[i]->Draw("SAME E2");
      }

      ROOTTools::PrintCanvas(&canvAllSpectra, outputDir + "/" + resonanceName + 
                             "_spectra_" + centralityName + "_all");

      legend.Clear();

      if (estimateFactors)
      {
         std::vector<TH1D *> distrRABsVsPTStatErr;
         std::vector<TH1D *> distrRABsVsPTSysErr;

         // scale for R_{AB} = c_{bias}/N_{coll}
         const double scaleRAB = centralityBin["bias_factor"].as<double>()/
                                 centralityBin["N_coll"].as<double>()*42.2;

         double maxRAB = 1e-31;

         for (unsigned int i = 0; i < spectrasVsPTStatErr.size(); i++)
         {
            distrRABsVsPTStatErr.
               emplace_back(static_cast<TH1D *>(spectrasVsPTStatErr[i]->Clone()));
            distrRABsVsPTSysErr.
               emplace_back(static_cast<TH1D *>(spectrasVsPTSysErr[i]->Clone()));

            distrRABsVsPTStatErr.back()->Divide(distrSpectraPPVsPTStatErr);
            distrRABsVsPTSysErr.back()->Divide(distrSpectraPPVsPTSysErr);

            distrRABsVsPTStatErr.back()->Scale(scaleRAB);
            distrRABsVsPTSysErr.back()->Scale(scaleRAB);

            for (unsigned int j = 1; j < pTNBins; j++)
            {
               if (distrRABsVsPTStatErr.back()->GetBinContent(j) < 1e-15) continue;

               maxRAB = CppTools::Maximum(maxRAB, distrRABsVsPTStatErr.back()->GetBinContent(j));
            }
         }

         maxRAB = CppTools::Maximum(maxRAB, 1.99);

         TCanvas canvAllRAB("all spectra canv", "", 800, 800);

         canvAllRAB.SetFillStyle(4000);
         canvAllRAB.SetFrameFillColor(0);
         canvAllRAB.SetFrameFillStyle(0);
         canvAllRAB.SetFrameBorderMode(0);

         gPad->SetRightMargin(0.002); gPad->SetTopMargin(0.002); 
         gPad->SetLeftMargin(0.1); gPad->SetBottomMargin(0.112);

         ROOTTools::DrawFrame(xMin - 0.1, 0., xMax + 0.1, maxRAB, 
                              "", "p_{T} [GeV/c]", "R_{AB}", 1., 0.95);

         if (maxRAB > 1.)
         {
            TLine line(xMin - 0.1, 1., xMax + 0.1, 1.);
            line.SetLineColorAlpha(kBlack, 0.5);
            line.SetLineStyle(2);
            line.SetLineWidth(4);
            line.Clone()->Draw();
         }

         for (unsigned int i = 0; i < distrRABsVsPTStatErr.size(); i++)
         {
            distrRABsVsPTStatErr[i]->SetLineColor(methodColors[i]);
            distrRABsVsPTStatErr[i]->Draw("SAME");

            distrRABsVsPTSysErr[i]->SetFillColorAlpha(methodColors[i], 0.3);
            distrRABsVsPTSysErr[i]->SetFillStyle(1001);
            distrRABsVsPTSysErr[i]->Draw("SAME E2");

            const std::string methodName = 
               inputYAMLResonance["pair_selection_methods"][i]["name"].as<std::string>();

            legend.AddEntry(distrRABsVsPTStatErr[i], methodName.c_str(), "PLC");
         }

         legend.Draw();

         TLatex tlText;

         tlText.SetTextFont(52);
         tlText.SetTextSize(0.05);

         tlText.DrawLatexNDC(0.15, 0.15, centralityNameTex.c_str());

         ROOTTools::PrintCanvas(&canvAllRAB, outputDir + "/" + resonanceName + 
                                "_RAB_all_methods_" + centralityName);

         TH1D *distrResultingRABVsPTStatErr = 
            static_cast<TH1D *>(distrResultingSpectraVsPTStatErr.Clone());
         TH1D *distrResultingRABVsPTSysErr = 
            static_cast<TH1D *>(distrResultingSpectraVsPTSysErr.Clone());

         distrResultingRABVsPTStatErr->Divide(distrSpectraPPVsPTStatErr);
         distrResultingRABVsPTSysErr->Divide(distrSpectraPPVsPTSysErr);

         distrResultingRABVsPTStatErr->Scale(scaleRAB);
         distrResultingRABVsPTSysErr->Scale(scaleRAB);

         for (int i = 1; i <= distrResultingRABVsPTStatErr->GetXaxis()->GetNbins(); i++)
         {
            if (distrResultingRABVsPTStatErr->GetBinContent(i) < 1e-7)
            {
               distrResultingRABVsPTStatErr->SetBinContent(i, -1e-7);
               distrResultingRABVsPTSysErr->SetBinContent(i, -1e-7);
            }
         }

         distrResultingRABVsPTStatErr->Write("RAB vs pT with stat errors");
         distrResultingRABVsPTSysErr->Write("RAB vs pT with sys errors");
      }
   }

   resultsOutputFile->Close();

   CppTools::PrintInfo("Results (spectra and RAB) were succesfully evaluated");
   CppTools::PrintInfo("Results were written in " + resultsOutputFileName);
   CppTools::PrintInfo("Pictures were written in " + outputDir + " directory");

   return 0;
}

#endif /* ESTIMATE_RESULTS_CPP */
