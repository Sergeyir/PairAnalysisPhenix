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

   inputRecEffFileName = "data/Parameters/RecEffResonance/" + runName + "/" + resonanceName + ".root";
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

   TH1D *distrSpectraPPVsPTWithStatErr = nullptr;
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
         distrSpectraPPVsPTWithStatErr = 
            static_cast<TH1D *>(spectraPPFile->Get("spectra vs pT with stat err"));
         //spectraPPVsPTWithSysErr = 
         //   static_cast<TH1D *>(spectraPPFile->Get("spectra vs pT with sys err"));
         if (!distrSpectraPPVsPTWithStatErr)
         {
            CppTools::PrintError("No p+p spectra was found in file " + spectraPPFileName);
         }
      }
   }

   for (const YAML::Node& centralityBin : inputYAMLResonance["centrality_bins"])
   {
      const std::string centralityName = centralityBin["name"].as<std::string>();
      const std::string centralityNameTex = centralityBin["name_tex"].as<std::string>();

      std::vector<TH1D *> spectrasVsPTWithStatErr;
      //std::vector<TH1D *> spectrasVsPTWithSysErrors;

      double yMin = 1e31;
      double yMax = 1e-31;

      TLegend legend(0.6, 0.7, 0.95, 0.95);
      legend.SetLineColorAlpha(0, 0.);
      legend.SetFillColorAlpha(0, 0.);

      for (const YAML::Node& method : inputYAMLResonance["pair_selection_methods"])
      {
         const std::string methodName = method["name"].as<std::string>();

         inputFileName = "data/RawYields/" + runName + "/Resonance/" + 
                         std::to_string(taxiNumber) + "_" + resonanceName + 
                         "_" + methodName + ".root";
         CppTools::CheckInputFile(inputFileName);
         inputFile = TFile::Open(inputFileName.c_str());

         TH1D *rawYieldVsPTWithStatErr = 
            static_cast<TH1D *>(inputFile->Get((centralityName + "/raw yield vs pT").c_str()));

         if (!rawYieldVsPTWithStatErr) 
         {
            CppTools::Print(methodName + "/" + centralityName + "/raw yield vs pT");
            CppTools::PrintError("No raw yield distribution was found in file " + inputFileName + 
                                 " for " + methodName + " in centrality " + centralityName);
         }

         TH1D *recEffVsPT = static_cast<TH1D *>
            (inputRecEffFile->Get((methodName + "/reconstruction efficiency "\
                                   "vs pT with stat errors").c_str()));

         if (!recEffVsPT)
         {
            CppTools::PrintError("No reconstruction efficiency was found in file " + 
                                 inputRecEffFileName + " for " + methodName);
         }

         rawYieldVsPTWithStatErr->Divide(recEffVsPT);
         rawYieldVsPTWithStatErr->Scale(1./spectraNorm);

         spectrasVsPTWithStatErr.emplace_back(rawYieldVsPTWithStatErr);

         for (unsigned int i = 1; i <= pTNBins; i++)
         {
            if (rawYieldVsPTWithStatErr->GetBinContent(i) < 1e-31) continue;

            yMin = CppTools::Minimum(yMin, rawYieldVsPTWithStatErr->GetBinContent(i));
            yMax = CppTools::Maximum(yMax, rawYieldVsPTWithStatErr->GetBinContent(i));
         }
      }

      resultsOutputFile->mkdir(centralityName.c_str());
      resultsOutputFile->cd(centralityName.c_str());

      TH1D distrResultingSpectraVsPTStatErr("spectra vs pT with stat errors", "", 
                                            pTNBins, &pTBinRanges[0]);

      distrResultingSpectraVsPTStatErr.SetLineColor(kBlack);
      distrResultingSpectraVsPTStatErr.SetLineWidth(4);

      // graph containing all points across different methods for better Tsallis fit
      TGraphErrors graphSpectraVsPTForTsallisFit;

      for (unsigned int i = 1; i <= pTNBins; i++)
      {
         double minErr = 1e31;
         double minStatErr = 1e31;
         //double minSysErr = 1e31;

         double valueOfMinErrHist = -1.;
         for (unsigned int j = 0; j < spectrasVsPTWithStatErr.size(); j++)
         {
            if (spectrasVsPTWithStatErr[j]->GetBinContent(i) < 1e-31) continue;

            const double value = spectrasVsPTWithStatErr[j]->GetBinContent(i);
            const double statErr = spectrasVsPTWithStatErr[j]->GetBinError(i);
            const double fullErr = statErr;

            if (fullErr < minErr)
            {
               valueOfMinErrHist = value;
               minStatErr = statErr;
               minErr = fullErr;
            }

            graphSpectraVsPTForTsallisFit.
               AddPoint((pTBinRanges[i - 1] + pTBinRanges[i])/2., value);
            graphSpectraVsPTForTsallisFit.
               SetPointError(graphSpectraVsPTForTsallisFit.GetN() - 1, 0., fullErr);
         }

         if (valueOfMinErrHist < 0.) continue;

         distrResultingSpectraVsPTStatErr.SetBinContent(i, valueOfMinErrHist);
         distrResultingSpectraVsPTStatErr.SetBinError(i, minStatErr);
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

         // clearing previos points so that corrected ones can be written
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
            distrResultingSpectraVsPTStatErr.
               SetBinError(j + 1, distrResultingSpectraVsPTStatErr.GetBinError(j + 1)/norm);

            for (unsigned int k = 0; k < spectrasVsPTWithStatErr.size(); k++)
            {
               if (spectrasVsPTWithStatErr[k]->GetBinContent(j + 1) < 1e-31) continue;

               spectrasVsPTWithStatErr[k]->
                  SetBinContent(j + 1, spectrasVsPTWithStatErr[k]->GetBinContent(j + 1)/norm);
               spectrasVsPTWithStatErr[k]->
                  SetBinError(j + 1, spectrasVsPTWithStatErr[k]->GetBinError(j + 1)/norm);

               graphSpectraVsPTForTsallisFit.
                  AddPoint((pTBinRanges[j] + pTBinRanges[j + 1])/2., 
                           spectrasVsPTWithStatErr[k]->GetBinContent(j + 1));
               graphSpectraVsPTForTsallisFit.
                  SetPointError(graphSpectraVsPTForTsallisFit.GetN() - 1, 0.,
                                spectrasVsPTWithStatErr[k]->GetBinError(j + 1));
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

      distrResultingSpectraVsPTStatErr.Clone()->Write();

      ROOTTools::PrintCanvas(&canvSpectra, outputDir + "/spectra_" + centralityName);

      std::vector<TH1D *> spectraRatiosVsPT;

      double ratioMin = 1e31;
      double ratioMax = 1e-31;

      for (unsigned int i = 0; i < spectrasVsPTWithStatErr.size(); i++)
      {
         spectraRatiosVsPT.
            emplace_back(static_cast<TH1D *>(spectrasVsPTWithStatErr[i]->Clone()));

         spectraRatiosVsPT[i]->SetLineWidth(4);

         spectraRatiosVsPT[i]->Divide(&tsallisFit);

         for (unsigned int j = 1; j < pTNBins; j++)
         {
            if (spectraRatiosVsPT[i]->GetBinContent(j) < 1e-15) continue;

            ratioMin = CppTools::Minimum(ratioMin, spectraRatiosVsPT[i]->GetBinContent(j));
            ratioMax = CppTools::Maximum(ratioMax, spectraRatiosVsPT[i]->GetBinContent(j));
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

      for (unsigned int i = 0; i < spectrasVsPTWithStatErr.size(); i++)
      {
         spectrasVsPTWithStatErr[i]->Draw("SAME PLC");

         const std::string methodName = 
            inputYAMLResonance["pair_selection_methods"][i]["name"].as<std::string>();

         legend.AddEntry(spectrasVsPTWithStatErr[i], methodName.c_str(), "PLC");

         //spectrasVsPTWithStatErr[i]->Write(("spectra vs pT with stat errors, " + 
         //                                      methodName).c_str());
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

      for (unsigned int i = 0; i < spectraRatiosVsPT.size(); i++)
      {
         spectraRatiosVsPT[i]->Draw("SAME PLC");
      }

      ROOTTools::PrintCanvas(&canvAllSpectra, outputDir + "/spectra_" + centralityName + "_all");

      legend.Clear();

      if (estimateFactors)
      {
         std::vector<TH1D *>distrRABsVsPTWithStatErr;
         //std::vector<TH1D *>distrRABsVsPTWithSysErr;

         // scale for R_{AB} = c_{bias}/N_{coll}
         const double scaleRAB = centralityBin["bias_factor"].as<double>()/
                                 centralityBin["N_coll"].as<double>()*42.2;

         double maxRAB = 1e-31;

         for (unsigned int i = 0; i < spectrasVsPTWithStatErr.size(); i++)
         {
            distrRABsVsPTWithStatErr.
               emplace_back(static_cast<TH1D *>(spectrasVsPTWithStatErr[i]->Clone()));

            distrRABsVsPTWithStatErr.back()->Divide(distrSpectraPPVsPTWithStatErr);
            distrRABsVsPTWithStatErr.back()->Scale(scaleRAB);

            for (unsigned int j = 1; j < pTNBins; j++)
            {
               if (distrRABsVsPTWithStatErr.back()->GetBinContent(j) < 1e-15) continue;

               maxRAB = CppTools::Maximum(maxRAB, distrRABsVsPTWithStatErr.back()->GetBinContent(j));
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

         for (unsigned int i = 0; i < distrRABsVsPTWithStatErr.size(); i++)
         {
            distrRABsVsPTWithStatErr[i]->Draw("SAME PLC");

            const std::string methodName = 
               inputYAMLResonance["pair_selection_methods"][i]["name"].as<std::string>();

            legend.AddEntry(distrRABsVsPTWithStatErr[i], methodName.c_str(), "PLC");
         }

         legend.Draw();

         TLatex tlText;

         tlText.SetTextFont(52);
         tlText.SetTextSize(0.05);

         tlText.DrawLatexNDC(0.15, 0.15, centralityNameTex.c_str());

         ROOTTools::PrintCanvas(&canvAllRAB, outputDir + "/RAB_all_methods_" + centralityName);

         TH1D *distrResultingRABVsPTStatErr = 
            static_cast<TH1D *>(distrResultingSpectraVsPTStatErr.Clone());

         distrResultingRABVsPTStatErr->Divide(distrSpectraPPVsPTWithStatErr);
         distrResultingRABVsPTStatErr->Scale(scaleRAB);

         distrResultingRABVsPTStatErr->Write("RAB vs pT with stat errors");
      }
   }

   resultsOutputFile->Close();

   CppTools::PrintInfo("Results (spectra and RAB) were succesfully evaluated");
   CppTools::PrintInfo("Results were written in " + resultsOutputFileName);
   CppTools::PrintInfo("Pictures were written in " + outputDir + " directory");

   return 0;
}

#endif /* ESTIMATE_RESULTS_CPP */
