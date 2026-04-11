/** 
 *  @file   EstimateSpectra.hpp 
 *  @brief  Contains realisations of functions and variables that are used for estimation of invariant pT spectra
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysis).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef ESTIMATE_SPECTRA_CPP
#define ESTIMATE_SPECTRA_CPP

#include "EstimateSpectra.hpp"

using namespace EstimateSpectra;

int main(int argc, char **argv)
{
   if (argc < 3) 
   {
      CppTools::PrintError("Expected 2 parameters while " + std::to_string(argc - 1) + " "\
                           "parameter(s) were provided \n Usage: bin/EstimateSpectra "\
                           "inputYAMLName taxiNumber");
   }


   gStyle->SetPalette(kGreenRedViolet);

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

   inputFileName = "data/RawYields/" + runName + "/Resonance/" + 
                   std::to_string(taxiNumber) + "_" + resonanceName + ".root";
   CppTools::CheckInputFile(inputFileName);
   inputFile = TFile::Open(inputFileName.c_str());

   inputRecEffFileName = "data/Parameters/ResonanceEff/" + runName + "/" + resonanceName + ".root";
   CppTools::CheckInputFile(inputRecEffFileName);
   inputRecEffFile = TFile::Open(inputRecEffFileName.c_str());

   pTNBins = inputYAMLResonance["pt_bins"].size();

   for (unsigned int i = 0; i < pTNBins; i++)
   {
      pTBinRanges.push_back(inputYAMLResonance["pt_bins"][i]["min"].as<double>());
   }
   pTBinRanges.push_back(inputYAMLResonance["pt_bins"][pTNBins - 1]["max"].as<double>());

   const std::string spectraOutputDir = "data/Spectra/" + runName;
   std::filesystem::create_directories(spectraOutputDir);

   const std::string outputDir = "output/Spectra/" + runName + "/" + std::to_string(taxiNumber);
   std::filesystem::create_directories(outputDir);

   double spectraNorm = inputYAMLResonance["branching_ratio"].as<double>();
   /*
   double relativeSpectraNormErr = inputYAMLResonance["branching_ratio_uncertainty"].as<double>()/
                                   spectraNorm;
                                   */

   if (inputYAMLResonance["has_antiparticle"]) spectraNorm *= 2.;

   TFile *spectraOutputFile = TFile::Open((spectraOutputDir + "/" + std::to_string(taxiNumber) + 
                                           "_" + resonanceName + ".root").c_str(),
                                          "RECREATE");

   for (const YAML::Node& centralityBin : inputYAMLResonance["centrality_bins"])
   {
      const std::string centralityName = centralityBin["name"].as<std::string>();

      std::vector<TH1D *> spectrasVsPTWithStatErrors;
      //std::vector<TH1D *> spectrasWithSysErrors;

      double yMin = 1e31;
      double yMax = 1e-31;

      TLegend legend(0.6, 0.7, 0.95, 0.95);
      legend.SetLineColorAlpha(0, 0.);
      legend.SetFillColorAlpha(0, 0.);

      for (const YAML::Node& method : inputYAMLResonance["pair_selection_methods"])
      {
         const std::string methodName = method["name"].as<std::string>();

         TH1D *rawYieldVsPTWithStatErrors = 
            static_cast<TH1D *>(inputFile->Get((methodName + "/" + centralityName + 
                                                "/raw yield vs pT").c_str()));

         if (!rawYieldVsPTWithStatErrors) 
         {
            CppTools::Print(methodName + "/" + centralityName + "/raw yield vs pT");
            CppTools::PrintError("No raw yield distribution was found in file " + inputFileName + 
                                 " for " + methodName + " in centrality " + centralityName);
         }

         TH1D *recEffVsPT = 
            static_cast<TH1D *>(inputRecEffFile->Get((methodName + 
                                                     "/reconstruction efficiency vs pT").c_str()));

         if (!recEffVsPT)
         {
            CppTools::PrintError("No reconstruction efficiency was found in file " + 
                                 inputRecEffFileName + " for " + methodName);
         }

         rawYieldVsPTWithStatErrors->Divide(recEffVsPT);
         rawYieldVsPTWithStatErrors->Scale(1./spectraNorm);

         spectrasVsPTWithStatErrors.emplace_back(rawYieldVsPTWithStatErrors);

         for (int i = 1; i <= rawYieldVsPTWithStatErrors->GetXaxis()->GetNbins(); i++)
         {
            if (rawYieldVsPTWithStatErrors->GetBinContent(i) < 1e-31) continue;

            yMin = CppTools::Minimum(yMin, rawYieldVsPTWithStatErrors->GetBinContent(i));
            yMax = CppTools::Maximum(yMax, rawYieldVsPTWithStatErrors->GetBinContent(i));
         }
      }

      spectraOutputFile->mkdir(centralityName.c_str());
      spectraOutputFile->cd(centralityName.c_str());

      TH1D resultingSpectraVsPT("spectra vs pT", "", pTNBins, &pTBinRanges[0]);

      resultingSpectraVsPT.SetLineColor(kBlack);
      resultingSpectraVsPT.SetLineWidth(4);

      for (int i = 1; i <= resultingSpectraVsPT.GetXaxis()->GetNbins(); i++)
      {
         double minErr = 1e31;
         double valueOfMinErrHist = -1.;
         for (unsigned int j = 0; j < spectrasVsPTWithStatErrors.size(); j++)
         {
            if (spectrasVsPTWithStatErrors[j]->GetBinContent(i) < 1e-31) continue;

            if (spectrasVsPTWithStatErrors[j]->GetBinError(i) < minErr)
            {
               valueOfMinErrHist = spectrasVsPTWithStatErrors[j]->GetBinContent(i);
               minErr = spectrasVsPTWithStatErrors[j]->GetBinError(i);
            }
         }

         if (valueOfMinErrHist < 0.) continue;

         resultingSpectraVsPT.SetBinContent(i, valueOfMinErrHist);
         resultingSpectraVsPT.SetBinError(i, minErr);
      }

      double xMin = pTBinRanges.front();
      double xMax = pTBinRanges.back();

      for (int i = 1; i <= resultingSpectraVsPT.GetXaxis()->GetNbins(); i++)
      {
         if (resultingSpectraVsPT.GetBinContent(i) < 1e-31) xMin = pTBinRanges[i - 1];
         else break;
      }

      for (int i = resultingSpectraVsPT.GetXaxis()->GetNbins(); i >= 1; i--)
      {
         if (resultingSpectraVsPT.GetBinContent(i) < 1e-31) xMax = pTBinRanges[i - 1];
         else break;
      }

      TF1 tsallisFit("tsallis", "0.5/pi*[0]*([1] - 1.)*([1] - 2.)/([2] + [3]*([1] - 1.))/"\
                                "([2] + [3])*([2] + sqrt(x^2 + [3]^2)/([2] + [3]))^(-[1])");
      tsallisFit.SetParameters(1., 1.1, 10.);
      tsallisFit.SetParLimits(1, 1.1, 30.);
      tsallisFit.FixParameter(3, resonanceMass);

      tsallisFit.SetRange(xMin + 0.1, xMax - 0.1);

      tsallisFit.SetLineStyle(2);
      tsallisFit.SetLineWidth(4);
      tsallisFit.SetLineColor(kRed - 3);

      for (unsigned int i = 0; i < fitNTries; i++)
      {
         resultingSpectraVsPT.Fit(&tsallisFit, "QBN");

         for (unsigned int j = 0; j < pTNBins; j++)
         {
            if (resultingSpectraVsPT.GetBinContent(j + 1) < 1e-31) continue;

            const double norm = tsallisFit.Integral(pTBinRanges[j], pTBinRanges[j + 1])/
                                tsallisFit.Eval((pTBinRanges[j] + pTBinRanges[j + 1])/2.)/
                                (pTBinRanges[j + 1] - pTBinRanges[j]);

            resultingSpectraVsPT.
               SetBinContent(j + 1, resultingSpectraVsPT.GetBinContent(j + 1)/norm);
            resultingSpectraVsPT.
               SetBinError(j + 1, resultingSpectraVsPT.GetBinError(j + 1)/norm);

            for (int k = 0; k < static_cast<int>(spectrasVsPTWithStatErrors.size()); k++)
            {
               if (spectrasVsPTWithStatErrors[k]->GetBinContent(j + 1) < 1e-31) continue;

               spectrasVsPTWithStatErrors[k]->
                  SetBinContent(j + 1, spectrasVsPTWithStatErrors[k]->GetBinContent(j + 1)/norm);
               spectrasVsPTWithStatErrors[k]->
                  SetBinError(j + 1, spectrasVsPTWithStatErrors[k]->GetBinError(j + 1)/norm);
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
      resultingSpectraVsPT.Draw("SAME");

      resultingSpectraVsPT.Clone()->Write();

      ROOTTools::PrintCanvas(&canvSpectra, outputDir + "/" + centralityName + "");

      std::vector<TH1D *> spectraRatiosVsPT;

      double ratioMin = 1e31;
      double ratioMax = 1e-31;

      for (unsigned int i = 0; i < spectrasVsPTWithStatErrors.size(); i++)
      {
         spectraRatiosVsPT.
            emplace_back(static_cast<TH1D *>(spectrasVsPTWithStatErrors[i]->Clone()));

         spectraRatiosVsPT[i]->SetLineWidth(4);

         spectraRatiosVsPT[i]->Divide(&tsallisFit);

         for (int j = 1; j < spectraRatiosVsPT[i]->GetXaxis()->GetNbins(); j++)
         {
            if (spectraRatiosVsPT[i]->GetBinContent(j) < 1e-15) continue;

            ratioMin = CppTools::Minimum(ratioMin, spectraRatiosVsPT[i]->GetBinContent(j));
            ratioMax = CppTools::Maximum(ratioMax, spectraRatiosVsPT[i]->GetBinContent(j));
         }
      }

      TCanvas canvAll("all spectra canv", "", 800, 1000);

      canvAll.SetFillStyle(4000);
      canvAll.SetFrameFillColor(0);
      canvAll.SetFrameFillStyle(0);
      canvAll.SetFrameBorderMode(0);

      canvAll.Divide(1, 2, 0., 0.);

      canvAll.cd(1);

      gPad->SetLogy();

      gPad->SetPad(0., 0.3, 1., 1.);
      gPad->SetRightMargin(0.002); gPad->SetTopMargin(0.002); 
      gPad->SetLeftMargin(0.152); gPad->SetBottomMargin(0.);

      ROOTTools::DrawFrame(xMin - 0.1, yMin/5., xMax + 0.1, yMax*5., "", "", 
                           "1/(2#pip_{T}) d^{2} N/dp_{T}/dy [(GeV/c)^{-2}]");

      tsallisFit.Draw("SAME");

      for (int i = 0; i < static_cast<int>(spectrasVsPTWithStatErrors.size()); i++)
      {
         spectrasVsPTWithStatErrors[i]->Draw("SAME PLC");

         const std::string methodName = 
            inputYAMLResonance["pair_selection_methods"][i]["name"].as<std::string>();

         legend.AddEntry(spectrasVsPTWithStatErrors[i], methodName.c_str(), "PLC");

         //spectrasVsPTWithStatErrors[i]->Write(("spectra vs pT, " + methodName).c_str());
      }

      legend.Draw();

      canvAll.cd(2);

      gPad->SetPad(0., 0., 1., 0.3);
      gPad->SetRightMargin(0.002); gPad->SetTopMargin(0.); 
      gPad->SetLeftMargin(0.152); gPad->SetBottomMargin(0.25);

      ROOTTools::DrawFrame(xMin - 0.1, ratioMin/1.1, xMax + 0.1, ratioMax*1.1, 
                           "", "p_{T} [GeV/c]", "Data/Fit", 1., 0.7, 0.11, 0.11);

      if (ratioMin < 1. && ratioMax > 1.)
      {
         TLine line(xMin - 0.1, 1., xMax + 0.1, 1.);
         line.SetLineColorAlpha(kBlack, 0.5);
         line.SetLineStyle(2);
         line.SetLineWidth(4);
         line.Clone()->Draw();
      }

      for (int i = 0; i < static_cast<int>(spectraRatiosVsPT.size()); i++)
      {
         spectraRatiosVsPT[i]->Draw("SAME PLC");
      }

      ROOTTools::PrintCanvas(&canvAll, outputDir + "/" + centralityName + "_all");
   }

   spectraOutputFile->Close();

   CppTools::PrintInfo("Spectra were succesfully evaluated");

   return 0;
}

#endif /* ESTIMATE_SPECTRA_CPP */
