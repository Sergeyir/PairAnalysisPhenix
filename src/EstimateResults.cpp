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
   if (argc < 2 || argc > 5) 
   {
      CppTools::PrintError("Expected 2 parameters while " + std::to_string(argc - 1) + " "\
                           "parameter(s) were provided \n Usage: bin/EstimateResults "\
                           "inputYAMLName taxiNumber=default* taxiNumberWithLoosenedCuts=default*"\
                           "taxiNumberWithTightenedCuts=default*"\
                           "* default taxi job numbers are defined in inputYAMLName; "\
                           "if default value in inputYAMLName is -9999 then you have to provide the "\
                           "taxi number manually or change the default value\n"\
                           "  default value of -9999 for cut variations will make the program skip "\
                           "cut variation analysis and systematics");
   }

   gStyle->SetPalette(kSouthWest);

   CppTools::CheckInputFile(argv[1]);

   inputYAMLResonance.OpenFile(argv[1]);
   inputYAMLResonance.CheckStatus("resonance");

   if (argc > 2) taxiNumber = std::stoi(argv[2]);
   else taxiNumber = inputYAMLResonance["taxi_job"].as<int>();
   if (argc > 4) 
   {
      taxiNumberLoosenedCuts = std::stoi(argv[3]);
      taxiNumberTightenedCuts = std::stoi(argv[4]);
   }
   else 
   {
      taxiNumberLoosenedCuts = inputYAMLResonance["taxi_job_loosened_cuts"].as<int>();
      taxiNumberTightenedCuts = inputYAMLResonance["taxi_job_tightened_cuts"].as<int>();
   }

   if (taxiNumber == -9999) CppTools::PrintError("Taxi job number was not defined");

   if (taxiNumberLoosenedCuts == -9999) 
   {
      CppTools::PrintWarning("Taxi job number for loosened cut variation was not defined");
      doCutsVarSys = false;
   }
   if (taxiNumberTightenedCuts == -9999) 
   {
      CppTools::PrintWarning("Taxi job number for tightened cut variation was not defined");
      doCutsVarSys = false;
   }

   if (doCutsVarSys == false)
   {
      CppTools::PrintWarning("Cuts variation systematics evaluation will be disabled");
   }

   runName = inputYAMLResonance["run_name"].as<std::string>();

   inputYAMLMain.OpenFile("input/" + runName + "/main.yaml");
   inputYAMLMain.CheckStatus("main");

   gStyle->SetOptStat(0);
   gErrorIgnoreLevel = kWarning;

   gROOT->SetBatch(kTRUE);
   gStyle->SetOptStat(kFALSE);

   TH1::AddDirectory(kFALSE);
   ROOT::EnableImplicitMT(std::thread::hardware_concurrency());

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

   const std::string outputDirResults = "output/Results/" + runName + "/" + std::to_string(taxiNumber);
   std::filesystem::create_directories(outputDirResults);

   const std::string outputDirSys = "output/Systematics/" + runName + "/" + std::to_string(taxiNumber);
   std::filesystem::create_directories(outputDirSys);

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
      // searching for pT min bin for RAB
      for (unsigned int i = 0; i < pTBinRanges.size(); i++)
      {
         bool pTBinMinChosen = false;
         for (int j = 1; j <= distrSpectraPPVsPTStatErr->GetXaxis()->GetNbins(); j++)
         {
            const double pT = distrSpectraPPVsPTStatErr->GetXaxis()->GetBinCenter(j);

            if (pT > pTBinRanges[i] && pT < pTBinRanges[i + 1])
            {
               pTBinMinRAB = i;
               pTBinMinChosen = true;
               break;
            }
         }
         if (pTBinMinChosen) break;
      }
      // searching for pT max bin for RAB
      for (int i = static_cast<int>(pTBinRanges.size()) - 1; i >= 0; i--)
      {
         bool pTBinMaxChosen = false;
         for (int j = distrSpectraPPVsPTStatErr->GetXaxis()->GetNbins(); j > 0; j--)
         {
            const double pT = distrSpectraPPVsPTStatErr->GetXaxis()->GetBinCenter(j);

            if (pT > pTBinRanges[i] && pT < pTBinRanges[i + 1])
            {
               pTBinMaxRAB = i;
               pTBinMaxChosen = true;
               break;
            }
         }
         if (pTBinMaxChosen) break;
      }
      if (pTBinMinRAB == -1 || pTBinMaxRAB == -1)
      {
         CppTools::PrintError("Could not determine the range for RAB");
      }
      // checking pT bins consistency
      for (int i = 1; i <= distrSpectraPPVsPTStatErr->GetXaxis()->GetNbins(); i++)
      {
         if (fabs(distrSpectraPPVsPTStatErr->GetXaxis()->GetBinLowEdge(i) - 
                  pTBinRanges[pTBinMinRAB + i - 1]) > 1e-7 ||
             fabs(distrSpectraPPVsPTStatErr->GetXaxis()->GetBinUpEdge(i) - 
                  pTBinRanges[pTBinMinRAB + i]) > 1e-7)
         {
            CppTools::Print(distrSpectraPPVsPTStatErr->GetXaxis()->GetBinLowEdge(i),
                  pTBinRanges[pTBinMinRAB + i - 1],
             distrSpectraPPVsPTStatErr->GetXaxis()->GetBinUpEdge(i),
                  pTBinRanges[pTBinMinRAB + i]);
            CppTools::PrintError("pT bins mismatch between p+p and A+B "\
                                 "spectra for pT bin " + std::to_string(i));
         }
      }
   }

   // counter for centrality bins
   int iC = 0;
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

         const double methodPTMin = pTBinRanges[method["centrality_bin_parameters"][iC]["pt_bin_min"].as<int>()];
         const double methodPTMax = pTBinRanges[method["centrality_bin_parameters"][iC]["pt_bin_max"].as<int>() + 1];

         const std::string inputFileName = 
            "data/RawYields/" + runName + "/Resonance/" + std::to_string(taxiNumber) + 
            "_" + resonanceName + "_" + methodName + ".root";
         CppTools::CheckInputFile(inputFileName);
         TFile *inputFile = TFile::Open(inputFileName.c_str());

         TFile *inputFileLoosenedCuts = nullptr;
         TFile *inputFileTightenedCuts = nullptr;

         if (doCutsVarSys)
         {
            const std::string inputFileLoosenedCutsName = 
               "data/RawYields/" + runName + "/Resonance/" + 
               std::to_string(taxiNumberLoosenedCuts) + 
               "_" + resonanceName + "_" + methodName + ".root";
            const std::string inputFileTightenedCutsName = 
               "data/RawYields/" + runName + "/Resonance/" + 
               std::to_string(taxiNumberTightenedCuts) +
               "_" + resonanceName + "_" + methodName + ".root";

            CppTools::CheckInputFile(inputFileLoosenedCutsName);
            CppTools::CheckInputFile(inputFileTightenedCutsName);

            inputFileLoosenedCuts = TFile::Open(inputFileLoosenedCutsName.c_str());
            inputFileTightenedCuts = TFile::Open(inputFileTightenedCutsName.c_str());
         }

         // contains systematic uncertainties of pT scale variation
         TH1D sysPTScale("pT scale sys", "", pTNBins, &pTBinRanges[0]);
         // contains systematic uncertainties of acceptance variation
         TH1D sysAccVar("acc var sys", "", pTNBins, &pTBinRanges[0]);
         // contains systematic uncertainties of cuts variation
         TH1D sysCutsVar("cut var sys", "", pTNBins, &pTBinRanges[0]);
         // contains systematic uncertainties of raw yield extraction
         TH1D sysYieldExtr("yield extr sys", "", pTNBins, &pTBinRanges[0]);
         // contains full systematic uncertainty (A+B)
         TH1D sysFull("full sys", "", pTNBins, &pTBinRanges[0]);

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

         TH1D *rawYieldVsPTLoosenedCuts = nullptr;
         TH1D *rawYieldVsPTTightenedCuts = nullptr;

         if (doCutsVarSys)
         {
            rawYieldVsPTLoosenedCuts = static_cast<TH1D *>(inputFileLoosenedCuts->
               Get((centralityName + "/raw yield vs pT with stat errors").c_str()));
            rawYieldVsPTTightenedCuts = static_cast<TH1D *>(inputFileTightenedCuts->
               Get((centralityName + "/raw yield vs pT with stat errors").c_str()));

            if (!rawYieldVsPTLoosenedCuts) 
            {
               CppTools::PrintError("No raw yield distribution for loosened cuts "\
                                    "was found in file " + inputFileName + " for " + methodName + 
                                    " in centrality " + centralityName);
            }
            if (!rawYieldVsPTTightenedCuts) 
            {
               CppTools::PrintError("No raw yield distribution for tightened cuts "\
                                    "was found in file " + inputFileName + " for " + methodName + 
                                    " in centrality " + centralityName);
            }
         }

         TH1D *recEffVsPTStatErr = static_cast<TH1D *>
            (inputRecEffFile->Get((methodName + "/reconstruction efficiency "\
                                   "vs pT with stat errors").c_str()));
         TH1D *recEffVsPTSysErrAltPT = static_cast<TH1D *>
            (inputRecEffFile->Get((methodName + "/reconstruction efficiency "\
                                   "vs pT with sys errors, alt pT scale").c_str()));
         TH1D *recEffVsPTSysErrAccVar = static_cast<TH1D *>
            (inputRecEffFile->Get((methodName + "/reconstruction efficiency "\
                                   "vs pT with sys errors, acceptance variation").c_str()));
         TH1D *recEffVsPTSysErr = static_cast<TH1D *>
            (inputRecEffFile->Get((methodName + "/reconstruction efficiency "\
                                   "vs pT with sys errors").c_str()));

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
         if (!recEffVsPTSysErrAltPT)
         {
            CppTools::PrintError("No reconstruction efficiency with systematic errors for "\
                                 " alternative pT scale was found in file " + 
                                 inputRecEffFileName + " for " + methodName);
         }
         if (!recEffVsPTSysErrAccVar)
         {
            CppTools::PrintError("No reconstruction efficiency with systematic errors for "\
                                 " acceptance variation was found in file " + 
                                 inputRecEffFileName + " for " + methodName);
         }

         TH1D *recEffVsPTLoosenedCuts = nullptr;
         TH1D *recEffVsPTTightenedCuts = nullptr;

         if (doCutsVarSys)
         {
            // no systematics needed for cut variations
            recEffVsPTLoosenedCuts = static_cast<TH1D *>
               (inputRecEffFile->Get((methodName + "/reconstruction efficiency "\
                                      "vs pT, loosened cuts").c_str()));
            recEffVsPTTightenedCuts = static_cast<TH1D *>
               (inputRecEffFile->Get((methodName + "/reconstruction efficiency "\
                                      "vs pT, tightened cuts").c_str()));
            if (!recEffVsPTLoosenedCuts)
            {
               CppTools::PrintError("No reconstruction efficiency for loosened cuts "\
                                    "analysis was found in file " + inputRecEffFileName + 
                                    " for " + methodName);
            }
            if (!recEffVsPTTightenedCuts)
            {
               CppTools::PrintError("No reconstruction efficiency for tightened cuts "\
                                    "analysis was found in file " + inputRecEffFileName + 
                                    " for " + methodName);
            }
         }

         rawYieldVsPTStatErr->Divide(recEffVsPTStatErr);
         rawYieldVsPTStatErr->Scale(1./spectraNorm);

         if (doCutsVarSys)
         {
            rawYieldVsPTLoosenedCuts->Divide(recEffVsPTLoosenedCuts);
            rawYieldVsPTTightenedCuts->Divide(recEffVsPTTightenedCuts);
            rawYieldVsPTLoosenedCuts->Scale(1./spectraNorm);
            rawYieldVsPTTightenedCuts->Scale(1./spectraNorm);
         }

         TH1D resultVarCutsVarL("results variation from loosened cuts", 
                                "", pTNBins, &pTBinRanges[0]);
         TH1D resultVarCutsVarT("results variation from tightened cuts", 
                                "", pTNBins, &pTBinRanges[0]);

         // setting all systematics (cuts variation and yield extraction are corrected lower)
         for (unsigned int i = 1; i <= pTNBins; i++)
         {
            if (rawYieldVsPTStatErr->GetBinContent(i) < 1e-31) continue;

            sysPTScale.SetBinContent(i, recEffVsPTSysErrAltPT->GetBinError(i)/
                                     recEffVsPTSysErrAltPT->GetBinContent(i));
            sysAccVar.SetBinContent(i, recEffVsPTSysErrAccVar->GetBinError(i)/
                                     recEffVsPTSysErrAccVar->GetBinContent(i));
            sysYieldExtr.SetBinContent(i, rawYieldVsPTSysErr->GetBinError(i)/
                                       rawYieldVsPTSysErr->GetBinContent(i));

            double cutsVarSys = 0.;
            if (doCutsVarSys)
            {
               // relative variations of loosened and tightened cuts 
               // and their statistical uncertainties
               const double varL = (rawYieldVsPTLoosenedCuts->GetBinContent(i) - 
                                    rawYieldVsPTStatErr->GetBinContent(i))/
                                   rawYieldVsPTStatErr->GetBinContent(i);
               const double varT = (rawYieldVsPTTightenedCuts->GetBinContent(i) - 
                                    rawYieldVsPTStatErr->GetBinContent(i))/
                                   rawYieldVsPTStatErr->GetBinContent(i);
               const double varLErr = 
                  CppTools::UncertaintyProp(rawYieldVsPTStatErr->GetBinError(i)/
                                            rawYieldVsPTStatErr->GetBinContent(i),
                                            rawYieldVsPTLoosenedCuts->GetBinError(i)/
                                            rawYieldVsPTLoosenedCuts->GetBinContent(i))/2.;
               const double varTErr = 
                  CppTools::UncertaintyProp(rawYieldVsPTStatErr->GetBinError(i)/
                                            rawYieldVsPTStatErr->GetBinContent(i),
                                            rawYieldVsPTTightenedCuts->GetBinError(i)/
                                            rawYieldVsPTTightenedCuts->GetBinContent(i))/2.;

               cutsVarSys = CppTools::RMS(varL, varT);

               sysCutsVar.SetBinContent(i, cutsVarSys);
               sysCutsVar.SetBinError(i, CppTools::UncertaintyProp(varLErr, varTErr)/2.);

               resultVarCutsVarL.SetBinContent(i, varL);
               resultVarCutsVarT.SetBinContent(i, varT);

               resultVarCutsVarL.SetBinError(i, varLErr*varL);
               resultVarCutsVarT.SetBinError(i, varTErr*varT);
            }
         }
         // Cuts variations: setting histograms and drawing
         { 
            double minY = CppTools::Minimum(resultVarCutsVarL.GetMinimum(),
                                            resultVarCutsVarT.GetMinimum());
            double maxY = CppTools::Maximum(resultVarCutsVarL.GetMaximum(),
                                            resultVarCutsVarT.GetMaximum());

            if (minY > 0) minY /= 1.3;
            else (minY) *= 1.3;

            if (maxY > 0) maxY *= 1.3;
            else maxY /= 1.3;
            
            resultVarCutsVarL.SetMinimum(minY);
            resultVarCutsVarL.SetMaximum(maxY);

            resultVarCutsVarL.GetXaxis()->
               SetRange(resultVarCutsVarL.GetXaxis()->FindBin(methodPTMin + 1e-7),
                        resultVarCutsVarL.GetXaxis()->FindBin(methodPTMax - 1e-7));

            resultVarCutsVarL.SetLineWidth(2);
            resultVarCutsVarT.SetLineWidth(2);

            resultVarCutsVarL.SetLineColorAlpha(kAzure - 3, 0.8);
            resultVarCutsVarT.SetLineColorAlpha(kRed - 3, 0.8);

            TCanvas canv("cuts var canv", "", 800, 800);

            canv.SetFillStyle(4000);
            canv.SetFrameFillColor(0);
            canv.SetFrameFillStyle(0);
            canv.SetFrameBorderMode(0);

            gPad->SetRightMargin(0.035); gPad->SetTopMargin(0.03); 
            gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.112);

            ROOTTools::DrawFrame(&resultVarCutsVarL, "", "#it{p}_{T} [GeV/#it{c}]", 
                                 "Var(#it{Y})", 1., 1.5);

            resultVarCutsVarT.Draw("SAME");

            TLegend legend(0.15, 0.85, 0.9, 0.95);
            legend.SetLineColorAlpha(0, 0.);
            legend.SetFillColorAlpha(0, 0.);
            legend.SetNColumns(2);

            legend.AddEntry(&resultVarCutsVarL, "Loosened cuts");
            legend.AddEntry(&resultVarCutsVarT, "Tightened cuts");

            legend.Draw();

            if (minY < 0. && maxY > 0.)
            {
               TLine line(methodPTMin, 0., methodPTMax, 0.);
               line.SetLineColorAlpha(kBlack, 0.5);
               line.SetLineStyle(2);
               line.SetLineWidth(4);
               line.Draw();
            }

            ROOTTools::PrintCanvas(&canv, outputDirSys + "/CutsVar_" +  resonanceName + 
                                   "_" + methodName + "_" + centralityName);
         }
         // Systematics of cuts variation: correcting systematics, setting histograms, fits, and drawing
         if (doCutsVarSys)
         { 
            sysCutsVar.GetXaxis()->
               SetRange(sysCutsVar.GetXaxis()->FindBin(methodPTMin + 1e-7),
                        sysCutsVar.GetXaxis()->FindBin(methodPTMax - 1e-7));

            sysCutsVar.SetLineWidth(2);

            sysCutsVar.SetLineColorAlpha(kBlack, 0.9);

            TF1 fit("cuts var sys fit", "pol2");

            fit.SetLineColor(kRed - 3);
            fit.SetLineWidth(4);
            fit.SetLineStyle(2);

            fit.SetRange(methodPTMin/1.05, methodPTMax*1.05);

            sysCutsVar.Fit(&fit, "QMN");

            TCanvas canv("cuts var systematics canv", "", 800, 800);

            canv.SetFillStyle(4000);
            canv.SetFrameFillColor(0);
            canv.SetFrameFillStyle(0);
            canv.SetFrameBorderMode(0);

            gPad->SetRightMargin(0.035); gPad->SetTopMargin(0.03); 
            gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.112);

            ROOTTools::DrawFrame(&sysCutsVar, "", "#it{p}_{T} [GeV/#it{c}]", 
                                 "Relative uncertainty", 1., 1.5);

            fit.Draw("SAME");

            ROOTTools::PrintCanvas(&canv, outputDirSys + "/CutsVarSys_" +  resonanceName + 
                                   "_" + methodName + "_" + centralityName);

            for (int i = 1; i < sysCutsVar.GetXaxis()->GetNbins(); i++)
            {
               const double val = fit.Eval(sysCutsVar.GetXaxis()->GetBinCenter(i));
               if (val > 0.)
               {
                  sysCutsVar.SetBinContent(i, fit.Eval(sysCutsVar.GetXaxis()->GetBinCenter(i)));
               }
               else sysCutsVar.SetBinContent(i, 0.);
            }
         }
         // Finishing setting systematics
         for (unsigned int i = 1; i <= pTNBins; i++)
         {
            if (rawYieldVsPTStatErr->GetBinContent(i) < 1e-31) continue;
            // full relative systematic uncertainty (A+B types)
            const double sys = 
               CppTools::UncertaintyProp(sysPTScale.GetBinContent(i), 
                                         sysAccVar.GetBinContent(i),
                                         sysCutsVar.GetBinContent(i),
                                         sysYieldExtr.GetBinContent(i));

            sysFull.SetBinContent(i, sys);

            rawYieldVsPTSysErr->SetBinContent(i, rawYieldVsPTStatErr->GetBinContent(i));
            rawYieldVsPTSysErr->SetBinError(i, sys*rawYieldVsPTSysErr->GetBinContent(i));

            yMin = CppTools::Minimum(yMin, rawYieldVsPTStatErr->GetBinContent(i));
            yMax = CppTools::Maximum(yMax, rawYieldVsPTStatErr->GetBinContent(i));
         }

         spectrasVsPTStatErr.emplace_back(rawYieldVsPTStatErr);
         spectrasVsPTSysErr.emplace_back(rawYieldVsPTSysErr);

         // All systematics: setting histograms and drawing
         { 
            sysFull.GetXaxis()->
               SetRange(sysFull.GetXaxis()->FindBin(methodPTMin + 1e-7),
                        sysFull.GetXaxis()->FindBin(methodPTMax - 1e-7));

            sysFull.SetMaximum(sysFull.GetMaximum()*1.2);
            sysFull.SetMinimum(0.001);

            sysPTScale.SetLineColorAlpha(kP6Yellow, 0.9);
            sysAccVar.SetLineColorAlpha(kP6Red, 0.9);
            sysYieldExtr.SetLineColorAlpha(kP6Blue, 0.9);
            sysFull.SetLineColor(kP6Gray);

            sysPTScale.SetLineWidth(2);
            sysAccVar.SetLineWidth(2);
            sysYieldExtr.SetLineWidth(2);
            sysFull.SetLineWidth(4);

            sysPTScale.SetLineStyle(2);
            sysAccVar.SetLineStyle(9);
            sysYieldExtr.SetLineStyle(8);

            TCanvas canv("sys canv", "", 800, 800);

            canv.SetFillStyle(4000);
            canv.SetFrameFillColor(0);
            canv.SetFrameFillStyle(0);
            canv.SetFrameBorderMode(0);

            gPad->SetRightMargin(0.035); gPad->SetTopMargin(0.03); 
            gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.112);

            ROOTTools::DrawFrame(&sysFull, "", "#it{p}_{T} [GeV/#it{c}]", "Relative uncertainty", 1., 1.5);

            sysPTScale.Draw("SAME");
            sysAccVar.Draw("SAME");
            sysYieldExtr.Draw("SAME");

            if (doCutsVarSys)
            {
               sysCutsVar.Sumw2(false);
               sysCutsVar.SetLineStyle(7);
               sysCutsVar.SetLineColorAlpha(kP6Violet, 0.9);
               sysCutsVar.SetLineWidth(2);
               sysCutsVar.Draw("SAME");
            }

            TLegend legend(0.15, 0.85, 0.9, 0.95);
            legend.SetLineColorAlpha(0, 0.);
            legend.SetFillColorAlpha(0, 0.);
            legend.SetNColumns(3);

            legend.AddEntry(&sysPTScale, "#it{p}_{T} scale");
            legend.AddEntry(&sysAccVar, "Acceptance");
            legend.AddEntry(&sysCutsVar, "Cuts");
            legend.AddEntry(&sysYieldExtr, "Yield extraction");
            legend.AddEntry(&sysFull, "Full");

            legend.Draw();

            ROOTTools::PrintCanvas(&canv, outputDirSys + "/" +  resonanceName + 
                                   "_" + methodName + "_" + centralityName);
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

      TF1 tsallisFit("MB spectra fit", 
                     "0.5/pi*[0]*([1] - 1.)*([1] - 2.)/([2] + [3]*([1] - 1.))/"\
                     "([2] + [3])*([2] + sqrt(x^2 + [3]^2)/([2] + [3]))^(-[1])");
      tsallisFit.SetParameters(1., 2.5, 10.);
      tsallisFit.SetParLimits(1, 2., 30.);
      tsallisFit.FixParameter(3, resonanceMass);

      tsallisFit.SetRange(xMin/1.05, xMax*1.05);

      tsallisFit.SetLineStyle(2);
      tsallisFit.SetLineWidth(4);
      tsallisFit.SetLineColorAlpha(kBlack, 0.5);

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

      ROOTTools::
         DrawFrame(xMin - 0.1, yMin/5., xMax + 0.1, yMax*5., "", "#it{p}_{T} [GeV/#it{c}]", 
                   "1/(2#pi#it{p}_{T}) #it{d}^{2} #it{N}/#it{dp}_{T}/#it{dy} [(GeV/#it{c})^{-2}]");

      canvSpectra.SetFillStyle(4000);
      canvSpectra.SetFrameFillColor(0);
      canvSpectra.SetFrameFillStyle(0);
      canvSpectra.SetFrameBorderMode(0);

      tsallisFit.Draw("SAME");
      distrResultingSpectraVsPTStatErr.Draw("SAME");
      distrResultingSpectraVsPTSysErr.Draw("SAME E2");

      distrResultingSpectraVsPTStatErr.Clone()->Write();
      distrResultingSpectraVsPTSysErr.Clone()->Write();

      ROOTTools::PrintCanvas(&canvSpectra, outputDirResults + "/" + resonanceName + 
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

            ratioMin = 
               CppTools::Minimum(ratioMin, spectraRatiosVsPTStatErr[i]->GetBinContent(j) - 
                                 CppTools::Maximum(spectraRatiosVsPTStatErr[i]->GetBinError(j),
                                                   spectraRatiosVsPTSysErr[i]->GetBinError(j)));
            ratioMax = 
               CppTools::Maximum(ratioMax, spectraRatiosVsPTStatErr[i]->GetBinContent(j) +
                                 CppTools::Maximum(spectraRatiosVsPTStatErr[i]->GetBinError(j),
                                                   spectraRatiosVsPTSysErr[i]->GetBinError(j)));
         }
      }

      tsallisFit.SetLineColor(kRed - 3);

      TCanvas canvAllSpectra("all spectra canv", "", 800, 1000);

      canvAllSpectra.SetFillStyle(4000);
      canvAllSpectra.SetFrameFillColor(0);
      canvAllSpectra.SetFrameFillStyle(0);
      canvAllSpectra.SetFrameBorderMode(0);

      canvAllSpectra.Divide(1, 2, 0., 0.);

      canvAllSpectra.cd(1);

      gPad->SetLogy();

      gPad->SetPad(0., 0.5, 1., 1.);
      gPad->SetRightMargin(0.002); gPad->SetTopMargin(0.002); 
      gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.);

      ROOTTools::
         DrawFrame(xMin - 0.1, yMin/5., xMax + 0.1, yMax*5., "", "", 
                   "1/(2#pi#it{p}_{T}) #it{d}^{2} #it{N}/#it{dp}_{T}/#it{dy} [(GeV/#it{c})^{-2}]", 0., 0.95, 0.07, 0.07);

      tsallisFit.Draw("SAME");

      TLatex tlText;

      tlText.SetTextFont(52);
      tlText.SetTextSize(0.1);
      tlText.DrawLatexNDC(0.2, 0.1, centralityNameTex.c_str());

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

      gPad->SetPad(0., 0., 1., 0.5);
      gPad->SetRightMargin(0.002); gPad->SetTopMargin(0.); 
      gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.16);

      ROOTTools::DrawFrame(xMin - 0.1, ratioMin/1.1, xMax + 0.1, ratioMax*1.1, 
                           "", "#it{p}_{T} [GeV/#it{c}]", "Data/Fit", 1., 0.95, 0.07, 0.07);

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

      ROOTTools::PrintCanvas(&canvAllSpectra, outputDirResults + "/" + resonanceName + 
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
               push_back(new TH1D(("rab stat " + std::to_string(i)).c_str(), "", 
                                   pTBinMaxRAB - pTBinMinRAB + 1,
                                   &pTBinRanges[pTBinMinRAB]));
            distrRABsVsPTSysErr.
               push_back(new TH1D(("rab sys " + std::to_string(i)).c_str(), "", 
                                   pTBinMaxRAB - pTBinMinRAB + 1,
                                   &pTBinRanges[pTBinMinRAB]));

            for (int j = 1; j <= distrRABsVsPTStatErr.back()->GetXaxis()->GetNbins(); j++)
            {
               distrRABsVsPTStatErr.back()->
                  SetBinContent(j, spectrasVsPTStatErr[i]->GetBinContent(pTBinMinRAB + j));
               distrRABsVsPTStatErr.back()->
                  SetBinError(j, spectrasVsPTStatErr[i]->GetBinError(pTBinMinRAB + j));
               distrRABsVsPTSysErr.back()->
                  SetBinContent(j, spectrasVsPTSysErr[i]->GetBinContent(pTBinMinRAB + j));
               distrRABsVsPTSysErr.back()->
                  SetBinError(j, spectrasVsPTSysErr[i]->GetBinError(pTBinMinRAB + j));
            }

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

         TCanvas canvAllRAB("rab canv", "", 800, 800);

         canvAllRAB.SetFillStyle(4000);
         canvAllRAB.SetFrameFillColor(0);
         canvAllRAB.SetFrameFillStyle(0);
         canvAllRAB.SetFrameBorderMode(0);

         gPad->SetRightMargin(0.002); gPad->SetTopMargin(0.002); 
         gPad->SetLeftMargin(0.1); gPad->SetBottomMargin(0.112);

         ROOTTools::DrawFrame(xMin - 0.1, 0., xMax + 0.1, maxRAB, 
                              "", "#it{p}_{T} [GeV/#it{c}]", "#it{R}_{AB}", 1., 0.95);

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

         ROOTTools::PrintCanvas(&canvAllRAB, outputDirResults + "/" + resonanceName + 
                                "_RAB_all_methods_" + centralityName);

         TH1D distrResultingRABVsPTStatErr("resulting rab stat", "", 
                                           pTBinMaxRAB - pTBinMinRAB + 1,
                                           &pTBinRanges[pTBinMinRAB]);
         TH1D distrResultingRABVsPTSysErr("resulting rab sys", "", 
                                           pTBinMaxRAB - pTBinMinRAB + 1,
                                           &pTBinRanges[pTBinMinRAB]);

         for (int i = 1; i <= distrResultingRABVsPTStatErr.GetXaxis()->GetNbins(); i++)
         {
            distrResultingRABVsPTStatErr.
               SetBinContent(i, distrResultingSpectraVsPTStatErr.GetBinContent(pTBinMinRAB + i));
            distrResultingRABVsPTStatErr.
               SetBinError(i, distrResultingSpectraVsPTStatErr.GetBinError(pTBinMinRAB + i));
            distrResultingRABVsPTSysErr.
               SetBinContent(i, distrResultingSpectraVsPTSysErr.GetBinContent(pTBinMinRAB + i));
            distrResultingRABVsPTSysErr.
               SetBinError(i, distrResultingSpectraVsPTSysErr.GetBinError(pTBinMinRAB + i));
         }

         distrResultingRABVsPTStatErr.Divide(distrSpectraPPVsPTStatErr);
         distrResultingRABVsPTSysErr.Divide(distrSpectraPPVsPTSysErr);

         distrResultingRABVsPTStatErr.Scale(scaleRAB);
         distrResultingRABVsPTSysErr.Scale(scaleRAB);

         for (int i = 1; i <= distrResultingRABVsPTStatErr.GetXaxis()->GetNbins(); i++)
         {
            if (distrResultingRABVsPTStatErr.GetBinContent(i) < 1e-7)
            {
               distrResultingRABVsPTStatErr.SetBinContent(i, -1e-7);
               distrResultingRABVsPTSysErr.SetBinContent(i, -1e-7);
            }
         }

         distrResultingRABVsPTStatErr.Write("RAB vs pT with stat errors");
         distrResultingRABVsPTSysErr.Write("RAB vs pT with sys errors");
      }
      iC++;

      if (centralityBin["is_mb"].as<bool>())
      {
         resultsOutputFile->cd();
         tsallisFit.Write();
      }
   }

   resultsOutputFile->Close();

   CppTools::PrintInfo("Results were succesfully evaluated");
   CppTools::PrintInfo("Spectras were written in " + resultsOutputFileName);
   CppTools::PrintInfo("RAB and spectra pictures were written in " + 
                       outputDirResults + " directory");
   CppTools::PrintInfo("Systematics were written in " + outputDirSys + " directory");

   return 0;
}

#endif /* ESTIMATE_RESULTS_CPP */
