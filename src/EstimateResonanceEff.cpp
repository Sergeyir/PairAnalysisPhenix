/** 
 *  @file   EstimateResonanceEff.cpp 
 *  @brief  Contains realisations of functions that are used for estimation of resonance reconstruction efficiency with he use of the data from MC
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysis).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef ESTIMATE_RESONANCE_EFF_CPP
#define ESTIMATE_RESONANCE_EFF_CPP

#include "EstimateResonanceEff.hpp"

using namespace EstimateResonanceEff;

int main(int argc, char **argv)
{
   if (argc < 2 || argc > 3) 
   {
      std::string errMsg = "Expected 1-2 parameters while " + std::to_string(argc - 1) + " ";
      errMsg += "parameter(s) were provided \n Usage: bin/EstimateResonanceEff ";
      errMsg += "inputYAMLName numberOfThreads=std::thread::hardware_concurrency()";
      CppTools::PrintError(errMsg);
   }
 
   CppTools::CheckInputFile(argv[1]);
 
   if (argc == 2) ROOT::EnableImplicitMT(std::thread::hardware_concurrency());
   else ROOT::EnableImplicitMT(std::stoi(argv[2]));

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

   SetGaussianBroadeningFunction();

   sigmalizedYieldExtractionRange = 
      inputYAMLResonance["sigmalized_yield_extraction_range"].as<double>();

   const std::string inputDir = "data/PostSim/" + runName + "/Resonance/";

   inputFileName = inputDir + resonanceName + ".root";
   CppTools::CheckInputFile(inputFileName);

   inputFile = TFile::Open(inputFileName.c_str(), "READ");

   for (auto &file : std::filesystem::directory_iterator(inputDir))
   {
      const std::string fileName = static_cast<std::string>(file.path());

      if (fileName == inputFileName) continue;

      if (std::regex_match(fileName, std::regex("(.*)" + resonanceName + 
                                                "_pTScale_([0-9\\.]*)\\.root")))
      {
         CppTools::PrintInfo("Found alternative simulation file " + fileName +   
                             " for pT scale systematic uncertainties evaluation");

         altPTScaleSimInputFileNames.emplace_back(fileName);
         altPTScaleSimInputFiles.emplace_back(TFile::Open(fileName.c_str()));
      }
   }

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

   // 17 different pairs selections methods
   numberOfIterations = pTNBins*inputYAMLResonance["pair_selection_methods"].size();

   const std::string parametersOutputDir = "data/Parameters/ResonanceEff/" + runName;
   outputFile = TFile::Open((parametersOutputDir + "/" + resonanceName + 
                             ".root").c_str(), "RECREATE");

   // performing fits for each pair selection method
   for (const auto& method : inputYAMLResonance["pair_selection_methods"])
   {
      PerformMInvFitsForMethod(method["name"].as<std::string>());
   }

   outputFile->Close();
   pBar.Finish();

   CppTools::PrintInfo("EstimateResonanceEff executable has finished running succesfully");

   return 0;
}

void EstimateResonanceEff::PerformMInvFitsForMethod(const std::string& methodName)
{
   const std::string outputDir = "output/ResonanceEff/" + runName + "/" + methodName;
   std::filesystem::create_directories(outputDir);

   outputFile->mkdir(methodName.c_str());
   outputFile->cd(methodName.c_str());

   text.SetTextAngle(270.);

   TH1D distrMeansVsPT("means vs pT", "", pTNBins, &pTBinRanges[0]);
   TH1D distrGammasVsPT("gammas vs pT", "", pTNBins, &pTBinRanges[0]);
   TH1D distrRecEffVsPTStatErr("reconstruction efficiency vs pT with stat errors", 
                               "", pTNBins, &pTBinRanges[0]);
   TH1D distrRecEffVsPTSysErr("reconstruction efficiency vs pT with sys errors", 
                              "", pTNBins, &pTBinRanges[0]);

   std::vector<TH1D> distrAltSimPTScaleMeansVsPT;
   std::vector<TH1D> distrAltSimPTScaleGammasVsPT;
   std::vector<TH1D> distrAltSimPTScaleRecEffVsPT;

   for (unsigned int i = 0; i < altPTScaleSimInputFiles.size(); i++)
   {
      distrAltSimPTScaleMeansVsPT.emplace_back(("means vs pT " + std::to_string(i)).c_str(), 
                                               "", pTNBins, &pTBinRanges[0]);
      distrAltSimPTScaleGammasVsPT.emplace_back(("gammas vs pT " + std::to_string(i)).c_str(), 
                                                 "", pTNBins, &pTBinRanges[0]);
      distrAltSimPTScaleRecEffVsPT.emplace_back(("rec eff vs pT " + std::to_string(i)).c_str(), 
                                                 "", pTNBins, &pTBinRanges[0]);
   }

   for (unsigned int i = 0; i < pTNBins; i++)
   {
      pBar.Print(static_cast<double>(numberOfCalls)/static_cast<double>(numberOfIterations));

      PerformMInvFit(i, methodName, inputFile, 
                     distrRecEffVsPTStatErr, distrMeansVsPT, distrGammasVsPT, 
                     outputDir + "/" + resonanceName + "_" + 
                     CppTools::DtoStr(pTBinRanges[i], 1) + "-" + 
                     CppTools::DtoStr(pTBinRanges[i + 1], 1));

      double recEffSysErr = 0.;

      for (unsigned j = 0; j < altPTScaleSimInputFiles.size(); j++)
      {
         PerformMInvFit(i, methodName, altPTScaleSimInputFiles[j], 
                        distrAltSimPTScaleRecEffVsPT[j], distrAltSimPTScaleMeansVsPT[j], 
                        distrAltSimPTScaleGammasVsPT[j], "");

         recEffSysErr += pow(distrRecEffVsPTStatErr.GetBinContent(i + 1) - 
                             distrAltSimPTScaleRecEffVsPT[j].GetBinContent(i + 1), 2);
      }

      if (altPTScaleSimInputFiles.size() != 0)
      {
         recEffSysErr /= static_cast<double>(altPTScaleSimInputFiles.size());
         recEffSysErr = sqrt(recEffSysErr);
      }

      distrRecEffVsPTSysErr.SetBinContent(i + 1, distrRecEffVsPTStatErr.GetBinContent(i + 1));
      distrRecEffVsPTSysErr.SetBinError(i + 1, recEffSysErr);

      numberOfCalls++;
   }

   text.SetTextAngle(0.);

   distrMeansVsPT.SetLineColor(kRed - 2);
   distrMeansVsPT.SetMarkerColorAlpha(0, 0.);
   distrMeansVsPT.SetLineWidth(4);

   distrGammasVsPT.SetLineColor(kRed - 2);
   distrGammasVsPT.SetMarkerColorAlpha(0, 0.);
   distrGammasVsPT.SetLineWidth(4);

   distrRecEffVsPTStatErr.SetLineColor(kRed - 2);
   distrRecEffVsPTStatErr.SetMarkerColorAlpha(0, 0.);
   distrRecEffVsPTStatErr.SetLineWidth(3);

   distrRecEffVsPTSysErr.SetFillStyle(1001);
   distrRecEffVsPTSysErr.SetFillColorAlpha(kRed - 2, 0.5);


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

   text.DrawTextNDC(0.38, 0.9, ("MC " + methodName).c_str());
   massResonancePDG.Draw();

   ROOTTools::PrintCanvas(&canvMeansVsPT, outputDir + "/" + resonanceName + "_means");

   TCanvas canvGammasVsPT("canv gammas vs pT", "", 800, 800);

   gPad->SetRightMargin(0.03);
   gPad->SetTopMargin(0.02);
   gPad->SetLeftMargin(0.172);
   gPad->SetBottomMargin(0.112);

   ROOTTools::DrawFrame(&distrGammasVsPT, "", "#it{p}_{T} [GeV/#it{c}]", 
                        "#it{#Gamma} [GeV/#it{c}^{2}]", 1., 1.82);

   text.DrawTextNDC(0.38, 0.9, ("MC " + methodName).c_str());
   gammaResonancePDG.Draw();

   ROOTTools::PrintCanvas(&canvGammasVsPT, outputDir + "/" + resonanceName + "_gammas");

   TCanvas canvRecEffVsPT("canv rec eff vs pT", "", 800, 800);

   gPad->SetLogy();

   gPad->SetRightMargin(0.03);
   gPad->SetTopMargin(0.02);
   gPad->SetLeftMargin(0.142);
   gPad->SetBottomMargin(0.112);

   ROOTTools::DrawFrame(&distrRecEffVsPTStatErr, "", "#it{p}_{T} [GeV/#it{c}]", 
                        "#it{#varepsilon}_{rec}");

   distrRecEffVsPTSysErr.Draw("SAME E2");

   text.DrawTextNDC(0.38, 0.9, ("MC " + methodName).c_str());

   ROOTTools::PrintCanvas(&canvRecEffVsPT, outputDir + "/" + resonanceName + "_rec_eff");

   distrMeansVsPT.Write();
   distrGammasVsPT.Write();
   distrRecEffVsPTStatErr.Write();
   distrRecEffVsPTSysErr.Write();

   pBar.Clear();
   CppTools::PrintInfo(methodName + " done");
   pBar.RePrint();
}

void EstimateResonanceEff::PerformMInvFit(const unsigned int pTBin, const std::string& methodName,
                                          TFile *file, TH1D& distrRecEffVsPT,
                                          TH1D& distrMeansVsPT, TH1D& distrGammasVsPT,
                                          const std::string& outputFileNameWithoutExt)
{
   TH1D *distrOrigUnscaledPT = static_cast<TH1D *>(file->Get("orig unscaled pT"));
   if (!distrOrigUnscaledPT) CppTools::PrintError("Original unscaled pT distribution was not found "\
                                                  " in file " + (std::string) file->GetName());
   TH1D *distrOrigPT = static_cast<TH1D *>(file->Get("orig pT"));
   if (!distrOrigPT) CppTools::PrintError("Original pT distribution was not found in file" + 
                                          (std::string) file->GetName());

   // number of generated resonance particles in the current pT bin
   // statistical uncertainty of this value is insignificant and therefore can be ignored
   const double numberOfGenerated = 
      distrOrigPT->Integral(distrOrigPT->GetXaxis()->FindBin(pTBinRanges[pTBin] + 1e-6), 
                            distrOrigPT->GetXaxis()->FindBin(pTBinRanges[pTBin + 1] - 1e-6));
   // relative uncertanty of the number of generated particles in the current pT bin
   const double numberOfGeneratedRelativeErr = 
      1./sqrt(distrOrigUnscaledPT->
              Integral(distrOrigPT->GetXaxis()->FindBin(pTBinRanges[pTBin] + 1e-6), 
                       distrOrigPT->GetXaxis()->FindBin(pTBinRanges[pTBin + 1] - 1e-6)));

   const std::string distrMInvVsPTName = "M_inv: " + methodName;
   TH2F *distrMInvVsPT = static_cast<TH2F *>(inputFile->Get(distrMInvVsPTName.c_str()));

   if (!distrMInvVsPT) CppTools::PrintError("Distribution named " + distrMInvVsPTName + "\" "\
                                            "was not found in file " + inputFileName);

   TH1D *distrMInv = distrMInvVsPT->
      ProjectionY("proj", distrMInvVsPT->GetXaxis()->FindBin(pTBinRanges[pTBin] + 1e-6),
                  distrMInvVsPT->GetXaxis()->FindBin(pTBinRanges[pTBin + 1] - 1e-6));

   if (distrMInv->GetEntries() == 0) return;

   // sigma of a gaus that is convoluted with Breit-Wigner
   const double gaussianBroadeningSigma = 
      gaussianBroadeningEstimatorFunc->Eval((pTBinRanges[pTBin] + pTBinRanges[pTBin + 1])/2.);

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
   fit.SetParLimits(2, gammaResonance/1.05, gammaResonance*1.05);
   //fit.FixParameter(2, gammaResonance);
   fit.SetParLimits(3, gaussianBroadeningSigma/1.05, gaussianBroadeningSigma*1.05);
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
      fit.SetParLimits(3, fit.GetParameter(3)/(1. + 0.05/static_cast<double>(j*j)),
                       fit.GetParameter(3)*(1. + 0.05/static_cast<double>(j*j)));
      fit.SetParLimits(5, fit.GetParameter(5)/(1. + 2./static_cast<double>(j*j)),
                       fit.GetParameter(5)*(1. + 2./static_cast<double>(j*j)));
      fit.SetParLimits(6, fit.GetParameter(6)/(1. + 1./static_cast<double>(j*j)),
                       fit.GetParameter(6)*(1. + 1./static_cast<double>(j*j)));

      fit.SetRange(fit.GetParameter(1) - (fit.GetParameter(2) + gaussianBroadeningSigma)*3., 
                   fit.GetParameter(1) + (fit.GetParameter(2) + gaussianBroadeningSigma)*3.);

      distrMInv->Fit(&fit, "RQMNBLC");
   }

   //distrMInv->Fit(&fit, "RQMNBLE");

   distrMeansVsPT.SetBinContent(pTBin + 1, fit.GetParameter(1));
   distrMeansVsPT.SetBinError(pTBin + 1, fit.GetParError(1));

   distrGammasVsPT.SetBinContent(pTBin + 1, fit.GetParameter(2));
   distrGammasVsPT.SetBinError(pTBin + 1, fit.GetParError(2));

   for (int j = 0; j < fitResonance.GetNpar(); j++)
   {
      fitResonance.SetParameter(j, fit.GetParameter(j));
   }

   for (int j = 0; j < fitBG.GetNpar(); j++)
   {
      fitBG.SetParameter(j, fit.GetParameter(j + fitResonance.GetNpar()));
   }

   const double lowIntegrationRange = 
      fit.GetParameter(1) - (fit.GetParameter(2) + fit.GetParameter(3))*
                             sigmalizedYieldExtractionRange;
   const double upIntegrationRange = 
      fit.GetParameter(1) + (fit.GetParameter(2) + fit.GetParameter(3))*
                             sigmalizedYieldExtractionRange;

   double recYieldErr;
   double recYield = GetYield(distrMInv, fitBG, lowIntegrationRange, 
                              upIntegrationRange, recYieldErr);

   distrRecEffVsPT.SetBinContent(pTBin + 1, recYield/numberOfGenerated);
   distrRecEffVsPT.SetBinError(pTBin + 1, CppTools::UncertaintyProp(recYieldErr/recYield, 
                                                                    numberOfGeneratedRelativeErr)*
                               recYield/numberOfGenerated);

   if (outputFileNameWithoutExt != "")
   {
      fitResonance.SetRange(fit.GetParameter(1) - 
                            (fit.GetParameter(2) + gaussianBroadeningSigma)*3., 
                            fit.GetParameter(1) + 
                            (fit.GetParameter(2) + gaussianBroadeningSigma)*3.);
      fitBG.SetRange(fit.GetParameter(1) - (fit.GetParameter(2) + gaussianBroadeningSigma)*3., 
                     fit.GetParameter(1) + (fit.GetParameter(2) + gaussianBroadeningSigma)*3.);

      distrMInv->GetXaxis()->
         SetRange(CppTools::Maximum(distrMInv->GetXaxis()->FindBin(fit.GetParameter(1) - 
                                                                   fit.GetParameter(2)*5. -
                                                                   gaussianBroadeningSigma*5.), 1), 
                  distrMInv->GetXaxis()->FindBin(fit.GetParameter(1) + fit.GetParameter(2)*5. +
                                                 gaussianBroadeningSigma*5.));

      fit.SetLineWidth(4);
      fitResonance.SetLineWidth(4);
      fitBG.SetLineWidth(4);
      distrMInv->SetLineWidth(2);

      fit.SetLineColorAlpha(kRed - 3, 0.8);
      fitResonance.SetLineColorAlpha(kAzure - 3, 0.8);
      fitBG.SetLineColorAlpha(kGreen - 3, 0.8);

      fitResonance.SetLineStyle(2);
      fitBG.SetLineStyle(7);

      distrMInv->SetLineColor(kBlack);
      distrMInv->SetMarkerColor(kBlack);

      TCanvas canvMInv("canv MInv", "", 800, 800);

      gPad->SetRightMargin(0.03);
      gPad->SetTopMargin(0.02);
      gPad->SetLeftMargin(0.173);
      gPad->SetBottomMargin(0.112);

      ROOTTools::DrawFrame(distrMInv, "", "#it{M}_{inv} [GeV/c^{2}]", "Weighted counts", 1., 1.9);

      text.DrawTextNDC(0.9, 0.95, ("MC " + methodName).c_str());
      texText.DrawLatexNDC(0.2, 0.9, (CppTools::DtoStr(pTBinRanges[pTBin], 1) + 
                                      " < #it{p}_{T} < " + 
                                      CppTools::DtoStr(pTBinRanges[pTBin + 1], 1)).c_str());
      texText.DrawLatexNDC(0.2, 0.83, 
                           ("#it{#chi}^{2}/NDF = " + 
                            CppTools::DtoStr(fit.GetChisquare()/fit.GetNDF(), 2)).c_str());

      fitBG.Draw("SAME");
      fitResonance.Draw("SAME");
      fit.Draw("SAME");

      ROOTTools::PrintCanvas(&canvMInv, outputFileNameWithoutExt);
   }
}

void EstimateResonanceEff::SetGaussianBroadeningFunction()
{
   const std::string inputFileName = "data/Parameters/GaussianBroadening/" + 
                                     runName + "/" + resonanceName + ".root";

   if (!std::filesystem::exists(inputFileName))
   {
      CppTools::PrintError(inputFileName + " does not exists. "\
                           "Run executable bin/EstimateGassianBroadening first");
   }

   gaussianBroadeningEstimatorFunc = 
      static_cast<TF1 *>(TFile::Open(inputFileName.c_str())->Get("gaussian broadening sigma fit"));
}

double EstimateResonanceEff::GetYield(TH1D *distrMInv, const TF1& funcBG, 
                                      const double xMin, const double xMax, double &err)
{
   // integral over the signal
   double integral = 0.;
   // integral over the background
   double integralBG = 0.;

   if (distrMInv->GetXaxis()->FindBin(xMin) < 1 || 
       distrMInv->GetXaxis()->FindBin(xMax) > distrMInv->GetXaxis()->GetNbins())
   {
      CppTools::PrintWarning("EstimateResonanceEff::GetYield: specified integration range is "\
                             "outside the histogram range; ignoring underflow and overflow bins");
   }

   for (int i = CppTools::Maximum(distrMInv->GetXaxis()->FindBin(xMin), 1); 
        i <= CppTools::Minimum(distrMInv->GetXaxis()->FindBin(xMax), 
                               distrMInv->GetXaxis()->GetNbins()); i++)
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

#endif /* ESTIMATE_RESONANCE_EFF_CPP */
