/** 
 *  @file   EstimateGaussianBroadening.cpp 
 *  @brief  Contains realisations of functions that are used for estimation gaussian broadening parameter sigma used as a parameter in gaus which is convoluted with breit-wigner
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysis).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef ESTIMATE_GAUSSIAN_BROADENING_CPP
#define ESTIMATE_GAUSSIAN_BROADENING_CPP

#include "EstimateGaussianBroadening.hpp"

using namespace EstimateGaussianBroadening;

int main(int argc, char **argv)
{
   if (argc < 2 || argc > 3) 
   {
      std::string errMsg = "Expected 1-2 parameters while " + std::to_string(argc) + " ";
      errMsg += "parameter(s) were provided \n Usage: bin/EstimateGaussianBroadening ";
      errMsg += "inputYAMLName numberOfThreads=std::thread::hardware_concurrency()";
      CppTools::PrintError(errMsg);
   }
 
   CppTools::CheckInputFile(argv[1]);
 
   if (argc == 2) ROOT::EnableImplicitMT(std::thread::hardware_concurrency());
   else ROOT::EnableImplicitMT(std::stoi(argv[2]));

   gStyle->SetOptStat(0);
   gErrorIgnoreLevel = kWarning;

   inputYAMLResonance.OpenFile(argv[1]);
   inputYAMLResonance.CheckStatus("resonance");

   runName = inputYAMLResonance["run_name"].as<std::string>();

   inputYAMLMain.OpenFile("input/" + runName + "/main.yaml");
   inputYAMLMain.CheckStatus("main");

   outputDir = "output/GaussianBroadening/" + runName;
   system(("mkdir -p " + outputDir).c_str());

   nameResonance = inputYAMLResonance["name"].as<std::string>();
   massResonance = inputYAMLResonance["mass"].as<double>();

   const std::string inputFileName = "data/PostSim/" + runName + "/WidthlessResonance/" + 
                                     nameResonance + ".root";

   CppTools::CheckInputFile(inputFileName);
   inputFile = TFile::Open(inputFileName.c_str(), "READ");

   distr2DMInv = static_cast<TH2F *>(inputFile->Get("M_inv: NoPID"));

   const unsigned int pTNBins = inputYAMLResonance["pt_bins"].size() - 1;
   const double pTMin = inputYAMLResonance["pt_bins"][0]["min"].as<double>()/1.1;
   const double pTMax = inputYAMLResonance["pt_bins"][pTNBins]["max"].as<double>()*1.1;

   //text.SetTextFont(43);
   //text.SetTextSize(50);
   texText.SetTextFont(43);
   texText.SetTextSize(50);

   ProgressBar pBar("FANCY", "", PBarColor::BOLD_CYAN);

   for (int i = distr2DMInv->GetXaxis()->FindBin(pTMin); 
        i < distr2DMInv->GetXaxis()->FindBin(pTMax); i += 2)
   {
      pBar.Print(static_cast<double>(i)/
                 static_cast<double>(distr2DMInv->GetXaxis()->FindBin(pTMax) - 
                                     distr2DMInv->GetXaxis()->FindBin(pTMin)));
      PerformMInvFit(i, i + 1);
   }

   pBar.Finish();

   TF1 fitSigmas("gaussian broadening sigma fit", fitSigmasFormula.c_str());
   fitSigmas.SetRange(pTMin, pTMax);

   fitSigmas.SetLineWidth(3);
   fitSigmas.SetLineStyle(2);
   fitSigmas.SetLineColor(kGray + 3);

   grSigmas.Fit(&fitSigmas, "RQMNB");

   grSigmas.SetMarkerStyle(20);
   grSigmas.SetMarkerSize(1.4);
   grSigmas.SetMarkerColor(kRed - 3);

   grSigmas.SetLineWidth(2);
   grSigmas.SetLineColor(kRed + 2);

   TCanvas canv("canv", "", 800, 800);

   gPad->SetRightMargin(0.01);
   gPad->SetTopMargin(0.02);
   gPad->SetLeftMargin(0.17);
   gPad->SetBottomMargin(0.112);

   ROOTTools::DrawFrame(pTMin, 0., pTMax, TMath::MaxElement(grSigmas.GetN(), grSigmas.GetY())*1.2,
                        "", "#it{p}_{T} [GeV/c]", "#it{#sigma} [GeV/c^{2}]", 1., 1.8);

   fitSigmas.Draw("SAME");
   grSigmas.Clone()->Draw("P");

   ROOTTools::PrintCanvas(&canv, outputDir + "/" + nameResonance + "_sigmas");

   const std::string parametersOutputDir = "data/Parameters/GaussianBroadening/" + runName;
   system(("mkdir -p " + parametersOutputDir).c_str());

   TFile parametersOutput((parametersOutputDir + "/" + nameResonance + ".root").c_str(), "RECREATE");
   parametersOutput.cd();

   fitSigmas.Write();

   CppTools::PrintInfo("EstimateGaussianBroadening executable has finished running succesfully");
}

void EstimateGaussianBroadening::PerformMInvFit(const int pTBinMin, const int pTBinMax)
{
   TH1D *distrMInv = distr2DMInv->ProjectionY("proj", pTBinMin, pTBinMax);

   // fit for resonance+bg approximation
   TF1 fit("resonance + bg fit", "gaus(0) + gaus(3)");
   // fit for resonance approximation
   TF1 fitResonance("resonance fit", "gaus");
   // fit for bg approximation
   TF1 fitBG("bg fit", "gaus");

   const double maxBinVal = distrMInv->GetBinContent(distrMInv->GetMaximumBin());

   fit.SetParameters(maxBinVal, massResonance, 5e-3, maxBinVal/20., massResonance, 0.2);

   fit.SetParLimits(0, maxBinVal/3., maxBinVal);
   fit.SetParLimits(1, massResonance/1.02, massResonance*1.02);
   fit.SetParLimits(2, 1e-3, 2e-2);
   fit.SetParLimits(3, 0., maxBinVal/5.);
   fit.SetParLimits(4, 0., massResonance*10.);
   fit.SetParLimits(5, 5e-2, 1.);

   fit.SetRange(massResonance - 1e-2, massResonance + 1e-2);
   distrMInv->Fit(&fit, "RQMNBL");

   for (unsigned int j = 1; j <= fitNTries; j++)
   {
      fit.SetParLimits(1, fit.GetParameter(1) - 1e-2/static_cast<double>(j*j*j), 
                       fit.GetParameter(1) + 1e-2/static_cast<double>(j*j*j));
      fit.SetParLimits(2, fit.GetParameter(2)/(1. + 2./static_cast<double>(j*j*j)),
                       fit.GetParameter(2)*(1. + 2./static_cast<double>(j*j*j)));
      fit.SetParLimits(4, fit.GetParameter(4)/(1. + 2./static_cast<double>(j*j*j)),
                       fit.GetParameter(4)*(1. + 2./static_cast<double>(j*j*j)));
      fit.SetParLimits(5, fit.GetParameter(2)*5.,
                       fit.GetParameter(5)*(1. + 2./static_cast<double>(j*j*j)));

      fit.SetRange(fit.GetParameter(1) - fit.GetParameter(2)*10., 
                  fit.GetParameter(1) + fit.GetParameter(2)*10.);

      distrMInv->Fit(&fit, "RQMNBL");
      if (j == fitNTries) distrMInv->Fit(&fit, "RQMNBLE");
   }

   for (int j = 0; j < fitResonance.GetNpar(); j++)
   {
      fitResonance.SetParameter(j, fit.GetParameter(j));
   }

   for (int j = 0; j < fitBG.GetNpar(); j++)
   {
      fitBG.SetParameter(j, fit.GetParameter(j + fitResonance.GetNpar()));
   }

   fitResonance.SetRange(fit.GetParameter(1) - fit.GetParameter(2)*10., 
                         fit.GetParameter(1) + fit.GetParameter(2)*10.);
   fitBG.SetRange(fit.GetParameter(1) - fit.GetParameter(2)*10., 
                  fit.GetParameter(1) + fit.GetParameter(2)*10.);
   distrMInv->GetXaxis()->
      SetRange(distrMInv->GetXaxis()->FindBin(fit.GetParameter(1) - fit.GetParameter(2)*10.), 
               distrMInv->GetXaxis()->FindBin(fit.GetParameter(1) + fit.GetParameter(2)*10.));

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

   TCanvas canv("canv", "", 800, 800);

   gPad->SetRightMargin(0.03);
   gPad->SetTopMargin(0.02);
   gPad->SetLeftMargin(0.142);
   gPad->SetBottomMargin(0.112);

   const double pTMin = distr2DMInv->GetXaxis()->GetBinLowEdge(pTBinMin);
   const double pTMax = distr2DMInv->GetXaxis()->GetBinUpEdge(pTBinMax);

   ROOTTools::DrawFrame(distrMInv, "", "#it{M}_{inv} [GeV/c^{2}]", "Weighted counts");

   texText.DrawLatexNDC(0.17, 0.9, (CppTools::DtoStr(pTMin, 1) + " < #it{p}_{T} < " + 
                        CppTools::DtoStr(pTMax, 1)).c_str());

   fitBG.Draw("SAME");
   fitResonance.Draw("SAME");
   fit.Draw("SAME");

   ROOTTools::PrintCanvas(&canv, outputDir + "/" + nameResonance + "_" + 
                          CppTools::DtoStr(pTMin, 1) + "-" + CppTools::DtoStr(pTMax, 1));

   grSigmas.AddPoint(CppTools::Average(pTMin, pTMax), fit.GetParameter(2));
   grSigmas.SetPointError(grSigmas.GetN() - 1, 0., fit.GetParError(2));
}

#endif /* ESTIMATE_ESTIMATE_GAUSSIAN_BROADENING_CPP */
