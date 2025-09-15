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

   distr2DInvM = static_cast<TH2F *>(inputFile->Get("M_inv: NoPID"));

   const unsigned int pTNBins = inputYAMLResonance["pt_bins"].size() - 1;
   const double pTMin = inputYAMLResonance["pt_bins"][0]["min"].as<double>()/1.1;
   const double pTMax = inputYAMLResonance["pt_bins"][pTNBins]["max"].as<double>()*1.1;

   ProgressBar pBar("FANCY", "", PBarColor::BOLD_CYAN);

   for (int i = distr2DInvM->GetXaxis()->FindBin(pTMin); 
        i < distr2DInvM->GetXaxis()->FindBin(pTMax); i += 2)
   {
      pBar.Print(static_cast<double>(i)/
                 static_cast<double>(distr2DInvM->GetXaxis()->FindBin(pTMax) - 
                                     distr2DInvM->GetXaxis()->FindBin(pTMin)));
      PerformInvMassFit(i, i + 1);
   }

   pBar.Finish();

   TF1 fitSigmas("fit for sigmas", fitSigmasFormula.c_str());
   fitSigmas.SetRange(pTMin, pTMax);

   fitSigmas.SetLineWidth(3);
   fitSigmas.SetLineStyle(2);
   fitSigmas.SetLineColor(kGray + 3);

   grSigmas.Fit(&fitSigmas, "RQMNB");

   grSigmas.SetMarkerStyle(20);
   grSigmas.SetMarkerSize(1.4);
   grSigmas.SetMarkerColor(kRed - 3);

   TCanvas canv("canv", "", 800, 800);

   gPad->SetRightMargin(0.01);
   gPad->SetTopMargin(0.01);
   gPad->SetLeftMargin(0.17);
   gPad->SetBottomMargin(0.112);

   ROOTTools::DrawFrame(pTMin, 0., pTMax, TMath::MaxElement(grSigmas.GetN(), grSigmas.GetY())*1.2,
                        "", "p_{T} [GeV/c]", "#sigma [GeV/c^{2}]", 1., 1.8);

   fitSigmas.Draw("SAME");
   grSigmas.Clone()->Draw("P");

   ROOTTools::PrintCanvas(&canv, outputDir + "/" + nameResonance + "_sigmas");

   system(("mkdir -p data/Parameters/GaussianBroadening/" + runName).c_str());

   std::ofstream parametersOutput("data/Parameters/GaussianBroadening/" + runName + 
                                  "/" + nameResonance + ".txt");

   parametersOutput << fitSigmasFormula << std::endl;
   parametersOutput << fitSigmas.GetNpar() << std::endl;
   for (int i = 0; i < fitSigmas.GetNpar() - 1; i++)
   {
      parametersOutput << fitSigmas.GetParameter(i) << " ";
   }
   parametersOutput << fitSigmas.GetParameter(fitSigmas.GetNpar() - 1);

   CppTools::PrintInfo("EstimateGaussianBroadening executable has finished running succesfully");
}

void EstimateGaussianBroadening::PerformInvMassFit(const int pTBinMin, const int pTBinMax)
{
   TH1D *distrInvM = distr2DInvM->ProjectionY("proj", pTBinMin, pTBinMax);

   // fit for resonance+bg approximation
   TF1 fit("resonance + bg fit", "gaus(0) + gaus(3)");
   // fit for resonance approximation
   TF1 fitResonance("resonance fit", "gaus");
   // fit for bg approximation
   TF1 fitBG("bg fit", "gaus");

   const double maxBinVal = distrInvM->GetBinContent(distrInvM->GetMaximumBin());

   fit.SetParameters(maxBinVal, massResonance, 5e-3, maxBinVal/20., massResonance, 0.2);

   fit.SetParLimits(0, maxBinVal/3., maxBinVal);
   fit.SetParLimits(1, massResonance/1.02, massResonance*1.02);
   fit.SetParLimits(2, 1e-3, 2e-2);
   fit.SetParLimits(3, 0., maxBinVal/5.);
   fit.SetParLimits(4, 0., massResonance*10.);
   fit.SetParLimits(5, 5e-2, 1.);

   fit.SetRange(massResonance - 1e-2, massResonance + 1e-2);
   distrInvM->Fit(&fit, "RQMNB");

   for (unsigned int i = 1; i <= fitNTries; i++)
   {
      fit.SetParLimits(1, fit.GetParameter(1) - 1e-2/static_cast<double>(i*i*i), 
                       fit.GetParameter(1) - 1e-2/static_cast<double>(i*i*i));
      fit.SetParLimits(2, fit.GetParameter(2)/(1. + 2./static_cast<double>(i*i*i)),
                       fit.GetParameter(2)*(1. + 2./static_cast<double>(i*i*i)));
      fit.SetParLimits(4, fit.GetParameter(4)/(1. + 2./static_cast<double>(i*i*i)),
                       fit.GetParameter(4)*(1. + 2./static_cast<double>(i*i*i)));
      fit.SetParLimits(5, fit.GetParameter(2)*5.,
                       fit.GetParameter(5)*(1. + 2./static_cast<double>(i*i*i)));

      fit.SetRange(fit.GetParameter(1) - fit.GetParameter(2)*10., 
                  fit.GetParameter(1) + fit.GetParameter(2)*10.);

      distrInvM->Fit(&fit, "RQMNB");
   }

   for (int i = 0; i < fitResonance.GetNpar(); i++)
   {
      fitResonance.SetParameter(i, fit.GetParameter(i));
   }

   for (int i = 0; i < fitBG.GetNpar(); i++)
   {
      fitBG.SetParameter(i, fit.GetParameter(i + fitResonance.GetNpar()));
   }
   fitResonance.SetRange(fit.GetParameter(1) - fit.GetParameter(2)*10., 
                         fit.GetParameter(1) + fit.GetParameter(2)*10.);
   fitBG.SetRange(fit.GetParameter(1) - fit.GetParameter(2)*10., 
                  fit.GetParameter(1) + fit.GetParameter(2)*10.);
   distrInvM->GetXaxis()->
      SetRange(distrInvM->GetXaxis()->FindBin(fit.GetParameter(1) - fit.GetParameter(2)*10.), 
               distrInvM->GetXaxis()->FindBin(fit.GetParameter(1) + fit.GetParameter(2)*10.));
   distrInvM->GetYaxis()->SetRange(1, distrInvM->GetYaxis()->GetNbins());

   fit.SetLineWidth(4);
   fitResonance.SetLineWidth(3);
   fitBG.SetLineWidth(3);

   fit.SetLineColorAlpha(kRed - 3, 0.8);
   fitResonance.SetLineColorAlpha(kAzure - 3, 0.8);
   fitBG.SetLineColorAlpha(kGreen - 3, 0.8);

   fitResonance.SetLineStyle(2);
   fitBG.SetLineStyle(3);

   distrInvM->SetLineColor(kBlack);
   distrInvM->SetMarkerColor(kBlack);

   TCanvas canv("canv", "", 800, 800);

   gPad->SetRightMargin(0.03);
   gPad->SetTopMargin(0.08);
   gPad->SetLeftMargin(0.14);
   gPad->SetBottomMargin(0.112);

   const double pTMin = distr2DInvM->GetXaxis()->GetBinLowEdge(pTBinMin);
   const double pTMax = distr2DInvM->GetXaxis()->GetBinUpEdge(pTBinMax);

   ROOTTools::DrawFrame(distrInvM, CppTools::DtoStr(pTMin, 1) + " < p_{T} < " + 
                        CppTools::DtoStr(pTMax, 1), "M_{inv} [GeV/c^{2}]", "Counts");

   fitBG.Draw("SAME");
   fitResonance.Draw("SAME");
   fit.Draw("SAME");

   ROOTTools::PrintCanvas(&canv, outputDir + "/" + nameResonance + "_" + 
                          CppTools::DtoStr(pTMin, 1) + "-" + CppTools::DtoStr(pTMax, 1));

   grSigmas.AddPoint(CppTools::Average(pTMin, pTMax), fit.GetParameter(2));
   grSigmas.SetPointError(grSigmas.GetN() - 1, 0., 0.0001);
}

#endif /* ESTIMATE_ESTIMATE_GAUSSIAN_BROADENING_CPP */
