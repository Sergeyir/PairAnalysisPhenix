/** 
 *  @file   EstimateResonanceEff.cpp 
 *  @brief  Contains realisations of functions that are used for estimation of resonance reconstruction efficiency witht he use of the data from MC
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
      std::string errMsg = "Expected 1-2 parameters while " + std::to_string(argc) + " ";
      errMsg += "parameter(s) were provided \n Usage: bin/EstimateResonanceEff ";
      errMsg += "inputYAMLName numberOfThreads=std::thread::hardware_concurrency()";
      CppTools::PrintError(errMsg);
   }
 
   CppTools::CheckInputFile(argv[1]);
 
   if (argc == 2) ROOT::EnableImplicitMT(std::thread::hardware_concurrency());
   else ROOT::EnableImplicitMT(std::stoi(argv[2]));

   gStyle->SetOptStat(0);
   gErrorIgnoreLevel = kWarning;

   TH1F::AddDirectory(kFALSE);
   TH1D::AddDirectory(kFALSE);
   TH2F::AddDirectory(kFALSE);

   inputYAMLResonance.OpenFile(argv[1]);
   inputYAMLResonance.CheckStatus("resonance");

   runName = inputYAMLResonance["run_name"].as<std::string>();

   inputYAMLMain.OpenFile("input/" + runName + "/main.yaml");
   inputYAMLMain.CheckStatus("main");

   nameResonance = inputYAMLResonance["name"].as<std::string>();
   massResonance = inputYAMLResonance["mass"].as<double>();
   gammaResonance = inputYAMLResonance["gamma"].as<double>();

   SetGaussianBroadeningParameters();

   inputFileName = "data/PostSim/" + runName + "/Resonance/" + nameResonance + ".root";

   CppTools::CheckInputFile(inputFileName);
   inputFile = TFile::Open(inputFileName.c_str(), "READ");

   distrOrigUnscaledPT = static_cast<TH1F *>(inputFile->Get("orig unscaled pT"));
   if (!distrOrigUnscaledPT) CppTools::PrintError("Original unscaled pT distribution was not found "\
                                                  " in file " + inputFileName);
   distrOrigPT = static_cast<TH1F *>(inputFile->Get("orig pT"));
   if (!distrOrigPT) CppTools::PrintError("Original pT distribution was not found in file" + 
                                          inputFileName);

   text.SetTextFont(43);
   text.SetTextSize(45);
   text.SetTextAngle(270.);
   texText.SetTextFont(43);
   texText.SetTextSize(45);

   pTNBins = inputYAMLResonance["pt_bins"].size();

   for (unsigned int i = 0; i < pTNBins; i++)
   {
      pTBinRanges.push_back(inputYAMLResonance["pt_bins"][i]["min"].as<double>());
   }
   pTBinRanges.push_back(inputYAMLResonance["pt_bins"][pTNBins - 1]["max"].as<double>());

   // 17 different pairs selections methods
   //numberOfIterations = pTNBins*17;
   numberOfIterations = pTNBins;

   // performing fits for each pair selection method
   PerformMInvFitsForMethod("DCPC1NoPID");
   /*
   PerformMInvFitsForMethod("NoPID");
   PerformMInvFitsForMethod("PC2NoPID");
   PerformMInvFitsForMethod("PC3NoPID");
   PerformMInvFitsForMethod("TOFeNoPID");
   PerformMInvFitsForMethod("TOFwNoPID");
   PerformMInvFitsForMethod("EMCalNoPID");
   PerformMInvFitsForMethod("DCPC11PID");
   PerformMInvFitsForMethod("1TOFDCPC11PID");
   PerformMInvFitsForMethod("1EMCalDCPC11PID");
   PerformMInvFitsForMethod("1PID");
   PerformMInvFitsForMethod("TOF1PID");
   PerformMInvFitsForMethod("EMCal1PID");
   PerformMInvFitsForMethod("2PID");
   PerformMInvFitsForMethod("TOF2PID");
   PerformMInvFitsForMethod("EMCal2PID");
   PerformMInvFitsForMethod("1TOF1EMCal2PID");
   */

   pBar.Finish();

   /*
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
   gPad->SetTopMargin(0.02);
   gPad->SetLeftMargin(0.17);
   gPad->SetBottomMargin(0.112);

   ROOTTools::DrawFrame(pTMin, 0., pTMax, TMath::MaxElement(grSigmas.GetN(), grSigmas.GetY())*1.2,
                        "", "p_{T} [GeV/c]", "#sigma [GeV/c^{2}]", 1., 1.8);

   fitSigmas.Draw("SAME");
   grSigmas.Clone()->Draw("P");

   ROOTTools::PrintCanvas(&canv, outputDir + "/" + nameResonance + "_sigmas");

   system(("mkdir -p data/Parameters/ResonanceEff/" + runName).c_str());

   std::ofstream parametersOutput("data/Parameters/ResonanceEff/" + runName + 
                                  "/" + nameResonance + ".txt");

   parametersOutput << fitSigmasFormula << std::endl;
   parametersOutput << fitSigmas.GetNpar() << std::endl;
   for (int i = 0; i < fitSigmas.GetNpar() - 1; i++)
   {
      parametersOutput << fitSigmas.GetParameter(i) << " ";
   }
   parametersOutput << fitSigmas.GetParameter(fitSigmas.GetNpar() - 1);

   */
   CppTools::PrintInfo("EstimateResonanceEff executable has finished running succesfully");
}

void EstimateResonanceEff::PerformMInvFitsForMethod(const std::string& methodName)
{
   const std::string distr2DMInvName = "M_inv: " + methodName;
   TH2F *distr2DMInv = static_cast<TH2F *>(inputFile->Get(distr2DMInvName.c_str()));
   if (!distr2DMInv) CppTools::PrintError("Distribution named " + distr2DMInvName + "\" "\
                                          "was not found in file " + inputFileName);

   const std::string outputDir = "output/ResonanceEff/" + runName + "/" + methodName;
   system(("mkdir -p " + outputDir).c_str());

   TH1D distrMeansVsPT("means vs pT", "", pTNBins, &pTBinRanges[0]);
   TH1D distrGammasVsPT("sigmas vs pT", "", pTNBins, &pTBinRanges[0]);
   TH1D distrRecEffVsPT("reconstruction efficiency vs pT", "", pTNBins, &pTBinRanges[0]);

   for (unsigned int i = 0; i < pTNBins; i++)
   {
      pBar.Print(static_cast<double>(numberOfCalls)/static_cast<double>(numberOfIterations));

      TH1D *distrMInv = distr2DMInv->
         ProjectionY("proj", distr2DMInv->GetXaxis()->FindBin(pTBinRanges[i] + 1e-6),
                     distr2DMInv->GetXaxis()->FindBin(pTBinRanges[i + 1] - 1e-6));

      // sigma of a gaus that is convoluted with Breit-Wigner
      const double gaussianBroadeningSigma = 
         gaussianBroadeningEstimatorFunc->Eval((pTBinRanges[i] + pTBinRanges[i + 1])/2.);

      // number of generated resonance particles in the current pT bin
      // statistical uncertainty of this value is insignificant and therefore can be ignored
      const double numberOfGenerated = 
         distrOrigPT->Integral(distrOrigPT->GetXaxis()->FindBin(pTBinRanges[i] + 1e-6), 
                               distrOrigPT->GetXaxis()->FindBin(pTBinRanges[i + 1] - 1e-6));
      // relative uncertanty of the number of generated particles in the current pT bin
      const double numberOfGeneratedRelativeErr = 
         1./sqrt(distrOrigUnscaledPT->
                 Integral(distrOrigPT->GetXaxis()->FindBin(pTBinRanges[i] + 1e-6), 
                          distrOrigPT->GetXaxis()->FindBin(pTBinRanges[i + 1] - 1e-6)));

      // for low entries inputs background is ignored since the value 
      // of the statistical uncertainty is significantly larger
      if (distrMInv->GetEntries() < 200.) 
      {
         const double reconstructedYield = distrMInv->
            Integral(distrMInv->GetXaxis()->FindBin(massResonance - 
                                                    (gammaResonance + gaussianBroadeningSigma)*3.),
                     distrMInv->GetXaxis()->FindBin(massResonance + 
                                                    (gammaResonance + gaussianBroadeningSigma)*3.));
         // reconstruction efficiency
         const double recEff = reconstructedYield/numberOfGenerated;
         // reconstructed yield relative statistical uncertainty
         // due to the spectra scaling statistical uncertainty is not tied to the integral 
         // but rather to the number of entries of the signal
         double recYieldRelativeErr = 
            1./sqrt(distrMInv->GetEntries()*reconstructedYield/
                    distrMInv->Integral(1, distrMInv->GetXaxis()->GetNbins()));
         const double recEffErr = CppTools::UncertaintyProp(numberOfGeneratedRelativeErr, 
                                                            recYieldRelativeErr)*recEff;

         distrRecEffVsPT.SetBinContent(i, recEff);
         distrRecEffVsPT.SetBinError(i, recEffErr);
      }

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

      fit.SetParLimits(0, maxBinVal/3., maxBinVal);
      fit.SetParLimits(1, massResonance/1.05, massResonance*1.05);
      fit.SetParLimits(2, gammaResonance/1.05, gammaResonance*1.05);
      fit.FixParameter(3, gaussianBroadeningSigma);
      fit.SetParLimits(4, 0., maxBinVal/3.);
      fit.SetParLimits(5, 0., massResonance*10.);
      fit.SetParLimits(6, gammaResonance*2., gammaResonance*20.);

      distrMInv->Fit(&fit, "RQMNB");

      for (unsigned int j = 1; j <= fitNTries; j++)
      {
         fit.SetParLimits(1, fit.GetParameter(1) - 0.02/static_cast<double>(j*j*j), 
                          fit.GetParameter(1) + 0.02/static_cast<double>(j*j*j));
         fit.SetParLimits(2, fit.GetParameter(2)/(1. + 0.2/static_cast<double>(j*j*j)),
                          fit.GetParameter(2)*(1. + 0.2/static_cast<double>(j*j*j)));
         fit.SetParLimits(5, fit.GetParameter(5)/(1. + 2./static_cast<double>(j*j*j)),
                          fit.GetParameter(5)*(1. + 2./static_cast<double>(j*j*j)));

         fit.SetRange(fit.GetParameter(1) - (fit.GetParameter(2) + gaussianBroadeningSigma)*3., 
                     fit.GetParameter(1) + (fit.GetParameter(2) + gaussianBroadeningSigma)*3.);

         distrMInv->Fit(&fit, "RQMNB");
      }

      distrMeansVsPT.SetBinContent(i + 1, fit.GetParameter(1));
      distrMeansVsPT.SetBinError(i + 1, fit.GetParError(1));

      distrGammasVsPT.SetBinContent(i + 1, fit.GetParameter(2));
      distrGammasVsPT.SetBinError(i + 1, fit.GetParError(2));

      double recYieldErr;
      const double recYield = 
         GetYield(distrMInv, fitBG, fit.GetParameter(1) - fit.GetParameter(2)*2., 
                  fit.GetParameter(1) + fit.GetParameter(2)*2., recYieldErr);

      distrRecEffVsPT.SetBinContent(i + 1, recYield/numberOfGenerated);
      distrRecEffVsPT.SetBinError(i + 1, CppTools::UncertaintyProp(recYieldErr/recYield, 
                                                               numberOfGeneratedRelativeErr)*
                                  recYield/numberOfGenerated);

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

      text.DrawTextNDC(0.9, 0.95, methodName.c_str());
      texText.DrawLatexNDC(0.2, 0.9, (CppTools::DtoStr(pTBinRanges[i], 1) + " < #it{p}_{T} < " + 
                           CppTools::DtoStr(pTBinRanges[i + 1], 1)).c_str());
      texText.DrawLatexNDC(0.2, 0.83, 
                           ("#it{#chi}^{2}/NDF = " + 
                            CppTools::DtoStr(fit.GetChisquare()/fit.GetNDF(), 2)).c_str());

      /*
      texText.DrawLatexNDC(0.6, 0.9, ("#it{#mu}=" + CppTools::DtoStr(fit.GetParameter(1)*1000., 1) + 
                                      " [MeV/c^{2}]").c_str());
      texText.DrawLatexNDC(0.6, 0.85, ("#it{#Gamma}=" + 
                                       CppTools::DtoStr(fit.GetParameter(2)*1000., 1) + 
                                       " [MeV/c^{2}]").c_str());
                                       */

      fitBG.Draw("SAME");
      fitResonance.Draw("SAME");
      fit.Draw("SAME");

      ROOTTools::PrintCanvas(&canvMInv, outputDir + "/" + nameResonance + "_" + 
                             CppTools::DtoStr(pTBinRanges[i], 1) + "-" + 
                             CppTools::DtoStr(pTBinRanges[i + 1], 1));

      numberOfCalls++;
   }

   distrMeansVsPT.SetLineColor(kBlack);
   distrMeansVsPT.SetMarkerColor(kBlack);
   distrMeansVsPT.SetLineWidth(3);

   distrGammasVsPT.SetLineColor(kBlack);
   distrGammasVsPT.SetMarkerColor(kBlack);
   distrGammasVsPT.SetLineWidth(3);

   distrRecEffVsPT.SetLineColor(kBlack);
   distrRecEffVsPT.SetMarkerColor(kBlack);
   distrRecEffVsPT.SetLineWidth(3);

   TLine massResonancePDG(pTBinRanges[0], massResonance, pTBinRanges[pTNBins], massResonance);
   massResonancePDG.SetLineColorAlpha(kBlack, 0.5);
   massResonancePDG.SetLineStyle(2);
   massResonancePDG.SetLineWidth(3);

   TLine gammaResonancePDG(pTBinRanges[0], gammaResonance, pTBinRanges[pTNBins], gammaResonance);
   gammaResonancePDG.SetLineColorAlpha(kBlack, 0.5);
   gammaResonancePDG.SetLineStyle(2);
   gammaResonancePDG.SetLineWidth(3);

   distrMeansVsPT.SetMaximum(CppTools::Maximum(distrMeansVsPT.GetMaximum(),
                                               massResonance)*1.01);
   distrMeansVsPT.SetMinimum(CppTools::Minimum(distrMeansVsPT.GetMinimum(),
                                               massResonance)/1.01);

   distrGammasVsPT.SetMaximum(CppTools::Maximum(distrGammasVsPT.GetMaximum(),
                                                gammaResonance)*1.01);
   distrGammasVsPT.SetMinimum(CppTools::Minimum(distrGammasVsPT.GetMinimum(),
                                                gammaResonance)/1.01);

   TCanvas canvMeansVsPT("canv means vs pT", "", 800, 800);

   gPad->SetRightMargin(0.03);
   gPad->SetTopMargin(0.02);
   gPad->SetLeftMargin(0.17);
   gPad->SetBottomMargin(0.112);

   ROOTTools::DrawFrame(&distrMeansVsPT, "", "#it{p}_{T} [GeV/c]", 
                        "#it{#mu} [GeV/c^{2}]", 1., 1.8);

   massResonancePDG.Draw();

   ROOTTools::PrintCanvas(&canvMeansVsPT, outputDir + "/" + nameResonance + "_means");

   TCanvas canvGammasVsPT("canv means vs pT", "", 800, 800);

   gPad->SetRightMargin(0.03);
   gPad->SetTopMargin(0.02);
   gPad->SetLeftMargin(0.17);
   gPad->SetBottomMargin(0.112);

   ROOTTools::DrawFrame(&distrGammasVsPT, "", "#it{p}_{T} [GeV/c]", 
                        "#it{#sigma} [GeV/c^{2}]", 1., 1.8);

   gammaResonancePDG.Draw();

   ROOTTools::PrintCanvas(&canvGammasVsPT, outputDir + "/" + nameResonance + "_sigmas");

   TCanvas canvRecEffVsPT("canv means vs pT", "", 800, 800);

   gPad->SetLogy();

   gPad->SetRightMargin(0.03);
   gPad->SetTopMargin(0.02);
   gPad->SetLeftMargin(0.142);
   gPad->SetBottomMargin(0.112);

   ROOTTools::DrawFrame(&distrRecEffVsPT, "", "#it{p}_{T} [GeV/c]", "#it{#varepsilon}_{rec}");

   ROOTTools::PrintCanvas(&canvRecEffVsPT, outputDir + "/" + nameResonance + "_rec_eff");
}

void EstimateResonanceEff::SetGaussianBroadeningParameters()
{
   const std::string inputFileName = "data/Parameters/GaussianBroadening/" + 
                                     runName + "/" + nameResonance + ".txt";

   if (!CppTools::FileExists(inputFileName))
   {
      CppTools::PrintError(inputFileName + " does not exists. "\
                           "Run executable bin/EstimateGassianBroadening first");
   }

   std::ifstream inputFile(inputFileName);

   int numberOfParameters;
   std::string function;

   if (!(inputFile >> function >> numberOfParameters)) 
   {
      CppTools::PrintError("Unexpected end of file " + inputFileName);
   }

   gaussianBroadeningEstimatorFunc = std::make_unique<TF1>("sigma estimator", function.c_str());

   for (int i = 0; i < numberOfParameters; i++)
   {
      double par;
      if (!(inputFile >> par))  
      {
         CppTools::PrintError("Unexpected end of file " + inputFileName);
      }
      gaussianBroadeningEstimatorFunc->SetParameter(i, par);
   }
}

double EstimateResonanceEff::GetYield(TH1D *distrMInv, const TF1& funcBG, 
                                      const double xMin, const double xMax, double &err)
{
   // integral over the signal
   double integral = 0.;
   // integral over the background
   double integralBG = 0.;

   for (int i = distrMInv->GetXaxis()->FindBin(xMin); 
        i <= distrMInv->GetXaxis()->FindBin(xMin); i++)
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
   err = 1./sqrt(distrMInv->GetEntries()*integral/
                 distrMInv->Integral(1, distrMInv->GetXaxis()->GetNbins()));
   // normalizing background integral by the number of integration steps
   integral -= integralBG/100.;

   return integral;
}

#endif /* ESTIMATE_ESTIMATE_RESONANCE_EFF_CPP */
