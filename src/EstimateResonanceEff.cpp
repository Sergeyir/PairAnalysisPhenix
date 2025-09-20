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

   inputYAMLResonance.OpenFile(argv[1]);
   inputYAMLResonance.CheckStatus("resonance");

   runName = inputYAMLResonance["run_name"].as<std::string>();

   inputYAMLMain.OpenFile("input/" + runName + "/main.yaml");
   inputYAMLMain.CheckStatus("main");

   const std::string inputFileNameGausBroadening = "data/Parameters/GaussianBroadening/" + 
                                                   runName + "/" + resonanceName + "txt";
   CppTools::FileExists(inputFile);

   nameResonance = inputYAMLResonance["name"].as<std::string>();
   massResonance = inputYAMLResonance["mass"].as<double>();
   gammaResonance = inputYAMLResonance["gamma"].as<double>();

   SetGaussianBroadeningParameters();

   const std::string inputFileName = "data/PostSim/" + runName + "/Resonance/" + 
                                     nameResonance + ".root";

   CppTools::CheckInputFile(inputFileName);
   inputFile = TFile::Open(inputFileName.c_str(), "READ");

   distrOrigVsPT = static_cast<TH1F *>(inputFile->Get("orig pT"));

   text.SetTextFont(43);
   text.SetTextSize(50);
   texText.SetTextFont(43);
   texText.SetTextSize(50);

   pTNBins = inputYAMLResonance["pt_bins"].size() - 1;

   for (unsigned int i = 0; i < pTNBins; i++)
   {
      pTBinRanges.push_back(inputYAMLResonance["pt_bins"][i]["min"].as<double>());
   }
   pTBinRanges.push_back(inputYAMLResonance["pt_bins"][i]["max"].as<double>());

   // 17 different pairs selections methods
   numberOfIterations = pTNBins*17;

   // performing fits for each pair selection method
   PerformInvMassFits("DCPC1NoPID");
   PerformInvMassFits("NoPID");
   PerformInvMassFits("PC2NoPID");
   PerformInvMassFits("PC3NoPID");
   PerformInvMassFits("TOFeNoPID");
   PerformInvMassFits("TOFwNoPID");
   PerformInvMassFits("EMCalNoPID");
   PerformInvMassFits("DCPC11PID");
   PerformInvMassFits("1TOFDCPC11PID");
   PerformInvMassFits("1EMCalDCPC11PID");
   PerformInvMassFits("1PID");
   PerformInvMassFits("TOF1PID");
   PerformInvMassFits("EMCal1PID");
   PerformInvMassFits("2PID");
   PerformInvMassFits("TOF2PID");
   PerformInvMassFits("EMCal2PID");
   PerformInvMassFits("1TOF1EMCal2PID");

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

void EstimateResonanceEff::PerformInvMassFitsForMethod(const std::string& methodName)
{
   TH2F *distrInvM distrInvM = static_cast<TH2F *>(inputFile->Get(("M_inv:" + methodName).c_str()));

   const std::string outputDir = "output/ResonanceEff/" + runName + "/" + methodName;
   system(("mkdir -p " + outputDir).c_str());

   TH1D distrMeansVsPT("means vs pT", "", pTNBins, &pTBinRanges[0]);
   TH1D distrSigmasVsPT("sigmas vs pT", "", pTNBins, &pTBinRanges[0]);
   TH1D distrRecEffVsPT("reconstruction efficiency vs pT", "", pTNBins, &pTBinRanges[0]);

   for (for unsigned int i = 0; i < pTNBins; i++)
   {
      pBar.Print(static_cast<double>(numberOfCalls)/static_cast<double>(numberOfIterations));

      TH1D *distrInvM = distr2DInvM->
         ProjectionY("proj", distr2DInvM->GetXaixs()->FindBin(pTBinRanges[i] + 1e-6),
                     distr2DInvM->GetXaixs()->FindBin(pTBinRanges[i + 1] - 1e-6));

      // sigma of a gaus that is convoluted with Breit-Wigner
      const double gaussBroadeningSigma = gaussianBroadeningEstimatorFunc.Eval();

      // number of generated resonance particles in the current pT bin
      const double numberOfGenerated = 
         distrOrig->Integral(distrOrig->GetXaxis()->FindBin(pTBinRanges[i] + 1e-6), 
                             distrOrig->GetXaxis()->FindBin(pTBinRanges[i] - 1e-6));

      // for low entries inputs background is ignored since the value 
      // of the statistical uncertainty is significantly larger
      if (distrInvM->GetEntries() < 200) 
      {
         const double reconstructedYield = distrInvMProj->
            Integral(distrInvM->GetXaxis()->FindBin(massResonance - gammaResonance*3. - 
                                                    gaussBroadeningSigma*3.),
                     distrInvM->GetXaxis()->FindBin(massResonance + gammaResonance*3. +
                                                    gaussBroadeningSigma*3.));
         // reconstruction efficiency and its statistical uncertainty (i.e. error)
         const double recEff = reconstructedYield/numberOfGenerated;
         const double recEffErr = CppTools::UncertaintyProp(1./sqrt(reconstructedYield), 
                                                            1./sqrt(numberOfGenerated))*recEff;
         distrRecEffVsPT->SetBinContent(i, recEff);
         distrRecEffVsPT->SetBinError(i, recEffErr);
      }

      // fit for resonance+bg approximation
      TF1 fit("resonance + bg fit", &RBWConvGausBGGaus, 
              massResonance - gammaResonance*3., massResonance + gammaResonance*3., 7);
      // fit for resonance approximation
      TF1 fitResonance("resonance fit", &RBWConvGaus, massResonance - gammaResonance*3., 
                       massResonance + gammaResonance*3., 7);
      // fit for bg approximation
      TF1 fitBG("bg fit", &Gaus, massResonance - gammaResonance*3., 
                massResonance + gammaResonance*3., 7);

      const double maxBinVal = distrInvM->GetBinContent(distrInvM->GetMaximumBin());

      fit.SetParameters(maxBinVal, massResonance, gammaResonance, gaussianBroadeningSigma,
                        maxBinVal/20., massResonance, 0.2);

      fit.SetParLimits(0, maxBinVal/3., maxBinVal);
      fit.SetParLimits(1, massResonance/1.05, massResonance*1.05);
      fit.SetParLimits(2, gammaResonance/1.05, gammaResonance*1.05);
      fit.FixParameter(3, gaussianBroadeningSigma);
      fit.SetParLimits(4, 0., maxBinVal/5.);
      fit.SetParLimits(5, 0., massResonance*10.);
      fit.SetParLimits(6, 5e-2, 1.);

      distrInvM->Fit(&fit, "RQMNB");

      for (unsigned int i = 1; i <= fitNTries; i++)
      {
         fit.SetParLimits(1, fit.GetParameter(1) - 0.02/static_cast<double>(i*i*i), 
                          fit.GetParameter(1) + 0.02/static_cast<double>(i*i*i));
         fit.SetParLimits(2, fit.GetParameter(2)/(1. + 0.2/static_cast<double>(i*i*i)),
                          fit.GetParameter(2)*(1. + 0.2/static_cast<double>(i*i*i)));
         fit.SetParLimits(5, fit.GetParameter(4)/(1. + 2./static_cast<double>(i*i*i)),
                          fit.GetParameter(4)*(1. + 2./static_cast<double>(i*i*i)));
         fit.SetParLimits(6, fit.GetParameter(2)*5.,
                          fit.GetParameter(5)*(1. + 2./static_cast<double>(i*i*i)));

         fit.SetRange(fit.GetParameter(1) - fit.GetParameter(2)*3. - gaussianBroadeningSigma, 
                     fit.GetParameter(1) + fit.GetParameter(2)*3. + gaussianBroadeningSigma);

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

      fitResonance.SetRange(fit.GetParameter(1) - fit.GetParameter(2)*3. - gaussianBroadeningSigma, 
                            fit.GetParameter(1) + fit.GetParameter(2)*3. + gaussianBroadeningSigma);
      fitBG.SetRange(fit.GetParameter(1) - fit.GetParameter(2)*3. - gaussianBroadeningSigma, 
                     fit.GetParameter(1) + fit.GetParameter(2)*3. + gaussianBroadeningSigma);

      distrInvM->GetXaxis()->
         SetRange(CppTools::Minimum(distrInvM->GetXaxis()->FindBin(fit.GetParameter(1) - 
                                                                   fit.GetParameter(2)*5. -
                                                                   gaussianBroadeningSigma*5.), 1), 
                  distrInvM->GetXaxis()->FindBin(fit.GetParameter(1) - fit.GetParameter(2)*5. -
                                                 gaussianBroadeningSigma*5.));

      fit.SetLineWidth(4);
      fitResonance.SetLineWidth(4);
      fitBG.SetLineWidth(4);
      distrInvM->SetLineWidth(2);

      fit.SetLineColorAlpha(kRed - 3, 0.8);
      fitResonance.SetLineColorAlpha(kAzure - 3, 0.8);
      fitBG.SetLineColorAlpha(kGreen - 3, 0.8);

      fitResonance.SetLineStyle(2);
      fitBG.SetLineStyle(7);

      distrInvM->SetLineColor(kBlack);
      distrInvM->SetMarkerColor(kBlack);

      TCanvas canv("canv", "", 800, 800);

      gPad->SetRightMargin(0.03);
      gPad->SetTopMargin(0.02);
      gPad->SetLeftMargin(0.142);
      gPad->SetBottomMargin(0.112);

      ROOTTools::DrawFrame(distrInvM, "", "#it{M}_{inv} [GeV/c^{2}]", "Weighted counts");

      texText.DrawLatexNDC(0.17, 0.9, (CppTools::DtoStr(pTBinRanges[i], 1) + " < #it{p}_{T} < " + 
                           CppTools::DtoStr(pTBinRanges[i + 1], 1)).c_str());

      fitBG.Draw("SAME");
      fitResonance.Draw("SAME");
      fit.Draw("SAME");

      ROOTTools::PrintCanvas(&canv, outputDir + "/" + nameResonance + "_" + 
                             CppTools::DtoStr(pTMin, 1) + "-" + CppTools::DtoStr(pTMax, 1));

      /*
      grSigmas.AddPoint(CppTools::Average(pTMin, pTMax), fit.GetParameter(2));
      grSigmas.SetPointError(grSigmas.GetN() - 1, 0., 0.0001);
      */
      numberOfIterations++;
   }
}

void EstimateResonanceEff::SetGaussianBroadeningParameters()
{
   const std::string inputFileName = "data/Parameters/GaussianBroadening/" + 
                                     runName + "/" + nameResonance + "txt";

   if (!CppTools::FileExists(inputFileName))
   {
      CppTools::PrintError(inputFileName " + does not exists. "\
                           "Run executable bin/EstimateGassianBroadening first");
   }

   std::ifstream inputFile(inputFileName);

   int numberOfParameters;
   std::string function;

   if (!(inputFile >> function >> numberOfParameters)) 
   {
      CppTools::PrintError("Unexpected end of file " + inputFileName);
   }

   gaussianBroadeningEstimationFunc = std::make_unique<TF1>("sigma estimator", function.c_str());

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

double EstimateResonanceEff::GetYield(TH1D *distrInvM, const TF1& funcBG, 
                                      const double xMin, const double xMax, double &err)
{
   double yield;
   for (int i = distrInvM->GetXaxis()->FindBin(xMin); 
        i <= distrInvM->GetXaxis()->FindBin(xMin); i++)
   {
      
   }
}

#endif /* ESTIMATE_ESTIMATE_RESONANCE_EFF_CPP */
