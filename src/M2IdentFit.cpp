/** 
 *  @file   M2IdentFit.cpp 
 *  @brief  Contains realistation of functions that are used for approximation of m2 distributions for the estimation of timing parameters that are used for identification via m2
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/CalPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef M2_IDENT_FIT_CPP
#define M2_IDENT_FIT_CPP

#include "M2IdentFit.hpp"

// this namespace is only used so that documentation will not become a mess
// so there is no need to enforce the contents inside of it 
// being accessed only via the scope resolution operator in this file
using namespace M2IdentFit;

int main(int argc, char **argv)
{
   if (argc < 2 || argc > 3) 
   {
      std::string errMsg = "Expected 1-2 parameters while " + std::to_string(argc) + " ";
      errMsg += "parameter(s) were provided \n Usage: bin/M2IdentFit ";
      errMsg += "inputYAMLName numberOfThreads=std::thread::hardware_concurrency()";
      CppTools::PrintError(errMsg);
   }
 
   CppTools::CheckInputFile(argv[1]);
 
   if (argc == 2) numberOfThreads = std::thread::hardware_concurrency();
   else numberOfThreads = std::stoi(argv[2]);

   ROOT::EnableImplicitMT(numberOfThreads);

   inputYAMLM2Id.OpenFile(argv[1], "m2id");
   inputYAMLM2Id.CheckStatus("m2id");

   runName = inputYAMLM2Id["run_name"].as<std::string>();

   inputYAMLMain.OpenFile("input/" + runName + "/main.yaml");
   inputYAMLMain.CheckStatus("main");
 
   collisionSystemName = inputYAMLMain["collision_system_name"].as<std::string>();

   const std::string inputDataFileName = 
      "data/Real/" + collisionSystemName + "/SingleTrack/sum.root";
   CppTools::CheckInputFile(inputDataFileName);
   inputDataFile = TFile::Open(inputDataFileName.c_str(), "READ");

   outputDir = "output/M2Id/" + runName;
   system(("mkdir -p " + outputDir).c_str());

   parametersDir = "data/Parameters/M2Id/" + runName;
   system(("mkdir -p " + parametersDir).c_str());

   rawYieldsDir = "data/RawYields/SingleTrack/" + runName;
   system(("mkdir -p " + rawYieldsDir).c_str());

   gROOT->ProcessLine("gErrorIngonreLevel = 1001;");
   gErrorIgnoreLevel = kWarning;

   gStyle->SetOptStat(0);
   gStyle->SetOptStat(0);

   for (const YAML::Node& detector : inputYAMLM2Id["detectors"])
   {
      PerformFitsForDetector(detector, inputYAMLMain["centrality_min"].as<double>(), 
                             inputYAMLMain["centrality_max"].as<double>());
   }
}

void M2IdentFit::PerformFitsForDetector(const YAML::Node& detector, 
                                        const double centralityMin, 
                                        const double centralityMax)
{
   // container that holds fit parameters and yields for pi^+
   FitParameters fitPiPlus("pi+", pow(0.139570, 2), detector);
   // container that holds fit parameters and yields for K^+
   FitParameters fitKPlus("k+", pow(0.493677, 2), detector);
   // container that holds fit parameters and yields for protons
   FitParameters fitP("p", pow(0.938272, 2), detector);
   // container that holds fit parameters and yields for pi^-
   FitParameters fitPiMinus("pi-", pow(0.139570, 2), detector);
   // container that holds fit parameters and yields for K^-
   FitParameters fitKMinus("k-", pow(0.493677, 2), detector);
   // container that holds fit parameters and yields for antiprotons
   FitParameters fitPBar("pbar", pow(0.938272, 2), detector);

   TH3F* m2DistrPos = inputDataFile->
      Get(("m2, " + detector["name"].as<std::string>() + ", charge>0".c_str());
   TH3F* m2DistrNeg = inputDataFile->
      Get(("m2, " + detector["name"].as<std::string>() + ", charge<0".c_str());

   // minimum pT in whole pT range
   const double pTMin = CppTools::Minimum(detector["pion_pt_min"].as<double>(), 
                                          detector["kaon_pt_min"].as<double>(),
                                          detector["proton_pt_min"].as<double>());
   // maximum pT in whole pT range
   const double pTMax = CppTools::Maximum(detector["pion_pt_max"].as<double>(), 
                                          detector["kaon_pt_max"].as<double>(),
                                          detector["proton_pt_max"].as<double>());

   for (const YAML::Node& pT : inputYAMLM2Id["pt_bins"])
   {
      // minium pT for the current pT bin
      const double binPTMin = pT["min"].as<double>();
      // maximum pT for the current pT bin
      const double binPTMax = pT["max"].as<double>();

      if (binPTMin + 1e-15 < pTMin) continue;
      if (binPTMax - 1e-15 > pTMax) break;

      TH1D *m2DistrPosProj = m2DistrPos->
         ProjectionY((m2DistrPos->GetName() + std::to_string((binPTMin + binPTMax)/2.)).c_str(),
                     m2DistrPos->GetXaxis()->FindBin(binPTMin + 1e-15),
                     m2DistrPos->GetXaxis()->FindBin(binPTMax - 1e-15),
                     m2DistrPos->GetXaxis()->FindBin(centralityMin + 1e-15),
                     m2DistrPos->GetXaxis()->FindBin(centralityMax - 1e-15));

      TH1D *m2DistrNegProj = m2DistrNeg->
         ProjectionY((m2DistrNeg->GetName() + std::to_string((binPTMin + binPTMax)/2.)).c_str(),
                     m2DistrNeg->GetXaxis()->FindBin(binPTMin + 1e-15),
                     m2DistrNeg->GetXaxis()->FindBin(binPTMax - 1e-15),
                     m2DistrNeg->GetXaxis()->FindBin(centralityMin + 1e-15),
                     m2DistrNeg->GetXaxis()->FindBin(centralityMax - 1e-15));

      if (binPTMin > detector["pion_pt_min"].as<double>() && 
          binPTMax < detector["pion_pt_max"].as<double>())
      {
         PerformSingleM2Fit(binPTMin, binPTMax, m2DistrPosProj, fitPiPlus, 
                            detector["pion_bg_func"].as<std::string>(), 
                            distance, sigmaAlpha, sigmaMS, sigmaT);
         PerformSingleM2Fit(binPTMin, binPTMax, m2DistrNegProj, fitPiMinus, 
                            detector["pion_bg_func"].as<std::string>(), 
                            distance, sigmaAlpha, sigmaMS, sigmaT);
      }
      if (binPTMin > detector["kaon_pt_min"].as<double>() && 
          binPTMax < detector["kaon_pt_max"].as<double>())
      {
         PerformSingleM2Fit(binPTMin, binPTMax, m2DistrPosProj, fitKPlus, 
                            detector["kaon_bg_func"].as<std::string>(), 
                            distance, sigmaAlpha, sigmaMS, sigmaT);
         PerformSingleM2Fit(binPTMin, binPTMax, m2DistrPosProj, fitKMinus, 
                            detector["kaon_bg_func"].as<std::string>(), 
                            distance, sigmaAlpha, sigmaMS, sigmaT);
      }
      if (binPTMin > detector["proton_pt_min"].as<double>() && 
          binPTMax < detector["proton_pt_max"].as<double>())
      {
         PerformSingleM2Fit(binPTMin, binPTMax, m2DistrPosProj, fitP,
                            detector["proton_bg_func"].as<std::string>(), 
                            distance, sigmaAlpha, sigmaMS, sigmaT);
         PerformSingleM2Fit(binPTMin, binPTMax, m2DistrPosProj, fitPBar,
                            detector["proton_bg_func"].as<std::string>(), 
                            distance, sigmaAlpha, sigmaMS, sigmaT);
      }
   }

   fitPiPlus.meansVsPTFit->SetRange(detector["pion_pt_min"].as<double>() - 0.05, 
                                    detector["pion_pt_max"].as<double>() + 0.05);
   fitKPlus.meansVsPTFit->SetRange(detector["kaon_pt_min"].as<double>() - 0.05, 
                                   detector["kaon_pt_max"].as<double>() + 0.05);
   fitP.meansVsPTFit->SetRange(detector["proton_pt_min"].as<double>() - 0.05, 
                               detector["proton_pt_max"].as<double>() + 0.05);
   fitPiMinus.meansVsPTFit->SetRange(detector["pion_pt_min"].as<double>() - 0.05, 
                                     detector["pion_pt_max"].as<double>() + 0.05);
   fitKMinus.meansVsPTFit->SetRange(detector["kaon_pt_min"].as<double>() - 0.05, 
                                    detector["kaon_pt_max"].as<double>() + 0.05);
   fitPBar.meansVsPTFit->SetRange(detector["proton_pt_min"].as<double>() - 0.05, 
                                  detector["proton_pt_max"].as<double>() + 0.05);
}

void M2IdentFit::PerformSingleM2Fit(const double pTMin, const double pTMax, TH1F *massDistr,
                                    FitParameters& fitPar, const std::string& funcFG)
{
   const double pT = (pTMin + pTMax)/2.

   fitPar.m2BGFit.push_back
      (std::make_unique<TF1>((name + ": " + std::to_string(pTMin) + "-" + 
                              std::to_string(pTMax)).c_str()), funcBG.c_str());
   fitPar.m2GausFit.push_back
      (std::make_unique<TF1>((name + ": " + std::to_string(pTMin) + "-" + 
                              std::to_string(pTMax)).c_str()), "gaus");
   fitPar.m2Fit.push_back
      (std::make_unique<TF1>((name + ": " + std::to_string(pTMin) + "-" + 
                              std::to_string(pTMax)).c_str()), (funcBG + "gaus(3)").c_str());

   fitPar.m2Fit.back()->
      SetParameter(3, massDistr->GetBinContent(massDistr->GetXaxis()->FindBin(fitPar.m2)));
   fitPar.m2Fit.back()->SetParameter(4, fitPar.m2);
   fitPar.m2Fit.back()->SetParameter(5, fitPar.sigmasVsPT->Eval(pT));

   fitPar.m2Fit.back()->SetParLimits(3, 0., massDistr->GetBinContent(massDistr->GetMaximumBin()));
   fitPar.m2Fit.back()->SetParLimits(4, fitPar.m2 - fitPar.sigmasVsPT->Eval(pT), 
                                     fitPar.m2 + fitPar.sigmasVsPT->Eval(pT));
   fitPar.m2Fit.back()->SetParLimits(5, fitPar.sigmasVsPT->Eval(pT)/1.2, 
                                     fitPar.sigmasVsPT->Eval(pT)*1.2);

   fitPar.m2FitGaus.back()->SetRange(m2 - fitPar.sigmasVsPT->Eval(pT)*3., 
                                     m2 + fitPar.sigmasVsPT->Eval(pT)*3.);

   massDistr->Fit(fitPar.m2FitGaus);
}

double M2IdentFit::GetYield(TH1F *hist, const double mean, const double sigma, 
                            const double sigmalizedYieldExtractionRange, TF1 *fitBG,
                            const double vetoLow, const double vetoHigh, double& err)
{
   double yield = 0; // yield of a signal
   double yieldNoBGSubtr = 0; // yield of a signal + bg i.e. without bg subtraction

   for (int i = hist->GetXaxis()->
           FindBin(CppTools::Maximum(mean - sigmalizedYieldExtractionRange*sigma, veto_low)); 
        i <= hist->GetXaxis()->
           FindBin(CppTools::Minimum(mean + sigmalizedYieldExtractionRange*sigma, veto_up)); i++)
   {
      // this is needed for bin shift correction
      double binAverageFitGaus1 = 0.;
      double binAverageFitGaus2 = 0.;
      double binAverageFitBG = 0.;
      // taking into account the bin shift correction
      for (double m2 = hist->GetXaxis()->GetBinLowEdge(i); 
           m2 < hist->GetXaxis()->GetBinUpEdge(i) + 1e-15; 
           m2 += hist->GetXaxis()->GetBinWidth(i)/99.)
      {
         binAverageFitBG = fitBG->Eval(m2);
      }

      // 100 iterations in the previous loop so the division is by 100
      binAverageFitBG /= 100.;

      yield += hist->GetBinContent(i) - binAverageFitBG;
      yieldNoBGSubtr += hist->GetBinContent(i);
   }

   // correction to the yield (see M2IdentFit::PerformSingleM2Fit definition)
   const double yieldCorrection = 
      (erf((hist->GetXaxis()->GetBinCenter(hist->GetXaxis()->FindBin(mean)) - 
       hist->GetXaxis()->GetBinLowEdge
       (hist->GetXaxis()->FindBin(Maximum(mean - sigmalizedYieldExtractionRange*sigma, veto_low))))*
       sigmalizedYieldExtractionRange/sqrt(2.)/sigma) +
       erf((hist->GetXaxis()->GetBinUpEdge
       (hist->GetXaxis()->FindBin(Minimum(mean + sigmalizedYieldExtractionRange*sigma, veto_up))) - 
        hist->GetXaxis()->GetBinCenter(hist->GetXaxis()->FindBin(mean)))*
       sigmalizedYieldExtractionRange/sqrt(2.)/sigma))/
      2./erf(sigmalizedYieldExtractionRange/sqrt(2.));

   err = sqrt(yield_nosubtr)/yield/yieldCorrection;

   return yield/yieldCorrection;
}

M2IdentFit::FitParameters::FitParameters(const std::string& particleName, const double massSquared,
                                         const YAML::Node& detector)
{
   name = particleName;
   m2 = massSquared;

   meansVsPTFit = std::make_unique<TF1>(("means fit " + name).c_str(), 
                                        inputYAMLM2Id["means_vs_pt_fit_func"].as<char *>());
   sigmasVsPTFit = std::make_unique<TF1>(("sigmas fit " + name).c_str(), 
                                         inputYAMLM2Id["sigmas_vs_pt_fit_func"].as<char *>());

   // expected parameters
   const double distance = detector["L"].as<double>();
   const double sigmaAlpha = detector["sigma_alpha"].as<double>();
   const double sigmaMS = detector["sigma_ms"].as<double>();
   const double sigmaT = detector["sigma_t"].as<double>();

   meansVsPTFit.SetParameter(0, m2);
   meansVsPTFit.SetParameter(1, 0.);
   meansVsPTFit.SetParameter(2, sigmaAlpha);
   meansVsPTFit.SetParameter(3, sigmaMS);
   meansVsPTFit.SetParameter(4, sigmaT);
   meansVsPTFit.FixParameter(5, detector["K1"].as<double>());
   meansVsPTFit.FixParameter(6, detector["L"].as<double>());

   meansVsPTFit.SetParLimits(2, sigmaAlpha/1.2, sigmaAlpha*1.2);
   meansVsPTFit.SetParLimits(3, sigmaMS/1.2, sigmaMS*1.2);
   meansVsPTFit.SetParLimits(4, sigmaT/1.2, sigmaT*1.2);

   rawYieldOutputFile.open(rawYieldsDir + "/" + name + ".txt");
}

#endif /* M2_IDENT_FIT_CPP */
