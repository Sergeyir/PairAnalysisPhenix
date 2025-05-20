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

   outputDir = "data/PostM2Id/" + runName + "/SingleTrack/";
   system(("mkdir -p " + outputDir).c_str());

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
   // containers that hold fit parameters and yields for pi^+, K^+, p, pi^-, K^-, \bar{p} 
   // i.e. charged pions and kaons, protons, and antiprotons
   FitParameters fitPi, fitK, fitP, fitAPi, fitAK, fitPBar;

   fitPi.meansVsPTFit = 
      std::make_unique<TF1>("means fit pi", inputYAMLM2Id["means_vs_pt_fit_func"].as<char *>());
   fitPi.sigmasVsPTFit = 
      std::make_unique<TF1>("sigmas fit pi", inputYAMLM2Id["sigmas_vs_pt_fit_func"].as<char *>());

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

   const double distance = detector["L"].as<double>();
   const double sigmaAlpha = detector["sigma_alpha"].as<double>();
   const double sigmaMS = detector["sigma_ms"].as<double>();
   const double sigmaT = detector["sigma_t"].as<double>();

   for (const YAML::Node& pT : inputYAMLM2Id["pt_bins"])
   {
      // minium pT for the current pT bin
      const double binPTMin = pT["min"].as<double>();
      // maximum pT for the current pT bin
      const double binPTMax = pT["max"].as<double>();

      if (binPTMin + 1e-15 < pTMin || binPTMax - 1e-15 > pTMax) continue;

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
         PerformSingleM2Fit(binPTMin, binPTMax, m2DistrPosProj, fitPi, 
                            detector["pion_bg_func"].as<std::string>(), 
                            distance, sigmaAlpha, sigmaMS, sigmaT);
         PerformSingleM2Fit(binPTMin, binPTMax, m2DistrNegProj, fitAPi, 
                            detector["pion_bg_func"].as<std::string>(), 
                            distance, sigmaAlpha, sigmaMS, sigmaT);
      }
      if (binPTMin > detector["kaon_pt_min"].as<double>() && 
          binPTMax < detector["kaon_pt_max"].as<double>())
      {
         PerformSingleM2Fit(binPTMin, binPTMax, m2DistrPosProj, fitK, 
                            detector["kaon_bg_func"].as<std::string>(), 
                            distance, sigmaAlpha, sigmaMS, sigmaT);
         PerformSingleM2Fit(binPTMin, binPTMax, m2DistrPosProj, fitAK, 
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
}

void M2IdentFit::PerformSingleM2Fit(const double pTMin, const double pTMax, TH1F *massProj,
                                    FitParameters& fitPar, const std::string& funcFG, 
                                    const double distance, const double sigmaAlpha, 
                                    const double sigmaMS, const double sigmaT)
{
}

double M2IdentFit::GetYield(TH1F *hist, const double mean, const double sigma, 
                            const double sigmalizedYieldExtractionRange,
                            TF1 *fitGaus1, TF1 *fitGaus2, TF1 *fitBG,
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
         binAverageFitGaus1 = fitGaus1->Eval(m2);
         binAverageFitGaus2 = fitGaus2->Eval(m2);
         binAverageFitBG = fitBG->Eval(m2);
      }

      // 100 iterations in the previous loop so the division is by 100
      binAverageFitGaus1 /= 100.; 
      binAverageFitGaus2 /= 100.;
      binAverageFitBG /= 100.;

      yield += hist->GetBinContent(i) - binAverageFitGaus1 - binAverageFitGaus2 - binAverageFitBG;
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

#endif /* M2_IDENT_FIT_CPP */
