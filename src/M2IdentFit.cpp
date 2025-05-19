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

   TH3F* m2DistrPos = inputDataFile->
      Get(("m2, " + detector["name"].as<std::string>() + ", charge>0".c_str());
   TH3F* m2DistrNeg = inputDataFile->
      Get(("m2, " + detector["name"].as<std::string>() + ", charge<0".c_str());

   for (const YAML::Node& pT : inputYAMLM2Id["pt_bins"])
   {
      // minium pT for the current pT bin
      const double pTMin = pT["min"].as<double>();
      // maximum pT for the current pT bin
      const double pTMax = pT["max"].as<double>();

      // continue here
      TH1F *m2DistrPosProj = static_cast<TH1F *>(m2DistrPos->ProjectionY());

      if (pTMin > detector["pion_pt_min"].as<double>() && 
          pTMax < detector["pion_pt_max"].as<double>())
      {
         PerformSingleM2Fit(pTMin, pTMax, ;
      }
   }
}

void M2IdentFit::PerformSingleM2Fit(const std::string& pTRangeName, 
                                    TH1F *massProj, FitParameters& fitPar)
{
}

double M2IdentFit::GetYield(TH1F *hist, const double mean, const double sigma, 
                            const double vetoLow, const double vetoHigh, double& err)
{
	double yield = 0;
	double yield_nosubtr = 0;
	
	for (int i = hist->GetXaxis()->FindBin(Maximum(mean-Par.yield_extr_srange*sigma, veto_low)); 
			i <= hist->GetXaxis()->FindBin(Minimum(mean+Par.yield_extr_srange*sigma, veto_up)); i++)
	{
		yield += hist->GetBinContent(i) - 
			gaus1->Eval(hist->GetXaxis()->GetBinCenter(i)) -
			gaus2->Eval(hist->GetXaxis()->GetBinCenter(i));
		yield_nosubtr += hist->GetBinContent(i);
	}
	
	const double norm = (
		erf((hist->GetXaxis()->GetBinCenter(hist->GetXaxis()->FindBin(mean)) - 
		hist->GetXaxis()->GetBinLowEdge(
		hist->GetXaxis()->FindBin(Maximum(mean - Par.yield_extr_srange*sigma, veto_low))))*
		Par.yield_extr_srange/sqrt(2.)/sigma) +
		erf((hist->GetXaxis()->GetBinUpEdge(
		hist->GetXaxis()->FindBin(Minimum(mean + Par.yield_extr_srange*sigma, veto_up))) - 
		hist->GetXaxis()->GetBinCenter(hist->GetXaxis()->FindBin(mean)))*
		Par.yield_extr_srange/sqrt(2.)/sigma))/2./erf(Par.yield_extr_srange/sqrt(2.));

	err = sqrt(yield_nosubtr)/yield/norm;
	
	return yield/norm;
}

#endif /* M2_IDENT_FIT_CPP */
