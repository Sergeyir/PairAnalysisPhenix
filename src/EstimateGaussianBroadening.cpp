/** 
 *  @file   EstimateGaussianBroadening.cpp 
 *  @brief  Contains realisations of functions that are used for estimation gaussian broadening parameter sigma that used as a parameter in gaus which is convoluted with breit-wigner
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysis).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef GAUSSIAN_BROADENING_CPP
#define GAUSSIAN_BROADENING_CPP

#include <thread>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TLine.h"
#include "TLegend.h"
#include "TMath.h"

#include "StrTools.hpp"
#include "IOTools.hpp"
#include "MathTools.hpp"

#include "TCanvasTools.hpp"

#include "PBar.hpp"

#include "InputYAMLReader.hpp"

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
 
   if (argc == 2) numberOfThreads = std::thread::hardware_concurrency();
   else numberOfThreads = std::stoi(argv[2]);

   ROOT::EnableImplicitMT(numberOfThreads);

   inputYAMLResonance.OpenFile(argv[1]);
   inputYAMLResonance.CheckStatus("resonance");

   runName = inputYAMLResonance["run_name"].as<std::string>();

   inputYAMLMain.OpenFile("input/" + runName + "/main.yaml");
   inputYAMLMain.CheckStatus("main");

   outputDir = "output/GaussianBroadening/" + runName;
   system(("mkdir -p " + outputDir).c_str());

   resonanceName = inputYAMLResonance["name"].as<std::string>();

   const std::string inputFileName = "data/PostSim/" + runName + "/WidthlessResonance/" + 
                                     resonanceName + ".root";

   CheckInputFile(inputFileName);
   inputFile = TFile::Open(inputFileName, "READ");

   ProgressBar pBar("FANCY", "", PBarColor::BOLD_CYAN);

   for (const auto& pTBin : inputYAMLResonance["pt_bins"])
   {
      PerformInvMassFit(pTBin["min"].as<double>(), pTBin["max"].as<double>());
   }
}

EstimateGaussianBroadening::PerformInvMassFit(const double pTMin, const double pTMax)
{
   TH1F *distrInvM = static_cast<TH1F *>(inputFile->Get("M_inv: NoPID"));
}

#endif /* ESTIMATE_GAUSSIAN_BROADENING_CPP */
