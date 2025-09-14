/** 
 *  @file   EstimateSingleTrackEff.cpp 
 *  @brief  Contains realisations of functions and variables that are used for estimation of registering/identification of charged tracks of pions, kaons, and protons in MC
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/CalPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef ESTIMATE_SINGLE_TRACK_EFF_CPP
#define ESTIMATE_SINGLE_TRACK_EFF_CPP

#include "EstimateSingleTrackEff.hpp"

int main(int argc, char **argv)
{
   if (argc != 2) 
   {
      std::string errMsg = "Expected 1 parameter while " + std::to_string(argc) + " ";
      errMsg += "parameter(s) were provided \n Usage: bin/EstimateSingleTrackEff ";
      errMsg += "inputYAMLName/path to .yaml input";
      CppTools::PrintError(errMsg);
   }
 
   inputYAMLSim.OpenFile(argv[1], "single_track_sim");
   inputYAMLSim.CheckStatus("single_track_sim");
 
   runName = inputYAMLSim["run_name"].as<std::string>();

   inputYAMLMain.OpenFile("input/" + runName "/main.yaml");
   inputYAMLMain.CheckStatus("main");

   pTMin = inputYAMLSim["pt_min"].as<double>();
   pTMax = inputYAMLSim["pt_max"].as<double>();

   gStyle->SetOptStat(0);
   gErrorIgnoreLevel = kWarning;

   outputDir = "output/SingleTrackEff/" + runName;
   system(("mkdir -p " + outputDir).c_str());

   const std::string inputDir = "data/PostSim/" + runName + "/SingleTrack/";

   const std::string inputDataFileNamePiPlus = inputDir + "pi+.root";
   const std::string inputDataFileNameKPlus = inputDir + "k+.root";
   const std::string inputDataFileNameP = inputDir + "p.root";
   const std::string inputDataFileNamePiMinus = inputDir + "pi-.root";
   const std::string inputDataFileNameKMinus = inputDir + "k-.root";
   const std::string inputDataFileNamePBar = inputDir + "pbar.root";

   CppTools::CheckInputFile(inputDataFileNamePiPlus);
   CppTools::CheckInputFile(inputDataFileNameKPlus);
   CppTools::CheckInputFile(inputDataFileNameP);
   CppTools::CheckInputFile(inputDataFileNamePiMinus);
   CppTools::CheckInputFile(inputDataFileNameKMinus);
   CppTools::CheckInputFile(inputDataFileNamePBar);

   inputDataFilePiPlus = TFile::Open(inputDataFileNamePiPlus.c_str());
   inputDataFileKPlus = TFile::Open(inputDataFileNameKPlus.c_str());
   inputDataFileP = TFile::Open(inputDataFileNameP.c_str());
   inputDataFilePiMinus = TFile::Open(inputDataFileNamePiMinus.c_str());
   inputDataFileKMinus = TFile::Open(inputDataFileNameKMinus.c_str());
   inputDataFilePBar = TFile::Open(inputDataFileNamePBar.c_str());

   distrOrigPTPiPlus = static_cast<TH1F *>(inputDataFilePiPlus.Get("orig pT"));
   distrOrigPTKPlus = static_cast<TH1F *>(inputDataFileKPlus.Get("orig pT"));
   distrOrigPTP = static_cast<TH1F *>(inputDataFileP.Get("orig pT"));
   distrOrigPTPiMinus = static_cast<TH1F *>(inputDataFilePiMinus.Get("orig pT"));
   distrOrigPTKMinus = static_cast<TH1F *>(inputDataFileKMinus.Get("orig pT"));
   distrOrigPTPBar = static_cast<TH1F *>(inputDataFilePBar.Get("orig pT"));

   EstimateEffForSingleDistribution("PC2");
   EstimateEffForSingleDistribution("PC3");
   EstimateEffForSingleDistribution("TOFe");
   EstimateEffForSingleDistribution("TOFw");
   EstimateEffForSingleDistribution("EMCale0");
   EstimateEffForSingleDistribution("EMCale1");
   EstimateEffForSingleDistribution("EMCale2");
   EstimateEffForSingleDistribution("EMCale3");
   EstimateEffForSingleDistribution("EMCalw0");
   EstimateEffForSingleDistribution("EMCalw1");
   EstimateEffForSingleDistribution("EMCalw2");
   EstimateEffForSingleDistribution("EMCalw3");
}

EstimateSingleTrackEff::EstimateEffForSingleDetector(const std::string& detectorName, 
                                                     const bools isIdentification)
{

}

#endif /* ESTIMATE_SINGLE_TRACK_EFF_CPP */
