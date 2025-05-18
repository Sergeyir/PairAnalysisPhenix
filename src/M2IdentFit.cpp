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

   inputYAMLSim.OpenFile(argv[1], "m2id");
   inputYAMLSim.CheckStatus("m2id");

   runName = inputYAMLSim["run_name"].as<std::string>();

   inputYAMLMain.OpenFile("input/" + runName + "/main.yaml");
   inputYAMLMain.CheckStatus("main");
 
   collisionSystemName = inputYAMLMain["collision_system_name"].as<std::string>();

   outputDir = "data/PostSim/" + runName + "/SingleTrack/";
   system(("mkdir -p " + outputDir).c_str());
}

double M2IdentFit::GetYield(TH1F *hist, const double mean, const double sigma, 
                            const double vetoLow, const double vetoHigh)
{
}
void M2IdentFit::PerformFitsForDetector(const YAML::Node& detector, 
                                        const double centralityMin, 
                                        const double centralityMax)
{
}
void M2IdentFit::PerformSingleM2Fit(const double pT, TH1F *massProj, FitParameters& fitPar)
{
}

#endif /* M2_IDENT_FIT_CPP */
