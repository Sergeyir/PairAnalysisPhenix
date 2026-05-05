/** 
 *  @file   MInvFit.cpp 
 *  @brief  Contains implementation for usage of ROOTTools::GUIFit for tweaking bad background approximations on PHENIX processed data
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysis).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#pragma once

#include "ErrorHandler.hpp"
#include "IOTools.hpp"
#include "StrTools.hpp"

#include "GUIFit.hpp"

#include "AnalyzeRealMInv.hpp"

/// Contains invariant mass distributions in different pT ranges in one centrality class.
std::vector<TH1D *> distrsMInv;
/// Contains pT value of invariant mass distributions
std::vector<double> pTsDistrMInv;
/// Contains names of invariant mass distributions to be displayed on canvas
std::vector<std::string> distrMInvNames;
/// Contains fits of type signal+background obtained from invariant mass distributions 
/// in different pT ranges in one centrality class. These fits are default and are 
/// the ones from which results are obtained. Therefore fits must represent an 
/// approximation that is good while chi^2/NDF value is as small as possible. 
/// Mean and Gamma are limited within 10% of PDG value but not fixed.
std::vector<TF1 *> fits;
/// Contains background fits obtained from invariant mass distributions in different pT 
/// ranges in one centrality class
std::vector<TF1 *> fitsBG;
/// Contains alternative fits of type signal+background obtained from invariant mass distributions 
/// in different pT ranges in one centrality class. These fits are same as the ones in "fits" vector
/// but have a different background approximation.
std::vector<TF1 *> altFitsAB;
/// Contains background fits obtained from invariant mass distributions in different pT 
/// ranges in one centrality class
std::vector<TF1 *> altFitsBGAB;
/// Contains background fits obtained from invariant mass distributions in different pT 
/// ranges in one centrality class
/// Contains alternative fits with Mean and Gamma as free parameters of type signal+background 
/// obtained from invariant mass distributions in different pT ranges in one centrality class
std::vector<TF1 *> altFitsFreeG;
/// Contains background fits for altFitsFreeG obtained from invariant 
/// mass distributions in different pT ranges in one centrality class
std::vector<TF1 *> altFitsBGFreeG;
/// Contains alternative fits with Mean and Gamma as fixed parameters of type signal+background 
/// obtained from invariant mass distributions in different pT ranges in one centrality class.
/// Mean and Gamma for these fits are equal to PDG values and are fixed.
std::vector<TF1 *> altFitsFixedG;
/// Contains background fits for altFitsFixedG obtained from invariant 
/// mass distributions in different pT ranges in one centrality class
std::vector<TF1 *> altFitsBGFixedG;

// shows whether the alternative fits will be performed and written
// this value needs to be true if you want to later use different fits in this run
// for systematic uncertainty evaluation
bool doAltFits = true;

using namespace AnalyzeRealMInv;

void MInvFit()
{
	gStyle->SetOptStat(0);
   gErrorIgnoreLevel = kWarning;

   TH1::AddDirectory(kFALSE);
   TH2::AddDirectory(kFALSE);

   ROOT::EnableImplicitMT();
   gROOT->SetBatch(true);

   CppTools::PrintInfo("List of directories in data/Real directory that contain "\
                       "non empty \"Resonance\" directory");
   system("find data/Real/ -name Resonance -type d -not -empty");

   CppTools::Print("Choose the run (part of a directory, example: Run14HeAu200) "\
                   "from the above directory list and type it in");
   std::cout << ">> ";
   std::cin >> runName;
   //runName = "Run14HeAu200"; // temporary for testing

   CppTools::PrintInfo("List of .root files in " + runName + " taxi directory");
   system(("find data/Real/" + runName + "/Resonance/ -type f -name *.root").c_str());

   CppTools::Print("Choose taxi number (part of .root file, example: 20292) "\
                   "from the above directory list and type it in");
   std::cout << ">> ";
   std::cin >> taxiNumber;

   //taxiNumber = 20025; // temporary for testing

   inputFileName = "data/Real/" + runName + "/Resonance/" + std::to_string(taxiNumber) + ".root";
   CppTools::CheckInputFile(inputFileName);

   /* Temporarily disable since only K*(892) is needed at the moment
   CppTools::Print("Choose the particle");
   std::string particleName;
   std::cout << ">> ";
   std::cin >> particleName;
   */
   const std::string particleName = "KStar892";

   inputYAMLResonance.OpenFile("input/" + runName + "/" + particleName + ".yaml");
   inputYAMLResonance.CheckStatus("resonance");

   InputYAMLReader inputYAMLMain("input/" + runName + "/main.yaml");
   inputYAMLMain.CheckStatus("main");

   int methodBinIndex;
   while (true) // ininite loop until exit or valid input is specified
   {
      CppTools::PrintInfo("List of pair selection method bins for " + 
                          particleName + " in " + runName);

      for (int i = 0; i < inputYAMLResonance["pair_selection_methods"].size(); i++)
      {
         CppTools::Print(i, inputYAMLResonance["pair_selection_methods"][i]
                                              ["name"].as<std::string>());
      }

      CppTools::Print("Choose the pair selection method bin index from the list above "\
                      "(typing in any text will exit the program)");
      std::cout << ">> ";
      if (!(std::cin >> methodBinIndex))
      {
         CppTools::PrintInfo("Exiting the program");
         exit(1);
      }

      if (methodBinIndex < inputYAMLResonance["pair_selection_methods"].size() && 
          methodBinIndex >= 0) break;
      else CppTools::PrintWarning("Chosen method bin is out of range");
   }

   int centralityBinIndex;
   while (true) // ininite loop until exit or valid input is specified
   {
      CppTools::PrintInfo("List of centrality bins");

      for (int i = 0; i < inputYAMLResonance["centrality_bins"].size(); i++)
      {
         CppTools::Print(i, inputYAMLResonance["centrality_bins"][i]["name"].as<std::string>());
      }

      CppTools::Print("Choose the centrality bin index from the list above "\
                      "(typing in any text will exit the program)");
      std::cout << ">> ";
      if (!(std::cin >> centralityBinIndex))
      {
         CppTools::PrintInfo("Exiting the program");
         exit(1);
      }

      if (centralityBinIndex < inputYAMLResonance["centrality_bins"].size() &&
          centralityBinIndex >= 0) break;
      else CppTools::PrintWarning("Chosen centrality bin is out of range");
   }

   /*
   CppTools::Print("Choose whether you want to perform and save alternative fits "\
                   "(any number if yes, 0 if no);\n (typing in any text will exit the program)");
   std::cout << ">> ";
   if (!(std::cin >> doAltFits))
   {
      CppTools::PrintInfo("Exiting the program");
      exit(1);
   }
   */

   resonanceName = inputYAMLResonance["name"].as<std::string>();
   massResonance = inputYAMLResonance["mass"].as<double>();
   gammaResonance = inputYAMLResonance["gamma"].as<double>();

   daughter1Id = inputYAMLResonance["daughter1_id"].as<int>();
   daughter2Id = inputYAMLResonance["daughter2_id"].as<int>();

   minMInv = inputYAMLResonance["m_inv_range_min"].as<double>();
   maxMInv = inputYAMLResonance["m_inv_range_max"].as<double>();

   SetGaussianBroadeningFunction();

   inputFile = TFile::Open(inputFileName.c_str(), "READ");

   pTNBins = inputYAMLResonance["pt_bins"].size();

   for (unsigned int i = 0; i < pTNBins; i++)
   {
      pTBinRanges.push_back(inputYAMLResonance["pt_bins"][i]["min"].as<double>());
   }
   pTBinRanges.push_back(inputYAMLResonance["pt_bins"][pTNBins - 1]["max"].as<double>());

   const YAML::Node method = inputYAMLResonance["pair_selection_methods"][methodBinIndex];

   const unsigned int pTBinFitMin = method["centrality_bin_parameters"][centralityBinIndex]
                                          ["pt_bin_min"].as<int>();
   const unsigned int pTBinFitMax = method["centrality_bin_parameters"][centralityBinIndex]
                                          ["pt_bin_max"].as<int>();

   numberOfIterations = pTBinFitMax - pTBinFitMin + 1;

   const std::string parametersOutputDir = "data/Parameters/BGFitResonance/" + runName;
   std::filesystem::create_directories(parametersOutputDir);

   PerformMInvFits(method, centralityBinIndex);

   gROOT->SetBatch(false);

	TCanvas *canv = new TCanvas("", "", 1080, 1080);

   const std::string methodName = 
      inputYAMLResonance["pair_selection_methods"][methodBinIndex]["name"].as<std::string>();
   const std::string centralityName = 
      inputYAMLResonance["centrality_bins"][centralityBinIndex]["name"].as<std::string>();
   
   const std::string fitOutputDir = "data/Parameters/BGFitResonance/" + runName + "/" + 
                                    std::to_string(taxiNumber) + "/";

   std::filesystem::create_directories(fitOutputDir);

   // comments show indices of fit types
   // 0
   GUIFit::AddFitType(fitOutputDir + methodName + "_" + 
                      centralityName + ".root", "Default");
   if (doAltFits)
   {
      // 1
      if (altFitsAB.size() != 0)
      {
         GUIFit::AddFitType(fitOutputDir + methodName + "_" + 
                            centralityName + "_AB.root", "alt BG");
      }
      // 2
      if (altFitsFreeG.size() != 0)
      {
         GUIFit::AddFitType(fitOutputDir + methodName + "_" + 
                            centralityName + "_FreeG.root", "Free #Gamma");
         }
      if (altFitsFixedG.size() != 0)
      {
         // 3
         GUIFit::AddFitType(fitOutputDir + methodName + "_" + 
                            centralityName + "_FixedG.root", "Fixed #Gamma");
      }
   }

   for (unsigned int i = 0; i < fits.size(); i++)
   {
      GUIFit::AddHistogram(distrsMInv[i], CppTools::DtoStr(pTsDistrMInv[i], 2), distrMInvNames[i]);

      GUIFit::AddFit(fits[i], fitsBG[i], 0, fits[i]->GetNpar() - fitsBG[i]->GetNpar());

      if (doAltFits)
      {
         GUIFit::AddFit(altFitsAB[i], altFitsBGAB[i], 1, 
                        altFitsAB[i]->GetNpar() - altFitsBGAB[i]->GetNpar());
         GUIFit::AddFit(altFitsFreeG[i], altFitsBGFreeG[i], 2, 
                        altFitsFreeG[i]->GetNpar() - altFitsBGFreeG[i]->GetNpar());
         GUIFit::AddFit(altFitsFixedG[i], altFitsBGFixedG[i], 3, 
                        altFitsFixedG[i]->GetNpar() - altFitsBGFixedG[i]->GetNpar());
      }
   }

   gPad->SetTopMargin(0.05);
   gPad->SetRightMargin(0.05);
   gPad->SetBottomMargin(0.05);
   gPad->AddExec("exec", "GUIFit::Exec()");
}

void AnalyzeRealMInv::PerformMInvFits(const YAML::Node& method, const unsigned int centralityBinIndex)
{
   const std::string methodName = method["name"].as<std::string>();

   const YAML::Node centrality = inputYAMLResonance["centrality_bins"][centralityBinIndex];

   const std::string centralityName = centrality["name"].as<std::string>();

   const double sigmalizedFitRange = method["sigmalized_fit_range"].as<double>();

   pBar.SetText("Preparing M_{inv}");

   // input files with fixed BG fits
   TFile *inputFileFitsBG = 
      SetFixedBGFile("data/Parameters/BGFitResonance/" + runName + "/" + 
                     std::to_string(taxiNumber) + "/" + methodName + "_" + 
                     centralityName + ".root", "default"); // for default fit
   TFile *inputFileFitsBGAB = nullptr;
   TFile *inputFileFitsBGFreeG = nullptr;
   TFile *inputFileFitsBGFixedG = nullptr;

   if (doAltFits)
   {
      inputFileFitsBGAB = 
         SetFixedBGFile("data/Parameters/BGFitResonance/" + runName + "/" + 
                        std::to_string(taxiNumber) + "/" + methodName + "_" + 
                        centralityName + "_AB.root", "alternative BG"); 
      inputFileFitsBGFreeG = 
         SetFixedBGFile("data/Parameters/BGFitResonance/" + runName + "/" + 
                        std::to_string(taxiNumber) + "/" + methodName + "_" + 
                        centralityName + "_FreeG.root", "free G"); 
      inputFileFitsBGFixedG = 
         SetFixedBGFile("data/Parameters/BGFitResonance/" + runName + "/" + 
                        std::to_string(taxiNumber) + "/" + methodName + "_" + 
                        centralityName + "_FixedG.root", "fixed BG"); 
   }

   const unsigned int pTBinFitMin = method["centrality_bin_parameters"][centralityBinIndex]
                                          ["pt_bin_min"].as<int>();
   const unsigned int pTBinFitMax = method["centrality_bin_parameters"][centralityBinIndex]
                                          ["pt_bin_max"].as<int>();

   for (unsigned int i = pTBinFitMin; i <= pTBinFitMax; i++)
   {
      pBar.Print(static_cast<double>(numberOfCalls)/static_cast<double>(numberOfIterations));

      TH1D *distrMInvFG = nullptr;
      TH1D *distrMInvBG = nullptr;
      TH1D *distrMInvFGLR = nullptr; // low resolution (for BG normalization)
      TH1D *distrMInvBGLR = nullptr; // low resolution (for BG normalization)

      std::string decayMode = ParticleMap::nameShort[daughter1Id] +
                              ParticleMap::nameShort[daughter2Id];

      double numberOfEvents = 0.;

      TH1D *distrMInv = 
         MInv::Merge(inputFile, methodName, decayMode, 
                     centrality["cb_c_min"].as<int>(), centrality["cb_c_max"].as<int>(),
                     0, inputYAMLResonance["cb_z_bins"].as<int>() - 1, 
                     0, inputYAMLResonance["cb_r_bins"].as<int>() - 1,
                     pTBinRanges[i], pTBinRanges[i + 1],
                     distrMInvFG, distrMInvBG, distrMInvFGLR, distrMInvBGLR, 
                     numberOfEvents);

      if (inputYAMLResonance["has_antiparticle"].as<bool>() && 
          !inputYAMLResonance["separate_antiparticle"].as<bool>())
      {
         decayMode = ParticleMap::nameShort[daughter2Id] +
                     ParticleMap::nameShort[daughter1Id];
         distrMInv->Add(MInv::Merge(inputFile, methodName, decayMode, 
                                    centrality["cb_c_min"].as<int>(), 
                                    centrality["cb_c_max"].as<int>(),
                                    0, inputYAMLResonance["cb_z_bins"].as<int>() - 1, 
                                    0, inputYAMLResonance["cb_r_bins"].as<int>() - 1,
                                    pTBinRanges[i], pTBinRanges[i + 1],
                                    distrMInvFG, distrMInvBG, distrMInvFGLR, distrMInvBGLR,
                                    numberOfEvents));
      }

      if (!distrMInv)
      {
         pBar.Clear();
         CppTools::PrintError("Resulting M_{inv} histogram could not be constructed for " + 
                              CppTools::DtoStr(pTBinRanges[i], 2) + "<pT<" + 
                              CppTools::DtoStr(pTBinRanges[i + 1], 2));
      }
      else if (distrMInv->Integral(1, distrMInv->GetXaxis()->GetNbins()) < 1e-7)
      {
         pBar.Clear();
         CppTools::PrintWarning("Resulting histogram is empty in " + 
                                centrality["name"].as<std::string>() + " " +
                                CppTools::DtoStr(pTBinRanges[i], 2) + "<pT<" + 
                                CppTools::DtoStr(pTBinRanges[i + 1], 2));
         pBar.RePrint();
         continue;
      }

      if (!distrMInvFG)
      {
         pBar.Clear();
         CppTools::PrintError("Resulting M_{inv} foreground histogram "\
                              "could not be constructed for " + 
                              centrality["name"].as<std::string>() + " " +
                              CppTools::DtoStr(pTBinRanges[i], 2) + "<pT<" + 
                              CppTools::DtoStr(pTBinRanges[i + 1], 2));
      }

      if (!distrMInvBG)
      {
         pBar.Clear();
         CppTools::PrintWarning("Resulting M_{inv} foreground histogram "\
                                "could not be constructed for " + 
                                centrality["name"].as<std::string>() + " " +
                                CppTools::DtoStr(pTBinRanges[i], 2) + "<pT<" + 
                                CppTools::DtoStr(pTBinRanges[i + 1], 2));
         pBar.RePrint();
      }
      else if (distrMInvBG->Integral(1, distrMInvBG->GetXaxis()->GetNbins()) < 1e-7)
      {
         pBar.Clear();
         CppTools::PrintWarning("Resulting background histogram is empty in " + 
                                centrality["name"].as<std::string>() + 
                                CppTools::DtoStr(pTBinRanges[i], 2) + "<pT<" + 
                                CppTools::DtoStr(pTBinRanges[i + 1], 2));
         pBar.RePrint();
      }

      if (rebinX != 1)
      {
         distrMInv->Rebin(rebinX);
         distrMInvFG->Rebin(rebinX);
         distrMInvBG->Rebin(rebinX);
      }

      // sigma of a gaus that is convoluted with Breit-Wigner
      const double gaussianBroadeningSigma = 
         gaussianBroadeningEstimatorFunc->Eval((pTBinRanges[i] + pTBinRanges[i + 1])/2.);

      const std::string bgFitFunc = method["bg_fit_func"].as<std::string>();
      if (bgFitFunc == "pol2")
      {
         fits.push_back(new TF1(("Default" + std::to_string(i)).c_str(), 
                                &FitFunc::RBWConvGausBGPol2, 
                                massResonance - gammaResonance*3., 
                                massResonance + gammaResonance*3., 7));
         fitsBG.push_back(new TF1(("Default BG" + std::to_string(i)).c_str(), 
                                  &FitFunc::Pol2, 
                                  massResonance - gammaResonance*3., 
                                  massResonance + gammaResonance*3., 3));

         if (doAltFits)
         {
            altFitsAB.push_back(new TF1(("AB" + std::to_string(i)).c_str(), 
                                        &FitFunc::RBWConvGausBGPol3, 
                                        massResonance - gammaResonance*3., 
                                        massResonance + gammaResonance*3., 8));
            altFitsBGAB.push_back(new TF1(("BGAB" + std::to_string(i)).c_str(), 
                                          &FitFunc::Pol3,
                                          massResonance - gammaResonance*3., 
                                          massResonance + gammaResonance*3., 4));

            altFitsFreeG.push_back(new TF1(("Free G" + std::to_string(i)).c_str(), 
                                            &FitFunc::RBWConvGausBGPol2, 
                                            massResonance - gammaResonance*3., 
                                            massResonance + gammaResonance*3., 7));
            altFitsBGFreeG.push_back(new TF1(("Free G BG" + std::to_string(i)).c_str(), 
                                              &FitFunc::Pol2, 
                                              massResonance - gammaResonance*3., 
                                              massResonance + gammaResonance*3., 3));

            altFitsFixedG.push_back(new TF1(("Fixed G" + std::to_string(i)).c_str(), 
                                             &FitFunc::RBWConvGausBGPol2, 
                                             massResonance - gammaResonance*3., 
                                             massResonance + gammaResonance*3., 7));
            altFitsBGFixedG.push_back(new TF1(("Fixed G BG" + std::to_string(i)).c_str(), 
                                               &FitFunc::Pol2, 
                                               massResonance - gammaResonance*3., 
                                               massResonance + gammaResonance*3., 3));
         }
      }
      else if (bgFitFunc == "pol3")
      {
         fits.push_back(new TF1(("Default" + std::to_string(i)).c_str(), 
                                &FitFunc::RBWConvGausBGPol3, 
                                massResonance - gammaResonance*3., 
                                massResonance + gammaResonance*3., 8));
         fitsBG.push_back(new TF1(("Default BG" + std::to_string(i)).c_str(), 
                                  &FitFunc::Pol3, 
                                  massResonance - gammaResonance*3., 
                                  massResonance + gammaResonance*3., 4));
         if (doAltFits)
         {
            altFitsAB.push_back(new TF1(("AB" + std::to_string(i)).c_str(), 
                                        &FitFunc::RBWConvGausBGPol2, 
                                        massResonance - gammaResonance*3., 
                                        massResonance + gammaResonance*3., 7));
            altFitsBGAB.push_back(new TF1(("BGAB" + std::to_string(i)).c_str(), 
                                          &FitFunc::Pol2,
                                          massResonance - gammaResonance*3., 
                                          massResonance + gammaResonance*3., 3));

            altFitsFreeG.push_back(new TF1(("Free G" + std::to_string(i)).c_str(), 
                                            &FitFunc::RBWConvGausBGPol3, 
                                            massResonance - gammaResonance*3., 
                                            massResonance + gammaResonance*3., 8));
            altFitsBGFreeG.push_back(new TF1(("Free G BG" + std::to_string(i)).c_str(), 
                                              &FitFunc::Pol3, 
                                              massResonance - gammaResonance*3., 
                                              massResonance + gammaResonance*3., 4));

            altFitsFixedG.push_back(new TF1(("Fixed G" + std::to_string(i)).c_str(), 
                                             &FitFunc::RBWConvGausBGPol3, 
                                             massResonance - gammaResonance*3., 
                                             massResonance + gammaResonance*3., 8));
            altFitsBGFixedG.push_back(new TF1(("Fixed G BG" + std::to_string(i)).c_str(), 
                                               &FitFunc::Pol3, 
                                               massResonance - gammaResonance*3., 
                                               massResonance + gammaResonance*3., 4));
         }
      }
      else if (bgFitFunc == "pol4")
      {
         fits.push_back(new TF1(("Default" + std::to_string(i)).c_str(), 
                                &FitFunc::RBWConvGausBGPol4, 
                                massResonance - gammaResonance*3., 
                                massResonance + gammaResonance*3., 9));
         fitsBG.push_back(new TF1(("Default BG" + std::to_string(i)).c_str(), 
                                  &FitFunc::Pol4, 
                                  massResonance - gammaResonance*3., 
                                  massResonance + gammaResonance*3., 5));
         if (doAltFits)
         {
            altFitsAB.push_back(new TF1(("AB" + std::to_string(i)).c_str(), 
                                        &FitFunc::RBWConvGausBGPol3, 
                                        massResonance - gammaResonance*3., 
                                        massResonance + gammaResonance*3., 8));
            altFitsBGAB.push_back(new TF1(("BGAB" + std::to_string(i)).c_str(), 
                                          &FitFunc::Pol3,
                                          massResonance - gammaResonance*3., 
                                          massResonance + gammaResonance*3., 4));

            altFitsFreeG.push_back(new TF1(("Free G" + std::to_string(i)).c_str(), 
                                            &FitFunc::RBWConvGausBGPol4, 
                                            massResonance - gammaResonance*3., 
                                            massResonance + gammaResonance*3., 9));
            altFitsBGFreeG.push_back(new TF1(("Free G BG" + std::to_string(i)).c_str(), 
                                              &FitFunc::Pol4, 
                                              massResonance - gammaResonance*3., 
                                              massResonance + gammaResonance*3., 5));

            altFitsFixedG.push_back(new TF1(("Fixed G" + std::to_string(i)).c_str(), 
                                             &FitFunc::RBWConvGausBGPol4, 
                                             massResonance - gammaResonance*3., 
                                             massResonance + gammaResonance*3., 9));
            altFitsBGFixedG.push_back(new TF1(("Fixed G BG" + std::to_string(i)).c_str(), 
                                               &FitFunc::Pol4, 
                                               massResonance - gammaResonance*3., 
                                               massResonance + gammaResonance*3., 5));
         }
      }
      else CppTools::PrintError("Unknown fit function specified in input file: " + bgFitFunc);

      const std::string pTBinRangeName =  
         CppTools::DtoStr(pTBinRanges[i], 2) + "<p_{T}<" + 
         CppTools::DtoStr(pTBinRanges[i + 1], 2);

      const bool isBGFixedForThisPT = SetBGFit(inputFileFitsBG, fitsBG.back(), pTBinRangeName);
      bool isBGFixedForThisPTAltFitAB = false;
      bool isBGFixedForThisPTAltFitFreeG = false;
      bool isBGFixedForThisPTAltFitFixedG = false;

      if (doAltFits)
      {
         isBGFixedForThisPTAltFitAB = 
            SetBGFit(inputFileFitsBGAB, altFitsBGAB.back(), pTBinRangeName);
         isBGFixedForThisPTAltFitFreeG = 
            SetBGFit(inputFileFitsBGFreeG, altFitsBGFreeG.back(), pTBinRangeName);
         isBGFixedForThisPTAltFitFixedG = 
            SetBGFit(inputFileFitsBGFixedG, altFitsBGFixedG.back(), pTBinRangeName);
      }

      const double maxBinVal = distrMInv->GetBinContent(distrMInv->GetMaximumBin());
      const double minBinVal = distrMInv->GetBinContent(distrMInv->GetMinimumBin());

      fits.back()->SetParameters(maxBinVal, massResonance, 
                                 gammaResonance, gaussianBroadeningSigma);

      fits.back()->SetParLimits(0, 1., maxBinVal - minBinVal);
      fits.back()->SetParLimits(1, massResonance/1.05, massResonance*1.05);
      fits.back()->SetParLimits(2, gammaResonance/1.10, gammaResonance*1.10);
      fits.back()->SetParLimits(3, gaussianBroadeningSigma/1.10, 
                                   gaussianBroadeningSigma*1.10);
      if (doAltFits)
      {
         altFitsAB.back()->SetParameters(maxBinVal, massResonance, 
                                         gammaResonance, gaussianBroadeningSigma);
         altFitsFreeG.back()->SetParameters(maxBinVal, massResonance, 
                                             gammaResonance, gaussianBroadeningSigma);
         altFitsFixedG.back()->SetParameters(maxBinVal, massResonance, 
                                              gammaResonance, gaussianBroadeningSigma);

         altFitsAB.back()->SetParLimits(0, 1., maxBinVal - minBinVal);
         altFitsAB.back()->SetParLimits(1, massResonance/1.05, massResonance*1.05);
         altFitsAB.back()->SetParLimits(2, gammaResonance/1.10, gammaResonance*1.10);
         altFitsAB.back()->SetParLimits(3, gaussianBroadeningSigma/1.10, 
                                           gaussianBroadeningSigma*1.10);

         altFitsFreeG.back()->SetParLimits(0, 1., maxBinVal - minBinVal);
         altFitsFreeG.back()->SetParLimits(1, massResonance/1.05, massResonance*1.05);
         altFitsFreeG.back()->SetParLimits(2, gammaResonance/2., gammaResonance*2.);
         altFitsFreeG.back()->SetParLimits(3, gaussianBroadeningSigma/1.10, 
                                            gaussianBroadeningSigma*1.10);

         altFitsFixedG.back()->SetParLimits(0, 1., maxBinVal - minBinVal);
         altFitsFixedG.back()->SetParLimits(1, massResonance/1.05, massResonance*1.05);
         altFitsFixedG.back()->FixParameter(2, gammaResonance);
         altFitsFixedG.back()->SetParLimits(3, gaussianBroadeningSigma/1.10, 
                                             gaussianBroadeningSigma*1.10);
      }

      if (isBGFixedForThisPT)
      {
         for (int i = fits.back()->GetNpar() - fitsBG.back()->GetNpar(); 
              i < fits.back()->GetNpar(); i++)
         {
            fits.back()->FixParameter(i, fitsBG.back()->GetParameter(i - fits.back()->GetNpar() + 
                                                                     fitsBG.back()->GetNpar()));
         }
      }

      if (isBGFixedForThisPTAltFitAB)
      {
         for (int i = altFitsAB.back()->GetNpar() - altFitsBGAB.back()->GetNpar(); 
              i < altFitsAB.back()->GetNpar(); i++)
         {
            altFitsAB.back()->FixParameter(i, altFitsBGAB.back()->
                                           GetParameter(i - altFitsAB.back()->GetNpar() + 
                                                        altFitsBGAB.back()->GetNpar()));
         }
      }

      if (isBGFixedForThisPTAltFitFreeG)
      {
         for (int i = altFitsFreeG.back()->GetNpar() - altFitsBGFreeG.back()->GetNpar(); 
              i < altFitsFreeG.back()->GetNpar(); i++)
         {
            altFitsFreeG.back()->FixParameter(i, altFitsBGFreeG.back()->
                                              GetParameter(i - altFitsFreeG.back()->GetNpar() + 
                                                           altFitsBGFreeG.back()->GetNpar()));
         }
      }

      if (isBGFixedForThisPTAltFitFixedG)
      {
         for (int i = altFitsFixedG.back()->GetNpar() - altFitsBGFixedG.back()->GetNpar(); 
              i < altFitsFixedG.back()->GetNpar(); i++)
         {
            altFitsFixedG.back()->FixParameter(i, altFitsBGFixedG.back()->
                                               GetParameter(i - altFitsFixedG.back()->GetNpar() + 
                                                            altFitsBGFixedG.back()->GetNpar()));
         }
      }

      distrMInv->Fit(fits.back(), "RQMNBLC");

      if (doAltFits)
      {
         distrMInv->Fit(altFitsAB.back(), "RQNMBLC");
         distrMInv->Fit(altFitsFreeG.back(), "RQNMBLC");
         distrMInv->Fit(altFitsFixedG.back(), "RQNMBLC");
      }

      if (!isBGFixedForThisPT)
      {
         for (int j = 0; j < fitsBG.back()->GetNpar(); j++)
         {
            fitsBG.back()->SetParameter(j, fits.back()->
                                        GetParameter(fits.back()->GetNpar() - 
                                                     fitsBG.back()->GetNpar() + j));
         }
      }

      if (doAltFits)
      {
         if (!isBGFixedForThisPTAltFitAB)
         {
            for (int j = 0; j < altFitsBGAB.back()->GetNpar(); j++)
            {
               altFitsBGAB.back()->
                  SetParameter(j, altFitsAB.back()->
                               GetParameter(altFitsAB.back()->GetNpar() - 
                                            altFitsBGAB.back()->GetNpar() + j));
            }
         }

         if (!isBGFixedForThisPTAltFitFreeG)
         {
            for (int j = 0; j < altFitsBGFreeG.back()->GetNpar(); j++)
            {
               altFitsBGFreeG.back()->
                  SetParameter(j, altFitsFreeG.back()->
                               GetParameter(altFitsFreeG.back()->GetNpar() - 
                                            altFitsBGFreeG.back()->GetNpar() + j));
            }
         }

         if (!isBGFixedForThisPTAltFitFixedG)
         {
            for (int j = 0; j < altFitsBGFreeG.back()->GetNpar(); j++)
            {
               altFitsBGFixedG.back()->
                  SetParameter(j, altFitsFixedG.back()->
                               GetParameter(altFitsFixedG.back()->GetNpar() - 
                                            altFitsBGFixedG.back()->GetNpar() + j));
            }
         }
      }

      fits.back()->SetRange(fits.back()->GetParameter(1) - 
                            (fits.back()->GetParameter(2) + gaussianBroadeningSigma)*
                            sigmalizedFitRange, 
                            fits.back()->GetParameter(1) + 
                            (fits.back()->GetParameter(2) + gaussianBroadeningSigma)*
                            sigmalizedFitRange);

      fitsBG.back()->SetRange(fits.back()->GetParameter(1) - 
                              (fits.back()->GetParameter(2) + gaussianBroadeningSigma)*
                              sigmalizedFitRange, 
                              fits.back()->GetParameter(1) + 
                              (fits.back()->GetParameter(2) + gaussianBroadeningSigma)*
                              sigmalizedFitRange);

      fits.back()->SetLineWidth(4);
      fitsBG.back()->SetLineWidth(4);
      fitsBG.back()->SetLineStyle(7);

      if (doAltFits)
      {
         altFitsAB.back()->SetRange(altFitsAB.back()->GetParameter(1) - 
                                    (altFitsAB.back()->GetParameter(2) + 
                                     gaussianBroadeningSigma)*sigmalizedFitRange, 
                                    altFitsAB.back()->GetParameter(1) + 
                                    (altFitsAB.back()->GetParameter(2) + gaussianBroadeningSigma)*
                                    sigmalizedFitRange);

         altFitsBGAB.back()->SetRange(altFitsAB.back()->GetParameter(1) - 
                                      (altFitsAB.back()->GetParameter(2) + 
                                       gaussianBroadeningSigma)*sigmalizedFitRange, 
                                      altFitsAB.back()->GetParameter(1) + 
                                      (altFitsAB.back()->GetParameter(2) + 
                                       gaussianBroadeningSigma)*sigmalizedFitRange);

         altFitsFreeG.back()->SetRange(altFitsFreeG.back()->GetParameter(1) - 
                                       (altFitsFreeG.back()->GetParameter(2) + 
                                        gaussianBroadeningSigma)*sigmalizedFitRange, 
                                       altFitsFreeG.back()->GetParameter(1) + 
                                       (altFitsFreeG.back()->GetParameter(2) + 
                                        gaussianBroadeningSigma)*sigmalizedFitRange);

         altFitsBGFreeG.back()->SetRange(altFitsFreeG.back()->GetParameter(1) - 
                                         (altFitsFreeG.back()->GetParameter(2) + 
                                          gaussianBroadeningSigma)*sigmalizedFitRange, 
                                         altFitsFreeG.back()->GetParameter(1) + 
                                         (altFitsFreeG.back()->GetParameter(2) + 
                                          gaussianBroadeningSigma)*sigmalizedFitRange);

         altFitsFixedG.back()->SetRange(altFitsFixedG.back()->GetParameter(1) - 
                                        (altFitsFixedG.back()->GetParameter(2) + 
                                         gaussianBroadeningSigma)*sigmalizedFitRange, 
                                        altFitsFixedG.back()->GetParameter(1) + 
                                        (altFitsFixedG.back()->GetParameter(2) + 
                                         gaussianBroadeningSigma)*sigmalizedFitRange);

         altFitsBGFixedG.back()->SetRange(altFitsFixedG.back()->GetParameter(1) - 
                                          (altFitsFixedG.back()->GetParameter(2) + 
                                           gaussianBroadeningSigma)*sigmalizedFitRange, 
                                          altFitsFixedG.back()->GetParameter(1) + 
                                          (altFitsFixedG.back()->GetParameter(2) + 
                                           gaussianBroadeningSigma)*sigmalizedFitRange);

         altFitsAB.back()->SetLineWidth(4);
         altFitsFreeG.back()->SetLineWidth(4);
         altFitsFixedG.back()->SetLineWidth(4);

         altFitsBGAB.back()->SetLineWidth(4);
         altFitsBGFreeG.back()->SetLineWidth(4);
         altFitsBGFixedG.back()->SetLineWidth(4);

         altFitsBGAB.back()->SetLineStyle(7);
         altFitsBGFreeG.back()->SetLineStyle(7);
         altFitsBGFixedG.back()->SetLineStyle(7);
      }

      distrMInv->SetMaximum(distrMInv->GetMaximum()*1.2);

      distrMInv->GetXaxis()->SetRange(distrMInv->GetXaxis()->FindBin(minMInv + 1e-7), 
                                      distrMInv->GetXaxis()->FindBin(maxMInv - 1e-7));
      distrMInvFG->GetXaxis()->SetRange(distrMInvFG->GetXaxis()->FindBin(minMInv + 1e-7), 
                                        distrMInvFG->GetXaxis()->FindBin(maxMInv - 1e-7));
      distrMInvBG->GetXaxis()->SetRange(distrMInvBG->GetXaxis()->FindBin(minMInv + 1e-7), 
                                        distrMInvBG->GetXaxis()->FindBin(maxMInv - 1e-7));

      distrMInv->SetLineWidth(2);
      distrMInvFG->SetLineWidth(3);
      distrMInvBG->SetLineWidth(2);

      distrMInv->SetLineColor(kGray + 3);
      distrMInv->SetMarkerColor(kGray + 3);

      distrMInv->SetMarkerStyle(20);
      distrMInv->SetMarkerSize(0.7);
      distrMInv->SetMarkerColor(kGray + 3);

      distrsMInv.emplace_back(static_cast<TH1D *>(distrMInv->Clone()));
      distrMInvNames.emplace_back((CppTools::DtoStr(pTBinRanges[i], 2) + "<p_{T}<" + 
                                   CppTools::DtoStr(pTBinRanges[i + 1], 2)));
      pTsDistrMInv.push_back((pTBinRanges[i] + pTBinRanges[i + 1])/2.);

      numberOfCalls++;
   }
}

// copied from AnalyzeRealMInv.cpp
void AnalyzeRealMInv::SetGaussianBroadeningFunction()
{
   const std::string inputFileName = "data/Parameters/GaussianBroadening/" + 
                                     runName + "/" + resonanceName + ".root";

   if (!std::filesystem::exists(inputFileName))
   {
      CppTools::PrintError(inputFileName + " does not exists. "\
                           "Run executable bin/EstimateGassianBroadening first");
   }

   gaussianBroadeningEstimatorFunc = 
      static_cast<TF1 *>(TFile::Open(inputFileName.c_str())->Get("gaussian broadening sigma fit"));
}

TFile *AnalyzeRealMInv::SetFixedBGFile(const std::string& inputFileName, 
                                       const std::string& fitTypeName)
{
   if (std::filesystem::exists(inputFileName))
   {
      pBar.Clear();
      CppTools::PrintInfo("Fixed BG fits file for " + fitTypeName + " fits was found");
      pBar.RePrint();
      return TFile::Open(inputFileName.c_str());
   }
   else
   {
      pBar.Clear();
      CppTools::PrintWarning("Fixed BG fits file for " + fitTypeName + " fits was not found; "\
                             "using free BG fit for these fits");
      pBar.RePrint();
      return nullptr;
   }
}

bool AnalyzeRealMInv::SetBGFit(TFile *&inputFile, TF1 *&fitBG, const std::string& identifier)
{
   if (inputFile)
   {
      TF1 *fitBGTmp = static_cast<TF1 *>(inputFile->Get(identifier.c_str()));

      if (!fitBGTmp)
      {
         pBar.Clear();
         CppTools::PrintWarning("Could not obtain BG approximation for + " + 
                                identifier + " from file " + std::string(inputFile->GetName()) + 
                                "; fixed BG will be disabled for this fit type and pT bin");
         pBar.RePrint();
      }
      else if (fitBGTmp->GetNpar() != fitBG->GetNpar())
      {
         pBar.Clear();
         CppTools::PrintWarning("Number of parameters for BG fit mismatch for " +
                                identifier + " from file " + std::string(inputFile->GetName()) + 
                                "; fixed BG will be disabled for this fit type and pT bin");
         pBar.RePrint();
      }
      else
      {
         for (int i = 0; i < fitBGTmp->GetNpar(); i++)
         {
            fitBG->SetParameter(i, fitBGTmp->GetParameter(i));
         }
         return true;
      }
   }

   return false;
}
