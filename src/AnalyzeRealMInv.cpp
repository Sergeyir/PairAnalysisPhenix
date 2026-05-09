/** 
 *  @file   AnalyzeRealMInv.cpp 
 *  @brief  Contains realisations of functions that are used for analyzis of expeirmental invariant mass distributions 
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysis).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef ANALYZE_REAL_M_INV_CPP
#define ANALYZE_REAL_M_INV_CPP

#include "AnalyzeRealMInv.hpp"

using namespace AnalyzeRealMInv;

int main(int argc, char **argv)
{
   if (argc < 3 || argc > 6) 
   {
      CppTools::PrintError("Expected 2-5 parameters while " + std::to_string(argc - 1) + " "\
                           "parameter(s) were provided \n Usage: bin/AnalyzeRealMInv "\
                           "inputYAMLName taxiNumber methodName=all rebinX=1 "\
                           "numberOfThreads=std::thread::hardware_concurrency()");
   }
 
   CppTools::CheckInputFile(argv[1]);

   taxiNumber = std::stoi(argv[2]);
 
   if (argc == 6) ROOT::EnableImplicitMT(std::stoi(argv[5]));
   else ROOT::EnableImplicitMT(std::thread::hardware_concurrency());

   if (argc == 5) 
   {
      rebinX = std::stoi(argv[4]);
      if (rebinX <= 0) CppTools::PrintError("Rebin value cannot be 0 or negative");
   }

   std::string methodToAnalyze;
   if (argc >= 4) methodToAnalyze = argv[3];
   else methodToAnalyze = "all";

   gStyle->SetOptStat(0);
   gErrorIgnoreLevel = kWarning;

   gStyle->SetPalette(kRainbow);

   TH1::AddDirectory(kFALSE);
   TH2::AddDirectory(kFALSE);

   inputYAMLResonance.OpenFile(argv[1]);
   inputYAMLResonance.CheckStatus("resonance");

   runName = inputYAMLResonance["run_name"].as<std::string>();

   inputYAMLMain.OpenFile("input/" + runName + "/main.yaml");
   inputYAMLMain.CheckStatus("main");

   resonanceName = inputYAMLResonance["name"].as<std::string>();
   massResonance = inputYAMLResonance["mass"].as<double>();
   gammaResonance = inputYAMLResonance["gamma"].as<double>();

   daughter1Id = inputYAMLResonance["daughter1_id"].as<int>();
   daughter2Id = inputYAMLResonance["daughter2_id"].as<int>();

   minMInv = inputYAMLResonance["m_inv_range_min"].as<double>();
   maxMInv = inputYAMLResonance["m_inv_range_max"].as<double>();

   SetGaussianBroadeningFunction();

   inputFileName = "data/Real/" + runName + "/Resonance/" + std::to_string(taxiNumber) + ".root";

   CppTools::CheckInputFile(inputFileName);
   inputFile = TFile::Open(inputFileName.c_str(), "READ");

   text.SetTextFont(43);
   text.SetTextSize(45);

   texText.SetTextFont(43);
   texText.SetTextSize(45);

   text.SetTextAngle(270.);

   pTNBins = inputYAMLResonance["pt_bins"].size();

   for (unsigned int i = 0; i < pTNBins; i++)
   {
      pTBinRanges.push_back(inputYAMLResonance["pt_bins"][i]["min"].as<double>());
   }
   pTBinRanges.push_back(inputYAMLResonance["pt_bins"][pTNBins - 1]["max"].as<double>());

   sigmalizedYieldExtractionRange = 
      inputYAMLResonance["sigmalized_yield_extraction_range"].as<double>();

   int methodIndex = -1;
   if (methodToAnalyze == "all")
   {
      for (const YAML::Node& method : inputYAMLResonance["pair_selection_methods"])
      {
         numberOfIterations += pTNBins*inputYAMLResonance["centrality_bins"].size();
      }
   }
   else
   {
      for (int i = 0; i < static_cast<int>(inputYAMLResonance["pair_selection_methods"].size()); i++)
      {
         if (inputYAMLResonance["pair_selection_methods"][i]["name"].as<std::string>() ==
             methodToAnalyze) methodIndex = i;
      }
      if (methodIndex == -1)
      {
         CppTools::PrintError("No method named " + methodToAnalyze + " in input file");
      }

      numberOfIterations = pTNBins*inputYAMLResonance["centrality_bins"].size();
   }

   parametersOutputDir = "data/RawYields/" + runName + "/Resonance";
   std::filesystem::create_directories(parametersOutputDir);

   // performing fits for specified pair selection methods
   if (methodToAnalyze == "all")
   {
      for (const YAML::Node& method : inputYAMLResonance["pair_selection_methods"])
      {
         PerformMInvFits(method);
      }
   }
   else 
   {
      PerformMInvFits(inputYAMLResonance["pair_selection_methods"][methodIndex]);
   }

   parametersOutputFile->Close();

   pBar.Finish();

   CppTools::PrintInfo("AnalyzeRealMInv executable has finished running succesfully");
}

void AnalyzeRealMInv::PerformMInvFits(const YAML::Node& method)
{
   const std::string methodName = method["name"].as<std::string>();

   parametersOutputFile = TFile::Open((parametersOutputDir + "/" + std::to_string(taxiNumber) + 
                                       "_" + resonanceName + "_" + methodName + 
                                       ".root").c_str(), "RECREATE");

   const std::string outputDir = "output/MInv/" + runName + "/" + 
                                 std::to_string(taxiNumber) + "/" + methodName;
   std::filesystem::create_directories(outputDir);

   const double sigmalizedFitRange = method["sigmalized_fit_range"].as<double>();

   for (unsigned int centralityBin = 0; 
        centralityBin < inputYAMLResonance["centrality_bins"].size(); centralityBin++)
   {
      const YAML::Node centrality = inputYAMLResonance["centrality_bins"][centralityBin];

      const std::string centralityName = centrality["name"].as<std::string>();

      TH1D distrMeansVsPT("means vs pT", "", pTNBins, &pTBinRanges[0]);
      TH1D distrGammasVsPT("gammas vs pT", "", pTNBins, &pTBinRanges[0]);
      TH1D distrRawYieldVsPTStatErr("raw yield vs pT with stat errors", 
                                    "", pTNBins, &pTBinRanges[0]);
      TH1D distrRawYieldVsPTSysErr("raw yield vs pT with sys errors", 
                                   "", pTNBins, &pTBinRanges[0]);

      std::vector<std::vector<double>> vecParBG;
      std::vector<double> vecPT;

      const std::string inputFileFitsBGName = "data/Parameters/BGFitResonance/" + runName + "/" + 
                                              std::to_string(taxiNumber) + "/" + methodName + 
                                              "_" + centralityName + ".root";

      // input file with BG fits for default fit
      TFile *inputFileFitsBG = 
         SetFixedBGFile("data/Parameters/BGFitResonance/" + runName + "/" + 
                        std::to_string(taxiNumber) + "/" + methodName + "_" + centralityName + 
                        ".root", "default " + methodName + " " + centralityName);

      TFile *inputFileFitsBGAB = 
         SetFixedBGFile("data/Parameters/BGFitResonance/" + runName + "/" + 
                        std::to_string(taxiNumber) + "/" + methodName + "_" + 
                        centralityName + "_AB.root", "alternative BG " + 
                        methodName + " " + centralityName, false); 
      TFile *inputFileFitsBGFreeG = 
         SetFixedBGFile("data/Parameters/BGFitResonance/" + runName + "/" + 
                        std::to_string(taxiNumber) + "/" + methodName + "_" + 
                        centralityName + "_FreeG.root", "free G " + 
                        methodName + " " + centralityName, false); 
      TFile *inputFileFitsBGFixedG = 
         SetFixedBGFile("data/Parameters/BGFitResonance/" + runName + "/" + 
                        std::to_string(taxiNumber) + "/" + methodName + "_" + 
                        centralityName + "_FixedG.root", "fixed BG " + 
                        methodName + " " + centralityName, false); 

      // only performing alternative fits if all parameters were found:
      // there are runs with alternative cuts for which no alternative fits are needed;
      // also alternative fits most of the time require tweaking with GUI/MInv.cpp
      // which writes the background approximation parameters
      const bool performAltFits = (inputFileFitsBGAB && 
                                   inputFileFitsBGFreeG && 
                                   inputFileFitsBGFixedG);

      const unsigned int pTBinFitMin = method["centrality_bin_parameters"][centralityBin]
                                             ["pt_bin_min"].as<int>();
      const unsigned int pTBinFitMax = method["centrality_bin_parameters"][centralityBin]
                                             ["pt_bin_max"].as<int>();

      for (unsigned int i = 0; i < pTNBins; i++)
      {
         pBar.Print(static_cast<double>(numberOfCalls)/static_cast<double>(numberOfIterations));

         TH1D *distrMInvFG = nullptr;
         TH1D *distrMInvBG = nullptr;
         TH1D *distrMInvFGLR = nullptr; // low resolution (for BG normalization)
         TH1D *distrMInvBGLR = nullptr; // low resolution (for BG normalization)

         // shows whether to perform approximations for the current bin
         // if false the histograms will not be approximated but will be printed anyways
         const bool performFit = (i >= pTBinFitMin && i <= pTBinFitMax);

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

            numberOfEvents /= 2.;
         }

         if (!distrMInv)
         {
            pBar.Clear();
            CppTools::PrintError("Resulting M_{inv} histogram could not be constructed for " + 
                                 centrality["name"].as<std::string>() + " " +
                                 CppTools::DtoStr(pTBinRanges[i], 2) + "<pT<" + 
                                 CppTools::DtoStr(pTBinRanges[i + 1], 2));
         }
         else if (performFit && distrMInv->Integral(1, distrMInv->GetXaxis()->GetNbins()) < 1e-7)
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
         else if (performFit && 
                  distrMInvBG->Integral(1, distrMInvBG->GetXaxis()->GetNbins()) < 1e-7)
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

         // fit for resonance+bg approximation
         TF1 *fit = nullptr;
         // fit for bg approximation
         TF1 *fitBG = nullptr;
         // alternatice fit for resonance+bg approximation with alternative BG
         TF1 *altFitAB = nullptr;
         TF1 *altFitBGAB = nullptr;
         // alternatice fit for resonance+bg approximation with free Gamma
         TF1 *altFitFreeG = nullptr;
         TF1 *altFitBGFreeG = nullptr;
         // alternatice fit for resonance+bg approximation with fixed Gamma
         TF1 *altFitFixedG = nullptr;
         TF1 *altFitBGFixedG = nullptr;

         // sigma of a gaus that is convoluted with Breit-Wigner
         const double gaussianBroadeningSigma = 
            gaussianBroadeningEstimatorFunc->Eval((pTBinRanges[i] + pTBinRanges[i + 1])/2.);

         double fitRangeMin = -1.;
         double fitRangeMax = -1.;
         double altFitABRangeMin = -1.;
         double altFitABRangeMax = -1.;
         double altFitFreeGRangeMin = -1.;
         double altFitFreeGRangeMax = -1.;
         double altFitFixedGRangeMin = -1.;
         double altFitFixedGRangeMax = -1.;

         double lowIntegrationRange = -1.;
         double upIntegrationRange = -1.;
         double lowIntegrationRangeAltFitAB = -1.;
         double upIntegrationRangeAltFitAB = -1.;
         double lowIntegrationRangeAltFitFreeG = -1.;
         double upIntegrationRangeAltFitFreeG = -1.;
         double lowIntegrationRangeAltFitFixedG = -1.;
         double upIntegrationRangeAltFitFixedG = -1.;

         if (performFit)
         {
            const std::string bgFitFunc = method["bg_fit_func"].as<std::string>();
            if (bgFitFunc == "pol2")
            {
               fit = new TF1("Default", &FitFunc::RBWConvGausBGPol2, 
                             massResonance - gammaResonance*3., 
                             massResonance + gammaResonance*3., 7);
               fitBG = new TF1("Default BG", &FitFunc::Pol2, 
                               massResonance - gammaResonance*3., 
                               massResonance + gammaResonance*3., 3);

               if (performAltFits)
               {
                  altFitAB = new TF1("AB", &FitFunc::RBWConvGausBGPol3, 
                                     massResonance - gammaResonance*3., 
                                     massResonance + gammaResonance*3., 8);
                  altFitBGAB = new TF1("AB BG", &FitFunc::Pol3, 
                                       massResonance - gammaResonance*3., 
                                       massResonance + gammaResonance*3., 4);
                  altFitFreeG = new TF1("Free #Gamma", &FitFunc::RBWConvGausBGPol2, 
                                        massResonance - gammaResonance*3., 
                                        massResonance + gammaResonance*3., 7);
                  altFitBGFreeG = new TF1("Free #Gamma BG", &FitFunc::Pol2, 
                                          massResonance - gammaResonance*3., 
                                          massResonance + gammaResonance*3., 3);
                  altFitFixedG = new TF1("Fixed #Gamma", &FitFunc::RBWConvGausBGPol2, 
                                         massResonance - gammaResonance*3., 
                                         massResonance + gammaResonance*3., 7);
                  altFitBGFixedG = new TF1("Fixed #Gamma BG", &FitFunc::Pol2, 
                                           massResonance - gammaResonance*3., 
                                           massResonance + gammaResonance*3., 3);
               }
            }
            else if (bgFitFunc == "pol3")
            {
               fit = new TF1("resonance + bg fit", &FitFunc::RBWConvGausBGPol3, 
                             massResonance - gammaResonance*3., 
                             massResonance + gammaResonance*3., 8);
               fitBG = new TF1("bg fit", &FitFunc::Pol3, 
                               massResonance - gammaResonance*3., 
                               massResonance + gammaResonance*3., 4);
               if (performAltFits)
               {
                  altFitAB = new TF1("AB", &FitFunc::RBWConvGausBGPol2, 
                                     massResonance - gammaResonance*3., 
                                     massResonance + gammaResonance*3., 7);
                  altFitBGAB = new TF1("AB BG", &FitFunc::Pol2, 
                                       massResonance - gammaResonance*3., 
                                       massResonance + gammaResonance*3., 3);
                  altFitFreeG = new TF1("FreeG", &FitFunc::RBWConvGausBGPol3, 
                                        massResonance - gammaResonance*3., 
                                        massResonance + gammaResonance*3., 8);
                  altFitBGFreeG = new TF1("FreeG BG", &FitFunc::Pol3, 
                                          massResonance - gammaResonance*3., 
                                          massResonance + gammaResonance*3., 4);
                  altFitFixedG = new TF1("FixedG", &FitFunc::RBWConvGausBGPol3, 
                                         massResonance - gammaResonance*3., 
                                         massResonance + gammaResonance*3., 8);
                  altFitBGFixedG = new TF1("FixedG BG", &FitFunc::Pol3, 
                                           massResonance - gammaResonance*3., 
                                           massResonance + gammaResonance*3., 4);
               }
            }
            else if (bgFitFunc == "pol4")
            {
               fit = new TF1("resonance + bg fit", &FitFunc::RBWConvGausBGPol4, 
                             massResonance - gammaResonance*3., 
                             massResonance + gammaResonance*3., 9);
               fitBG = new TF1("bg fit", &FitFunc::Pol4, 
                               massResonance - gammaResonance*3., 
                               massResonance + gammaResonance*3., 5);
               if (performAltFits)
               {
                  altFitAB = new TF1("AB", &FitFunc::RBWConvGausBGPol3, 
                                     massResonance - gammaResonance*3., 
                                     massResonance + gammaResonance*3., 8);
                  altFitBGAB = new TF1("AB BG", &FitFunc::Pol3, 
                                       massResonance - gammaResonance*3., 
                                       massResonance + gammaResonance*3., 4);
                  altFitFreeG = new TF1("Free #Gamma", &FitFunc::RBWConvGausBGPol4, 
                                        massResonance - gammaResonance*3., 
                                        massResonance + gammaResonance*3., 9);
                  altFitBGFreeG = new TF1("Free #Gamma BG", &FitFunc::Pol4, 
                                          massResonance - gammaResonance*3., 
                                          massResonance + gammaResonance*3., 5);
                  altFitFixedG = new TF1("Fixed #Gamma", &FitFunc::RBWConvGausBGPol4, 
                                         massResonance - gammaResonance*3., 
                                         massResonance + gammaResonance*3., 9);
                  altFitBGFixedG = new TF1("Fixed #Gamma BG", &FitFunc::Pol4, 
                                           massResonance - gammaResonance*3., 
                                           massResonance + gammaResonance*3., 5);
               }
            }
            else CppTools::PrintError("Unknown fit function specified in input file: " + bgFitFunc);

            const std::string pTBinRangeName =  
               CppTools::DtoStr(pTBinRanges[i], 2) + "<p_{T}<" + 
               CppTools::DtoStr(pTBinRanges[i + 1], 2);

            const bool isBGFixedForThisPT = 
               SetBGFit(inputFileFitsBG, fitBG, pTBinRangeName);
            // only performing alternative fits if all parameters were succesfully read
            const bool isBGFixedForThisPTAltFit = 
               (SetBGFit(inputFileFitsBGAB, altFitBGAB, pTBinRangeName) &&
                SetBGFit(inputFileFitsBGFreeG, altFitBGFreeG, pTBinRangeName) &&
                SetBGFit(inputFileFitsBGFixedG, altFitBGFixedG, pTBinRangeName));

            const double maxBinVal = distrMInv->GetBinContent(distrMInv->GetMaximumBin());
            const double minBinVal = distrMInv->GetBinContent(distrMInv->GetMinimumBin());

            fit->SetParameters(maxBinVal, massResonance, gammaResonance, gaussianBroadeningSigma);
            fit->SetParLimits(0, 1., maxBinVal - minBinVal);
            fit->SetParLimits(1, massResonance/1.05, massResonance*1.05);
            fit->SetParLimits(2, gammaResonance/1.10, gammaResonance*1.10);
            fit->SetParLimits(3, gaussianBroadeningSigma/1.10, gaussianBroadeningSigma*1.10);

            if (isBGFixedForThisPTAltFit)
            {
               altFitAB->SetParameters(maxBinVal, massResonance, 
                                       gammaResonance, gaussianBroadeningSigma);
               altFitAB->SetParLimits(0, 1., maxBinVal - minBinVal);
               altFitAB->SetParLimits(1, massResonance/1.05, massResonance*1.05);
               altFitAB->SetParLimits(2, gammaResonance/1.10, gammaResonance*1.10);
               altFitAB->SetParLimits(3, gaussianBroadeningSigma/1.10, 
                                      gaussianBroadeningSigma*1.10);

               altFitFreeG->SetParameters(maxBinVal, massResonance, 
                                          gammaResonance, gaussianBroadeningSigma);
               altFitFreeG->SetParLimits(0, 1., maxBinVal - minBinVal);
               altFitFreeG->SetParLimits(1, massResonance/1.05, massResonance*1.05);
               altFitFreeG->SetParLimits(2, gammaResonance/2., gammaResonance*2.);
               altFitFreeG->SetParLimits(3, gaussianBroadeningSigma/1.10, 
                                         gaussianBroadeningSigma*1.10);

               altFitFixedG->SetParameters(maxBinVal, massResonance, 
                                           gammaResonance, gaussianBroadeningSigma);
               altFitFixedG->SetParLimits(0, 1., maxBinVal - minBinVal);
               altFitFixedG->SetParLimits(1, massResonance/1.05, massResonance*1.05);
               altFitFixedG->FixParameter(2, gammaResonance);
               altFitFixedG->FixParameter(3, gaussianBroadeningSigma);
            }

            if (isBGFixedForThisPT)
            {
               for (int i = fit->GetNpar() - fitBG->GetNpar(); i < fit->GetNpar(); i++)
               {
                  fit->FixParameter(i, fitBG->GetParameter(i - fit->GetNpar() + fitBG->GetNpar()));
               }
            }

            if (isBGFixedForThisPTAltFit)
            {
               for (int i = altFitAB->GetNpar() - altFitBGAB->GetNpar(); 
                    i < altFitAB->GetNpar(); i++)
               {
                  altFitAB->FixParameter(i, altFitBGAB->
                                         GetParameter(i - altFitAB->GetNpar() + 
                                                      altFitBGAB->GetNpar()));
               }
               for (int i = altFitFreeG->GetNpar() - altFitBGFreeG->GetNpar(); 
                    i < altFitFreeG->GetNpar(); i++)
               {
                  altFitFreeG->FixParameter(i, altFitBGFreeG->
                                            GetParameter(i - altFitFreeG->GetNpar() + 
                                                         altFitBGFreeG->GetNpar()));
               }
               for (int i = altFitFixedG->GetNpar() - altFitBGFixedG->GetNpar(); 
                    i < altFitFixedG->GetNpar(); i++)
               {
                  altFitFixedG->FixParameter(i, altFitBGFixedG->
                                             GetParameter(i - altFitFixedG->GetNpar() + 
                                                          altFitBGFixedG->GetNpar()));
               }
            }

            distrMInv->Fit(fit, "RQMNBLC");

            if (isBGFixedForThisPTAltFit)
            {
               distrMInv->Fit(altFitAB, "RQMNBLC");
               distrMInv->Fit(altFitFreeG, "RQMNBLC");
               distrMInv->Fit(altFitFixedG, "RQMNBLC");
            }

            for (unsigned int j = 1; j <= fitNTries; j++)
            {
               fitRangeMin = fit->GetParameter(1) - (fit->GetParameter(2) + 
                                                     fit->GetParameter(3))*sigmalizedFitRange;
               fitRangeMax = fit->GetParameter(1) + (fit->GetParameter(2) + 
                                                     fit->GetParameter(3))*sigmalizedFitRange;

               fit->SetParLimits(1, fit->GetParameter(1)/(1. + 0.05/static_cast<double>(j*j)), 
                                 fit->GetParameter(1)*(1. + 0.05/static_cast<double>(j*j)));
               fit->SetParLimits(2, fit->GetParameter(2)/(1. + 0.05/static_cast<double>(j*j)),
                                 fit->GetParameter(2)*(1. + 0.1/static_cast<double>(j*j)));

               fit->SetRange(fitRangeMin, fitRangeMax);

               distrMInv->Fit(fit, "RQMNBLC");

               if (isBGFixedForThisPTAltFit)
               {
                  altFitABRangeMin = altFitAB->GetParameter(1) - 
                                     (altFitAB->GetParameter(2) + 
                                      altFitAB->GetParameter(3))*sigmalizedFitRange;
                  altFitABRangeMax = altFitAB->GetParameter(1) + 
                                     (altFitAB->GetParameter(2) + 
                                      altFitAB->GetParameter(3))*sigmalizedFitRange;
                  altFitFreeGRangeMin = altFitFreeG->GetParameter(1) - 
                                        (altFitFreeG->GetParameter(2) + 
                                         altFitFreeG->GetParameter(3))*sigmalizedFitRange;
                  altFitFreeGRangeMax = altFitFreeG->GetParameter(1) + 
                                        (altFitFreeG->GetParameter(2) + 
                                         altFitFreeG->GetParameter(3))*sigmalizedFitRange;
                  altFitFixedGRangeMin = altFitFixedG->GetParameter(1) - 
                                         (altFitFixedG->GetParameter(2) + 
                                          altFitFixedG->GetParameter(3))*sigmalizedFitRange;
                  altFitFixedGRangeMax = altFitFixedG->GetParameter(1) + 
                                         (altFitFixedG->GetParameter(2) + 
                                          altFitFixedG->GetParameter(3))*sigmalizedFitRange;

                  altFitAB->SetParLimits(1, altFitAB->GetParameter(1)/
                                         (1. + 0.05/static_cast<double>(j*j)), 
                                         altFitAB->GetParameter(1)*
                                         (1. + 0.05/static_cast<double>(j*j)));
                  altFitAB->SetParLimits(2, altFitAB->GetParameter(2)/
                                         (1. + 0.05/static_cast<double>(j*j)),
                                         altFitAB->GetParameter(2)*
                                         (1. + 0.1/static_cast<double>(j*j)));
                  altFitFreeG->SetParLimits(1, altFitFreeG->GetParameter(1)/
                                            (1. + 0.05/static_cast<double>(j*j)), 
                                            altFitFreeG->GetParameter(1)*
                                            (1. + 0.05/static_cast<double>(j*j)));
                  altFitFreeG->SetParLimits(2, altFitFreeG->GetParameter(2)/
                                            (1. + 0.05/static_cast<double>(j*j)),
                                            altFitFreeG->GetParameter(2)*
                                            (1. + 0.1/static_cast<double>(j*j)));
                  altFitFixedG->SetParLimits(1, altFitFixedG->GetParameter(1)/
                                             (1. + 0.05/static_cast<double>(j*j)), 
                                             altFitFixedG->GetParameter(1)*
                                             (1. + 0.05/static_cast<double>(j*j)));
                  altFitFixedG->SetParLimits(2, altFitFixedG->GetParameter(2)/
                                             (1. + 0.05/static_cast<double>(j*j)),
                                             altFitFixedG->GetParameter(2)*
                                             (1. + 0.1/static_cast<double>(j*j)));

                  altFitAB->SetRange(altFitABRangeMin, altFitABRangeMax);
                  altFitFreeG->SetRange(altFitFreeGRangeMin, altFitFreeGRangeMax);
                  altFitFixedG->SetRange(altFitFixedGRangeMin, altFitFixedGRangeMax);

                  distrMInv->Fit(altFitAB, "RQMNBLC");
                  distrMInv->Fit(altFitFreeG, "RQMNBLC");
                  distrMInv->Fit(altFitFixedG, "RQMNBLC");
               }
            }
            //distrMInv->Fit(&fit, "RQMNBLE");


            fitRangeMin = fit->GetParameter(1) - (fit->GetParameter(2) + 
                                                  fit->GetParameter(3))*sigmalizedFitRange;
            fitRangeMax = fit->GetParameter(1) + (fit->GetParameter(2) + 
                                                  fit->GetParameter(3))*sigmalizedFitRange;

            if (isBGFixedForThisPTAltFit)
            {
               altFitABRangeMin = altFitAB->GetParameter(1) - 
                                  (altFitAB->GetParameter(2) + 
                                   altFitAB->GetParameter(3))*sigmalizedFitRange;
               altFitABRangeMax = altFitAB->GetParameter(1) + 
                                  (altFitAB->GetParameter(2) + 
                                   altFitAB->GetParameter(3))*sigmalizedFitRange;
               altFitFreeGRangeMin = altFitFreeG->GetParameter(1) - 
                                     (altFitFreeG->GetParameter(2) + 
                                      altFitFreeG->GetParameter(3))*sigmalizedFitRange;
               altFitFreeGRangeMax = altFitFreeG->GetParameter(1) + 
                                     (altFitFreeG->GetParameter(2) + 
                                      altFitFreeG->GetParameter(3))*sigmalizedFitRange;
               altFitFixedGRangeMin = altFitFixedG->GetParameter(1) - 
                                      (altFitFixedG->GetParameter(2) + 
                                       altFitFixedG->GetParameter(3))*sigmalizedFitRange;
               altFitFixedGRangeMax = altFitFixedG->GetParameter(1) + 
                                      (altFitFixedG->GetParameter(2) + 
                                       altFitFixedG->GetParameter(3))*sigmalizedFitRange;
            }

            distrMeansVsPT.SetBinContent(i + 1, fit->GetParameter(1));
            distrMeansVsPT.SetBinError(i + 1, fit->GetParError(1));

            distrGammasVsPT.SetBinContent(i + 1, fit->GetParameter(2));
            distrGammasVsPT.SetBinError(i + 1, fit->GetParError(2));

            if (!isBGFixedForThisPT)
            {
               for (int j = 0; j < fitBG->GetNpar(); j++)
               {
                  fitBG->SetParameter(j, fit->GetParameter(fit->GetNpar() - fitBG->GetNpar() + j));
               }
            }

            if (performAltFits && !isBGFixedForThisPTAltFit)
            {
               for (int j = 0; j < altFitBGAB->GetNpar(); j++)
               {
                  altFitBGAB->
                     SetParameter(j, altFitAB->GetParameter(altFitAB->GetNpar() - 
                                                            altFitBGAB->GetNpar() + j));
               }
               for (int j = 0; j < altFitBGFreeG->GetNpar(); j++)
               {
                  altFitBGFreeG->
                     SetParameter(j, altFitFreeG->GetParameter(altFitFreeG->GetNpar() - 
                                                               altFitBGFreeG->GetNpar() + j));
               }
               for (int j = 0; j < altFitBGFixedG->GetNpar(); j++)
               {
                  altFitBGFixedG->
                     SetParameter(j, altFitFixedG->GetParameter(altFitFixedG->GetNpar() - 
                                                                altFitBGFixedG->GetNpar() + j));
               }
            }

            fitBG->SetRange(fitRangeMin, fitRangeMax);

            if (isBGFixedForThisPTAltFit)
            {
               altFitBGAB->SetRange(altFitABRangeMin, altFitABRangeMax);
               altFitBGFreeG->SetRange(altFitFreeGRangeMin, altFitFreeGRangeMax);
               altFitBGFixedG->SetRange(altFitFixedGRangeMin, altFitFixedGRangeMax);
            }

            lowIntegrationRange = fit->GetParameter(1) - 
                                  (fit->GetParameter(2) + 
                                   fit->GetParameter(3))*sigmalizedYieldExtractionRange;
            upIntegrationRange = fit->GetParameter(1) + 
                                 (fit->GetParameter(2) + 
                                  fit->GetParameter(3))*sigmalizedYieldExtractionRange;

            if (isBGFixedForThisPTAltFit)
            {
               lowIntegrationRangeAltFitAB = altFitAB->GetParameter(1) - 
                                             (altFitAB->GetParameter(2) + 
                                              altFitAB->GetParameter(3))*
                                             sigmalizedYieldExtractionRange;
               upIntegrationRangeAltFitAB = altFitAB->GetParameter(1) + 
                                            (altFitAB->GetParameter(2) + 
                                             altFitAB->GetParameter(3))*
                                            sigmalizedYieldExtractionRange;
               lowIntegrationRangeAltFitFreeG = altFitFreeG->GetParameter(1) - 
                                                (altFitFreeG->GetParameter(2) + 
                                                 altFitFreeG->GetParameter(3))*
                                                sigmalizedYieldExtractionRange;
               upIntegrationRangeAltFitFreeG = altFitFreeG->GetParameter(1) + 
                                               (altFitFreeG->GetParameter(2) + 
                                                altFitFreeG->GetParameter(3))*
                                               sigmalizedYieldExtractionRange;
               lowIntegrationRangeAltFitFixedG = altFitFixedG->GetParameter(1) - 
                                                 (altFitFixedG->GetParameter(2) + 
                                                  altFitFixedG->GetParameter(3))*
                                                 sigmalizedYieldExtractionRange;
               upIntegrationRangeAltFitFixedG = altFitFixedG->GetParameter(1) + 
                                                (altFitFixedG->GetParameter(2) + 
                                                 altFitFixedG->GetParameter(3))*
                                                sigmalizedYieldExtractionRange;
            }

            double rawYield = GetYield(distrMInv, fitBG, lowIntegrationRange, upIntegrationRange);
            double rawYieldSysErr = 1e-15;

            if (isBGFixedForThisPTAltFit)
            {
               const double rawYieldAltFitAB = GetYield(distrMInv, altFitBGAB, 
                                                        lowIntegrationRangeAltFitAB, 
                                                        upIntegrationRangeAltFitAB);
               const double rawYieldAltFitFreeG = GetYield(distrMInv, altFitBGFreeG, 
                                                           lowIntegrationRangeAltFitFreeG, 
                                                           upIntegrationRangeAltFitFreeG);
               const double rawYieldAltFitFixedG = GetYield(distrMInv, altFitBGFixedG, 
                                                            lowIntegrationRangeAltFitFixedG, 
                                                            upIntegrationRangeAltFitFixedG);
               rawYieldSysErr = CppTools::RMS(rawYield - rawYieldAltFitAB,
                                              rawYield - rawYieldAltFitFreeG,
                                              rawYield - rawYieldAltFitFixedG);
            }

            double rawYieldStatErr = 
               sqrt(distrMInvFG->Integral(distrMInvFG->GetXaxis()->FindBin(lowIntegrationRange),
                                          distrMInvFG->GetXaxis()->FindBin(upIntegrationRange)));

            // 2*pi*pT*dpT*N_{evt}
            const double rawYieldNorm = 2.*M_PI*(pTBinRanges[i] + pTBinRanges[i + 1])/2.*
                                        (pTBinRanges[i + 1] - pTBinRanges[i])*numberOfEvents;
            rawYield /= rawYieldNorm;
            rawYieldStatErr /= rawYieldNorm;
            rawYieldSysErr /= rawYieldNorm;

            distrRawYieldVsPTStatErr.SetBinContent(i + 1, rawYield);
            distrRawYieldVsPTStatErr.SetBinError(i + 1, rawYieldStatErr);

            distrRawYieldVsPTSysErr.SetBinContent(i + 1, rawYield);
            distrRawYieldVsPTSysErr.SetBinError(i + 1, rawYieldSysErr);

            fit->SetLineWidth(4);
            fitBG->SetLineWidth(4);

            fit->SetLineColorAlpha(kRed - 3, 0.8);
            fitBG->SetLineColorAlpha(kGray + 2, 0.8);

            fitBG->SetLineStyle(7);

            if (isBGFixedForThisPTAltFit)
            {
               altFitAB->SetLineWidth(3);
               altFitBGAB->SetLineWidth(3);
               altFitFreeG->SetLineWidth(3);
               altFitBGFreeG->SetLineWidth(3);
               altFitFixedG->SetLineWidth(3);
               altFitBGFixedG->SetLineWidth(3);

               altFitAB->SetLineColorAlpha(kP6Blue, 0.8);
               altFitBGAB->SetLineColorAlpha(kP6Blue, 0.6);
               altFitFreeG->SetLineColorAlpha(kP6Red, 0.8);
               altFitBGFreeG->SetLineColorAlpha(kP6Red, 0.6);
               altFitFixedG->SetLineColorAlpha(kP6Gray, 0.8);
               altFitBGFixedG->SetLineColorAlpha(kP6Gray, 0.6);

               altFitBGAB->SetLineStyle(7);
               altFitBGFreeG->SetLineStyle(7);
               altFitBGFixedG->SetLineStyle(7);
            }
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
         distrMInvFGLR->SetLineWidth(3);
         distrMInvBGLR->SetLineWidth(2);

         distrMInv->SetLineColor(kGray + 3);

         int lastNonZeroBinLR = 1;
         for (int i = distrMInvFGLR->GetXaxis()->GetNbins(); i >= 1; i--)
         {
            if (distrMInvFGLR->GetBinContent(i) > 1e-7 || distrMInvBGLR->GetBinContent(i) > 1e-7) 
            {
               lastNonZeroBinLR = i;
               break;
            }
         }

         distrMInvFGLR->GetXaxis()->SetRange(1, lastNonZeroBinLR + 2);
         distrMInvBGLR->GetXaxis()->SetRange(1, lastNonZeroBinLR + 2);

         { /* canvas with invariant mass distribution with subtracted background only */
            TCanvas canvMInv("canv MInv", "", 800, 800);

            gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05); 
            gPad->SetLeftMargin(0.173); gPad->SetBottomMargin(0.112);

            ROOTTools::DrawFrame(distrMInv, "", "#it{M}_{inv} [GeV/#it{c}^{2}]", "Counts", 
                                 1., 1.9, 0.05, 0.05, true, false);

            text.DrawTextNDC(0.9, 0.93, (methodName).c_str());
            texText.DrawLatexNDC(0.2, 0.88, (CppTools::DtoStr(pTBinRanges[i], 1) + 
                                 " < #it{p}_{T} < " + 
                                 CppTools::DtoStr(pTBinRanges[i + 1], 1)).c_str());
            text.DrawTextNDC(0.85, 0.93, centrality["name_tex"].as<std::string>().c_str());
            if (performFit)
            {
               /*
               texText.DrawLatexNDC(0.2, 0.81, 
                                    ("#it{#chi}^{2}/NDF=" + 
                                     CppTools::DtoStr(fit->GetChisquare()/
                                                      fit->GetNDF(), 1)).c_str());
                                     */
               fitBG->Draw("SAME");
               fit->Draw("SAME");
            }

            distrMInv->Draw("SAME");

            ROOTTools::PrintCanvas(&canvMInv, outputDir + "/" + resonanceName + "_" + 
                                   centrality["name"].as<std::string>() + "_" +
                                   CppTools::DtoStr(pTBinRanges[i], 1) + "-" + 
                                   CppTools::DtoStr(pTBinRanges[i + 1], 1));
         } /* canvas with invariant mass distribution with subtracted background only */


         { /* summary canvas: MInv, FGBG, FG/BG */
            TCanvas canvMInvSummary("canv MInv summary", "", 1920, 1080);

            distrMInv->SetFillStyle(3003);
            distrMInvFG->SetFillStyle(3005);
            distrMInvBG->SetFillStyle(3004);
            distrMInvFGLR->SetFillStyle(3005);
            distrMInvBGLR->SetFillStyle(3004);

            canvMInvSummary.Divide(3, 2);
            
            canvMInvSummary.cd(1);

            gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05); 
            gPad->SetLeftMargin(0.173); gPad->SetBottomMargin(0.112);

            ROOTTools::DrawFrame(distrMInv, "", "#it{M}_{inv} [GeV/#it{c}^{2}]", "Counts", 
                                 1., 1.7, 0.05, 0.05, true, false);

            texText.DrawLatexNDC(0.2, 0.85, (CppTools::DtoStr(pTBinRanges[i], 1) + 
                                 " < #it{p}_{T} < " + 
                                 CppTools::DtoStr(pTBinRanges[i + 1], 1)).c_str());
            text.DrawTextNDC(0.8, 0.92, centrality["name_tex"].as<std::string>().c_str());
            text.DrawTextNDC(0.88, 0.92, (methodName).c_str());

            if (performFit)
            {
               /*
               texText.DrawLatexNDC(0.2, 0.74, 
                                    ("#it{#chi}^{2}/NDF = " + 
                                     CppTools::DtoStr(fit->GetChisquare()/
                                                      fit->GetNDF(), 1)).c_str());
                                     */
               /*

               texText.DrawLatexNDC(0.6, 0.9, ("#it{#mu}=" + CppTools::DtoStr(fit->GetParameter(1)*
                                                                              1000., 1) + 
                                               " [MeV/c^{2}]").c_str());
               texText.DrawLatexNDC(0.6, 0.85, ("#it{#Gamma}=" + 
                                                CppTools::DtoStr(fit->GetParameter(2)*1000., 1) + 
                                                " [MeV/c^{2}]").c_str());
                                                */

               fitBG->Draw("SAME");
               fit->Draw("SAME");

               TLine lowYieldExtrRangeLine(lowIntegrationRange, 
                                           distrMInv->GetBinContent(distrMInv->GetMinimumBin()), 
                                           lowIntegrationRange, 
                                           distrMInv->GetBinContent(distrMInv->GetMaximumBin()));

               TLine upYieldExtrRangeLine(upIntegrationRange, 
                                          distrMInv->GetBinContent(distrMInv->GetMinimumBin()), 
                                          upIntegrationRange, 
                                          distrMInv->GetBinContent(distrMInv->GetMaximumBin()));

               lowYieldExtrRangeLine.SetLineColorAlpha(kGray + 1, 0.5);
               lowYieldExtrRangeLine.SetLineStyle(3);
               lowYieldExtrRangeLine.SetLineWidth(3);

               upYieldExtrRangeLine.SetLineColorAlpha(kGray + 1, 0.5);
               upYieldExtrRangeLine.SetLineStyle(3);
               upYieldExtrRangeLine.SetLineWidth(3);

               lowYieldExtrRangeLine.Clone()->Draw();
               upYieldExtrRangeLine.Clone()->Draw();
            }

            distrMInv->Clone()->Draw("SAME");

            canvMInvSummary.cd(4);

            gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05); 
            gPad->SetLeftMargin(0.173); gPad->SetBottomMargin(0.112);

            ROOTTools::DrawFrame(static_cast<TH1D *>(distrMInv->Clone()), 
                                 "", "#it{M}_{inv} [GeV/#it{c}^{2}]", "Counts", 1., 1.7);

            if (performFit && performAltFits)
            {
               altFitBGAB->Draw("SAME");
               altFitBGFreeG->Draw("SAME");
               altFitBGFixedG->Draw("SAME");

               altFitAB->Draw("SAME");
               altFitFreeG->Draw("SAME");
               altFitFixedG->Draw("SAME");

               TLegend legend(0.75, 0.6, 0.95, 0.95);

               legend.SetLineColorAlpha(0, 0.);
               legend.SetFillColorAlpha(0, 0.);

               legend.AddEntry(altFitAB, "Alt BG", "L");
               legend.AddEntry(altFitFreeG, "Free #Gamma", "L");
               legend.AddEntry(altFitFixedG, "Fixed #Gamma", "L");

               legend.Clone()->Draw();
            }

            canvMInvSummary.cd(5);

            gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05); 
            gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.112);

            if (distrMInvBG->GetEntries() > 1e-3 && distrMInvFG->GetEntries() > 1e-3)
            {
               TH1D *ratioFGBG = static_cast<TH1D *>(distrMInvFG->Clone());
               ratioFGBG->Divide(distrMInvBG);

               ratioFGBG->SetLineWidth(2);
               ratioFGBG->SetLineColor(kBlack);

               // draws frame with histogram ratioFGBG points
               ROOTTools::DrawFrame(ratioFGBG, "", "#it{M}_{inv} [GeV/#it{c}^{2}]", 
                                    "FG/BG", 1., 1.7, 0.05, 0.05, true, false);

               if (performFit)
               {
                  TLine lowYieldExtrRangeLine(lowIntegrationRange, 
                                              ratioFGBG->GetBinContent(ratioFGBG->GetMinimumBin()), 
                                              lowIntegrationRange, 
                                              ratioFGBG->GetBinContent(ratioFGBG->GetMaximumBin()));

                  TLine upYieldExtrRangeLine(upIntegrationRange, 
                                             ratioFGBG->GetBinContent(ratioFGBG->GetMinimumBin()), 
                                             upIntegrationRange, 
                                             ratioFGBG->GetBinContent(ratioFGBG->GetMaximumBin()));

                  lowYieldExtrRangeLine.SetLineColorAlpha(kGray + 1, 0.5);
                  lowYieldExtrRangeLine.SetLineStyle(3);
                  lowYieldExtrRangeLine.SetLineWidth(3);

                  upYieldExtrRangeLine.SetLineColorAlpha(kGray + 1, 0.5);
                  upYieldExtrRangeLine.SetLineStyle(3);
                  upYieldExtrRangeLine.SetLineWidth(3);

                  lowYieldExtrRangeLine.Clone()->Draw();
                  upYieldExtrRangeLine.Clone()->Draw();

                  // graph that contains fit to correlated BG ratio
                  TGraph fitToBGRatio;
                  // graph that contains background fit to correlated BG ratio
                  TGraph fitBGToBGRatio;

                  for (int i = CppTools::Maximum(ratioFGBG->GetXaxis()->FindBin(fitRangeMin), 1);
                       i <= CppTools::Minimum(ratioFGBG->GetXaxis()->FindBin(fitRangeMax), 
                                              ratioFGBG->GetXaxis()->GetNbins()); i++)
                  {
                     const double xVal = ratioFGBG->GetXaxis()->GetBinCenter(i);
                     fitToBGRatio.
                        AddPoint(xVal, (fit->Eval(xVal) + distrMInvBG->GetBinContent(i))/
                                       distrMInvBG->GetBinContent(i));
                     fitBGToBGRatio.
                        AddPoint(xVal, (fitBG->Eval(xVal) + distrMInvBG->GetBinContent(i))/
                                       distrMInvBG->GetBinContent(i));
                  }

                  fitBGToBGRatio.SetLineStyle(2);

                  fitToBGRatio.SetLineColorAlpha(kRed - 3, 0.8);
                  fitBGToBGRatio.SetLineColorAlpha(kGray + 2, 0.8);

                  fitToBGRatio.SetLineWidth(4);
                  fitBGToBGRatio.SetLineWidth(4);

                  fitToBGRatio.Clone()->Draw("L");
                  fitBGToBGRatio.Clone()->Draw("L");
               }

               ratioFGBG->Draw("SAME");
            }
            else
            {
               text.DrawTextNDC(0.88, 0.92, "No data");
            }

            canvMInvSummary.cd(6);

            gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05); 
            gPad->SetLeftMargin(0.148); gPad->SetBottomMargin(0.112);

            if (distrMInvBGLR->GetEntries() > 1e-3 && distrMInvFGLR->GetEntries() > 1e-3)
            {
               TH1D *ratioFGBGLR = static_cast<TH1D *>(distrMInvFGLR->Clone());
               ratioFGBGLR->Divide(distrMInvBGLR);

               ratioFGBGLR->SetLineWidth(2);
               ratioFGBGLR->SetLineColor(kBlack);

               ROOTTools::DrawFrame(ratioFGBGLR, "", "#it{M}_{inv} [GeV/#it{c}^{2}]", 
                                    "FG/BG", 1., 1.7);
            }
            else
            {
               text.DrawTextNDC(0.88, 0.92, "No data");
            }

            distrMInvFG->SetMaximum(distrMInvFG->GetMaximum()*1.15);
            distrMInvFGLR->SetMaximum(distrMInvFGLR->GetMaximum()*3.);

            canvMInvSummary.cd(2);

            gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05); 
            gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.112);

            distrMInv->Sumw2(false);
            distrMInvBG->Sumw2(false);
            distrMInvFG->Sumw2(false);

            if (distrMInvFG->GetEntries() > distrMInv->GetEntries())
            {
               distrMInvFG->SetMinimum(0.);
               ROOTTools::DrawFrame(distrMInvFG, "", "#it{M}_{inv} [GeV/#it{c}^{2}]", 
                                    "Counts", 1., 1.7, 0.05, 0.05, true, false);
            }
            else
            {
               ROOTTools::DrawFrame(distrMInv, "", "#it{M}_{inv} [GeV/#it{c}^{2}]", 
                                    "Counts", 1., 1.7, 0.05, 0.05, true, false);
            }

            if (distrMInvFG->GetEntries() > 1e-3) 
            {
               distrMInvFG->SetLineColorAlpha(kAzure + 2, 0.8);
               distrMInvFG->Draw("SAME PFC");
            }
            else text.DrawTextNDC(0.88, 0.9, "No data on foreground");

            if (distrMInvBG->GetEntries() > 1e-3) 
            {
               distrMInvBG->SetLineColorAlpha(kGreen + 2, 0.8);
               distrMInvBG->Draw("SAME PFC");
            }
            else text.DrawTextNDC(0.78, 0.9, "No data on background");

            distrMInv->SetLineColorAlpha(kRed + 2, 0.8);
            distrMInv->Draw("SAME PFC");

            canvMInvSummary.cd(3);

            gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05); 
            gPad->SetLeftMargin(0.148); gPad->SetBottomMargin(0.112);

            gPad->SetLogy();

            distrMInvBGLR->Sumw2(false);

            if (distrMInvFGLR->GetEntries() > 1e-3)
            {
               ROOTTools::DrawFrame(distrMInvFGLR, "", "#it{M}_{inv} [GeV/#it{c}^{2}]", 
                                    "Counts", 1., 1.7, 0.05, 0.05, true, false);
               distrMInvFGLR->SetLineColorAlpha(kAzure + 2, 0.8);
               distrMInvFGLR->Draw("SAME PFC");

               if (distrMInvBGLR->GetEntries() > 1e-3) 
               {
                  distrMInvBGLR->SetLineColorAlpha(kRed + 2, 0.8);
                  distrMInvBGLR->Draw("SAME PFC");
               }
               else text.DrawTextNDC(0.78, 0.9, "No data on background");
            }
            else text.DrawTextNDC(0.88, 0.9, "No data on foreground");

            ROOTTools::PrintCanvas(&canvMInvSummary, outputDir + "/Summary_" + resonanceName + "_" + 
                                   centrality["name"].as<std::string>() + "_" +
                                   CppTools::DtoStr(pTBinRanges[i], 1) + "-" + 
                                   CppTools::DtoStr(pTBinRanges[i + 1], 1), true, false);
         } /* summary canvas: MInv, FGBG, FG/BG */

         { /* FG, BG, and signal on the same canvas */
            TCanvas canvMInvFGBG("canv MInv FGBG", "", 800, 800);

            gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05); 
            gPad->SetLeftMargin(0.173); gPad->SetBottomMargin(0.112);

            distrMInv->Sumw2(false);
            distrMInvBG->Sumw2(false);
            distrMInvFG->Sumw2(false);

            if (distrMInvFG->GetEntries() > distrMInv->GetEntries())
            {
               distrMInvFG->SetMinimum(0.);
               ROOTTools::DrawFrame(distrMInvFG, "", "#it{M}_{inv} [GeV/#it{c}^{2}]", 
                                    "Counts", 1., 1.7, 0.05, 0.05, true, false);
            }
            else
            {
               ROOTTools::DrawFrame(distrMInv, "", "#it{M}_{inv} [GeV/#it{c}^{2}]", 
                                    "Counts", 1., 1.7, 0.05, 0.05, true, false);
            }

            if (distrMInvFG->GetEntries() > 1e-3) 
            {
               distrMInvFG->SetLineColorAlpha(kAzure + 2, 0.8);
               distrMInvFG->Draw("SAME PFC");
            }
            else text.DrawTextNDC(0.88, 0.9, "No data on foreground");

            if (distrMInvBG->GetEntries() > 1e-3) 
            {
               distrMInvBG->SetLineColorAlpha(kGreen + 2, 0.8);
               distrMInvBG->Draw("SAME PFC");
            }
            else text.DrawTextNDC(0.78, 0.9, "No data on background");

            distrMInv->SetLineColorAlpha(kRed + 2, 0.8);
            distrMInv->Draw("SAME PFC");

            TLegend legend(0.3, 0.86, 0.98, 0.94);
            legend.SetLineColorAlpha(0, 0.);
            legend.SetFillColorAlpha(0, 0.);
            legend.SetNColumns(3);

            legend.AddEntry(distrMInvFG, "FG", "PFC");
            legend.AddEntry(distrMInvBG, "BG", "PFC");
            legend.AddEntry(distrMInv, "FG-BG", "PFC");

            legend.Draw();

            ROOTTools::PrintCanvas(&canvMInvFGBG, outputDir + "/FGBG_" + resonanceName + "_" + 
                                   centrality["name"].as<std::string>() + "_" +
                                   CppTools::DtoStr(pTBinRanges[i], 1) + "-" + 
                                   CppTools::DtoStr(pTBinRanges[i + 1], 1), false);
         } /* FG, BG, and signal on the same canvas */

         numberOfCalls++;
      }

      distrMeansVsPT.SetLineColor(kRed - 2);
      distrMeansVsPT.SetMarkerColor(kRed - 2);
      distrMeansVsPT.SetLineWidth(4);

      distrGammasVsPT.SetLineColor(kRed - 2);
      distrGammasVsPT.SetMarkerColor(kRed - 2);
      distrGammasVsPT.SetLineWidth(4);

      distrRawYieldVsPTStatErr.SetLineColor(kRed - 2);
      distrRawYieldVsPTStatErr.SetMarkerColor(kRed - 2);
      distrRawYieldVsPTStatErr.SetLineWidth(4);

      distrRawYieldVsPTSysErr.SetLineColor(kRed - 2);
      distrRawYieldVsPTSysErr.SetFillStyle(1001);
      distrRawYieldVsPTSysErr.SetFillColorAlpha(kRed - 2, 0.5);

      TLine massResonancePDG(pTBinRanges[0], massResonance, pTBinRanges[pTNBins], massResonance);
      massResonancePDG.SetLineColorAlpha(kBlack, 0.5);
      massResonancePDG.SetLineStyle(2);
      massResonancePDG.SetLineWidth(4);

      TLine gammaResonancePDG(pTBinRanges[0], gammaResonance, pTBinRanges[pTNBins], gammaResonance);
      gammaResonancePDG.SetLineColorAlpha(kBlack, 0.5);
      gammaResonancePDG.SetLineStyle(2);
      gammaResonancePDG.SetLineWidth(4);

      distrMeansVsPT.SetMaximum(massResonance*1.05);
      distrMeansVsPT.SetMinimum(massResonance*0.95);

      distrGammasVsPT.SetMaximum(gammaResonance*1.5);
      distrGammasVsPT.SetMinimum(gammaResonance/2.);

      TCanvas canvMeansVsPT("canv means vs pT", "", 800, 800);

      gPad->SetRightMargin(0.03);
      gPad->SetTopMargin(0.02);
      gPad->SetLeftMargin(0.172);
      gPad->SetBottomMargin(0.112);

      ROOTTools::DrawFrame(&distrMeansVsPT, "", "#it{p}_{T} [GeV/#it{c}]", 
                           "#it{#mu} [GeV/#it{c}^{2}]", 1., 1.82);

      text.DrawTextNDC(0.9, 0.95, (methodName).c_str());
      massResonancePDG.Draw();

      ROOTTools::PrintCanvas(&canvMeansVsPT, outputDir + "/" + resonanceName + 
                             "_means_" + centralityName);

      TCanvas canvGammasVsPT("canv gammas vs pT", "", 800, 800);

      gPad->SetRightMargin(0.03);
      gPad->SetTopMargin(0.02);
      gPad->SetLeftMargin(0.172);
      gPad->SetBottomMargin(0.112);

      ROOTTools::DrawFrame(&distrGammasVsPT, "", "#it{p}_{T} [GeV/#it{c}]", 
                           "#it{#Gamma} [GeV/#it{c}^{2}]", 1., 1.82);

      text.DrawTextNDC(0.9, 0.95, (methodName).c_str());
      gammaResonancePDG.Draw();

      ROOTTools::PrintCanvas(&canvGammasVsPT, outputDir + "/" + resonanceName + 
                             "_gammas_" + centralityName);

      TCanvas canvRawYieldVsPT("canv raw yield vs pT", "", 800, 800);

      gPad->SetLogy();

      gPad->SetRightMargin(0.03);
      gPad->SetTopMargin(0.02);
      gPad->SetLeftMargin(0.141);
      gPad->SetBottomMargin(0.112);

      ROOTTools::DrawFrame(&distrRawYieldVsPTStatErr, "", "#it{p}_{T} [GeV/#it{c}]", 
                           "#it{dY}_{raw}/#it{dp}_{T} [(GeV/#it{c})^{-1}]", 1., 1.35);
      distrRawYieldVsPTSysErr.Draw("SAME E2");

      text.DrawTextNDC(0.9, 0.95, (methodName).c_str());

      ROOTTools::PrintCanvas(&canvRawYieldVsPT, outputDir + "/" + resonanceName + 
                             "_raw_yield_" + centralityName);

      parametersOutputFile->mkdir(centralityName.c_str());
      parametersOutputFile->cd(centralityName.c_str());

      distrMeansVsPT.Write();
      distrGammasVsPT.Write();
      distrRawYieldVsPTStatErr.Write();
      distrRawYieldVsPTSysErr.Write();
   }

   parametersOutputFile->Close();
}

void AnalyzeRealMInv::SetGaussianBroadeningFunction()
{
   const std::string inputFileName = "data/Parameters/GaussianBroadening/" + 
                                     runName + "/" + resonanceName + ".root";

   if (!std::filesystem::exists(inputFileName))
   {
      pBar.Clear();
      CppTools::PrintWarning(inputFileName + " does not exists. Using rough estimate sigma = 5 MeV "\
                             "Run executable bin/EstimateGassianBroadening first to better "\
                             "estimate the gaussian broadening parameter");
      pBar.RePrint();

      gaussianBroadeningEstimatorFunc = new TF1("gaussian broadening rough estimate", "5e-3");
   }
   else
   {
      gaussianBroadeningEstimatorFunc = 
         static_cast<TF1 *>(TFile::Open(inputFileName.c_str())->Get("gaussian broadening sigma fit"));
   }
}

double AnalyzeRealMInv::GetYield(TH1D *distrMInv, TF1 *funcBG, const double xMin, const double xMax)
{
   // integral over the signal
   double integral = 0.;
   // integral over the background
   double integralBG = 0.;

   if (distrMInv->GetXaxis()->FindBin(xMin) < 1 || 
       distrMInv->GetXaxis()->FindBin(xMax) > distrMInv->GetXaxis()->GetNbins())
   {
      pBar.Clear();
      CppTools::PrintWarning("AnalyzeRealMInv::GetYield: specified integration range is "\
                             "outside the histogram range; ignoring underflow and overflow bins");
      pBar.RePrint();
   }

   for (int i = CppTools::Maximum(distrMInv->GetXaxis()->FindBin(xMin), 1); 
        i <= CppTools::Minimum(distrMInv->GetXaxis()->FindBin(xMax), 
                               distrMInv->GetXaxis()->GetNbins()); i++)
   {
      integral += distrMInv->GetBinContent(i);

      // integrating background over the single bin for better estimation
      for (double m = distrMInv->GetXaxis()->GetBinLowEdge(i); 
           m <= distrMInv->GetXaxis()->GetBinUpEdge(i); 
           m += distrMInv->GetXaxis()->GetBinWidth(i)/100.)
      {
         integralBG += funcBG->Eval(m);
      }
   }

   // normalizing background integral by the number of integration steps
   integral -= integralBG/101.;

   return integral;
}

TFile *AnalyzeRealMInv::SetFixedBGFile(const std::string& inputFileName, 
                                       const std::string& fitTypeName,
                                       const bool printFreeFitWarning)
{
   if (std::filesystem::exists(inputFileName))
   {
      pBar.Clear();
      CppTools::PrintInfo("Fixed BG fits file for " + fitTypeName + " fits was found");
      pBar.RePrint();
      return TFile::Open(inputFileName.c_str());
   }
   else if (printFreeFitWarning)
   {
      pBar.Clear();
      CppTools::PrintWarning("Fixed BG fits file for " + fitTypeName + " fits was not found; "\
                             "using free BG fit for these fits");
      pBar.RePrint();
      return nullptr;
   }
   pBar.Clear();
   CppTools::PrintWarning("Fixed BG fits file for " + fitTypeName + " fits was not found; "\
                          "no fits will be performed for this fit type");
   pBar.RePrint();
   return nullptr;
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

#endif /* ANALYZE_REAL_M_INV_CPP */
