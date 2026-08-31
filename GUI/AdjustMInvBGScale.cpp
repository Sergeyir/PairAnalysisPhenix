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
#include "Box.hpp"

#include "TCanvasTools.hpp"

#include "GUIFit.hpp"

#include "AnalyzeRealMInv.hpp"

void PrintHelp();
void Draw();
void Exec();

/// Contains invariant mass foreground distributions in different pT ranges in one centrality class
std::vector<TH1D *> distrsMInvFG;
/// Contains invariant mass background distributions in different pT ranges in one centrality class
std::vector<TH1D *> distrsMInvBG;
/// Contains low resolution invariant mass foreground distributions in different pT ranges in one centrality class
std::vector<TH1D *> distrsMInvFGLR;
/// Contains low resolution invariant mass background distributions in different pT ranges in one centrality class
std::vector<TH1D *> distrsMInvBGLR;
/// Contains names of invariant mass distributions to be displayed on canvas
std::vector<std::string> distrMInvNames;
/// Contains MInv BG rescale values
std::vector<double> rescaleMInvBG;

std::string outputFileName;

TCanvas *canv;

bool parametersSet = false;

// shows the current pT bin index of an invariant mass distribution to display
int currentPTBinIndex = 0;

YAML::Node method;

using namespace AnalyzeRealMInv;

void AdjustMInvBGScale()
{
	gStyle->SetOptStat(0);
   gErrorIgnoreLevel = kWarning;

   TH1::AddDirectory(kFALSE);
   TH2::AddDirectory(kFALSE);

   ROOT::EnableImplicitMT(4);
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

   resonanceName = inputYAMLResonance["name"].as<std::string>();

   daughter1Id = inputYAMLResonance["daughter1_id"].as<int>();
   daughter2Id = inputYAMLResonance["daughter2_id"].as<int>();

   minMInv = inputYAMLResonance["m_inv_range_min"].as<double>();
   maxMInv = inputYAMLResonance["m_inv_range_max"].as<double>();

   inputFile = TFile::Open(inputFileName.c_str(), "READ");

   pTNBins = inputYAMLResonance["pt_bins"].size();

   for (unsigned int i = 0; i < pTNBins; i++)
   {
      pTBinRanges.push_back(inputYAMLResonance["pt_bins"][i]["min"].as<double>());
   }
   pTBinRanges.push_back(inputYAMLResonance["pt_bins"][pTNBins - 1]["max"].as<double>());

   method = inputYAMLResonance["pair_selection_methods"][methodBinIndex];

   const unsigned int pTBinFitMin = method["centrality_bin_parameters"][centralityBinIndex]
                                          ["pt_bin_min"].as<int>();
   const unsigned int pTBinFitMax = method["centrality_bin_parameters"][centralityBinIndex]
                                          ["pt_bin_max"].as<int>();

   numberOfIterations = pTBinFitMax - pTBinFitMin + 1;

   parametersOutputDir = "data/Parameters/MInvBGRescale/" + runName;
   std::filesystem::create_directories(parametersOutputDir);

   const std::string methodName = 
      inputYAMLResonance["pair_selection_methods"][methodBinIndex]["name"].as<std::string>();
   const std::string centralityName = 
      inputYAMLResonance["centrality_bins"][centralityBinIndex]["name"].as<std::string>();

   outputFileName = 
      parametersOutputDir + "/" + resonanceName + "_" + methodName + "_" + centralityName + ".txt";

   if (std::filesystem::exists(outputFileName))
   {
      CppTools::PrintWarning("File " + outputFileName + " already exists; the old file "\
                             "will be renamed to " + outputFileName + ".backup");

      std::ifstream inputFile(outputFileName);

      double val;
      while (inputFile >> val)
      {
         rescaleMInvBG.push_back(val);
      }

      parametersSet = true;
   }

   PerformMInvFits(method, centralityBinIndex);

   if (parametersSet)
   {
      if (rescaleMInvBG.size() != distrsMInvFG.size())
      {
         CppTools::PrintError("Number of rescale values mismatch the number of pT bins");
      }
      
      static_cast<void>(system(("cp " + outputFileName + " " + outputFileName + ".backup").c_str()));
   }

   gROOT->SetBatch(false);
 
   canv = new TCanvas("", "", 1600, 400);
   canv->Divide(4);
   Draw();
   PrintHelp();

   gPad->AddExec("exec", "Exec()");
}

void Draw()
{
   TH1D *distrMInv = static_cast<TH1D *>(distrsMInvFG[currentPTBinIndex]->Clone());
   TH1D *distrMInvFG = static_cast<TH1D *>(distrsMInvFG[currentPTBinIndex]->Clone());
   TH1D *distrMInvBG = static_cast<TH1D *>(distrsMInvBG[currentPTBinIndex]->Clone());
   TH1D *distrMInvFGLR = static_cast<TH1D *>(distrsMInvFGLR[currentPTBinIndex]->Clone());
   TH1D *distrMInvBGLR = static_cast<TH1D *>(distrsMInvBGLR[currentPTBinIndex]->Clone());

   distrMInv->SetLineColor(kRed - 3);
   distrMInv->SetLineWidth(2);

   distrMInvBG->Scale(rescaleMInvBG[currentPTBinIndex]);
   distrMInvBGLR->Scale(rescaleMInvBG[currentPTBinIndex]);

   distrMInv->Add(distrMInvBG, -1.);

   canv->cd(1);

   ROOTTools::DrawFrame(distrMInvFG->GetXaxis()->GetBinLowEdge(1), 0.,
                        distrMInvFG->GetXaxis()->GetBinUpEdge(distrMInvFG->GetXaxis()->GetNbins()),
                        distrMInvFG->GetMaximum()*1.1, distrMInvNames[currentPTBinIndex], "", "");

   distrMInvFG->Clone()->Draw("SAME");
   distrMInvBG->Clone()->Draw("SAME");
   distrMInv->Clone()->Draw("SAME");

   gPad->Modified();
   gPad->Update();

   distrMInvFG->Divide(distrMInvBG);
   distrMInvFGLR->Divide(distrMInvBGLR);

   distrMInv->SetLineColor(kBlack);
   distrMInvFG->SetLineColor(kBlack);

   canv->cd(2);

   distrMInv->Draw();

   gPad->Modified();
   gPad->Update();

   canv->cd(3);

   distrMInvFG->Draw();

   gPad->Modified();
   gPad->Update();

   canv->cd(4);

   distrMInvFGLR->Draw();

   gPad->Modified();
   gPad->Update();

   canv->Modified();
   canv->Update();
   CppTools::Print("Current pT bin scale:", rescaleMInvBG[currentPTBinIndex]);

}

void Write()
{
   std::ofstream outputFile(outputFileName);

   for (int i = 0; i < rescaleMInvBG.size(); i++)
   {
      outputFile << rescaleMInvBG[i];
      if (i < rescaleMInvBG.size() - 1) outputFile << std::endl;
   }
   CppTools::PrintInfo("Rescale values were written in " + outputFileName);

   outputFile.close();
}

void PrintHelp()
{
   CppTools::PrintInfo("Move your mouse towards the right histogram. You can use the "\
                       "following keys to direct the rescale process. After each key "\
                       "press root GUI will only detect the next key press after "\
                       "the mouse movement over the right histogram.");

   CppTools::Box box("Keyboard configuration");
   box.AddEntry("Next pT bin", "j");
   box.AddEntry("Previos pT bin", "k");
   box.AddEntry("Increase the BG scale by 0.001", "i");
   box.AddEntry("Decrease the BG scale by 0.001", "d");
   box.AddEntry("Print results", "p");
   box.AddEntry("Print help", "h");
   box.Print();
}

void Exec()
{
   const int px = gPad->GetEventX();
   const int py = gPad->GetEventY();

   switch (px)
   {
      case 'j':
      {
         if (distrsMInvFG.size() == 1)
         {
            CppTools::PrintWarning("Cannot switch between histograms "\
                                   "since only one was added");
         }
         else if (currentPTBinIndex < static_cast<int>(distrsMInvFG.size() - 1)) currentPTBinIndex++;
         else currentPTBinIndex = 0;

         Draw();

         break;
      }
      case 'k':
      {
         if (distrsMInvFG.size() == 1)
         {
         CppTools::PrintWarning("Cannot switch between histograms "\
                                "since only one was added");
         }
         else if (currentPTBinIndex > 0) currentPTBinIndex--;
         else currentPTBinIndex = distrsMInvFG.size() - 1;

         Draw();

         break;
      }
      case 'p':
      {
         Write();
         break;
      }
      case 'i':
      {
         rescaleMInvBG[currentPTBinIndex] += 0.005;
         Draw();
         break;
      }
      case 'd':
      {
         rescaleMInvBG[currentPTBinIndex] -= 0.005;
         Draw();
         break;
      }
      case 'h':
      {
         PrintHelp();
         break;
      }
   }
}

void AnalyzeRealMInv::PerformMInvFits(const YAML::Node& method, const unsigned int centralityBinIndex)
{
   const std::string methodName = method["name"].as<std::string>();

   const YAML::Node centrality = inputYAMLResonance["centrality_bins"][centralityBinIndex];

   const std::string centralityName = centrality["name"].as<std::string>();

   pBar.SetText("Preparing M_{inv}");

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

      distrMInvFG->SetLineColor(kAzure - 3);
      distrMInvBG->SetLineColor(kGreen - 3);

      distrMInvFGLR->SetLineColor(kBlack);
      distrMInvFGLR->SetLineWidth(2);

      distrMInvFG->SetLineWidth(2);
      distrMInvBG->SetLineWidth(2);
      distrMInv->SetLineWidth(2);

      distrsMInvFG.emplace_back(static_cast<TH1D *>(distrMInvFG->Clone()));
      distrsMInvBG.emplace_back(static_cast<TH1D *>(distrMInvBG->Clone()));
      distrsMInvFGLR.emplace_back(static_cast<TH1D *>(distrMInvFGLR->Clone()));
      distrsMInvBGLR.emplace_back(static_cast<TH1D *>(distrMInvBGLR->Clone()));
      distrMInvNames.emplace_back((CppTools::DtoStr(pTBinRanges[i], 2) + "<p_{T}<" + 
                                   CppTools::DtoStr(pTBinRanges[i + 1], 2)));
      if (!parametersSet)
      {
         rescaleMInvBG.push_back(1.);
      }

      numberOfCalls++;
   }
   pBar.Finish();
}
