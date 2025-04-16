/** 
 *  @file  EMCalDM.cpp
 *  @brief Contains function that can be used to obtain deadmaps of EMCal. More accurate deadmaps should also include fiducial cuts with the use of GUI/DM.cpp
 *
 * Usage: from the root of the repository run
 *
 * @code
 * root -q cling/EMCalDM.cpp
 * @endcode
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysisPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#pragma once

#include "IOTools.hpp"

void AnalyzeSector(const std::string& sectorName);

std::unique_ptr<TFile> inputFile;
std::string runName;

void EMCalDM()
{
   CppTools::PrintInfo("This is a program that help you obtain the EMCal deadmaps by "\
                       "excluding EMCal towers that have too small or too little signal "\
                       "compared to the average (more than 3 sigmas)");

   CppTools::PrintInfo("Available runs:");
   system("ls data/PostSim");

   CppTools::Print("Type in the run name you want to process:");
   std::cout << ">> ";
   
   std::cin >> runName;

   CppTools::Print();

   const std::string inputFileName = "data/Real/" + runName + "/SingleTrack/sum.root";
   CppTools::CheckInputFile(inputFileName);
   inputFile = std::unique_ptr<TFile>(TFile::Open(inputFileName.c_str(), "READ"));

   for (int i = 0; i < 4; i++)
   {
      AnalyzeSector("EMCale" + std::to_string(i));
      AnalyzeSector("EMCalw" + std::to_string(i));
   }
}

void AnalyzeSector(const std::string& sectorName)
{
   const std::string parametersOutputFileName = "data/Parameters/Deadmaps/" + 
                                                runName + "/" + sectorName + ".txt";

   if (CppTools::FileExists(parametersOutputFileName))
   {
      CppTools::PrintInfo("File " + parametersOutputFileName + " already exists. "\
                          "Skipping sector " + sectorName);
      return;
   }

   TH2F *heatmap = static_cast<TH2F *>(inputFile->Get(("Heatmap: " + sectorName).c_str()));

   TH1F distrTowerSignals(sectorName.c_str(), "", 100, 1., -1.);

   const double towersSignalIntegral = heatmap->Integral(1, heatmap->GetXaxis()->GetNbins(), 
                                                         1, heatmap->GetYaxis()->GetNbins());

   // Obtaining the distribution of the signals in every tower
   for (int i = 1; i <= heatmap->GetXaxis()->GetNbins(); i++)
   {
      for (int j = 1; j <= heatmap->GetYaxis()->GetNbins(); j++)
      {
         if (heatmap->GetBinContent(i, j) < 1e-15) continue;

         distrTowerSignals.Fill(heatmap->GetBinContent(i, j));
      }
   }

   TF1 towerSignalsFit("tower signal fit", "gaus(0)");
   towerSignalsFit.
      SetParameters(distrTowerSignals.GetMaximum(), 
                    distrTowerSignals.GetBinCenter(distrTowerSignals.GetXaxis()->GetNbins()/2),
                    distrTowerSignals.GetBinCenter(distrTowerSignals.GetXaxis()->GetNbins()/2));

   // approximating the tower signals distribution to get the sigma of the said distribution
   distrTowerSignals.Fit(&towerSignalsFit, "QMBN");

   // file for writing the deadmap
   std::ofstream parametersOutputFile(parametersOutputFileName);

   parametersOutputFile << heatmap->GetXaxis()->GetNbins() << " ";
   parametersOutputFile << heatmap->GetXaxis()->GetBinLowEdge(1) << " ";
   parametersOutputFile << heatmap->GetXaxis()->GetBinUpEdge(heatmap->GetXaxis()->GetNbins()) << " ";
   parametersOutputFile << heatmap->GetYaxis()->GetNbins() << " ";
   parametersOutputFile << heatmap->GetYaxis()->GetBinLowEdge(1) << " ";
   parametersOutputFile << heatmap->GetYaxis()->GetBinUpEdge(heatmap->GetYaxis()->GetNbins());
   parametersOutputFile << std::endl;

   TH2F *cutHeatmap = static_cast<TH2F *>(heatmap->Clone());

   const double mean = towerSignalsFit.GetParameter(1);
   const double sigma = towerSignalsFit.GetParameter(2);

   // cutting towers that are beyond 3 sigmas from the mean of the towers distribution
   for (int i = 1; i <= heatmap->GetYaxis()->GetNbins(); i++)
   {
      for (int j = 1; j <= heatmap->GetXaxis()->GetNbins(); j++)
      {
         if (heatmap->GetBinContent(j, i) < 1e-15 ||
             heatmap->GetBinContent(j, i) < mean - 3.*sigma ||
             heatmap->GetBinContent(j, i) > mean + 3.*sigma)
         {
            cutHeatmap->SetBinContent(j, i, 0.);
            parametersOutputFile << 1;
         }
         else
         {
            parametersOutputFile << 0;
         }
         if (j < heatmap->GetXaxis()->GetNbins())
         {
            parametersOutputFile << " ";
         }
         else if (i < heatmap->GetYaxis()->GetNbins())
         {
            parametersOutputFile << std::endl;
         }
      }
   }
   parametersOutputFile.close();
   CppTools::PrintInfo("Deadmap of " + sectorName + " were written in " + parametersOutputFileName);
}
