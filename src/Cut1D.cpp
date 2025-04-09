/** 
 *  @file   Cut1D.cpp 
 *  @brief  Contains realisations of functions to cut bad/dead areas from 1D distributions (slat, strip)
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/CalPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef CUT_1D_CPP
#define CUT_1D_CPP

#include "../include/Cut1D.hpp"

void Cut1D::CutDistribution(TH1D *distr, const std::string& name, const double minValueFromAverage)
{
   const int nBins = distr->GetXaxis()->GetNbins();
   /*
   if (nBins != static_cast<int>(distr->GetXaxis()->GetBinUpEdge(nBins) - 
                                 distr->GetXaxis()->GetBinLowEdge(1)))
   {
      CppTools::PrintError("Number of bins is inconsistent witht the number of " + name + "s");
   }
   */

   TH1D *cutDistr = static_cast<TH1D *>(distr->Clone());

   const double fullIntegral = distr->Integral(1, nBins);
   const double averageYValue = fullIntegral/static_cast<double>(nBins);

   system("mkdir -p par/Deadmaps/Run14HeAu200");
   std::ofstream cutOutputFile("par/Deadmaps/Run14HeAu200/" + name + ".txt");

   cutOutputFile << nBins << " " << distr->GetXaxis()->GetBinLowEdge(1) << " " << 
                    distr->GetXaxis()->GetBinUpEdge(nBins) << std::endl;

   for (int i = 1; i <= nBins; i++)
   {
      if (distr->GetBinContent(i) < averageYValue*minValueFromAverage)
      {
         cutOutputFile << 1;
         cutDistr->SetBinContent(i, 0);
      }
      else 
      {
         cutOutputFile << 0;
      }
      if (i < nBins) cutOutputFile << " ";
   }

   const double cutIntegral = cutDistr->Integral(1, nBins);

   CppTools::PrintInfo("Data lost for " + name + ": " + 
                       std::to_string((1. - cutIntegral/fullIntegral)*100.) + "%");


   distr->SetLineWidth(1);
   cutDistr->SetLineWidth(1);

   distr->SetLineColor(kRed - 3);
   cutDistr->SetLineColor(kAzure - 3);

   TCanvas canv(name.c_str(), "", 1800, 600);

   distr->Draw();
   cutDistr->Draw("SAME");

   system("mkdir -p output/Deadmaps/Run14HeAu200");
   ROOTTools::PrintCanvas(&canv, "output/Deadmaps/Run14HeAu200/" + name);
}

int main()
{
   gErrorIgnoreLevel = kWarning;
   gStyle->SetOptStat(0);

   TFile inputFile("data/Real/Run14HeAu200/SingleTrack/sum.root");

   Cut1D::CutDistribution(static_cast<TH1D *>(inputFile.Get("slat: TOFe")), "Slat", 0.4);
   Cut1D::CutDistribution(static_cast<TH1D *>(inputFile.Get("strip: TOFw")), "Strip", 0.4);
}

#endif /* CUT_1D_CPP */
