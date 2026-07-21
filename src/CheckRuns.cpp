/** 
 *  @file   CheckRuns.cpp
 *  @brief  Contains realisations of functions and variables that are used for the determination of good runs 
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysisPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef CHECK_RUNS_CPP
#define CHECK_RUNS_CPP

#include "CheckRuns.hpp"

using namespace CheckRuns;

int main(int argc, char **argv)
{
   if (argc < 2 || argc > 4) 
   {
      CppTools::PrintError("Expected 1-3 parameters while " + std::to_string(argc - 1) + " "\
                           "parameter(s) were provided \n Usage: bin/CheckRuns "\
                           "inputYAMLMain multThreshold=2. chi2NDFThreshold=3.");
   }

   CppTools::CheckInputFile(argv[1]);

   if (argc > 2) 
   {
      multThreshold = std::stof(argv[2]);
      if (argc > 3) 
      {
         chi2NDFThreshold = std::stof(argv[3]);
      }
   }

   gStyle->SetOptStat(0);
   gErrorIgnoreLevel = kWarning;

   inputYAMLMain.OpenFile(argv[1], "main");
   inputYAMLMain.CheckStatus("main");

   runName = inputYAMLMain["run_name"].as<std::string>();

   dmCutter.Initialize(runName, inputYAMLMain["detectors_configuration"].as<std::string>());

   sumFileName = "data/Real/" + runName + "/SingleTrack/sum.root";

   CppTools::CheckInputFile(sumFileName);

   inputDir = "data/Real/" + runName + "/SingleTrack";

   outputDir = "output/Runs/" + runName;
   std::filesystem::create_directories(outputDir);

   if (!std::filesystem::exists(inputDir))
   {
      CppTools::PrintError("No such directory" + inputDir);
   }

   for (const auto& file : std::filesystem::directory_iterator(inputDir))
   {
      const std::string fileName = static_cast<std::string>(file.path());
      if (!std::regex_match(fileName, std::regex("(.*)se-[0-9]{6}\\.root"))) continue;

      std::smatch m;
      std::regex_search(fileName, m, std::regex("[0-9]{6}"));
      runs.emplace_back(std::atoi(m.str().c_str()));
   }

   // At first assuming all runs to be good. Bad runs will be removed from this list in later checks
   goodRuns = runs;
   CheckRunsByMultiplicity();

   CheckRunsByDCBoard();

   CppTools::PrintInfo("Pictures and .root files were written in " + outputDir);

   std::filesystem::create_directories("data/GoodRunsList/");

   std::ofstream goodRunsListFile("data/GoodRunsList/" + runName + ".txt");
   // for se-xxxxxx.root files
   std::ofstream goodRunsListFileForSE("data/GoodRunsList/" + runName + "_for_se.txt");
   // for cabana boy output
   std::ofstream goodRunsListFileForCC("data/GoodRunsList/" + runName + "_for_cc.txt");

   for (const int& run : goodRuns)
   {
      goodRunsListFile << run << std::endl;
      goodRunsListFileForSE << "se-" << run << ".root" << std::endl;
      goodRunsListFileForCC << run << ".root" << std::endl;
   }
}

void CheckRuns::CheckRunsByMultiplicity()
{
   double averageMult = 0.;
   double averageChargeRatio = 0.;

   for (const int &run : goodRuns)
   {
      const std::string fileName = inputDir + "/se-" + std::to_string(run) + ".root";
      TFile *inputFile = TFile::Open(fileName.c_str());

      TH1D *distrMult = static_cast<TH1D *>(inputFile->Get("multiplicity"));
      TH1D *distrCentrality = static_cast<TH1D *>(inputFile->Get("centrality"));

      if (distrMult->GetBinContent(1) < 1e-6 || distrCentrality->Integral() < 1e-6)
      {
         CppTools::PrintWarning("Run " + std::to_string(run) + " is possibly empty");
         continue;
      }

      const double mult = distrMult->GetBinContent(1)/distrCentrality->Integral();

      averageMult += mult;
      averageChargeRatio += distrMult->GetBinContent(2)/distrMult->GetBinContent(3);

      inputFile->Close();
   }

   averageMult /= static_cast<double>(goodRuns.size());
   averageChargeRatio /= static_cast<double>(goodRuns.size());

   // graphs for multiplicity and charge ratio
   TGraph grMult;
   TGraph grChargeRatio;

   std::vector<int> passedRuns;

   for (const int &run : goodRuns)
   {
      const std::string fileName = inputDir + "/se-" + std::to_string(run) + ".root";
      TFile *inputFile = TFile::Open(fileName.c_str());

      TH1D *distrMult = static_cast<TH1D *>(inputFile->Get("multiplicity"));
      TH1D *distrCentrality = static_cast<TH1D *>(inputFile->Get("centrality"));

      const double mult = distrMult->GetBinContent(1)/distrCentrality->Integral();
      const double chargeRatio = distrMult->GetBinContent(2)/distrMult->GetBinContent(3);

      grMult.AddPoint(static_cast<double>(run), mult);
      grChargeRatio.AddPoint(static_cast<double>(run), chargeRatio);

      if (!std::isnan(mult) && !std::isnan(chargeRatio) && 
          averageMult/mult < multThreshold && averageMult/mult > 1./multThreshold &&
          averageChargeRatio/chargeRatio < multThreshold && 
          averageChargeRatio/chargeRatio > 1./multThreshold)
      {
         passedRuns.push_back(run);
      }
      else
      {
         badRuns.push_back(run);
      }
      inputFile->Close();
   }

   CppTools::PrintInfo(std::to_string(passedRuns.size()) + " runs out of " + 
                       std::to_string(goodRuns.size()) + " passed multiplicity check");

   goodRuns = passedRuns;

   grMult.SetMarkerStyle(7);
   grChargeRatio.SetMarkerStyle(7);

   TCanvas canv("canv", "", 1000, 400);
   canv.Divide(1, 2, 0., 0.);

   canv.cd(1);

   gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.01); 
   gPad->SetLeftMargin(0.065); gPad->SetBottomMargin(0.);

   ROOTTools::DrawFrame(grMult.GetXaxis()->GetBinLowEdge(1), 
                        grMult.GetYaxis()->GetBinLowEdge(1),
                        grMult.GetXaxis()->GetBinUpEdge(grMult.GetXaxis()->GetNbins()),
                        grMult.GetYaxis()->GetBinUpEdge(grMult.GetYaxis()->GetNbins()), 
                        "", "", "N_{charged}/N_{evt}", 0., 0.4, 0., 0.08);

   ROOTTools::DrawLine(grMult.GetXaxis()->GetBinLowEdge(1), averageMult, 
                       grMult.GetXaxis()->GetBinUpEdge(grMult.GetXaxis()->GetNbins()), averageMult,
                       kGray + 2, 0.5, 2, 2);
   grMult.Draw("P");

   canv.cd(2);

   gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.); 
   gPad->SetLeftMargin(0.065); gPad->SetBottomMargin(0.14);

   ROOTTools::DrawFrame(grChargeRatio.GetXaxis()->GetBinLowEdge(1), 
                        grChargeRatio.GetYaxis()->GetBinLowEdge(1),
                        grChargeRatio.GetXaxis()->GetBinUpEdge(grChargeRatio.GetXaxis()->GetNbins()),
                        grChargeRatio.GetYaxis()->GetBinUpEdge(grChargeRatio.GetYaxis()->GetNbins()), 
                        "", "run index", "N_{charged}^{+}/N_{charged}^{-}", 0.7, 0.4, 0.08, 0.08);

   ROOTTools::DrawLine(grChargeRatio.GetXaxis()->GetBinLowEdge(1), averageChargeRatio, 
                       grChargeRatio.GetXaxis()->GetBinUpEdge(grChargeRatio.GetXaxis()->GetNbins()), 
                       averageChargeRatio, kGray + 2, 0.5, 2, 2);

   grChargeRatio.Draw("P");

   ROOTTools::PrintCanvas(&canv, outputDir + "/mult");
}

void CheckRuns::CheckRunsByDCBoard()
{
   TFile *sumFile = TFile::Open(sumFileName.c_str());

   TH2F *heatmapSumDCe0 = static_cast<TH2F *>(sumFile->Get("_Heatmap: DCe, zDC>=0"));
   TH2F *heatmapSumDCe1 = static_cast<TH2F *>(sumFile->Get("_Heatmap: DCe, zDC<0"));
   TH2F *heatmapSumDCw0 = static_cast<TH2F *>(sumFile->Get("_Heatmap: DCw, zDC>=0"));
   TH2F *heatmapSumDCw1 = static_cast<TH2F *>(sumFile->Get("_Heatmap: DCw, zDC<0"));

   for (int i = 1; i <= heatmapSumDCe0->GetXaxis()->GetNbins(); i++)
   {
      for (int j = 1; j <= heatmapSumDCe0->GetYaxis()->GetNbins(); j++)
      {
         if (dmCutter.IsDeadDC(0, 1., heatmapSumDCe0->GetXaxis()->GetBinCenter(i),
                               heatmapSumDCe0->GetYaxis()->GetBinCenter(j)))
         {
            heatmapSumDCe0->SetBinContent(i, j, 0.);
         }
      }
   }

   for (int i = 1; i <= heatmapSumDCe1->GetXaxis()->GetNbins(); i++)
   {
      for (int j = 1; j <= heatmapSumDCe1->GetYaxis()->GetNbins(); j++)
      {
         if (dmCutter.IsDeadDC(0, -1., heatmapSumDCe1->GetXaxis()->GetBinCenter(i),
                               heatmapSumDCe1->GetYaxis()->GetBinCenter(j)))
         {
            heatmapSumDCe1->SetBinContent(i, j, 0.);
         }
      }
   }

   for (int i = 1; i <= heatmapSumDCw0->GetXaxis()->GetNbins(); i++)
   {
      for (int j = 1; j <= heatmapSumDCw0->GetYaxis()->GetNbins(); j++)
      {
         if (dmCutter.IsDeadDC(1, 1., heatmapSumDCw0->GetXaxis()->GetBinCenter(i),
                               heatmapSumDCw0->GetYaxis()->GetBinCenter(j)))
         {
            heatmapSumDCw0->SetBinContent(i, j, 0.);
         }
      }
   }

   for (int i = 1; i <= heatmapSumDCw1->GetXaxis()->GetNbins(); i++)
   {
      for (int j = 1; j <= heatmapSumDCw1->GetYaxis()->GetNbins(); j++)
      {
         if (dmCutter.IsDeadDC(1, -1., heatmapSumDCw1->GetXaxis()->GetBinCenter(i),
                               heatmapSumDCw1->GetYaxis()->GetBinCenter(j)))
         {
            heatmapSumDCw1->SetBinContent(i, j, 0.);
         }
      }
   }

   TH1D *projSumDCe0Board = heatmapSumDCe0->
      ProjectionY("sum DCe0 board", 1, heatmapSumDCe0->GetXaxis()->GetNbins());
   TH1D *projSumDCe1Board = heatmapSumDCe1->
      ProjectionY("sum DCe1 board", 1, heatmapSumDCe1->GetXaxis()->GetNbins());
   TH1D *projSumDCw0Board = heatmapSumDCw0->
      ProjectionY("sum DCw0 board", 1, heatmapSumDCw0->GetXaxis()->GetNbins());
   TH1D *projSumDCw1Board = heatmapSumDCw1->
      ProjectionY("sum DCw1 board", 1, heatmapSumDCw1->GetXaxis()->GetNbins());

   projSumDCe0Board->Scale(1./projSumDCe0Board->Integral());
   projSumDCe1Board->Scale(1./projSumDCe1Board->Integral());
   projSumDCw0Board->Scale(1./projSumDCw0Board->Integral());
   projSumDCw1Board->Scale(1./projSumDCw1Board->Integral());

   // graphs for chi2/NDF
   TGraph grChi2NDFDCe0Board;
   TGraph grChi2NDFDCe1Board;
   TGraph grChi2NDFDCw0Board;
   TGraph grChi2NDFDCw1Board;

   // graphs for constant fit parameter
   TGraph grConstParDCe0Board;
   TGraph grConstParDCe1Board;
   TGraph grConstParDCw0Board;
   TGraph grConstParDCw1Board;

   TFile outputFile((outputDir + "/DC.root").c_str(), "RECREATE");

   std::vector<int> passedRuns;

   double averageConstParDCe0Board = 0.;
   double averageConstParDCe1Board = 0.;
   double averageConstParDCw0Board = 0.;
   double averageConstParDCw1Board = 0.;

   for (const int &run : goodRuns)
   {
      const std::string fileName = inputDir + "/se-" + std::to_string(run) + ".root";
      TFile *inputFile = TFile::Open(fileName.c_str());

      TH2F *heatmapDCe0 = static_cast<TH2F *>(inputFile->Get("_Heatmap: DCe, zDC>=0"));
      TH2F *heatmapDCe1 = static_cast<TH2F *>(inputFile->Get("_Heatmap: DCe, zDC<0"));
      TH2F *heatmapDCw0 = static_cast<TH2F *>(inputFile->Get("_Heatmap: DCw, zDC>=0"));
      TH2F *heatmapDCw1 = static_cast<TH2F *>(inputFile->Get("_Heatmap: DCw, zDC<0"));

      for (int i = 1; i <= heatmapDCe0->GetXaxis()->GetNbins(); i++)
      {
         for (int j = 1; j <= heatmapDCe0->GetYaxis()->GetNbins(); j++)
         {
            if (dmCutter.IsDeadDC(0, 1., heatmapDCe0->GetXaxis()->GetBinCenter(i),
                                  heatmapDCe0->GetYaxis()->GetBinCenter(j)))
            {
               heatmapDCe0->SetBinContent(i, j, 0.);
            }
         }
      }

      for (int i = 1; i <= heatmapDCe1->GetXaxis()->GetNbins(); i++)
      {
         for (int j = 1; j <= heatmapDCe1->GetYaxis()->GetNbins(); j++)
         {
            if (dmCutter.IsDeadDC(0, -1., heatmapDCe1->GetXaxis()->GetBinCenter(i),
                                  heatmapDCe1->GetYaxis()->GetBinCenter(j)))
            {
               heatmapDCe1->SetBinContent(i, j, 0.);
            }
         }
      }

      for (int i = 1; i <= heatmapDCw0->GetXaxis()->GetNbins(); i++)
      {
         for (int j = 1; j <= heatmapDCw0->GetYaxis()->GetNbins(); j++)
         {
            if (dmCutter.IsDeadDC(1, 1., heatmapDCw0->GetXaxis()->GetBinCenter(i),
                                  heatmapDCw0->GetYaxis()->GetBinCenter(j)))
            {
               heatmapDCw0->SetBinContent(i, j, 0.);
            }
         }
      }

      for (int i = 1; i <= heatmapDCw1->GetXaxis()->GetNbins(); i++)
      {
         for (int j = 1; j <= heatmapDCw1->GetYaxis()->GetNbins(); j++)
         {
            if (dmCutter.IsDeadDC(1, -1., heatmapDCw1->GetXaxis()->GetBinCenter(i),
                                  heatmapDCw1->GetYaxis()->GetBinCenter(j)))
            {
               heatmapDCw1->SetBinContent(i, j, 0.);
            }
         }
      }

      TH1D *projDCe0Board = heatmapDCe0->
         ProjectionY("proj DCe0 board to ref ratio", 1, heatmapDCe0->GetXaxis()->GetNbins());
      TH1D *projDCe1Board = heatmapDCe1->
         ProjectionY("proj DCe1 board to ref ratio", 1, heatmapDCe1->GetXaxis()->GetNbins());
      TH1D *projDCw0Board = heatmapDCw0->
         ProjectionY("proj DCw0 board to ref ratio", 1, heatmapDCw0->GetXaxis()->GetNbins());
      TH1D *projDCw1Board = heatmapDCw1->
         ProjectionY("proj DCw1 board to ref ratio", 1, heatmapDCw1->GetXaxis()->GetNbins());

      projDCe0Board->Scale(1./projDCe0Board->Integral());
      projDCe1Board->Scale(1./projDCe1Board->Integral());
      projDCw0Board->Scale(1./projDCw0Board->Integral());
      projDCw1Board->Scale(1./projDCw1Board->Integral());

      projDCe0Board->Divide(projSumDCe0Board);
      projDCe1Board->Divide(projSumDCe1Board);
      projDCw0Board->Divide(projSumDCw0Board);
      projDCw1Board->Divide(projSumDCw1Board);

      outputFile.mkdir(std::to_string(run).c_str());
      outputFile.cd(std::to_string(run).c_str());

      TF1 fit("const fit", "pol0");

      fit.SetParameter(0, 1.);
      //fit.SetParLimits(0, 0., 2.);

      fit.FixParameter(0, GetYWeightedAverage(projDCe0Board));
      projDCe0Board->Fit(&fit, "QBP");
      const double constParDCe0Board = fit.GetParameter(0);
      const double chi2NDFDCe0Board = GetChi2NDF(projDCe0Board, &fit);
      projDCe0Board->Write();

      fit.FixParameter(0, GetYWeightedAverage(projDCe1Board));
      projDCe1Board->Fit(&fit, "QBP");
      const double constParDCe1Board = fit.GetParameter(0);
      const double chi2NDFDCe1Board = GetChi2NDF(projDCe1Board, &fit);
      projDCe1Board->Write();

      fit.FixParameter(0, GetYWeightedAverage(projDCw0Board));
      projDCw0Board->Fit(&fit, "QBP");
      const double constParDCw0Board = fit.GetParameter(0);
      const double chi2NDFDCw0Board = GetChi2NDF(projDCw0Board, &fit);
      projDCw0Board->Write();

      fit.FixParameter(0, GetYWeightedAverage(projDCw1Board));
      projDCw1Board->Fit(&fit, "QBP");
      const double constParDCw1Board = fit.GetParameter(0);
      const double chi2NDFDCw1Board = GetChi2NDF(projDCw1Board, &fit);
      projDCw1Board->Write();

      grChi2NDFDCe0Board.AddPoint(run, chi2NDFDCe0Board);
      grChi2NDFDCe1Board.AddPoint(run, chi2NDFDCe1Board);
      grChi2NDFDCw0Board.AddPoint(run, chi2NDFDCw0Board);
      grChi2NDFDCw1Board.AddPoint(run, chi2NDFDCw1Board);

      grConstParDCe0Board.AddPoint(run, constParDCe0Board);
      grConstParDCe1Board.AddPoint(run, constParDCe1Board);
      grConstParDCw0Board.AddPoint(run, constParDCw0Board);
      grConstParDCw1Board.AddPoint(run, constParDCw1Board);

      averageConstParDCe0Board += constParDCe0Board;
      averageConstParDCe1Board += constParDCe1Board;
      averageConstParDCw0Board += constParDCw0Board;
      averageConstParDCw1Board += constParDCw1Board;

      if (chi2NDFDCe0Board < chi2NDFThreshold && chi2NDFDCe1Board < chi2NDFThreshold &&
          chi2NDFDCw0Board < chi2NDFThreshold && chi2NDFDCw1Board < chi2NDFThreshold)
      {
         passedRuns.push_back(run);
      }
      else
      {
         badRuns.push_back(run);
      }
      inputFile->Close();
   }

   averageConstParDCe0Board /= static_cast<double>(goodRuns.size());
   averageConstParDCe1Board /= static_cast<double>(goodRuns.size());
   averageConstParDCw0Board /= static_cast<double>(goodRuns.size());
   averageConstParDCw1Board /= static_cast<double>(goodRuns.size());

   CppTools::PrintInfo(std::to_string(passedRuns.size()) + " runs out of " + 
                       std::to_string(goodRuns.size()) + " passed DC board check");

   goodRuns = passedRuns;

   grConstParDCe0Board.SetMarkerStyle(7);
   grConstParDCe1Board.SetMarkerStyle(7);
   grConstParDCw0Board.SetMarkerStyle(7);
   grConstParDCw1Board.SetMarkerStyle(7);

   grChi2NDFDCe0Board.SetMarkerStyle(7);
   grChi2NDFDCe1Board.SetMarkerStyle(7);
   grChi2NDFDCw0Board.SetMarkerStyle(7);
   grChi2NDFDCw1Board.SetMarkerStyle(7);

   grConstParDCe0Board.SetMarkerColorAlpha(kP6Red, 0.8);
   grConstParDCe1Board.SetMarkerColorAlpha(kP6Blue, 0.8);
   grConstParDCw0Board.SetMarkerColorAlpha(kP6Violet, 0.8);
   grConstParDCw1Board.SetMarkerColorAlpha(kP6Gray, 0.8);

   grChi2NDFDCe0Board.SetMarkerColorAlpha(kP6Red, 0.8);
   grChi2NDFDCe1Board.SetMarkerColorAlpha(kP6Blue, 0.8);
   grChi2NDFDCw0Board.SetMarkerColorAlpha(kP6Violet, 0.8);
   grChi2NDFDCw1Board.SetMarkerColorAlpha(kP6Gray, 0.8);

   TCanvas canv("canv", "", 1000, 400);
   canv.Divide(1, 2, 0., 0.);

   canv.cd(1);

   gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.01); 
   gPad->SetLeftMargin(0.065); gPad->SetBottomMargin(0.);

   const int nXBins = grConstParDCe0Board.GetXaxis()->GetNbins();
   const int nYBins = grConstParDCe0Board.GetYaxis()->GetNbins();
   const double xMin = grConstParDCe0Board.GetXaxis()->GetBinLowEdge(1);
   const double xMax = grConstParDCe0Board.GetXaxis()->GetBinUpEdge(nXBins);
   const double yMin = CppTools::Minimum(grConstParDCe0Board.GetYaxis()->GetBinLowEdge(1),
                                         grConstParDCe1Board.GetYaxis()->GetBinLowEdge(1),
                                         grConstParDCw0Board.GetYaxis()->GetBinLowEdge(1),
                                         grConstParDCw1Board.GetYaxis()->GetBinLowEdge(1));
   const double yMax = CppTools::Maximum(grConstParDCe0Board.GetYaxis()->GetBinUpEdge(nYBins),
                                         grConstParDCe1Board.GetYaxis()->GetBinUpEdge(nYBins),
                                         grConstParDCw0Board.GetYaxis()->GetBinUpEdge(nYBins),
                                         grConstParDCw1Board.GetYaxis()->GetBinUpEdge(nYBins));

   ROOTTools::DrawFrame(xMin, yMin, xMax, yMax, "", "run index", "const", 0., 0.4, 0., 0.08);

   ROOTTools::DrawLine(xMin, averageConstParDCe0Board, xMax, averageConstParDCe0Board, 
                       kP6Red, 0.5, 2, 2);
   ROOTTools::DrawLine(xMin, averageConstParDCe1Board, xMax, averageConstParDCe1Board, 
                       kP6Blue, 0.5, 2, 2);
   ROOTTools::DrawLine(xMin, averageConstParDCw0Board, xMax, averageConstParDCw0Board, 
                       kP6Violet, 0.5, 2, 2);
   ROOTTools::DrawLine(xMin, averageConstParDCw1Board, xMax, averageConstParDCw1Board, 
                       kP6Gray, 0.5, 2, 2);

   grConstParDCe0Board.Draw("P");
   grConstParDCe1Board.Draw("P");
   grConstParDCw0Board.Draw("P");
   grConstParDCw1Board.Draw("P");

   canv.cd(2);

   gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.); 
   gPad->SetLeftMargin(0.065); gPad->SetBottomMargin(0.14);

   const double chi2NDFMax = 
      CppTools::Maximum(grChi2NDFDCe0Board.GetYaxis()->GetBinUpEdge(nYBins),
                        grChi2NDFDCe1Board.GetYaxis()->GetBinUpEdge(nYBins),
                        grChi2NDFDCw0Board.GetYaxis()->GetBinUpEdge(nYBins),
                        grChi2NDFDCw1Board.GetYaxis()->GetBinUpEdge(nYBins), 3.);

   ROOTTools::DrawFrame(xMin, 0., xMax, chi2NDFMax, 
                        "", "run index", "#chi^{2}/NDF", 0.7, 0.4, 0.08, 0.08);

   grChi2NDFDCe0Board.Draw("P");
   grChi2NDFDCe1Board.Draw("P");
   grChi2NDFDCw0Board.Draw("P");
   grChi2NDFDCw1Board.Draw("P");

   ROOTTools::PrintCanvas(&canv, outputDir + "/DC");
}

double CheckRuns::GetYWeightedAverage(TH1D *hist)
{
   double result = 0.;
   double weightSum = 0.;

   for (int i = 1; i <= hist->GetXaxis()->GetNbins(); i++)
   {
      if (hist->GetBinContent(i) < 1e-7) continue;

      const double weight = 1./pow(hist->GetBinError(i), 2);

      result += hist->GetBinContent(i)*weight;
      weightSum += weight;
   }

   return result/weightSum;
}

double CheckRuns::GetChi2NDF(TH1D *hist, TF1 *fit)
{
   double result = 0.;
   unsigned int numberOfNonZeroPoints = 0;

   for (int i = 1; i <= hist->GetXaxis()->GetNbins(); i++)
   {
      if (hist->GetBinContent(i) < 1e-7) continue;

      result += pow((hist->GetBinContent(i) - 
                     fit->Eval(hist->GetXaxis()->GetBinCenter(i))), 2)/hist->GetBinError(i);
      numberOfNonZeroPoints++;
   }

   return result/static_cast<double>(numberOfNonZeroPoints - fit->GetNpar());
}

#endif /* CHECK_RUNS_CPP */
