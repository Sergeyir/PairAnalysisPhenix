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

   simFileName = "data/PostSim/" + runName + "/SingleTrack/all.root";

   CppTools::CheckInputFile(simFileName);

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

      const double mult = distrMult->GetBinContent(1)/distrCentrality->Integral();

      averageMult += mult;
      averageChargeRatio += distrMult->GetBinContent(2)/distrMult->GetBinContent(3);
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

      if (averageMult/mult < multThreshold && averageMult/mult > 1./multThreshold &&
          averageChargeRatio/chargeRatio < multThreshold && 
          averageChargeRatio/chargeRatio > 1./multThreshold)
      {
         passedRuns.push_back(run);
      }
      else
      {
         badRuns.push_back(run);
      }
   }

   CppTools::PrintInfo(std::to_string(passedRuns.size()) + " runs out of " + 
                       std::to_string(goodRuns.size()) + " passed multiplicity check");

   goodRuns = passedRuns;

   grMult.SetMarkerStyle(7);
   grChargeRatio.SetMarkerStyle(7);

   TCanvas canv("canv", "", 1000, 600);
   canv.Divide(1, 2, 0., 0.);

   canv.cd(1);

   gPad->SetRightMargin(0.045); gPad->SetTopMargin(0.02); 
   gPad->SetLeftMargin(0.08); gPad->SetBottomMargin(0.);

   ROOTTools::DrawFrame(grMult.GetXaxis()->GetBinLowEdge(1), 
                        grMult.GetYaxis()->GetBinLowEdge(1),
                        grMult.GetXaxis()->GetBinUpEdge(grMult.GetXaxis()->GetNbins()),
                        grMult.GetYaxis()->GetBinUpEdge(grMult.GetYaxis()->GetNbins()), 
                        "", "", "N_{charged}/N_{evt}", 0., 0.6, 0., 0.07);

   ROOTTools::DrawLine(grMult.GetXaxis()->GetBinLowEdge(1), averageMult, 
                       grMult.GetXaxis()->GetBinUpEdge(grMult.GetXaxis()->GetNbins()), averageMult,
                       kGray + 2, 0.5, 2, 2);
   grMult.Draw("P");

   canv.cd(2);

   gPad->SetRightMargin(0.045); gPad->SetTopMargin(0.); 
   gPad->SetLeftMargin(0.08); gPad->SetBottomMargin(0.12);

   ROOTTools::DrawFrame(grChargeRatio.GetXaxis()->GetBinLowEdge(1), 
                        grChargeRatio.GetYaxis()->GetBinLowEdge(1),
                        grChargeRatio.GetXaxis()->GetBinUpEdge(grChargeRatio.GetXaxis()->GetNbins()),
                        grChargeRatio.GetYaxis()->GetBinUpEdge(grChargeRatio.GetYaxis()->GetNbins()), 
                        "", "run index", "N_{charged}^{+}/N_{charged}^{-}", 0.7, 0.6, 0.07, 0.07);

   ROOTTools::DrawLine(grChargeRatio.GetXaxis()->GetBinLowEdge(1), averageChargeRatio, 
                       grChargeRatio.GetXaxis()->GetBinUpEdge(grChargeRatio.GetXaxis()->GetNbins()), 
                       averageChargeRatio, kGray + 2, 0.5, 2, 2);

   grChargeRatio.Draw("P");

   ROOTTools::PrintCanvas(&canv, outputDir + "/mult");
}

void CheckRuns::CheckRunsByDCBoard()
{
   TFile *simFile = TFile::Open(simFileName.c_str());

   TH2F *heatmapSimDCe0 = static_cast<TH2F *>(simFile->Get("_Heatmap: DCe, zDC>=0"));
   TH2F *heatmapSimDCe1 = static_cast<TH2F *>(simFile->Get("_Heatmap: DCe, zDC<0"));
   TH2F *heatmapSimDCw0 = static_cast<TH2F *>(simFile->Get("_Heatmap: DCw, zDC>=0"));
   TH2F *heatmapSimDCw1 = static_cast<TH2F *>(simFile->Get("_Heatmap: DCw, zDC<0"));

   for (int i = 1; i <= heatmapSimDCe0->GetXaxis()->GetNbins(); i++)
   {
      for (int j = 1; j <= heatmapSimDCe0->GetYaxis()->GetNbins(); j++)
      {
         if (dmCutter.IsDeadDC(0, 1., heatmapSimDCe0->GetXaxis()->GetBinCenter(i),
                               heatmapSimDCe0->GetYaxis()->GetBinCenter(j)))
         {
            heatmapSimDCe0->SetBinContent(i, j, 0.);
         }
      }
   }

   for (int i = 1; i <= heatmapSimDCe1->GetXaxis()->GetNbins(); i++)
   {
      for (int j = 1; j <= heatmapSimDCe1->GetYaxis()->GetNbins(); j++)
      {
         if (dmCutter.IsDeadDC(0, -1., heatmapSimDCe1->GetXaxis()->GetBinCenter(i),
                               heatmapSimDCe1->GetYaxis()->GetBinCenter(j)))
         {
            heatmapSimDCe1->SetBinContent(i, j, 0.);
         }
      }
   }

   for (int i = 1; i <= heatmapSimDCw0->GetXaxis()->GetNbins(); i++)
   {
      for (int j = 1; j <= heatmapSimDCw0->GetYaxis()->GetNbins(); j++)
      {
         if (dmCutter.IsDeadDC(1, 1., heatmapSimDCw0->GetXaxis()->GetBinCenter(i),
                               heatmapSimDCw0->GetYaxis()->GetBinCenter(j)))
         {
            heatmapSimDCw0->SetBinContent(i, j, 0.);
         }
      }
   }

   for (int i = 1; i <= heatmapSimDCw1->GetXaxis()->GetNbins(); i++)
   {
      for (int j = 1; j <= heatmapSimDCw1->GetYaxis()->GetNbins(); j++)
      {
         if (dmCutter.IsDeadDC(1, -1., heatmapSimDCw1->GetXaxis()->GetBinCenter(i),
                               heatmapSimDCw1->GetYaxis()->GetBinCenter(j)))
         {
            heatmapSimDCw1->SetBinContent(i, j, 0.);
         }
      }
   }

   TH1D *projSimDCe0Board = heatmapSimDCe0->
      ProjectionY("sim DCe0 board", 1, heatmapSimDCe0->GetXaxis()->GetNbins());
   TH1D *projSimDCe1Board = heatmapSimDCe1->
      ProjectionY("sim DCe1 board", 1, heatmapSimDCe1->GetXaxis()->GetNbins());
   TH1D *projSimDCw0Board = heatmapSimDCw0->
      ProjectionY("sim DCw0 board", 1, heatmapSimDCw0->GetXaxis()->GetNbins());
   TH1D *projSimDCw1Board = heatmapSimDCw1->
      ProjectionY("sim DCw1 board", 1, heatmapSimDCw1->GetXaxis()->GetNbins());

   projSimDCe0Board->Scale(1./projSimDCe0Board->Integral());
   projSimDCe1Board->Scale(1./projSimDCe1Board->Integral());
   projSimDCw0Board->Scale(1./projSimDCw0Board->Integral());
   projSimDCw1Board->Scale(1./projSimDCw1Board->Integral());

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

      projDCe0Board->Divide(projSimDCe0Board);
      projDCe1Board->Divide(projSimDCe1Board);
      projDCw0Board->Divide(projSimDCw0Board);
      projDCw1Board->Divide(projSimDCw1Board);

      outputFile.mkdir(std::to_string(run).c_str());
      outputFile.cd(std::to_string(run).c_str());

      TF1 fit("const fit", "pol0");

      fit.SetParameter(0, 1.);
      fit.SetParLimits(0, 0., 2.);
      
      projDCe0Board->Fit(&fit, "QBMG");
      const double constParDCe0Board = fit.GetParameter(0);
      const double chi2NDFDCe0Board = fit.GetChisquare()/fit.GetNDF();
      projDCe0Board->Write();

      projDCe1Board->Fit(&fit, "QBMG");
      const double constParDCe1Board = fit.GetParameter(0);
      const double chi2NDFDCe1Board = fit.GetChisquare()/fit.GetNDF();
      projDCe1Board->Write();

      projDCw0Board->Fit(&fit, "QBMG");
      const double constParDCw0Board = fit.GetParameter(0);
      const double chi2NDFDCw0Board = fit.GetChisquare()/fit.GetNDF();
      projDCw0Board->Write();

      projDCw1Board->Fit(&fit, "QBMG");
      const double constParDCw1Board = fit.GetParameter(0);
      const double chi2NDFDCw1Board = fit.GetChisquare()/fit.GetNDF();
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
   }

   averageConstParDCe0Board /= static_cast<double>(goodRuns.size());
   averageConstParDCe1Board /= static_cast<double>(goodRuns.size());
   averageConstParDCw0Board /= static_cast<double>(goodRuns.size());
   averageConstParDCw1Board /= static_cast<double>(goodRuns.size());

   CppTools::PrintInfo(std::to_string(passedRuns.size()) + " runs out of " + 
                       std::to_string(goodRuns.size()) + " passed DC board check");

   goodRuns = passedRuns;

   grChi2NDFDCe0Board.SetMarkerStyle(7);
   grChi2NDFDCe1Board.SetMarkerStyle(7);
   grChi2NDFDCw0Board.SetMarkerStyle(7);
   grChi2NDFDCw1Board.SetMarkerStyle(7);

   grConstParDCe0Board.SetMarkerStyle(7);
   grConstParDCe1Board.SetMarkerStyle(7);
   grConstParDCw0Board.SetMarkerStyle(7);
   grConstParDCw1Board.SetMarkerStyle(7);

   TCanvas canv("canv", "", 1000, 600);
   canv.Divide(1, 2, 0., 0.);

   canv.cd(1);

   gPad->SetRightMargin(0.045); gPad->SetTopMargin(0.02); 
   gPad->SetLeftMargin(0.08); gPad->SetBottomMargin(0.);

   ROOTTools::DrawFrame(grConstParDCe0Board.GetXaxis()->GetBinLowEdge(1), 
                        grConstParDCe0Board.GetYaxis()->GetBinLowEdge(1),
                        grConstParDCe0Board.GetXaxis()->
                        GetBinUpEdge(grConstParDCe0Board.GetXaxis()->GetNbins()),
                        grConstParDCe0Board.GetYaxis()->
                        GetBinUpEdge(grConstParDCe0Board.GetYaxis()->GetNbins()), 
                        "", "run index", "const", 0., 0.6, 0., 0.07);

   ROOTTools::DrawLine(grConstParDCe0Board.GetXaxis()->GetBinLowEdge(1), averageConstParDCe0Board, 
                       grConstParDCe0Board.GetXaxis()->GetBinUpEdge(grConstParDCe0Board.GetXaxis()->
                                                                    GetNbins()), 
                       averageConstParDCe0Board, kGray + 2, 0.5, 2, 2);

   grConstParDCe0Board.Draw("P");

   canv.cd(2);

   gPad->SetRightMargin(0.045); gPad->SetTopMargin(0.); 
   gPad->SetLeftMargin(0.08); gPad->SetBottomMargin(0.12);

   ROOTTools::DrawFrame(grChi2NDFDCe0Board.GetXaxis()->GetBinLowEdge(1), 
                        grChi2NDFDCe0Board.GetYaxis()->GetBinLowEdge(1),
                        grChi2NDFDCe0Board.GetXaxis()->
                        GetBinUpEdge(grChi2NDFDCe0Board.GetXaxis()->GetNbins()),
                        grChi2NDFDCe0Board.GetYaxis()->
                        GetBinUpEdge(grChi2NDFDCe0Board.GetYaxis()->GetNbins()), 
                        "", "run index", "#chi^{2}/NDF", 0.7, 0.6, 0.07, 0.07);

   grChi2NDFDCe0Board.Draw("P");

   ROOTTools::PrintCanvas(&canv, outputDir + "/DCe0");

   canv.Clear();

   canv.Divide(1, 2, 0., 0.);

   canv.cd(1);

   gPad->SetRightMargin(0.045); gPad->SetTopMargin(0.02); 
   gPad->SetLeftMargin(0.08); gPad->SetBottomMargin(0.);

   ROOTTools::DrawFrame(grConstParDCe1Board.GetXaxis()->GetBinLowEdge(1), 
                        grConstParDCe1Board.GetYaxis()->GetBinLowEdge(1),
                        grConstParDCe1Board.GetXaxis()->
                        GetBinUpEdge(grConstParDCe1Board.GetXaxis()->GetNbins()),
                        grConstParDCe1Board.GetYaxis()->
                        GetBinUpEdge(grConstParDCe1Board.GetYaxis()->GetNbins()), 
                        "", "run index", "const", 0., 0.6, 0., 0.07);

   ROOTTools::DrawLine(grConstParDCe1Board.GetXaxis()->GetBinLowEdge(1), averageConstParDCe1Board, 
                       grConstParDCe1Board.GetXaxis()->GetBinUpEdge(grConstParDCe1Board.GetXaxis()->
                                                                    GetNbins()), 
                       averageConstParDCe1Board, kGray + 2, 0.5, 2, 2);

   grConstParDCe1Board.Draw("P");

   canv.cd(2);

   gPad->SetRightMargin(0.045); gPad->SetTopMargin(0.); 
   gPad->SetLeftMargin(0.08); gPad->SetBottomMargin(0.12);

   ROOTTools::DrawFrame(grChi2NDFDCe1Board.GetXaxis()->GetBinLowEdge(1), 
                        grChi2NDFDCe1Board.GetYaxis()->GetBinLowEdge(1),
                        grChi2NDFDCe1Board.GetXaxis()->
                        GetBinUpEdge(grChi2NDFDCe1Board.GetXaxis()->GetNbins()),
                        grChi2NDFDCe1Board.GetYaxis()->
                        GetBinUpEdge(grChi2NDFDCe1Board.GetYaxis()->GetNbins()), 
                        "", "run index", "#chi^{2}/NDF", 0.7, 0.6, 0.07, 0.07);

   grChi2NDFDCe1Board.Draw("P");

   ROOTTools::PrintCanvas(&canv, outputDir + "/DCe1");

   canv.Clear();

   canv.Divide(1, 2, 0., 0.);

   canv.cd(1);

   gPad->SetRightMargin(0.045); gPad->SetTopMargin(0.02); 
   gPad->SetLeftMargin(0.08); gPad->SetBottomMargin(0.);

   ROOTTools::DrawFrame(grConstParDCw0Board.GetXaxis()->GetBinLowEdge(1), 
                        grConstParDCw0Board.GetYaxis()->GetBinLowEdge(1),
                        grConstParDCw0Board.GetXaxis()->
                        GetBinUpEdge(grConstParDCw0Board.GetXaxis()->GetNbins()),
                        grConstParDCw0Board.GetYaxis()->
                        GetBinUpEdge(grConstParDCw0Board.GetYaxis()->GetNbins()), 
                        "", "run index", "const", 0., 0.6, 0., 0.07);

   ROOTTools::DrawLine(grConstParDCw0Board.GetXaxis()->GetBinLowEdge(1), averageConstParDCw0Board, 
                       grConstParDCw0Board.GetXaxis()->GetBinUpEdge(grConstParDCw0Board.GetXaxis()->
                                                                    GetNbins()), 
                       averageConstParDCw0Board, kGray + 2, 0.5, 2, 2);

   grConstParDCw0Board.Draw("P");

   canv.cd(2);

   gPad->SetRightMargin(0.045); gPad->SetTopMargin(0.); 
   gPad->SetLeftMargin(0.08); gPad->SetBottomMargin(0.12);

   ROOTTools::DrawFrame(grChi2NDFDCw0Board.GetXaxis()->GetBinLowEdge(1), 
                        grChi2NDFDCw0Board.GetYaxis()->GetBinLowEdge(1),
                        grChi2NDFDCw0Board.GetXaxis()->
                        GetBinUpEdge(grChi2NDFDCw0Board.GetXaxis()->GetNbins()),
                        grChi2NDFDCw0Board.GetYaxis()->
                        GetBinUpEdge(grChi2NDFDCw0Board.GetYaxis()->GetNbins()), 
                        "", "run index", "#chi^{2}/NDF", 0.7, 0.6, 0.07, 0.07);

   grChi2NDFDCw0Board.Draw("P");

   ROOTTools::PrintCanvas(&canv, outputDir + "/DCw0");

   canv.Clear();

   canv.Divide(1, 2, 0., 0.);

   canv.cd(1);

   gPad->SetRightMargin(0.045); gPad->SetTopMargin(0.02); 
   gPad->SetLeftMargin(0.08); gPad->SetBottomMargin(0.);

   ROOTTools::DrawFrame(grConstParDCw1Board.GetXaxis()->GetBinLowEdge(1), 
                        grConstParDCw1Board.GetYaxis()->GetBinLowEdge(1),
                        grConstParDCw1Board.GetXaxis()->
                        GetBinUpEdge(grConstParDCw1Board.GetXaxis()->GetNbins()),
                        grConstParDCw1Board.GetYaxis()->
                        GetBinUpEdge(grConstParDCw1Board.GetYaxis()->GetNbins()), 
                        "", "run index", "const", 0., 0.6, 0., 0.07);

   ROOTTools::DrawLine(grConstParDCw1Board.GetXaxis()->GetBinLowEdge(1), averageConstParDCw1Board, 
                       grConstParDCw1Board.GetXaxis()->GetBinUpEdge(grConstParDCw1Board.GetXaxis()->
                                                                    GetNbins()), 
                       averageConstParDCw1Board, kGray + 2, 0.5, 2, 2);

   grConstParDCw1Board.Draw("P");

   canv.cd(2);

   gPad->SetRightMargin(0.045); gPad->SetTopMargin(0.); 
   gPad->SetLeftMargin(0.08); gPad->SetBottomMargin(0.12);

   ROOTTools::DrawFrame(grChi2NDFDCw1Board.GetXaxis()->GetBinLowEdge(1), 
                        grChi2NDFDCw1Board.GetYaxis()->GetBinLowEdge(1),
                        grChi2NDFDCw1Board.GetXaxis()->
                        GetBinUpEdge(grChi2NDFDCw1Board.GetXaxis()->GetNbins()),
                        grChi2NDFDCw1Board.GetYaxis()->
                        GetBinUpEdge(grChi2NDFDCw1Board.GetYaxis()->GetNbins()), 
                        "", "run index", "#chi^{2}/NDF", 0.7, 0.6, 0.07, 0.07);

   grChi2NDFDCw1Board.Draw("P");

   ROOTTools::PrintCanvas(&canv, outputDir + "/DCw1");
}

/*
void CheckRuns::CheckRunsByPC1()
{
   TFile *simFile = TFile::Open(simFileName.c_str());

   TH2F *heatmapSimPC1e = static_cast<TH2F *>(simFile->Get("Heatmap: PC1e"));
   TH2F *heatmapSimPC1w = static_cast<TH2F *>(simFile->Get("Heatmap: PC1e"));

   for (int i = 1; i <= heatmapPC1e->GetXaxis()->GetNbins(); i++)
   {
      for (int j = 1; j <= heatmapPC1e->GetYaxis()->GetNbins(); j++)
      {
         if (dmCutter.IsDeadPC1(0, heatmapSimPC1e->GetXaxis()->GetBinCenter(i),
                                heatmapSimPC1e->GetYaxis()->GetBinCenter(j)))
         {
            heatmapSimPC1e->SetBinContent(i, j, 0.);
         }
      }
   }

   for (int i = 1; i <= heatmapPC1w->GetXaxis()->GetNbins(); i++)
   {
      for (int j = 1; j <= heatmapPC1w->GetYaxis()->GetNbins(); j++)
      {
         if (dmCutter.IsDeadPC1(0, heatmapSimPC1w->GetXaxis()->GetBinCenter(i),
                                heatmapSimPC1w->GetYaxis()->GetBinCenter(j)))
         {
            heatmapSimPC1w->SetBinContent(i, j, 0.);
         }
      }
   }

}
*/

#endif /* CHECK_RUNS_CPP */
