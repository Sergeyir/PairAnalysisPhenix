/** 
 *  @file   FindGoodRuns.cpp
 *  @brief  Contains realisations of functions and variables that are used for the determination of good runs 
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysisPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef FIND_GOOD_RUNS_CPP
#define FIND_GOOD_RUNS_CPP

#include "FindGoodRuns.hpp"

using namespace FindGoodRuns;

int main(int argc, char **argv)
{
   if (argc < 2 || argc > 4) 
   {
      CppTools::PrintError("Expected 1-3 parameters while " + std::to_string(argc - 1) + " "\
                           "parameter(s) were provided \n Usage: bin/FindReferenceRun "\
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
   CheckRunsByMultiplicityAndFindReferenceRun();

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

void FindGoodRuns::CheckRunsByMultiplicityAndFindReferenceRun()
{
   double averageMult = 0.;
   double averageChargeRatio = 0.;

   double referenceRunMult = 0.;

   for (const int &run : goodRuns)
   {
      const std::string fileName = inputDir + "/se-" + std::to_string(run) + ".root";
      TFile *inputFile = TFile::Open(fileName.c_str());

      TH1D *distrMult = static_cast<TH1D *>(inputFile->Get("multiplicity"));
      TH1D *distrCentrality = static_cast<TH1D *>(inputFile->Get("centrality"));

      const double mult = distrMult->GetBinContent(1)/distrCentrality->Integral();

      averageMult += mult;
      averageChargeRatio += distrMult->GetBinContent(2)/distrMult->GetBinContent(3);

      if (mult > referenceRunMult)
      {
         referenceRunMult = mult;
         referenceRun = run;
      }
   }
   CppTools::PrintInfo("Reference run : " + std::to_string(referenceRun) + 
                       " with mult/event value " + std::to_string(referenceRunMult));

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

      if (referenceRunMult/mult < multThreshold && referenceRunMult/mult > 1./multThreshold &&
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

void FindGoodRuns::CheckRunsByDCBoard()
{
   TFile *referenceRunFile = TFile::Open((inputDir + "/se-" + 
                                          std::to_string(referenceRun) + ".root").c_str());

   TH2F *heatmapReferenceDCe0 = static_cast<TH2F *>(referenceRunFile->Get("_Heatmap: DCe, zDC>=0"));
   TH2F *heatmapReferenceDCe1 = static_cast<TH2F *>(referenceRunFile->Get("_Heatmap: DCe, zDC<0"));
   TH2F *heatmapReferenceDCw0 = static_cast<TH2F *>(referenceRunFile->Get("_Heatmap: DCw, zDC>=0"));
   TH2F *heatmapReferenceDCw1 = static_cast<TH2F *>(referenceRunFile->Get("_Heatmap: DCw, zDC<0"));

   for (int i = 1; i <= heatmapReferenceDCe0->GetXaxis()->GetNbins(); i++)
   {
      for (int j = 1; j <= heatmapReferenceDCe0->GetYaxis()->GetNbins(); j++)
      {
         if (dmCutter.IsDeadDC(0, 1., heatmapReferenceDCe0->GetXaxis()->GetBinCenter(i),
                               heatmapReferenceDCe0->GetYaxis()->GetBinCenter(j)))
         {
            heatmapReferenceDCe0->SetBinContent(i, j, 0.);
         }
      }
   }

   for (int i = 1; i <= heatmapReferenceDCe1->GetXaxis()->GetNbins(); i++)
   {
      for (int j = 1; j <= heatmapReferenceDCe1->GetYaxis()->GetNbins(); j++)
      {
         if (dmCutter.IsDeadDC(0, -1., heatmapReferenceDCe1->GetXaxis()->GetBinCenter(i),
                               heatmapReferenceDCe1->GetYaxis()->GetBinCenter(j)))
         {
            heatmapReferenceDCe1->SetBinContent(i, j, 0.);
         }
      }
   }

   for (int i = 1; i <= heatmapReferenceDCw0->GetXaxis()->GetNbins(); i++)
   {
      for (int j = 1; j <= heatmapReferenceDCw0->GetYaxis()->GetNbins(); j++)
      {
         if (dmCutter.IsDeadDC(1, 1., heatmapReferenceDCw0->GetXaxis()->GetBinCenter(i),
                               heatmapReferenceDCw0->GetYaxis()->GetBinCenter(j)))
         {
            heatmapReferenceDCw0->SetBinContent(i, j, 0.);
         }
      }
   }

   for (int i = 1; i <= heatmapReferenceDCw1->GetXaxis()->GetNbins(); i++)
   {
      for (int j = 1; j <= heatmapReferenceDCw1->GetYaxis()->GetNbins(); j++)
      {
         if (dmCutter.IsDeadDC(1, -1., heatmapReferenceDCw1->GetXaxis()->GetBinCenter(i),
                               heatmapReferenceDCw1->GetYaxis()->GetBinCenter(j)))
         {
            heatmapReferenceDCw1->SetBinContent(i, j, 0.);
         }
      }
   }

   TH1D *referenceDCe0Board = heatmapReferenceDCe0->
      ProjectionY("reference DCe0 board", 1, heatmapReferenceDCe0->GetXaxis()->GetNbins());
   TH1D *referenceDCe1Board = heatmapReferenceDCe1->
      ProjectionY("reference DCe1 board", 1, heatmapReferenceDCe1->GetXaxis()->GetNbins());
   TH1D *referenceDCw0Board = heatmapReferenceDCw0->
      ProjectionY("reference DCw0 board", 1, heatmapReferenceDCw0->GetXaxis()->GetNbins());
   TH1D *referenceDCw1Board = heatmapReferenceDCw1->
      ProjectionY("reference DCw1 board", 1, heatmapReferenceDCw1->GetXaxis()->GetNbins());

   const double nEvtRef = static_cast<TH1D *>(referenceRunFile->Get("centrality"))->Integral();

   referenceDCe0Board->Scale(1./nEvtRef);
   referenceDCe1Board->Scale(1./nEvtRef);
   referenceDCw0Board->Scale(1./nEvtRef);
   referenceDCw1Board->Scale(1./nEvtRef);

   // graphs for chi2/NDF
   TGraph grChi2NDFDCe0Board;
   TGraph grChi2NDFDCe1Board;
   TGraph grChi2NDFDCw0Board;
   TGraph grChi2NDFDCw1Board;

   TFile outputFile((outputDir + "/DC.root").c_str(), "RECREATE");

   std::vector<int> passedRuns;

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

      const double nEvt = static_cast<TH1D *>(inputFile->Get("centrality"))->Integral();

      projDCe0Board->Scale(1./nEvt);
      projDCe1Board->Scale(1./nEvt);
      projDCw0Board->Scale(1./nEvt);
      projDCw1Board->Scale(1./nEvt);

      projDCe0Board->Divide(referenceDCe0Board);
      projDCe1Board->Divide(referenceDCe1Board);
      projDCw0Board->Divide(referenceDCw0Board);
      projDCw1Board->Divide(referenceDCw1Board);

      outputFile.mkdir(std::to_string(run).c_str());
      outputFile.cd(std::to_string(run).c_str());

      TF1 fit("const fit", "pol0");
      
      projDCe0Board->Fit(&fit, "QP");
      const double chi2NDFDCe0Board = fit.GetChisquare()/fit.GetNDF();
      projDCe0Board->Write();

      projDCe1Board->Fit(&fit, "QP");
      const double chi2NDFDCe1Board = fit.GetChisquare()/fit.GetNDF();
      projDCe1Board->Write();

      projDCw0Board->Fit(&fit, "QP");
      const double chi2NDFDCw0Board = fit.GetChisquare()/fit.GetNDF();
      projDCw0Board->Write();

      projDCw1Board->Fit(&fit, "QP");
      const double chi2NDFDCw1Board = fit.GetChisquare()/fit.GetNDF();
      projDCw1Board->Write();

      grChi2NDFDCe0Board.AddPoint(run, chi2NDFDCe0Board);
      grChi2NDFDCe1Board.AddPoint(run, chi2NDFDCe1Board);
      grChi2NDFDCw0Board.AddPoint(run, chi2NDFDCw0Board);
      grChi2NDFDCw1Board.AddPoint(run, chi2NDFDCw1Board);

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

   CppTools::PrintInfo(std::to_string(passedRuns.size()) + " runs out of " + 
                       std::to_string(goodRuns.size()) + " passed DC board check");

   goodRuns = passedRuns;

   grChi2NDFDCe0Board.SetMarkerStyle(7);
   grChi2NDFDCe1Board.SetMarkerStyle(7);
   grChi2NDFDCw0Board.SetMarkerStyle(7);
   grChi2NDFDCw1Board.SetMarkerStyle(7);

   TCanvas canv("canv", "", 1000, 600);

   gPad->SetRightMargin(0.06); gPad->SetTopMargin(0.02); 
   gPad->SetLeftMargin(0.11); gPad->SetBottomMargin(0.1);

   ROOTTools::DrawFrame(grChi2NDFDCe0Board.GetXaxis()->GetBinLowEdge(1), 
                        grChi2NDFDCe0Board.GetYaxis()->GetBinLowEdge(1),
                        grChi2NDFDCe0Board.GetXaxis()->
                        GetBinUpEdge(grChi2NDFDCe0Board.GetXaxis()->GetNbins()),
                        grChi2NDFDCe0Board.GetYaxis()->
                        GetBinUpEdge(grChi2NDFDCe0Board.GetYaxis()->GetNbins()), 
                        "", "run index", "#chi^{2}/NDF", 0.9, 1.2);

   grChi2NDFDCe0Board.Draw("P");

   ROOTTools::PrintCanvas(&canv, outputDir + "/Chi2NDFDCe0");

   canv.Clear();

   ROOTTools::DrawFrame(grChi2NDFDCe1Board.GetXaxis()->GetBinLowEdge(1), 
                        grChi2NDFDCe1Board.GetYaxis()->GetBinLowEdge(1),
                        grChi2NDFDCe1Board.GetXaxis()->
                        GetBinUpEdge(grChi2NDFDCe1Board.GetXaxis()->GetNbins()),
                        grChi2NDFDCe1Board.GetYaxis()->
                        GetBinUpEdge(grChi2NDFDCe1Board.GetYaxis()->GetNbins()), 
                        "", "run index", "#chi^{2}/NDF", 0.9, 1.2);

   grChi2NDFDCe1Board.Draw("P");

   ROOTTools::PrintCanvas(&canv, outputDir + "/Chi2NDFDCe1");

   canv.Clear();

   ROOTTools::DrawFrame(grChi2NDFDCw0Board.GetXaxis()->GetBinLowEdge(1), 
                        grChi2NDFDCw0Board.GetYaxis()->GetBinLowEdge(1),
                        grChi2NDFDCw0Board.GetXaxis()->
                        GetBinUpEdge(grChi2NDFDCw0Board.GetXaxis()->GetNbins()),
                        grChi2NDFDCw0Board.GetYaxis()->
                        GetBinUpEdge(grChi2NDFDCw0Board.GetYaxis()->GetNbins()), 
                        "", "run index", "#chi^{2}/NDF", 0.9, 1.2);

   grChi2NDFDCw0Board.Draw("P");

   ROOTTools::PrintCanvas(&canv, outputDir + "/Chi2NDFDCw0");

   canv.Clear();

   ROOTTools::DrawFrame(grChi2NDFDCw1Board.GetXaxis()->GetBinLowEdge(1), 
                        grChi2NDFDCw1Board.GetYaxis()->GetBinLowEdge(1),
                        grChi2NDFDCw1Board.GetXaxis()->
                        GetBinUpEdge(grChi2NDFDCw1Board.GetXaxis()->GetNbins()),
                        grChi2NDFDCw1Board.GetYaxis()->
                        GetBinUpEdge(grChi2NDFDCw1Board.GetYaxis()->GetNbins()), 
                        "", "run index", "#chi^{2}/NDF", 0.9, 1.2);

   grChi2NDFDCw1Board.Draw("P");

   ROOTTools::PrintCanvas(&canv, outputDir + "/Chi2NDFDCw1");
}

#endif /* FIND_GOOD_RUNS_CPP */
