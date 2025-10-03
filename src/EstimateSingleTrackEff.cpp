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

// this namespace is only used so that documentation does not become a mess
// so there is no need to enforce the contents inside of it 
// being accessed only via the scope resolution operator in this file
using namespace EstimateSingleTrackEff;

int main(int argc, char **argv)
{
   if (argc != 2) 
   {
      std::string errMsg = "Expected 1 parameter while " + std::to_string(argc - 1) + " ";
      errMsg += "parameter(s) were provided \n Usage: bin/EstimateSingleTrackEff ";
      errMsg += "inputYAMLName/path to .yaml input";
      CppTools::PrintError(errMsg);
   }
 
   inputYAMLSim.OpenFile(argv[1], "single_track_sim");
   inputYAMLSim.CheckStatus("single_track_sim");
 
   runName = inputYAMLSim["run_name"].as<std::string>();

   inputYAMLMain.OpenFile("input/" + runName + "/main.yaml");
   inputYAMLMain.CheckStatus("main");

   pTMin = inputYAMLSim["pt_min"].as<double>();
   pTMax = inputYAMLSim["pt_max"].as<double>();

   gStyle->SetOptStat(0);
   gErrorIgnoreLevel = kWarning;

   outputDir = "output/SingleTrackEff/" + runName;
   void(system(("mkdir -p " + outputDir).c_str()));

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

   text.SetTextFont(43);
   text.SetTextSize(60);

   inputDataFilePiPlus = TFile::Open(inputDataFileNamePiPlus.c_str());
   inputDataFileKPlus = TFile::Open(inputDataFileNameKPlus.c_str());
   inputDataFileP = TFile::Open(inputDataFileNameP.c_str());
   inputDataFilePiMinus = TFile::Open(inputDataFileNamePiMinus.c_str());
   inputDataFileKMinus = TFile::Open(inputDataFileNameKMinus.c_str());
   inputDataFilePBar = TFile::Open(inputDataFileNamePBar.c_str());

   distrOrigPTPiPlus = static_cast<TH1F *>(inputDataFilePiPlus->Get("orig pT"));
   distrOrigPTKPlus = static_cast<TH1F *>(inputDataFileKPlus->Get("orig pT"));
   distrOrigPTP = static_cast<TH1F *>(inputDataFileP->Get("orig pT"));
   distrOrigPTPiMinus = static_cast<TH1F *>(inputDataFilePiMinus->Get("orig pT"));
   distrOrigPTKMinus = static_cast<TH1F *>(inputDataFileKMinus->Get("orig pT"));
   distrOrigPTPBar = static_cast<TH1F *>(inputDataFilePBar->Get("orig pT"));

   EstimateEffForSingleDetector("DC-PC1");
   EstimateEffForSingleDetector("PC2");
   EstimateEffForSingleDetector("PC3");
   EstimateEffForSingleDetector("TOFe");
   EstimateEffForSingleDetector("TOFw");
   EstimateEffForSingleDetector("EMCale0");
   EstimateEffForSingleDetector("EMCale1");
   EstimateEffForSingleDetector("EMCale2");
   EstimateEffForSingleDetector("EMCale3");
   EstimateEffForSingleDetector("EMCalw0");
   EstimateEffForSingleDetector("EMCalw1");
   EstimateEffForSingleDetector("EMCalw2");
   EstimateEffForSingleDetector("EMCalw3");
   EstimateEffForSingleDetector("TOFe", true);
   EstimateEffForSingleDetector("TOFw", true);
   EstimateEffForSingleDetector("EMCale2", true);
   EstimateEffForSingleDetector("EMCale3", true);
   EstimateEffForSingleDetector("EMCalw0", true);
   EstimateEffForSingleDetector("EMCalw1", true);
   EstimateEffForSingleDetector("EMCalw2", true);
   EstimateEffForSingleDetector("EMCalw3", true);

   return 1;
}

void EstimateSingleTrackEff::EstimateEffForSingleDetector(const std::string& detectorName, 
                                                          const bool isIdentification)
{
   std::string nameAux = "";
   if (isIdentification) nameAux = "id ";

   TH1F* distrRecPTPiPlus = static_cast<TH1F *>
      (inputDataFilePiPlus->Get(("rec " + nameAux + "pT: " + detectorName).c_str()));
   TH1F* distrRecPTKPlus = static_cast<TH1F *>
      (inputDataFileKPlus->Get(("rec " + nameAux + "pT: " + detectorName).c_str()));
   TH1F* distrRecPTP = static_cast<TH1F *>
      (inputDataFileP->Get(("rec " + nameAux + "pT: " + detectorName).c_str()));
   TH1F* distrRecPTPiMinus = static_cast<TH1F *>
      (inputDataFilePiMinus->Get(("rec " + nameAux + "pT: " + detectorName).c_str()));
   TH1F* distrRecPTKMinus = static_cast<TH1F *>
      (inputDataFileKMinus->Get(("rec " + nameAux + "pT: " + detectorName).c_str()));
   TH1F* distrRecPTPBar = static_cast<TH1F *>
      (inputDataFilePBar->Get(("rec " + nameAux + "pT: " + detectorName).c_str()));

   distrRecPTPiPlus->Divide(distrOrigPTPiPlus);
   distrRecPTKPlus->Divide(distrOrigPTKPlus);
   distrRecPTP->Divide(distrOrigPTP);
   distrRecPTPiMinus->Divide(distrOrigPTPiMinus);
   distrRecPTKMinus->Divide(distrOrigPTKMinus);
   distrRecPTPBar->Divide(distrOrigPTPBar);

   distrRecPTPiPlus->SetLineColor(kP6Violet);
   distrRecPTKPlus->SetLineColor(kP6Gray);
   distrRecPTP->SetLineColor(kP6Grape);
   distrRecPTPiMinus->SetLineColor(kP6Red);
   distrRecPTKMinus->SetLineColor(kP6Yellow);
   distrRecPTPBar->SetLineColor(kP6Blue);

   distrRecPTPiPlus->SetMarkerColor(kP6Violet);
   distrRecPTKPlus->SetMarkerColor(kP6Gray);
   distrRecPTP->SetMarkerColor(kP6Grape);
   distrRecPTPiMinus->SetMarkerColor(kP6Red);
   distrRecPTKMinus->SetMarkerColor(kP6Yellow);
   distrRecPTPBar->SetMarkerColor(kP6Blue);

   distrRecPTPiPlus->SetLineWidth(8);
   distrRecPTKPlus->SetLineWidth(7);
   distrRecPTP->SetLineWidth(6);
   distrRecPTPiMinus->SetLineWidth(5);
   distrRecPTKMinus->SetLineWidth(4);
   distrRecPTPBar->SetLineWidth(3);

   TLegend legend(0.65, 0.15, 0.95, 0.3);
   legend.SetLineColorAlpha(0, 0.);
   legend.SetFillColorAlpha(0, 0.);
   legend.SetNColumns(3);

   TLine line(pTMin, 1., pTMax, 1.);
   line.SetLineWidth(4);
   line.SetLineStyle(2);
   line.SetLineColorAlpha(kBlack, 0.5);

   TCanvas canv(detectorName.c_str(), "", 800, 800);

   gPad->SetLogy();
   gPad->SetLeftMargin(0.14);
   gPad->SetBottomMargin(0.115);
   gPad->SetRightMargin(0.03);
   gPad->SetTopMargin(0.02);

   // first bin in pT range (pTMin, pTMax)
   const int firstBinInRange = distrRecPTP->GetXaxis()->FindBin(pTMin + 1e-6);

   const double yMin = 
      CppTools::Maximum(CppTools::Minimum(distrRecPTPiPlus->GetBinContent(firstBinInRange),
                                          distrRecPTKPlus->GetBinContent(firstBinInRange),
                                          distrRecPTP->GetBinContent(firstBinInRange),
                                          distrRecPTPiMinus->GetBinContent(firstBinInRange),
                                          distrRecPTPiMinus->GetBinContent(firstBinInRange),
                                          distrRecPTPBar->GetBinContent(firstBinInRange))/20., 
                        1e-8)/2.;

   ROOTTools::DrawFrame(pTMin, yMin, pTMax, 2., "", "p_{T} [GeV/c]", "#varepsilon_{" + 
                        static_cast<std::string>((isIdentification) ? "id" : "reg") + "}");

   line.Draw();

   distrRecPTPiPlus->Draw("SAME");
   distrRecPTKPlus->Draw("SAME");
   distrRecPTP->Draw("SAME");
   distrRecPTPiMinus->Draw("SAME");
   distrRecPTKMinus->Draw("SAME");
   distrRecPTPBar->Draw("SAME");

   legend.AddEntry(distrRecPTPiPlus, "#pi^{+}");
   legend.AddEntry(distrRecPTKPlus, "K^{+}");
   legend.AddEntry(distrRecPTP, "p");
   legend.AddEntry(distrRecPTPiMinus, "#pi^{-}");
   legend.AddEntry(distrRecPTKMinus, "K^{-}");
   legend.AddEntry(distrRecPTPBar, "#bar{p}");

   text.DrawTextNDC(0.3, 0.15, detectorName.c_str());

   legend.Draw();

   ROOTTools::PrintCanvas(&canv, outputDir + "/" + detectorName + 
                          ((isIdentification) ? "IdEff" : "RecEff"));
}

#endif /* ESTIMATE_SINGLE_TRACK_EFF_CPP */
