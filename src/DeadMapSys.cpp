/** 
 *  @file   DeadMapSys.cpp 
 *  @brief  Contains realisations of functions that are used for estimation of systematic uncertainties of heatmaps (real vs simulated) 
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysis).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/

#ifndef DEAD_MAP_SYS_CPP
#define DEAD_MAP_SYS_CPP

#include "../include/DeadMapSys.hpp"

// this namespace is only used so that documentation will not become a mess
// so there is no need to enforce the contents inside of it 
// being accessed only via the scope resolution operator in this file
using namespace DeadMapSys;

void DeadMapSys::CheckHists(const TH2F *histReal, const TH2F *histSim, const std::string& name)
{
   if (!histReal) 
   {
      CppTools::PrintError("Real histogram named " + name + " does not exists");
   }
   if (!histSim) 
   {
      CppTools::PrintError("Simulated histogram named" + name + " does not exists");
   }

   if (histReal->GetXaxis()->GetNbins() != histSim->GetXaxis()->GetNbins() ||
       fabs(histReal->GetXaxis()->GetBinLowEdge(1) - 
            histReal->GetXaxis()->GetBinLowEdge(1)) > 1e-6 ||
       fabs(histReal->GetXaxis()->GetBinUpEdge(histReal->GetXaxis()->GetNbins()) - 
            histReal->GetXaxis()->GetBinUpEdge(histSim->GetXaxis()->GetNbins())) > 1e-6) 
   {
      CppTools::PrintError("X axis bins are inconsistent for histograms named " + name);
   }
   if (histReal->GetYaxis()->GetNbins() != histSim->GetYaxis()->GetNbins() ||
       fabs(histReal->GetYaxis()->GetBinLowEdge(1) - 
            histReal->GetYaxis()->GetBinLowEdge(1)) > 1e-6 ||
       fabs(histReal->GetYaxis()->GetBinUpEdge(histReal->GetYaxis()->GetNbins()) - 
            histReal->GetYaxis()->GetBinUpEdge(histSim->GetYaxis()->GetNbins())) > 1e-6) 
   {
      CppTools::PrintError("Y axis bins are inconsistent for histograms named " + name);
   }
}

double DeadMapSys::GetNormRatio(const double ratio)
{
   if (ratio > 1.) return ratio;
   return 1./ratio;
}

double DeadMapSys::GetUncertainty(TH2F *&realDistr, TH2F *&simDistr, 
                                  TH2F *&realCutDistr, TH2F *&simCutDistr, 
                                  const int numberOfHeatmapDivisions,
                                  const std::string& detectorName, const std::string& title, 
                                  const std::string& xTitle, const std::string& yTitle,
                                  const int rebinX, const int rebinY, const bool drawYProj)
{	
	TH1D *realCutDistrProjX = realCutDistr->
      ProjectionX((title + "_real_x").c_str(), 1, realCutDistr->GetYaxis()->GetNbins());
	TH1D *simCutDistrProjX = simCutDistr->
      ProjectionX((title + "_sim_x").c_str(), 1, simCutDistr->GetYaxis()->GetNbins());

   realCutDistr->SetMinimum(0.5);
   simCutDistr->SetMinimum(1e-15);

   realCutDistrProjX->RebinX(rebinX);
   simCutDistrProjX->RebinX(rebinX);

	realCutDistrProjX->SetFillColorAlpha(kOrange - 4, 0.5);
	
	realCutDistrProjX->Scale(1./realCutDistrProjX->Integral(1, 
                            realCutDistrProjX->GetXaxis()->GetNbins()), "nosw2");
	simCutDistrProjX->Scale(1./simCutDistrProjX->Integral(1, 
                           simCutDistrProjX->GetXaxis()->GetNbins()), "nosw2");

	TLegend projXLegend = TLegend(0.4, 0.78, 0.9, 0.88);
	
	projXLegend.SetNColumns(2);
	
	projXLegend.SetLineColorAlpha(0, 0);
	projXLegend.SetFillColorAlpha(0, 0);
	
	realCutDistrProjX->SetLineColor(kRed - 3);
	simCutDistrProjX->SetLineColor(kAzure - 3);

	projXLegend.AddEntry(realCutDistrProjX, "Data");
	projXLegend.AddEntry(simCutDistrProjX, "MC");
	
	realCutDistrProjX->SetTitle((title + " X projections data vs MC").c_str());
	realCutDistrProjX->GetXaxis()->SetTitle(xTitle.c_str());

	realCutDistrProjX->SetMaximum(CppTools::Maximum(realCutDistrProjX->GetMaximum()*1.3, 
                                                   simCutDistrProjX->GetMaximum()*1.3));

	TCanvas canv = TCanvas("canv", "canv", 900, 450);

	canv.Divide(2);

	canv.cd(1);
   gPad->SetRightMargin(0.14); gPad->SetTopMargin(0.07); 
   gPad->SetLeftMargin(0.125); gPad->SetBottomMargin(0.105);

   ROOTTools::DrawFrame(realDistr, title, xTitle, yTitle, 
                        1., 1.2, 0.05, 0.05, true, true, "COLZ");
	
	canv.cd(2);
   gPad->SetRightMargin(0.14); gPad->SetTopMargin(0.07); 
   gPad->SetLeftMargin(0.125); gPad->SetBottomMargin(0.105);

   ROOTTools::DrawFrame(realCutDistr, "Cut " + title, xTitle, yTitle, 
                        1., 1.2, 0.05, 0.05, true, true, "COLZ");

   ROOTTools::PrintCanvas(&canv, "output/Deadmaps/" + runName + "/" + detectorName);

   if (drawYProj)
   {
      TH1D *realCutDistrProjY = realCutDistr->
         ProjectionY((title + "_real_y").c_str(), 1, realCutDistr->GetXaxis()->GetNbins());
      TH1D *simCutDistrProjY = simCutDistr->
         ProjectionY((title + "_sim_y").c_str(), 1, simCutDistr->GetXaxis()->GetNbins());

      realCutDistrProjY->RebinX(rebinY);
      simCutDistrProjY->RebinX(rebinY);

      realCutDistrProjY->SetFillColorAlpha(kOrange - 4, 0.5);
      
      realCutDistrProjY->Scale(1./realCutDistrProjY->Integral(1, 
                               realCutDistrProjY->GetXaxis()->GetNbins()), "nosw2");
      simCutDistrProjY->Scale(1./simCutDistrProjY->Integral(1, 
                              simCutDistrProjY->GetXaxis()->GetNbins()), "nosw2");

      TLegend projYLegend = TLegend(0.4, 0.78, 0.9, 0.88);

      projYLegend.SetNColumns(2);

      projYLegend.SetLineColorAlpha(0, 0);
      projYLegend.SetFillColorAlpha(0, 0);

      projYLegend.AddEntry(realCutDistrProjY, "data");
      projYLegend.AddEntry(simCutDistrProjY, "MC");

      realCutDistrProjY->SetTitle((title + " Y projections data vs MC").c_str());
      realCutDistrProjY->GetXaxis()->SetTitle(yTitle.c_str());

      realCutDistrProjY->SetMaximum(CppTools::Maximum(realCutDistrProjY->GetMaximum()*1.3, 
                                                      simCutDistrProjY->GetMaximum()*1.3));

      realCutDistrProjY->GetXaxis()->SetLabelSize(0.05);
      realCutDistrProjY->GetYaxis()->SetLabelSize(0.05);

      realCutDistrProjY->SetLineColor(kRed-3);
      simCutDistrProjY->SetLineColor(kAzure-3);

      TCanvas projCanv = TCanvas("projCanv", "canv", 1800, 900);
      
      projCanv.cd();
      gPad->Divide(2, 2);
      
      projCanv.cd(1);
      gPad->SetRightMargin(0.1); gPad->SetTopMargin(0.07); 
      gPad->SetLeftMargin(0.07); gPad->SetBottomMargin(0.1);

      ROOTTools::DrawFrame(realCutDistr, "Cut " + title, xTitle, yTitle, 
                           0.9, 0.7, 0.05, 0.05, true, true, "COLZ");
      
      projCanv.cd(3);
      gPad->SetRightMargin(0.1); gPad->SetTopMargin(0.07); 
      gPad->SetLeftMargin(0.07); gPad->SetBottomMargin(0.1);

      ROOTTools::DrawFrame(simCutDistr, "MC cut " + title, xTitle, yTitle, 
                           0.9, 0.7, 0.05, 0.05, true, true, "COLZ");
  
      projCanv.cd(2);
      gPad->SetRightMargin(0.02); gPad->SetTopMargin(0.07); 
      gPad->SetLeftMargin(0.1); gPad->SetBottomMargin(0.1);

      ROOTTools::DrawFrame(realCutDistrProjX, "X projections of cut " + title, xTitle, "Counts", 
                           0.9, 1.0, 0.05, 0.05, true, true);

      simCutDistrProjX->Draw("SAME HIST");
      projXLegend.Draw();

      projCanv.cd(4);
      gPad->SetRightMargin(0.02); gPad->SetTopMargin(0.07); 
      gPad->SetLeftMargin(0.1); gPad->SetBottomMargin(0.1);

      ROOTTools::DrawFrame(realCutDistrProjY, "Y projections of cut " + title, yTitle, "Counts", 
                           0.9, 1.0, 0.05, 0.05, true, true);
      simCutDistrProjY->Draw("SAME HIST");
      projYLegend.Draw();

      ROOTTools::PrintCanvas(&projCanv, "output/Systematics/" + runName + "/" + detectorName);
   }
   else
   {
      TCanvas projCanv = TCanvas("projCanv", "canv", 1800, 600);
      
      projCanv.cd();
      gPad->Divide(3, 1);
      
      projCanv.cd(1);
      gPad->SetRightMargin(0.145); gPad->SetTopMargin(0.07); 
      gPad->SetLeftMargin(0.09); gPad->SetBottomMargin(0.09);

      ROOTTools::DrawFrame(realCutDistr, "Cut " + title, xTitle, yTitle, 
                           0.8, 1.0, 0.05, 0.05, true, true, "COLZ");
      
      projCanv.cd(2);
      gPad->SetRightMargin(0.145); gPad->SetTopMargin(0.07); 
      gPad->SetLeftMargin(0.09); gPad->SetBottomMargin(0.09);
      simCutDistr->Draw("colz");

      ROOTTools::DrawFrame(simCutDistr, "Cut MC " + title, xTitle, yTitle, 
                           0.8, 1.0, 0.05, 0.05, true, true, "COLZ");
      
      projCanv.cd(3);
      gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.07); 
      gPad->SetLeftMargin(0.145); gPad->SetBottomMargin(0.09);

      ROOTTools::DrawFrame(realCutDistrProjX, "X projections of cut " + title, xTitle, "Counts", 
                           0.8, 1.65, 0.05, 0.05, true, true);
      simCutDistrProjX->Draw("SAME HIST");
      projXLegend.Draw();

      ROOTTools::PrintCanvas(&projCanv, "output/Systematics/" + runName + "/" + detectorName);
   }

   const double realDataLost = (1. - realCutDistr->Integral()/realDistr->Integral())*100.;
   const double simDataLost = (1. - simCutDistr->Integral()/simDistr->Integral())*100.;

   realCutDistr->Scale(1./realCutDistr->Integral(1, realCutDistr->GetXaxis()->GetNbins(),
                                                 1, realCutDistr->GetYaxis()->GetNbins()));
   simCutDistr->Scale(1./simCutDistr->Integral(1, simCutDistr->GetXaxis()->GetNbins(),
                                               1, simCutDistr->GetYaxis()->GetNbins()));

   std::vector<double> divisionRealIntegral;
   std::vector<double> divisionSimIntegral;

   divisionRealIntegral.resize(numberOfHeatmapDivisions);
   divisionSimIntegral.resize(numberOfHeatmapDivisions);

   // bins for each division are chosen uniformly
   for (int i = 1; i <= realCutDistr->GetXaxis()->GetNbins(); i++)
   {
      for (int j = 1; j <= realCutDistr->GetYaxis()->GetNbins(); j++)
      {
         divisionRealIntegral[realCutDistr->GetBin(i, j) % numberOfHeatmapDivisions] += 
            realCutDistr->GetBinContent(i, j);
         divisionSimIntegral[realCutDistr->GetBin(i, j) % numberOfHeatmapDivisions] += 
            simCutDistr->GetBinContent(i, j);
      }
   }

   // relative uncertainty
   double uncertainty = 0.;

   for (int i = 0; i < numberOfHeatmapDivisions; i++)
   {
      if (divisionRealIntegral[i] < 1e-3 || divisionSimIntegral[i] < 1e-3) continue;

      const double ratio = divisionRealIntegral[i]/divisionSimIntegral[i];
      uncertainty += (1. - ratio)*(1. - ratio);
   }

   uncertainty = sqrt(uncertainty/static_cast<double>(numberOfHeatmapDivisions));

   table.PrintRow(detectorName, CppTools::DtoStr(uncertainty*100.) + " %", 
                  CppTools::DtoStr(realDataLost, 3) + "%", 
                  CppTools::DtoStr(simDataLost, 3) + "%");

	return uncertainty;
}

int main(int argc, char **argv)
{
   if (argc < 2 || argc > 3) 
   {
      std::string errMsg = "Expected 1 parameter while " + std::to_string(argc - 1) + " ";
      errMsg += "parameter(s) were provided \n Usage: bin/DeadMapSys ";
      errMsg += "inputYAMLName numberOfThreads=std::thread::hardware_concurrency()";
      CppTools::PrintError(errMsg);
   }
 
   CppTools::CheckInputFile(argv[1]);

   inputYAMLMain.OpenFile(argv[1], "main");
   inputYAMLMain.CheckStatus("main");

   runName = inputYAMLMain["run_name"].as<std::string>();

   gErrorIgnoreLevel = kWarning;
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);

   outputDirDM = "output/Deadmaps/" + runName + "/";
   outputDirSys = "output/Systematics/" + runName + "/";
   outputDirParameters = "data/Parameters/Systematics/" + runName + "/";

   std::filesystem::create_directories(outputDirDM);
   std::filesystem::create_directories(outputDirSys);
   std::filesystem::create_directories(outputDirParameters);

   const std::string inputRealDataFileName = "data/Real/" + runName + "/SingleTrack/sum.root";
   const std::string inputSimDataFileName = "data/PostSim/" + runName + "/SingleTrack/all.root";

   CppTools::CheckInputFile(inputRealDataFileName);
   CppTools::CheckInputFile(inputSimDataFileName);

   inputRealDataFile = TFile::Open(inputRealDataFileName.c_str());
   inputSimDataFile = TFile::Open(inputSimDataFileName.c_str());

   std::ofstream systematicsOutputFile(outputDirParameters + "Acceptance.txt");

   const std::string detectorsConfiguration = 
      inputYAMLMain["detectors_configuration"].as<std::string>();

   dmCutter.Initialize(runName, detectorsConfiguration);

	table.Begin("Heatmap acceptance info");
	table.PrintHeader("detector", "uncertainty", "data lost", "MC lost");

   if (detectorsConfiguration[0] == '1') // DC
   {
      TH2F *realHeatmapDCe0 = static_cast<TH2F *>(inputRealDataFile->Get("_Heatmap: DCe, zDC>=0"));
      TH2F *realHeatmapDCe1 = static_cast<TH2F *>(inputRealDataFile->Get("_Heatmap: DCe, zDC<0"));
      TH2F *realHeatmapDCw0 = static_cast<TH2F *>(inputRealDataFile->Get("_Heatmap: DCw, zDC>=0"));
      TH2F *realHeatmapDCw1 = static_cast<TH2F *>(inputRealDataFile->Get("_Heatmap: DCw, zDC<0"));

      TH2F *simHeatmapDCe0 = static_cast<TH2F *>(inputSimDataFile->Get("_Heatmap: DCe, zDC>=0"));
      TH2F *simHeatmapDCe1 = static_cast<TH2F *>(inputSimDataFile->Get("_Heatmap: DCe, zDC<0"));
      TH2F *simHeatmapDCw0 = static_cast<TH2F *>(inputSimDataFile->Get("_Heatmap: DCw, zDC>=0"));
      TH2F *simHeatmapDCw1 = static_cast<TH2F *>(inputSimDataFile->Get("_Heatmap: DCw, zDC<0"));

      CheckHists(realHeatmapDCe0, simHeatmapDCe0, "_Heatmap: DCe, zDC>=0");
      CheckHists(realHeatmapDCe1, simHeatmapDCe1, "_Heatmap: DCe, zDC<0");
      CheckHists(realHeatmapDCw0, simHeatmapDCw0, "_Heatmap: DCw, zDC>=0");
      CheckHists(realHeatmapDCw1, simHeatmapDCw1, "_Heatmap: DCw, zDC<0");

      TH2F *realCutHeatmapDCe0 = static_cast<TH2F *>(realHeatmapDCe0->Clone("DCe, zDC>=0 real"));
      TH2F *realCutHeatmapDCe1 = static_cast<TH2F *>(realHeatmapDCe1->Clone("DCe, zDC<0 real"));
      TH2F *realCutHeatmapDCw0 = static_cast<TH2F *>(realHeatmapDCw0->Clone("DCw, zDC>=0 real"));
      TH2F *realCutHeatmapDCw1 = static_cast<TH2F *>(realHeatmapDCw1->Clone("DCw, zDC<0 real"));

      TH2F *simCutHeatmapDCe0 = static_cast<TH2F *>(simHeatmapDCe0->Clone("DCe, zDC>=0 sim"));
      TH2F *simCutHeatmapDCe1 = static_cast<TH2F *>(simHeatmapDCe1->Clone("DCe, zDC<0 sim"));
      TH2F *simCutHeatmapDCw0 = static_cast<TH2F *>(simHeatmapDCw0->Clone("DCw, zDC>=0 sim"));
      TH2F *simCutHeatmapDCw1 = static_cast<TH2F *>(simHeatmapDCw1->Clone("DCw, zDC<0 sim"));

      for (int i = 1; i <= realHeatmapDCe0->GetXaxis()->GetNbins(); i++)
      {
         for (int j = 1; j <= realHeatmapDCe0->GetYaxis()->GetNbins(); j++)
         {
            if (dmCutter.IsDeadDC(0, 1., realHeatmapDCe0->GetXaxis()->GetBinCenter(i),
                                  realHeatmapDCe0->GetYaxis()->GetBinCenter(j)))
            {
               realCutHeatmapDCe0->SetBinContent(i, j, 0.);
               simCutHeatmapDCe0->SetBinContent(i, j, 0.);
            }
         }
      }
      for (int i = 1; i <= realHeatmapDCe1->GetXaxis()->GetNbins(); i++)
      {
         for (int j = 1; j <= realHeatmapDCe1->GetYaxis()->GetNbins(); j++)
         {
            if (dmCutter.IsDeadDC(0, -1., realHeatmapDCe1->GetXaxis()->GetBinCenter(i),
                                  realHeatmapDCe1->GetYaxis()->GetBinCenter(j)))
            {
               realCutHeatmapDCe1->SetBinContent(i, j, 0.);
               simCutHeatmapDCe1->SetBinContent(i, j, 0.);
            }
         }
      }
      for (int i = 1; i <= realHeatmapDCw0->GetXaxis()->GetNbins(); i++)
      {
         for (int j = 1; j <= realHeatmapDCw0->GetYaxis()->GetNbins(); j++)
         {
            if (dmCutter.IsDeadDC(1, 1., realHeatmapDCw0->GetXaxis()->GetBinCenter(i),
                                  realHeatmapDCw0->GetYaxis()->GetBinCenter(j)))
            {
               realCutHeatmapDCw0->SetBinContent(i, j, 0.);
               simCutHeatmapDCw0->SetBinContent(i, j, 0.);
            }
         }
      }
      for (int i = 1; i <= realHeatmapDCw1->GetXaxis()->GetNbins(); i++)
      {
         for (int j = 1; j <= realHeatmapDCw1->GetYaxis()->GetNbins(); j++)
         {
            if (dmCutter.IsDeadDC(1, -1., realHeatmapDCw1->GetXaxis()->GetBinCenter(i),
                                  realHeatmapDCw1->GetYaxis()->GetBinCenter(j)))
            {
               realCutHeatmapDCw1->SetBinContent(i, j, 0.);
               simCutHeatmapDCw1->SetBinContent(i, j, 0.);
            }
         }
      }

      systematicsOutputFile << 
         GetUncertainty(realHeatmapDCe0, simHeatmapDCe0, 
                        realCutHeatmapDCe0, simCutHeatmapDCe0, 10,
                        "DCe0", "DC east, zDC>=0", "board", "#alpha", 3, 1, false) << " " <<
         GetUncertainty(realHeatmapDCe1, simHeatmapDCe1,
                        realCutHeatmapDCe1, simCutHeatmapDCe1, 10,
                        "DCe1", "DC east, zDC<0", "board", "#alpha", 3, 1, false) << " " <<
         GetUncertainty(realHeatmapDCw0, simHeatmapDCw0,
                        realCutHeatmapDCw0, simCutHeatmapDCw0, 10,
                        "DCw0", "DC west, zDC>=0", "board", "#alpha", 3, 1, false) << " " <<
         GetUncertainty(realHeatmapDCw1, simHeatmapDCw1,
                        realCutHeatmapDCw1, simCutHeatmapDCw1, 10,
                        "DCw1", "DC west, zDC<0", "board", "#alpha", 3, 1, false) << std::endl;
   }
   else
   {
      systematicsOutputFile << 0 << " " << 0 << " " << 0 << " " << 0 << std::endl;
   }

   if (detectorsConfiguration[1] == '1') // PC1
   {
      TH2F *realHeatmapPC1e = static_cast<TH2F *>(inputRealDataFile->Get("Heatmap: PC1e"));
      TH2F *realHeatmapPC1w = static_cast<TH2F *>(inputRealDataFile->Get("Heatmap: PC1w"));

      TH2F *simHeatmapPC1e = static_cast<TH2F *>(inputSimDataFile->Get("Heatmap: PC1e"));
      TH2F *simHeatmapPC1w = static_cast<TH2F *>(inputSimDataFile->Get("Heatmap: PC1w"));

      CheckHists(realHeatmapPC1e, simHeatmapPC1e, "Heatmap: PC1e");
      CheckHists(realHeatmapPC1w, simHeatmapPC1w, "Heatmap: PC1w");

      TH2F *realCutHeatmapPC1e = static_cast<TH2F *>(realHeatmapPC1e->Clone());
      TH2F *realCutHeatmapPC1w = static_cast<TH2F *>(realHeatmapPC1w->Clone());

      TH2F *simCutHeatmapPC1e = static_cast<TH2F *>(simHeatmapPC1e->Clone());
      TH2F *simCutHeatmapPC1w = static_cast<TH2F *>(simHeatmapPC1w->Clone());

      for (int i = 1; i <= realHeatmapPC1e->GetXaxis()->GetNbins(); i++)
      {
         for (int j = 1; j <= realHeatmapPC1e->GetYaxis()->GetNbins(); j++)
         {
            if (dmCutter.IsDeadPC1(0, realHeatmapPC1e->GetXaxis()->GetBinCenter(i),
                                   realHeatmapPC1e->GetYaxis()->GetBinCenter(j)))
            {
               realCutHeatmapPC1e->SetBinContent(i, j, 0.);
               simCutHeatmapPC1e->SetBinContent(i, j, 0.);
            }
         }
      }
      for (int i = 1; i <= realHeatmapPC1w->GetXaxis()->GetNbins(); i++)
      {
         for (int j = 1; j <= realHeatmapPC1w->GetYaxis()->GetNbins(); j++)
         {
            if (dmCutter.IsDeadPC1(1, realHeatmapPC1w->GetXaxis()->GetBinCenter(i),
                                   realHeatmapPC1w->GetYaxis()->GetBinCenter(j)))
            {
               realCutHeatmapPC1w->SetBinContent(i, j, 0.);
               simCutHeatmapPC1w->SetBinContent(i, j, 0.);
            }
         }
      }

      systematicsOutputFile << 
         GetUncertainty(realHeatmapPC1e, simHeatmapPC1e, 
                        realCutHeatmapPC1e, simCutHeatmapPC1e, 10,
                        "PC1e", "PC1 east", "z_{PC1}", "y_{PC1}", 2) << " " <<
         GetUncertainty(realHeatmapPC1w, simHeatmapPC1w, 
                        realCutHeatmapPC1w, simCutHeatmapPC1w, 10,
                        "PC1w", "PC1 west", "z_{PC1}", "y_{PC1}", 2) << std::endl;
   }
   else
   {
      systematicsOutputFile << 0 << " " << 0 << std::endl;
   }

   if (detectorsConfiguration[2] == '1') // PC2
   {
      TH2F *realHeatmapPC2 = static_cast<TH2F *>(inputRealDataFile->Get("Heatmap: PC2"));

      TH2F *simHeatmapPC2 = static_cast<TH2F *>(inputSimDataFile->Get("Heatmap: PC2"));

      CheckHists(realHeatmapPC2, simHeatmapPC2, "Heatmap: PC2");

      TH2F *realCutHeatmapPC2 = static_cast<TH2F *>(realHeatmapPC2->Clone());

      TH2F *simCutHeatmapPC2 = static_cast<TH2F *>(simHeatmapPC2->Clone());

      for (int i = 1; i <= realHeatmapPC2->GetXaxis()->GetNbins(); i++)
      {
         for (int j = 1; j <= realHeatmapPC2->GetYaxis()->GetNbins(); j++)
         {
            if (dmCutter.IsDeadPC2(realHeatmapPC2->GetXaxis()->GetBinCenter(i),
                                   realHeatmapPC2->GetYaxis()->GetBinCenter(j)))
            {
               realCutHeatmapPC2->SetBinContent(i, j, 0.);
               simCutHeatmapPC2->SetBinContent(i, j, 0.);
            }
         }
      }

      systematicsOutputFile << 
         GetUncertainty(realHeatmapPC2, simHeatmapPC2, 
                        realCutHeatmapPC2, simCutHeatmapPC2, 10,
                        "PC2", "PC2", "z_{PC2}", "y_{PC2}") << std::endl;
   }
   else
   {
      systematicsOutputFile << 0 << std::endl;
   }

   if (detectorsConfiguration[3] == '1') // PC3
   {
      TH2F *realHeatmapPC3e = static_cast<TH2F *>(inputRealDataFile->Get("Heatmap: PC3e"));
      TH2F *realHeatmapPC3w = static_cast<TH2F *>(inputRealDataFile->Get("Heatmap: PC3w"));

      TH2F *simHeatmapPC3e = static_cast<TH2F *>(inputSimDataFile->Get("Heatmap: PC3e"));
      TH2F *simHeatmapPC3w = static_cast<TH2F *>(inputSimDataFile->Get("Heatmap: PC3w"));

      CheckHists(realHeatmapPC3e, simHeatmapPC3e, "Heatmap: PC3e");
      CheckHists(realHeatmapPC3w, simHeatmapPC3w, "Heatmap: PC3w");

      TH2F *realCutHeatmapPC3e = static_cast<TH2F *>(realHeatmapPC3e->Clone());
      TH2F *realCutHeatmapPC3w = static_cast<TH2F *>(realHeatmapPC3w->Clone());

      TH2F *simCutHeatmapPC3e = static_cast<TH2F *>(simHeatmapPC3e->Clone());
      TH2F *simCutHeatmapPC3w = static_cast<TH2F *>(simHeatmapPC3w->Clone());

      for (int i = 1; i <= realHeatmapPC3e->GetXaxis()->GetNbins(); i++)
      {
         for (int j = 1; j <= realHeatmapPC3e->GetYaxis()->GetNbins(); j++)
         {
            if (dmCutter.IsDeadPC3(0, realHeatmapPC3e->GetXaxis()->GetBinCenter(i),
                                   realHeatmapPC3e->GetYaxis()->GetBinCenter(j)))
            {
               realCutHeatmapPC3e->SetBinContent(i, j, 0.);
               simCutHeatmapPC3e->SetBinContent(i, j, 0.);
            }
         }
      }
      for (int i = 1; i <= realHeatmapPC3w->GetXaxis()->GetNbins(); i++)
      {
         for (int j = 1; j <= realHeatmapPC3w->GetYaxis()->GetNbins(); j++)
         {
            if (dmCutter.IsDeadPC3(1, realHeatmapPC3w->GetXaxis()->GetBinCenter(i),
                                   realHeatmapPC3w->GetYaxis()->GetBinCenter(j)))
            {
               realCutHeatmapPC3w->SetBinContent(i, j, 0.);
               simCutHeatmapPC3w->SetBinContent(i, j, 0.);
            }
         }
      }

      systematicsOutputFile << 
         GetUncertainty(realHeatmapPC3e, simHeatmapPC3e, 
                        realCutHeatmapPC3e, simCutHeatmapPC3e, 10,
                        "PC3e", "PC3 east", "z_{PC3}", "y_{PC3}", 2) << " " <<
         GetUncertainty(realHeatmapPC3w, simHeatmapPC3w, 
                        realCutHeatmapPC3w, simCutHeatmapPC3w, 10,
                        "PC3w", "PC3 west", "z_{PC3}", "y_{PC3}", 2) << std::endl;
   }
   else
   {
      systematicsOutputFile << 0 << " " << 0 << std::endl;
   }

   if (detectorsConfiguration[4] == '1') // TOFe
   {
      TH2F *realHeatmapTOFe = static_cast<TH2F *>(inputRealDataFile->Get("Heatmap: TOFe"));

      TH2F *simHeatmapTOFe = static_cast<TH2F *>(inputSimDataFile->Get("Heatmap: TOFe"));

      CheckHists(realHeatmapTOFe, simHeatmapTOFe, "Heatmap: TOFe");

      TH2F *realCutHeatmapTOFe = static_cast<TH2F *>(realHeatmapTOFe->Clone());

      TH2F *simCutHeatmapTOFe = static_cast<TH2F *>(simHeatmapTOFe->Clone());

      for (int i = 1; i <= realHeatmapTOFe->GetXaxis()->GetNbins(); i++)
      {
         for (int j = 1; j <= realHeatmapTOFe->GetYaxis()->GetNbins(); j++)
         {
            if (dmCutter.IsDeadTOFe(i - 1, j - 1))
            {
               realCutHeatmapTOFe->SetBinContent(i, j, 0.);
               simCutHeatmapTOFe->SetBinContent(i, j, 0.);
            }
         }
      }

      systematicsOutputFile << 
         GetUncertainty(realHeatmapTOFe, simHeatmapTOFe, 
                        realCutHeatmapTOFe, simCutHeatmapTOFe, 5,
                        "TOFe", "TOFe", "chamber", "slat") << std::endl;
   }
   else
   {
      systematicsOutputFile << 0 << std::endl;
   }

   if (detectorsConfiguration[5] == '1') // TOFw
   {
      TH2F *realHeatmapTOFw = static_cast<TH2F *>(inputRealDataFile->Get("Heatmap: TOFw"));

      TH2F *simHeatmapTOFw = static_cast<TH2F *>(inputSimDataFile->Get("Heatmap: TOFw"));

      CheckHists(realHeatmapTOFw, simHeatmapTOFw, "Heatmap: TOFw");

      TH2F *realCutHeatmapTOFw = static_cast<TH2F *>(realHeatmapTOFw->Clone());
      TH2F *simCutHeatmapTOFw = static_cast<TH2F *>(simHeatmapTOFw->Clone());

      for (int i = 1; i <= realHeatmapTOFw->GetXaxis()->GetNbins(); i++)
      {
         for (int j = 1; j <= realHeatmapTOFw->GetYaxis()->GetNbins(); j++)
         {
            if (dmCutter.IsDeadTOFw(i - 1, j - 1))
            {
               realCutHeatmapTOFw->SetBinContent(i, j, 0.);
               simCutHeatmapTOFw->SetBinContent(i, j, 0.);
            }
         }
      }

      systematicsOutputFile << 
         CppTools::RMS(GetUncertainty(realHeatmapTOFw, simHeatmapTOFw, 
                                      realCutHeatmapTOFw, simCutHeatmapTOFw, 4,
                                      "TOFw", "TOFw", "chamber", "strip")) << std::endl;
   }
   else
   {
      systematicsOutputFile << 0 << std::endl;
   }

   if (detectorsConfiguration[6] == '1') // EMCal
   {
      for (int i = 0; i < 4; i++)
      {
         TH2F *realHeatmapEMCale = static_cast<TH2F *>
            (inputRealDataFile->Get(("Heatmap: EMCale" + std::to_string(i)).c_str()));

         TH2F *simHeatmapEMCale = static_cast<TH2F *>
            (inputSimDataFile->Get(("Heatmap: EMCale" + std::to_string(i)).c_str()));

         CheckHists(realHeatmapEMCale, simHeatmapEMCale, "Heatmap: EMCale" + std::to_string(i));

         TH2F *realCutHeatmapEMCale = static_cast<TH2F *>(realHeatmapEMCale->Clone());

         TH2F *simCutHeatmapEMCale = static_cast<TH2F *>(simHeatmapEMCale->Clone());

         for (int j = 1; j <= realHeatmapEMCale->GetXaxis()->GetNbins(); j++)
         {
            for (int k = 1; k <= realHeatmapEMCale->GetYaxis()->GetNbins(); k++)
            {
               if (dmCutter.IsDeadEMCal(0, i, j - 1, k - 1))
               {
                  realCutHeatmapEMCale->SetBinContent(j, k, 0.);
                  simCutHeatmapEMCale->SetBinContent(j, k, 0.);
               }
            }
         }

         systematicsOutputFile << 
            GetUncertainty(realHeatmapEMCale, simHeatmapEMCale, 
                           realCutHeatmapEMCale, simCutHeatmapEMCale, 8,
                           "EMCale" + std::to_string(i), "EMCale" + std::to_string(i), 
                           "y_{tower}", "z_{tower}");
      }
      systematicsOutputFile << std::endl;
      for (int i = 0; i < 4; i++)
      {
         TH2F *realHeatmapEMCalw = static_cast<TH2F *>
            (inputRealDataFile->Get(("Heatmap: EMCalw" + std::to_string(i)).c_str()));

         TH2F *simHeatmapEMCalw = static_cast<TH2F *>
            (inputSimDataFile->Get(("Heatmap: EMCalw" + std::to_string(i)).c_str()));

         CheckHists(realHeatmapEMCalw, simHeatmapEMCalw, "Heatmap: EMCalw" + std::to_string(i));

         TH2F *realCutHeatmapEMCalw = static_cast<TH2F *>(realHeatmapEMCalw->Clone());

         TH2F *simCutHeatmapEMCalw = static_cast<TH2F *>(simHeatmapEMCalw->Clone());

         for (int j = 1; j <= realHeatmapEMCalw->GetXaxis()->GetNbins(); j++)
         {
            for (int k = 1; k <= realHeatmapEMCalw->GetYaxis()->GetNbins(); k++)
            {
               if (dmCutter.IsDeadEMCal(1, i, j - 1, k - 1))
               {
                  realCutHeatmapEMCalw->SetBinContent(j, k, 0.);
                  simCutHeatmapEMCalw->SetBinContent(j, k, 0.);
               }
            }
         }

         systematicsOutputFile << 
            GetUncertainty(realHeatmapEMCalw, simHeatmapEMCalw, 
                           realCutHeatmapEMCalw, simCutHeatmapEMCalw, 8,
                           "EMCalw" + std::to_string(i), "EMCalw" + std::to_string(i), 
                           "y_{tower}", "z_{tower}");
         if (i < 3) systematicsOutputFile << " ";
      }
      systematicsOutputFile << std::endl;
   }
   else
   {
      std::ofstream systematicsOutputFile(outputDirParameters + "EMCal.txt");
      for (int i = 0; i < 4; i++)
      {
         systematicsOutputFile << 0 << " ";
      }
      systematicsOutputFile << std::endl;
      for (int i = 0; i < 4; i++)
      {
         systematicsOutputFile << 0;
         if (i < 3) systematicsOutputFile << " ";
      }
      systematicsOutputFile << std::endl;
   }

	table.End();
}

#endif /* DEAD_MAP_SYS_CPP */
