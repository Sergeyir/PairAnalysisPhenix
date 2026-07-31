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

void DeadMapSys::DrawDeadmap(TH2F *&heatmap, TH2F *&cutHeatmap,
                             const std::string& detectorName, const std::string& title,
                             const std::string& xTitle, const std::string& yTitle)
{	
   heatmap->SetMinimum(0.5);
   cutHeatmap->SetMinimum(0.5);

	TCanvas canv = TCanvas("canv", "canv", 900, 450);

	canv.Divide(2);

	canv.cd(1);
   gPad->SetRightMargin(0.14); gPad->SetTopMargin(0.07); 
   gPad->SetLeftMargin(0.125); gPad->SetBottomMargin(0.105);

   ROOTTools::DrawFrame(heatmap, title, xTitle, yTitle, 
                        1., 1.2, 0.05, 0.05, true, true, "COLZ");
	
	canv.cd(2);
   gPad->SetRightMargin(0.14); gPad->SetTopMargin(0.07); 
   gPad->SetLeftMargin(0.125); gPad->SetBottomMargin(0.105);

   ROOTTools::DrawFrame(cutHeatmap, "Cut " + title, xTitle, yTitle, 
                        1., 1.2, 0.05, 0.05, true, true, "COLZ");

   ROOTTools::PrintCanvas(&canv, "output/Deadmaps/" + runName + "/" + detectorName);
}

double DeadMapSys::GetUncertainty(TH2F *&realHeatmap, TH2F *&simHeatmap, 
                                  TH2F *&realCutHeatmap, TH2F *&simCutHeatmap, 
                                  const int numberOfHeatmapDivisions,
                                  const std::string& detectorName, const std::string& title, 
                                  const std::string& xTitle, const std::string& yTitle,
                                  const int rebinX, const int rebinY, const bool drawYProj)
{	
	TH1D *realCutHeatmapProjX = realCutHeatmap->
      ProjectionX((title + "_real_x").c_str(), 1, realCutHeatmap->GetYaxis()->GetNbins());
	TH1D *simCutHeatmapProjX = simCutHeatmap->
      ProjectionX((title + "_sim_x").c_str(), 1, simCutHeatmap->GetYaxis()->GetNbins());

   realCutHeatmap->SetMinimum(0.5);
   simCutHeatmap->SetMinimum(1e-15);

   realCutHeatmapProjX->RebinX(rebinX);
   simCutHeatmapProjX->RebinX(rebinX);

	realCutHeatmapProjX->SetFillColorAlpha(kOrange - 4, 0.5);
	
	realCutHeatmapProjX->Scale(1./realCutHeatmapProjX->Integral(1, 
                            realCutHeatmapProjX->GetXaxis()->GetNbins()), "nosw2");
	simCutHeatmapProjX->Scale(1./simCutHeatmapProjX->Integral(1, 
                           simCutHeatmapProjX->GetXaxis()->GetNbins()), "nosw2");

	TLegend projXLegend = TLegend(0.4, 0.78, 0.9, 0.88);
	
	projXLegend.SetNColumns(2);
	
	projXLegend.SetLineColorAlpha(0, 0);
	projXLegend.SetFillColorAlpha(0, 0);
	
	realCutHeatmapProjX->SetLineColor(kRed - 3);
	simCutHeatmapProjX->SetLineColor(kAzure - 3);

	realCutHeatmapProjX->SetLineWidth(2);

	projXLegend.AddEntry(realCutHeatmapProjX, "Data");
	projXLegend.AddEntry(simCutHeatmapProjX, "MC");
	
	realCutHeatmapProjX->SetTitle((title + " X projections data vs MC").c_str());
	realCutHeatmapProjX->GetXaxis()->SetTitle(xTitle.c_str());

	realCutHeatmapProjX->SetMaximum(CppTools::Maximum(realCutHeatmapProjX->GetMaximum()*1.3, 
                                                   simCutHeatmapProjX->GetMaximum()*1.3));

   if (drawYProj)
   {
      TH1D *realCutHeatmapProjY = realCutHeatmap->
         ProjectionY((title + "_real_y").c_str(), 1, realCutHeatmap->GetXaxis()->GetNbins());
      TH1D *simCutHeatmapProjY = simCutHeatmap->
         ProjectionY((title + "_sim_y").c_str(), 1, simCutHeatmap->GetXaxis()->GetNbins());

      realCutHeatmapProjY->RebinX(rebinY);
      simCutHeatmapProjY->RebinX(rebinY);

      realCutHeatmapProjY->SetFillColorAlpha(kOrange - 4, 0.5);
      
      realCutHeatmapProjY->Scale(1./realCutHeatmapProjY->Integral(1, 
                               realCutHeatmapProjY->GetXaxis()->GetNbins()), "nosw2");
      simCutHeatmapProjY->Scale(1./simCutHeatmapProjY->Integral(1, 
                              simCutHeatmapProjY->GetXaxis()->GetNbins()), "nosw2");

      TLegend projYLegend = TLegend(0.4, 0.78, 0.9, 0.88);

      projYLegend.SetNColumns(2);

      projYLegend.SetLineColorAlpha(0, 0);
      projYLegend.SetFillColorAlpha(0, 0);

      projYLegend.AddEntry(realCutHeatmapProjY, "data");
      projYLegend.AddEntry(simCutHeatmapProjY, "MC");

      realCutHeatmapProjY->SetTitle((title + " Y projections data vs MC").c_str());
      realCutHeatmapProjY->GetXaxis()->SetTitle(yTitle.c_str());

      realCutHeatmapProjY->SetMaximum(CppTools::Maximum(realCutHeatmapProjY->GetMaximum()*1.3, 
                                                      simCutHeatmapProjY->GetMaximum()*1.3));

      realCutHeatmapProjY->GetXaxis()->SetLabelSize(0.05);
      realCutHeatmapProjY->GetYaxis()->SetLabelSize(0.05);

      realCutHeatmapProjY->SetLineColor(kRed-3);
      simCutHeatmapProjY->SetLineColor(kAzure-3);

      TCanvas projCanv = TCanvas("projCanv", "canv", 1400, 1200);
      
      projCanv.cd();
      gPad->Divide(2, 2);
      
      projCanv.cd(1);
      gPad->SetRightMargin(0.13); gPad->SetTopMargin(0.07); 
      gPad->SetLeftMargin(0.1); gPad->SetBottomMargin(0.1);

      ROOTTools::DrawFrame(realCutHeatmap, "Cut " + title, xTitle, yTitle, 
                           0.95, 1.1, 0.045, 0.045, true, true, "COLZ");
      
      projCanv.cd(3);
      gPad->SetRightMargin(0.13); gPad->SetTopMargin(0.07); 
      gPad->SetLeftMargin(0.1); gPad->SetBottomMargin(0.1);

      ROOTTools::DrawFrame(simCutHeatmap, "MC cut " + title, xTitle, yTitle, 
                           0.95, 1.1, 0.045, 0.045, true, true, "COLZ");
  
      projCanv.cd(2);
      gPad->SetRightMargin(0.01); gPad->SetTopMargin(0.07); 
      gPad->SetLeftMargin(0.11); gPad->SetBottomMargin(0.1);

      realCutHeatmapProjX->SetMinimum(CppTools::Minimum(realCutHeatmapProjX->GetMinimum(), 
                                                        simCutHeatmapProjX->GetMinimum())/1.1);
      realCutHeatmapProjX->SetMaximum(CppTools::Maximum(realCutHeatmapProjX->GetMaximum(), 
                                                        simCutHeatmapProjX->GetMaximum())*1.1);

      ROOTTools::DrawFrame(realCutHeatmapProjX, "X projections of cut " + title, xTitle, "", 
                           0.9, 0., 0.045, 0.045, true, true);

      simCutHeatmapProjX->Draw("SAME HIST");
      projXLegend.Draw();

      projCanv.cd(4);
      gPad->SetRightMargin(0.01); gPad->SetTopMargin(0.07); 
      gPad->SetLeftMargin(0.11); gPad->SetBottomMargin(0.1);

      realCutHeatmapProjY->SetMinimum(CppTools::Minimum(realCutHeatmapProjY->GetMinimum(), 
                                                        simCutHeatmapProjY->GetMinimum()));
      realCutHeatmapProjY->SetMaximum(CppTools::Maximum(realCutHeatmapProjY->GetMaximum(), 
                                                        simCutHeatmapProjY->GetMaximum()));

      ROOTTools::DrawFrame(realCutHeatmapProjY, "Y projections of cut " + title, yTitle, "", 
                           0.9, 0., 0.045, 0.045, true, true);
      simCutHeatmapProjY->Draw("SAME HIST");
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

      ROOTTools::DrawFrame(realCutHeatmap, "Cut " + title, xTitle, yTitle, 
                           0.8, 1.0, 0.05, 0.05, true, true, "COLZ");
      
      projCanv.cd(2);
      gPad->SetRightMargin(0.145); gPad->SetTopMargin(0.07); 
      gPad->SetLeftMargin(0.09); gPad->SetBottomMargin(0.09);
      simCutHeatmap->Draw("colz");

      ROOTTools::DrawFrame(simCutHeatmap, "Cut MC " + title, xTitle, yTitle, 
                           0.8, 1.0, 0.05, 0.05, true, true, "COLZ");
      
      projCanv.cd(3);
      gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.07); 
      gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.09);

      realCutHeatmapProjX->SetMinimum(CppTools::Minimum(realCutHeatmapProjX->GetMinimum(), 
                                                        simCutHeatmapProjX->GetMinimum()));
      realCutHeatmapProjX->SetMaximum(CppTools::Maximum(realCutHeatmapProjX->GetMaximum(), 
                                                        simCutHeatmapProjX->GetMaximum()));

      ROOTTools::DrawFrame(realCutHeatmapProjX, "X projections of cut " + title, xTitle, "", 
                           0.8, 0., 0.05, 0.05, true, true);
      simCutHeatmapProjX->Draw("SAME HIST");
      projXLegend.Draw();

      ROOTTools::PrintCanvas(&projCanv, "output/Systematics/" + runName + "/" + detectorName);
   }

   const double realDataLost = (1. - realCutHeatmap->Integral()/realHeatmap->Integral())*100.;
   const double simDataLost = (1. - simCutHeatmap->Integral()/simHeatmap->Integral())*100.;

   const double realHeatmapScaling = 
      realCutHeatmap->Integral(1, realCutHeatmap->GetXaxis()->GetNbins(),
                               1, realCutHeatmap->GetYaxis()->GetNbins());
   const double simHeatmapScaling  = 
      simCutHeatmap->Integral(1, simCutHeatmap->GetXaxis()->GetNbins(),
                              1, simCutHeatmap->GetYaxis()->GetNbins());

   std::vector<double> divisionRealIntegral;
   std::vector<double> divisionSimIntegral;

   divisionRealIntegral.resize(numberOfHeatmapDivisions);
   divisionSimIntegral.resize(numberOfHeatmapDivisions);

   // bins for each division are chosen uniformly
   for (int i = 1; i <= realCutHeatmap->GetXaxis()->GetNbins(); i++)
   {
      for (int j = 1; j <= realCutHeatmap->GetYaxis()->GetNbins(); j++)
      {
         divisionRealIntegral[realCutHeatmap->GetBin(i, j) % numberOfHeatmapDivisions] += 
            realCutHeatmap->GetBinContent(i, j)/realHeatmapScaling;
         divisionSimIntegral[realCutHeatmap->GetBin(i, j) % numberOfHeatmapDivisions] += 
            simCutHeatmap->GetBinContent(i, j)/simHeatmapScaling;
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
   if (argc < 2) 
   {
      std::string errMsg = "Expected 1 parameter while " + std::to_string(argc - 1) + " ";
      errMsg += "parameter(s) were provided \n Usage: bin/DeadMapSys inputYAMLName";
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
   dmCutterMC.Initialize(runName, detectorsConfiguration, "data/Parameters/SimDeadmaps");

	table.Begin("Heatmap acceptance info");
	table.PrintHeader("detector", "uncertainty", "data lost", "MC lost");

   std::array<double, 2> reweightDCe{1., 1.};
   std::array<double, 2> reweightDCw{1., 1.};
   double reweightPC1e = 1.;
   double reweightPC1w = 1.;
   double reweightPC2 = 1.;
   double reweightPC3e = 1.;
   double reweightPC3w = 1.;
   double reweightTOFe = 1.;
   double reweightTimingTOFe = 1.;
   double reweightTOFw = 1.;
   double reweightTimingTOFw = 1.;
   std::array<double, 4> reweightEMCale{1., 1., 1., 1.};
   std::array<double, 4> reweightEMCalw{1., 1., 1., 1.};

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

      // for this project X and Y axis for different DC heatmaps are the same
      for (int i = 1; i <= realHeatmapDCe0->GetXaxis()->GetNbins(); i++)
      {
         const double x = realHeatmapDCe0->GetXaxis()->GetBinCenter(i);
         for (int j = 1; j <= realHeatmapDCe0->GetYaxis()->GetNbins(); j++)
         {
            const double y = realHeatmapDCe0->GetYaxis()->GetBinCenter(j);

            if (dmCutter.IsDeadDC(0, 1., x, y))
            {
               realCutHeatmapDCe0->SetBinContent(i, j, 0.);
               simCutHeatmapDCe0->SetBinContent(i, j, 0.);
            }
            if (dmCutter.IsDeadDC(0, -1., x, y))
            {
               realCutHeatmapDCe1->SetBinContent(i, j, 0.);
               simCutHeatmapDCe1->SetBinContent(i, j, 0.);
            }
            if (dmCutter.IsDeadDC(1, 1., x, y))
            {
               realCutHeatmapDCw0->SetBinContent(i, j, 0.);
               simCutHeatmapDCw0->SetBinContent(i, j, 0.);
            }
            if (dmCutter.IsDeadDC(1, -1., x, y))
            {
               realCutHeatmapDCw1->SetBinContent(i, j, 0.);
               simCutHeatmapDCw1->SetBinContent(i, j, 0.);
            }
         }
      }

      DrawDeadmap(realHeatmapDCe0, realCutHeatmapDCe0,
                  "DCe0", "DC east, #it{z}_{DC}#geq0", "board", "#it{#alpha}");
      DrawDeadmap(realHeatmapDCe1, realCutHeatmapDCe1,
                  "DCe1", "DC east, #it{z}_{DC}<0", "board", "#it{#alpha}");
      DrawDeadmap(realHeatmapDCw0, realCutHeatmapDCw0,
                  "DCw0", "DC west, #it{z}_{DC}#geq0", "board", "#it{#alpha}");
      DrawDeadmap(realHeatmapDCw1, realCutHeatmapDCw1,
                  "DCw1", "DC west, #it{z}_{DC}<0", "board", "#it{#alpha}");

      systematicsOutputFile << 
         GetUncertainty(realHeatmapDCe0, simHeatmapDCe0, realCutHeatmapDCe0, simCutHeatmapDCe0, 10,
                        "DCe0", "DC east, #it{z}_{DC}#geq0", 
                        "board", "#it{#alpha}", 3, 1, false) << " " <<
         GetUncertainty(realHeatmapDCe1, simHeatmapDCe1, realCutHeatmapDCe1, simCutHeatmapDCe1, 10,
                        "DCe1", "DC east, #it{z}_{DC}<0", 
                        "board", "#it{#alpha}", 3, 1, false) << " " <<
         GetUncertainty(realHeatmapDCw0, simHeatmapDCw0, realCutHeatmapDCw0, simCutHeatmapDCw0, 10,
                        "DCw0", "DC west, #it{z}_{DC}#geq0", 
                        "board", "#it{#alpha}", 3, 1, false) << " " <<
         GetUncertainty(realHeatmapDCw1, simHeatmapDCw1, realCutHeatmapDCw1, simCutHeatmapDCw1, 10,
                        "DCw1", "DC west, #it{z}_{DC}<0", "board", "#it{#alpha}", 3, 1, false) <<
         std::endl;

      /* unused and not needed for now
      // setting uncut heatmaps to be the heatmaps after fiducial cuts on real data
      // before starting additional fiducial cuts on MC
      simHeatmapDCe0 = static_cast<TH2F *>(simCutHeatmapDCe0->Clone());
      simHeatmapDCe1 = static_cast<TH2F *>(simCutHeatmapDCe1->Clone());
      simHeatmapDCw0 = static_cast<TH2F *>(simCutHeatmapDCw0->Clone());
      simHeatmapDCw1 = static_cast<TH2F *>(simCutHeatmapDCw1->Clone());

      realHeatmapDCe0 = static_cast<TH2F *>(realCutHeatmapDCe0->Clone());
      realHeatmapDCe1 = static_cast<TH2F *>(realCutHeatmapDCe1->Clone());
      realHeatmapDCw0 = static_cast<TH2F *>(realCutHeatmapDCw0->Clone());
      realHeatmapDCw1 = static_cast<TH2F *>(realCutHeatmapDCw1->Clone());

      // for this project X and Y axis for different DC heatmaps are the same
      for (int i = 1; i <= realHeatmapDCe0->GetXaxis()->GetNbins(); i++)
      {
         const double x = realHeatmapDCe0->GetXaxis()->GetBinCenter(i);
         for (int j = 1; j <= realHeatmapDCe0->GetYaxis()->GetNbins(); j++)
         {
            const double y = realHeatmapDCe0->GetYaxis()->GetBinCenter(j);

            if (dmCutterMC.IsDeadDC(0, 1., x, y))
            {
               realCutHeatmapDCe0->SetBinContent(i, j, 0.);
               simCutHeatmapDCe0->SetBinContent(i, j, 0.);
            }
            if (dmCutterMC.IsDeadDC(0, -1., x, y))
            {
               realCutHeatmapDCe1->SetBinContent(i, j, 0.);
               simCutHeatmapDCe1->SetBinContent(i, j, 0.);
            }
            if (dmCutterMC.IsDeadDC(1, 1., x, y))
            {
               realCutHeatmapDCw0->SetBinContent(i, j, 0.);
               simCutHeatmapDCw0->SetBinContent(i, j, 0.);
            }
            if (dmCutterMC.IsDeadDC(1, -1., x, y))
            {
               realCutHeatmapDCw1->SetBinContent(i, j, 0.);
               simCutHeatmapDCw1->SetBinContent(i, j, 0.);
            }
         }
      }

      DrawDeadmap(simHeatmapDCe0, simCutHeatmapDCe0,
                  "DCe0_MC", "DC east, #it{z}_{DC}#geq0", "board", "#it{#alpha}");
      DrawDeadmap(simHeatmapDCe1, simCutHeatmapDCe1,
                  "DCe1_MC", "DC east, #it{z}_{DC}<0", "board", "#it{#alpha}");
      DrawDeadmap(simHeatmapDCw0, simCutHeatmapDCw0,
                  "DCw0_MC", "DC west, #it{z}_{DC}#geq0", "board", "#it{#alpha}");
      DrawDeadmap(simHeatmapDCw1, simCutHeatmapDCw1,
                  "DCw1_MC", "DC west, #it{z}_{DC}<0", "board", "#it{#alpha}");

      systematicsOutputFile << 
         GetUncertainty(realHeatmapDCe0, simHeatmapDCe0, realCutHeatmapDCe0, simCutHeatmapDCe0, 10,
                        "DCe0_MC", "DC east, #it{z}_{DC}#geq0", 
                        "board", "#it{#alpha}", 3, 1, false) << " " <<
         GetUncertainty(realHeatmapDCe1, simHeatmapDCe1, realCutHeatmapDCe1, simCutHeatmapDCe1, 10,
                        "DCe1_MC", "DC east, #it{z}_{DC}<0", 
                        "board", "#it{#alpha}", 3, 1, false) << " " <<
         GetUncertainty(realHeatmapDCw0, simHeatmapDCw0, realCutHeatmapDCw0, simCutHeatmapDCw0, 10,
                        "DCw0_MC", "DC west, #it{z}_{DC}#geq0", 
                        "board", "#it{#alpha}", 3, 1, false) << " " <<
         GetUncertainty(realHeatmapDCw1, simHeatmapDCw1, realCutHeatmapDCw1, simCutHeatmapDCw1, 10,
                        "DCw1_MC", "DC west, #it{z}_{DC}<0", 
                        "board", "#it{#alpha}", 3, 1, false) << std::endl;

      reweightDCe[0] = realHeatmapDCe0->Integral(1, realHeatmapDCe0->GetXaxis()->GetNbins(),
                                                 1, realHeatmapDCe0->GetYaxis()->GetNbins())/
                       realCutHeatmapDCe0->Integral(1, realHeatmapDCe0->GetXaxis()->GetNbins(),
                                                    1, realHeatmapDCe0->GetYaxis()->GetNbins());
      reweightDCe[1] = realHeatmapDCe1->Integral(1, realHeatmapDCe1->GetXaxis()->GetNbins(),
                                                 1, realHeatmapDCe1->GetYaxis()->GetNbins())/
                       realCutHeatmapDCe1->Integral(1, realHeatmapDCe1->GetXaxis()->GetNbins(),
                                                    1, realHeatmapDCe1->GetYaxis()->GetNbins());
      reweightDCw[0] = realHeatmapDCw0->Integral(1, realHeatmapDCw0->GetXaxis()->GetNbins(),
                                                 1, realHeatmapDCw0->GetYaxis()->GetNbins())/
                       realCutHeatmapDCw0->Integral(1, realHeatmapDCw0->GetXaxis()->GetNbins(),
                                                    1, realHeatmapDCw0->GetYaxis()->GetNbins());
      reweightDCw[1] = realHeatmapDCw1->Integral(1, realHeatmapDCw1->GetXaxis()->GetNbins(),
                                                 1, realHeatmapDCw1->GetYaxis()->GetNbins())/
                       realCutHeatmapDCw1->Integral(1, realHeatmapDCw1->GetXaxis()->GetNbins(),
                                                    1, realHeatmapDCw1->GetYaxis()->GetNbins());
      */
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
         const double x = realHeatmapPC1e->GetXaxis()->GetBinCenter(i);
         for (int j = 1; j <= realHeatmapPC1e->GetYaxis()->GetNbins(); j++)
         {
            const double y = realHeatmapPC1e->GetYaxis()->GetBinCenter(j);

            if (dmCutter.IsDeadPC1(0, x, y))
            {
               realCutHeatmapPC1e->SetBinContent(i, j, 0.);
               simCutHeatmapPC1e->SetBinContent(i, j, 0.);
            }
         }
      }
      for (int i = 1; i <= realHeatmapPC1w->GetXaxis()->GetNbins(); i++)
      {
         const double x = realHeatmapPC1w->GetXaxis()->GetBinCenter(i);
         for (int j = 1; j <= realHeatmapPC1w->GetYaxis()->GetNbins(); j++)
         {
            const double y = realHeatmapPC1w->GetYaxis()->GetBinCenter(j);

            if (dmCutter.IsDeadPC1(1, x, y))
            {
               realCutHeatmapPC1w->SetBinContent(i, j, 0.);
               simCutHeatmapPC1w->SetBinContent(i, j, 0.);
            }
         }
      }

      DrawDeadmap(realHeatmapPC1e, realCutHeatmapPC1e,
                  "PC1e", "PC1 east", "#it{z}_{PC1}", "#it{#varphi}_{PC1}");
      DrawDeadmap(realHeatmapPC1w, realCutHeatmapPC1w,
                  "PC1w", "PC1 west", "#it{z}_{PC1}", "#it{#varphi}_{PC1}");

      GetUncertainty(realHeatmapPC1e, simHeatmapPC1e, realCutHeatmapPC1e, simCutHeatmapPC1e, 10,
                     "PC1e", "PC1 east", "#it{z}_{PC1}", "#it{#varphi}_{PC1}", 2);
      GetUncertainty(realHeatmapPC1w, simHeatmapPC1w, realCutHeatmapPC1w, simCutHeatmapPC1w, 10,
                     "PC1w", "PC1 west", "#it{z}_{PC1}", "#it{#varphi}_{PC1}", 2);

      // setting uncut heatmaps to be the heatmaps after fiducial cuts on real data
      // before starting additional fiducial cuts on MC

      realHeatmapPC1e = static_cast<TH2F *>(realCutHeatmapPC1e->Clone());
      realHeatmapPC1w = static_cast<TH2F *>(realCutHeatmapPC1w->Clone());
      simHeatmapPC1e = static_cast<TH2F *>(simCutHeatmapPC1e->Clone());
      simHeatmapPC1w = static_cast<TH2F *>(simCutHeatmapPC1w->Clone());

      for (int i = 1; i <= realHeatmapPC1e->GetXaxis()->GetNbins(); i++)
      {
         const double x = realHeatmapPC1e->GetXaxis()->GetBinCenter(i);
         for (int j = 1; j <= realHeatmapPC1e->GetYaxis()->GetNbins(); j++)
         {
            const double y = realHeatmapPC1e->GetYaxis()->GetBinCenter(j);

            if (dmCutterMC.IsDeadPC1(0, x, y))
            {
               realCutHeatmapPC1e->SetBinContent(i, j, 0.);
               simCutHeatmapPC1e->SetBinContent(i, j, 0.);
            }
         }
      }
      for (int i = 1; i <= realHeatmapPC1w->GetXaxis()->GetNbins(); i++)
      {
         const double x = realHeatmapPC1w->GetXaxis()->GetBinCenter(i);
         for (int j = 1; j <= realHeatmapPC1w->GetYaxis()->GetNbins(); j++)
         {
            const double y = realHeatmapPC1w->GetYaxis()->GetBinCenter(j);

            if (dmCutterMC.IsDeadPC1(1, x, y))
            {
               realCutHeatmapPC1w->SetBinContent(i, j, 0.);
               simCutHeatmapPC1w->SetBinContent(i, j, 0.);
            }
         }
      }

      DrawDeadmap(simHeatmapPC1e, simCutHeatmapPC1e,
                  "PC1e_MC", "PC1 east", "#it{z}_{PC1}", "#it{#varphi}_{PC1}");
      DrawDeadmap(simHeatmapPC1w, simCutHeatmapPC1w,
                  "PC1w_MC", "PC1 west", "#it{z}_{PC1}", "#it{#varphi}_{PC1}");

      systematicsOutputFile << 
         GetUncertainty(realHeatmapPC1e, simHeatmapPC1e, realCutHeatmapPC1e, simCutHeatmapPC1e, 10,
                        "PC1e_MC", "PC1 east", "#it{z}_{PC1}", "#it{#varphi}_{PC1}", 2) << " " <<
         GetUncertainty(realHeatmapPC1w, simHeatmapPC1w, realCutHeatmapPC1w, simCutHeatmapPC1w, 10,
                        "PC1w_MC", "PC1 west", "#it{z}_{PC1}", "#it{#varphi}_{PC1}", 2) << 
         std::endl;

      reweightPC1e = realHeatmapPC1e->Integral(1, realHeatmapPC1e->GetXaxis()->GetNbins(),
                                             1, realHeatmapPC1e->GetYaxis()->GetNbins())/
                     realCutHeatmapPC1e->Integral(1, realHeatmapPC1e->GetXaxis()->GetNbins(),
                                                1, realHeatmapPC1e->GetYaxis()->GetNbins());
      reweightPC1w = realHeatmapPC1w->Integral(1, realHeatmapPC1w->GetXaxis()->GetNbins(),
                                             1, realHeatmapPC1w->GetYaxis()->GetNbins())/
                     realCutHeatmapPC1w->Integral(1, realHeatmapPC1w->GetXaxis()->GetNbins(),
                                                1, realHeatmapPC1w->GetYaxis()->GetNbins());
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
         const double x = realHeatmapPC2->GetXaxis()->GetBinCenter(i);
         for (int j = 1; j <= realHeatmapPC2->GetYaxis()->GetNbins(); j++)
         {
            const double y = realHeatmapPC2->GetYaxis()->GetBinCenter(j);

            if (dmCutter.IsDeadPC2(x, y))
            {
               realCutHeatmapPC2->SetBinContent(i, j, 0.);
               simCutHeatmapPC2->SetBinContent(i, j, 0.);
            }
         }
      }

      DrawDeadmap(realHeatmapPC2, realCutHeatmapPC2,
                  "PC2", "PC2", "#it{z}_{PC2}", "#it{#varphi}_{PC2}");

      GetUncertainty(realHeatmapPC2, simHeatmapPC2, realCutHeatmapPC2, simCutHeatmapPC2, 10,
                     "PC2", "PC2", "#it{z}_{PC2}", "#it{#varphi}_{PC2}");

      // setting uncut heatmaps to be the heatmaps after fiducial cuts on real data
      // before starting additional fiducial cuts on MC
      realHeatmapPC2 = static_cast<TH2F *>(realCutHeatmapPC2->Clone());
      simHeatmapPC2 = static_cast<TH2F *>(simCutHeatmapPC2->Clone());

      for (int i = 1; i <= realHeatmapPC2->GetXaxis()->GetNbins(); i++)
      {
         const double x = realHeatmapPC2->GetXaxis()->GetBinCenter(i);
         for (int j = 1; j <= realHeatmapPC2->GetYaxis()->GetNbins(); j++)
         {
            const double y = realHeatmapPC2->GetYaxis()->GetBinCenter(j);

            if (dmCutterMC.IsDeadPC2(x, y))
            {
               realCutHeatmapPC2->SetBinContent(i, j, 0.);
               simCutHeatmapPC2->SetBinContent(i, j, 0.);
            }
         }
      }

      DrawDeadmap(simHeatmapPC2, simCutHeatmapPC2,
                  "PC2_MC", "PC2", "#it{z}_{PC2}", "#it{#varphi}_{PC2}");

      systematicsOutputFile << 
         GetUncertainty(realHeatmapPC2, simHeatmapPC2, 
                        realCutHeatmapPC2, simCutHeatmapPC2, 10,
                        "PC2_MC", "PC2", "#it{z}_{PC2}", "#it{#varphi}_{PC2}") << std::endl;

      reweightPC2 = realHeatmapPC2->Integral(1, realHeatmapPC2->GetXaxis()->GetNbins(),
                                           1, realHeatmapPC2->GetYaxis()->GetNbins())/
                     realCutHeatmapPC2->Integral(1, realHeatmapPC2->GetXaxis()->GetNbins(),
                                               1, realHeatmapPC2->GetYaxis()->GetNbins());
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
         const double x = realHeatmapPC3e->GetXaxis()->GetBinCenter(i);
         for (int j = 1; j <= realHeatmapPC3e->GetYaxis()->GetNbins(); j++)
         {
            const double y = realHeatmapPC3e->GetYaxis()->GetBinCenter(j);

            if (dmCutter.IsDeadPC3(0, x, y))
            {
               realCutHeatmapPC3e->SetBinContent(i, j, 0.);
               simCutHeatmapPC3e->SetBinContent(i, j, 0.);
            }
         }
      }
      for (int i = 1; i <= realHeatmapPC3w->GetXaxis()->GetNbins(); i++)
      {
         const double x = realHeatmapPC3w->GetXaxis()->GetBinCenter(i);
         for (int j = 1; j <= realHeatmapPC3w->GetYaxis()->GetNbins(); j++)
         {
            const double y = realHeatmapPC3w->GetYaxis()->GetBinCenter(j);

            if (dmCutter.IsDeadPC3(1, x, y))
            {
               realCutHeatmapPC3w->SetBinContent(i, j, 0.);
               simCutHeatmapPC3w->SetBinContent(i, j, 0.);
            }
         }
      }

      DrawDeadmap(realHeatmapPC3e, realCutHeatmapPC3e,
                  "PC3e", "PC3 east", "#it{z}_{PC3}", "#it{#varphi}_{PC3}");
      DrawDeadmap(realHeatmapPC3w, realCutHeatmapPC3w,
                  "PC3w", "PC3 west", "#it{z}_{PC3}", "#it{#varphi}_{PC3}");

      systematicsOutputFile << 
         GetUncertainty(realHeatmapPC3e, simHeatmapPC3e, 
                        realCutHeatmapPC3e, simCutHeatmapPC3e, 10,
                        "PC3e", "PC3 east", "#it{z}_{PC3}", "#it{#varphi}_{PC3}", 2) << " " <<
         GetUncertainty(realHeatmapPC3w, simHeatmapPC3w, 
                        realCutHeatmapPC3w, simCutHeatmapPC3w, 10,
                        "PC3w", "PC3 west", "#it{z}_{PC3}", "#it{#varphi}_{PC3}", 2) << std::endl;

      // setting uncut heatmaps to be the heatmaps after fiducial cuts on real data
      // before starting additional fiducial cuts on MC

      realHeatmapPC3e = static_cast<TH2F *>(realCutHeatmapPC3e->Clone());
      realHeatmapPC3w = static_cast<TH2F *>(realCutHeatmapPC3w->Clone());
      simHeatmapPC3e = static_cast<TH2F *>(simCutHeatmapPC3e->Clone());
      simHeatmapPC3w = static_cast<TH2F *>(simCutHeatmapPC3w->Clone());

      for (int i = 1; i <= realHeatmapPC3e->GetXaxis()->GetNbins(); i++)
      {
         const double x = realHeatmapPC3e->GetXaxis()->GetBinCenter(i);
         for (int j = 1; j <= realHeatmapPC3e->GetYaxis()->GetNbins(); j++)
         {
            const double y = realHeatmapPC3e->GetYaxis()->GetBinCenter(j);

            if (dmCutterMC.IsDeadPC3(0, x, y))
            {
               realCutHeatmapPC3e->SetBinContent(i, j, 0.);
               simCutHeatmapPC3e->SetBinContent(i, j, 0.);
            }
         }
      }
      for (int i = 1; i <= realHeatmapPC3w->GetXaxis()->GetNbins(); i++)
      {
         const double x = realHeatmapPC3w->GetXaxis()->GetBinCenter(i);
         for (int j = 1; j <= realHeatmapPC3w->GetYaxis()->GetNbins(); j++)
         {
            const double y = realHeatmapPC3w->GetYaxis()->GetBinCenter(j);

            if (dmCutterMC.IsDeadPC3(1, x, y))
            {
               realCutHeatmapPC3w->SetBinContent(i, j, 0.);
               simCutHeatmapPC3w->SetBinContent(i, j, 0.);
            }
         }
      }

      DrawDeadmap(simHeatmapPC3e, simCutHeatmapPC3e,
                  "PC3e_MC", "PC3 east", "#it{z}_{PC3}", "#it{#varphi}_{PC3}");
      DrawDeadmap(simHeatmapPC3w, simCutHeatmapPC3w,
                  "PC3w_MC", "PC3 west", "#it{z}_{PC3}", "#it{#varphi}_{PC3}");

      systematicsOutputFile << 
         GetUncertainty(realHeatmapPC3e, simHeatmapPC3e, 
                        realCutHeatmapPC3e, simCutHeatmapPC3e, 10,
                        "PC3e_MC", "PC3 east", "#it{z}_{PC3}", "#it{#varphi}_{PC3}", 2) << " " <<
         GetUncertainty(realHeatmapPC3w, simHeatmapPC3w, 
                        realCutHeatmapPC3w, simCutHeatmapPC3w, 10,
                        "PC3w_MC", "PC3 west", "#it{z}_{PC3}", "#it{#varphi}_{PC3}", 2) << 
         std::endl;

      reweightPC3e = realHeatmapPC3e->Integral(1, realHeatmapPC3e->GetXaxis()->GetNbins(),
                                             1, realHeatmapPC3e->GetYaxis()->GetNbins())/
                     realCutHeatmapPC3e->Integral(1, realHeatmapPC3e->GetXaxis()->GetNbins(),
                                                1, realHeatmapPC3e->GetYaxis()->GetNbins());
      reweightPC3w = realHeatmapPC3w->Integral(1, realHeatmapPC3w->GetXaxis()->GetNbins(),
                                             1, realHeatmapPC3w->GetYaxis()->GetNbins())/
                     realCutHeatmapPC3w->Integral(1, realHeatmapPC3w->GetXaxis()->GetNbins(),
                                                1, realHeatmapPC3w->GetYaxis()->GetNbins());
   }
   else
   {
      systematicsOutputFile << 0 << " " << 0 << std::endl;
   }

   if (detectorsConfiguration[4] == '1') // TOFe
   {
      TH2F *realHeatmapTOFe = static_cast<TH2F *>(inputRealDataFile->Get("Heatmap: TOFe"));
      TH2F *simCutHeatmapTOFe = static_cast<TH2F *>(inputSimDataFile->Get("Heatmap: TOFe"));

      CheckHists(realHeatmapTOFe, simCutHeatmapTOFe, "Heatmap: TOFe");

      // hetmaps without eloss weighting for systematics
      TH2F *realHeatmapTOFeSys = static_cast<TH2F *>(inputRealDataFile->Get("_Heatmap: TOFe hit"));
      TH2F *simHeatmapTOFeSys= static_cast<TH2F *>(inputSimDataFile->Get("_Heatmap: TOFe hit"));

      CheckHists(realHeatmapTOFeSys, simHeatmapTOFeSys, "_Heatmap: TOFe hit");

      TH2F *realCutHeatmapTOFe = static_cast<TH2F *>(realHeatmapTOFe->Clone());
      TH2F *simMCCutHeatmapTOFe = static_cast<TH2F *>(simCutHeatmapTOFe->Clone());

      TH2F *realCutHeatmapTOFeSys = static_cast<TH2F *>(realHeatmapTOFeSys->Clone());
      TH2F *simCutHeatmapTOFeSys = static_cast<TH2F *>(simHeatmapTOFeSys->Clone());

      TH2F *realMCCutHeatmapTOFeSys = static_cast<TH2F *>(realCutHeatmapTOFeSys->Clone());
      TH2F *simMCCutHeatmapTOFeSys = static_cast<TH2F *>(simCutHeatmapTOFeSys->Clone());

      for (int i = 1; i <= realHeatmapTOFe->GetXaxis()->GetNbins(); i++)
      {
         for (int j = 1; j <= realHeatmapTOFe->GetYaxis()->GetNbins(); j++)
         {
            if (dmCutter.IsDeadTOFe(i - 1, j - 1))
            {
               realCutHeatmapTOFe->SetBinContent(i, j, 0.);
               simCutHeatmapTOFe->SetBinContent(i, j, 0.);

               simMCCutHeatmapTOFe->SetBinContent(i, j, 0.);

               realCutHeatmapTOFeSys->SetBinContent(i, j, 0.);
               simCutHeatmapTOFeSys->SetBinContent(i, j, 0.);

               realMCCutHeatmapTOFeSys->SetBinContent(i, j, 0.);
               simMCCutHeatmapTOFeSys->SetBinContent(i, j, 0.);
            }
            if (dmCutterMC.IsDeadTOFe(i - 1, j - 1))
            {
               simMCCutHeatmapTOFe->SetBinContent(i, j, 0.);

               realMCCutHeatmapTOFeSys->SetBinContent(i, j, 0.);
               simMCCutHeatmapTOFeSys->SetBinContent(i, j, 0.);
            }
         }
      }

      DrawDeadmap(realHeatmapTOFe, realCutHeatmapTOFe,
                  "TOFe", "TOFe", "#it{Z}_{slat}", "#it{Y}_{slat}");
      DrawDeadmap(simCutHeatmapTOFe, simMCCutHeatmapTOFe,
                  "TOFe_MC", "TOFe", "#it{Z}_{slat}", "#it{Y}_{slat}");

      GetUncertainty(realHeatmapTOFeSys, simHeatmapTOFeSys, 
                     realCutHeatmapTOFeSys, simCutHeatmapTOFeSys, 5,
                     "TOFe", "TOFe", "#it{Z}_{slat}", "#it{Y}_{slat}");

      systematicsOutputFile << 
         GetUncertainty(realHeatmapTOFeSys, simHeatmapTOFeSys, 
                        realMCCutHeatmapTOFeSys, simMCCutHeatmapTOFeSys, 5,
                        "TOFe_MC", "TOFe", "#it{Z}_{slat}", "#it{Y}_{slat}") << " ";

      reweightTOFe = 
         realCutHeatmapTOFeSys->Integral(1, realHeatmapTOFe->GetXaxis()->GetNbins(),
                                         1, realHeatmapTOFe->GetYaxis()->GetNbins())/
         realMCCutHeatmapTOFeSys->Integral(1, realHeatmapTOFe->GetXaxis()->GetNbins(),
                                           1, realHeatmapTOFe->GetYaxis()->GetNbins());

      for (int i = 1; i <= realHeatmapTOFe->GetXaxis()->GetNbins(); i++)
      {
         for (int j = 1; j <= realHeatmapTOFe->GetYaxis()->GetNbins(); j++)
         {
            if (dmCutter.IsDeadTimingTOFe(i - 1, j - 1))
            {
               realCutHeatmapTOFe->SetBinContent(i, j, 0.);
               simCutHeatmapTOFe->SetBinContent(i, j, 0.);

               realCutHeatmapTOFeSys->SetBinContent(i, j, 0.);
               simCutHeatmapTOFeSys->SetBinContent(i, j, 0.);

               realMCCutHeatmapTOFeSys->SetBinContent(i, j, 0.);
               simMCCutHeatmapTOFeSys->SetBinContent(i, j, 0.);
            }
            if (dmCutterMC.IsDeadTimingTOFe(i - 1, j - 1))
            {
               realMCCutHeatmapTOFeSys->SetBinContent(i, j, 0.);
               simMCCutHeatmapTOFeSys->SetBinContent(i, j, 0.);
            }
         }
      }

      DrawDeadmap(realCutHeatmapTOFe, realCutHeatmapTOFe,
                  "TimingTOFe", "TOFe", "#it{Z}_{slat}", "#it{Y}_{slat}");
      DrawDeadmap(simCutHeatmapTOFe, simMCCutHeatmapTOFe,
                  "TimingTOFe_MC", "TOFe", "#it{Z}_{slat}", "#it{Y}_{slat}");

      GetUncertainty(realHeatmapTOFeSys, simHeatmapTOFeSys, 
                     realCutHeatmapTOFeSys, simCutHeatmapTOFeSys, 5,
                     "TimingTOFe", "TOFe", "#it{Z}_{slat}", "#it{Y}_{slat}");

      systematicsOutputFile << 
         GetUncertainty(realHeatmapTOFeSys, simHeatmapTOFeSys, 
                        realMCCutHeatmapTOFeSys, simMCCutHeatmapTOFeSys, 5,
                        "TimingTOFe_MC", "TOFe", "#it{Z}_{slat}", "#it{Y}_{slat}") << std::endl;

      reweightTimingTOFe = 
         realCutHeatmapTOFeSys->Integral(1, realHeatmapTOFe->GetXaxis()->GetNbins(),
                                         1, realHeatmapTOFe->GetYaxis()->GetNbins())/
         realMCCutHeatmapTOFeSys->Integral(1, realHeatmapTOFe->GetXaxis()->GetNbins(),
                                           1, realHeatmapTOFe->GetYaxis()->GetNbins());
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

      TH2F *realMCCutHeatmapTOFw = static_cast<TH2F *>(realHeatmapTOFw->Clone());
      TH2F *simMCCutHeatmapTOFw = static_cast<TH2F *>(simHeatmapTOFw->Clone());

      for (int i = 1; i <= realHeatmapTOFw->GetXaxis()->GetNbins(); i++)
      {
         for (int j = 1; j <= realHeatmapTOFw->GetYaxis()->GetNbins(); j++)
         {
            if (dmCutter.IsDeadTOFw(i - 1, j - 1))
            {
               realCutHeatmapTOFw->SetBinContent(i, j, 0.);
               simCutHeatmapTOFw->SetBinContent(i, j, 0.);

               realMCCutHeatmapTOFw->SetBinContent(i, j, 0.);
               simMCCutHeatmapTOFw->SetBinContent(i, j, 0.);
            }
            if (dmCutterMC.IsDeadTOFw(i - 1, j - 1))
            {
               realMCCutHeatmapTOFw->SetBinContent(i, j, 0.);
               simMCCutHeatmapTOFw->SetBinContent(i, j, 0.);
            }
         }
      }

      DrawDeadmap(realHeatmapTOFw, realCutHeatmapTOFw, "TOFw", "TOFw", "#it{Z}_{strip}", "#it{Y}_{strip}");
      DrawDeadmap(simCutHeatmapTOFw, simMCCutHeatmapTOFw, "TOFw_MC", "TOFw", "#it{Z}_{strip}", "#it{Y}_{strip}");

      GetUncertainty(realHeatmapTOFw, simHeatmapTOFw, realCutHeatmapTOFw, simCutHeatmapTOFw, 4,
                     "TOFw", "TOFw", "#it{Z}_{strip}", "#it{Y}_{strip}");

      systematicsOutputFile << 
         GetUncertainty(realHeatmapTOFw, simHeatmapTOFw, 
                        realMCCutHeatmapTOFw, simMCCutHeatmapTOFw, 4,
                        "TOFw_MC", "TOFw", "#it{Z}_{strip}", "#it{Y}_{strip}") << " ";

      reweightTOFw = 
         realCutHeatmapTOFw->Integral(1, realHeatmapTOFw->GetXaxis()->GetNbins(),
                                      1, realHeatmapTOFw->GetYaxis()->GetNbins())/
         realMCCutHeatmapTOFw->Integral(1, realHeatmapTOFw->GetXaxis()->GetNbins(),
                                        1, realHeatmapTOFw->GetYaxis()->GetNbins());

      for (int i = 1; i <= realHeatmapTOFw->GetXaxis()->GetNbins(); i++)
      {
         for (int j = 1; j <= realHeatmapTOFw->GetYaxis()->GetNbins(); j++)
         {
            if (dmCutter.IsDeadTimingTOFw(i - 1, j - 1))
            {
               realCutHeatmapTOFw->SetBinContent(i, j, 0.);
               simCutHeatmapTOFw->SetBinContent(i, j, 0.);

               realMCCutHeatmapTOFw->SetBinContent(i, j, 0.);
               simMCCutHeatmapTOFw->SetBinContent(i, j, 0.);
            }
            if (dmCutterMC.IsDeadTimingTOFw(i - 1, j - 1))
            {
               realMCCutHeatmapTOFw->SetBinContent(i, j, 0.);
               simMCCutHeatmapTOFw->SetBinContent(i, j, 0.);
            }
         }
      }

      DrawDeadmap(realHeatmapTOFw, realCutHeatmapTOFw, "TimingTOFw", "TOFw", "#it{Z}_{strip}", "#it{Y}_{strip}");
      DrawDeadmap(simCutHeatmapTOFw, simMCCutHeatmapTOFw, 
                  "TimingTOFw_MC", "TOFw", "#it{Z}_{strip}", "#it{Y}_{strip}");

     GetUncertainty(realHeatmapTOFw, simHeatmapTOFw, realCutHeatmapTOFw, simCutHeatmapTOFw, 4,
                    "TimingTOFw", "TOFw", "#it{Z}_{strip}", "#it{Y}_{strip}");
      systematicsOutputFile << 
        GetUncertainty(realHeatmapTOFw, simHeatmapTOFw, 
                       realMCCutHeatmapTOFw, simMCCutHeatmapTOFw, 4,
                       "TimingTOFw_MC", "TOFw", "#it{Z}_{strip}", "#it{Y}_{strip}") << std::endl;

      reweightTimingTOFw = 
         realCutHeatmapTOFw->Integral(1, realHeatmapTOFw->GetXaxis()->GetNbins(),
                                      1, realHeatmapTOFw->GetYaxis()->GetNbins())/
         realMCCutHeatmapTOFw->Integral(1, realHeatmapTOFw->GetXaxis()->GetNbins(),
                                        1, realHeatmapTOFw->GetYaxis()->GetNbins());
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
         TH2F *simCutHeatmapEMCale = static_cast<TH2F *>
            (inputSimDataFile->Get(("Heatmap: EMCale" + std::to_string(i)).c_str()));

         CheckHists(realHeatmapEMCale, simCutHeatmapEMCale, "Heatmap: EMCale" + std::to_string(i));

         // hetmaps without ecore weighting for systematics
         TH2F *realHeatmapEMCaleSys = static_cast<TH2F *>
            (inputRealDataFile->Get(("_Heatmap: EMCale hit" + std::to_string(i)).c_str()));
         TH2F *simHeatmapEMCaleSys = static_cast<TH2F *>
            (inputSimDataFile->Get(("_Heatmap: EMCale" + std::to_string(i) + " hit").c_str()));

         CheckHists(realHeatmapEMCaleSys, simHeatmapEMCaleSys, 
                    "_Heatmap: EMCale" + std::to_string(i) + " hit");

         TH2F *realCutHeatmapEMCale = static_cast<TH2F *>(realHeatmapEMCale->Clone());
         TH2F *simMCCutHeatmapEMCale = static_cast<TH2F *>(simCutHeatmapEMCale->Clone());

         TH2F *realCutHeatmapEMCaleSys = static_cast<TH2F *>(realHeatmapEMCaleSys->Clone());
         TH2F *simCutHeatmapEMCaleSys = static_cast<TH2F *>(simHeatmapEMCaleSys->Clone());

         TH2F *realMCCutHeatmapEMCaleSys = static_cast<TH2F *>(realHeatmapEMCaleSys->Clone());
         TH2F *simMCCutHeatmapEMCaleSys = static_cast<TH2F *>(simHeatmapEMCaleSys->Clone());

         for (int j = 1; j <= realHeatmapEMCale->GetXaxis()->GetNbins(); j++)
         {
            for (int k = 1; k <= realHeatmapEMCale->GetYaxis()->GetNbins(); k++)
            {
               if (dmCutter.IsDeadEMCal(0, i, j - 1, k - 1))
               {
                  realCutHeatmapEMCale->SetBinContent(j, k, 0.);
                  simCutHeatmapEMCale->SetBinContent(j, k, 0.);
                  simMCCutHeatmapEMCale->SetBinContent(j, k, 0.);

                  realCutHeatmapEMCaleSys->SetBinContent(j, k, 0.);
                  simCutHeatmapEMCaleSys->SetBinContent(j, k, 0.);

                  realMCCutHeatmapEMCaleSys->SetBinContent(j, k, 0.);
                  simMCCutHeatmapEMCaleSys->SetBinContent(j, k, 0.);
               }
               if (dmCutterMC.IsDeadEMCal(0, i, j - 1, k - 1))
               {
                  simMCCutHeatmapEMCale->SetBinContent(j, k, 0.);

                  realMCCutHeatmapEMCaleSys->SetBinContent(j, k, 0.);
                  simMCCutHeatmapEMCaleSys->SetBinContent(j, k, 0.);
               }
            }
         }

         DrawDeadmap(realHeatmapEMCale, realCutHeatmapEMCale,
                     "EMCale" + std::to_string(i), "EMCale" + std::to_string(i), 
                     "#it{Y}_{tower}", "#it{Z}_{tower}");

         DrawDeadmap(simCutHeatmapEMCale, simMCCutHeatmapEMCale,
                     "EMCale" + std::to_string(i) + "_MC", "EMCale" + std::to_string(i), 
                     "#it{Y}_{tower}", "#it{Z}_{tower}");

         GetUncertainty(realHeatmapEMCaleSys, simHeatmapEMCaleSys, 
                        realCutHeatmapEMCaleSys, simCutHeatmapEMCaleSys, 8,
                        "EMCale" + std::to_string(i), "EMCale" + std::to_string(i), 
                        "#it{Y}_{tower}", "#it{Z}_{tower}");

         systematicsOutputFile << 
            GetUncertainty(realHeatmapEMCaleSys, simHeatmapEMCaleSys, 
                           realMCCutHeatmapEMCaleSys, simMCCutHeatmapEMCaleSys, 8,
                           "EMCale" + std::to_string(i) + "_MC", "EMCale" + std::to_string(i), 
                           "#it{Y}_{tower}", "#it{Z}_{tower}");

         if (i < 3) systematicsOutputFile << " ";

         reweightEMCale[i] = 
            realCutHeatmapEMCaleSys->Integral(1, realHeatmapEMCale->GetXaxis()->GetNbins(),
                                              1, realHeatmapEMCale->GetYaxis()->GetNbins())/
            realMCCutHeatmapEMCaleSys->Integral(1, realHeatmapEMCale->GetXaxis()->GetNbins(),
                                                1, realHeatmapEMCale->GetYaxis()->GetNbins());
      }
      systematicsOutputFile << std::endl;

      for (int i = 0; i < 4; i++)
      {
         TH2F *realHeatmapEMCalw = static_cast<TH2F *>
            (inputRealDataFile->Get(("Heatmap: EMCalw" + std::to_string(i)).c_str()));
         TH2F *simCutHeatmapEMCalw = static_cast<TH2F *>
            (inputSimDataFile->Get(("Heatmap: EMCalw" + std::to_string(i)).c_str()));

         CheckHists(realHeatmapEMCalw, simCutHeatmapEMCalw, 
                    "Heatmap: EMCalw" + std::to_string(i));

         // hetmaps without ecore weighting for systematics
         TH2F *realHeatmapEMCalwSys = static_cast<TH2F *>
            (inputRealDataFile->Get(("_Heatmap: EMCalw hit" + std::to_string(i)).c_str()));
         TH2F *simHeatmapEMCalwSys = static_cast<TH2F *>
            (inputSimDataFile->Get(("_Heatmap: EMCalw" + std::to_string(i) + " hit").c_str()));

         CheckHists(realHeatmapEMCalwSys, simHeatmapEMCalwSys, 
                    "_Heatmap: EMCalw" + std::to_string(i) + " hit");

         TH2F *realCutHeatmapEMCalw = static_cast<TH2F *>(realHeatmapEMCalw->Clone());
         TH2F *simMCCutHeatmapEMCalw = static_cast<TH2F *>(simCutHeatmapEMCalw->Clone());

         TH2F *realCutHeatmapEMCalwSys = static_cast<TH2F *>(realHeatmapEMCalwSys->Clone());
         TH2F *simCutHeatmapEMCalwSys = static_cast<TH2F *>(simHeatmapEMCalwSys->Clone());

         TH2F *realMCCutHeatmapEMCalwSys = static_cast<TH2F *>(realHeatmapEMCalwSys->Clone());
         TH2F *simMCCutHeatmapEMCalwSys = static_cast<TH2F *>(simHeatmapEMCalwSys->Clone());

         for (int j = 1; j <= realHeatmapEMCalw->GetXaxis()->GetNbins(); j++)
         {
            for (int k = 1; k <= realHeatmapEMCalw->GetYaxis()->GetNbins(); k++)
            {
               if (dmCutter.IsDeadEMCal(1, i, j - 1, k - 1))
               {
                  realCutHeatmapEMCalw->SetBinContent(j, k, 0.);
                  simCutHeatmapEMCalw->SetBinContent(j, k, 0.);
                  simMCCutHeatmapEMCalw->SetBinContent(j, k, 0.);

                  realCutHeatmapEMCalwSys->SetBinContent(j, k, 0.);
                  simCutHeatmapEMCalwSys->SetBinContent(j, k, 0.);

                  realMCCutHeatmapEMCalwSys->SetBinContent(j, k, 0.);
                  simMCCutHeatmapEMCalwSys->SetBinContent(j, k, 0.);
               }
               if (dmCutterMC.IsDeadEMCal(1, i, j - 1, k - 1))
               {
                  simMCCutHeatmapEMCalw->SetBinContent(j, k, 0.);

                  realMCCutHeatmapEMCalwSys->SetBinContent(j, k, 0.);
                  simMCCutHeatmapEMCalwSys->SetBinContent(j, k, 0.);
               }
            }
         }

         DrawDeadmap(realHeatmapEMCalw, realCutHeatmapEMCalw,
                     "EMCalw" + std::to_string(i), "EMCalw" + std::to_string(i), 
                     "#it{Y}_{tower}", "#it{Z}_{tower}");

         DrawDeadmap(simCutHeatmapEMCalw, simMCCutHeatmapEMCalw,
                     "EMCalw" + std::to_string(i) + "_MC", "EMCalw" + std::to_string(i), 
                     "#it{Y}_{tower}", "#it{Z}_{tower}");

         GetUncertainty(realHeatmapEMCalwSys, simHeatmapEMCalwSys, 
                        realCutHeatmapEMCalwSys, simCutHeatmapEMCalwSys, 8,
                        "EMCalw" + std::to_string(i), "EMCalw" + std::to_string(i), 
                        "#it{Y}_{tower}", "#it{Z}_{tower}");

         systematicsOutputFile << 
            GetUncertainty(realHeatmapEMCalwSys, simHeatmapEMCalwSys, 
                           realMCCutHeatmapEMCalwSys, simMCCutHeatmapEMCalwSys, 8,
                           "EMCalw" + std::to_string(i) + "_MC", "EMCalw" + std::to_string(i), 
                           "#it{Y}_{tower}", "#it{Z}_{tower}");
         if (i < 3) systematicsOutputFile << " ";

         reweightEMCalw[i] = 
            realCutHeatmapEMCalwSys->Integral(1, realHeatmapEMCalw->GetXaxis()->GetNbins(),
                                              1, realHeatmapEMCalw->GetYaxis()->GetNbins())/
            realMCCutHeatmapEMCalwSys->Integral(1, realHeatmapEMCalw->GetXaxis()->GetNbins(),
                                                1, realHeatmapEMCalw->GetYaxis()->GetNbins());
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

   const std::string reweightsDir = "data/Parameters/MCReweights/" + runName;
   std::filesystem::create_directories(reweightsDir);
   std::ofstream reweightsFile(reweightsDir + "/ConstantScale.txt");

   reweightsFile << reweightDCe[0] << " " << reweightDCe[1] << " " <<
                    reweightDCw[0] << " " << reweightDCw[1] << std::endl <<
                    reweightPC1e << " " << reweightPC1w << std::endl <<
                    reweightPC2 << std::endl <<
                    reweightPC3e << " " << reweightPC3w << std::endl <<
                    reweightTOFe << " " << reweightTimingTOFe << std::endl <<
                    reweightTOFw << " " << reweightTimingTOFw << std::endl <<
                    reweightEMCale[0] << " " << reweightEMCale[1] << " " <<
                    reweightEMCale[2] << " " << reweightEMCale[3] << std::endl <<
                    reweightEMCalw[0] << " " << reweightEMCalw[1] << " " <<
                    reweightEMCalw[2] << " " << reweightEMCalw[3];
}

#endif /* DEAD_MAP_SYS_CPP */
