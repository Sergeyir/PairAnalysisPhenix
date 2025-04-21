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

void DeadMapSys::SetHistStyle(TH2F *hist, const std::string& title, 
                              const std::string& xTitle, const std::string &yTitle)
{
	hist->SetTitle(title.c_str());
	hist->GetXaxis()->SetTitle(xTitle.c_str());
	hist->GetYaxis()->SetTitle(yTitle.c_str());
	
	hist->SetTitleSize(0.05, "X");
	hist->SetTitleSize(0.05, "Y");

	hist->GetXaxis()->SetLabelSize(0.05);
	hist->GetYaxis()->SetLabelSize(0.05);

	hist->GetXaxis()->SetTitleOffset(0.9);
	hist->GetYaxis()->SetTitleOffset(1.25);
}

double DeadMapSys::GetNormRatio(const double ratio)
{
   if (ratio > 1.) return ratio;
   return 1./ratio;
}

double DeadMapSys::CalculateUncertaintyFromProj(TH1D *dataCutDistr, TH1D *simCutDistr)
{
	std::vector<double> dataNormIntegrals;
	std::vector<double> simNormIntegrals;
	std::vector<double> ratios;

	for (int i = 0; i < numberOfProjectionDivisions; i++)
	{
		dataNormIntegrals.push_back(0.);
		simNormIntegrals.push_back(0.);
		for (int j = 1; j < dataCutDistr->GetXaxis()->GetNbins(); j++)
		{
			if (j % 10 != i) continue;
			dataNormIntegrals.back() += dataCutDistr->GetBinContent(j);
			simNormIntegrals.back() += simCutDistr->GetBinContent(j);
		}
		
		simCutDistr->Scale(dataNormIntegrals[i]/simNormIntegrals[i]);
		ratios.push_back(GetNormRatio(dataCutDistr->Integral()/simCutDistr->Integral()));	
	}
	
	for (int i = 0; i < numberOfProjectionDivisions; i++)
	{
		const double dataIntegral = dataCutDistr->Integral(
			static_cast<int>(dataCutDistr->GetXaxis()->GetNbins()*i/numberOfProjectionDivisions), 
			 static_cast<int>(dataCutDistr->GetXaxis()->GetNbins()*(i+1)/numberOfProjectionDivisions));
		
		const double simIntegral = simCutDistr->Integral(
			static_cast<int>(dataCutDistr->GetXaxis()->GetNbins()*i/numberOfProjectionDivisions), 
			 static_cast<int>(dataCutDistr->GetXaxis()->GetNbins()*(i+1)/numberOfProjectionDivisions));
		if (simIntegral > 0 && dataIntegral > 0)
		{
			simCutDistr->Scale(dataIntegral/simIntegral);
			ratios.push_back(GetNormRatio(dataCutDistr->Integral()/simCutDistr->Integral()));	
		}
	}
	return CppTools::RMSFromCArray(&ratios[0], ratios.size()) - 1.;
}

double DeadMapSys::GetUncertaintyFromXProj(TH2F *realDistr, TH2F *simDistr, 
                                           TH2F *realCutDistr, TH2F* simCutDistr, 
                                           const std::string& detectorName, 
                                           const std::string& title, 
                                           const std::string& xTitle, const std::string& yTitle,
                                           const int rebinX)
{
	SetHistStyle(realDistr, title, xTitle, yTitle);
	SetHistStyle(simDistr, "MC " + title, xTitle, yTitle);
	SetHistStyle(realCutDistr, "Cut " + title, xTitle, yTitle);
	SetHistStyle(simCutDistr, "Cut MC " + title, xTitle, yTitle);
	
	TH1D *realCutDistrProj = realCutDistr->
      ProjectionX(((std::string) realCutDistr->GetName() + "_projX").c_str(), 
                  1, realCutDistr->GetYaxis()->GetNbins());
	TH1D *simCutDistrProj = simCutDistr->
      ProjectionX(((std::string) simCutDistr->GetName() + "_projX").c_str(), 
                  1, simCutDistr->GetYaxis()->GetNbins());

   realCutDistrProj->RebinX(rebinX);
   simCutDistrProj->RebinX(rebinX);

	realCutDistrProj->
      Scale(1./realCutDistrProj->Integral(1, realCutDistrProj->GetXaxis()->GetNbins()), "nosw2");
	simCutDistrProj->
      Scale(1./simCutDistrProj->Integral(1, simCutDistrProj->GetXaxis()->GetNbins()), "nosw2");

	realCutDistrProj->SetFillColorAlpha(kOrange - 4, 0.5);
	
	TLegend projLegend = TLegend(0.4, 0.78, 0.9, 0.88);
	projLegend.SetNColumns(2);
	
	projLegend.SetLineColorAlpha(0, 0);
	projLegend.SetFillColorAlpha(0, 0);
	
	realCutDistrProj->SetLineColor(kRed - 3);
	simCutDistrProj->SetLineColor(kAzure - 3);

	realCutDistrProj->SetLineWidth(2);
	simCutDistrProj->SetLineWidth(2);

	projLegend.AddEntry(realCutDistrProj, "real");
	projLegend.AddEntry(simCutDistrProj, "MC");
	
	realCutDistrProj->SetTitle((title + " projections real vs MC").c_str());
	realCutDistrProj->GetXaxis()->SetTitle(xTitle.c_str());

	realCutDistrProj->SetMaximum(realCutDistrProj->GetMaximum()*1.3);
	
	realCutDistrProj->GetXaxis()->SetLabelSize(0.05);
	realCutDistrProj->GetYaxis()->SetLabelSize(0.05);
	
	TCanvas canv("canv", "canv", 900, 450);
	canv.Divide(2);
	
	canv.cd(1);
	
	gPad->cd(1);
	gPad->SetLeftMargin(0.12);
	gPad->SetRightMargin(0.12);
	realDistr->Clone()->Draw("colz");

	canv.cd(2);
	
	gPad->cd(1);
	gPad->SetLeftMargin(0.12);
	gPad->SetRightMargin(0.12);
	realCutDistr->Clone()->Draw("colz");

   ROOTTools::PrintCanvas(&canv, outputDirDM + detectorName);

	TCanvas projCanv = TCanvas("projCanv", "canv", 1800, 600);

	projCanv.Divide(3, 1);

	projCanv.cd(1);
	gPad->SetLeftMargin(0.12);
	gPad->SetRightMargin(0.12);
	realCutDistr->Clone()->Draw("colz");

	projCanv.cd(2);
	gPad->SetLeftMargin(0.12);
	gPad->SetRightMargin(0.12);
	simCutDistr->Clone()->Draw("colz");

	projCanv.cd(3);
	gPad->SetLeftMargin(0.13);
	realCutDistrProj->Draw();
	simCutDistrProj->Draw("SAME HIST");

	projLegend.Draw();

   ROOTTools::PrintCanvas(&projCanv, outputDirSys + detectorName);
	
	const double uncertainty = CalculateUncertaintyFromProj(realCutDistrProj, simCutDistrProj);
	
	const double realLost = (1. - realCutDistr->Integral()/realDistr->Integral())*100.;
	const double simLost = (1. - simCutDistr->Integral()/simDistr->Integral())*100.;
	
	table.PrintRow(detectorName, 
		            CppTools::DtoStr(uncertainty*100., 3), 
		            "-", "-", 
		            CppTools::DtoStr(realLost, 3) + "%",
		            CppTools::DtoStr(simLost, 3) + "%");
	
	return uncertainty;
}

double DeadMapSys::GetUncertaintyFromXYProj(TH2F *realDistr, TH2F *simDistr, 
                                            TH2F *realCutDistr, TH2F* simCutDistr, 
                                            const std::string& detectorName, 
                                            const std::string& title, 
                                            const std::string& xTitle, const std::string& yTitle,
                                            const int rebinX, const int rebinY)
{
	SetHistStyle(realDistr, title, xTitle, yTitle);
	SetHistStyle(simDistr, "MC " + title, xTitle, yTitle);
	SetHistStyle(realCutDistr, "Cut " + title, xTitle, yTitle);
	SetHistStyle(simCutDistr, "Cut MC " + title, xTitle, yTitle);
	
	TH1D *realCutDistrProjX = realCutDistr->
      ProjectionX((title + "_real_x").c_str(), 1, realCutDistr->GetYaxis()->GetNbins());
	TH1D *realCutDistrProjY = realCutDistr->
      ProjectionY((title + "_real_y").c_str(), 1, realCutDistr->GetXaxis()->GetNbins());
	TH1D *simCutDistrProjX = simCutDistr->
      ProjectionX((title + "_sim_x").c_str(), 1, simCutDistr->GetYaxis()->GetNbins());
	TH1D *simCutDistrProjY = simCutDistr->
      ProjectionY((title + "_sim_y").c_str(), 1, simCutDistr->GetXaxis()->GetNbins());

   realCutDistrProjX->RebinX(rebinX);
   simCutDistrProjX->RebinX(rebinX);
   realCutDistrProjY->RebinX(rebinY);
   simCutDistrProjY->RebinX(rebinY);

	realCutDistrProjX->SetFillColorAlpha(kOrange-4, 0.5);
	realCutDistrProjY->SetFillColorAlpha(kOrange-4, 0.5);
	
	realCutDistrProjX->Scale(1./realCutDistrProjX->Integral(1, realCutDistrProjX->GetXaxis()->GetNbins()), "nosw2");
	simCutDistrProjX->Scale(1./simCutDistrProjX->Integral(1, simCutDistrProjX->GetXaxis()->GetNbins()), "nosw2");
	
	realCutDistrProjY->Scale(1./realCutDistrProjY->Integral(1, realCutDistrProjY->GetXaxis()->GetNbins()), "nosw2");
	simCutDistrProjY->Scale(1./simCutDistrProjY->Integral(1, simCutDistrProjY->GetXaxis()->GetNbins()), "nosw2");

	TLegend projXLegend = TLegend(0.4, 0.78, 0.9, 0.88);
	TLegend projYLegend = TLegend(0.4, 0.78, 0.9, 0.88);
	
	projXLegend.SetNColumns(2);
	projYLegend.SetNColumns(2);
	
	projXLegend.SetLineColorAlpha(0, 0);
	projXLegend.SetFillColorAlpha(0, 0);
	projYLegend.SetLineColorAlpha(0, 0);
	projYLegend.SetFillColorAlpha(0, 0);
	
	realCutDistrProjX->SetLineColor(kRed-3);
	realCutDistrProjY->SetLineColor(kRed-3);
	simCutDistrProjX->SetLineColor(kAzure-3);
	simCutDistrProjY->SetLineColor(kAzure-3);

	realCutDistrProjX->SetLineWidth(2);
	realCutDistrProjY->SetLineWidth(2);
	simCutDistrProjX->SetLineWidth(2);
	simCutDistrProjY->SetLineWidth(2);

	projXLegend.AddEntry(realCutDistrProjX, "real");
	projXLegend.AddEntry(simCutDistrProjX, "MC");

	projYLegend.AddEntry(realCutDistrProjY, "real");
	projYLegend.AddEntry(simCutDistrProjY, "MC");
	
	realCutDistrProjX->SetTitle((title + " X projections real vs MC").c_str());
	realCutDistrProjX->GetXaxis()->SetTitle(xTitle.c_str());

	realCutDistrProjY->SetTitle((title + " Y projections real vs MC").c_str());
	realCutDistrProjY->GetXaxis()->SetTitle(yTitle.c_str());

	realCutDistrProjX->SetMaximum(realCutDistrProjX->GetMaximum()*1.3);
	realCutDistrProjY->SetMaximum(realCutDistrProjY->GetMaximum()*1.3);

	realCutDistrProjX->GetXaxis()->SetLabelSize(0.05);
	realCutDistrProjX->GetYaxis()->SetLabelSize(0.05);

	realCutDistrProjY->GetXaxis()->SetLabelSize(0.05);
	realCutDistrProjY->GetYaxis()->SetLabelSize(0.05);
	
	TCanvas canv = TCanvas("canv", "canv", 900, 450);
	canv.Divide(2);

	canv.cd(1);
	gPad->SetLeftMargin(0.125);
	gPad->SetRightMargin(0.12);
	realDistr->Clone()->Draw("colz");
	
	canv.cd(2);
	gPad->SetLeftMargin(0.125);
	gPad->SetRightMargin(0.12);
	realCutDistr->Clone()->Draw("colz");

   ROOTTools::PrintCanvas(&canv, "output/Deadmaps/" + runName + "/" + detectorName);

	TCanvas projCanv = TCanvas("projCanv", "canv", 1800, 900);
	
	projCanv.cd();
	gPad->Divide(2, 2);
	
	projCanv.cd(1);
	gPad->SetLeftMargin(0.12);
	gPad->SetRightMargin(0.12);
	realCutDistr->Draw("colz");
	
	projCanv.cd(3);
	gPad->SetLeftMargin(0.12);
	gPad->SetRightMargin(0.12);
	simCutDistr->Draw("colz");
	
	projCanv.cd(2);
	gPad->SetLeftMargin(0.12);
	gPad->SetBottomMargin(0.11);
	realCutDistrProjX->Draw();
	simCutDistrProjX->Draw("SAME HIST");
	projXLegend.Draw();

	projCanv.cd(4);
	gPad->SetLeftMargin(0.12);
	gPad->SetBottomMargin(0.11);
	realCutDistrProjY->Draw();
	simCutDistrProjY->Draw("SAME HIST");
	projYLegend.Draw();

   ROOTTools::PrintCanvas(&projCanv, "output/Systematics/" + runName + "/" + detectorName);

	const double uncertaintyX = CalculateUncertaintyFromProj(realCutDistrProjX, simCutDistrProjX);
	const double uncertaintyY = CalculateUncertaintyFromProj(realCutDistrProjY, simCutDistrProjY);
	const double uncertainty = CppTools::RMS(uncertaintyX, uncertaintyY);

	const double realLost = (1. - realCutDistr->Integral()/realDistr->Integral())*100.;
	const double simLost = (1. - simCutDistr->Integral()/simDistr->Integral())*100.;

	table.PrintRow(detectorName, 
		CppTools::DtoStr(uncertainty*100., 3), 
		CppTools::DtoStr(uncertaintyX*100., 3), 
		CppTools::DtoStr(uncertaintyY*100., 3), 
		CppTools::DtoStr(realLost, 3) + "%",
		CppTools::DtoStr(simLost, 3) + "%");

	return uncertainty;
}

int main(int argc, char **argv)
{
   if (argc < 2 || argc > 3) 
   {
      std::string errMsg = "Expected 1 parameter while " + std::to_string(argc) + " ";
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

   system(("mkdir -p " + outputDirDM + " " + outputDirSys + " " + outputDirParameters).c_str());

   const std::string inputRealDataFileName = "data/Real/" + runName + "/SingleTrack/sum.root";
   const std::string inputSimDataFileName = "data/PostSim/" + runName + "/SingleTrack/all.root";

   CppTools::CheckInputFile(inputRealDataFileName);
   CppTools::CheckInputFile(inputSimDataFileName);

   inputRealDataFile = TFile::Open(inputRealDataFileName.c_str());
   inputSimDataFile = TFile::Open(inputSimDataFileName.c_str());

   const std::string detectorsConfiguration = 
      inputYAMLMain["detectors_configuration"].as<std::string>();

   dmCutter.Initialize(runName, detectorsConfiguration);

	table.Begin("Acceptance uncertainty");
	table.PrintHeader("detector", "sys", "sys x", "sys y", "data lost", "MC lost");

   if (detectorsConfiguration[0] == '1') // DC
   {
      TH2F *realHeatmapDCe0 = static_cast<TH2F *>(inputRealDataFile->Get("Heatmap: DCe, zDC>=0"));
      TH2F *realHeatmapDCe1 = static_cast<TH2F *>(inputRealDataFile->Get("Heatmap: DCe, zDC<0"));
      TH2F *realHeatmapDCw0 = static_cast<TH2F *>(inputRealDataFile->Get("Heatmap: DCw, zDC>=0"));
      TH2F *realHeatmapDCw1 = static_cast<TH2F *>(inputRealDataFile->Get("Heatmap: DCw, zDC<0"));

      TH2F *simHeatmapDCe0 = static_cast<TH2F *>(inputSimDataFile->Get("Heatmap: DCe, zDC>=0"));
      TH2F *simHeatmapDCe1 = static_cast<TH2F *>(inputSimDataFile->Get("Heatmap: DCe, zDC<0"));
      TH2F *simHeatmapDCw0 = static_cast<TH2F *>(inputSimDataFile->Get("Heatmap: DCw, zDC>=0"));
      TH2F *simHeatmapDCw1 = static_cast<TH2F *>(inputSimDataFile->Get("Heatmap: DCw, zDC<0"));

      CheckHists(realHeatmapDCe0, simHeatmapDCe0, "Heatmap: DCe, zDC>=0");
      CheckHists(realHeatmapDCe1, simHeatmapDCe1, "Heatmap: DCe, zDC<0");
      CheckHists(realHeatmapDCw0, simHeatmapDCw0, "Heatmap: DCw, zDC>=0");
      CheckHists(realHeatmapDCw1, simHeatmapDCw1, "Heatmap: DCw, zDC<0");

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
            if (dmCutter.IsDeadDC(1, 1., realHeatmapDCe0->GetXaxis()->GetBinCenter(i),
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
            if (dmCutter.IsDeadDC(1, -1., realHeatmapDCe1->GetXaxis()->GetBinCenter(i),
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
            if (dmCutter.IsDeadDC(0, 1., realHeatmapDCw0->GetXaxis()->GetBinCenter(i),
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
            if (dmCutter.IsDeadDC(0, -1., realHeatmapDCw1->GetXaxis()->GetBinCenter(i),
                                  realHeatmapDCw1->GetYaxis()->GetBinCenter(j)))
            {
               realCutHeatmapDCw1->SetBinContent(i, j, 0.);
               simCutHeatmapDCw1->SetBinContent(i, j, 0.);
            }
         }
      }

      std::ofstream systematicsOutputFile(outputDirParameters + "DC.txt");
      systematicsOutputFile << 
         GetUncertaintyFromXProj(realHeatmapDCe0, simHeatmapDCe0, 
                                 realCutHeatmapDCe0, simCutHeatmapDCe0, 
                                 "DCe0", "DC east, zDC>=0", "board", "#alpha", 2) << " " <<
         GetUncertaintyFromXProj(realHeatmapDCe1, simHeatmapDCe1, 
                                 realCutHeatmapDCe1, simCutHeatmapDCe1, 
                                 "DCe1", "DC east, zDC<0", "board", "#alpha", 2) << " " <<
         GetUncertaintyFromXProj(realHeatmapDCw0, simHeatmapDCw0, 
                                 realCutHeatmapDCw0, simCutHeatmapDCw0, 
                                 "DCw0", "DC west, zDC>=0", "board", "#alpha", 2) << " " <<
         GetUncertaintyFromXProj(realHeatmapDCw1, simHeatmapDCw1, 
                                 realCutHeatmapDCw1, simCutHeatmapDCw1, 
                                 "DCw1", "DC west, zDC<0", "board", "#alpha", 2);
   }
   else
   {
      std::ofstream systematicsOutputFile(outputDirParameters + "DC.txt");
      systematicsOutputFile << 0 << " " << 0 << " " << 0 << " " << 0;
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
            if (dmCutter.IsDeadPC1(1, realHeatmapPC1e->GetXaxis()->GetBinCenter(i),
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
            if (dmCutter.IsDeadPC1(0, realHeatmapPC1w->GetXaxis()->GetBinCenter(i),
                                   realHeatmapPC1w->GetYaxis()->GetBinCenter(j)))
            {
               realCutHeatmapPC1w->SetBinContent(i, j, 0.);
               simCutHeatmapPC1w->SetBinContent(i, j, 0.);
            }
         }
      }

      std::ofstream systematicsOutputFile(outputDirParameters + "PC1.txt");
      systematicsOutputFile << 
         GetUncertaintyFromXYProj(realHeatmapPC1e, simHeatmapPC1e, 
                                  realCutHeatmapPC1e, simCutHeatmapPC1e, 
                                  "PC1e", "PC1 east", "z_{PC1}", "y_{PC1}", 2) << " " <<
         GetUncertaintyFromXYProj(realHeatmapPC1w, simHeatmapPC1w, 
                                  realCutHeatmapPC1w, simCutHeatmapPC1w, 
                                  "PC1w", "PC1 west", "z_{PC1}", "y_{PC1}", 2);
   }
   else
   {
      std::ofstream systematicsOutputFile(outputDirParameters + "PC1.txt");
      systematicsOutputFile << 0 << " " << 0;
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

      std::ofstream systematicsOutputFile(outputDirParameters + "PC2.txt");
      systematicsOutputFile << 
         GetUncertaintyFromXYProj(realHeatmapPC2, simHeatmapPC2, 
                                  realCutHeatmapPC2, simCutHeatmapPC2, 
                                  "PC2", "PC2", "z_{PC2}", "y_{PC2}");
   }
   else
   {
      std::ofstream systematicsOutputFile(outputDirParameters + "PC2.txt");
      systematicsOutputFile << 0;
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
            if (dmCutter.IsDeadPC3(1, realHeatmapPC3e->GetXaxis()->GetBinCenter(i),
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
            if (dmCutter.IsDeadPC3(0, realHeatmapPC3w->GetXaxis()->GetBinCenter(i),
                                   realHeatmapPC3w->GetYaxis()->GetBinCenter(j)))
            {
               realCutHeatmapPC3w->SetBinContent(i, j, 0.);
               simCutHeatmapPC3w->SetBinContent(i, j, 0.);
            }
         }
      }

      std::ofstream systematicsOutputFile(outputDirParameters + "PC3.txt");
      systematicsOutputFile << 
         GetUncertaintyFromXYProj(realHeatmapPC3e, simHeatmapPC3e, 
                                  realCutHeatmapPC3e, simCutHeatmapPC3e, 
                                  "PC3e", "PC3 east", "z_{PC3}", "y_{PC3}", 2) << " " <<
         GetUncertaintyFromXYProj(realHeatmapPC3w, simHeatmapPC3w, 
                                  realCutHeatmapPC3w, simCutHeatmapPC3w, 
                                  "PC3w", "PC3 west", "z_{PC3}", "y_{PC3}", 2);
   }
   else
   {
      std::ofstream systematicsOutputFile(outputDirParameters + "PC3.txt");
      systematicsOutputFile << 0 << " " << 0;
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
            if (dmCutter.IsDeadTOFe(realHeatmapTOFe->GetXaxis()->GetBinCenter(i),
                                   realHeatmapTOFe->GetYaxis()->GetBinCenter(j)))
            {
               realCutHeatmapTOFe->SetBinContent(i, j, 0.);
               simCutHeatmapTOFe->SetBinContent(i, j, 0.);
            }
         }
      }

      std::ofstream systematicsOutputFile(outputDirParameters + "TOFe.txt");
      systematicsOutputFile << 
         GetUncertaintyFromXYProj(realHeatmapTOFe, simHeatmapTOFe, 
                                  realCutHeatmapTOFe, simCutHeatmapTOFe, 
                                  "TOFe", "TOFe", "y_{TOFe}", "z_{TOFe}");
   }
   else
   {
      std::ofstream systematicsOutputFile(outputDirParameters + "TOFe.txt");
      systematicsOutputFile << 0;
   }

   if (detectorsConfiguration[5] == '1') // TOFw
   {
      TH2F *realHeatmapTOFw0 = 
         static_cast<TH2F *>(inputRealDataFile->Get("Heatmap: TOFw, ptofy<100"));
      TH2F *realHeatmapTOFw1 = 
         static_cast<TH2F *>(inputRealDataFile->Get("Heatmap: TOFw, ptofy>100"));

      TH2F *simHeatmapTOFw0 = 
         static_cast<TH2F *>(inputSimDataFile->Get("Heatmap: TOFw, ptofy<100"));
      TH2F *simHeatmapTOFw1 = 
         static_cast<TH2F *>(inputSimDataFile->Get("Heatmap: TOFw, ptofy>100"));

      CheckHists(realHeatmapTOFw0, simHeatmapTOFw0, "Heatmap: TOFw, ptofy<100");
      CheckHists(realHeatmapTOFw1, simHeatmapTOFw1, "Heatmap: TOFw, ptofy>100");

      TH2F *realCutHeatmapTOFw0 = static_cast<TH2F *>(realHeatmapTOFw0->Clone());
      TH2F *realCutHeatmapTOFw1 = static_cast<TH2F *>(realHeatmapTOFw1->Clone());

      TH2F *simCutHeatmapTOFw0 = static_cast<TH2F *>(simHeatmapTOFw0->Clone());
      TH2F *simCutHeatmapTOFw1 = static_cast<TH2F *>(simHeatmapTOFw1->Clone());

      for (int i = 1; i <= realHeatmapTOFw0->GetXaxis()->GetNbins(); i++)
      {
         for (int j = 1; j <= realHeatmapTOFw0->GetYaxis()->GetNbins(); j++)
         {
            if (dmCutter.IsDeadTOFw(realHeatmapTOFw0->GetXaxis()->GetBinCenter(i),
                                    realHeatmapTOFw0->GetYaxis()->GetBinCenter(j)))
            {
               realCutHeatmapTOFw0->SetBinContent(i, j, 0.);
               simCutHeatmapTOFw0->SetBinContent(i, j, 0.);
            }
         }
      }
      for (int i = 1; i <= realHeatmapTOFw1->GetXaxis()->GetNbins(); i++)
      {
         for (int j = 1; j <= realHeatmapTOFw1->GetYaxis()->GetNbins(); j++)
         {
            if (dmCutter.IsDeadTOFw(realHeatmapTOFw1->GetXaxis()->GetBinCenter(i),
                                    realHeatmapTOFw1->GetYaxis()->GetBinCenter(j)))
            {
               realCutHeatmapTOFw1->SetBinContent(i, j, 0.);
               simCutHeatmapTOFw1->SetBinContent(i, j, 0.);
            }
         }
      }

      std::ofstream systematicsOutputFile(outputDirParameters + "TOFw.txt");
      systematicsOutputFile << 
         CppTools::RMS(GetUncertaintyFromXYProj(realHeatmapTOFw0, simHeatmapTOFw0, 
                                                realCutHeatmapTOFw0, simCutHeatmapTOFw0, 
                                                "TOFw0", "TOFw, y_{TOFw}<100", 
                                                "y_{TOFw}", "z_{TOFw}"),
                       GetUncertaintyFromXYProj(realHeatmapTOFw1, simHeatmapTOFw1, 
                                                realCutHeatmapTOFw1, simCutHeatmapTOFw1, 
                                                "TOFw1", "TOFw, y_{TOFw}>100", 
                                                "y_{TOFw}", "z_{TOFw}"));
   }
   else
   {
      std::ofstream systematicsOutputFile(outputDirParameters + "TOFw.txt");
      systematicsOutputFile << 0;
   }

   if (detectorsConfiguration[6] == '1') // EMCal
   {
      std::ofstream systematicsOutputFile(outputDirParameters + "EMCal.txt");

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
               if (dmCutter.IsDeadEMCal(1, i, realHeatmapEMCale->GetXaxis()->GetBinCenter(j),
                                        realHeatmapEMCale->GetYaxis()->GetBinCenter(k)))
               {
                  realCutHeatmapEMCale->SetBinContent(j, k, 0.);
                  simCutHeatmapEMCale->SetBinContent(j, k, 0.);
               }
            }
         }

         systematicsOutputFile << 
            GetUncertaintyFromXYProj(realHeatmapEMCale, simHeatmapEMCale, 
                                     realCutHeatmapEMCale, simCutHeatmapEMCale, 
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
               if (dmCutter.IsDeadEMCal(0, i, realHeatmapEMCalw->GetXaxis()->GetBinCenter(j),
                                        realHeatmapEMCalw->GetYaxis()->GetBinCenter(k)))
               {
                  realCutHeatmapEMCalw->SetBinContent(j, k, 0.);
                  simCutHeatmapEMCalw->SetBinContent(j, k, 0.);
               }
            }
         }

         systematicsOutputFile << 
            GetUncertaintyFromXYProj(realHeatmapEMCalw, simHeatmapEMCalw, 
                                     realCutHeatmapEMCalw, simCutHeatmapEMCalw, 
                                     "EMCalw" + std::to_string(i), "EMCalw" + std::to_string(i), 
                                     "y_{tower}", "z_{tower}");
         if (i < 3) systematicsOutputFile << " ";
      }
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
   }

	table.End();
}

#endif /* DEAD_MAP_SYS_CPP */
