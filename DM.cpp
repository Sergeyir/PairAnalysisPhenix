#include "DM.hpp"

using namespace Run14HeAu200Cuts;

struct
{
	std::string runName = "Run14HeAu200";
	Table table{6};

	TFile *dataFile = TFile::Open(("data/Real/" + runName + "/SingleTrack/sum.root").c_str());
	TFile *simFile = TFile::Open(("data/PostSim/" + runName + "/Heatmaps/all.root").c_str());

	const int sysNDiv = 4;
} Par;

void CheckHists(TH2F *hist1, TH2F *hist2)
{
   if (hist1->GetXaxis()->GetNbins() != hist2->GetXaxis()->GetNbins() ||
       fabs(hist1->GetXaxis()->GetBinLowEdge(1) - hist1->GetXaxis()->GetBinLowEdge(1)) > 1e-6 ||
       fabs(hist1->GetXaxis()->GetBinUpEdge(hist1->GetXaxis()->GetNbins()) - 
            hist1->GetXaxis()->GetBinUpEdge(hist2->GetXaxis()->GetNbins())) > 1e-6) 
   {
      PrintError("X axis bins are inconsistent for histograms named" + 
                 static_cast<std::string>(hist1->GetName()));
   }
   if (hist1->GetYaxis()->GetNbins() != hist2->GetYaxis()->GetNbins() ||
       fabs(hist1->GetYaxis()->GetBinLowEdge(1) - hist1->GetYaxis()->GetBinLowEdge(1)) > 1e-6 ||
       fabs(hist1->GetYaxis()->GetBinUpEdge(hist1->GetYaxis()->GetNbins()) - 
            hist1->GetYaxis()->GetBinUpEdge(hist2->GetYaxis()->GetNbins())) > 1e-6) 
   {
      PrintError("Y axis bins are inconsistent for histograms named" + 
                 static_cast<std::string>(hist1->GetName()));
   }
}

void SetHistStyle(TH2F *hist, const std::string& title, 
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

double GetNormRatio(const double ratio)
{
   if (ratio > 1.) return ratio;
   return 1./ratio;
}

double CalculateUncertaintyFromProj(TH1F *dataCutDistr, TH1F *simCutDistr)
{
	std::vector<double> dataNormIntegrals;
	std::vector<double> simNormIntegrals;
	std::vector<double> ratios;

	for (int i = 0; i < 10; i++)
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
	
	for (int i = 0; i < Par.sysNDiv; i++)
	{
		const double dataIntegral = dataCutDistr->Integral(
			static_cast<int>(dataCutDistr->GetXaxis()->GetNbins()*i/Par.sysNDiv), 
			 static_cast<int>(dataCutDistr->GetXaxis()->GetNbins()*(i+1)/Par.sysNDiv));
		
		const double simIntegral = simCutDistr->Integral(
			static_cast<int>(dataCutDistr->GetXaxis()->GetNbins()*i/Par.sysNDiv), 
			 static_cast<int>(dataCutDistr->GetXaxis()->GetNbins()*(i+1)/Par.sysNDiv));
		if (simIntegral > 0 && dataIntegral > 0)
		{
			simCutDistr->Scale(dataIntegral/simIntegral);
			ratios.push_back(GetNormRatio(dataCutDistr->Integral()/simCutDistr->Integral()));	
		}
	}
	return RMSv(ratios) - 1.;
}

double CalculateUncertainty1Proj(TH2F *dataDistr, TH2F *simDistr, 
                                 TH2F *dataCutDistr, TH2F* simCutDistr, 
                                 const std::string& detectorName, const std::string& title, 
                                 const std::string& xTitle, const std::string& yTitle,
                                 const int rebinX)
{   
	SetHistStyle(dataDistr, title, xTitle, yTitle);
	SetHistStyle(simDistr, "MC " + title, xTitle, yTitle);
	SetHistStyle(dataCutDistr, "Cut " + title, xTitle, yTitle);
	SetHistStyle(simCutDistr, "Cut MC " + title, xTitle, yTitle);
	
	TH1F *dataCutDistrProj = (TH1F*) dataCutDistr->ProjectionX(((std::string) dataCutDistr->GetName() + "_data").c_str(), 1, dataCutDistr->GetYaxis()->GetNbins());
	TH1F *simCutDistrProj = (TH1F*) simCutDistr->ProjectionX(((std::string) simCutDistr->GetName() + "_sim").c_str(), 1, simCutDistr->GetYaxis()->GetNbins());

   dataCutDistrProj->RebinX(rebinX);
   simCutDistrProj->RebinX(rebinX);

	dataCutDistrProj->Scale(1./dataCutDistrProj->Integral(
		1, dataCutDistrProj->GetXaxis()->GetNbins()), "nosw2");
	simCutDistrProj->Scale(1./simCutDistrProj->Integral(
		1, simCutDistrProj->GetXaxis()->GetNbins()), "nosw2");

	dataCutDistrProj->SetFillColorAlpha(kOrange-4, 0.5);
	
	TLegend projLegend = TLegend(0.4, 0.78, 0.9, 0.88);
	projLegend.SetNColumns(2);
	
	projLegend.SetLineColorAlpha(0, 0);
	projLegend.SetFillColorAlpha(0, 0);
	
	dataCutDistrProj->SetLineColor(kRed-3);
	simCutDistrProj->SetLineColor(kAzure-3);

	dataCutDistrProj->SetLineWidth(2);
	simCutDistrProj->SetLineWidth(2);

	projLegend.AddEntry(dataCutDistrProj, "data");
	projLegend.AddEntry(simCutDistrProj, "MC");
	
	dataCutDistrProj->SetTitle((title + " projections data vs MC").c_str());
	dataCutDistrProj->GetXaxis()->SetTitle(xTitle.c_str());

	dataCutDistrProj->SetMaximum(dataCutDistrProj->GetMaximum()*1.3);
	
	dataCutDistrProj->GetXaxis()->SetLabelSize(0.05);
	dataCutDistrProj->GetYaxis()->SetLabelSize(0.05);
	
	TCanvas canv = TCanvas("canv", "canv", 900, 450);
	canv.Divide(2);
	
	canv.cd(1);
	
	gPad->cd(1);
	gPad->SetLeftMargin(0.12);
	gPad->SetRightMargin(0.12);
	dataDistr->Clone()->Draw("colz");

	canv.cd(2);
	
	gPad->cd(1);
	gPad->SetLeftMargin(0.12);
	gPad->SetRightMargin(0.12);
	dataCutDistr->Clone()->Draw("colz");

	PrintCanvas(&canv, "output/Deadmaps/" + Par.runName + "/" + detectorName);

	TCanvas projCanv = TCanvas("projCanv", "canv", 1800, 600);

	projCanv.Divide(3, 1);

	projCanv.cd(1);
	gPad->SetLeftMargin(0.12);
	gPad->SetRightMargin(0.12);
	dataCutDistr->Clone()->Draw("colz");

	projCanv.cd(2);
	gPad->SetLeftMargin(0.12);
	gPad->SetRightMargin(0.12);
	simCutDistr->Clone()->Draw("colz");

	projCanv.cd(3);
	gPad->SetLeftMargin(0.13);
	dataCutDistrProj->Draw();
	simCutDistrProj->Draw("SAME HIST");

	projLegend.Draw();

	PrintCanvas(&projCanv, "output/Systematics/" + Par.runName + "/" + detectorName);
	
	const double uncertainty = CalculateUncertaintyFromProj(dataCutDistrProj, simCutDistrProj);
	
	const double data_lost = (1. - dataCutDistr->Integral()/dataDistr->Integral())*100.;
	const double sim_lost = (1. - simCutDistr->Integral()/simDistr->Integral())*100.;
	
	Par.table.PrintRow(detectorName, 
		DtoStr(uncertainty*100., 3), 
		"-", "-", 
		DtoStr(data_lost, 3) + "%",
		DtoStr(sim_lost, 3) + "%");
	
	return uncertainty;
}

double CalculateUncertainty2Proj(TH2F *dataDistr, TH2F *simDistr, 
                                 TH2F *dataCutDistr, TH2F* simCutDistr, 
                                 const std::string& detectorName, const std::string& title, 
                                 const std::string& xTitle, const std::string& yTitle,
                                 const int rebinX, const int rebinY)
{
	SetHistStyle(dataDistr, title, xTitle, yTitle);
	SetHistStyle(simDistr, "MC " + title, xTitle, yTitle);
	SetHistStyle(dataCutDistr, "Cut " + title, xTitle, yTitle);
	SetHistStyle(simCutDistr, "Cut MC " + title, xTitle, yTitle);
	
	TH1F *dataCutDistrProjX = (TH1F*) dataCutDistr->ProjectionX(
		(title + "_data_x").c_str(), 1, dataCutDistr->GetXaxis()->GetNbins());
	TH1F *dataCutDistrProjY = (TH1F*) dataCutDistr->ProjectionY(
		(title + "_data_y").c_str(), 1, dataCutDistr->GetYaxis()->GetNbins());
	TH1F *simCutDistrProjX = (TH1F*) simCutDistr->ProjectionX(
		(title + "_sim_x").c_str(), 1, simCutDistr->GetXaxis()->GetNbins());
	TH1F *simCutDistrProjY = (TH1F*) simCutDistr->ProjectionY(
		(title + "_sim_y").c_str(), 1, simCutDistr->GetYaxis()->GetNbins());

   dataCutDistrProjX->RebinX(rebinX);
   simCutDistrProjX->RebinX(rebinX);
   dataCutDistrProjY->RebinX(rebinY);
   simCutDistrProjY->RebinX(rebinY);

	dataCutDistrProjX->SetFillColorAlpha(kOrange-4, 0.5);
	dataCutDistrProjY->SetFillColorAlpha(kOrange-4, 0.5);
	
	dataCutDistrProjX->Scale(1./dataCutDistrProjX->Integral(1, dataCutDistrProjX->GetXaxis()->GetNbins()), "nosw2");
	simCutDistrProjX->Scale(1./simCutDistrProjX->Integral(1, simCutDistrProjX->GetXaxis()->GetNbins()), "nosw2");
	
	dataCutDistrProjY->Scale(1./dataCutDistrProjY->Integral(1, dataCutDistrProjY->GetXaxis()->GetNbins()), "nosw2");
	simCutDistrProjY->Scale(1./simCutDistrProjY->Integral(1, simCutDistrProjY->GetXaxis()->GetNbins()), "nosw2");

	TLegend projXLegend = TLegend(0.4, 0.78, 0.9, 0.88);
	TLegend projYLegend = TLegend(0.4, 0.78, 0.9, 0.88);
	
	projXLegend.SetNColumns(2);
	projYLegend.SetNColumns(2);
	
	projXLegend.SetLineColorAlpha(0, 0);
	projXLegend.SetFillColorAlpha(0, 0);
	projYLegend.SetLineColorAlpha(0, 0);
	projYLegend.SetFillColorAlpha(0, 0);
	
	dataCutDistrProjX->SetLineColor(kRed-3);
	dataCutDistrProjY->SetLineColor(kRed-3);
	simCutDistrProjX->SetLineColor(kAzure-3);
	simCutDistrProjY->SetLineColor(kAzure-3);

	dataCutDistrProjX->SetLineWidth(2);
	dataCutDistrProjY->SetLineWidth(2);
	simCutDistrProjX->SetLineWidth(2);
	simCutDistrProjY->SetLineWidth(2);

	projXLegend.AddEntry(dataCutDistrProjX, "data");
	projXLegend.AddEntry(simCutDistrProjX, "MC");

	projYLegend.AddEntry(dataCutDistrProjY, "data");
	projYLegend.AddEntry(simCutDistrProjY, "MC");
	
	dataCutDistrProjX->SetTitle((title + " X projections data vs MC").c_str());
	dataCutDistrProjX->GetXaxis()->SetTitle(xTitle.c_str());

	dataCutDistrProjY->SetTitle((title + " Y projections data vs MC").c_str());
	dataCutDistrProjY->GetXaxis()->SetTitle(yTitle.c_str());

	dataCutDistrProjX->SetMaximum(dataCutDistrProjX->GetMaximum()*1.3);
	dataCutDistrProjY->SetMaximum(dataCutDistrProjY->GetMaximum()*1.3);

	dataCutDistrProjX->GetXaxis()->SetLabelSize(0.05);
	dataCutDistrProjX->GetYaxis()->SetLabelSize(0.05);

	dataCutDistrProjY->GetXaxis()->SetLabelSize(0.05);
	dataCutDistrProjY->GetYaxis()->SetLabelSize(0.05);
	
	TCanvas canv = TCanvas("canv", "canv", 900, 450);
	canv.Divide(2);

	canv.cd(1);
	gPad->SetLeftMargin(0.125);
	gPad->SetRightMargin(0.12);
	dataDistr->Clone()->Draw("colz");
	
	canv.cd(2);
	gPad->SetLeftMargin(0.125);
	gPad->SetRightMargin(0.12);
	dataCutDistr->Clone()->Draw("colz");

	PrintCanvas(&canv, "output/Deadmaps/" + Par.runName + "/" + detectorName);

	TCanvas projCanv = TCanvas("projCanv", "canv", 1800, 900);
	
	projCanv.cd();
	gPad->Divide(2, 2);
	
	projCanv.cd(1);
	gPad->SetLeftMargin(0.12);
	gPad->SetRightMargin(0.12);
	dataCutDistr->Draw("colz");
	
	projCanv.cd(3);
	gPad->SetLeftMargin(0.12);
	gPad->SetRightMargin(0.12);
	simCutDistr->Draw("colz");
	
	projCanv.cd(2);
	gPad->SetLeftMargin(0.12);
	gPad->SetBottomMargin(0.11);
	dataCutDistrProjX->Draw();
	simCutDistrProjX->Draw("SAME HIST");
	projXLegend.Draw();

	projCanv.cd(4);
	gPad->SetLeftMargin(0.12);
	gPad->SetBottomMargin(0.11);
	dataCutDistrProjY->Draw();
	simCutDistrProjY->Draw("SAME HIST");
	projYLegend.Draw();

	PrintCanvas(&projCanv, "output/Systematics/" + Par.runName + "/" + detectorName);

	const double uncertaintyX = CalculateUncertaintyFromProj(dataCutDistrProjX, simCutDistrProjX);
	const double uncertaintyY = CalculateUncertaintyFromProj(dataCutDistrProjY, simCutDistrProjY);
	const double uncertainty = RMS(uncertaintyX, uncertaintyY);

	const double data_lost = (1. - dataCutDistr->Integral()/dataDistr->Integral())*100.;
	const double sim_lost = (1. - simCutDistr->Integral()/simDistr->Integral())*100.;

	Par.table.PrintRow(detectorName, 
		DtoStr(uncertainty*100., 3), 
		DtoStr(uncertaintyX*100., 3), 
		DtoStr(uncertaintyY*100., 3), 
		DtoStr(data_lost, 3) + "%",
		DtoStr(sim_lost, 3) + "%");

	return uncertainty;
}

double GetUncertainty(const std::string& histName, const std::string &detectorName, 
                      const std::string &title, const std::string xTitle, 
                      const std::string yTitle, bool sysCalc2Proj, 
                      bool (*cut_func)(const double, const double),
                      const int rebinX, const int rebinY)
{
   TH2F *dataDistr = static_cast<TH2F *>(Par.dataFile->Get(histName.c_str()));
   TH2F *simDistr = static_cast<TH2F *>(Par.simFile->Get(histName.c_str()));

   if (!dataDistr || dataDistr == NULL) PrintError(histName + " does not exist in real data file");
   if (!dataDistr || simDistr == NULL) PrintError(histName + " does not exist in sim data file");
   
   CheckHists(dataDistr, simDistr);
   
   TH2F *dataCutDistr = static_cast<TH2F *>(dataDistr->Clone());
   TH2F *simCutDistr = static_cast<TH2F *>(simDistr->Clone());

	for (int i = 1; i <= dataDistr->GetXaxis()->GetNbins(); i++)
	{
		for (int j = 1; j <= dataDistr->GetYaxis()->GetNbins(); j++)
		{
			const double xVal = dataDistr->GetXaxis()->GetBinCenter(i);
			const double yVal = dataDistr->GetYaxis()->GetBinCenter(j);
         
			if (cut_func(xVal, yVal)) 
         {
            dataCutDistr->SetBinContent(i, j, 0.);
            simCutDistr->SetBinContent(i, j, 0.);
         }
		}
	}
   if (!sysCalc2Proj) 
   {
      return CalculateUncertainty1Proj(dataDistr, simDistr, dataCutDistr, simCutDistr, 
                                       detectorName, title, xTitle, yTitle, rebinX);
   }
   return CalculateUncertainty2Proj(dataDistr, simDistr, dataCutDistr, simCutDistr, 
                                    detectorName, title, xTitle, yTitle, rebinX, rebinY);
}

double  GetUncertainty(const std::string& histName, const std::string &detectorName, 
                       const std::string &title, const std::string xTitle, 
                       const std::string yTitle, bool sysCalc2Proj, 
                       bool (*cut_func)(const double, const double, const double), 
                       const double auxVal, const int rebinX, const int rebinY)
{
   TH2F *dataDistr = static_cast<TH2F *>(Par.dataFile->Get(histName.c_str()));
   TH2F *simDistr = static_cast<TH2F *>(Par.simFile->Get(histName.c_str()));
   
   if (!dataDistr || dataDistr == NULL) PrintError(histName + " does not exist in real data file");
   if (!dataDistr || simDistr == NULL) PrintError(histName + " does not exist in sim data file");
   
   CheckHists(dataDistr, simDistr);
   
   TH2F *dataCutDistr = static_cast<TH2F *>(dataDistr->Clone());
   TH2F *simCutDistr = static_cast<TH2F *>(simDistr->Clone());

	for (int i = 1; i <= dataDistr->GetXaxis()->GetNbins(); i++)
	{
		for (int j = 1; j <= dataDistr->GetYaxis()->GetNbins(); j++)
		{
			const double xVal = dataDistr->GetXaxis()->GetBinCenter(i);
			const double yVal = dataDistr->GetYaxis()->GetBinCenter(j);
			
			if (cut_func(auxVal, xVal, yVal)) 
			{
            dataCutDistr->SetBinContent(i, j, 0.);
            simCutDistr->SetBinContent(i, j, 0.);
			}
		}
	}
   if (!sysCalc2Proj) 
   {
      return CalculateUncertainty1Proj(dataDistr, simDistr, dataCutDistr, simCutDistr, 
                                       detectorName, title, xTitle, yTitle, rebinX);
   }
   return CalculateUncertainty2Proj(dataDistr, simDistr, dataCutDistr, simCutDistr, 
                                    detectorName, title, xTitle, yTitle, rebinX, rebinY);
}

double GetUncertainty(const std::string& histName, const std::string &detectorName, 
                      const std::string &title, const std::string xTitle, 
                      const std::string yTitle, bool sysCalc2Proj, 
                      bool (*cut_func)(const double, const double, const double, const double), 
                      const double auxVal1, const double auxVal2, 
                      const int rebinX, const int rebinY)
{
   TH2F *dataDistr = static_cast<TH2F *>(Par.dataFile->Get(histName.c_str()));
   TH2F *simDistr = static_cast<TH2F *>(Par.simFile->Get(histName.c_str()));

   if (!dataDistr || dataDistr == NULL) PrintError(histName + " does not exist in real data file");
   if (!dataDistr || simDistr == NULL) PrintError(histName + " does not exist in sim data file");
   
   CheckHists(dataDistr, simDistr);
   
   TH2F *dataCutDistr = static_cast<TH2F *>(dataDistr->Clone());
   TH2F *simCutDistr = static_cast<TH2F *>(simDistr->Clone());

	for (int i = 1; i <= dataDistr->GetXaxis()->GetNbins(); i++)
	{
		for (int j = 1; j <= dataDistr->GetYaxis()->GetNbins(); j++)
		{
			const double xVal = dataDistr->GetXaxis()->GetBinCenter(i);
			const double yVal = dataDistr->GetYaxis()->GetBinCenter(j);
			
			if (cut_func(auxVal1, auxVal2, xVal, yVal)) 
			{
            dataCutDistr->SetBinContent(i, j, 0.);
            simCutDistr->SetBinContent(i, j, 0.);
			}
		}
	}
   if (!sysCalc2Proj) 
   {
      return CalculateUncertainty1Proj(dataDistr, simDistr, dataCutDistr, simCutDistr, 
                                       detectorName, title, xTitle, yTitle, rebinX);
   }
   return CalculateUncertainty2Proj(dataDistr, simDistr, dataCutDistr, simCutDistr, 
                                    detectorName, title, xTitle, yTitle, rebinX, rebinY);
}

double GetUncertainty(const std::string& histName, const std::string &detectorName, 
                      const std::string &title, const std::string xTitle, 
                      const std::string yTitle, bool sysCalc2Proj, 
                      bool (*cut_func)(const double, const double, const int, 
                                       const double, const double), 
                      const double auxVal1, const double auxVal2, const int auxVal3,
                      const int rebinX, const int rebinY)
{
   TH2F *dataDistr = static_cast<TH2F *>(Par.dataFile->Get(histName.c_str()));
   TH2F *simDistr = static_cast<TH2F *>(Par.simFile->Get(histName.c_str()));

   if (!dataDistr || dataDistr == NULL) PrintError(histName + " does not exist in real data file");
   if (!dataDistr || simDistr == NULL) PrintError(histName + " does not exist in sim data file");
   
   CheckHists(dataDistr, simDistr);
   
   TH2F *dataCutDistr = static_cast<TH2F *>(dataDistr->Clone());
   TH2F *simCutDistr = static_cast<TH2F *>(simDistr->Clone());

	for (int i = 1; i <= dataDistr->GetXaxis()->GetNbins(); i++)
	{
		for (int j = 1; j <= dataDistr->GetYaxis()->GetNbins(); j++)
		{
			const double xVal = dataDistr->GetXaxis()->GetBinCenter(i);
			const double yVal = dataDistr->GetYaxis()->GetBinCenter(j);
			
			if (cut_func(auxVal1, auxVal2, auxVal3, xVal, yVal)) 
			{
            dataCutDistr->SetBinContent(i, j, 0.);
            simCutDistr->SetBinContent(i, j, 0.);
			}
		}
	}
   
   if (!sysCalc2Proj) 
   {
      return CalculateUncertainty1Proj(dataDistr, simDistr, dataCutDistr, simCutDistr, 
                                       detectorName, title, xTitle, yTitle, rebinX);
   }
   return CalculateUncertainty2Proj(dataDistr, simDistr, dataCutDistr, simCutDistr, 
                                    detectorName, title, xTitle, yTitle, rebinX, rebinY);
}

void DM()
{
	gROOT->SetBatch(kTRUE);
	gStyle->SetOptStat(0);
	
	gErrorIgnoreLevel = kWarning;

	CheckInputFile("data/Real/" + Par.runName + "/SingleTrack/sum.root");
	CheckInputFile("data/PostSim/" + Par.runName + "/Heatmaps/");

	system(("mkdir -p data/Deadmaps/" + Par.runName).c_str());
	system(("mkdir -p data/Systematics/" + Par.runName).c_str());
	system(("mkdir -p output/Deadmaps/" + Par.runName).c_str());
	system(("mkdir -p output/Systematics/" + Par.runName).c_str());

	Par.table.Begin("Acceptance uncertainty");
	Par.table.PrintHeader("detector", "sys", "sys x", "sys y", "data lost", "MC lost");

   const double sysDCe0 = 
      GetUncertainty("Heatmap: DCe, zDC>=0", "DCe0", "DCe, z#geq0", "board", "#alpha", false,
                     &IsDeadDC, 2., 1., 2);
   const double sysDCe1 = 
      GetUncertainty("Heatmap: DCe, zDC<0", "DCe1", "DCe, z<0", "board", "#alpha", false,
                     &IsDeadDC, 2., -1., 2);
   const double sysDCw0 = 
      GetUncertainty("Heatmap: DCw, zDC>=0", "DCw0", "DCw, z#geq0", "board", "#alpha", false,
                     &IsDeadDC, 1., 1., 2);
   const double sysDCw1 = 
      GetUncertainty("Heatmap: DCw, zDC<0", "DCw1", "DCw, z<0", "board", "#alpha", false,
                     &IsDeadDC, 1., -1., 2);

   /*
	TH2F *dceast0_data_cut = CutDeadAreas((TH2F *) dceast0_data->Clone(), &IsDeadDC, 2., 1.);
	TH2F *dceast1_data_cut = CutDeadAreas((TH2F *) dceast1_data->Clone(), &IsDeadDC, 2., -1.);
	TH2F *dcwest0_data_cut = CutDeadAreas((TH2F *) dcwest0_data->Clone(), &IsDeadDC, 1., 1.);
	TH2F *dcwest1_data_cut = CutDeadAreas((TH2F *) dcwest1_data->Clone(), &IsDeadDC, 1., -1.);

	//PC1
	TH2F *pc1e_data_cut = CutDeadAreas((TH2F *) pc1e_data->Clone(), &IsDeadPC1, 2.);
	TH2F *pc1w_data_cut = CutDeadAreas((TH2F *) pc1w_data->Clone(), &IsDeadPC1, 1.);

	//PC2
	TH2F *pc2_data_cut = CutDeadAreas((TH2F *) pc2_data->Clone(), &IsDeadPC2);

	//PC3
	TH2F *pc3e_data_cut = CutDeadAreas((TH2F *) pc3e_data->Clone(), &IsDeadPC3, 2.);
	TH2F *pc3w_data_cut = CutDeadAreas((TH2F *) pc3w_data->Clone(), &IsDeadPC3, 1.);

	//EMCal
	TH2F *emcale0_pos_data_cut = CutDeadAreas((TH2F *) emcale0_pos_data->Clone(), &IsDeadEMCal, 2., 1., 0);
	TH2F *emcale1_pos_data_cut = CutDeadAreas((TH2F *) emcale1_pos_data->Clone(), &IsDeadEMCal, 2., 1., 1);
	TH2F *emcale2_pos_data_cut = CutDeadAreas((TH2F *) emcale2_pos_data->Clone(), &IsDeadEMCal, 2., 1., 2);
	TH2F *emcale3_pos_data_cut = CutDeadAreas((TH2F *) emcale3_pos_data->Clone(), &IsDeadEMCal, 2., 1., 3);
	TH2F *emcale0_neg_data_cut = CutDeadAreas((TH2F *) emcale0_neg_data->Clone(), &IsDeadEMCal, 2., -1., 0);
	TH2F *emcale1_neg_data_cut = CutDeadAreas((TH2F *) emcale1_neg_data->Clone(), &IsDeadEMCal, 2., -1., 1);
	TH2F *emcale2_neg_data_cut = CutDeadAreas((TH2F *) emcale2_neg_data->Clone(), &IsDeadEMCal, 2., -1., 2);
	TH2F *emcale3_neg_data_cut = CutDeadAreas((TH2F *) emcale3_neg_data->Clone(), &IsDeadEMCal, 2., -1., 3);
	TH2F *emcalw0_pos_data_cut = CutDeadAreas((TH2F *) emcalw0_pos_data->Clone(), &IsDeadEMCal, 1., 1., 0);
	TH2F *emcalw1_pos_data_cut = CutDeadAreas((TH2F *) emcalw1_pos_data->Clone(), &IsDeadEMCal, 1., 1., 1);
	TH2F *emcalw2_pos_data_cut = CutDeadAreas((TH2F *) emcalw2_pos_data->Clone(), &IsDeadEMCal, 1., 1., 2);
	TH2F *emcalw3_pos_data_cut = CutDeadAreas((TH2F *) emcalw3_pos_data->Clone(), &IsDeadEMCal, 1., 1., 3);
	TH2F *emcalw0_neg_data_cut = CutDeadAreas((TH2F *) emcalw0_neg_data->Clone(), &IsDeadEMCal, 1., -1., 0);
	TH2F *emcalw1_neg_data_cut = CutDeadAreas((TH2F *) emcalw1_neg_data->Clone(), &IsDeadEMCal, 1., -1., 1);
	TH2F *emcalw2_neg_data_cut = CutDeadAreas((TH2F *) emcalw2_neg_data->Clone(), &IsDeadEMCal, 1., -1., 2);
	TH2F *emcalw3_neg_data_cut = CutDeadAreas((TH2F *) emcalw3_neg_data->Clone(), &IsDeadEMCal, 1., -1., 3);

	//TOFe
	TH2F *tofe0_data_cut = CutDeadAreas((TH2F *) tofe0_data->Clone(), &IsDeadTOFe, 1.);
	TH2F *tofe1_data_cut = CutDeadAreas((TH2F *) tofe1_data->Clone(), &IsDeadTOFe, -1.);
	
	//TOFw
	TH2F *tofw0_data_cut = CutDeadAreas((TH2F *) tofw0_data->Clone(), &IsDeadTOFw, 1.);
	TH2F *tofw1_data_cut = CutDeadAreas((TH2F *) tofw1_data->Clone(), &IsDeadTOFw, -1.);
   */
   
	Par.table.End();

	//system(("mkdir -p ../../sim_analysis/input/Systematics/" + Par.runName).c_str());
			
   /*
	WriteFile("../../sim_analysis/input/Systematics/" + Par.runName + "/acceptance.txt", 
		dceast0_sys, dceast1_sys, dcwest0_sys, dcwest1_sys,
		pc1e_sys, pc1w_sys, pc2_sys, pc3e_sys, pc3w_sys,
		tofe0_sys, tofe1_sys, tofw0_sys, tofw1_sys,
		emcale0_pos_sys, emcale0_neg_sys, emcale1_pos_sys, emcale1_neg_sys, 
		emcale2_pos_sys, emcale2_neg_sys, emcale3_pos_sys, emcale3_neg_sys,
		emcalw0_pos_sys, emcalw0_neg_sys, emcalw1_pos_sys, emcalw1_neg_sys,
		emcalw2_pos_sys, emcalw2_neg_sys, emcalw3_pos_sys, emcalw3_neg_sys);
      */
}
