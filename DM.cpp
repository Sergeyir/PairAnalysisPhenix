#include "DM.hpp"

using namespace Run14HeAu200Cuts;

struct
{
	std::string runName = "Run14HeAu200";
	Table table{6};

	TFile *dataFile = TFile::Open(("data/Real/" + runName + "/sum.root").c_str());
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
                                 const std::string& xTitle, const std::string& yTitle)
{   
	SetHistStyle(dataDistr, title, xTitle, yTitle);
	SetHistStyle(simDistr, "MC " + title, xTitle, yTitle);
	SetHistStyle(dataCutDistr, "Cut " + title, xTitle, yTitle);
	SetHistStyle(simCutDistr, "Cut MC " + title, xTitle, yTitle);
	
	TH1F *dataCutDistr_proj = (TH1F*) dataCutDistr->ProjectionX(((std::string) dataCutDistr->GetName() + "_data").c_str(), 1, dataCutDistr->GetYaxis()->GetNbins());
	TH1F *simCutDistr_proj = (TH1F*) simCutDistr->ProjectionX(((std::string) simCutDistr->GetName() + "_sim").c_str(), 1, simCutDistr->GetYaxis()->GetNbins());

	dataCutDistr_proj->Scale(1./dataCutDistr_proj->Integral(
		1, dataCutDistr_proj->GetXaxis()->GetNbins()), "nosw2");
	simCutDistr_proj->Scale(1./simCutDistr_proj->Integral(
		1, simCutDistr_proj->GetXaxis()->GetNbins()), "nosw2");

	dataCutDistr_proj->SetFillColorAlpha(kOrange-4, 0.5);
	
	TLegend proj_legend = TLegend(0.4, 0.78, 0.9, 0.88);
	proj_legend.SetNColumns(2);
	
	proj_legend.SetLineColorAlpha(0, 0);
	proj_legend.SetFillColorAlpha(0, 0);
	
	dataCutDistr_proj->SetLineColor(kRed-3);
	simCutDistr_proj->SetLineColor(kAzure-3);

	dataCutDistr_proj->SetLineWidth(2);
	simCutDistr_proj->SetLineWidth(2);

	proj_legend.AddEntry(dataCutDistr_proj, "data");
	proj_legend.AddEntry(simCutDistr_proj, "MC");
	
	dataCutDistr_proj->SetTitle((title + " projections data vs MC").c_str());
	dataCutDistr_proj->GetXaxis()->SetTitle(xTitle.c_str());

	dataCutDistr_proj->SetMaximum(dataCutDistr_proj->GetMaximum()*1.3);
	
	dataCutDistr_proj->GetXaxis()->SetLabelSize(0.05);
	dataCutDistr_proj->GetYaxis()->SetLabelSize(0.05);
	
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

	TCanvas proj_canv = TCanvas("proj_canv", "canv", 1800, 600);

	proj_canv.Divide(3, 1);

	proj_canv.cd(1);
	gPad->SetLeftMargin(0.12);
	gPad->SetRightMargin(0.12);
	dataCutDistr->Clone()->Draw("colz");

	proj_canv.cd(2);
	gPad->SetLeftMargin(0.12);
	gPad->SetRightMargin(0.12);
	simCutDistr->Clone()->Draw("colz");

	proj_canv.cd(3);
	gPad->SetLeftMargin(0.13);
	dataCutDistr_proj->Draw();
	simCutDistr_proj->Draw("SAME HIST");

	proj_legend.Draw();

	PrintCanvas(&proj_canv, "output/Systematics/" + Par.runName + "/" + detectorName);
	
	const double uncertainty = CalculateUncertaintyFromProj(dataCutDistr_proj, simCutDistr_proj);
	
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
                                 const std::string& xTitle, const std::string& yTitle)
{
	SetHistStyle(dataDistr, title, xTitle, yTitle);
	SetHistStyle(simDistr, "MC " + title, xTitle, yTitle);
	SetHistStyle(dataCutDistr, "Cut " + title, xTitle, yTitle);
	SetHistStyle(simCutDistr, "Cut MC " + title, xTitle, yTitle);
	
	TH1F *dataCutDistr_proj_x = (TH1F*) dataCutDistr->ProjectionX(
		(title + "_data_x").c_str(), 1, dataCutDistr->GetXaxis()->GetNbins());
	TH1F *dataCutDistr_proj_y = (TH1F*) dataCutDistr->ProjectionY(
		(title + "_data_y").c_str(), 1, dataCutDistr->GetYaxis()->GetNbins());
	TH1F *simCutDistr_proj_x = (TH1F*) simCutDistr->ProjectionX(
		(title + "_sim_x").c_str(), 1, simCutDistr->GetXaxis()->GetNbins());
	TH1F *simCutDistr_proj_y = (TH1F*) simCutDistr->ProjectionY(
		(title + "_sim_y").c_str(), 1, simCutDistr->GetYaxis()->GetNbins());

	dataCutDistr_proj_x->SetFillColorAlpha(kOrange-4, 0.5);
	dataCutDistr_proj_y->SetFillColorAlpha(kOrange-4, 0.5);
	
	dataCutDistr_proj_x->Scale(1./dataCutDistr_proj_x->Integral(1, dataCutDistr_proj_x->GetXaxis()->GetNbins()), "nosw2");
	simCutDistr_proj_x->Scale(1./simCutDistr_proj_x->Integral(1, simCutDistr_proj_x->GetXaxis()->GetNbins()), "nosw2");
	
	dataCutDistr_proj_y->Scale(1./dataCutDistr_proj_y->Integral(1, dataCutDistr_proj_y->GetXaxis()->GetNbins()), "nosw2");
	simCutDistr_proj_y->Scale(1./simCutDistr_proj_y->Integral(1, simCutDistr_proj_y->GetXaxis()->GetNbins()), "nosw2");

	TLegend proj_x_legend = TLegend(0.4, 0.78, 0.9, 0.88);
	TLegend proj_y_legend = TLegend(0.4, 0.78, 0.9, 0.88);
	
	proj_x_legend.SetNColumns(2);
	proj_y_legend.SetNColumns(2);
	
	proj_x_legend.SetLineColorAlpha(0, 0);
	proj_x_legend.SetFillColorAlpha(0, 0);
	proj_y_legend.SetLineColorAlpha(0, 0);
	proj_y_legend.SetFillColorAlpha(0, 0);
	
	dataCutDistr_proj_x->SetLineColor(kRed-3);
	dataCutDistr_proj_y->SetLineColor(kRed-3);
	simCutDistr_proj_x->SetLineColor(kAzure-3);
	simCutDistr_proj_y->SetLineColor(kAzure-3);

	dataCutDistr_proj_x->SetLineWidth(2);
	dataCutDistr_proj_y->SetLineWidth(2);
	simCutDistr_proj_x->SetLineWidth(2);
	simCutDistr_proj_y->SetLineWidth(2);

	proj_x_legend.AddEntry(dataCutDistr_proj_x, "data");
	proj_x_legend.AddEntry(simCutDistr_proj_x, "MC");

	proj_y_legend.AddEntry(dataCutDistr_proj_y, "data");
	proj_y_legend.AddEntry(simCutDistr_proj_y, "MC");
	
	dataCutDistr_proj_x->SetTitle((title + " X projections data vs MC").c_str());
	dataCutDistr_proj_x->GetXaxis()->SetTitle(xTitle.c_str());

	dataCutDistr_proj_y->SetTitle((title + " Y projections data vs MC").c_str());
	dataCutDistr_proj_y->GetXaxis()->SetTitle(yTitle.c_str());

	dataCutDistr_proj_x->SetMaximum(dataCutDistr_proj_x->GetMaximum()*1.3);
	dataCutDistr_proj_y->SetMaximum(dataCutDistr_proj_y->GetMaximum()*1.3);

	dataCutDistr_proj_x->GetXaxis()->SetLabelSize(0.05);
	dataCutDistr_proj_x->GetYaxis()->SetLabelSize(0.05);

	dataCutDistr_proj_y->GetXaxis()->SetLabelSize(0.05);
	dataCutDistr_proj_y->GetYaxis()->SetLabelSize(0.05);
	
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

	TCanvas proj_canv = TCanvas("proj_canv", "canv", 1800, 900);
	
	proj_canv.cd();
	gPad->Divide(2, 2);
	
	proj_canv.cd(1);
	gPad->SetLeftMargin(0.12);
	gPad->SetRightMargin(0.12);
	dataCutDistr->Draw("colz");
	
	proj_canv.cd(3);
	gPad->SetLeftMargin(0.12);
	gPad->SetRightMargin(0.12);
	simCutDistr->Draw("colz");
	
	proj_canv.cd(2);
	gPad->SetLeftMargin(0.12);
	gPad->SetBottomMargin(0.11);
	dataCutDistr_proj_x->Draw();
	simCutDistr_proj_x->Draw("SAME HIST");
	proj_x_legend.Draw();

	proj_canv.cd(4);
	gPad->SetLeftMargin(0.12);
	gPad->SetBottomMargin(0.11);
	dataCutDistr_proj_y->Draw();
	simCutDistr_proj_y->Draw("SAME HIST");
	proj_y_legend.Draw();

	PrintCanvas(&proj_canv, "output/Systematics/" + Par.runName + "/" + detectorName);

	const double uncertainty_x = CalculateUncertaintyFromProj(dataCutDistr_proj_x, simCutDistr_proj_x);
	const double uncertainty_y = CalculateUncertaintyFromProj(dataCutDistr_proj_y, simCutDistr_proj_y);
	const double uncertainty = RMS(uncertainty_x, uncertainty_y);

	const double data_lost = (1. - dataCutDistr->Integral()/dataDistr->Integral())*100.;
	const double sim_lost = (1. - simCutDistr->Integral()/simDistr->Integral())*100.;

	Par.table.PrintRow(detectorName, 
		DtoStr(uncertainty*100., 3), 
		DtoStr(uncertainty_x*100., 3), 
		DtoStr(uncertainty_y*100., 3), 
		DtoStr(data_lost, 3) + "%",
		DtoStr(sim_lost, 3) + "%");

	return uncertainty;
}

double GetUncertainty(const std::string& histName, const std::string &detectorName, 
                      const std::string &title, const std::string xTitle, 
                      const std::string yTitle, bool sysCalc2Proj, 
                      bool (*cut_func)(const double, const double))
{
   TH2F *dataDistr = static_cast<TH2F *>(Par.dataFile->Get(histName.c_str()));
   TH2F *simDistr = static_cast<TH2F *>(Par.simFile->Get(histName.c_str()));
   
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
                                       detectorName, title, xTitle, yTitle);
   }
   return CalculateUncertainty2Proj(dataDistr, simDistr, dataCutDistr, simCutDistr, 
                                    detectorName, title, xTitle, yTitle);
}

double  GetUncertainty(const std::string& histName, const std::string &detectorName, 
                       const std::string &title, const std::string xTitle, 
                       const std::string yTitle, bool sysCalc2Proj, 
                       bool (*cut_func)(const double, const double, const double), 
                       const double auxVal)
{
   TH2F *dataDistr = static_cast<TH2F *>(Par.dataFile->Get(histName.c_str()));
   TH2F *simDistr = static_cast<TH2F *>(Par.simFile->Get(histName.c_str()));
   
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
                                       detectorName, title, xTitle, yTitle);
   }
   return CalculateUncertainty2Proj(dataDistr, simDistr, dataCutDistr, simCutDistr, 
                                    detectorName, title, xTitle, yTitle);
}

double GetUncertainty(const std::string& histName, const std::string &detectorName, 
                      const std::string &title, const std::string xTitle, 
                      const std::string yTitle, bool sysCalc2Proj, 
                      bool (*cut_func)(const double, const double, const double, const double), 
                      const double auxVal1, const double auxVal2)
{
   TH2F *dataDistr = static_cast<TH2F *>(Par.dataFile->Get(histName.c_str()));
   TH2F *simDistr = static_cast<TH2F *>(Par.simFile->Get(histName.c_str()));
   
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
                                       detectorName, title, xTitle, yTitle);
   }
   return CalculateUncertainty2Proj(dataDistr, simDistr, dataCutDistr, simCutDistr, 
                                    detectorName, title, xTitle, yTitle);
}

double GetUncertainty(const std::string& histName, const std::string &detectorName, 
                      const std::string &title, const std::string xTitle, 
                      const std::string yTitle, bool sysCalc2Proj, 
                      bool (*cut_func)(const double, const double, const int, 
                                       const double, const double), 
                      const double auxVal1, const double auxVal2, const int auxVal3)
{
   TH2F *dataDistr = static_cast<TH2F *>(Par.dataFile->Get(histName.c_str()));
   TH2F *simDistr = static_cast<TH2F *>(Par.simFile->Get(histName.c_str()));
   
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
                                       detectorName, title, xTitle, yTitle);
   }
   return CalculateUncertainty2Proj(dataDistr, simDistr, dataCutDistr, simCutDistr, 
                                    detectorName, title, xTitle, yTitle);
}

void DM()
{
	gROOT->SetBatch(kTRUE);
	gStyle->SetOptStat(0);
	
	gErrorIgnoreLevel = kWarning;

	CheckInputFile("data/Real/" + Par.runName + "/sum.root");
	CheckInputFile("data/PostSim/" + Par.runName + "/Heatmaps/");

	system(("mkdir -p data/Deadmaps/" + Par.runName).c_str());
	system(("mkdir -p data/Systematics/" + Par.runName).c_str());
	system(("mkdir -p output/Deadmaps/" + Par.runName).c_str());
	system(("mkdir -p output/Systematics/" + Par.runName).c_str());

	Par.table.Begin("Acceptance uncertainty");
	Par.table.PrintHeader("detector", "sys", "sys x", "sys y", "data lost", "MC lost");

   const double sysDCe0 = 
      GetUncertainty("Heatmap: DCe, zed>=0", "DCe0", "DCe, z#geq0", "board", "#alpha", false,
                     &IsDeadDC, 2., 1.);
   const double sysDCe1 = 
      GetUncertainty("Heatmap: DCe, zed<0", "DCe1", "DCe, z<0", "board", "#alpha", false,
                     &IsDeadDC, 2., -1.);
   const double sysDCw0 = 
      GetUncertainty("Heatmap: DCw, zed>=0", "DCw0", "DCw, z#geq0", "board", "#alpha", false,
                     &IsDeadDC, 1., 1.);
   const double sysDCw1 = 
      GetUncertainty("Heatmap: DCw, zed<0", "DCw1", "DCw, z<0", "board", "#alpha", false,
                     &IsDeadDC, 1., -1.);

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
