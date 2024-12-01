#include "IOTools.hpp"
#include "DeadAreasCuts.h"
#include "Table.hpp"
#include "StrTools.hpp"
#include "TCanvasPrinter.hpp"

using namespace CutsRun7AuAu200MinBias;

struct
{
	std::string run_name = "Run7AuAu200";
	std::string data_file_name = "sum.root";
	std::string sim_file_name = "heatmaps.root";
	const int sys_ndiv = 4;

	Table<std::string, std::string, std::string, std::string, std::string, std::string> table;
} Par;

double GetSystUncertainty(TH1F *data_hist, TH1F *sim_hist)
{
	std::vector<double> data_norm_integrals;
	std::vector<double> sim_norm_integrals;
	std::vector<double> ratios;

	for (int i = 0; i < 10; i++)
	{
		data_norm_integrals.push_back(0.);
		sim_norm_integrals.push_back(0.);
		for (int j = 1; j < data_hist->GetXaxis()->GetNbins(); j++)
		{
			if (j % 10 != i) continue;
			data_norm_integrals.back() += data_hist->GetBinContent(j);
			sim_norm_integrals.back() += sim_hist->GetBinContent(j);
		}
		
		sim_hist->Scale(data_norm_integrals[i]/sim_norm_integrals[i]);
		ratios.push_back(GetNormRatio(data_hist->Integral()/sim_hist->Integral()));	
	}
	
	for (int i = 0; i < Par.sys_ndiv; i++)
	{
		const double data_int = data_hist->Integral(
			static_cast<int>(data_hist->GetXaxis()->GetNbins()*i/Par.sys_ndiv), 
			 static_cast<int>(data_hist->GetXaxis()->GetNbins()*(i+1)/Par.sys_ndiv));
		
		const double sim_int = sim_hist->Integral(
			static_cast<int>(data_hist->GetXaxis()->GetNbins()*i/Par.sys_ndiv), 
			 static_cast<int>(data_hist->GetXaxis()->GetNbins()*(i+1)/Par.sys_ndiv));
		if (sim_int > 0 && data_int > 0)
		{
			sim_hist->Scale(data_int/sim_int);
			ratios.push_back(GetNormRatio(data_hist->Integral()/sim_hist->Integral()));	
		}
	}
	return RMSv(ratios)-1.;
}

TH2F* CutDeadAreas(TH2F *hist, bool (*cut_func)(const double, const double))
{
	for (int i = 1; i <= hist->GetXaxis()->GetNbins(); i++)
	{
		for (int j = 1; j <= hist->GetYaxis()->GetNbins(); j++)
		{
			const double val_x = hist->GetXaxis()->GetBinCenter(i);
			const double val_y = hist->GetYaxis()->GetBinCenter(j);
			
			if (cut_func(val_x, val_y)) 
			{
				hist->SetBinContent(i, j, 0.);
			}
		}
	}
	return hist;
}

TH2F* CutDeadAreas(TH2F *hist, bool (*cut_func)(const double, const double, const double), const double aux_val)
{
	for (int i = 1; i <= hist->GetXaxis()->GetNbins(); i++)
	{
		for (int j = 1; j <= hist->GetYaxis()->GetNbins(); j++)
		{
			const double val_x = hist->GetXaxis()->GetBinCenter(i);
			const double val_y = hist->GetYaxis()->GetBinCenter(j);
			
			if (cut_func(aux_val, val_x, val_y)) 
			{
				hist->SetBinContent(i, j, 0.);
			}
		}
	}
	return hist;
}

TH2F* CutDeadAreas(TH2F *hist, bool (*cut_func)(const double, const double, const double, const double), const double aux_val1, const double aux_val2)
{
	for (int i = 1; i <= hist->GetXaxis()->GetNbins(); i++)
	{
		for (int j = 1; j <= hist->GetYaxis()->GetNbins(); j++)
		{
			const double val_x = hist->GetXaxis()->GetBinCenter(i);
			const double val_y = hist->GetYaxis()->GetBinCenter(j);
			
			if (cut_func(aux_val1, aux_val2, val_x, val_y)) 
			{
				hist->SetBinContent(i, j, 0.);
			}
		}
	}
	return hist;
}

TH2F* CutDeadAreas(TH2F *hist, bool (*cut_func)(const double, const double, const int, const double, const double), const double aux_val1, const double aux_val2, const int aux_val3)
{
	for (int i = 1; i <= hist->GetXaxis()->GetNbins(); i++)
	{
		for (int j = 1; j <= hist->GetYaxis()->GetNbins(); j++)
		{
			const double val_x = hist->GetXaxis()->GetBinCenter(i);
			const double val_y = hist->GetYaxis()->GetBinCenter(j);
			
			if (cut_func(aux_val1, aux_val2, aux_val3, val_x, val_y)) 
			{
				hist->SetBinContent(i, j, 0.);
			}
		}
	}
	return hist;
}

void SetHistStyle(TH2F *hist, std::string title, std::string x_title, std::string y_title)
{
	hist->SetTitle(title.c_str());
	hist->GetXaxis()->SetTitle(x_title.c_str());
	hist->GetYaxis()->SetTitle(y_title.c_str());
	
	hist->SetTitleSize(0.05, "X");
	hist->SetTitleSize(0.05, "Y");

	hist->GetXaxis()->SetLabelSize(0.05);
	hist->GetYaxis()->SetLabelSize(0.05);

	hist->GetXaxis()->SetTitleOffset(0.9);
	hist->GetYaxis()->SetTitleOffset(1.25);
}

double GetUncertainty1Proj(TH2F *data_raw_hist, TH2F *sim_raw_hist, TH2F *data_hist, TH2F* sim_hist, std::string name, std::string title, std::string x_title, std::string y_title, const int rebin_x = 1)
{
	SetHistStyle(data_raw_hist, title, x_title, y_title);
	SetHistStyle(sim_raw_hist, "MC " + title, x_title, y_title);
	SetHistStyle(data_hist, "Cut " + title, x_title, y_title);
	SetHistStyle(sim_hist, "Cut MC " + title, x_title, y_title);
	
	TH1F *data_hist_proj = (TH1F*) data_hist->ProjectionX((name + "_data").c_str(), 1, data_hist->GetYaxis()->GetNbins());
	TH1F *sim_hist_proj = (TH1F*) sim_hist->ProjectionX((name + "_sim").c_str(), 1, sim_hist->GetYaxis()->GetNbins());

	data_hist_proj->Scale(1./data_hist_proj->Integral(
		1, data_hist_proj->GetXaxis()->GetNbins()), "nosw2");
	sim_hist_proj->Scale(1./sim_hist_proj->Integral(
		1, sim_hist_proj->GetXaxis()->GetNbins()), "nosw2");

	data_hist_proj->Rebin(rebin_x);
	sim_hist_proj->Rebin(rebin_x);

	data_hist_proj->SetFillColorAlpha(kOrange-4, 0.5);
	
	TLegend proj_legend = TLegend(0.4, 0.78, 0.9, 0.88);
	proj_legend.SetNColumns(2);
	
	proj_legend.SetLineColorAlpha(0, 0);
	proj_legend.SetFillColorAlpha(0, 0);
	
	data_hist_proj->SetLineColor(kRed-3);
	sim_hist_proj->SetLineColor(kAzure-3);

	data_hist_proj->SetLineWidth(2);
	sim_hist_proj->SetLineWidth(2);

	proj_legend.AddEntry(data_hist_proj, "data");
	proj_legend.AddEntry(sim_hist_proj, "MC");
	
	data_hist_proj->SetTitle((title + " projections data vs MC").c_str());
	data_hist_proj->GetXaxis()->SetTitle(x_title.c_str());

	data_hist_proj->SetMaximum(data_hist_proj->GetMaximum()*1.3);
	
	data_hist_proj->GetXaxis()->SetLabelSize(0.05);
	data_hist_proj->GetYaxis()->SetLabelSize(0.05);
	
	TCanvas canv = TCanvas("canv", "canv", 900, 450);
	canv.Divide(2);
	
	canv.cd(1);
	
	gPad->cd(1);
	gPad->SetLeftMargin(0.12);
	gPad->SetRightMargin(0.12);
	data_raw_hist->Clone()->Draw("colz");

	canv.cd(2);
	
	gPad->cd(1);
	gPad->SetLeftMargin(0.12);
	gPad->SetRightMargin(0.12);
	data_hist->Clone()->Draw("colz");

	PrintCanvas(&canv, "../output/deadmaps/" + Par.run_name + "/" + name);

	TCanvas proj_canv = TCanvas("proj_canv", "canv", 1800, 900);

	proj_canv.Divide(3, 1);

	proj_canv.cd(1);
	gPad->SetLeftMargin(0.12);
	gPad->SetRightMargin(0.12);
	data_hist->Clone()->Draw("colz");

	proj_canv.cd(2);
	gPad->SetLeftMargin(0.12);
	gPad->SetRightMargin(0.12);
	sim_hist->Clone()->Draw("colz");

	proj_canv.cd(3);
	gPad->SetLeftMargin(0.13);
	data_hist_proj->Draw();
	sim_hist_proj->Draw("SAME HIST");

	proj_legend.Draw();

	PrintCanvas(&proj_canv, "../output/systematics/" + Par.run_name + "/" + name);
	
	const double uncertainty = GetSystUncertainty(data_hist_proj, sim_hist_proj);
	
	const double data_lost = (1. - data_hist->Integral()/data_raw_hist->Integral())*100.;
	const double sim_lost = (1. - sim_hist->Integral()/sim_raw_hist->Integral())*100.;
	
	Par.table.PrintRow(name, 
		DtoStr(uncertainty*100., 3), 
		"-", "-", 
		DtoStr(data_lost, 3) + "%",
		DtoStr(sim_lost, 3) + "%");
	
	return uncertainty;
}

double GetUncertainty2Proj(TH2F *data_raw_hist, TH2F *sim_raw_hist, TH2F *data_hist, TH2F* sim_hist, std::string name, std::string title, std::string x_title, std::string y_title, const int rebin_x = 1, const int rebin_y = 1)
{
	SetHistStyle(data_raw_hist, title, x_title, y_title);
	SetHistStyle(sim_raw_hist, "MC " + title, x_title, y_title);
	SetHistStyle(data_hist, "Cut " + title, x_title, y_title);
	SetHistStyle(sim_hist, "Cut MC " + title, x_title, y_title);
	
	TH1F *data_hist_proj_x = (TH1F*) data_hist->ProjectionX(
		(title + "_data_x").c_str(), 1, data_hist->GetXaxis()->GetNbins());
	TH1F *data_hist_proj_y = (TH1F*) data_hist->ProjectionY(
		(title + "_data_y").c_str(), 1, data_hist->GetYaxis()->GetNbins());
	TH1F *sim_hist_proj_x = (TH1F*) sim_hist->ProjectionX(
		(title + "_sim_x").c_str(), 1, sim_hist->GetXaxis()->GetNbins());
	TH1F *sim_hist_proj_y = (TH1F*) sim_hist->ProjectionY(
		(title + "_sim_y").c_str(), 1, sim_hist->GetYaxis()->GetNbins());

	data_hist_proj_x->Rebin(rebin_x);
	sim_hist_proj_x->Rebin(rebin_x);
	data_hist_proj_y->Rebin(rebin_y);
	sim_hist_proj_y->Rebin(rebin_y);

	data_hist_proj_x->SetFillColorAlpha(kOrange-4, 0.5);
	data_hist_proj_y->SetFillColorAlpha(kOrange-4, 0.5);
	
	data_hist_proj_x->Scale(1./data_hist_proj_x->Integral(1, data_hist_proj_x->GetXaxis()->GetNbins()), "nosw2");
	sim_hist_proj_x->Scale(1./sim_hist_proj_x->Integral(1, sim_hist_proj_x->GetXaxis()->GetNbins()), "nosw2");
	
	data_hist_proj_y->Scale(1./data_hist_proj_y->Integral(1, data_hist_proj_y->GetXaxis()->GetNbins()), "nosw2");
	sim_hist_proj_y->Scale(1./sim_hist_proj_y->Integral(1, sim_hist_proj_y->GetXaxis()->GetNbins()), "nosw2");

	TLegend proj_x_legend = TLegend(0.4, 0.78, 0.9, 0.88);
	TLegend proj_y_legend = TLegend(0.4, 0.78, 0.9, 0.88);
	
	proj_x_legend.SetNColumns(2);
	proj_y_legend.SetNColumns(2);
	
	proj_x_legend.SetLineColorAlpha(0, 0);
	proj_x_legend.SetFillColorAlpha(0, 0);
	proj_y_legend.SetLineColorAlpha(0, 0);
	proj_y_legend.SetFillColorAlpha(0, 0);
	
	data_hist_proj_x->SetLineColor(kRed-3);
	data_hist_proj_y->SetLineColor(kRed-3);
	sim_hist_proj_x->SetLineColor(kAzure-3);
	sim_hist_proj_y->SetLineColor(kAzure-3);

	data_hist_proj_x->SetLineWidth(2);
	data_hist_proj_y->SetLineWidth(2);
	sim_hist_proj_x->SetLineWidth(2);
	sim_hist_proj_y->SetLineWidth(2);

	proj_x_legend.AddEntry(data_hist_proj_x, "data");
	proj_x_legend.AddEntry(sim_hist_proj_x, "MC");

	proj_y_legend.AddEntry(data_hist_proj_y, "data");
	proj_y_legend.AddEntry(sim_hist_proj_y, "MC");
	
	data_hist_proj_x->SetTitle((title + " X projections data vs MC").c_str());
	data_hist_proj_x->GetXaxis()->SetTitle(x_title.c_str());

	data_hist_proj_y->SetTitle((title + " Y projections data vs MC").c_str());
	data_hist_proj_y->GetXaxis()->SetTitle(y_title.c_str());

	data_hist_proj_x->Rebin(rebin_x);
	data_hist_proj_y->Rebin(rebin_y);

	sim_hist_proj_x->Rebin(rebin_x);
	sim_hist_proj_y->Rebin(rebin_y);
	
	data_hist_proj_x->SetMaximum(data_hist_proj_x->GetMaximum()*1.3);
	data_hist_proj_y->SetMaximum(data_hist_proj_y->GetMaximum()*1.3);

	data_hist_proj_x->GetXaxis()->SetLabelSize(0.05);
	data_hist_proj_x->GetYaxis()->SetLabelSize(0.05);

	data_hist_proj_y->GetXaxis()->SetLabelSize(0.05);
	data_hist_proj_y->GetYaxis()->SetLabelSize(0.05);
	
	TCanvas canv = TCanvas("canv", "canv", 900, 450);
	canv.Divide(2);

	canv.cd(1);
	gPad->SetLeftMargin(0.125);
	gPad->SetRightMargin(0.12);
	data_raw_hist->Clone()->Draw("colz");
	
	canv.cd(2);
	gPad->SetLeftMargin(0.125);
	gPad->SetRightMargin(0.12);
	data_hist->Clone()->Draw("colz");

	PrintCanvas(&canv, "../output/deadmaps/" + Par.run_name + "/" + name);

	TCanvas proj_canv = TCanvas("proj_canv", "canv", 1800, 900);
	
	proj_canv.cd();
	gPad->Divide(2, 2);
	
	proj_canv.cd(1);
	gPad->SetLeftMargin(0.12);
	gPad->SetRightMargin(0.12);
	data_hist->Draw("colz");
	
	proj_canv.cd(3);
	gPad->SetLeftMargin(0.12);
	gPad->SetRightMargin(0.12);
	sim_hist->Draw("colz");
	
	proj_canv.cd(2);
	gPad->SetLeftMargin(0.12);
	gPad->SetBottomMargin(0.11);
	data_hist_proj_x->Draw();
	sim_hist_proj_x->Draw("SAME HIST");
	proj_x_legend.Draw();

	proj_canv.cd(4);
	gPad->SetLeftMargin(0.12);
	gPad->SetBottomMargin(0.11);
	data_hist_proj_y->Draw();
	sim_hist_proj_y->Draw("SAME HIST");
	proj_y_legend.Draw();

	PrintCanvas(&proj_canv, "../output/systematics/" + Par.run_name + "/" + name);

	const double uncertainty_x = GetSystUncertainty(data_hist_proj_x, sim_hist_proj_x);
	const double uncertainty_y = GetSystUncertainty(data_hist_proj_y, sim_hist_proj_y);
	const double uncertainty = RMS(uncertainty_x, uncertainty_y);

	const double data_lost = (1. - data_hist->Integral()/data_raw_hist->Integral())*100.;
	const double sim_lost = (1. - sim_hist->Integral()/sim_raw_hist->Integral())*100.;

	Par.table.PrintRow(name, 
		DtoStr(uncertainty*100., 3), 
		DtoStr(uncertainty_x*100., 3), 
		DtoStr(uncertainty_y*100., 3), 
		DtoStr(data_lost, 3) + "%",
		DtoStr(sim_lost, 3) + "%");

	return uncertainty;
}

void DM()
{
	gROOT->SetBatch(kTRUE);
	gStyle->SetOptStat(0);
	
	gErrorIgnoreLevel = kWarning;

	std::string data_input_file_name = "../data/" + Par.run_name + "/" + Par.data_file_name;
	std::string sim_input_file_name = "../data/phenix_sim/" + Par.run_name + "/" + Par.sim_file_name;
	
	CheckInputFile(data_input_file_name);
	CheckInputFile(sim_input_file_name);
	
	TFile data_file = TFile(data_input_file_name.c_str());
	TFile sim_file = TFile(sim_input_file_name.c_str());

	system(("mkdir -p ../data/deadmaps/" + Par.run_name).c_str());
	system(("mkdir -p ../data/systematics/" + Par.run_name).c_str());
	system(("mkdir -p ../output/deadmaps/" + Par.run_name).c_str());
	system(("mkdir -p ../output/systematics/" + Par.run_name).c_str());

	Par.table.Begin("Acceptance uncertainty");
	Par.table.PrintHeader("detector", "sys", "sys x", "sys y", "data lost", "MC lost");

	//Reading hists from files
	
	//DC
	TH2F *dceast0_data = (TH2F *) data_file.Get("dceast0");
	TH2F *dceast1_data = (TH2F *) data_file.Get("dceast1");
	TH2F *dcwest0_data = (TH2F *) data_file.Get("dcwest0");
	TH2F *dcwest1_data = (TH2F *) data_file.Get("dcwest1");
	
	TH2F *dceast0_sim = (TH2F *) sim_file.Get("dceast0");
	TH2F *dceast1_sim = (TH2F *) sim_file.Get("dceast1");
	TH2F *dcwest0_sim = (TH2F *) sim_file.Get("dcwest0");
	TH2F *dcwest1_sim = (TH2F *) sim_file.Get("dcwest1");

	//PC1
	TH2F *pc1e_data = (TH2F *) data_file.Get("pc1e_z_vs_phi");
	TH2F *pc1w_data = (TH2F *) data_file.Get("pc1w_z_vs_phi");

	TH2F *pc1e_sim = (TH2F *) sim_file.Get("pc1e_z_vs_phi");
	TH2F *pc1w_sim = (TH2F *) sim_file.Get("pc1w_z_vs_phi");
	
	//PC2
	TH2F *pc2_data = (TH2F *) data_file.Get("pc2_z_vs_phi");
	TH2F *pc2_sim = (TH2F *) sim_file.Get("pc2_z_vs_phi");

	//PC3
	TH2F *pc3e_data = (TH2F *) data_file.Get("pc3e_z_vs_phi");
	TH2F *pc3w_data = (TH2F *) data_file.Get("pc3w_z_vs_phi");
	
	TH2F *pc3e_sim = (TH2F *) sim_file.Get("pc3e_z_vs_phi");
	TH2F *pc3w_sim = (TH2F *) sim_file.Get("pc3w_z_vs_phi");

	//TOF
	TH2F *tofe0_data = (TH2F *) data_file.Get("tofe0");
	TH2F *tofe1_data = (TH2F *) data_file.Get("tofe1");
	
	TH2F *tofe0_sim = (TH2F *) sim_file.Get("tofe0");
	TH2F *tofe1_sim = (TH2F *) sim_file.Get("tofe1");

	TH2F *tofw0_data = (TH2F *) data_file.Get("tofw0");
	TH2F *tofw1_data = (TH2F *) data_file.Get("tofw1");
	
	TH2F *tofw0_sim = (TH2F *) sim_file.Get("tofw0");
	TH2F *tofw1_sim = (TH2F *) sim_file.Get("tofw1");

	//EMCal
	TH2F *emcale0_pos_data = (TH2F *) data_file.Get("emcale0_pos");
	TH2F *emcale1_pos_data = (TH2F *) data_file.Get("emcale1_pos");
	TH2F *emcale2_pos_data = (TH2F *) data_file.Get("emcale2_pos");
	TH2F *emcale3_pos_data = (TH2F *) data_file.Get("emcale3_pos");
	TH2F *emcale0_neg_data = (TH2F *) data_file.Get("emcale0_neg");
	TH2F *emcale1_neg_data = (TH2F *) data_file.Get("emcale1_neg");
	TH2F *emcale2_neg_data = (TH2F *) data_file.Get("emcale2_neg");
	TH2F *emcale3_neg_data = (TH2F *) data_file.Get("emcale3_neg");
	TH2F *emcalw0_pos_data = (TH2F *) data_file.Get("emcalw0_pos");
	TH2F *emcalw1_pos_data = (TH2F *) data_file.Get("emcalw1_pos");
	TH2F *emcalw2_pos_data = (TH2F *) data_file.Get("emcalw2_pos");
	TH2F *emcalw3_pos_data = (TH2F *) data_file.Get("emcalw3_pos");
	TH2F *emcalw0_neg_data = (TH2F *) data_file.Get("emcalw0_neg");
	TH2F *emcalw1_neg_data = (TH2F *) data_file.Get("emcalw1_neg");
	TH2F *emcalw2_neg_data = (TH2F *) data_file.Get("emcalw2_neg");
	TH2F *emcalw3_neg_data = (TH2F *) data_file.Get("emcalw3_neg");
	
	TH2F *emcale0_pos_sim = (TH2F *) sim_file.Get("emcale0_pos");
	TH2F *emcale1_pos_sim = (TH2F *) sim_file.Get("emcale1_pos");
	TH2F *emcale2_pos_sim = (TH2F *) sim_file.Get("emcale2_pos");
	TH2F *emcale3_pos_sim = (TH2F *) sim_file.Get("emcale3_pos");
	TH2F *emcale0_neg_sim = (TH2F *) sim_file.Get("emcale0_neg");
	TH2F *emcale1_neg_sim = (TH2F *) sim_file.Get("emcale1_neg");
	TH2F *emcale2_neg_sim = (TH2F *) sim_file.Get("emcale2_neg");
	TH2F *emcale3_neg_sim = (TH2F *) sim_file.Get("emcale3_neg");
	TH2F *emcalw0_pos_sim = (TH2F *) sim_file.Get("emcalw0_pos");
	TH2F *emcalw1_pos_sim = (TH2F *) sim_file.Get("emcalw1_pos");
	TH2F *emcalw2_pos_sim = (TH2F *) sim_file.Get("emcalw2_pos");
	TH2F *emcalw3_pos_sim = (TH2F *) sim_file.Get("emcalw3_pos");
	TH2F *emcalw0_neg_sim = (TH2F *) sim_file.Get("emcalw0_neg");
	TH2F *emcalw1_neg_sim = (TH2F *) sim_file.Get("emcalw1_neg");
	TH2F *emcalw2_neg_sim = (TH2F *) sim_file.Get("emcalw2_neg");
	TH2F *emcalw3_neg_sim = (TH2F *) sim_file.Get("emcalw3_neg");

	//Dead areas cuts
	
	//DC
	TH2F *dceast0_data_cut = CutDeadAreas((TH2F *) dceast0_data->Clone(), &IsDeadDC, 2., 1.);
	TH2F *dceast1_data_cut = CutDeadAreas((TH2F *) dceast1_data->Clone(), &IsDeadDC, 2., -1.);
	TH2F *dcwest0_data_cut = CutDeadAreas((TH2F *) dcwest0_data->Clone(), &IsDeadDC, 1., 1.);
	TH2F *dcwest1_data_cut = CutDeadAreas((TH2F *) dcwest1_data->Clone(), &IsDeadDC, 1., -1.);

	TH2F *dceast0_sim_cut = CutDeadAreas((TH2F *) dceast0_sim->Clone(), &IsDeadDC, 2., 1.);
	TH2F *dceast1_sim_cut = CutDeadAreas((TH2F *) dceast1_sim->Clone(), &IsDeadDC, 2., -1.);
	TH2F *dcwest0_sim_cut = CutDeadAreas((TH2F *) dcwest0_sim->Clone(), &IsDeadDC, 1., 1.);
	TH2F *dcwest1_sim_cut = CutDeadAreas((TH2F *) dcwest1_sim->Clone(), &IsDeadDC, 1., -1.);

	//PC1
	TH2F *pc1e_data_cut = CutDeadAreas((TH2F *) pc1e_data->Clone(), &IsDeadPC1, 2.);
	TH2F *pc1w_data_cut = CutDeadAreas((TH2F *) pc1w_data->Clone(), &IsDeadPC1, 1.);

	TH2F *pc1e_sim_cut = CutDeadAreas((TH2F *) pc1e_sim->Clone(), &IsDeadPC1, 2.);
	TH2F *pc1w_sim_cut = CutDeadAreas((TH2F *) pc1w_sim->Clone(), &IsDeadPC1, 1.);
	
	//PC2
	TH2F *pc2_data_cut = CutDeadAreas((TH2F *) pc2_data->Clone(), &IsDeadPC2);
	TH2F *pc2_sim_cut = CutDeadAreas((TH2F *) pc2_sim->Clone(), &IsDeadPC2);

	//PC3
	TH2F *pc3e_data_cut = CutDeadAreas((TH2F *) pc3e_data->Clone(), &IsDeadPC3, 2.);
	TH2F *pc3w_data_cut = CutDeadAreas((TH2F *) pc3w_data->Clone(), &IsDeadPC3, 1.);

	TH2F *pc3e_sim_cut = CutDeadAreas((TH2F *) pc3e_sim->Clone(), &IsDeadPC3, 2.);
	TH2F *pc3w_sim_cut = CutDeadAreas((TH2F *) pc3w_sim->Clone(), &IsDeadPC3, 1.);

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

	TH2F *emcale0_pos_sim_cut = CutDeadAreas((TH2F *) emcale0_pos_sim->Clone(), &IsDeadEMCal, 2., 1., 0);
	TH2F *emcale1_pos_sim_cut = CutDeadAreas((TH2F *) emcale1_pos_sim->Clone(), &IsDeadEMCal, 2., 1., 1);
	TH2F *emcale2_pos_sim_cut = CutDeadAreas((TH2F *) emcale2_pos_sim->Clone(), &IsDeadEMCal, 2., 1., 2);
	TH2F *emcale3_pos_sim_cut = CutDeadAreas((TH2F *) emcale3_pos_sim->Clone(), &IsDeadEMCal, 2., 1., 3);
	TH2F *emcale0_neg_sim_cut = CutDeadAreas((TH2F *) emcale0_neg_sim->Clone(), &IsDeadEMCal, 2., -1., 0);
	TH2F *emcale1_neg_sim_cut = CutDeadAreas((TH2F *) emcale1_neg_sim->Clone(), &IsDeadEMCal, 2., -1., 1);
	TH2F *emcale2_neg_sim_cut = CutDeadAreas((TH2F *) emcale2_neg_sim->Clone(), &IsDeadEMCal, 2., -1., 2);
	TH2F *emcale3_neg_sim_cut = CutDeadAreas((TH2F *) emcale3_neg_sim->Clone(), &IsDeadEMCal, 2., -1., 3);
	TH2F *emcalw0_pos_sim_cut = CutDeadAreas((TH2F *) emcalw0_pos_sim->Clone(), &IsDeadEMCal, 1., 1., 0);
	TH2F *emcalw1_pos_sim_cut = CutDeadAreas((TH2F *) emcalw1_pos_sim->Clone(), &IsDeadEMCal, 1., 1., 1);
	TH2F *emcalw2_pos_sim_cut = CutDeadAreas((TH2F *) emcalw2_pos_sim->Clone(), &IsDeadEMCal, 1., 1., 2);
	TH2F *emcalw3_pos_sim_cut = CutDeadAreas((TH2F *) emcalw3_pos_sim->Clone(), &IsDeadEMCal, 1., 1., 3);
	TH2F *emcalw0_neg_sim_cut = CutDeadAreas((TH2F *) emcalw0_neg_sim->Clone(), &IsDeadEMCal, 1., -1., 0);
	TH2F *emcalw1_neg_sim_cut = CutDeadAreas((TH2F *) emcalw1_neg_sim->Clone(), &IsDeadEMCal, 1., -1., 1);
	TH2F *emcalw2_neg_sim_cut = CutDeadAreas((TH2F *) emcalw2_neg_sim->Clone(), &IsDeadEMCal, 1., -1., 2);
	TH2F *emcalw3_neg_sim_cut = CutDeadAreas((TH2F *) emcalw3_neg_sim->Clone(), &IsDeadEMCal, 1., -1., 3);
	
	//TOFe
	TH2F *tofe0_data_cut = CutDeadAreas((TH2F *) tofe0_data->Clone(), &IsDeadTOFe, 1.);
	TH2F *tofe1_data_cut = CutDeadAreas((TH2F *) tofe1_data->Clone(), &IsDeadTOFe, -1.);

	TH2F *tofe0_sim_cut = CutDeadAreas((TH2F *) tofe0_sim->Clone(), &IsDeadTOFe, 1.);
	TH2F *tofe1_sim_cut = CutDeadAreas((TH2F *) tofe1_sim->Clone(), &IsDeadTOFe, -1.);
	
	//TOFw
	TH2F *tofw0_data_cut = CutDeadAreas((TH2F *) tofw0_data->Clone(), &IsDeadTOFw, 1.);
	TH2F *tofw1_data_cut = CutDeadAreas((TH2F *) tofw1_data->Clone(), &IsDeadTOFw, -1.);

	TH2F *tofw0_sim_cut = CutDeadAreas((TH2F *) tofw0_sim->Clone(), &IsDeadTOFw, 1.);
	TH2F *tofw1_sim_cut = CutDeadAreas((TH2F *) tofw1_sim->Clone(), &IsDeadTOFw, -1.);
	
	//DC
	const double dceast0_sys = GetUncertainty1Proj(
		dceast0_data, dceast0_sim, dceast0_data_cut, dceast0_sim_cut,
		"dceast0", "DC east (z #geq 0)", "board_{DC}", "#alpha_{DC}");
	const double dceast1_sys = GetUncertainty1Proj(
		dceast1_data, dceast1_sim, dceast1_data_cut, dceast1_sim_cut, 
		"dceast1", "DC east (z < 0)", "board_{DC}", "#alpha_{DC}");
	const double dcwest0_sys = GetUncertainty1Proj(
		dcwest0_data, dcwest0_sim, dcwest0_data_cut, dcwest0_sim_cut, 
		"dcwest0", "DC west (z #geq 0)", "board_{DC}", "#alpha_{DC}");
	const double dcwest1_sys = GetUncertainty1Proj(
		dcwest1_data, dcwest1_sim, dcwest1_data_cut, dcwest1_sim_cut, 
		"dcwest1", "DC west (z < 0)", "board_{DC}", "#alpha_{DC}");

	//PC1
	const double pc1e_sys = GetUncertainty2Proj(
		pc1e_data, pc1e_sim, pc1e_data_cut, pc1e_sim_cut, 
		"pc1e", "PC1e", "Z_{PC1}", "#varphi_{PC1}");
	const double pc1w_sys = GetUncertainty2Proj(
		pc1w_data, pc1w_sim, pc1w_data_cut, pc1w_sim_cut, 
		"pc1w", "PC1w", "Z_{PC1}", "#varphi_{PC1}");
	
	//PC2
	const double pc2_sys = GetUncertainty2Proj(
		pc2_data, pc2_sim, pc2_data_cut, pc2_sim_cut, 
		"pc2", "PC2", "Z_{PC2}", "#varphi_{PC2}");

	//PC3
	const double pc3e_sys = GetUncertainty2Proj(
		pc3e_data, pc3e_sim, pc3e_data_cut, pc3e_sim_cut, 
		"pc3e", "PC3e", "Z_{PC3}", "#varphi_{PC3}");
	const double pc3w_sys = GetUncertainty2Proj(
		pc3w_data, pc3w_sim, pc3w_data_cut, pc3w_sim_cut, 
		"pc3w", "PC3w", "Z_{PC3}", "#varphi_{PC3}");

	//TOF
	const double tofe0_sys = GetUncertainty2Proj(
		tofe0_data, tofe0_sim, tofe0_data_cut, tofe0_sim_cut, 
		"tofe0", "TOFe, z #geq 0", "Y_{TOF}", "Z_{TOF}");
	const double tofe1_sys = GetUncertainty2Proj(
		tofe1_data, tofe1_sim, tofe1_data_cut, tofe1_sim_cut, 
		"tofe1", "TOFe, z < 0", "Y_{TOF}", "Z_{TOF}");
	const double tofw0_sys = GetUncertainty1Proj(
		tofw0_data, tofw0_sim, tofw0_data_cut, tofw0_sim_cut, 
		"tofw0", "TOFw, z #geq 0", "board_{DC}", "#alpha_{DC}");
	const double tofw1_sys = GetUncertainty1Proj(
		tofw1_data, tofw1_sim, tofw1_data_cut, tofw1_sim_cut, 
		"tofw1", "TOFw, z < 0", "board_{DC}", "#alpha_{DC}");

	//EMCal
	const double emcale0_pos_sys = GetUncertainty2Proj(
		emcale0_pos_data, emcale0_pos_sim, emcale0_pos_data_cut, emcale0_pos_sim_cut, 
		"emcale0_pos", "EMCale0, z #geq 0", "Y_{EMC}", "Z_{EMC}", 1, 2);
	const double emcale1_pos_sys = GetUncertainty2Proj(
		emcale1_pos_data, emcale1_pos_sim, emcale1_pos_data_cut, emcale1_pos_sim_cut, 
		"emcale1_pos", "EMCale1, z #geq 0", "Y_{EMC}", "Z_{EMC}", 1, 2);
	const double emcale2_pos_sys = GetUncertainty2Proj(
		emcale2_pos_data, emcale2_pos_sim, emcale2_pos_data_cut, emcale2_pos_sim_cut, 
		"emcale2_pos", "EMCale2, z #geq 0", "Y_{EMC}", "Z_{EMC}", 1, 2);
	const double emcale3_pos_sys = GetUncertainty2Proj(
		emcale3_pos_data, emcale3_pos_sim, emcale3_pos_data_cut, emcale3_pos_sim_cut, 
		"emcale3_pos", "EMCale3, z #geq 0", "Y_{EMC}", "Z_{EMC}", 1, 2);

	const double emcale0_neg_sys = GetUncertainty2Proj(
		emcale0_neg_data, emcale0_neg_sim, emcale0_neg_data_cut, emcale0_neg_sim_cut, 
		"emcale0_neg", "EMCale0, z < 0", "Y_{EMC}", "Z_{EMC}", 1, 2);
	const double emcale1_neg_sys = GetUncertainty2Proj(
		emcale1_neg_data, emcale1_neg_sim, emcale1_neg_data_cut, emcale1_neg_sim_cut, 
		"emcale1_neg", "EMCale1, z < 0", "Y_{EMC}", "Z_{EMC}", 1, 2);
	const double emcale2_neg_sys = GetUncertainty2Proj(
		emcale2_neg_data, emcale2_neg_sim, emcale2_neg_data_cut, emcale2_neg_sim_cut, 
		"emcale2_neg", "EMCale2, z < 0", "Y_{EMC}", "Z_{EMC}", 1, 2);
	const double emcale3_neg_sys = GetUncertainty2Proj(
		emcale3_neg_data, emcale3_neg_sim, emcale3_neg_data_cut, emcale3_neg_sim_cut, 
		"emcale3_neg", "EMCale3, z < 0", "Y_{EMC}", "Z_{EMC}", 1, 2);

	const double emcalw0_pos_sys = GetUncertainty2Proj(
		emcalw0_pos_data, emcalw0_pos_sim, emcalw0_pos_data_cut, emcalw0_pos_sim_cut, 
		"emcalw0_pos", "EMCalw0, z #geq 0", "Y_{EMC}", "Z_{EMC}", 1, 2);
	const double emcalw1_pos_sys = GetUncertainty2Proj(
		emcalw1_pos_data, emcalw1_pos_sim, emcalw1_pos_data_cut, emcalw1_pos_sim_cut, 
		"emcalw1_pos", "EMCalw1, z #geq 0", "Y_{EMC}", "Z_{EMC}", 1, 2);
	const double emcalw2_pos_sys = GetUncertainty2Proj(
		emcalw2_pos_data, emcalw2_pos_sim, emcalw2_pos_data_cut, emcalw2_pos_sim_cut, 
		"emcalw2_pos", "EMCalw2, z #geq 0", "Y_{EMC}", "Z_{EMC}", 1, 2);
	const double emcalw3_pos_sys = GetUncertainty2Proj(
		emcalw3_pos_data, emcalw3_pos_sim, emcalw3_pos_data_cut, emcalw3_pos_sim_cut, 
		"emcalw3_pos", "EMCalw3, z #geq 0", "Y_{EMC}", "Z_{EMC}", 1, 2);

	const double emcalw0_neg_sys = GetUncertainty2Proj(
		emcalw0_neg_data, emcalw0_neg_sim, emcalw0_neg_data_cut, emcalw0_neg_sim_cut, 
		"emcalw0_neg", "EMCalw0, z < 0", "Y_{EMC}", "Z_{EMC}", 1, 2);
	const double emcalw1_neg_sys = GetUncertainty2Proj(
		emcalw1_neg_data, emcalw1_neg_sim, emcalw1_neg_data_cut, emcalw1_neg_sim_cut, 
		"emcalw1_neg", "EMCalw1, z < 0", "Y_{EMC}", "Z_{EMC}", 1, 2);
	const double emcalw2_neg_sys = GetUncertainty2Proj(
		emcalw2_neg_data, emcalw2_neg_sim, emcalw2_neg_data_cut, emcalw2_neg_sim_cut, 
		"emcalw2_neg", "EMCalw2, z < 0", "Y_{EMC}", "Z_{EMC}", 1, 2);
	const double emcalw3_neg_sys = GetUncertainty2Proj(
		emcalw3_neg_data, emcalw3_neg_sim, emcalw3_neg_data_cut, emcalw3_neg_sim_cut, 
		"emcalw3_neg", "EMCalw3, z < 0", "Y_{EMC}", "Z_{EMC}", 1, 2);

	Par.table.End();

	system(("mkdir -p ../../sim_analysis/input/Systematics/" + Par.run_name).c_str());
			
	WriteFile("../../sim_analysis/input/Systematics/" + Par.run_name + "/acceptance.txt", 
		dceast0_sys, dceast1_sys, dcwest0_sys, dcwest1_sys,
		pc1e_sys, pc1w_sys, pc2_sys, pc3e_sys, pc3w_sys,
		tofe0_sys, tofe1_sys, tofw0_sys, tofw1_sys,
		emcale0_pos_sys, emcale0_neg_sys, emcale1_pos_sys, emcale1_neg_sys, 
		emcale2_pos_sys, emcale2_neg_sys, emcale3_pos_sys, emcale3_neg_sys,
		emcalw0_pos_sys, emcalw0_neg_sys, emcalw1_pos_sys, emcalw1_neg_sys,
		emcalw2_pos_sys, emcalw2_neg_sys, emcalw3_pos_sys, emcalw3_neg_sys);
}
