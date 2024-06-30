#include <vector>
#include <cmath>
#include <string>

#include "TH2.h"
#include "TH1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TLatex.h"

#include "../lib/ErrorHandler.h"
#include "../lib/Particles.h"
#include "../lib/StrTool.h"
#include "../lib/OutputTool.h"
#include "../lib/InputTool.h"
#include "../lib/Tools.h"
#include "../lib/CentralityTypes.h"
#include "../lib/TCanvasPrinter.h"

#include "../lib/BW.h"

struct pmm_noPID
{
	std::string name = "noPID";
	std::string title = "noPID";
	std::vector<double> ptmin = 
		{2.3, 2.6, 2.9, 3.4, 4.0, 4.5, 5.0};
	std::vector<double> ptmax = 
		{2.6, 2.9, 3.4, 4.0, 4.5, 5.0, 6.5};
	bool contains_emcal_id = false;
};

struct pmm_1PID
{
	std::string name = "1PID";
	std::string title = "1PID";
	std::vector<double> ptmin = 
		{1.7, 1.9, 2.1, 2.3, 2.6, 2.9, 3.4, 4.0};
	std::vector<double> ptmax = 
		{1.9, 2.1, 2.3, 2.6, 2.9, 3.4, 4.0, 4.5};
	bool contains_emcal_id = false;
};

struct pmm_2PID
{
	std::string name = "2PID";
	std::string title = "2PID";
	std::vector<double> ptmin = 
		{0.9, 1.1, 1.4, 1.7, 1.9, 2.1, 2.3, 2.6, 2.9, 3.4};
	std::vector<double> ptmax = 
		{1.1, 1.4, 1.7, 1.9, 2.1, 2.3, 2.6, 2.9, 3.4, 4.0};
	bool contains_emcal_id = true;
};

struct pmm_EMC2PID
{
	std::string name = "EMC2PID";
	std::string title = "EMC2PID";
	std::vector<double> ptmin = 
		{0.9, 1.1, 1.4, 1.7, 1.9, 2.1};
	std::vector<double> ptmax = 
		{1.1, 1.4, 1.7, 1.9, 2.1, 2.3};
	bool contains_emcal_id = true;
};

struct pmm_TOF2PID
{
	std::string name = "TOF2PID";
	std::string title = "TOF2PID";
	std::vector<double> ptmin = 
		{1.9, 2.1, 2.3, 2.6, 2.9};
	std::vector<double> ptmax = 
		{2.1, 2.3, 2.6, 2.9, 3.4};
	bool contains_emcal_id = false;
};

struct
{
	//pmm_noPID method;
	//pmm_1PID method;
	//pmm_2PID method;
	pmm_EMC2PID method;
	//pmm_TOF2PID method;
	
	std::string run_name = "Run7AuAu200";

	KStar892 particle;

	AuAu200CTypeMB4 CType;	

	bool print_result = 0;

	const double pt_decreased = 0.995;
	const double pt_increased = 1.005;

	TF1 gdf_sigma = TF1("gdf_sigma", "pol1");

	double min_y = 1e31;
	double max_y = -1e31;

	const double srange = 1;
	const double yield_extr_srange = 1;

	std::vector<TGraphErrors *> VGr;
	
	TLegend legend = TLegend(0.55, 0.1, 0.9, 0.4);

	const int rebin = 1;
	const double min_minv_range = 0.632;
} Par;

double GetYield(TH1F *hist, const double pt, double &err, std::string name)
{
	TCanvas canv = TCanvas(name.c_str(), "canv", 900, 900);
	
	TF1 fit_func = TF1("fit_RBW_GDF_BG", &RBW_GDF_POL2, 0., 10., 7);
	TF1 bg_func = TF1("bg", "pol2");
	TF1 bw_func = TF1("bw_func", &RBW_GDF_conv, 0., 10., 7);
	
	fit_func.SetParameters(
		hist->GetBinContent(hist->FindBin(Par.particle.mass)), 
		Par.particle.mass, Par.particle.gamma);
	
	fit_func.SetParLimits(0, hist->GetBinContent(hist->FindBin(Par.particle.mass))/1.05, 
		hist->GetBinContent(hist->FindBin(Par.particle.mass))*10.);
	fit_func.SetParLimits(1, Par.particle.mass - Par.particle.gamma/5., 
		Par.particle.mass + Par.particle.gamma/5.);
	fit_func.SetParLimits(2, Par.particle.gamma/1.02, Par.particle.gamma*1.02);
	
	fit_func.FixParameter(3, Par.gdf_sigma.Eval(pt));
	
	//background is negligible
	fit_func.FixParameter(4., 0.);
	fit_func.FixParameter(5., 0.);
	fit_func.FixParameter(6., 0.);

	if (pt < 2.)
	{
		fit_func.SetRange(Par.particle.mass - Par.particle.gamma*Par.srange*5., 
			Par.particle.mass + Par.particle.gamma*Par.srange*1.);
	}
	else
	{
		fit_func.SetRange(Par.particle.mass - Par.particle.gamma*Par.srange*3., 
			Par.particle.mass + Par.particle.gamma*Par.srange*3.);
	}
	
	bw_func.SetRange(Par.particle.mass - Par.particle.gamma*Par.srange*5., 
		Par.particle.mass + Par.particle.gamma*Par.srange*5.);

	bw_func.SetLineStyle(2);
	bw_func.SetLineColor(kAzure);

	hist->SetMarkerStyle(8);
	hist->SetMarkerSize(1.);
	hist->SetLineWidth(2);
	hist->SetLineColor(kBlack);
	hist->SetMarkerColor(kBlack);

	hist->Fit(&fit_func, "RQMBN");

	for (int i = 0; i < 4; i++)
	{
		bw_func.FixParameter(i, fit_func.GetParameter(i));
	}

	fit_func.SetRange(Par.particle.mass - Par.particle.gamma*Par.srange*5., 
		Par.particle.mass + Par.particle.gamma*Par.srange*5.);

	double yield = 0;

	double mean = bw_func.GetParameter(1);
	double sigma = bw_func.GetParameter(2);
	
	for (int i = hist->FindBin(mean - sigma*Par.yield_extr_srange);
		i <= hist->FindBin(mean + sigma*Par.yield_extr_srange); i++)
	{
		yield += hist->GetBinContent(i);
	}

	err = sqrt(hist->Integral(
		hist->FindBin(Par.particle.mass - Par.particle.gamma*Par.srange), 
		hist->FindBin(Par.particle.mass + Par.particle.gamma*Par.srange)))/yield;
		
	hist->GetXaxis()->SetRange(
		hist->GetXaxis()->FindBin(Par.min_minv_range),
		hist->GetXaxis()->FindBin(Par.particle.mass + Par.particle.gamma*10.));

	TLine lleft = TLine(mean - sigma*Par.srange, 
		0., mean - sigma*Par.srange, 
		hist->GetBinContent(hist->GetMaximumBin()));

	TLine lright = TLine(mean + sigma*Par.srange, 
		0., mean + sigma*Par.srange, 
		hist->GetBinContent(hist->GetMaximumBin()));
	
	lleft.SetLineColor(kGray);
	lright.SetLineColor(kGray);

	lleft.SetLineStyle(2);
	lright.SetLineStyle(2);

	lleft.SetLineWidth(2);
	lright.SetLineWidth(2);
	
	hist->Draw();

	lleft.Draw();
	lright.Draw();
	
	bg_func.Draw("SAME");
	fit_func.Draw("SAME");
	bw_func.Draw("SAME");	
	
	PrintCanvas(&canv, name, false);
	
	return yield;
}

void AddEffGraph(TFile *input_file, TFile *input_file_pt_decreased, TFile *input_file_pt_increased, const int cnum, Color_t color, Style_t marker_style)
{
	TH2F *hist2d = (TH2F *) input_file->Get((Par.CType.cname_nop[cnum] + "/" + 
		Par.method.name + "_" + Par.CType.cname_nop[cnum]).c_str());
	TH2F *hist2d_pt_decreased = (TH2F *) input_file_pt_decreased->Get((Par.CType.cname_nop[cnum] + "/" + 
		Par.method.name + "_" + Par.CType.cname_nop[cnum]).c_str());
	TH2F *hist2d_pt_increased = (TH2F *) input_file_pt_increased->Get((Par.CType.cname_nop[cnum] + "/" + 
		Par.method.name + "_" + Par.CType.cname_nop[cnum]).c_str());
	TH2F *hist2d_acc_decreased = (TH2F *) input_file->Get((Par.CType.cname_nop[cnum] + "/" +
		Par.method.name + "_acc_decreased_" + Par.CType.cname_nop[cnum]).c_str());
	TH2F *hist2d_acc_increased = (TH2F *) input_file->Get((Par.CType.cname_nop[cnum] + "/" +
		Par.method.name + "_acc_increased_" + Par.CType.cname_nop[cnum]).c_str());

	TH2F *hist2d_m2_eff_decreased, *hist2d_m2_eff_increased;

	if (Par.method.contains_emcal_id)
	{
		hist2d_m2_eff_decreased = (TH2F *) input_file->Get((Par.CType.cname_nop[cnum] + "/" +
			Par.method.name + "_m2_eff_decreased_" + Par.CType.cname_nop[cnum]).c_str());
		hist2d_m2_eff_increased = (TH2F *) input_file->Get((Par.CType.cname_nop[cnum] + "/" +
			Par.method.name + "_m2_eff_increased_" + Par.CType.cname_nop[cnum]).c_str());
	}

	TH1F *orig = (TH1F *) input_file->Get("orig");

	TGraphErrors graph = TGraphErrors();
	
	graph.SetMarkerStyle(marker_style);
	graph.SetMarkerSize(2);
	graph.SetMarkerColor(color-3);
	graph.SetLineColor(color-3);
	graph.SetLineWidth(2);
	graph.SetFillColor(color+9);
	graph.SetFillStyle(3001);
	
	TF1 gausn = TF1("gausn", "gausn");
	gausn.SetParameters(1, 0, 1);

	std::string output_name = "";
	if (Par.print_result) output_name = "../par/Efficiency/" + Par.run_name + "/" + 
		Par.particle.name_nl + "/" + Par.method.name + "_" + Par.CType.cname_nop[cnum] + ".txt";
	
	system(("mkdir -p ../par/Efficiency/" + Par.run_name + "/" + Par.particle.name_nl).c_str());
	ofstream output(output_name.c_str());
	
	system(("mkdir -p ../output/Efficiency/" + Par.run_name + "/" + 
	Par.method.name + "_" + Par.particle.name_nl + 
	"/" + Par.CType.cname_nop[cnum]).c_str());

	Print("pt", "eff", "err", "pair_uncertainty", "pt_scale_uncertainty", 
		"acc_uncertainty", "m2_eff_uncertainty");
	for (int i = 0; i < Par.method.ptmin.size(); i++)
	{
		const double orig_yield = orig->Integral(orig->FindBin(Par.method.ptmin[i]+0.01), 
			orig->FindBin(Par.method.ptmax[i]-0.01));
		
		TH1F *hist1d = (TH1F *) hist2d->ProjectionY(Par.method.name.c_str(), 
			hist2d->GetXaxis()->FindBin(Par.method.ptmin[i]+0.01), 
			hist2d->GetXaxis()->FindBin(Par.method.ptmax[i]-0.01));

		TH1F *hist1d_pt_decreased = (TH1F *) hist2d_pt_decreased->ProjectionY(
			(Par.method.name+"_pt_decreased" + Par.CType.cname_nop[cnum]).c_str(), 
			hist2d_pt_decreased->GetXaxis()->FindBin(Par.method.ptmin[i]+0.01), 
			hist2d_pt_decreased->GetXaxis()->FindBin(Par.method.ptmax[i]-0.01));

		TH1F *hist1d_pt_increased = (TH1F *) hist2d_pt_increased->ProjectionY(
			(Par.method.name+"_pt_increased" + Par.CType.cname_nop[cnum]).c_str(), 
			hist2d_pt_increased->GetXaxis()->FindBin(Par.method.ptmin[i]+0.01), 
			hist2d_pt_increased->GetXaxis()->FindBin(Par.method.ptmax[i]-0.01));

		TH1F *hist1d_acc_decreased = (TH1F *) hist2d_acc_decreased->ProjectionY(
			(Par.method.name+"_acc_decreased" + Par.CType.cname_nop[cnum]).c_str(), 
			hist2d_pt_decreased->GetXaxis()->FindBin(Par.method.ptmin[i]+0.01), 
			hist2d_pt_decreased->GetXaxis()->FindBin(Par.method.ptmax[i]-0.01));

		TH1F *hist1d_acc_increased = (TH1F *) hist2d_acc_increased->ProjectionY(
			(Par.method.name+"_acc_increased" + Par.CType.cname_nop[cnum]).c_str(), 
			hist2d_pt_increased->GetXaxis()->FindBin(Par.method.ptmin[i]+0.01), 
			hist2d_pt_increased->GetXaxis()->FindBin(Par.method.ptmax[i]-0.01));

		TH1F *hist1d_m2_eff_decreased, *hist1d_m2_eff_increased;

		if (Par.method.contains_emcal_id)
		{
			hist1d_m2_eff_decreased = (TH1F *) hist2d_m2_eff_decreased->ProjectionY(
				(Par.method.name+"_m2_eff_decreased").c_str(), 
				hist2d_pt_decreased->GetXaxis()->FindBin(Par.method.ptmin[i]+0.01), 
				hist2d_pt_decreased->GetXaxis()->FindBin(Par.method.ptmax[i]-0.01));

			hist1d_m2_eff_increased = (TH1F *) hist2d_m2_eff_increased->ProjectionY(
				(Par.method.name+"_m2_eff_increased").c_str(), 
				hist2d_pt_increased->GetXaxis()->FindBin(Par.method.ptmin[i]+0.01), 
				hist2d_pt_increased->GetXaxis()->FindBin(Par.method.ptmax[i]-0.01));
		}

		hist1d->Rebin(Par.rebin);
		hist1d_pt_decreased->Rebin(Par.rebin);
		hist1d_pt_increased->Rebin(Par.rebin);

		double eff, err;
		double tmp;

		const double pt = (Par.method.ptmin[i] + Par.method.ptmax[i])/2.;

		const double yield = GetYield(hist1d, pt, err, 
			"../output/Efficiency/" + Par.run_name + "/" + Par.method.name + 
			"_" + Par.particle.name_nl + "/" + Par.CType.cname_nop[cnum] + "/" + 
			DtoStr(Par.method.ptmin[i]) + "-" + DtoStr(Par.method.ptmax[i]));

		const double yield_pt_decreased_ratio = 
			GetNormRatio(GetYield(hist1d_pt_decreased, pt, tmp, 
			"../output/Efficiency/" + Par.run_name + "/" + Par.method.name + 
			"_" + Par.particle.name_nl + "/" + Par.CType.cname_nop[cnum] + "/pT_decreased_" + 
			DtoStr(Par.method.ptmin[i]) + "-" + DtoStr(Par.method.ptmax[i]))/yield);

		const double yield_pt_increased_ratio = 
			GetNormRatio(GetYield(hist1d_pt_increased, pt, tmp, 
			"../output/Efficiency/" + Par.run_name + "/" + Par.method.name + 
			"_" + Par.particle.name_nl + "/" + Par.CType.cname_nop[cnum] + "/pT_increased_" + 
			DtoStr(Par.method.ptmin[i]) + "-" + DtoStr(Par.method.ptmax[i]))/yield);

		const double yield_acc_decreased_ratio = 
			GetNormRatio(GetYield(hist1d_acc_decreased, pt, tmp, 
			"../output/Efficiency/" + Par.run_name + "/" + Par.method.name + 
			"_" + Par.particle.name_nl + "/" + Par.CType.cname_nop[cnum] + "/acc_decreased_" + 
			DtoStr(Par.method.ptmin[i]) + "-" + DtoStr(Par.method.ptmax[i]))/yield);

		const double yield_acc_increased_ratio = 
			GetNormRatio(GetYield(hist1d_acc_increased, pt, tmp, 
			"../output/Efficiency/" + Par.run_name + "/" + Par.method.name + 
			"_" + Par.particle.name_nl + "/" + Par.CType.cname_nop[cnum] + "/acc_increased_" + 
			DtoStr(Par.method.ptmin[i]) + "-" + DtoStr(Par.method.ptmax[i]))/yield);

		double yield_m2_eff_decreased_ratio = 1.;
		double yield_m2_eff_increased_ratio = 1.;

		if (Par.method.contains_emcal_id)
		{
			yield_m2_eff_decreased_ratio = 
				GetNormRatio(GetYield(hist1d_m2_eff_decreased, pt, tmp, 
				"../output/Efficiency/" + Par.run_name + "/" + Par.method.name + 
				"_" + Par.particle.name_nl + "/" + Par.CType.cname_nop[cnum] + "/m2_eff_decreased_" + 
				DtoStr(Par.method.ptmin[i]) + "-" + DtoStr(Par.method.ptmax[i]))/yield);

			yield_m2_eff_increased_ratio = 
				GetNormRatio(GetYield(hist1d_m2_eff_increased, pt, tmp, 
				"../output/Efficiency/" + Par.run_name + "/" + Par.method.name + 
				"_" + Par.particle.name_nl + "/" + Par.CType.cname_nop[cnum] + "/m2_eff_increased_" + 
				DtoStr(Par.method.ptmin[i]) + "-" + DtoStr(Par.method.ptmax[i]))/yield);
		}
		
		const double pt_scale_uncertainty = 
			RMS(yield_pt_decreased_ratio, yield_pt_increased_ratio) - 1.;

		const double acc_uncertainty = 
			RMS(yield_acc_decreased_ratio, yield_acc_increased_ratio) - 1.;

		const double m2_eff_uncertainty = 
			RMS(yield_m2_eff_decreased_ratio, yield_m2_eff_increased_ratio) - 1.;

		const double pair_uncertainty = sqrt(
			pow(pt_scale_uncertainty, 2) + pow(acc_uncertainty, 2) + pow(m2_eff_uncertainty, 2));
			
		eff = yield/orig_yield;
		
		Print(pt, eff, err*100., pair_uncertainty*100., pt_scale_uncertainty*100., 
			acc_uncertainty*100., m2_eff_uncertainty*100.);
		
		if (eff <= 0.) continue;
	
		graph.AddPoint(pt, eff);
		graph.SetPointError(graph.GetN()-1, 0, eff*err);

		Par.min_y = Minimum(Par.min_y, Maximum(eff - eff*err, 1e-6));
		Par.max_y = Maximum(Par.max_y, eff + eff*err);

		output << eff << " " << err << " " << 
			pt_scale_uncertainty << " " << acc_uncertainty << " " << m2_eff_uncertainty << std::endl;
	}

	output.close();
	if (Par.print_result) PrintInfo("File " + output_name + " was written");
	
	Par.legend.AddEntry((TGraphErrors *) graph.Clone(), (Par.CType.cname[cnum]).c_str(), "P");
	
	Par.VGr.push_back((TGraphErrors *) graph.Clone());
}

void WritePtRange()
{
	std::string output_name = "../par/PtRange/" + Par.run_name + "/" + 
		Par.particle.name_nl + "/" + Par.method.name + ".txt";
	
	system(("mkdir -p ../par/PtRange/" + Par.run_name + "/" + Par.particle.name_nl).c_str());
	ofstream output(output_name.c_str());

	for (int i = 0; i < Par.method.ptmin.size(); i++)
	{
		output << Par.method.ptmin[i] << " " << Par.method.ptmax[i] << std::endl;
	}

	PrintInfo("File " + output_name + " was written");
	output.close();
}

void PairEff()
{
	gROOT->SetBatch(kTRUE);
	gStyle->SetOptStat(0);
	
	gErrorIgnoreLevel = kWarning;

	ROOT::EnableImplicitMT();

	system(("mkdir -p ../output/systematics/" + Par.run_name).c_str());
	
	std::string gdf_par_file_name = "../par/GDF/" + Par.run_name + "/" + Par.particle.name_nl + ".txt";
	CheckInputFile(gdf_par_file_name);
	
	CheckInputFile("../data/phenix_sim/" + Par.run_name + "/" + Par.particle.name_nl + ".root");
	CheckInputFile("../data/phenix_sim/" + Par.run_name + "/" + Par.particle.name_nl +
		 "_pt" + DtoStr(Par.pt_decreased, 3) + ".root");
	CheckInputFile("../data/phenix_sim/" + Par.run_name + "/" + Par.particle.name_nl +
		 "_pt" + DtoStr(Par.pt_increased, 3) + ".root");

	WritePtRange();
	
	TFile sim_input_file = TFile(
		("../data/phenix_sim/" + Par.run_name + "/" + Par.particle.name_nl + ".root").c_str());

	TFile sim_input_file_pt_decreased = TFile(
		("../data/phenix_sim/" + Par.run_name + "/" + Par.particle.name_nl +
		 "_pt" + DtoStr(Par.pt_decreased, 3) + ".root").c_str());
	
	TFile sim_input_file_pt_increased = TFile(
		("../data/phenix_sim/" + Par.run_name + "/" + Par.particle.name_nl +
		 "_pt" + DtoStr(Par.pt_increased, 3) + ".root").c_str());

	double *gdf_par = ReadFileIntoArray(gdf_par_file_name, 2);
	Par.gdf_sigma.SetParameters(gdf_par);
	
	for (int i = 0; i < Par.CType.size; i++)
	{
		AddEffGraph(&sim_input_file, &sim_input_file_pt_decreased, &sim_input_file_pt_increased, 
			i, Par.CType.color[i], Par.CType.marker_style[i]);
	}
	
	TH1F hr = TH1F("range_hist", 
		"", 10, Par.method.ptmin.front()/1.1, Par.method.ptmax.back());
	
	hr.GetXaxis()->SetTitle("p_{T}, GeV/c");
	hr.GetYaxis()->SetTitle("#epsilon");
	
	hr.GetXaxis()->SetTitleSize(0.05);
	hr.GetYaxis()->SetTitleSize(0.05);

	hr.GetXaxis()->SetLabelSize(0.05);
	hr.GetYaxis()->SetLabelSize(0.05);
	
	hr.SetMinimum(Par.min_y/1.5);
	hr.SetMaximum(Par.max_y*1.5);

	TCanvas canv = TCanvas("canv", "canv", 900, 600);
	gPad->SetLogy();

	gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.11);

	hr.Draw();
	hr.Draw("SAME AXIS X+ Y+");

	for (TGraphErrors *graph : Par.VGr) graph->Draw("P");
	
	Par.legend.SetNColumns(2);
	Par.legend.SetLineColorAlpha(0, 0);
	Par.legend.SetFillColorAlpha(0, 0);
	Par.legend.Draw();	

	PrintCanvas(&canv, 
		"../output/Efficiency/" + Par.run_name + "/PairEff_" 
		+ Par.method.name + "_" + Par.particle.name_nl);
}
