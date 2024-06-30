#include <fstream>
#include <string>
#include <vector>
#include <chrono>

#include "cmath"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TStyle.h"

#include "PBar.hpp"

#include "../lib/OutputTool.h"
#include "../lib/ErrorHandler.h"
#include "../lib/Tools.h"
#include "../lib/DeadAreasCuts.h"
#include "../lib/TCanvasPrinter.h"

using namespace CutsRun7AuAu200MinBias;

struct
{
	std::string run_type = "Run7AuAu200";
	std::string ref_file_name1 = "se-231211.root";
	std::string ref_file_name2 = "se-239461.root";
	const double run_transition = 236100;
	const double chi2ndf_threshold = 2.;
	const double min_ratio = 0.9;
	const double max_ratio = 1.1;
	const double min_mult = 0.8;
	const double max_mult = 1.2;
	const double mult_min_ratio = 1.;
	const double mult_max_ratio = 1.3;
} Par;

TH1D* GetNormProj(TH2D *hist2d1, TH2D *hist2d2, const double nevents)
{
	TH1D *hist1d = (TH1D*) hist2d1->ProjectionY("proj", -1, -1);
	hist1d->Add((TH1D*) hist2d2->ProjectionY("proj", -1, -1));
	
	hist1d->Scale(100./nevents);

	return hist1d;
}

double GetRatio(TH1D *hist1d, TH1D *ref_hist1d)
{
	TH1D *ratio_hist = (TH1D *) hist1d->Clone();
	ratio_hist->Divide(ref_hist1d);
	
	TF1 const_fit = TF1("const", "pol0");
	
	const_fit.SetParameter(0, 1.);
	ratio_hist->Fit(&const_fit, "QNC");
	return const_fit.GetParameter(0.);
}

double GetChi2NDF(TH1D *hist1d, TH1D *ref_hist1d, const double mean)
{	
	double chi2 = 0;
	double ndf = 1.;
	
	TH1D *ratio_hist = (TH1D *) hist1d->Clone();
	ratio_hist->Divide(ref_hist1d);

	TF1 const_fit = TF1("const", "pol0");
	
	const_fit.SetParameter(0, 1.);
	ratio_hist->Fit(&const_fit, "QNP");
	
	return const_fit.GetChisquare()/const_fit.GetNDF();
}

void ApplyDM(TH2D *map, const double phi, const double zed)
{
	double board, alpha;

	int x_entries = map->GetNbinsX();
	int y_entries = map->GetNbinsY();
	int num = map->GetEntries();
	
	//vectors for storing parameters of each bin
	
	vector <bool> pos;

	for (double i = 0; i < x_entries; i++)
	{
		for (int j = 0; j < y_entries; j++)
		{
			board = map->GetXaxis()->GetBinCenter(i+1);
			alpha = map->GetYaxis()->GetBinCenter(j+1);

			if (IsDeadDC(phi, zed, board, alpha)) map->SetBinContent(i, j, 0);
			if (IsDeadDC(phi, zed, board, alpha)) map->SetBinContent(i, j, 0);
		}
	}
}

void GoodRunsSelector()
{
	gROOT->SetBatch(kTRUE);
	std::string input_dir  = "../data/" + Par.run_type + "/";
	
	CheckInputFile(input_dir + Par.ref_file_name1);
	CheckInputFile(input_dir + Par.ref_file_name2);
	
	std::string ls_com = "cd " + input_dir + " && ls se* > runs.txt";
	system(ls_com.c_str());

	std::string runs_file_name = input_dir + "runs.txt";
	ifstream runs_file(runs_file_name.c_str());
	
	vector<std::string> runs;
		
	TFile ref_file1 = TFile((input_dir + Par.ref_file_name1).c_str());
	TFile ref_file2 = TFile((input_dir + Par.ref_file_name2).c_str());

	TH2D *ref_east0_2d1 = (TH2D*) ref_file1.Get("dceast0");
	TH2D *ref_east1_2d1 = (TH2D*) ref_file1.Get("dceast1");
	TH2D *ref_west0_2d1 = (TH2D*) ref_file1.Get("dcwest0");
	TH2D *ref_west1_2d1 = (TH2D*) ref_file1.Get("dcwest1");

	TH2D *ref_east0_2d2 = (TH2D*) ref_file2.Get("dceast0");
	TH2D *ref_east1_2d2 = (TH2D*) ref_file2.Get("dceast1");
	TH2D *ref_west0_2d2 = (TH2D*) ref_file2.Get("dcwest0");
	TH2D *ref_west1_2d2 = (TH2D*) ref_file2.Get("dcwest1");

	ApplyDM(ref_east0_2d1, 1, 2);
	ApplyDM(ref_east1_2d1, -1, 2);
	ApplyDM(ref_west0_2d1, 1, 1);
	ApplyDM(ref_west1_2d1, -1, 1);

	ApplyDM(ref_east0_2d2, 1, 2);
	ApplyDM(ref_east1_2d2, -1, 2);
	ApplyDM(ref_west0_2d2, 1, 1);
	ApplyDM(ref_west1_2d2, -1, 1);

	std::vector<double> bad_runs; 
	std::vector<std::string> good_runs, good_runs_name;

	TH1D *ref_centr_hist1 = (TH1D*) ref_file1.Get("central_bin");
	const double ref_nevents1 = ref_centr_hist1->Integral(1, ref_centr_hist1->GetXaxis()->GetNbins());
	
	TH1D *ref_centr_hist2 = (TH1D*) ref_file2.Get("central_bin");
	const double ref_nevents2 = ref_centr_hist2->Integral(1, ref_centr_hist2->GetXaxis()->GetNbins());
	
	TH1D *ref_east_1d1 = GetNormProj(ref_east0_2d1, ref_east1_2d1, ref_nevents1);
	TH1D *ref_west_1d1 = GetNormProj(ref_west0_2d1, ref_west1_2d1, ref_nevents1);
	TH1D *ref_east_1d2 = GetNormProj(ref_east0_2d2, ref_east1_2d2, ref_nevents2);
	TH1D *ref_west_1d2 = GetNormProj(ref_west0_2d2, ref_west1_2d2, ref_nevents2);
	
	double min_x = 1e9, max_x = -1;
	
	double average_mult = 0.;
	while (!runs_file.eof())
	{
		std::string run_name;
		runs_file >> run_name;
		if (run_name == "") break;

		runs.push_back(run_name);
		
		TFile se_file = TFile((input_dir + run_name).c_str());
		TH1D *mult_hist = (TH1D*) se_file.Get("multiplicity");
		TH1D *centr_hist = (TH1D*) se_file.Get("central_bin");
		const double nevents = centr_hist->Integral(1, centr_hist->GetXaxis()->GetNbins());
		average_mult += mult_hist->GetBinContent(1)/nevents;

		std::string to_erase = "se-.root";

		for (char c : to_erase)
		{
			run_name.erase(std::remove(run_name.begin(), run_name.end(), c), run_name.end());
		}
		
		double run_num = atoi(run_name.c_str());

		min_x = Minimum(min_x, run_num);
		max_x = Maximum(max_x, run_num);
	}

	average_mult /= static_cast<double>(runs.size());

	TGraph mult_graph = TGraph();
	TGraph mult_ratio_graph = TGraph();

	TGraph chi2ndf_east_graph = TGraph();
	TGraph chi2ndf_west_graph = TGraph();
	
	TGraph ratio_east_graph = TGraph();
	TGraph ratio_west_graph = TGraph();
	
	ProgressBar pbar("FANCY", OutputColor::bold_cyan);

	for (unsigned int i = 0; i < runs.size(); i++)
	{
		std::string run_name = runs[i];

		pbar.Print((double) (i+1.)/(runs.size()-1.));

		std::string se_file_name = input_dir + run_name;
		TFile se_file = TFile(se_file_name.c_str());

		TH1D *mult_hist = (TH1D*) se_file.Get("multiplicity");
		TH1D *centr_hist = (TH1D*) se_file.Get("central_bin");

		const double nevents = centr_hist->Integral(1, centr_hist->GetXaxis()->GetNbins());
		const double mult = mult_hist->GetBinContent(1)/nevents/average_mult;
		const double mult_ratio = mult_hist->GetBinContent(2)/mult_hist->GetBinContent(3);
		
		std::string to_erase = "se-.root";
		for (char c : to_erase)
		{
			run_name.erase(std::remove(run_name.begin(), run_name.end(), c), run_name.end());
		}
		
		double run_num = atoi(run_name.c_str());
		
		mult_graph.AddPoint(run_num, mult);
		mult_ratio_graph.AddPoint(run_num, mult_ratio);
		
		TH2D *east0_2d = (TH2D*) se_file.Get("dceast0");
		TH2D *east1_2d = (TH2D*) se_file.Get("dceast1");
		TH2D *west0_2d = (TH2D*) se_file.Get("dcwest0");
		TH2D *west1_2d = (TH2D*) se_file.Get("dcwest1");
		
		ApplyDM(east0_2d, 1, 2);
		ApplyDM(east1_2d, -1, 2);
		ApplyDM(west0_2d, 1, 1);
		ApplyDM(west1_2d, -1, 1);
		
		TH1D *east_1d = GetNormProj(east0_2d, east1_2d, nevents);
		TH1D *west_1d = GetNormProj(west0_2d, west1_2d, nevents);
		
		double chi2ndf_east, chi2ndf_west;
		double ratio_east, ratio_west;
		
		if (run_num < Par.run_transition)
		{
			ratio_east = GetRatio(east_1d, ref_east_1d1);
			ratio_west = GetRatio(west_1d, ref_west_1d1);
			chi2ndf_east = GetChi2NDF(east_1d, ref_east_1d1, ratio_east);
			chi2ndf_west = GetChi2NDF(west_1d, ref_west_1d1, ratio_west);
		}
		else
		{
			ratio_east = GetRatio(east_1d, ref_east_1d2);
			ratio_west = GetRatio(west_1d, ref_west_1d2);
			chi2ndf_east = GetChi2NDF(east_1d, ref_east_1d2, ratio_east);
			chi2ndf_west = GetChi2NDF(west_1d, ref_west_1d2, ratio_west);
		}
		
		//chi2/ndf, const, and multiplicity check
		if (chi2ndf_east > Par.chi2ndf_threshold ||
			chi2ndf_west > Par.chi2ndf_threshold ||
			ratio_east > Par.max_ratio ||
			ratio_west > Par.max_ratio ||
			ratio_east < Par.min_ratio ||
			ratio_west < Par.min_ratio ||
			mult < Par.min_mult || 
			mult > Par.max_mult ||
			mult_ratio < Par.mult_min_ratio || 
			mult_ratio > Par.mult_max_ratio) bad_runs.push_back(run_num);
		else 
		{	
			good_runs_name.push_back(run_name + ".root");
			good_runs.push_back("se-" + run_name + ".root");
		}
		
		chi2ndf_east_graph.AddPoint(run_num, chi2ndf_east);
		chi2ndf_west_graph.AddPoint(run_num, chi2ndf_west);
		
		ratio_east_graph.AddPoint(run_num, ratio_east);
		ratio_west_graph.AddPoint(run_num, ratio_west);
	}

	std::ofstream good_runs_file(input_dir + "good_se_runs.txt");
	std::ofstream good_runs_name_file(input_dir + "good_runs.txt");

	for (auto run : good_runs) good_runs_file << run << " " << std::endl;
	for (auto run : good_runs_name) good_runs_name_file << run << " " << std::endl;

	PrintInfo("File " + input_dir + "good_se_runs.txt was written");
	PrintInfo("File " + input_dir + "good_runs.txt was written");

	PrintSeparator(" Bad runs found in DC ");

	for (double bad_run : bad_runs) std::cout << "se-" << bad_run << ".root ";
	Print();
	for (double bad_run : bad_runs) std::cout << bad_run << ".root ";
	Print();
	for (double bad_run : bad_runs) std::cout << bad_run << " ";
	Print();

	Print(bad_runs.size(), "runs were found to be problematic");

	const int marker_style = 33;
	Color_t marker_color = kBlue+2;

	mult_graph.SetMarkerStyle(marker_style);
	mult_ratio_graph.SetMarkerStyle(marker_style);
	
	mult_graph.SetMarkerColor(marker_color);
	mult_ratio_graph.SetMarkerColor(marker_color);
	
	chi2ndf_east_graph.SetMarkerStyle(marker_style);
	chi2ndf_west_graph.SetMarkerStyle(marker_style);

	chi2ndf_east_graph.SetMarkerColor(marker_color);
	chi2ndf_west_graph.SetMarkerColor(marker_color);
	
	ratio_east_graph.SetMarkerStyle(marker_style);
	ratio_west_graph.SetMarkerStyle(marker_style);
	
	ratio_east_graph.SetMarkerColor(marker_color);
	ratio_west_graph.SetMarkerColor(marker_color);
	
	gStyle->SetOptStat(0);

	TH1D range_hist = TH1D("range hist", "", max_x-min_x, min_x, max_x);
	range_hist.GetXaxis()->SetTitle("run index");	

	range_hist.SetTitleSize(0.05, "X");
	range_hist.SetTitleSize(0.05, "Y");

	range_hist.GetXaxis()->SetLabelSize(0.05);
	range_hist.GetYaxis()->SetLabelSize(0.05);

	range_hist.SetMinimum(0.);
	
	std::string output_dir = "../output/deadmaps/" + Par.run_type + "/";
	system(("mkdir -p " + output_dir).c_str());

	range_hist.SetMaximum(2.);

	TCanvas mult_canv("multcanv", "canv", 1200, 600);
	
	mult_canv.Divide(2);

	mult_canv.cd(1);
	gPad->SetLeftMargin(0.14);
	gPad->SetBottomMargin(0.13);
	range_hist.GetYaxis()->SetTitle("Mean multiplicity per event");
	range_hist.DrawClone();

	mult_graph.Draw("P");
	
	mult_canv.cd(2);
	gPad->SetLeftMargin(0.14);
	gPad->SetBottomMargin(0.13);
	range_hist.GetYaxis()->SetTitle("Charge ratio (+/-)");
	range_hist.DrawClone();

	mult_ratio_graph.Draw("P");

	PrintCanvas(&mult_canv, output_dir + "multiplicity");
	
	range_hist.SetMaximum(Par.chi2ndf_threshold*3.);
	range_hist.GetYaxis()->SetTitle("#chi^{2}/ndf");
	
	TCanvas chi2ndf_canv("chi2ndfcanv", "canv", 1200, 400);
	chi2ndf_canv.Divide(2);
	
	chi2ndf_canv.cd(1);
	gPad->SetBottomMargin(0.13);
	gPad->SetLeftMargin(0.11);
	
	range_hist.SetTitle("DC east");
	range_hist.DrawClone();
	chi2ndf_east_graph.Draw("P");
	
	chi2ndf_canv.cd(2);
	gPad->SetBottomMargin(0.13);
	gPad->SetLeftMargin(0.11);
	
	range_hist.SetTitle("DC west");
	range_hist.DrawClone();
	chi2ndf_west_graph.Draw("P");
	
	PrintCanvas(&chi2ndf_canv, output_dir + "DCChi2Ndf");

	range_hist.SetMaximum(2.);
	range_hist.GetYaxis()->SetTitle("const");

	TCanvas ratio_canv("ratiocanv", "canv", 1200, 400);
	ratio_canv.Divide(2);

	ratio_canv.cd(1);
	gPad->SetBottomMargin(0.13);
	gPad->SetLeftMargin(0.11);
	
	range_hist.SetTitle("DC east");
	range_hist.DrawClone();
	ratio_east_graph.Draw("P");
	
	ratio_canv.cd(2);
	gPad->SetBottomMargin(0.13);
	gPad->SetLeftMargin(0.11);

	range_hist.SetTitle("DC west");
	range_hist.DrawClone();
	ratio_west_graph.Draw("P");

	PrintCanvas(&ratio_canv, output_dir + "Ratio");
}
