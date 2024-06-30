#include "../lib/ErrorHandler.h"
#include "../lib/OutputTool.h"
#include "../lib/StrTool.h"
#include "../lib/TCanvasPrinter.h"

struct TOF
{
	std::string name = "tof";
	std::string title = "TOF";
	std::array<std::string, 2> side = {"e", "w"};
	std::array<std::string, 1> sect = {""};
};

struct EMCal
{
	std::string name = "emcal";
	std::string title = "EMCal";
	std::array<std::string, 2> side = {"e", "w"};
	std::array<std::string, 4> sect = {"0", "1", "2", "3"};
};

struct
{
	std::string data_dir = "../data";
	std::string run_name = "Run7AuAu200";
	std::array<std::string, 3> magf = {"", "+-", "-+"};

	EMCal Det;
	
	/*
	std::vector<double> ptmin = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
		1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
		2.2, 2.4, 2.6, 2.8, 3.0, 3.5};

	std::vector<double> ptmax = {0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 
		1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 
		2.4, 2.6, 2.8, 3.0, 3.5, 4.};
	*/

	std::vector<double> ptmin = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
		1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9};

	std::vector<double> ptmax = {0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 
		1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0};

	std::vector<TGraphErrors *> VGr;

	TLegend *legend;
	TH1F hr;
	TLine line = TLine(ptmin.front(), 1., ptmax.back(), 1.);

	double min_y;
	double max_y;
} Par;

void CheckMinY(double y)
{
	if (y < Par.min_y) Par.min_y = y;
}

void CheckMaxY(double y)
{
	if (y > Par.max_y) Par.max_y = y;
}

void AddEffGraph(std::string magf, std::string side, std::string sect, std::string particle, std::string legend_entry, Color_t color, Style_t marker_style)
{
	std::string input_file_name = Par.data_dir + "/phenix_sim/" + 
		Par.run_name + "/" + particle + magf + ".root";

	CheckInputFile(input_file_name);
	
	TFile input_file = TFile(input_file_name.c_str());

	
	TGraphErrors graph = TGraphErrors();

	TH1F *orig_hist = (TH1F*) input_file.Get("orig");
	TH1F *reg_hist = (TH1F*) input_file.Get(("reg_" + Par.Det.name + side + sect).c_str());
	
	graph.SetMarkerColor(color);
	graph.SetLineColor(color);
	graph.SetLineWidth(2);
	
	graph.SetMarkerStyle(marker_style);
	graph.SetMarkerSize(2);

	std::string outfile_name = Par.data_dir + "/Efficiency/" + Par.run_name +
		"/Reg/"  + Par.Det.title + side + sect + "/" + particle + magf + ".txt";
	ofstream output_file(outfile_name.c_str());

	for (int i = 0; i < Par.ptmin.size(); i++)
	{
		double reg_int_err = 0;
		double orig_int_err = 0;
		
		double reg_int = 
			reg_hist->IntegralAndError(reg_hist->FindBin(Par.ptmin[i]+0.01), 
			reg_hist->FindBin(Par.ptmax[i]-0.01), reg_int_err);
			
		double orig_int;
		
		orig_int = 
			orig_hist->IntegralAndError(orig_hist->FindBin(Par.ptmin[i]+0.01), 
			orig_hist->FindBin(Par.ptmax[i]-0.01), orig_int_err);

		const double pt = (Par.ptmin[i]+Par.ptmax[i])/2.;
		
		if (reg_int <= 0) PrintError("Integral of registering histogram "\
			"is insufficient at pt = " + DtoStr(pt, 2) + " GeV : " + to_string(reg_int));
		if (orig_int <= 0) PrintError("Integral of original histogram "\
			"is insufficient at pt = " + DtoStr(pt, 2) + " GeV : " + to_string(orig_int));

		const double eff = reg_int/orig_int;
		const double err = sqrt(
			pow(reg_int_err/reg_int, 2) +
			pow(orig_int_err/orig_int, 2))*eff;

		CheckMinY(eff-err);
		CheckMaxY(eff+err);

		output_file << pt << " " << eff << " " << err/eff << std::endl;

		graph.AddPoint(pt, eff);
		graph.SetPointError(graph.GetN()-1, 0, err);
	}

	PrintInfo("File " + outfile_name + " was written");
	output_file.close();

	Par.legend->AddEntry(graph.Clone(), legend_entry.c_str(), "P");
	Par.VGr.push_back((TGraphErrors *) graph.Clone());
}

void PerformSingleRegEff(std::string magf, std::string side, std::string sect = "")
{
	Par.VGr.clear();
	
	TCanvas canv = TCanvas("canv", "canv", 900, 450);
	Par.legend = new TLegend(0.65, 0.15, 0.95, 0.4);
	Par.hr = TH1F((magf + side + sect).c_str(), 
		"", 10, Par.ptmin.front(), Par.ptmax.back());
			
	Par.min_y = 1e31;
	Par.max_y = -1e31;

	gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.16);
	
	system(("mkdir -p " + Par.data_dir + "/Efficiency/" + 
		Par.run_name + "/Reg/" + Par.Det.title + side + sect).c_str());

	AddEffGraph(magf, side, sect, "pion", "#pi^{+}", kRed-3, 55);
	AddEffGraph(magf, side, sect, "kaon", "K^{+}", kOrange-3, 53);
	AddEffGraph(magf, side, sect, "apion", "#pi^{-}", kAzure-3, 59);
	AddEffGraph(magf, side, sect, "akaon", "K^{-}", kGreen-3, 54);

	Par.hr.GetXaxis()->SetTitle("p_{T}, GeV/c");
	Par.hr.GetYaxis()->SetTitle("#epsilon_{reg}");

	Par.hr.SetTitleSize(0.07, "X");
	Par.hr.SetTitleSize(0.07, "Y");

	Par.hr.GetXaxis()->SetLabelSize(0.07);
	Par.hr.GetYaxis()->SetLabelSize(0.07);

	Par.hr.GetXaxis()->SetTitleOffset(1.);
	Par.hr.GetYaxis()->SetTitleOffset(1.);

	Par.hr.SetMaximum(Par.max_y*1.1);
	Par.hr.SetMinimum(Par.min_y/1.2);

	Par.hr.Clone()->Draw();
	Par.hr.Clone()->Draw("SAME AXIS X+ Y+");

	Par.line.Clone()->Draw();

	for (TGraphErrors *graph : Par.VGr) graph->Draw("P");

	Par.legend->SetNColumns(2);
	Par.legend->SetLineColorAlpha(0., 0.);
	Par.legend->SetFillColorAlpha(0., 0.);

	Par.legend->Draw();

	system(("mkdir -p  ../output/Efficiency/" + Par.run_name).c_str());
	PrintCanvas(&canv, "../output/Efficiency/" + Par.run_name + 
		"/Reg" + magf + "_" + Par.Det.title + side + sect + ".png");
	delete Par.legend;
}

void RegEff()
{
	gStyle->SetOptStat(0);
	gROOT->SetBatch(kTRUE);
 	gErrorIgnoreLevel = kWarning;

	Par.line.SetLineColor(kGray+1);
	Par.line.SetLineStyle(2);
	Par.line.SetLineWidth(2);

	for (int i = 0; i < Par.magf.size(); i++)
	{
		for (int j = 0; j < Par.Det.side.size(); j++)
		{
			for (int k = 0; k < Par.Det.sect.size(); k++)
			{
				PerformSingleRegEff(Par.magf[i], Par.Det.side[j], Par.Det.sect[k]);
			}
		}
	}
}
