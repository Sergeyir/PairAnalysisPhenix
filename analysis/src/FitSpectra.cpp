#include "../lib/ErrorHandler.h"
#include "../lib/OutputTool.h"

struct
{
	std::string system_name = "AuAu200";
	std::string part_name = "proton";
	std::string input_dir = "../ext/Spectra/" + system_name;
	//std::string input_dir = "../data/Spectra/Run7AuAu200/19079/";
	//const double mass = 0.139570;
	//const double mass = 0.493677;
	const double mass = 0.938272;
	std::string centr = "0093";
	const double kappa = 2.;
} Par;

void FitSpectra()
{
	gROOT->SetBatch(kTRUE);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	
	std::string input_file_name = Par.input_dir + "/" + Par.part_name + "_" + Par.centr + ".txt";

	CheckInputFile(input_file_name);

	TGraphErrors gr = TGraphErrors(input_file_name.c_str(), "%lg %lg %lg");
	//tsallis distribution
	//https://iopscience.iop.org/article/10.1088/1742-6596/878/1/012016/pdf
	TF1 tsallis = TF1("tsallis", "[0]*x^([5]-1)*(1 + ([1] - 1)*(sqrt([4] + x^2)^([5]) - [2])/[3])^(-1./([1]-1.))");
	
	tsallis.SetParameter(0, 10); //scale
	tsallis.SetParameter(1, 1.01); //q
	tsallis.SetParameter(2, 0.2); //mu
	tsallis.SetParameter(3, 1000); //T
	
	tsallis.SetParLimits(1, 1.00001, 2.);
	tsallis.FixParameter(4, Par.mass*Par.mass); //m2
	tsallis.FixParameter(5, Par.kappa); //kappa

	gr.Fit(&tsallis, "QB");
	
	TCanvas *canv = new TCanvas("canv", "canv", 600, 800);
	gPad->SetLeftMargin(0.12);
	TH1F *pad = gPad->DrawFrame(0.3, 
		gr.GetPointY(gr.GetN()-1)/2., 
		gr.GetPointX(gr.GetN()-1)*1.1, 
		gr.Eval(0.2)*2., 
		(Par.part_name + " spectra").c_str());

	pad->Draw("SAME AXIS X+ Y+");

	pad->GetXaxis()->SetTitle("p_{T}");
	pad->GetYaxis()->SetTitle("d^{2} N / dp_{T} dy");

	gPad->SetLogy();

	gr.SetMarkerStyle(20);
	gr.SetMarkerColor(kRed-3);
	gr.SetLineColor(kRed-3);

	tsallis.SetLineColor(kRed+1);

	gr.Draw("P");

	system(("mkdir -p ../output/fits/" + Par.system_name).c_str());
	canv->SaveAs(("../output/fits/" + Par.system_name + "/" + Par.part_name + ".png").c_str());

	PrintInfo("Tsallis fit parameters: " + 
		to_string(tsallis.GetParameter(0)) + " " +
		to_string(tsallis.GetParameter(1)) + " " +
		to_string(tsallis.GetParameter(2)) + " " +
		to_string(tsallis.GetParameter(3)) + " " +
		to_string(tsallis.GetParameter(4)) + " " +
		to_string(tsallis.GetParameter(5)));
}
