#include "../lib/Particles.h"
#include "../lib/ErrorHandler.h"
#include "../lib/OutputTool.h"
#include "../lib/StrTool.h"
#include "../lib/Tools.h"
#include "../lib/TCanvasPrinter.h"

struct
{
	std::string run_name = "Run7AuAu200";
	KStar892 part;
	const double ptmin = 0.9;
	const double ptmax = 7.5;
	const double pt_step = 0.4;
	
	const double sigma_min = 1e-4;
	const double sigma_max = 15e-3;

	std::string sigma_fit_func = "pol1";

	const int nsteps = 5;
} Par;

void Width()
{
	gROOT->SetBatch(kTRUE);
	gStyle->SetOptStat(0);
	gErrorIgnoreLevel = kWarning;
	
	std::string input_file_name = "../data/phenix_sim/" + Par.run_name + "/Widthless_" + Par.part.name_nl + ".root";
	CheckInputFile(input_file_name);

	system(("mkdir -p ../output/InvM/" + Par.run_name + "/Widthless").c_str());
	
	TFile input_file = TFile(input_file_name.c_str());
	TH2F *invm = (TH2F *) input_file.Get("InvM");

	TGraph sigma_gr = TGraph();
	sigma_gr.SetMarkerStyle(71);
	sigma_gr.SetMarkerSize(2);
	sigma_gr.SetMarkerColor(kRed-3);

	double max_sigma = 0.;

	for (double pt = Par.ptmin + Par.pt_step/2.; pt < Par.ptmax; pt += Par.pt_step)
	{
		TH1F *invm_proj = (TH1F *) invm->ProjectionY(DtoStr(pt, 1).c_str(),
			invm->GetXaxis()->FindBin(pt - Par.pt_step/2. + 0.01), 
			invm->GetXaxis()->FindBin(pt + Par.pt_step/2. - 0.01));
		
		TF1 gaus_fit = TF1("gaus", "gaus");
		gaus_fit.SetParameter(0, invm_proj->GetBinContent(invm_proj->GetXaxis()->FindBin(Par.part.mass)));
		gaus_fit.SetParameter(1, Par.part.mass);
		gaus_fit.SetParLimits(1, Par.part.mass - 5e-3, Par.part.mass + 5e-3);
		gaus_fit.SetParameter(2, 1.);
		gaus_fit.SetParLimits(2, Par.sigma_min, Par.sigma_max);

		gaus_fit.SetRange(Par.part.mass - Par.sigma_max*2., Par.part.mass + Par.sigma_max*2.);
		
		TCanvas canv = TCanvas("canv", "canv", 800, 800);
		for (int j = 0; j < Par.nsteps; j++)
		{
			invm_proj->Fit(&gaus_fit, "RQMBN");
			gaus_fit.SetRange(
				Par.part.mass - 2.5*gaus_fit.GetParameter(2),
				Par.part.mass + 2.5*gaus_fit.GetParameter(2));
		}
		
		invm_proj->GetXaxis()->SetRange(
			invm_proj->GetXaxis()->FindBin(Par.part.mass - 5.*gaus_fit.GetParameter(2)),
			invm_proj->GetXaxis()->FindBin(Par.part.mass + 5.*gaus_fit.GetParameter(2)));

		invm_proj->SetMaximum(invm_proj->GetBinContent(invm_proj->GetXaxis()->FindBin(Par.part.mass))*2.);
		
		invm_proj->Fit(&gaus_fit, "RQMB");
		sigma_gr.AddPoint(pt, gaus_fit.GetParameter(2));

		PrintCanvas(&canv, "../output/InvM/" + Par.run_name + "/Widthless/" + 
			DtoStr(pt - Par.pt_step/2., 1) + "-" + DtoStr(pt + Par.pt_step/2., 1), false);

		max_sigma = Maximum(max_sigma, gaus_fit.GetParameter(2));
	}

	TF1 sigma_fit = TF1("sigma_fit", Par.sigma_fit_func.c_str());
	sigma_fit.SetLineColor(kRed+2);
	sigma_fit.SetLineWidth(3);
	sigma_fit.SetLineStyle(2);

	TCanvas canv = TCanvas("canv", "canv", 800, 600);
	gPad->SetLeftMargin(0.135);
	gPad->SetBottomMargin(0.13);

	TH1F *hr = gPad->DrawFrame(Par.ptmin, 0., Par.ptmax, max_sigma*1.5);
	hr->GetXaxis()->SetTitle("p_{T}, GeV/c");
	hr->GetYaxis()->SetTitle("#sigma, (GeV/c^{2})^{2}");
	hr->GetXaxis()->SetTitleOffset(1.);
	hr->GetYaxis()->SetTitleOffset(1.4);
	hr->GetXaxis()->SetLabelSize(0.05);
	hr->GetYaxis()->SetLabelSize(0.05);
	hr->SetTitleSize(0.05, "X");
	hr->SetTitleSize(0.05, "Y");
	
	sigma_gr.Fit(&sigma_fit, "QMB");
	sigma_gr.DrawClone("P");
	
	PrintCanvas(&canv, "../output/InvM/" + Par.run_name + "/Widthless/sigmas");
	
	Print(sigma_fit.GetParameter(0), sigma_fit.GetParameter(1));
	system(("mkdir -p ../par/GDF/" + Par.run_name).c_str());
	
	WriteFile("../par/GDF/" + Par.run_name + "/" + Par.part.name_nl + ".txt", 
		sigma_fit.GetParameter(0), sigma_fit.GetParameter(1));
}
