#include <cmath>
#include <fstream>
#include <algorithm>

#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "../lib/BW.h"

#include "ParInv.h"

#include "StrTool.h"
#include "ErrorHandler.h"
#include "OutputTool.h"
#include "Tools.h"
#include "Box.h"
#include "GuiFit.h"
#include "TCanvasPrinter.h"

struct
{
	std::vector<double> pt;
	std::vector<double> raw_yield;
	std::vector<double> yield;
	std::vector<double> stat_err;
	std::vector<double> sys_err;
	std::vector<double> ptmax, ptmin;
} Spectra;

void SetEff()
{	
	if (Par.use_eff)
	{
		std::string ptrange_name = "../par/PtRange/" + Par.run + "/" + 
			Par.particle.name_nl + "/" + Par.method_name + ".txt";
		
		CheckInputFile(ptrange_name);
		ifstream ptrange_file(ptrange_name.c_str());
		
		double ptmin, ptmax;
		while (ptrange_file >> ptmin >> ptmax)
		{
			Par.ptmin.push_back(ptmin);
			Par.ptmax.push_back(ptmax);
		}
		
		for (int i = 0; i < Par.CType.size; i++)
		{
			Par.eff.resize(Par.CType.size);
			Par.eff_err.resize(Par.CType.size);
			Par.eff_pt_sc_sys_err.resize(Par.CType.size);
			Par.eff_acc_sys_err.resize(Par.CType.size);
			Par.eff_m2_eff_sys_err.resize(Par.CType.size);
	
			std::string method_name = "../par/Efficiency/" + Par.run + "/" + 
				Par.particle.name_nl + "/" + Par.method_name + "_" +
				Par.CType.cname_nop[i] + ".txt";
			
			CheckInputFile(method_name.c_str());
			ifstream eff_file(method_name.c_str());
			
			double eff, eff_err, pt_sc_unc, acc_unc, m2_eff_unc;
			
			while (eff_file >> eff >> eff_err >> pt_sc_unc >> acc_unc >> m2_eff_unc)
			{
				Par.eff[i].push_back(eff);
				Par.eff_err[i].push_back(eff_err);
				Par.eff_pt_sc_sys_err[i].push_back(pt_sc_unc);
				Par.eff_acc_sys_err[i].push_back(acc_unc);
				Par.eff_m2_eff_sys_err[i].push_back(m2_eff_unc);
			}
		}
	}
	else
	{
		Par.do_draw_spectra = false;
		Par.do_draw_rcp = false;
		Par.do_draw_rab = false;

		for (double ptmin : DefaultPt.ptmin) Par.ptmin.push_back(ptmin);
		for (double ptmax : DefaultPt.ptmax) Par.ptmax.push_back(ptmax);
	}
}

void SetFrame(TH1F *frame, std::string x_label, std::string y_label, const double size = 0.05)
{
	frame->GetXaxis()->SetTitle(x_label.c_str());
	frame->GetYaxis()->SetTitle(y_label.c_str());
	frame->GetXaxis()->SetTitleSize(size);
	frame->GetYaxis()->SetTitleSize(size);
	frame->GetXaxis()->SetLabelSize(size);
	frame->GetYaxis()->SetLabelSize(size);
}

TGraph GetFGBGGraph(TF1 *func, TH1F *bg_hist, Color_t color, Style_t line_style, const int line_width)
{
	TGraph gr = TGraph();

	for (double m2 = Par.particle.mass - Par.particle.gamma*2.; m2 < Par.particle.mass + Par.particle.gamma*2.; m2 += Par.particle.gamma/50.)
	{
		const double bg_val = bg_hist->GetBinContent(bg_hist->FindBin(m2));
		gr.AddPoint(m2, (func->Eval(m2) + bg_val)/bg_val);
	}

	gr.SetLineStyle(line_style);
	gr.SetMarkerStyle(0);
	gr.SetLineWidth(line_width);
	gr.SetLineColor(color);
	
	return gr;
}

std::string GetTitle(const double ptmin, const double ptmax, std::string cname)
{
	std::string title = DtoStr(ptmin, 1) + " < p_{T} < " + DtoStr(ptmax, 1) + " [GeV/c], " + cname;
	return title;
}

double GetYield(TH1F *hist, TF1 *pol_fit, const double mean, const double gamma)
{
	double yield = 0;
	for (int i = hist->FindBin(mean - gamma*2.); i <= hist->FindBin(mean + gamma*2.); i++)
	{
		yield += hist->GetBinContent(i) - pol_fit->Eval(hist->GetBinCenter(i));
	}
	return yield;
}

double GetYieldErr(TH1F *fg_hist, const double yield, const double mean, const double gamma)
{
	double yield_err = abs(sqrt(fg_hist->Integral(fg_hist->FindBin(mean - 2.*gamma), fg_hist->FindBin(mean + 2.*gamma)))*100./yield);
	return yield_err;
}

void SetHistStyle(TH1F *hist, std::string xtitle, Color_t color, bool do_set_marker = 1, int marker_style = 24, int marker_size = 2)
{
	hist->SetLineColor(color);
	
	hist->GetXaxis()->SetTitle(xtitle.c_str());
	
	hist->SetTitleSize(0.05, "X");
	hist->SetTitleSize(0.05, "Y");

	hist->GetXaxis()->SetLabelSize(0.05);
	hist->GetYaxis()->SetLabelSize(0.05);

	hist->SetLineWidth(2);

	if (do_set_marker == true)
	{
		hist->SetMarkerStyle(marker_style);
		hist->SetMarkerSize(marker_size);
		hist->SetMarkerColor(color);
	}
}

void DrawPartOut(TH1F *hist, const double yield, const double yield_err, const double x = 0.14, const double y = 0.15, Color_t color = kBlack)
{
	TLatex tltext;
	tltext.SetTextFont(52);
	tltext.SetTextSize(0.06);
	tltext.SetTextColor(color);
	tltext.DrawLatexNDC(x, y, (DtoStr(yield, 3) + "/1e6 #pm " + DtoStr(yield_err, 3) + "%").c_str());
}

void DrawFitPar(TF1 *bw_fit)
{
   const double mean = bw_fit->GetParameter(1);
   const double gamma = bw_fit->GetParameter(2);
   const double chi2ndf = bw_fit->GetChisquare()/bw_fit->GetNDF();
      
	TLatex tltext;
	tltext.SetTextFont(52);
	tltext.SetTextSize(0.05);
	tltext.DrawLatexNDC(0.16, 0.82, ("M = " + DtoStr(mean*1e3, 1) + " MeV/c^{2}").c_str());
	tltext.DrawLatexNDC(0.16, 0.78, ("#Gamma = " + DtoStr(gamma*1e3, 2) + " MeV/c^{2}").c_str());
	tltext.DrawLatexNDC(0.58, 0.82, ("#chi^{2}/NDF = " + DtoStr(chi2ndf, 2)).c_str());
}

void SetFitPar(TF1* fg_fit, TF1* bg_fit, const int ptc, std::string bg_par_file_name = "")
{
	fg_fit->SetParameters(10, Par.particle.gamma, Par.particle.mass);
	
	fg_fit->SetParName(0, "K");
	fg_fit->SetParName(1, "M");
	fg_fit->SetParName(2, "#Gamma");
	fg_fit->SetParName(3, "#sigma_{GDF}");

	fg_fit->SetParLimits(1, Par.particle.mass*0.99, Par.particle.mass*1.01);
	fg_fit->FixParameter(3, Par.gdf_sigma.Eval(Average(Par.ptmin[ptc], Par.ptmax[ptc])));

	bg_fit->SetLineStyle(5);

	for (int i = 0; i < bg_fit->GetNpar(); i++)
	{
		bg_fit->SetParameter(i, 1.);
	}

	if (Par.user_defined_bg_fit)
	{
		double pt_tmp[2];
		double ypoints[4];

		ifstream bg_par_file(bg_par_file_name);
		
		pt_tmp[0] = 0;
		bool file_ended = false;
		while (pt_tmp[0] < Par.ptmin[ptc]-0.01 && !file_ended)
		{
			if (!(bg_par_file >> pt_tmp[0] >> pt_tmp[1] >> 
				ypoints[0] >> ypoints[1] >> ypoints[2] >> ypoints[3])) file_ended = true;
		}
		
		bg_fit->SetParameters(ypoints);
		
		for (int i = 0; i < bg_fit->GetNpar(); i++) 
		{
			fg_fit->FixParameter(i+4, bg_fit->GetParameter(i));
		}
	}
}

double PerformBWFit(TH1F *hist, TF1 *bw_fit, TF1 *pol_fit, Color_t color = kRed-3, bool free_Gamma = false, const double alpha = 1.)
{
	bw_fit->SetLineColorAlpha(color, alpha);
	pol_fit->SetLineColorAlpha(color, alpha);

	bw_fit->SetLineWidth(3);
	pol_fit->SetLineWidth(2);
	
	if (!free_Gamma) bw_fit->SetParLimits(2, 
		Par.particle.gamma*Par.min_ssigma, Par.particle.gamma*Par.max_ssigma);
	else bw_fit->SetParLimits(2, Par.particle.gamma/10., Par.particle.gamma*10.);

	if (!Par.user_defined_bg_fit)
	{
		hist->Fit(pol_fit, "RQNCF");
		
		for (int i = 0; i < pol_fit->GetNpar(); i++)
		{
			bw_fit->SetParameter(i+4, pol_fit->GetParameter(i));
		}
	}
	
	hist->Fit(bw_fit, "RQMN");
	
	for (int i = 4; i < bw_fit->GetNpar(); i++) 
	{
		pol_fit->FixParameter(i-4, bw_fit->GetParameter(i));
	}
	
	hist->Fit(pol_fit, "RQBNC");
	
	const double mean = bw_fit->GetParameter(1);
	return mean;
}

void MergeFG(int ptc, int c, int z, int r, std::string channel, TFile *input, TH1F *merged_fg_hist)
{	
	std::string tf_dir = "c0" + to_string(c) + "_z0" + 
		to_string(z) + "_r0" + to_string(r);

	TH2F *fg_hist_2d = (TH2F*) input->Get((tf_dir + "/" + 
		Par.method_name + "_" + channel + Par.addit_name + "_FG12").c_str());

	if (std::isnan(fg_hist_2d->Integral()))
	{
		PrintWarning("Histogram with the name " + 
			static_cast<std::string>(fg_hist_2d->GetName()) + " is corrupted at c = " +
			to_string(c) + ", z = " + to_string(z) + ", r = " + to_string(r));
	}
	
	TH1F *fg_hist = (TH1F*) fg_hist_2d->ProjectionY("FGInvM_proj",
		fg_hist_2d->GetXaxis()->FindBin(Par.ptmin[ptc]+0.01),
		fg_hist_2d->GetXaxis()->FindBin(Par.ptmax[ptc]-0.01));

	fg_hist->Rebin(Par.rebin_num);

	if (merged_fg_hist->GetNbinsX() == 1)
	{
		const int nbins = fg_hist->GetNbinsX();
		const double xmin = fg_hist->GetBinLowEdge(1);
		const double xmax = fg_hist->GetBinLowEdge(nbins) + fg_hist->GetBinWidth(1);
		merged_fg_hist->SetBins(nbins, xmin, xmax);
	}

	merged_fg_hist->Add(fg_hist);
}

void PerformBGSubtr(int ptc, int c, int z, int r, std::string channel, TFile *input,
	TH1F *merged_hist, TH1F *merged_fg_hist, TH1F *merged_bg_hist, const int first_subtr_bin)
{
	std::string tf_dir = "c0" + to_string(c) + "_z0" + 
		to_string(z) + "_r0" + to_string(r);

	TH2F *fg_hist_2d = (TH2F*) input->Get((tf_dir + "/" + 
		Par.method_name + "_" + channel + Par.addit_name + "_FG12").c_str());
	TH2F *bg_hist_2d = (TH2F*) input->Get((tf_dir + "/" + 
		Par.method_name + "_" + channel + Par.addit_name + "_BG12").c_str());

	if (std::isnan(bg_hist_2d->Integral()))
	{
		PrintWarning("Histogram with the name " +
			static_cast<std::string>(bg_hist_2d->GetName()) + " is corrupted at c = " +
			to_string(c) + ", z = " + to_string(z) + ", r = " + to_string(r));
	}

	TH1F *fg_hist = (TH1F*) fg_hist_2d->ProjectionY("FGInvM_proj",
		fg_hist_2d->GetXaxis()->FindBin(Par.ptmin[ptc]+0.01),
		fg_hist_2d->GetXaxis()->FindBin(Par.ptmax[ptc]-0.01));

	TH1F *bg_hist = (TH1F*) bg_hist_2d->ProjectionY("BGInvM_proj",
		fg_hist_2d->GetXaxis()->FindBin(Par.ptmin[ptc]+0.01), 
		fg_hist_2d->GetXaxis()->FindBin(Par.ptmax[ptc]-0.01));

	fg_hist->Rebin(Par.rebin_num);
	bg_hist->Rebin(Par.rebin_num);

	if (merged_fg_hist->GetNbinsX() == 1)
	{
		const int nbins = fg_hist->GetNbinsX();
		const double xmin = fg_hist->GetBinLowEdge(1);
		const double xmax = fg_hist->GetBinLowEdge(nbins) + fg_hist->GetBinWidth(1);
		merged_hist->SetBins(nbins, xmin, xmax);
		merged_fg_hist->SetBins(nbins, xmin, xmax);
		merged_bg_hist->SetBins(nbins, xmin, xmax);
	}
	
	TH1F *hist = (TH1F*) fg_hist->Clone("invM");
	
	bg_hist->Scale(
		fg_hist->Integral(first_subtr_bin, fg_hist->GetXaxis()->GetNbins())/
		bg_hist->Integral(first_subtr_bin, bg_hist->GetXaxis()->GetNbins()));

	hist->Add(bg_hist, -1.);

	const double fg_full_integral = fg_hist->Integral(1, fg_hist->GetXaxis()->GetNbins());
	double fg_corr_integral = 0.;
	double min_val = 1e15;
	int min_val_bin = 0;
	
	//recalculating the scaling for overscaled bg hists
	for (int i = 1; fg_corr_integral < fg_full_integral*Par.fg_corr_part_integral_threshold; i++)
	{
		fg_corr_integral += fg_hist->GetBinContent(i);
		if (hist->GetBinContent(i) < min_val && bg_hist->GetBinContent(i) > 0. && fg_hist->GetBinContent(i) > 0.)
		{
			min_val = hist->GetBinContent(i);
			min_val_bin = i;
		}
	}
	
	if (min_val < 0)
	{
		bg_hist->Scale(fg_hist->GetBinContent(min_val_bin)/bg_hist->GetBinContent(min_val_bin));
		
		hist = (TH1F*) fg_hist->Clone("invM");
		hist->Add(bg_hist, -1.);
	}
	
	//underscaling can be negated by setting bigger Par.fg_norm_part_integral_threshold
	
	merged_hist->Add(hist);
	merged_fg_hist->Add(fg_hist);
	merged_bg_hist->Add(bg_hist);
}

void WriteSpectra(const int cnum)
{
	ofstream spectra_outfile;
	spectra_outfile.open(("../data/Spectra/" + Par.run + "/" + 
	Par.runnum + "/" + Par.method_name + "_" + Par.particle.name_nl + "_" + 
	Par.CType.cname_nop[cnum] + ".txt").c_str());
	
	for (int i = 0; i < Spectra.pt.size(); i++)
	{
		spectra_outfile << Spectra.pt[i] << "	" << 
			Spectra.yield[i] << "	" <<
			Spectra.stat_err[i] << "	" << 
			Spectra.sys_err[i]<< std::endl;
	}
	
	Spectra.pt.clear();
	Spectra.yield.clear();
	Spectra.stat_err.clear();
	Spectra.sys_err.clear();
	Spectra.ptmin.clear();
	Spectra.ptmax.clear();

	spectra_outfile.close();
	PrintInfo("File ../data/Spectra/" + Par.run + "/" + 
	Par.runnum + "/" + Par.method_name + "_" + Par.particle.name_nl + "_" + 
	Par.CType.cname_nop[cnum] + ".txt was written");
}

void WriteRawYields(const int cnum)
{
	ofstream raw_yield_outfile(("../data/RawYields/" + Par.run + "/" + 
	Par.runnum + "/" + Par.method_name + "_" + Par.particle.name_nl + "_" + 
	Par.CType.cname_nop[cnum] + ".txt").c_str());

	for (int i = 0; i < Spectra.pt.size(); i++)
	{
		raw_yield_outfile << Spectra.pt[i] << "	" << 
			Spectra.raw_yield[i] << "	" << 
			Spectra.stat_err[i] << " " << 
			Spectra.sys_err[i] << std::endl;
	}
	raw_yield_outfile.close();
	PrintInfo("File ../data/RawYields/" + Par.run + "/" + 
	Par.runnum + "/" + Par.method_name + "_" + Par.particle.name_nl + "_" + 
	Par.CType.cname_nop[cnum] + ".txt was written");
}

void PerformInvM(TFile *input, const int i, const int ptcmin, const int ptcmax, Box *box)
{
	table.Begin(Par.particle.name_nl + " in " + 
		Par.method_name + "_" + Par.channel + " " + Par.CType.cname[i]);
	table.PrintHeader("pt, GeV", "Y_raw, 10^6", "Y_inv, 10^6", "err, %", "sys, %");

	double nevents = 0;			
	
	for (int c = Par.CType.cmin[i]; c <= Par.CType.cmax[i]; c++)
	{
		for (int z = Par.zmin; z <= Par.zmax; z++)
		{
			for (int r = Par.rmin; r <= Par.rmax; r++)
			{
				std::string tf_dir = "c0" + to_string(c) + "_z0" + 
					to_string(z) + "_r0" + to_string(r);
				
				TH1F *pool_stat = (TH1F*) input->Get((tf_dir + "/PoolStatistics").c_str());
				
				if (std::isnan(pool_stat->GetBinContent(2))) 
				{
					PrintWarning("PoolStatistics histogram is corrupted at c = " + 
					to_string(c) + ", z = " + to_string(z) + ", r = " + to_string(r));
				}
				
				nevents += pool_stat->GetBinContent(2);
			}
		}
	}

	const std::string bg_file_name = 
			"../par/BgFit/" + Par.run + "/" + Par.runnum + "/" +
			Par.particle.name_nl + "/" + Par.method_name + "_" +
			Par.CType.cname_nop[i] + ".txt";
	const std::string bg_file_name_free = 
			"../par/BgFit/" + Par.run + "/" + Par.runnum + "/" +
			Par.particle.name_nl + "/" + Par.method_name + "_" +
			Par.CType.cname_nop[i] + "_free.txt";
	const std::string bg_file_name_wide = 
			"../par/BgFit/" + Par.run + "/" + Par.runnum + "/" +
			Par.particle.name_nl + "/" + Par.method_name + "_" +
			Par.CType.cname_nop[i] + "_wide.txt";
	const std::string bg_file_name_wide_free = 
			"../par/BgFit/" + Par.run + "/" + Par.runnum + "/" +
			Par.particle.name_nl + "/" + Par.method_name + "_" +
			Par.CType.cname_nop[i] + "_wide_free.txt";

	ofstream bg_file;
	ofstream bg_file_free;
	ofstream bg_file_wide;
	ofstream bg_file_wide_free;

	if (Par.user_defined_bg_fit)
	{
		CheckInputFile(bg_file_name.c_str());
		CheckInputFile(bg_file_name_free.c_str());
		CheckInputFile(bg_file_name_wide.c_str());
		CheckInputFile(bg_file_name_wide_free.c_str());
	}         
	else if (Par.write_bg_par)
	{
		system(("mkdir -p ../par/BgFit/" + Par.run + "/" + 
			Par.runnum + "/" + Par.particle.name_nl + "/").c_str());

		bg_file.open(bg_file_name);
		bg_file_free.open(bg_file_name_free);
		bg_file_wide.open(bg_file_name_wide);
		bg_file_wide_free.open(bg_file_name_wide_free);
	}

	box->AddEntry(Par.channel + " " + Par.method_name + " " + Par.CType.cname[i], to_string(nevents/1E6));

	TGraphErrors sys_graph = TGraphErrors();
	std::array<TGraphErrors, 3> sys_ratio = {TGraphErrors(), TGraphErrors(), TGraphErrors()};

	const int nbins = Par.ptmin.size()+1;
	float bins_ranges[nbins];
	
	for (int j = 0; j < nbins; j++)
	{
		bins_ranges[j] = Par.ptmin[j];
	}
	bins_ranges[nbins-1] = Par.ptmax.back();

	std::vector<TH1F> sys_hists = 
	{
		TH1F(("Acceptance" + Par.CType.cname[i]).c_str(), 
			"Acceptance", nbins-1, bins_ranges), 
		TH1F(("p_{T} scale" + Par.CType.cname[i]).c_str(), 
			"p_{T} scale", nbins-1, bins_ranges), 
		TH1F(("m2 eff" + Par.CType.cname[i]).c_str(), 
			"#epsilon_{m^{2}}+#epsilon_{id}", nbins-1, bins_ranges), 
		TH1F(("Raw yield" + Par.CType.cname[i]).c_str(), 
			"Raw yield", nbins-1, bins_ranges),
		TH1F(("BR" + Par.CType.cname[i]).c_str(), 
			"BR", nbins-1, bins_ranges),
		TH1F(("Total" + Par.CType.cname[i]).c_str(), 
			"Total", nbins-1, bins_ranges),
	};

	for (TH1F &hist : sys_hists) hist.SetLineWidth(3);
	
	sys_hists[0].SetLineStyle(2);
	sys_hists[1].SetLineStyle(7);
	sys_hists[2].SetLineWidth(3);
	sys_hists[3].SetLineStyle(2);
	sys_hists[4].SetLineStyle(5);

	sys_hists[3].SetLineWidth(4);
	sys_hists[5].SetLineWidth(4);
	
	sys_graph.SetMarkerStyle(53);
	sys_graph.SetMarkerSize(2);
	sys_graph.SetMarkerColor(kBlack);
	sys_graph.SetLineWidth(2);

	for (int j = 0; j < sys_ratio.size(); j++)
	{
		sys_ratio[j].SetMarkerStyle(53+j);
		sys_ratio[j].SetMarkerSize(2);
		sys_ratio[j].SetLineWidth(2);
	}

	sys_ratio[0].SetMarkerColor(kRed-3);
	sys_ratio[1].SetMarkerColor(kAzure-3);
	sys_ratio[2].SetMarkerColor(kGreen-3);

	sys_ratio[0].SetLineColorAlpha(kRed-3, 0.5);
	sys_ratio[1].SetLineColorAlpha(kAzure-3, 0.5);
	sys_ratio[2].SetLineColorAlpha(kGreen-3, 0.5);
	
	for (int ptc = ptcmin; ptc < ptcmax; ptc++)
	{
		std::string title = GetTitle(Par.ptmin[ptc], Par.ptmax[ptc], Par.CType.cname[i]);

		TH1F *merged_hist = new TH1F((Par.method_name + 
			title).c_str(), title.c_str(), 1, 0, 1);
		TH1F *merged_fg_hist = new TH1F((Par.method_name + 
			title + "fg").c_str(), title.c_str(), 1, 0, 1);
		TH1F *merged_bg_hist = new TH1F((Par.method_name + 
			title + "bg").c_str(), title.c_str(), 1, 0, 1);
		
		const double pt = Average(Par.ptmin[ptc], Par.ptmax[ptc]);
		
		for (int c = Par.CType.cmin[i]; c <= Par.CType.cmax[i]; c++)
		{
			TH1F *current_centr_merged_fg_hist = new TH1F((Par.method_name + 
				title + "fg_current_centr").c_str(), title.c_str(), 1, 0, 1);
			
			//merging fg hists into one for the current centrality class
			for (int z = Par.zmin; z <= Par.zmax; z++)
			{
				for (int r = Par.rmin; r <= Par.rmax; r++)
				{
					MergeFG(ptc, c, z, r, Par.channel, 
						input, current_centr_merged_fg_hist);
					if (Par.particle.has_antipart) 
					{
						MergeFG(ptc, c, z, r, Par.particle.dname2 + Par.particle.dname1, 
							input, current_centr_merged_fg_hist);
					}
				}
			}

			//determining the range of scaling of bg histograms
			int first_subtr_bin = current_centr_merged_fg_hist->GetXaxis()->GetNbins();
			double fg_part_integral = 0.;
			double fg_full_integral = current_centr_merged_fg_hist->Integral(1, first_subtr_bin);
			
			for (;fg_part_integral < fg_full_integral*Par.fg_norm_part_integral_threshold; first_subtr_bin--)
			{
				fg_part_integral += current_centr_merged_fg_hist->GetBinContent(first_subtr_bin);
			}

			delete current_centr_merged_fg_hist;
			
			//subtracting the background
			for (int z = Par.zmin; z <= Par.zmax; z++)
			{
				for (int r = Par.rmin; r <= Par.rmax; r++)
				{
					PerformBGSubtr(ptc, c, z, r, 
						Par.particle.dname2 + Par.particle.dname1, input, 
						merged_hist, merged_fg_hist, merged_bg_hist, first_subtr_bin);
					
					if (Par.particle.has_antipart) 
					{
						PerformBGSubtr(ptc, c, z, r,
							Par.particle.dname2 + Par.particle.dname1, input,
							merged_hist, merged_fg_hist, merged_bg_hist, first_subtr_bin);
					}
				}
			}	
		}
		
		merged_hist->GetXaxis()->SetRange(
			merged_hist->FindBin(Par.particle.zoom_min), 
			merged_hist->FindBin(Par.particle.zoom_max));
		
		TCanvas sum_canv = TCanvas("Summarized_canvas", "Summarized canvas", 1000, 1000);
		
		TH1F *fgbg_hist = (TH1F*) merged_fg_hist->Clone();
		merged_fg_hist->SetTitle(((std::string) "FG, BG, " + merged_fg_hist->GetTitle()).c_str());
		fgbg_hist->SetTitle(((std::string) "FG/BG, " + fgbg_hist->GetTitle()).c_str());
		fgbg_hist->Divide(merged_bg_hist);
		
		fgbg_hist->GetXaxis()->SetRange(1, fgbg_hist->FindBin(Maximum(Par.particle.int_upp, Par.particle.zoom_max)));
		
		TH1F *zoomed_fgbg_hist = (TH1F*) fgbg_hist->Clone();
		zoomed_fgbg_hist->GetXaxis()->SetRange(zoomed_fgbg_hist->FindBin(Par.particle.zoom_min), zoomed_fgbg_hist->FindBin(Par.particle.zoom_max));
		
		SetHistStyle(fgbg_hist, "M_{inv}, GeV/c^{2}", kBlue+2, 1, 20, 1);
		SetHistStyle(merged_fg_hist, "M_{inv}, GeV/c^{2}", kBlue+2, 0);
		SetHistStyle(merged_bg_hist, "M_{inv}, GeV/c^{2}", kRed-3, 0);
		SetHistStyle(zoomed_fgbg_hist, "M_{inv}, GeV/c^{2}", kBlue+2, 1, 20, 1);
      SetHistStyle(merged_hist, "M_{inv}, GeV/c^{2}", kBlue-3, 1, 20, 1);
		
		merged_fg_hist->GetXaxis()->SetRange(merged_fg_hist->FindBin(Par.particle.zoom_min), merged_fg_hist->FindBin(Maximum(Par.particle.int_upp, Par.particle.zoom_max)));
		merged_fg_hist->SetMaximum(merged_fg_hist->GetMaximum()*1.2);
		
		TH1F *hist = (TH1F*) merged_hist->Clone();
		hist->GetXaxis()->SetRange(1, hist->GetNbinsX());
		
		merged_fg_hist->SetFillStyle(3001);
		merged_bg_hist->SetFillStyle(3001);
		hist->SetFillStyle(3001);
		
		hist->SetLineColor(kRed-3);
		merged_bg_hist->SetLineColor(kGreen-3);
		
		hist->SetLineWidth(2);
		
		TLegend legend = TLegend(0.65, 0.65, 0.9, 0.9);
		
		sum_canv.Divide(2);
		
		sum_canv.cd(2);
		gPad->Divide(1, 3);
		gPad->cd(2);
		
		gPad->SetLeftMargin(0.13);
		gPad->SetBottomMargin(0.13);
		gPad->SetGrid();
		gPad->SetLogy();
		fgbg_hist->Draw();
		fgbg_hist->SetMaximum(fgbg_hist->GetBinContent(fgbg_hist->FindBin(Par.particle.mass))*1.2);
		fgbg_hist->SetMinimum(0.95);
		
		sum_canv.cd(2);
		gPad->cd(1);
		
		gPad->SetLeftMargin(0.13);
		gPad->SetBottomMargin(0.13);
		gPad->SetGrid();

		merged_fg_hist->SetMinimum(0);
		merged_fg_hist->Draw("PFC");
		merged_bg_hist->Sumw2(false);
		merged_bg_hist->Draw("PFC SAME");
		hist->Sumw2(false);
		hist->Draw("PFC SAME");

		legend.SetLineColorAlpha(0, 0);
		legend.SetFillColorAlpha(0, 0);
		legend.AddEntry(merged_fg_hist, "FG");
		legend.AddEntry(merged_bg_hist, "BG");
		legend.AddEntry(hist, "FG-BG");
		legend.Draw();
	
		if (merged_hist->Integral(merged_hist->FindBin(Par.particle.mass - Par.particle.gamma*2.), 
			merged_hist->FindBin(Par.particle.mass + Par.particle.gamma*2.)) == 0)
		{
			table.PrintRow(DtoStr(Par.ptmin[ptc], Par.nf) + " - " +
				DtoStr(Par.ptmax[ptc], Par.nf), "empty", "-", "-", "-");
		}
		else
		{	
			double range_min = Par.particle.mass - Par.particle.gamma*2.5;
			double range_max = Par.particle.mass + Par.particle.gamma*2.5;
		
			TF1 bw_fit = TF1("RBW*GDF basic", &RBW_GDF_POL3, 
				Par.particle.mass-Par.particle.gamma*Par.srange, 
				Par.particle.mass+Par.particle.gamma*Par.srange, 8);
			TF1 bw_fit_free = TF1("RBW*GDF free #Gamma", &RBW_GDF_POL3, 
				Par.particle.mass-Par.particle.gamma*Par.srange, 
				Par.particle.mass+Par.particle.gamma*Par.srange, 8);
			TF1 bw_fit_wide = TF1("RBW*GDF wide", &RBW_GDF_POL3, 
				Par.particle.mass-Par.particle.gamma*Par.wide_srange, 
				Par.particle.mass+Par.particle.gamma*Par.wide_srange, 8);
			TF1 bw_fit_wide_free = TF1("RBW*GDF wide free #Gamma", &RBW_GDF_POL3, 
				Par.particle.mass-Par.particle.gamma*Par.wide_srange, 
				Par.particle.mass+Par.particle.gamma*Par.wide_srange, 8);
			
			TF1 pol_fit = TF1("fitPol", "pol3", 
				Par.particle.mass-Par.particle.gamma*Par.srange, 
				Par.particle.mass+Par.particle.gamma*Par.srange);
			TF1 pol_fit_free = TF1("fitPol_free", "pol3", 
				Par.particle.mass-Par.particle.gamma*Par.srange, 
				Par.particle.mass+Par.particle.gamma*Par.srange);
			TF1 pol_fit_wide = TF1("fitPol_wide", "pol3", 
				Par.particle.mass-Par.particle.gamma*Par.wide_srange, 
				Par.particle.mass+Par.particle.gamma*Par.wide_srange);
			TF1 pol_fit_wide_free = TF1("fitPol_wide_free", "pol3", 
				Par.particle.mass-Par.particle.gamma*Par.wide_srange, 
				Par.particle.mass+Par.particle.gamma*Par.wide_srange);

			std::vector<std::thread> set_fits;

			set_fits.push_back(std::thread(SetFitPar, 
				&bw_fit, &pol_fit, ptc, bg_file_name));
			set_fits.push_back(std::thread(SetFitPar, 
				&bw_fit_free, &pol_fit_free, ptc, bg_file_name_free));
			set_fits.push_back(std::thread(SetFitPar, 
				&bw_fit_wide, &pol_fit_wide, ptc, bg_file_name_wide));
			set_fits.push_back(std::thread(SetFitPar, 
				&bw_fit_wide_free, &pol_fit_wide_free, ptc, bg_file_name_wide_free));
			
			while (!set_fits.empty())
			{
				set_fits.back().join();
				set_fits.pop_back();
			}

			const double mean = PerformBWFit(merged_hist,
				&bw_fit, &pol_fit, kRed+1, false, 1);
			const double mean_free = PerformBWFit(merged_hist,
				&bw_fit_free, &pol_fit_free, kRed+1, true, 0.8);
			const double mean_wide = PerformBWFit(merged_hist,
				&bw_fit_wide, &pol_fit_wide, kAzure+1, false, 0.8);
			const double mean_wide_free = PerformBWFit(merged_hist,
				&bw_fit_wide_free, &pol_fit_wide_free, kGreen+1, true, 0.8);
			
			if (!Par.user_defined_bg_fit && Par.write_bg_par)
			{
				bg_file << Par.ptmin[ptc] << " " << Par.ptmax[ptc] << " ";
				for (int i = 0; i < pol_fit.GetNpar()-1; i++)
				{
					bg_file << pol_fit.GetParameter(i) << " ";
				}
				bg_file << pol_fit.GetParameter(pol_fit.GetNpar()-1) << std::endl;

				bg_file_free << Par.ptmin[ptc] << " " << Par.ptmax[ptc] << " ";
				for (int i = 0; i < pol_fit_free.GetNpar()-1; i++)
				{
					bg_file_free << pol_fit_free.GetParameter(i) << " ";
				}
				bg_file_free << pol_fit_free.GetParameter(pol_fit_free.GetNpar()-1) << std::endl;

				bg_file_wide << Par.ptmin[ptc] << " " << Par.ptmax[ptc] << " ";
				for (int i = 0; i < pol_fit_wide.GetNpar()-1; i++)
				{
					bg_file_wide << pol_fit_wide.GetParameter(i) << " ";
				}
				bg_file_wide << pol_fit_wide.GetParameter(pol_fit_wide.GetNpar()-1) << std::endl;

				bg_file_wide_free << Par.ptmin[ptc] << " " << Par.ptmax[ptc] << " ";
				for (int i = 0; i < pol_fit_wide_free.GetNpar()-1; i++)
				{
					bg_file_wide_free << pol_fit_wide_free.GetParameter(i) << " ";
				}
				bg_file_wide_free << pol_fit_wide_free.GetParameter(pol_fit_wide_free.GetNpar()-1) << std::endl;
			}
			
			double yield = GetYield(merged_hist, 
				&pol_fit, mean, Par.particle.gamma);
			double yield_free = GetYield(merged_hist, 
				&pol_fit_free, mean_free, Par.particle.gamma);
			double yield_wide = GetYield(merged_hist, 
				&pol_fit_wide, mean_wide, Par.particle.gamma);
			double yield_wide_free = GetYield(merged_hist, 
				&pol_fit_wide_free, mean_wide_free, Par.particle.gamma);

			double yield_err = GetYieldErr(merged_fg_hist, 
				yield, mean, Par.particle.gamma);
			double yield_free_err = GetYieldErr(merged_fg_hist, 
				yield_free, mean_free, Par.particle.gamma);
			double yield_wide_err = GetYieldErr(merged_fg_hist, 
				yield_wide, mean_wide, Par.particle.gamma);
			double yield_wide_free_err = GetYieldErr(merged_fg_hist, 
				yield_wide_free, mean_wide_free, Par.particle.gamma);
	
			double fit_uncertainty = RMS(
				GetNormRatio(yield/yield_free), 
				GetNormRatio(yield/yield_wide), 
				GetNormRatio(yield/yield_wide_free))-1.;
			
			sys_graph.AddPoint((Par.ptmin[ptc] + Par.ptmax[ptc])/2., fit_uncertainty*100.);
			sys_graph.SetPointError(sys_graph.GetN()-1, 0., yield_err);

			if (Par.draw_sys)
			{
				sys_ratio[0].AddPoint((Par.ptmin[ptc] + Par.ptmax[ptc])/2., yield_free/yield);
				sys_ratio[0].SetPointError(sys_ratio[0].GetN()-1, 0., 
					ErrPropagation(yield_err, yield_free_err)/100.*yield_free/yield);

				sys_ratio[1].AddPoint((Par.ptmin[ptc] + Par.ptmax[ptc])/2., yield_wide/yield);
				sys_ratio[1].SetPointError(sys_ratio[1].GetN()-1, 0., 
					ErrPropagation(yield_err, yield_wide_err)/100.*yield_wide/yield);

				sys_ratio[2].AddPoint((Par.ptmin[ptc] + Par.ptmax[ptc])/2., yield_wide_free/yield);
				sys_ratio[2].SetPointError(sys_ratio[2].GetN()-1, 0., 
					ErrPropagation(yield_err, yield_wide_free_err)/100.*yield_wide_free/yield);
			}

			yield /= nevents;
			
			TCanvas canv = TCanvas("canv", "canv", 800, 800);
			gPad->SetLeftMargin(0.13);
			gPad->SetBottomMargin(0.12);
			gPad->SetTopMargin(0.125);
			
			merged_hist->SetMaximum(merged_hist->GetMaximum()*1.1);
			merged_hist->Draw("P");
         
			if (Par.draw_fit)
			{
				bw_fit.Draw("SAME");
				pol_fit.Draw("SAME");
			}
			
			canv.Update();

			//DrawPartOut(merged_hist, yield*1e6, yield_err);
			DrawFitPar(&bw_fit);

			PrintCanvas(&canv, 
				"../output/InvM/" + Par.run + "/" + Par.runnum + "/" + Par.method_name + 
				"_" + Par.channel + Par.CType.cname_nop[i] +
				"_pt" + DtoStr(Par.ptmin[ptc], Par.nf) + "-" + 
				DtoStr(Par.ptmax[ptc], Par.nf));

			TGraph bw_gr = GetFGBGGraph(&bw_fit, merged_bg_hist, kRed-3, 1, 3);
			TGraph pol_gr = GetFGBGGraph(&pol_fit, merged_bg_hist, kRed-3, 2, 2);

			TText text;
			text.SetTextFont(43);
			text.SetTextSize(30);

			sum_canv.cd(2);
			gPad->cd(3);
			gPad->SetLeftMargin(0.13);
			gPad->SetBottomMargin(0.13);
			zoomed_fgbg_hist->Draw();

			if (Par.draw_fit)
			{
				bw_gr.Clone()->Draw("SAME C");
				pol_gr.Clone()->Draw("SAME C");
			}

			text.DrawText(Par.particle.mass + Par.particle.gamma, 
				zoomed_fgbg_hist->GetMaximum() - (zoomed_fgbg_hist->GetMaximum() - 
				zoomed_fgbg_hist->GetMinimum())*0.1, "Basic fit");
			
			sum_canv.cd(1);
			gPad->Divide(1, 2);
			gPad->cd(2);
			gPad->SetBottomMargin(0.13);
			gPad->SetLeftMargin(0.13);

			merged_hist->Draw("PLC PMC");
			if (Par.draw_fit)
			{
				bw_fit.Draw("SAME");
				pol_fit.Draw("SAME");
			}

			text.DrawText(Par.particle.mass + Par.particle.gamma, 
				merged_hist->GetMaximum()*0.92, "Basic fit");
			
			DrawPartOut(merged_hist, yield*1e6, yield_err);
			
			sum_canv.cd(1);
			gPad->cd(1);
			gPad->SetBottomMargin(0.13);
			gPad->SetLeftMargin(0.13);

			merged_hist->Draw("PLC PMC");

			TLegend fit_legend = TLegend(0.5, 0.65, 0.92, 0.92);
			fit_legend.SetLineColorAlpha(0., 0.);
			fit_legend.SetFillColorAlpha(0., 0.);
			
			bw_fit_free.Draw("SAME");
			pol_fit_free.Draw("SAME");
			bw_fit_wide.Draw("SAME");
			pol_fit_wide.Draw("SAME");
			bw_fit_wide_free.Draw("SAME");
			pol_fit_wide_free.Draw("SAME");

			fit_legend.AddEntry(&bw_fit_free, "Free #Gamma", "L");
			fit_legend.AddEntry(&bw_fit_wide, "Wide", "L");
			fit_legend.AddEntry(&bw_fit_wide_free, "Wide + free #Gamma", "L");

			fit_legend.Draw();

			DrawPartOut(merged_hist, yield_free*1E6/nevents, 
				yield_free_err, 0.14, 0.23, kRed+2);
			DrawPartOut(merged_hist, yield_wide*1E6/nevents, 
				yield_wide_err, 0.14, 0.19, kAzure+2);
			DrawPartOut(merged_hist, yield_wide_free*1E6/nevents, 
				yield_wide_free_err, 0.14, 0.15, kGreen+2);
			
			PrintCanvas(&sum_canv, 
				"../output/InvM/" + Par.run + 
				"/" + Par.runnum + "/SUM_" + Par.method_name + "_" + Par.channel + 
				"_" + Par.CType.cname_nop[i] +
				"_pt" + DtoStr(Par.ptmin[ptc], Par.nf) + "-" +
				DtoStr(Par.ptmax[ptc], Par.nf) + ".png");
				
			if (Par.gui_fit_mode)
			{	
				GuiFitPar.AddHist((TH1F *) merged_hist->Clone(), 
					DtoStr(Par.ptmin[ptc], Par.nf) + " " +
					DtoStr(Par.ptmax[ptc], Par.nf));
				
				GuiFitPar.AddFit(&bw_fit, &pol_fit, &RBW_GDF_POL3, 4);
				GuiFitPar.AddFit(&bw_fit_free, &pol_fit_free, &RBW_GDF_POL3, 4);
				GuiFitPar.AddFit(&bw_fit_wide, &pol_fit_wide, &RBW_GDF_POL3, 4);
				GuiFitPar.AddFit(&bw_fit_wide_free, &pol_fit_wide_free, &RBW_GDF_POL3, 4);
			}
			
			double spectra_norm = 1.;

			Spectra.raw_yield.push_back(yield);
			
			if (Par.use_eff)
			{
				spectra_norm = 2.*pi*(Par.ptmax[ptc]-Par.ptmin[ptc])*pt*
					Par.eff[i][ptc]*Par.particle.BR;
				if (Par.particle.has_antipart == 1) spectra_norm *= 2.;
			}

			double yield_sys_err = 0.;
			if (Par.use_eff)
			{
				yield_sys_err = ErrPropagation(
					fit_uncertainty,
					Par.eff_acc_sys_err[i][ptc], 
					Par.eff_m2_eff_sys_err[i][ptc], 
					Par.eff_pt_sc_sys_err[i][ptc],
					Par.particle.BR_uncertainty);
			}

			table.PrintRow(DtoStr(Par.ptmin[ptc], Par.nf) + " - " +
				DtoStr(Par.ptmax[ptc], Par.nf), to_string(yield*1E6), 
				to_string(yield/spectra_norm), to_string(yield_err), 
				to_string(yield_sys_err*100.));

			Spectra.pt.push_back(pt);
			Spectra.yield.push_back(yield/spectra_norm);
			Spectra.stat_err.push_back(yield_err);
		}	
		delete merged_fg_hist;
		delete merged_bg_hist;
		delete merged_hist;
	}
	
	if (Par.draw_sys)
	{
		TCanvas ratio_canv = TCanvas("ratio", "sys", 1200, 400);
		
		ratio_canv.Divide(2);

		TLine line = TLine(Par.ptmin.front(), 1., Par.ptmax.back(), 1.);
		line.SetLineColor(kGray-1);
		line.SetLineWidth(2);
		line.SetLineStyle(2);
		
		ratio_canv.cd(1);
		gPad->cd(1);
		gPad->SetLeftMargin(0.13);
		gPad->SetBottomMargin(0.13);
		
		TF1 sys_fit = TF1("sys_fit", "expo(0)+expo(2)");
		sys_fit.SetLineColorAlpha(kRed-3, 0.8);
		sys_fit.SetLineStyle(2);
		sys_fit.SetLineWidth(2);
		sys_fit.SetRange(
			Par.ptmin.front() + (Par.ptmax.front() - Par.ptmin.front())/4., 
			Par.ptmax.back() - (Par.ptmax.back() - Par.ptmin.back())/4.);
		
		sys_graph.Fit(&sys_fit, "RQN");

		double sys_graph_max_y = 0.;

		for (int ptc = 0; ptc < Par.ptmin.size(); ptc++)
		{
			double yield_sys_err = 0;
			const double pt = Average(Par.ptmin[ptc], Par.ptmax[ptc]);

			//setting uncertainties
			sys_hists[0].SetBinContent(ptc+1, Par.eff_acc_sys_err[i][ptc]*100.);
			sys_hists[1].SetBinContent(ptc+1, Par.eff_pt_sc_sys_err[i][ptc]*100.);
			sys_hists[2].SetBinContent(ptc+1, Par.eff_m2_eff_sys_err[i][ptc]*100.);
			sys_hists[4].SetBinContent(ptc+1, Par.particle.BR_uncertainty*100.);
			
			yield_sys_err = ErrPropagation(
				sys_fit.Eval(pt)/100.,
				Par.eff_acc_sys_err[i][ptc], 
				Par.eff_m2_eff_sys_err[i][ptc], 
				Par.eff_pt_sc_sys_err[i][ptc],
				Par.particle.BR_uncertainty);
			
			Spectra.sys_err.push_back(yield_sys_err*100.);
			Spectra.ptmin.push_back(Par.ptmin[ptc]);
			Spectra.ptmax.push_back(Par.ptmax[ptc]);

			sys_graph_max_y = Maximum(sys_graph_max_y, sys_fit.Eval(pt), sys_graph.GetPointY(ptc));

			sys_hists[3].SetBinContent(ptc+1, sys_fit.Eval(pt));
			sys_hists[5].SetBinContent(ptc+1, Spectra.sys_err.back());
		}

		TH1F *frame = gPad->DrawFrame(Par.ptmin.front(), 0., Par.ptmax.back(), sys_graph_max_y*1.5, 
			("Raw yield uncertainty"));
		SetFrame(frame, "p_{T}", "Uncertainty, \%", 0.08);

		gPad->SetLeftMargin(0.13);
		gPad->SetBottomMargin(0.18);
		frame->GetYaxis()->SetTitleOffset(0.75);
		
		sys_fit.Draw("SAME");
		sys_graph.Draw("P");

		ratio_canv.cd(2);
		gPad->SetLeftMargin(0.13);
		gPad->SetBottomMargin(0.18);
		
		frame = gPad->DrawFrame(Par.ptmin.front(), 0., Par.ptmax.back(), 3., 
			("Varied/Basic yields ratio"));
		SetFrame(frame, "p_{T}", "Y_{Varied}/Y_{Basic}", 0.08);

      line.DrawClone();

		frame->GetYaxis()->SetTitleOffset(0.7);
		
		TLegend ratio_legend = TLegend(0.6, 0.6, 0.92, 0.92);
		ratio_legend.SetLineColorAlpha(0., 0.);
		ratio_legend.SetFillColorAlpha(0., 0.);
		for (TGraphErrors &gr : sys_ratio) 
		{
			gr.Draw("P");
		}

		ratio_legend.AddEntry(&sys_ratio[0], "Free #Gamma",  "P");
		ratio_legend.AddEntry(&sys_ratio[1], "Wide",  "P");
		ratio_legend.AddEntry(&sys_ratio[2], "Wide + free #Gamma",  "P");

		ratio_legend.Draw();
		
		line.Draw();

		PrintCanvas(&ratio_canv, 
			"../output/InvM/" + Par.run + "/" + 
			Par.runnum + "/RawYieldUnc_" + Par.method_name + "_" + 
			Par.CType.cname_nop[i]);

		TCanvas sys_canv = TCanvas("sys_canv", "sys", 600, 600);
		
		gPad->SetLeftMargin(0.145);
		gPad->SetBottomMargin(0.13);
		
		frame = gPad->DrawFrame(Par.ptmin.front(), 0., Par.ptmax.back(), 
			sys_hists[5].GetMaximum()*1.5,
			("Systematic uncertainties, " + Par.CType.cname[i]).c_str());
		SetFrame(frame, "p_{T}", "Uncertainty, \%");

		TLegend sys_legend(0.3, 0.65, 0.9, 0.88);
		sys_legend.SetNColumns(2);
		
		sys_legend.SetLineColorAlpha(0, 0);
		sys_legend.SetFillColorAlpha(0, 0);

		for (TH1F &hist : sys_hists)
		{
			if (hist.Integral() != 0)
			{
				hist.Draw("SAME PLC ][");
				sys_legend.AddEntry(&hist, hist.GetTitle(), "L");
			}
		}

		sys_legend.Draw();
		
		PrintCanvas(&sys_canv, 
			"../output/InvM/" + Par.run + "/" + 
			Par.runnum + "/SYS_" + Par.method_name + "_" + 
			Par.CType.cname_nop[i]);
	}

	table.End();
	Print();
	
	if (Par.write_yields) 
	{
		WriteRawYields(i);
		WriteSpectra(i);
	}
	
	if (Par.write_bg_par)
	{
		PrintInfo("File " + bg_file_name + " was written");
		PrintInfo("File " + bg_file_name_free + " was written");
		PrintInfo("File " + bg_file_name_wide + " was written");
		PrintInfo("File " + bg_file_name_wide_free + " was written");
	}
}

void PrintParameters()
{
	Box box("Parameters");
	
	box.AddEntry("Taxi directory", Par.taxi_dir);
	box.AddEntry("Run", Par.run);
	box.AddEntry("Run number", Par.runnum);
	box.AddEntry("File", Par.file_name);
	
	box.AddEntry("Particle", Par.particle.name_nl);
	box.AddEntry("Channel", Par.channel);
	
	box.AddEntry("Rebin", to_string(Par.rebin_num));
	box.AddEntry("Sigma normalized yield calculation range", Par.srange);
	box.AddEntry("Sigma normalized range for systematic uncertainty", Par.srange);
		
	box.AddEntry("Method name", Par.method_name);
	box.AddEntry("z range", to_string(Par.zmin) + " - " + to_string(Par.zmax));
	box.AddEntry("r range", to_string(Par.zmin) + " - " + to_string(Par.zmax));
	
	box.AddEntry("Draw spectra", Par.do_draw_spectra);
	box.AddEntry("Draw rcp", to_string(Par.do_draw_rcp));
	box.AddEntry("Draw rab", to_string(Par.do_draw_rab));

	box.Print();
}

void SimpleCheck()
{
	if (Par.ptmin.size() != Par.ptmax.size()) PrintError("ptmin and ptmax sizes are different: " + to_string(Par.ptmin.size()) + " vs " + to_string(Par.ptmax.size()));
	
	for (int i = 0; i < Par.ptmin.size(); i++)
	{
		if (Par.ptmin[i] > Par.ptmax[i]) PrintError("ptmin is bigger than ptmax at position " + to_string(i));
	}
}

unsigned int CalcLines(std::string input_dir)
{
	unsigned int nlines = 0;

	ifstream runs_file("../data/runs.txt");

	std::string subrun_name;

	while (runs_file >> subrun_name)
	{
		CheckInputFile(input_dir + "/" + subrun_name);
		if (subrun_name != (std::string) "sum.root") nlines++;
	}

	if (nlines == 0) PrintError("directory " + input_dir + " is empty or does not exist");

	return nlines;
}
