#pragma once

#include <vector>
#include <array>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TColor.h"
#include "TVirtualPad.h"

#include "ErrorHandler.h"
#include "OutputTool.h"
#include "Tools.h"
#include "StrTool.h"
#include "TCanvasPrinter.h"

using namespace std;

class RABPainter
{
	private:

	const double syst_err_width = 0.05;
	
	std::unique_ptr<TH1F> hr;
	TCanvas canv;
	std::unique_ptr<TLine> tline;
	TLegend legend = TLegend(0.3, 0.78, 1., 0.95);

	TLatex tltext = TLatex();

	std::vector<TGraphErrors> egraphs;
	std::vector<TGraphErrors> sgraphs;

	double syst_x;

	bool new_canv = false;

	std::vector<std::vector<double>> vec_pt;
	std::vector<std::vector<double>> vec_y;
	std::vector<std::vector<double>> vec_y_stat;
	std::vector<std::vector<double>> vec_y_syst;

	public:

   bool printInfo = false;
	
	RABPainter(std::string title, const double minX, const double maxX, std::string XaxisTitle, std::string YaxisTitle, bool create_new_canv = false, double syst_width = 0.05)
	{
		syst_x = syst_width;
		
		if (create_new_canv)
		{
			canv.Constructor(title.c_str(), title.c_str(), 1200, 900);
			new_canv = true;
		}
		
		hr = std::unique_ptr<TH1F>(new TH1F("RAB_range_hist", title.c_str(), 10, minX, maxX));
		tline  = std::unique_ptr<TLine>(new TLine(minX, 1., maxX, 1.));

		gStyle->SetOptStat(0);

		hr->SetTitle(title.c_str());
		hr->GetXaxis()->SetLabelSize(0.05);
		hr->GetYaxis()->SetLabelSize(0.05);
		
		hr->GetXaxis()->SetTitle(XaxisTitle.c_str());
		hr->GetYaxis()->SetTitle(YaxisTitle.c_str());

		hr->GetXaxis()->SetTitleSize(0.06);
		hr->GetYaxis()->SetTitleSize(0.06);

      hr->GetYaxis()->SetTitleOffset(0.75);

		legend.SetNColumns(1);
		legend.SetLineColorAlpha(0., 0.);
		legend.SetFillColorAlpha(0., 0.);

		legend.SetTextFont(43);
		legend.SetTextSize(18);

		tline->SetLineColor(kGray);
		tline->SetLineStyle(2);
		tline->SetLineWidth(3);

		tltext.SetTextFont(53);
		tltext.SetTextColor(kGray+3);
		tltext.SetTextSize(30);

		gPad->SetLeftMargin(0.15);
	}

	void AddGraph(std::string file_name, std::string legend_entry, Color_t color, Style_t marker_style, const double minX = -1, const double maxX = 1e31, bool use_stat = true, bool use_syst = true)
	{
		CheckInputFile(file_name);

		ifstream file(file_name.c_str());
		
      if (use_stat) egraphs.push_back(TGraphErrors());
		if (use_syst) sgraphs.push_back(TGraphErrors());

		double x, y, y_stat_err, y_syst_err;
	
		if (printInfo) Print("From file:", file_name);

		while (file >> x >> y)
		{
         if (use_stat && !(file >> y_stat_err)) break;
         if (use_syst && !(file >> y_syst_err)) break;
         
			if (printInfo) Print(x, y, y_stat_err, y_syst_err);

			if (x < minX) 
			{
				if (printInfo) PrintInfo("x value is less than a minimum threshold : " 
					+ DtoStr(minX, 2) + " ; skipping line");
				continue;
			}
			if (x > maxX) 
			{
				if (printInfo) PrintInfo("x value is greater than a maximum threshold : " 
					+ DtoStr(maxX, 2) + " ; skipping file");
				break;
			}
			if (y <= 0) 
			{
				if (printInfo) PrintWarning("y <= 0; skipping line");
				continue;
			}

			egraphs.back().AddPoint(x, y);
			if (use_stat) egraphs.back().SetPointError(egraphs.back().GetN()-1, 0., y_stat_err);
         
         if (use_syst)
         {
            sgraphs.back().AddPoint(x, y);
            sgraphs.back().SetPointError(sgraphs.back().GetN()-1, syst_x, y_syst_err);
         }
		}

		if (egraphs.back().GetN() == 0) PrintError("Graph is empty");
		
		egraphs.back().SetMarkerStyle(marker_style);
		egraphs.back().SetMarkerColor(color);
		egraphs.back().SetMarkerSize(1.4);
		egraphs.back().SetLineWidth(2);
		egraphs.back().SetLineColor(color);

      if (use_syst)
      {
         sgraphs.back().SetLineColorAlpha(color, 0.5);
         sgraphs.back().SetFillStyle(1001);
         sgraphs.back().SetFillColorAlpha(color, 0.3);
      }
		
		legend.AddEntry(egraphs.back().Clone(), legend_entry.c_str(), "p");
	}

	void AddScalingUncertainty(double x, double syst_y, Color_t color, const double alpha = 0.5, const double width = 0.05)
	{
		sgraphs.push_back(TGraphErrors());

		sgraphs.back().AddPoint(x, 1.);
		sgraphs.back().SetPointError(sgraphs.back().GetN()-1, width, syst_y);

		sgraphs.back().SetLineColorAlpha(0, 0);
		sgraphs.back().SetFillStyle(1001);
		sgraphs.back().SetFillColorAlpha(color, alpha);
	}

	void AddGraph(std::string ab_file_name, std::string pp_file_name, std::string legend_entry, Color_t color, Style_t marker_style, const double ncolls, const double mult, const double minX = -1, const double maxX = 1e31)
	{
		CheckInputFile(ab_file_name);
		CheckInputFile(pp_file_name);

		ifstream ab_file(ab_file_name.c_str());
		ifstream pp_file(pp_file_name.c_str());
		
		egraphs.push_back(TGraphErrors());
		sgraphs.push_back(TGraphErrors());

		double ab_x, ab_y, ab_y_stat_err, ab_y_syst_err;
		double pp_x, pp_y, pp_y_stat_err, pp_y_syst_err;
	
		if (printInfo) Print("From files", ab_file_name, pp_file_name);

		vec_pt.resize(vec_pt.size() + 1);
		vec_y.resize(vec_y.size() + 1);
		vec_y_stat.resize(vec_y_stat.size() + 1);
		vec_y_syst.resize(vec_y_syst.size() + 1);

		pp_file >> pp_x >> pp_y >> pp_y_stat_err >> pp_y_syst_err;

		while (ab_file >> ab_x >> ab_y >> ab_y_stat_err >> ab_y_syst_err)
		{
			if (printInfo) Print("AB:", ab_x, ab_y, ab_y_stat_err, ab_y_syst_err);

			if (ab_x < minX) 
			{
				if (printInfo) PrintInfo("x value is less than a minimum threshold : " 
					+ DtoStr(minX, 2) + " ; skipping AB file line");
				continue;
			}
			if (ab_x > maxX) 
			{
				if (printInfo) PrintInfo("x value is greater than a maximum threshold : " 
					+ DtoStr(maxX, 2) + " ; skipping files");
				break;
			}
			if (ab_y <= 0) 
			{
				if (printInfo) PrintWarning("y <= 0; skipping AB file line");
				continue;
			}

			if (printInfo) Print("pp:", pp_x, pp_y, pp_y_stat_err, pp_y_syst_err);
			
			if (ab_x > pp_x + 0.01)
			{
				while (ab_x > pp_x + 0.01 && pp_file >> pp_x >> pp_y >> pp_y_stat_err >> pp_y_syst_err)
				{	
					if (printInfo) PrintInfo("x from AB file is less than x from pp file; skipping pp file line");
					if (printInfo) Print("pp:", pp_x, pp_y, pp_y_stat_err, pp_y_syst_err);
				}
			}
			
			if (ab_x < pp_x - 0.01)
			{
				if (printInfo) Print(ab_x - pp_x);
				if (printInfo) PrintInfo("x from AB file is greater than x from pp file; skipping AB file line");
				continue;
			}

			const double x = ab_x;
			const double y = ab_y/pp_y/ncolls*mult;
			const double y_stat_err = sqrt(pow(ab_y_stat_err, 2) + pow(pp_y_stat_err, 2))*y/100.;
			const double y_syst_err = sqrt(pow(pp_y_syst_err/100., 2) + pow(ab_y_syst_err/100., 2))*y;
         
         if (printInfo) 
         {
            PrintInfo("Calculating RAB at pt_ab = " + DtoStr(ab_x, 2) + ", pt_pp = " + DtoStr(pp_x, 2) + ": " + to_string(y));
            if (y < 0. || y > 2.) PrintWarning("RAB at pt = " + DtoStr(ab_x, 2) + " is outside of y axis range");
         }

			egraphs.back().AddPoint(x, y);
			egraphs.back().SetPointError(egraphs.back().GetN()-1, 0., y_stat_err);

			sgraphs.back().AddPoint(x, y);
			sgraphs.back().SetPointError(sgraphs.back().GetN()-1, syst_x, y_syst_err);

			vec_pt.back().push_back(x);
			vec_y.back().push_back(y);
			vec_y_stat.back().push_back(y_stat_err*100.);
			vec_y_syst.back().push_back(y_syst_err*100.);
		}

		if (egraphs.back().GetN() == 0) PrintError("Graph is empty");
		
		egraphs.back().SetMarkerStyle(marker_style);
		egraphs.back().SetMarkerColor(color);
		egraphs.back().SetMarkerSize(1.4);
		egraphs.back().SetLineWidth(2);
		egraphs.back().SetLineColor(color);

		sgraphs.back().SetLineColorAlpha(color, 0.5);
		sgraphs.back().SetFillStyle(1001);
		sgraphs.back().SetFillColorAlpha(color, 0.1);
		
		legend.AddEntry(egraphs.back().Clone(), legend_entry.c_str(), "p");
	}

	void ExtendGraph(std::string ab_file_name, std::string pp_file_name, const double ncolls, const double mult, const double minX = -1, const double maxX = 1e31)
	{
		CheckInputFile(ab_file_name);
		CheckInputFile(pp_file_name);

		ifstream ab_file(ab_file_name.c_str());
		ifstream pp_file(pp_file_name.c_str());
		
		double ab_x, ab_y, ab_y_stat_err, ab_y_syst_err;
		double pp_x, pp_y, pp_y_stat_err, pp_y_syst_err;
	
		if (printInfo) Print("From files", ab_file_name, pp_file_name);

		pp_file >> pp_x >> pp_y >> pp_y_stat_err >> pp_y_syst_err;

		while (ab_file >> ab_x >> ab_y >> ab_y_stat_err >> ab_y_syst_err)
		{
			if (printInfo) Print("AB:", ab_x, ab_y, ab_y_stat_err, ab_y_syst_err);

			if (ab_x < minX) 
			{
				if (printInfo) PrintInfo("x value is less than a minimum threshold : " 
					+ DtoStr(minX, 2) + " ; skipping AB file line");
				continue;
			}
			if (ab_x > maxX) 
			{
				if (printInfo) PrintInfo("x value is greater than a maximum threshold : " 
					+ DtoStr(maxX, 2) + " ; skipping files");
				break;
			}
			if (ab_y <= 0) 
			{
				if (printInfo) PrintWarning("y <= 0; skipping AB file line");
				continue;
			}

			if (printInfo) Print("pp:", pp_x, pp_y, pp_y_stat_err, pp_y_syst_err);
			
			if (ab_x > pp_x + 0.01)
			{
				while (ab_x > pp_x + 0.01 && pp_file >> pp_x >> pp_y >> pp_y_stat_err >> pp_y_syst_err)
				{	
               if (printInfo)
               {
                  PrintInfo("x from AB file is less than x from pp file; skipping pp file line");
                  Print("pp:", pp_x, pp_y, pp_y_stat_err, pp_y_syst_err);
               }
				}
			}
			
			if (ab_x < pp_x - 0.01)
			{
            if (printInfo)
            {
               Print(ab_x - pp_x);
               PrintInfo("x from AB file is greater than x from pp file; skipping AB file line");
            }
				continue;
			}

			if (printInfo) PrintInfo("Calculating RAB at pt_ab = " + DtoStr(ab_x, 2) + ", pt_pp = " + DtoStr(pp_x, 2));

			const double x = ab_x;
			const double y = ab_y/pp_y/ncolls*mult;
			const double y_stat_err = sqrt(pow(ab_y_stat_err, 2) + pow(pp_y_stat_err, 2))*y/100.;
			const double y_syst_err = sqrt(pow(pp_y_syst_err/100., 2) + pow(ab_y_syst_err/100., 2))*y;

			egraphs.back().AddPoint(ab_x, y);
			egraphs.back().SetPointError(egraphs.back().GetN()-1, 0., y_stat_err);

			sgraphs.back().AddPoint(ab_x, y);
			sgraphs.back().SetPointError(sgraphs.back().GetN()-1, syst_x, y_syst_err);

			vec_pt.back().push_back(x);
			vec_y.back().push_back(y);
			vec_y_stat.back().push_back(y_stat_err*100.);
			vec_y_syst.back().push_back(y_syst_err*100.);
		}
	}
	
	void Draw(std::string name, bool draw_opposite_axis = true, bool draw_x_labels = true, bool draw_y_labels = true, bool draw_legend = true, bool draw_syst = true)
	{
		hr->SetMinimum(0.);
		hr->SetMaximum(2.);

      TH1F *current_hr = (TH1F*) hr->Clone();

      if (!draw_x_labels) 
      {
         current_hr->GetXaxis()->SetTitle("");
         current_hr->GetXaxis()->SetLabelSize(0);
      }

      if (!draw_y_labels) 
      {
         current_hr->GetYaxis()->SetTitle("");
         current_hr->GetYaxis()->SetLabelSize(0);
      }

		current_hr->Draw("AXIS");
		if (draw_opposite_axis) current_hr->DrawClone("SAME AXIS X+ Y+");	

		if (egraphs.size() != 0) tline->DrawClone();

		for (TGraphErrors gr : egraphs) gr.DrawClone("P");
		if (draw_syst) for (TGraphErrors gr : sgraphs) gr.DrawClone("5");
		if (draw_legend) legend.DrawClone();
		if (!new_canv)
		{
			if (egraphs.size() != 0) tltext.DrawLatex(current_hr->GetBinLowEdge(2), 0.2, name.c_str());
			else tltext.DrawLatex(current_hr->GetBinLowEdge(2), 1.6, name.c_str());
		}

		egraphs.clear();
		sgraphs.clear();

		legend.Clear();
		
		if (new_canv) PrintCanvas(&canv, name);
	}

	void PrintTable(std::string output_file_name)
	{
		CheckOutputFile(output_file_name);
		ofstream output(output_file_name);
		
		std::cout.precision(6);
		output << "\\begin{table}[h!]" << std::endl;
		output << "\\centering" << std::endl;
		output << "\\begin{tabular}{c|c|c|c}" << std::endl;
		output << "$p_T$, GeV/c" << " & " << "$R_{AA}$" << " & " 
			<< "Statistical uncertainty, \\%" << " & " << "Systematic uncertainty, \\% \\\\" << std::endl;
		for (int j = 0; j < vec_pt.back().size(); j++)
		{
			output << "	" << vec_pt.back()[j] << " & " << 
				std::scientific << vec_y.back()[j] << std::defaultfloat << " & " 
				<< vec_y_stat.back()[j] << " & " << vec_y_syst.back()[j] << "\\\\" << std::endl;
		}
		output << "\\end{tabular}" << std::endl;
		output << "\\caption{\\Kstar meson $R_{AA}$ in Au+Au at \\snn}" << std::endl;
		output << "\\end{table}" << std::endl;
		
		output.close();
		PrintInfo("File " + output_file_name + " was written");

		vec_pt.pop_back();
		vec_y.pop_back();
		vec_y_stat.pop_back();
		vec_y_syst.pop_back();
	}
};
