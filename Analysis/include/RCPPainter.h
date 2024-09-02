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

using namespace std;

class RCPPainter
{
	private:

	const double sys_err_width = 0.05;
	
	std::unique_ptr<TH1F> hr;
	TCanvas canv;
	std::unique_ptr<TLine> tline;
	TLegend legend = TLegend(0.15, 0.65, 0.85, 0.85);

	TLatex tltext = TLatex();

	std::vector<TGraphErrors> egraphs;
	std::vector<TGraphErrors> sgraphs;

	double sys_x;

	bool new_canv = false;

	std::vector<double> vec_pt;
	std::vector<double> vec_y;
	std::vector<double> vec_y_stat;
	std::vector<double> vec_y_sys;

   bool printInfo = false;

	public:
	
	RCPPainter(std::string title, const double minX, const double maxX, std::string XaxisTitle, std::string YaxisTitle, bool create_new_canv = false, double sys_width = 0.03)
	{
		sys_x = sys_width;
		
		if (create_new_canv)
		{
			canv.Constructor(title.c_str(), title.c_str(), 1200, 900);
			new_canv = true;
		}
		
		hr = std::unique_ptr<TH1F>(new TH1F("RCP_range_hist", title.c_str(), 10, minX, maxX));
		tline  = std::unique_ptr<TLine>(new TLine(minX, 1., maxX, 1));

		gStyle->SetOptStat(0);

		hr->SetTitle(title.c_str());
		hr->GetXaxis()->SetLabelSize(0.05);
		hr->GetYaxis()->SetLabelSize(0.05);
		
		hr->GetXaxis()->SetTitle(XaxisTitle.c_str());
		hr->GetYaxis()->SetTitle(YaxisTitle.c_str());

		hr->GetXaxis()->SetTitleSize(0.06);
		hr->GetYaxis()->SetTitleSize(0.06);

		legend.SetLineColorAlpha(0., 0.);
		legend.SetFillColorAlpha(0., 0.);

		legend.SetTextFont(43);
		legend.SetTextSize(20);

		tline->SetLineColor(kGray);
		tline->SetLineStyle(2);
		tline->SetLineWidth(3);

		tltext.SetTextFont(53);
		tltext.SetTextColor(kGray+3);
		tltext.SetTextSize(30);

		gPad->SetLeftMargin(0.15);
	}

	void AddScalingUncertainty(double x, double sys_y, Color_t color, const double alpha = 0.5, const double width = 0.05)
	{
		sgraphs.push_back(TGraphErrors());

		sgraphs.back().AddPoint(x, 1.);
		sgraphs.back().SetPointError(sgraphs.back().GetN()-1, width, sys_y);

		sgraphs.back().SetLineColorAlpha(0, 0);
		sgraphs.back().SetFillStyle(1001);
		sgraphs.back().SetFillColorAlpha(color, alpha);
	}

	//only for pions, kaons, and protons
	void AddGraph(std::string file_name, std::string legend_entry, Color_t color, Style_t marker_style, const double minX = -1, const double maxX = 1e31)
	{
		CheckInputFile(file_name);

		ifstream file(file_name.c_str());
		
		egraphs.push_back(TGraphErrors());
		sgraphs.push_back(TGraphErrors());

		double x, p_y, p_y_stat_err, a_y, a_y_stat_err, single_y_sys_err;
	
		if (printInfo) Print("From file", file_name);

		while (file >> x >> p_y >> p_y_stat_err  >> a_y >> a_y_stat_err >> single_y_sys_err)
		{
			if (printInfo) Print(x, p_y, p_y_stat_err, a_y, a_y_stat_err);
			
			const double y = (a_y+p_y)/2.;
			const double y_stat_err = sqrt(pow(a_y_stat_err/a_y, 2) + pow(p_y_stat_err/p_y, 2))*y;
			const double y_sys_err = single_y_sys_err/100.*y*sqrt(2.);

			egraphs.back().AddPoint(x, y);
			egraphs.back().SetPointError(egraphs.back().GetN()-1, 0., y_stat_err);

			sgraphs.back().AddPoint(x, y);
			sgraphs.back().SetPointError(sgraphs.back().GetN()-1, sys_x, y_sys_err);
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

	void AddGraph(std::string c_file_name, std::string p_file_name, std::string legend_entry, Color_t color, Style_t marker_style, const double c_ncolls, const double p_ncolls, const double minX = -1, const double maxX = 1e31, bool has_systematic = true)
	{
		CheckInputFile(c_file_name);
		CheckInputFile(p_file_name);

		ifstream c_file(c_file_name.c_str());
		ifstream p_file(p_file_name.c_str());
		
		egraphs.push_back(TGraphErrors());
		sgraphs.push_back(TGraphErrors());

		double c_x, c_y, c_y_stat_err, c_y_sys_err;
		double p_x, p_y, p_y_stat_err, p_y_sys_err;
	
		if (printInfo) Print("From files", c_file_name, p_file_name);

		vec_pt.clear();
		vec_y.clear();
		vec_y_stat.clear();
		vec_y_sys.clear();

		p_file >> p_x >> p_y >> p_y_stat_err; 
		if (has_systematic) p_file >> p_y_sys_err;
		else p_y_sys_err = 0;

		while (c_file >> c_x >> c_y >> c_y_stat_err)
		{
			if (has_systematic) c_file >> c_y_sys_err;
			else c_y_sys_err = 0;
			if (printInfo) Print("central:", c_x, c_y, c_y_stat_err, c_y_sys_err);

			if (c_x < minX) 
			{
				if (printInfo) PrintInfo("x value is less than a minimum threshold : " 
					+ DtoStr(minX, 2) + " ; skiping central file line");
				continue;
			}
			if (c_x > maxX) 
			{
				if (printInfo) PrintInfo("x value is greater than a maximum threshold : " 
					+ DtoStr(maxX, 2) + " ; skiping files");
				break;
			}
			if (c_y <= 0) 
			{
				if (printInfo) PrintWarning("y <= 0; skiping central file line");
				continue;
			}

			if (printInfo) Print("p:", p_x, p_y, p_y_stat_err, p_y_sys_err);
			
			if (c_x > p_x + 0.01)
			{
				while (c_x > p_x + 0.01 && p_file >> p_x >> p_y >> p_y_stat_err >> p_y_sys_err)
				{	
               if (printInfo)
               {
                  PrintInfo("x from central file is less than x from p file; skiping p file line");
                  Print("p:", p_x, p_y, p_y_stat_err, p_y_sys_err);
               }
				}
			}
			
         if (c_x < p_x - 0.01)
         {
            if (printInfo)
            {
               Print(c_x - p_x);
               PrintInfo("x from central file is greater than x from p file; skiping central file line");
            }
            continue;
         }

			const double x = c_x;
			const double y = c_y/p_y/c_ncolls*p_ncolls;
			const double y_stat_err = sqrt(pow(c_y_stat_err, 2) + pow(p_y_stat_err, 2))*y/100.;
			const double y_sys_err = sqrt(pow(p_y_sys_err/100., 2) + pow(c_y_sys_err/100., 2))*y;

         if (printInfo)
         {
            PrintInfo("Calculating central at pt_c = " + DtoStr(c_x, 2) + ", pt_p = " + DtoStr(p_x, 2) + ": " + to_string(y));
            if (y < 0. || y > 2.) PrintWarning("central at pt = " + DtoStr(c_x, 2) + " is outside of y axis range");
         }

			egraphs.back().AddPoint(x, y);
			egraphs.back().SetPointError(egraphs.back().GetN()-1, 0., y_stat_err);

			sgraphs.back().AddPoint(x, y);
			sgraphs.back().SetPointError(sgraphs.back().GetN()-1, sys_x, y_sys_err);

			vec_pt.push_back(x);
			vec_y.push_back(y);
			vec_y_stat.push_back(y_stat_err*100.);
			vec_y_sys.push_back(y_sys_err*100.);
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

	void ExtendGraph(std::string c_file_name, std::string p_file_name, const double c_ncolls, const double p_ncolls, const double minX = -1, const double maxX = 1e31)
	{
		CheckInputFile(c_file_name);
		CheckInputFile(p_file_name);

		ifstream c_file(c_file_name.c_str());
		ifstream p_file(p_file_name.c_str());
		
		double c_x, c_y, c_y_stat_err, c_y_sys_err;
		double p_x, p_y, p_y_stat_err, p_y_sys_err;
	
		if (printInfo) Print("From files", c_file_name, p_file_name);

		p_file >> p_x >> p_y >> p_y_stat_err >> p_y_sys_err;

		while (c_file >> c_x >> c_y >> c_y_stat_err >> c_y_sys_err)
		{
			if (printInfo) Print("AB:", c_x, c_y, c_y_stat_err, c_y_sys_err);

			if (c_x < minX) 
			{
				if (printInfo) PrintInfo("x value is less than a minimum threshold : " 
					+ DtoStr(minX, 2) + " ; skiping AB file line");
				continue;
			}
			if (c_x > maxX) 
			{
				if (printInfo) PrintInfo("x value is greater than a maximum threshold : " 
					+ DtoStr(maxX, 2) + " ; skiping files");
				break;
			}
			if (c_y <= 0) 
			{
				if (printInfo) PrintWarning("y <= 0; skiping AB file line");
				continue;
			}

			if (printInfo) Print("p:", p_x, p_y, p_y_stat_err, p_y_sys_err);
			
			if (c_x > p_x + 0.01)
			{
				while (c_x > p_x + 0.01 && p_file >> p_x >> p_y >> p_y_stat_err >> p_y_sys_err)
				{	
               if (printInfo)
               {
                  PrintInfo("x from AB file is less than x from p file; skiping p file line");
                  Print("p:", p_x, p_y, p_y_stat_err, p_y_sys_err);
               }
				}
			}
			
			if (c_x < p_x - 0.01)
			{
            if (printInfo)
            {
               Print(c_x - p_x);
               PrintInfo("x from AB file is greater than x from p file; skiping AB file line");
            }
				continue;
			}

			if (printInfo) PrintInfo("Calculating central at pt_c = " + DtoStr(c_x, 2) + ", pt_p = " + DtoStr(p_x, 2));

			const double x = c_x;
			const double y = c_y/p_y/c_ncolls*p_ncolls;
			const double y_stat_err = sqrt(pow(c_y_stat_err, 2) + pow(p_y_stat_err, 2))*y/100.;
			const double y_sys_err = sqrt(pow(p_y_sys_err/100., 2) + pow(c_y_sys_err/100., 2))*y;

			egraphs.back().AddPoint(c_x, y);
			egraphs.back().SetPointError(egraphs.back().GetN()-1, 0., y_stat_err);

			sgraphs.back().AddPoint(c_x, y);
			sgraphs.back().SetPointError(sgraphs.back().GetN()-1, sys_x, y_sys_err);

			vec_pt.push_back(x);
			vec_y.push_back(y);
			vec_y_stat.push_back(y_stat_err*100.);
			vec_y_sys.push_back(y_sys_err*100.);
		}
	}
	
	void Draw(std::string name, bool draw_legend = true, bool draw_sys = true)
	{
		hr->SetMinimum(0.);
		hr->SetMaximum(2.);

		hr->Clone()->Draw();	
		hr->Clone()->Draw("SAME AXIS X+ Y+");	

		if (egraphs.size() != 0) tline->Clone()->Draw();

		for (TGraphErrors gr : egraphs) gr.Clone()->Draw("P");
		if (draw_sys) for (TGraphErrors gr : sgraphs) gr.Clone()->Draw("5");
		if (draw_legend) legend.DrawClone();
		if (!new_canv)
		{
			if (egraphs.size() != 0) tltext.DrawLatex(hr->GetBinLowEdge(2), 0.2, name.c_str());
			else tltext.DrawLatex(hr->GetBinLowEdge(2), 1.6, name.c_str());
		}

		egraphs.clear();
		sgraphs.clear();

		legend.Clear();
		
		if (new_canv) canv.SaveAs(name.c_str());
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
		for (int j = 0; j < vec_pt.size(); j++)
		{
			output << "	" << vec_pt[j] << " & " << 
				std::scientific << vec_y[j] << std::defaultfloat << " & " 
				<< vec_y_stat[j] << " & " << vec_y_sys[j] << "\\\\" << std::endl;
		}
		output << "\\end{tabular}" << std::endl;
		output << "\\caption{\\Kstar meson $R_{CP}$ in Au+Au at \\snn}" << std::endl;
		output << "\\end{table}" << std::endl;
		
		output.close();
		PrintInfo("File " + output_file_name + " was written");

		vec_pt.clear();
		vec_y.clear();
		vec_y_stat.clear();
		vec_y_sys.clear();
	}
};
