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

class SpectraPainter
{
	private:

	const double sys_err_width = 0.05;
	
	double min_x;
	double max_x;
	double max_y{-1e30};
	double min_y{1e30};

	std::unique_ptr<TLine> ratio_line;
	std::unique_ptr<TH1F> hr;
	TCanvas canv;
	TLegend legend = TLegend(0.15, 0.12, 0.75, 0.38);

	std::vector<std::string> legend_entries;

	std::vector<TGraphErrors> egraphs;
	std::vector<TGraphErrors> sgraphs;
		
	std::vector<TGraphErrors> ratio_egraphs;
	std::vector<TGraphErrors> ratio_sgraphs;
	std::vector<double> mult_factor;

	std::vector<TF1> levy_funcs;
	std::vector<TF1> exp_funcs;
	
	std::vector<std::vector<double>> vec_pt;
	std::vector<std::vector<double>> vec_y;
	std::vector<std::vector<double>> vec_y_stat;
	std::vector<std::vector<double>> vec_y_sys;

   bool printInfo = false;
   
   public:
	
	SpectraPainter(std::string title, const double minX, const double maxX, std::string XaxisTitle, std::string YaxisTitle, const int canv_size_x = 700, const int canv_size_y = 900)
	{
		canv.Constructor(title.c_str(), title.c_str(), canv_size_x, canv_size_y);
		hr = std::make_unique<TH1F>(TH1F("Spectra_range_hist", title.c_str(), 10, minX, maxX));
		ratio_line = std::make_unique<TLine>(TLine(minX, 1., maxX, 1.));
		min_x = minX;
		max_x = maxX;

		canv.SetLogy();
		gStyle->SetOptStat(0);

		hr->SetTitle(title.c_str());
		hr->GetXaxis()->SetLabelSize(0.05);
		hr->GetYaxis()->SetLabelSize(0.05);
		
		hr->GetXaxis()->SetTitle(XaxisTitle.c_str());
		hr->GetYaxis()->SetTitle(YaxisTitle.c_str());

		hr->GetXaxis()->SetTitleSize(0.05);
		hr->GetYaxis()->SetTitleSize(0.05);

		legend.SetLineColorAlpha(0., 0.);
		legend.SetFillColorAlpha(0., 0.);
		
		ratio_line->SetLineColor(kGray);
		ratio_line->SetLineStyle(2);
		ratio_line->SetLineWidth(3);

		gPad->SetLeftMargin(0.17);
		gPad->SetBottomMargin(0.11);
	}

	TF1 *AddGraph(std::string spectra_file_name, std::string legend_entry, const int pow10, double mult, const double particle_mass, Style_t marker_style, Color_t color, const double minX = -1, const double maxX = 1e31, bool add_percent = false)
	{
		CheckInputFile(spectra_file_name);

		ifstream input_file(spectra_file_name.c_str());
		
		egraphs.push_back(TGraphErrors());
		sgraphs.push_back(TGraphErrors());

		ratio_egraphs.push_back(TGraphErrors());
		ratio_sgraphs.push_back(TGraphErrors());

		double x, y, y_stat_err, y_sys_err;
	
		mult *= pow(10., pow10);
		if (printInfo) PrintInfo("Mult is " + to_string(mult));

		mult_factor.push_back(mult);
		
		if (printInfo) Print("From file:", spectra_file_name);

		vec_pt.resize(vec_pt.size() + 1);
		vec_y.resize(vec_y.size() + 1);
		vec_y_stat.resize(vec_y_stat.size() + 1);
		vec_y_sys.resize(vec_y_sys.size() + 1);
		
		while (input_file >> x >> y >> y_stat_err >> y_sys_err)
		{
			if (printInfo) Print(x, y, y_stat_err, y_sys_err);
			if (x < minX) 
			{
				if (printInfo) PrintInfo("x value is less than a minimum threshold : " 
					+ DtoStr(minX, 2) + " ; skipping the line");
				continue;
			}
			if (x > maxX) 
			{
				if (printInfo) PrintInfo("x value is greater than a maximum threshold : " 
					+ DtoStr(maxX, 2) + " ; skipping the file");
				break;
			}
			if (y <= 0) 
			{
				if (printInfo) PrintWarning("y <= 0; skipping the line");
				continue;
			}
	
			egraphs.back().AddPoint(x, y*mult);
			egraphs.back().SetPointError(egraphs.back().GetN()-1, 0., y_stat_err/100.*y*mult);

			sgraphs.back().AddPoint(x, y*mult);
			sgraphs.back().SetPointError(sgraphs.back().GetN()-1, 0.1, y_sys_err/100.*y*mult);

			ratio_egraphs.back().AddPoint(x, y*mult);
			ratio_egraphs.back().SetPointError(ratio_egraphs.back().GetN()-1, 0., ErrPropagation(y_stat_err/100., y_sys_err/100.)*y*mult);

			ratio_sgraphs.back().AddPoint(x, y*mult);
			ratio_sgraphs.back().SetPointError(ratio_sgraphs.back().GetN()-1, 0.1, y_sys_err/100.*y*mult);
			
			vec_pt.back().push_back(x);
			vec_y.back().push_back(y);
			vec_y_stat.back().push_back(y_stat_err);
			vec_y_sys.back().push_back(y_sys_err);
			
			min_y = Minimum(min_y, y*mult);
			max_y = Maximum(max_y, y*mult);
		}

		if (egraphs.back().GetN() == 0) PrintError("Graph is empty");
		
		egraphs.back().SetMarkerStyle(marker_style);
		egraphs.back().SetMarkerColorAlpha(color, 1.);
		egraphs.back().SetMarkerSize(2.);
		egraphs.back().SetLineWidth(2);
		egraphs.back().SetLineColorAlpha(color, 0.8);

		sgraphs.back().SetMarkerStyle(0);
		sgraphs.back().SetFillColorAlpha(color, 0.4);
		sgraphs.back().SetLineColorAlpha(color, 0.8);

		ratio_egraphs.back().SetMarkerStyle(marker_style);
		ratio_egraphs.back().SetMarkerColorAlpha(color, 1.);
		ratio_egraphs.back().SetMarkerSize(2.);
		ratio_egraphs.back().SetLineWidth(2);
		ratio_egraphs.back().SetLineColorAlpha(color, 0.8);

		ratio_sgraphs.back().SetMarkerStyle(0);
		ratio_sgraphs.back().SetFillColorAlpha(0, 0);
		ratio_sgraphs.back().SetLineColorAlpha(color, 0.8);
			
		legend_entries.push_back(legend_entry);

		if (add_percent) legend_entry += "\\%";
		if (pow10 != 0) legend_entry += "#times10^{" + to_string(pow10) + "}";
		
		legend.AddEntry(egraphs.back().Clone(), legend_entry.c_str(), "P");

		levy_funcs.push_back(TF1("levy", "0.5/pi*[0]*([1] - 1)*([1] - 2)/([2] + [3]*([1]-1))/([2] + [3])*([2] + sqrt(x^2 + [3]^2)/([2] + [3]))^(-[1])*[4]"));
	
		levy_funcs.back().SetParameters(1, 1.1, 1.);
		levy_funcs.back().SetParLimits(1, 2.0001, 30.);
		levy_funcs.back().FixParameter(3, Par.particle.mass);
		levy_funcs.back().FixParameter(4, mult);
		
		exp_funcs.push_back(TF1("exp", "0.5/pi*[0]/[1]^2*exp(-x/[1])*[2]"));

		exp_funcs.back().SetParameters(0., 1.);
		exp_funcs.back().SetParameters(1., 1.);
      exp_funcs.back().FixParameter(2, mult);

		levy_funcs.back().SetLineWidth(3);
		exp_funcs.back().SetLineWidth(3);
		
		levy_funcs.back().SetLineColorAlpha(color, 0.7);
		exp_funcs.back().SetLineColorAlpha(kBlack, 0.5);

		exp_funcs.back().SetLineStyle(2);

		levy_funcs.back().SetRange(min_x*1.05, max_x/1.05);
		exp_funcs.back().SetRange(min_x*1.05, max_x/1.05);

		egraphs.back().Fit(&levy_funcs.back(), "RQMNB");
		return (TF1 *) levy_funcs.back().Clone();
	}

	TF1 *ExtendGraph(std::string spectra_file_name, const double minX, const double maxX = 1e31)
	{
		CheckInputFile(spectra_file_name);

		ifstream input_file(spectra_file_name.c_str());
		
		double x, y, y_stat_err, y_sys_err;
	
		if (printInfo) Print("From file:", spectra_file_name);

		while (input_file >> x >> y >> y_stat_err >> y_sys_err)
		{
			if (printInfo) Print(x, y, y_stat_err, y_sys_err);
			
			if (x < minX) 
			{
				if (printInfo) PrintInfo("x value is less than a minimum threshold : " 
					+ DtoStr(minX, 2) + " ; skipping the line");
				continue;
			}
			if (x > maxX) 
			{
				if (printInfo) PrintInfo("x value is greater than a maximum threshold : " 
					+ DtoStr(maxX, 2) + " ; skipping the file");
				break;
			}
			if (y <= 0) 
			{
				if (printInfo) PrintWarning("y <= 0; skipping the line");
				continue;
			}

			egraphs.back().AddPoint(x, y*mult_factor.back());
			egraphs.back().SetPointError(egraphs.back().GetN()-1, 0., y_stat_err/100.*y*mult_factor.back());

			sgraphs.back().AddPoint(x, y*mult_factor.back());
			sgraphs.back().SetPointError(sgraphs.back().GetN()-1, 0., y_sys_err/100.*y*mult_factor.back());

			ratio_egraphs.back().AddPoint(x, y*mult_factor.back());
			ratio_egraphs.back().SetPointError(ratio_egraphs.back().GetN()-1, 0., ErrPropagation(y_stat_err/100., y_sys_err/100.)*y*mult_factor.back());

			ratio_sgraphs.back().AddPoint(x, y*mult_factor.back());
			ratio_sgraphs.back().SetPointError(ratio_sgraphs.back().GetN()-1, 0., y_sys_err/100.*y*mult_factor.back());

			vec_pt.back().push_back(x);
			vec_y.back().push_back(y*mult_factor.back());
			vec_y_stat.back().push_back(y_stat_err);
			vec_y_sys.back().push_back(y_sys_err);
			
			min_y = Minimum(min_y, y*mult_factor.back());
			max_y = Maximum(max_y, y*mult_factor.back());
		}

		egraphs.back().Fit(&levy_funcs.back(), "RQMNB");
		return (TF1 *) levy_funcs.back().Clone();
	}

	void DrawDataToFitRatio(std::string output_file_name, TF1 *fit_func, std::string title = "", bool draw_sys = true)
	{
		TCanvas ratio_canv("ratio_data_to_fit", "ratio", 800, 600);

		gPad->SetLeftMargin(0.17);
		gPad->SetBottomMargin(0.11);

		TH1 *frame = gPad->DrawFrame(min_x, 0., max_x, 2.);

		frame->SetTitle(title.c_str());

		frame->GetYaxis()->SetTitleOffset(0.8);

		frame->GetXaxis()->SetLabelSize(0.05);
		frame->GetYaxis()->SetLabelSize(0.05);
		
		frame->GetXaxis()->SetTitle("p_{T}, GeV/c");
		frame->GetYaxis()->SetTitle("Data to fit ratio");

		frame->GetXaxis()->SetTitleSize(0.05);
		frame->GetYaxis()->SetTitleSize(0.05);
		
		frame->DrawClone("SAME AXIS X+ Y+");	

		ratio_line->DrawClone();

		TLegend ratio_legend(0.15, 0.8, 0.9, 0.9);
		ratio_legend.SetNColumns(5);

		ratio_legend.SetLineColorAlpha(0., 0.);
		ratio_legend.SetFillColorAlpha(0., 0.);

		for (int i = 0; i < egraphs.size(); i++) 
		{
			for (int j = 0; j < ratio_egraphs[i].GetN(); j++)
			{
				ratio_egraphs[i].SetPointY(j, ratio_egraphs[i].GetPointY(j)/
					fit_func->Eval(ratio_egraphs[i].GetPointX(j)));
				ratio_sgraphs[i].SetPointY(j, ratio_sgraphs[i].GetPointY(j)/
					fit_func->Eval(ratio_sgraphs[i].GetPointX(j)));

				ratio_egraphs[i].SetPointError(j, 0., ratio_egraphs[i].GetErrorY(j)/
					fit_func->Eval(ratio_egraphs[i].GetPointX(j)));

				ratio_sgraphs[i].SetPointError(j, 0.05, ratio_sgraphs[i].GetErrorY(j)/
					fit_func->Eval(ratio_sgraphs[i].GetPointX(j)));
			}
			
			ratio_egraphs[i].Clone()->Draw("P");
			//ratio_sgraphs[i].Clone()->Draw("5");
			ratio_legend.AddEntry(&ratio_egraphs[i], legend_entries[i].c_str(), "P");
		}
		ratio_legend.Clone()->Draw();
		
		PrintCanvas(&ratio_canv, output_file_name);
	}

	void Draw(bool draw_fits = true, bool draw_sys = false)
	{
		canv.cd();
		
		hr->SetMinimum(min_y/1e1);
		hr->SetMaximum(max_y*5.);
		
		hr->Clone()->Draw();	
		hr->Clone()->Draw("SAME AXIS X+ Y+");	

		for (int i = 0; i < egraphs.size(); i++) 
		{
			egraphs[i].Fit(&exp_funcs[i], "RQMNB");

			if (draw_fits)
			{
				exp_funcs[i].Draw("same");
				levy_funcs[i].Draw("same");
			}

			egraphs[i].Clone()->Draw("P");
			if (draw_sys) sgraphs[i].Clone()->Draw();
			
			if (printInfo) PrintInfo("Function approximation parameters for " + legend_entries[i] + ":");
			if (printInfo) Print("Levy",
				levy_funcs[i].GetParameter(0), 
				levy_funcs[i].GetParameter(1), 
				levy_funcs[i].GetParameter(2), 
				levy_funcs[i].GetParameter(3), 
				levy_funcs[i].GetParameter(4));
			
			if (printInfo) Print("expo", exp_funcs[i].GetParameter(0), exp_funcs[i].GetParameter(1));
		}
		legend.Clone()->Draw();
	}

   void Write(std::string output_file_name)
   {
		PrintCanvas(&canv, output_file_name);
   }

   void Draw(std::string output_file_name, bool draw_fits = true, bool draw_sys = false)
   {
      Draw(draw_fits, draw_sys);
      Write(output_file_name);
   }

	void PrintFitFuncsParameters(std::string output_dir)
	{
      std::string output_file_name = output_dir + "/exp_fit.tex";
      
		CheckOutputFile(output_file_name);
		ofstream output(output_file_name);
		
		std::cout.precision(6);
      
		output << "\\begin{table}[h!]" << std::endl;
		output << "\\centering" << std::endl;
		output << "\\begin{tabular}{c|c|c}" << std::endl;
		output << "Centrality class" << " & " << "dN/dy" << " & " << "T, GeV \\\\" << std::endl;
      output << "\\hline" << std::endl;
		for (int i = 0; i < exp_funcs.size(); i++)
		{
			output << " " << legend_entries[i] + "\\%" << " & " << 
				std::scientific << exp_funcs[i].GetParameter(0) << " & " << 
            exp_funcs[i].GetParameter(1) << std::defaultfloat << " \\\\ " << std::endl;
		}
		output << "\\end{tabular}" << std::endl;
		output << "\\caption{Parameters of exponential approximations of invariant $p_T$ spectra of \\Kstar meson in Au+Au at \\snn}" << std::endl;
		output << "\\end{table}" << std::endl;
		
		output.close();
		PrintInfo("File " + output_file_name + " was written");
      
      output_file_name = output_dir + "/levy_fit.tex";
      output.open(output_file_name);

		output << "\\begin{table}[h!]" << std::endl;
		output << "\\centering" << std::endl;
		output << "\\begin{tabular}{c|c|c|c|c}" << std::endl;
		output << "Centrality class" << " & " << "dN/dy" << " & " << " n & " << "$\\Lambda$, GeV/$c^2$ \\\\" << std::endl;
      output << "\\hline" << std::endl;
		for (int i = 0; i < levy_funcs.size(); i++)
		{
			output << " " << legend_entries[i] + "\\%" << " & " << 
				std::scientific << levy_funcs[i].GetParameter(0) << " & " << 
				std::scientific << levy_funcs[i].GetParameter(1) << " & " << 
				std::scientific << levy_funcs[i].GetParameter(2) << " & " << 
            std::defaultfloat << " \\\\ " << std::endl;
		}
		output << "\\end{tabular}" << std::endl;
		output << "\\caption{Parameters of levy approximations of invariant $p_T$ spectra of \\Kstar meson in Au+Au at \\snn}" << std::endl;
		output << "\\end{table}" << std::endl;
      
		output.close();
		PrintInfo("File " + output_file_name + " was written");
	}

	void PrintSpectra(std::string output_file_name)
	{
		CheckOutputFile(output_file_name);
		ofstream output(output_file_name);
		
		std::cout.precision(6);
		output << "\\begin{table}[h!]" << std::endl;
		output << "\\centering" << std::endl;
		output << "\\begin{tabular}{c|c|c|c}" << std::endl;
		output << "$p_T$, GeV/c" << " & " << "$1/(2\\pi p_T) d^2N/dp_T dy$" << " & " 
			<< "Statistical uncertainty, \\%" << " & " << "Systematic uncertainty, \\% \\\\" << std::endl;
		for (int i = 0; i < vec_pt.back().size(); i++)
		{
			output << "	" << vec_pt.back()[i] << " & " << 
				std::scientific << vec_y.back()[i] << std::defaultfloat << " & " 
				<< vec_y_stat.back()[i] << " & " << vec_y_sys.back()[i] << "\\\\" << std::endl;
		}
		output << "\\end{tabular}" << std::endl;
		output << "\\caption{Invariant $p_T$ spectra of \\Kstar meson in Au+Au at \\snn}" << std::endl;
		output << "\\end{table}" << std::endl;
		
		output.close();
		PrintInfo("File " + output_file_name + " was written");

		vec_pt.pop_back();
		vec_y.pop_back();
		vec_y_stat.pop_back();
		vec_y_sys.pop_back();
	}
};
