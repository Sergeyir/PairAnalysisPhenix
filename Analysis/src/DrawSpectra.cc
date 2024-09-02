#include "../lib/ParInv.h"
#include "../lib/SpectraPainter.h"
#include "../lib/LogoDrawer.h"

struct
{
	std::array<std::string, 5> methods = {"EMC2PID", "2PID", "TOF2PID", "1PID", "noPID"};
	std::array<double, 5> methods_y = {};

	std::array<Color_t, 5> methods_color = {kAzure-3, kRed-3, kGreen-3, kViolet-3, kOrange-3};
	std::array<int, 5> methods_marker = {59, 55, 57, 67, 53};

	std::array<std::array<double, 5>, Par.CType.size> methods_ptmin = {{
		{0.9, 0.9, 1.9, 1.7, 2.3},
		{0.9, 0.9, 1.9, 1.7, 2.3},
		{0.9, 0.9, 1.9, 1.7, 2.3},
		{0.9, 0.9, 1.9, 1.7, 2.3},
		{0.9, 0.9, 1.9, 1.7, 2.3}}};

	std::array<std::array<double, 5>, Par.CType.size> methods_ptmax = {{
		{2.3, 4.0, 3.4, 4.5, 6.5},
		{2.3, 4.0, 3.4, 4.5, 6.5},
		{2.3, 4.0, 3.4, 4.5, 6.5},
		{2.3, 4.0, 3.4, 4.5, 6.5},
		{2.3, 4.0, 3.4, 4.5, 6.5}}};
	
   /*
	std::array<std::array<double, 5>, Par.CType.size> merge_methods_ptmin = {{
		{-9., 0.9, -9., 2.9, 4.0},
		{-9., 0.9, -9., 2.6, 3.4},
		{-9., 0.9, -9., 2.6, 4.0},
		{-9., 0.9, -9., 1.9, 3.4},
		{-9., 0.9, -9., 2.9, 4.0}}};
	
	std::array<std::array<double, 5>, Par.CType.size> merge_methods_ptmax = {{
		{-9., 2.9, -9., 4.0, 6.5},
		{-9., 2.6, -9., 3.4, 6.5},
		{-9., 2.6, -9., 4.0, 6.5},
		{-9., 1.9, -9., 3.4, 6.5},
		{-9., 2.9, -9., 4.0, 6.5}}};
   */

	std::array<std::array<double, 5>, Par.CType.size> merge_methods_ptmin = {{
		{-9., 0.9, -9., 2.6, 4.0},
		{-9., 0.9, -9., 2.6, 4.0},
		{-9., 0.9, -9., 2.6, 4.0},
		{-9., 0.9, -9., 2.6, 4.0},
		{-9., 0.9, -9., 2.6, 4.0}}};
	
	std::array<std::array<double, 5>, Par.CType.size> merge_methods_ptmax = {{
		{-9., 2.6, -9., 4.0, 6.5},
		{-9., 2.6, -9., 4.0, 6.5},
		{-9., 2.6, -9., 4.0, 6.5},
		{-9., 2.6, -9., 4.0, 6.5},
		{-9., 2.6, -9., 4.0, 6.5}}};
	
	//for bin shift correction
	std::array<std::vector<double>, Par.CType.size> pt, ptmin, ptmax, y, y_stat, y_sys;

	void ApplyBinShiftCorr(const int cnum)
	{
		//tsallis
		/*
		TF1 fit_func("func", "[0]*x^([5]-1)*(1 + ([1] - 1)*(sqrt([4] + x^2)^([5]) - [2])/[3])^(-1./([1]-1.))");
		fit_func.SetParameters(10., 1.01, 0.2, 1000.); //scale, q, mu, and T

		fit_func.SetParLimits(1, 1.00001, 2.);
		fit_func.FixParameter(4, pow(Par.particle.mass, 2)); //m2
		fit_func.FixParameter(5, 1.); //kappa
		*/
		
		//Levy
		TF1 fit_func("func", "0.5/pi*[0]*([1] - 1)*([1] - 2)/([2] + [3]*([1]-1))/([2] + [3])*([2] + sqrt(x^2 + [3]^2)/([2] + [3]))^(-[1])");

		fit_func.SetParameters(1, 1.1, 1.);
		fit_func.SetParLimits(0, 0., 10.); //dN/dy
		fit_func.SetParLimits(1, 2.0001, 30.); //n
		fit_func.SetParLimits(2, 1., 10.); //Lambda
		fit_func.FixParameter(3, Par.particle.mass); //m
		
		for (int i = 0; i < 5; i++)
		{
			TGraphErrors gr;
			for (int j = 0; j < pt[cnum].size(); j++)
			{
				gr.AddPoint(pt[cnum][j], y[cnum][j]);
				gr.SetPointError(gr.GetN()-1, 0., 
					y_stat[cnum][j]/100.*y[cnum][j]);
			}
			gr.Fit(&fit_func, "QB");
			for (int j = 0; j < pt[cnum].size(); j++)
			{
				const double r = fit_func.Integral(
					ptmin[cnum][j], ptmax[cnum][j])/
					fit_func.Eval(pt[cnum][j])/
					(ptmax[cnum][j] - ptmin[cnum][j]);
				
				y[cnum][j] = y[cnum][j]/r;
			}
			PrintInfo("Hagedorn parameters for " + Par.CType.cname_nop[cnum] + ":");
			Print("	", fit_func.GetParameter(0), fit_func.GetParameter(1), fit_func.GetParameter(2));
		}
	}
} Vec;

void DrawSpectra() 
{
	gROOT->SetBatch(kTRUE);
	gStyle->SetOptStat(kFALSE);

	for (int i = 0; i < Par.CType.size; i++)
	{
		for (int j = 0; j < Vec.methods.size(); j++)
		{
			std::string input_file_name = "../data/Spectra/" + Par.run + "/" + 
				Par.runnum + "/" + Vec.methods[j] + "_" + Par.particle.name_nl + "_" + 
				Par.CType.cname_nop[i] + ".txt";
			CheckInputFile(input_file_name);
			
			ifstream input_file(input_file_name);
			
			double x, y, stat, sys;
			while (input_file >> x >> y >> stat >> sys)
			{
				if (x > Vec.merge_methods_ptmin[i][j] && x < Vec.merge_methods_ptmax[i][j])
				{
					Vec.pt[i].push_back(x);
					Vec.y[i].push_back(y);
					Vec.y_stat[i].push_back(stat);
					Vec.y_sys[i].push_back(sys);
				}
			}
		}

		Vec.ptmin[i].resize(Vec.pt[i].size());
		Vec.ptmax[i].resize(Vec.pt[i].size());

		for (int j = 0; j < Vec.pt[i].size(); j++)
		{
			for (int k = 0; k < DefaultPt.ptmin.size(); k++)
			{
				if (Vec.pt[i][j] > DefaultPt.ptmin[k] && Vec.pt[i][j] < DefaultPt.ptmax[k])
				{
					Vec.ptmin[i][j] = DefaultPt.ptmin[k];
					Vec.ptmax[i][j] = DefaultPt.ptmax[k];
					break;
				}
			}
		}

		Vec.ApplyBinShiftCorr(i);

		std::string spectra_file_name = "../data/Spectra/" + Par.run + "/" + 
			Par.runnum + "/" + Par.particle.name_nl + "_" + Par.CType.cname_nop[i] + ".txt";

		ofstream spectra_file(spectra_file_name);
		for (int j = 0; j < Vec.pt[i].size(); j++)
		{
			spectra_file << Vec.pt[i][j] << " " << Vec.y[i][j] << 
				" " << Vec.y_stat[i][j] << " " << Vec.y_sys[i][j] << std::endl;
		}
	}

	std::array<TF1 *, Par.CType.size> levy_fit;

   //signle collision system spectra only
	SpectraPainter single_spectra = SpectraPainter(Par.particle.name + " " + Par.CType.name, 
		0.7, 6.1, "p_{T} [GeV/c]", 
		"1/(2#pip_{T}) d^{2} N/dp_{T}/dy [(GeV/c)^{-2}]");

	const int pow10 = Par.CType.size-1;
	std::string legend_entry = Par.CType.cname_nop[Par.CType.size-1];
	
	std::string input_file = "../data/Spectra/" + Par.run + "/" + 
		Par.runnum + "/" + Par.particle.name_nl + "_" + 
		Par.CType.cname_nop[Par.CType.size-1] + ".txt";

	system(("mkdir -p ../output/Tables/Spectra/" + Par.run + "/" + Par.runnum).c_str());

	levy_fit[Par.CType.size-1] = single_spectra.AddGraph(
		input_file, legend_entry, pow10, 1., Par.particle.mass,
		Par.CType.marker_style[Par.CType.size-1], Par.CType.color[Par.CType.size-1],
		Vec.ptmin.back().front(), Vec.ptmax.back().back());
	
	single_spectra.PrintSpectra("../output/Tables/Spectra/" + 
		Par.run + "/" + Par.runnum + "/" + Par.particle.name_nl + "_" + 
		Par.CType.cname_nop.back() + ".tex");
	
	for (int i = 0; i < Par.CType.size-1; i++)
	{
		const int pow10 = Par.CType.size-i-2;
		std::string legend_entry = Par.CType.cname_nop[i];
		
		std::string input_file = "../data/Spectra/" + Par.run + "/" + 
			Par.runnum + "/" + Par.particle.name_nl + "_" + 
			Par.CType.cname_nop[i] + ".txt";
		
		levy_fit[i] = single_spectra.AddGraph(input_file,
			legend_entry, pow10, 1., Par.particle.mass, Par.CType.marker_style[i], Par.CType.color[i],
			Vec.ptmin.back().front(), Vec.ptmax.back().back());

		system(("mkdir -p ../output/Tables/Spectra/" + Par.run + "/" + Par.runnum).c_str());
		
		single_spectra.PrintSpectra("../output/Tables/Spectra/" + 
			Par.run + "/" + Par.runnum + "/" + Par.particle.name_nl + "_" + 
			Par.CType.cname_nop[i] + ".tex");
	}
   	
   single_spectra.PrintFitFuncsParameters("../output/Tables/Spectra/" + Par.run + "/" + Par.runnum);

	single_spectra.Draw();

   DrawPHENIXLogoPreliminary(0.6, 0.7, 1.8);
	
   single_spectra.Write("../output/InvM/" + Par.run + "/" + 
		Par.runnum + "/SingleSpectra_" + Par.particle.name_nl);
	
   //all spectra
	SpectraPainter spectra = SpectraPainter(Par.particle.name + " " + Par.CType.name, 
		0.7, 6.1, "p_{T} [GeV/c]", 
		"1/(2#pip_{T}) d^{2} N/dp_{T}/dy [(GeV/c)^{-2}]");

	legend_entry = Par.CType.cname_nop[Par.CType.size-1];
	
	input_file = "../data/Spectra/" + Par.run + "/" + 
		Par.runnum + "/" + Par.particle.name_nl + "_" + 
		Par.CType.cname_nop[Par.CType.size-1] + ".txt";

	levy_fit[Par.CType.size-1] = spectra.AddGraph(
		input_file, legend_entry, pow10, 1., Par.particle.mass,
		Par.CType.marker_style[Par.CType.size-1], Par.CType.color[Par.CType.size-1],
		Vec.ptmin.back().front(), Vec.ptmax.back().back());
	
	for (int i = 0; i < Par.CType.size-1; i++)
	{
		const int pow10 = Par.CType.size-i-2;
		std::string legend_entry = Par.CType.cname_nop[i];
		
		std::string input_file = "../data/Spectra/" + Par.run + "/" + 
			Par.runnum + "/" + Par.particle.name_nl + "_" + 
			Par.CType.cname_nop[i] + ".txt";
		
		levy_fit[i] = spectra.AddGraph(input_file,
			legend_entry, pow10, 1., Par.particle.mass, Par.CType.marker_style[i], Par.CType.color[i],
			Vec.ptmin.back().front(), Vec.ptmax.back().back());
	}
   	
	spectra.AddGraph("../ext/Spectra/pp200/KStar892.txt",
		//"p+p (AN911)", 0, 1./Par.pp_sigma, Par.particle.mass, 57, kGray+2);
		"p+p, PRC90, 054905", 0, 1./Par.pp_sigma, Par.particle.mass, 57, kGray+2);

	spectra.Draw();

   DrawPHENIXLogoPreliminary(0.6, 0.7, 1.8);
	
   spectra.Write("../output/InvM/" + Par.run + "/" + 
		Par.runnum + "/Spectra_" + Par.particle.name_nl);
   
	for (int i = 0; i < Par.CType.size; i++)
	{
		levy_fit[i]->SetParameter(4, 1.);
		SpectraPainter spectra_comp = SpectraPainter(Par.particle.name + " " + Par.CType.name, 
			0.7, 6.1, "p_{T} [GeV/c]", 
			"1/(2 #pi) 1/p_{T} d^{2} N/dp_{T}/dy [(GeV/c)^{-2}]", 400, 500);

		for (int j = Vec.methods.size()-1; j >= 0; j--)
		{
			spectra_comp.AddGraph("../data/Spectra/" + Par.run + "/" + 
				Par.runnum + "/" + Vec.methods[j] + "_" + Par.particle.name_nl + "_" + 
				Par.CType.cname_nop[i] + ".txt",
				Vec.methods[j], 0, 1., Par.particle.mass, 
				Vec.methods_marker[j], Vec.methods_color[j],
				0.9, 6.5, false);
		}

		spectra_comp.Draw("../output/InvM/" + Par.run + "/" + 
			Par.runnum + "/Spectra_" + Par.particle.name_nl + 
			"_" + Par.CType.cname_nop[i] + "_comp", 0);
		
		spectra_comp.DrawDataToFitRatio("../output/InvM/" + Par.run + "/" + 
			Par.runnum + "/Spectra_" + Par.particle.name_nl + 
			"_" + Par.CType.cname_nop[i] + "_ratio", (TF1 *) levy_fit[i]->Clone());
	}
}	
