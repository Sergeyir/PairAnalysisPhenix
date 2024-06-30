#include "../lib/InvM.h"
#include "../lib/SpectraPainter.h"
#include "../lib/RABPainter.h"
#include "../lib/InputTool.h"

void InvM(bool do_bg_subtr = true)
{
	ROOT::EnableImplicitMT();
	
	gStyle->SetPalette(kRainBow);

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);

	gStyle->SetTitleFontSize(0.065);
	
	gROOT->SetBatch(kTRUE);
	gROOT->ProcessLine("gErrorIgnoreLevel = 1001;");

	CheckInputFile(Par.taxi_dir + "/" + Par.run + "/" + Par.runnum + "/" + Par.file_name);

	SetEff();

	if (Par.use_gdf)
	{
		std::string gdf_par_file_name = "../par/GDF/" + Par.run + "/" + Par.particle.name_nl + ".txt";
		CheckInputFile(gdf_par_file_name);
		double *gdf_par = ReadFileIntoArray(gdf_par_file_name, 2);
		Par.gdf_sigma.SetParameters(gdf_par);
	}
	else
	{
		Par.gdf_sigma.SetParameters(0., 0.);
	}

	SimpleCheck();
	
	system(("mkdir -p ../output/InvM/" + Par.run + "/" + Par.runnum).c_str());
	system(("mkdir -p ../data/Spectra/" + Par.run + "/" + Par.runnum).c_str());
	system(("mkdir -p ../data/RawYields/" + Par.run + "/" + Par.runnum).c_str());

	TFile input = TFile((Par.taxi_dir + "/" + Par.run + 
		"/" + Par.runnum + "/" + Par.file_name).c_str());
	
	PrintParameters();
	
	PrintBigSeparator("Invariant mass subtracting and fitting");

	Box box("Events*10^-6");
	
	//for one centrality class only
	//PerformInvM(&input, 1, 0, Par.ptmin.size(), &box);
	
	//for all centrality classes in CType
	for (int i = 0; i < Par.CType.cmin.size(); i++)
	{
		PerformInvM(&input, i, 0, Par.ptmin.size(), &box);
	}

	box.Print();
	Print();

	if (Par.use_eff)
	{
		if (Par.do_draw_spectra) system("root -q -b -l DrawSpectra.cc");
		if (Par.do_draw_rcp) system("root -q -b -l DrawRCP.cc");
		if (Par.do_draw_rab) system("root -q -b -l DrawRAB.cc");
	}
}
