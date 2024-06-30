#include "../lib/InvM.h"
#include "../lib/SpectraPainter.h"
#include "../lib/RABPainter.h"
#include "../lib/InputTool.h"

void GuiInvM()
{
	//centrality class number
	const int cnum = 0;
	
	Par.user_defined_bg_fit = true;
	Par.write_bg_par = false;
	
	Par.gui_fit_mode = true;
	
	ROOT::EnableImplicitMT();
	gStyle->SetPalette(kRainBow);

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);

	gStyle->SetTitleFontSize(0.065);
	
	gROOT->SetBatch(kTRUE);
	gROOT->ProcessLine("gErrorIgnoreLevel = 1001;");

	CheckInputFile(Par.taxi_dir + "/" + Par.run + "/" + Par.runnum + "/" + Par.file_name);

	SetEff();

	TH1::AddDirectory(false);

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

	TFile input = TFile((Par.taxi_dir + "/" + Par.run + 
		"/" + Par.runnum + "/" + Par.file_name).c_str());
	
	Box box("Events*10^-6");

	GuiFitPar.AddFitType("../par/BgFit/Run7AuAu200/" + 
		Par.runnum + "/" + Par.particle.name_nl + "/" + Par.method_name + 
		"_" + Par.CType.cname_nop[cnum] + ".txt");
	GuiFitPar.AddFitType("../par/BgFit/Run7AuAu200/" + 
		Par.runnum + "/" + Par.particle.name_nl + "/" + Par.method_name + 
		"_" + Par.CType.cname_nop[cnum] + "_free.txt");
	GuiFitPar.AddFitType("../par/BgFit/Run7AuAu200/" + 
		Par.runnum + "/" + Par.particle.name_nl + "/" + Par.method_name + 
		"_" + Par.CType.cname_nop[cnum] + "_wide.txt");
	GuiFitPar.AddFitType("../par/BgFit/Run7AuAu200/" + 
		Par.runnum + "/" + Par.particle.name_nl + "/" + Par.method_name + 
		"_" + Par.CType.cname_nop[cnum] + "_wide_free.txt");

	RBW_GDF_conv_par.sigma_conv_range = 2.;
	RBW_GDF_conv_par.number_of_iter = 10;
	
	PerformInvM(&input, cnum, 0, Par.ptmin.size(), &box);

	box.Print();

	gROOT->SetBatch(kFALSE);
	
	TCanvas *gui_canv = new TCanvas("GuiFit canv", "GuiFit", 1000., 1000.);
	gui_canv->cd();
	
	gPad->SetBottomMargin(0.13);
	gPad->SetLeftMargin(0.13);
		
	GuiFitPar.Draw();

	gPad->Modified();
	gPad->Update();

	gPad->AddExec("GuiFit", "GuiFit()");

	Print();
}
