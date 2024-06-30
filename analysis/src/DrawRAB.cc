#include "../lib/ParInv.h"
#include "../lib/RABPainter.h"
#include "../lib/TCanvasPrinter.h"

void DrawRAB() 
{
	gROOT->SetBatch(kTRUE);
	gStyle->SetOptStat(kFALSE);

	const double ptmin = 0.9;
	const double ptmax = 6.5;

	TCanvas canv("canv", "canv", 1200, 900);
	canv.Divide(3, 2, 0., 0.);

	RABPainter rab = RABPainter("", 
		ptmin/1.1, ptmax*1.1, "p_{T}, GeV/c", "R_{AA}", false, 0.05);
	
	system(("mkdir -p ../output/Tables/RAB/" + Par.run + "/" + Par.runnum).c_str());
	
	for (int i = 0; i < Par.CType.cmin.size(); i++)
	{
		if (i >= 2) canv.cd(i+2);
		else canv.cd(i+1);
		if (i == 0 || i == 2) gPad->SetLeftMargin(0.16);
		if (i < 2) gPad->SetTopMargin(0.02);
		else gPad->SetBottomMargin(0.14);
		
		rab.AddGraph("../ext/RAB/AuAu200/phi1020_" + 
			Par.CType.cname_nop[i] + ".txt",	
			//"#varphi(1020), PPG096", kAzure+2, 54,
			"#varphi(1020)", kAzure+2, 54,
			ptmin, ptmax);

		rab.AddGraph("../ext/RAB/AuAu200/pi0_" + 
			Par.CType.cname_nop[i] + ".txt",	
			//"#pi^{0}, PPG014", kViolet-8, 55,
			"#pi^{0}", kViolet-8, 55,
			ptmin, ptmax);

		rab.AddGraph("../ext/RAB/AuAu200/proton_" + 
			Par.CType.cname_nop[i] + ".txt",	
			//"p+#bar{p}, PPG146", kSpring+2, 59,
			"p+#bar{p}", kSpring+2, 59,
			ptmin, ptmax);
		
		rab.AddGraph("../data/Spectra/" + Par.run + "/" + 
			Par.runnum + "/" + Par.particle.name_nl + "_" + 
			Par.CType.cname_nop[i] + ".txt",
			"../ext/Spectra/pp200/KStar892.txt",
			Par.particle.name, kRed+1, 53, Par.CType.ncolls[i], Par.pp_sigma,
			ptmin, ptmax);

		rab.PrintTable("../output/Tables/RAB/" + Par.run + "/" + 
			Par.runnum + "/KStar892_" + Par.CType.cname_nop[i] + ".tex");

		rab.AddScalingUncertainty(ptmax*1.05, 
			Par.CType.ncolls_uncertainty[i]/Par.CType.ncolls[i], kGray+2, 0.5, 0.1);
	
		rab.Draw(Par.CType.cname[i].c_str(), 1);
	}

	canv.cd(3);
	gPad->SetTopMargin(0.02);
	rab.Draw(Par.CType.name);

	PrintCanvas(&canv, "../output/InvM/" + Par.run + "/" + 
		Par.runnum + "/RAB_" + Par.particle.name_nl);
}	
