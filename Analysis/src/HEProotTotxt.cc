#include "../lib/OutputTool.h"
#include "../lib/ErrorHandler.h"

struct
{
	std::string input_dir = "../data/HEPData/";
	std::string output_dir = "../ext/Spectra/";
	std::string system = "AuAu200";
	std::string file_name = "phi1020";
	std::string addit_file_name = "_60-93";
	std::string tfile_dir = "Figure 3";
	std::string hist_name = "Hist1D_y8";
} Par;

void HEProotTotxt()
{
	std::string input_file_name = Par.input_dir + Par.file_name + ".root";
	std::string output_file_name = Par.output_dir + Par.system + "/" + Par.file_name + Par.addit_file_name + ".txt";

	CheckInputFile(input_file_name);

	TFile input = TFile(input_file_name.c_str());
	TH1F *data_points = (TH1F *) input.Get((Par.tfile_dir + "/" + Par.hist_name).c_str());
	TH1F *stat_err_points = (TH1F *) input.Get((Par.tfile_dir + "/" + Par.hist_name + "_e1").c_str());
	TH1F *syst_err_points = (TH1F *) input.Get((Par.tfile_dir + "/" + Par.hist_name + "_e2").c_str());

	Print(Par.tfile_dir + Par.hist_name);

	PrintInfo("File " + input_file_name + " was opened for reading");
	
	ofstream output(output_file_name);

	for (int i = 1; i <= data_points->GetXaxis()->GetNbins(); i++)
	{
		double x = data_points->GetXaxis()->GetBinCenter(i);
		double y = data_points->GetBinContent(i);
		double y_stat = stat_err_points->GetBinContent(i);
		double y_sys = syst_err_points->GetBinContent(i);

		output << x << "	" << y << "	" << y_stat/y*100. << "	" << y_sys/y*100. << std::endl;
	}

	output.close();
	PrintInfo("File " + output_file_name + " was written");
}
