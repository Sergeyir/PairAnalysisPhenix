#include "../lib/ErrorHandler.h"
#include "../lib/OutputTool.h"
#include "../lib/Tool.h"

using namespace Tool;

struct
{
	std::string run_name = "Run7AuAu200";
	std::string particle_name = "kaon";
	std::unique_ptr<TFile> input_file;
	
	std::vector<std::string> centr_name = {"00-20", "20-40", "40-60", "60-93", "00-93"};
	std::vector<std::string> det_name = {"dc_pc1", "pc3", "pc2", "tofe", "tofw", "emcale0", "emcale1", "emcale2", "emcale3", "emcalw0", "emcalw1", "emcalw2", "emcalw3"};
} Par;

void PrintEmb(std::string det_name)
{	
	TH1F *sim_hist = (TH1F *) Par.input_file->Get(("reg_" + det_name + "_sim").c_str());
	TH1F *real_hist = (TH1F *) Par.input_file->Get(("reg_" + det_name + "_real").c_str());
		
	std::vector<double> emb_err;
	if (det_name != "dc_pc1") std::cout << "DC-PC1 + ";
	else std::cout << "DC-PC1 ";

	const std::string output_file_name = "../../sim_analysis/input/Embedding/" + 
		Par.run_name + "/" + Par.particle_name + "_" + det_name + ".txt";
	ofstream output_file(output_file_name);
	
	for (int i = 1; i <= Par.centr_name.size(); i++)
	{
		const double emb_val = real_hist->GetBinContent(i)/sim_hist->GetBinContent(i);
		emb_err.push_back(emb_val*ErrPropagation(
			real_hist->GetBinError(i)/real_hist->GetBinContent(i), 
			sim_hist->GetBinError(i)/sim_hist->GetBinContent(i)));
		if (i != Par.det_name.size() - 1) std::cout << "& ";
		std::cout << emb_val;
		if (i != Par.det_name.size() - 1) std::cout << " $\\pm$ ";

		output_file << emb_val << " ";
	}
	std::cout << "\\\\" << std::endl;

	output_file.close();
	
	if (det_name == "dc_pc1") std::cout << " ";
	else if (det_name == "pc3") std::cout << "PC3 ";
	else if (det_name == "pc2") std::cout << "PC2 ";
	else if (det_name == "tofe") std::cout << "TOFe ";
	else if (det_name == "tofw") std::cout << "TOFw ";
	else if (det_name == "emcale0") std::cout << "EMCale0 ";
	else if (det_name == "emcale1") std::cout << "EMCale1 ";
	else if (det_name == "emcale2") std::cout << "EMCale2 ";
	else if (det_name == "emcale3") std::cout << "EMCale3 ";
	else if (det_name == "emcalw0") std::cout << "EMCalw0 ";
	else if (det_name == "emcalw1") std::cout << "EMCalw1 ";
	else if (det_name == "emcalw2") std::cout << "EMCalw2 ";
	else if (det_name == "emcalw3") std::cout << "EMCalw3 ";
	else PrintError("Unknown detector name: " + det_name);
	
	for (int i = 0; i < Par.centr_name.size(); i++)
	{
		if (i != Par.det_name.size() - 1) std::cout << "& ";
		std::cout << emb_err[i];
		if (i != Par.det_name.size() - 1) std::cout << " ";
	}
	std::cout << "\\\\" << std::endl;
	std::cout << "\\hline" << std::endl;
}

void EmbEff()
{
	std::string input_file_name = "../data/phenix_sim/" + 
		Par.run_name + "/Emb_" + Par.particle_name + ".root";
	CheckInputFile(input_file_name);
	Par.input_file = std::unique_ptr<TFile>(new TFile(input_file_name.c_str()));

	std::cout << "	 ";
	for (std::string centr_name : Par.centr_name)
	{
		std::cout << centr_name << " ";
	}
	std::cout << std::endl;

	system(("mkdir -p ../../sim_analysis/input/Embedding/" + Par.run_name).c_str());
	PrintInfo("Output files will be saved in ../../sim_analysis/input/Embedding/" + Par.run_name);

	for (std::string det_name : Par.det_name)
	{
		PrintEmb(det_name);
	}
}
