#pragma once

#include <thread>

#include "TFile.h"
#include "TH1.h"

#include "Particles.h"

#include "../../analysis/lib/CentralityTypes.h"
#include "../../analysis/lib/DeadAreasCuts.h"

using namespace CutsRun7AuAu200MinBias;

struct
{
	const std::string system = "AuAu200";
	const std::string run = "Run7";
	const bool do_use_weight_func = true;
	bool do_reweight_alpha = true;

	std::vector<std::string> part_queue = {"pion", "kaon", "apion", "akaon", "proton", "aproton"};
	std::vector<std::string> magf_queue = {"+-", "-+"};
	std::vector<std::string> aux_name_queue = {"_lpt", "_hpt"};

	const double ptmin = 0.3;
	const double ptmax = 8.;

	const int pt_nbins = static_cast<int>((ptmax - ptmin)*10.);

	const int nthreads = std::thread::hardware_concurrency();

	const std::string run_name = run + system;
	
	std::unique_ptr<TFile> data_input_file;
	std::unique_ptr<TFile> sim_input_file;
	
	TH1F *dce0_alpha_reweight;
	TH1F *dce1_alpha_reweight;
	TH1F *dcw0_alpha_reweight;
	TH1F *dcw1_alpha_reweight;
} Par;
