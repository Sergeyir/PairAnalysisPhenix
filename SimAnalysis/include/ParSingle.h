#include <thread>

#include "TFile.h"

#include "Particles.h"
#include "../../analysis/lib/DeadAreasCuts.h"

using namespace CutsRun7AuAu200MinBias;

struct
{	
	const std::string system = "AuAu200";
	const std::string run = "Run7";
	bool do_use_weight_func = true;
	bool do_reweight_alpha = false;
	
	std::vector<std::string> magf_queue = {"+-", "-+"};
	std::vector<std::string> part_queue = {"pion", "apion", "kaon", "akaon"};
	std::vector<std::string> aux_name_queue = {"_lpt", "_hpt"};
	
	const double ptmin = 0.3;
	const double ptmax = 8.;
	
	const int nthreads = std::thread::hardware_concurrency();

	const int pt_nbins = static_cast<int>((ptmax - ptmin)*10.);
	
	const std::string run_name = run + system;

	std::unique_ptr<TFile> data_input_file;
   std::unique_ptr<TFile> sim_input_file;
} Par;
