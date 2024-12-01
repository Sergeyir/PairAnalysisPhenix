#pragma once

#include <thread>

#include "Particles.h"

#include "../../analysis/lib/DeadAreasCuts.h"

using namespace CutsRun7AuAu200MinBias;

struct
{
	const std::string system = "AuAu200";
	const std::string run = "Run7";

	std::string part = "KStar892";
	std::vector<std::string> magf_queue = {"+-", "-+"};

	std::vector<std::string> daughter1_queue = {"pion", "kaon"};
	std::vector<std::string> daughter2_queue = {"akaon", "apion"};

	bool do_use_weight_func = true;
	
	const double ptmin = 0.3;
	const double ptmax = 8.;

	const int nthreads = std::thread::hardware_concurrency();

	const int pt_nbins = static_cast<int>((ptmax - ptmin)*10.);

	const std::string run_name = run + system;
} Par;
