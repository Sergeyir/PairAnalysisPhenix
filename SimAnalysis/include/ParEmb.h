#include <thread>

#include "Particles.h"
#include "../../analysis/lib/DeadAreasCuts.h"

using namespace CutsRun7AuAu200MinBias;

struct
{	
	const std::string system = "AuAu200";
	const std::string run = "Run7";

	std::vector<std::string> magf_queue = {"+-", "-+"};
	std::vector<std::string> part_queue = {"pion", "kaon"};
	std::vector<std::string> centr_queue = {"00-20", "20-40", "40-60", "60-93", "00-93"};
	std::vector<std::string> part_charge_queue = {"", "a"};

	const double ptmin = 0.3;
	const double ptmax = 8.;

	const int nthreads = std::thread::hardware_concurrency();

	const unsigned int centr_nbins = centr_queue.size();

	const std::string run_name = run + system;
} Par;
