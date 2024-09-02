#include <array>
#include <string>
#include "ErrorHandler.h"
#include "CentralityTypes.h"

struct Part
{
	double mean, sigma;
	double ptmin, ptmax;

	Color_t color;
};

struct TOFe
{
	std::string hist_name = "m2_tofe";
	std::string name = "TOFe";
	const double ptmax_mult = 2.;
	const int rebin = 1;
	
	const double max_ms_dev = 2.;

	const double L = 5.2;

	const double sigma_alpha = 0.9;
	const double sigma_ms = 1.;
	const double sigma_t = 120;
	
	bool separate_protons = false;
	std::string pi_bg_func = "0*[0]";
	std::string k_bg_func = "0*[0]";
	std::string p_bg_func = "0*[0]";
};

struct TOFw
{
	std::string hist_name = "m2_tofw";
	std::string name = "TOFw";
	const double ptmax_mult = 1.8;
	const int rebin = 1;
	
	const double max_ms_dev = 1.8;

	const double L = 4.95;

	const double sigma_alpha = 0.9;
	const double sigma_ms = 1.;
	const double sigma_t = 84;
	
	bool separate_protons = false;
	std::string pi_bg_func = "0*[0]";
	std::string k_bg_func = "0*[0]";
	std::string p_bg_func = "0*[0]";
};

struct EMCale2
{
	std::string hist_name = "m2_emcale2";
	std::string name = "EMCale2";
	const double ptmax_mult = 1.2;
	const int rebin = 1;
	
	const double max_ms_dev = 10.;

	//distance from the vertex to a detector
	const double L = 5.22;

	const double sigma_alpha = 0.835;
	const double sigma_ms = 2.3;
	const double sigma_t = 500;
	
	bool separate_protons = false;
	std::string pi_bg_func = "0*[0]";
	std::string k_bg_func = "exp(([0]+exp([1]+[2]*x)))";
	//std::string k_bg_func = "expo(0)";
	std::string p_bg_func = "exp(pol1(0))";
};

struct EMCale3
{
	std::string hist_name = "m2_emcale3";
	std::string name = "EMCale3";
	const double ptmax_mult = 1.2;
	const int rebin = 1;
	
	const double max_ms_dev = 10.;

	//distance from the vertex to a detector
	const double L = 5.22;

	const double sigma_alpha = 0.835;
	const double sigma_ms = 2.4;
	const double sigma_t = 520;
	
	bool separate_protons = false;
	std::string pi_bg_func = "0*[0]";
	std::string k_bg_func = "exp(([0]+exp([1]+[2]*x)))";
	//std::string k_bg_func = "expo(0)";
	std::string p_bg_func = "exp(pol1(0))";
};

struct EMCalw0
{
	std::string hist_name = "m2_emcalw0";
	std::string name = "EMCalw0";
	const double ptmax_mult = 1.2;
	const int rebin = 1;
	
	const double max_ms_dev = 10.;

	//distance from the vertex to a detector
	const double L = 5.22;

	const double sigma_alpha = 0.9;
	const double sigma_ms = 1.8;
	const double sigma_t = 480;
	
	bool separate_protons = false;
	std::string pi_bg_func = "0*[0]";
	std::string k_bg_func = "exp(([0]+exp([1]+[2]*x)))";
	//std::string k_bg_func = "expo(0)";
	std::string p_bg_func = "exp(pol1(0))";
};

struct EMCalw1
{
	std::string hist_name = "m2_emcalw1";
	std::string name = "EMCalw1";
	const double ptmax_mult = 1.2;
	const int rebin = 1;
	
	const double max_ms_dev = 10.;

	//distance from the vertex to a detector
	const double L = 5.22;

	const double sigma_alpha = 0.9;
	const double sigma_ms = 1.7;
	const double sigma_t = 480;
	
	bool separate_protons = false;
	std::string pi_bg_func = "0*[0]";
	std::string k_bg_func = "exp(([0]+exp([1]+[2]*x)))";
	//std::string k_bg_func = "expo(0)";
	std::string p_bg_func = "exp(pol1(0))";
};

struct EMCalw2
{
	std::string hist_name = "m2_emcalw2";
	std::string name = "EMCalw2";
	const double ptmax_mult = 1.2;
	const int rebin = 1;
	
	const double max_ms_dev = 10.;

	//distance from the vertex to a detector
	const double L = 5.22;

	const double sigma_alpha = 0.9;
	const double sigma_ms = 1.5;
	const double sigma_t = 480;
	
	bool separate_protons = false;
	std::string pi_bg_func = "0*[0]";
	std::string k_bg_func = "exp([0]+exp([1]+[2]*x))";
	//std::string k_bg_func = "expo(0)";
	std::string p_bg_func = "exp(pol1(0))";
};

struct EMCalw3
{
	std::string hist_name = "m2_emcalw3";
	std::string name = "EMCalw3";
	const double ptmax_mult = 1.2;
	const int rebin = 1;
	
	const double max_ms_dev = 10.;

	//distance from the vertex to a detector
	const double L = 5.22;

	const double sigma_alpha = 0.9;
	const double sigma_ms = 2.;
	const double sigma_t = 500;
	
	bool separate_protons = false;
	std::string pi_bg_func = "0*[0]";
	std::string k_bg_func = "exp(([0]+exp([1]+[2]*x)))";
	//std::string k_bg_func = "expo(0)";
	std::string p_bg_func = "exp(pol1(0))";
};

struct
{
	std::vector<double> ptmin = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 
		1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 
		2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5};
	
	std::vector<double> ptmax = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 
		1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 
		2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0};

	const double write_yield_ptmin = 0.5;

	std::string mean_fit = "pol1(0)";
	std::string sigma_fit = "sqrt(pow([2]/[5], 2)*(4.*pow(pol1(0)*x, 2)) + pow([3]/[5], 2)*(4.*pow(pol1(0), 2)*(1+pol1(0)/pow(x, 2))) + pow([4]*2.9972e-4/[6], 2)*(4.*pow(x, 2)*(pol1(0) + pow(x, 2))))";
	
	Part pion, kaon, proton;
	Part apion, akaon, aproton;
	
	std::string col_sys = "AuAu200";
	std::string run_dir = "../data/";
	std::string run_name = "Run7AuAu200";
	std::string file_name = "sum";

	bool draw_log = 1;

	double land_m2min = -0.6;
	double land_m2max = 1.2;
	
	double srange = 2.;
	double draw_srange = 3.;
	double yield_extr_srange = 1.;
	double veto_srange = 3.;

	AuAu200CTypeMB4 CType;
	
	std::string magf = "";
	EMCale2 Det;
	
	//fixes sigma parameters
	bool do_fix_sigma_par = 1;
	bool do_fix_mean_par = 0;

	//field integral
	const double K1 = 104;
	
	const int iter = 20;
	
	const int cnum = 4;
	
	const double pi_m2_min = 0.015;
	const double pi_m2_max = 0.1;
	
	const double k_m2_min = 0.23;
	const double k_m2_max = 0.4;
	
	const double p_m2_min = 0.85;
	const double p_m2_max = 0.95;
	
	std::array<std::unique_ptr<TF1>, 2> def_pion_sigma = 
		{std::unique_ptr<TF1>(new TF1("pion", sigma_fit.c_str())),
		std::unique_ptr<TF1>(new TF1("apion", sigma_fit.c_str()))};
	
	std::array<std::unique_ptr<TF1>, 2> def_kaon_sigma = 
		{std::unique_ptr<TF1>(new TF1("kaon", sigma_fit.c_str())), 
		std::unique_ptr<TF1>(new TF1("kaon", sigma_fit.c_str()))};
	std::array<std::unique_ptr<TF1>, 2> def_proton_sigma = {
		std::unique_ptr<TF1>(new TF1("proton", sigma_fit.c_str())),
		std::unique_ptr<TF1>(new TF1("proton", sigma_fit.c_str()))};
} Par;

void SetPartParCommon()
{
	Par.pion.mean = pow(0.13957, 2);
	Par.kaon.mean = pow(0.493677, 2);
	Par.proton.mean = pow(0.938272, 2);

	Par.apion.mean = pow(0.13957, 2);
	Par.akaon.mean = pow(0.493677, 2);
	Par.aproton.mean = pow(0.938272, 2);
	
	Par.pion.sigma = 0.02;
	Par.kaon.sigma = 0.02;
	Par.proton.sigma = 0.1;

	Par.apion.sigma = 0.02;
	Par.akaon.sigma = 0.02;
	Par.aproton.sigma = 0.1;

	Par.pion.color = kRed;
	Par.kaon.color = kGreen;
	Par.proton.color = kBlue;

	Par.apion.color = Par.pion.color;
	Par.akaon.color = Par.kaon.color;
	Par.aproton.color = Par.proton.color;
}

void SetPartPar()
{
	SetPartParCommon();

	Par.def_pion_sigma[0]->SetParameters(Par.pion.mean, 0., 
		Par.Det.sigma_alpha, Par.Det.sigma_ms, Par.Det.sigma_t, Par.K1, Par.Det.L);
	Par.def_kaon_sigma[0]->SetParameters(Par.kaon.mean, 0., 
		Par.Det.sigma_alpha, Par.Det.sigma_ms, Par.Det.sigma_t, Par.K1, Par.Det.L);
	Par.def_proton_sigma[0]->SetParameters(Par.proton.mean, 0., 
		Par.Det.sigma_alpha, Par.Det.sigma_ms, Par.Det.sigma_t, Par.K1, Par.Det.L);

	Par.def_pion_sigma[1]->SetParameters(Par.pion.mean, 0., 
		Par.Det.sigma_alpha, Par.Det.sigma_ms, Par.Det.sigma_t, Par.K1, Par.Det.L);
	Par.def_kaon_sigma[1]->SetParameters(Par.kaon.mean, 0., 
		Par.Det.sigma_alpha, Par.Det.sigma_ms, Par.Det.sigma_t, Par.K1, Par.Det.L);
	Par.def_proton_sigma[1]->SetParameters(Par.proton.mean, 0., 
		Par.Det.sigma_alpha, Par.Det.sigma_ms, Par.Det.sigma_t, Par.K1, Par.Det.L);

	if (Par.Det.name == "TOFe")
	{
		Par.pion.ptmin = 0.4;
		Par.kaon.ptmin = 0.3;
		Par.proton.ptmin = 0.4;

		Par.pion.ptmax = 1.6;
		Par.kaon.ptmax = 1.6;
		Par.proton.ptmax = 1.6;

		Par.apion.ptmin = Par.pion.ptmin;
		Par.akaon.ptmin = Par.kaon.ptmin;
		Par.aproton.ptmin = Par.proton.ptmin;

		Par.apion.ptmax = Par.pion.ptmax;
		Par.akaon.ptmax = Par.kaon.ptmax;
		Par.aproton.ptmax = Par.proton.ptmax;
		
	}
	else if (Par.Det.name == "TOFw")
	{
		Par.pion.ptmin = 0.3;
		Par.kaon.ptmin = 0.3;
		Par.proton.ptmin = 0.4;

		Par.pion.ptmax = 1.6;
		Par.kaon.ptmax = 1.6;
		Par.proton.ptmax = 1.6;

		Par.apion.ptmin = Par.pion.ptmin;
		Par.akaon.ptmin = Par.kaon.ptmin;
		Par.aproton.ptmin = Par.proton.ptmin;

		Par.apion.ptmax = Par.pion.ptmax;
		Par.akaon.ptmax = Par.kaon.ptmax;
		Par.aproton.ptmax = Par.proton.ptmax;
	}
	else if (Par.Det.name == "EMCale2" || Par.Det.name == "EMCalw0" || Par.Det.name == "EMCalw1" || Par.Det.name == "EMCalw2" || Par.Det.name == "EMCalw3")
	{
		Par.pion.ptmin = 0.3;
		Par.kaon.ptmin = 0.4;
		Par.proton.ptmin = 0.6;

		Par.pion.ptmax = 1.2;
		Par.kaon.ptmax = 1.2;
		Par.proton.ptmax = 1.6;

		Par.apion.ptmin = 0.3;
		Par.akaon.ptmin = 0.3;
		Par.aproton.ptmin = 0.4;

		Par.apion.ptmax = Par.pion.ptmax;
		Par.akaon.ptmax = Par.kaon.ptmax;
		Par.aproton.ptmax = Par.proton.ptmax;
	}
	else if (Par.Det.name == "EMCale3")
	{
		Par.pion.ptmin = 0.3;
		Par.kaon.ptmin = 0.4;
		Par.proton.ptmin = 0.6;

		Par.pion.ptmax = 1.2;
		Par.kaon.ptmax = 1.2;
		Par.proton.ptmax = 1.6;

		Par.apion.ptmin = 0.3;
		Par.akaon.ptmin = 0.3;
		Par.aproton.ptmin = 0.6;

		Par.apion.ptmax = Par.pion.ptmax;
		Par.akaon.ptmax = Par.kaon.ptmax;
		Par.aproton.ptmax = Par.proton.ptmax;
	}
	else PrintError("There is no detector named " + Par.Det.name);
}
