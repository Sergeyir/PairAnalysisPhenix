#pragma once

#include "CentralityTypes.h"
#include "Particles.h"
#include "Table.h"

const double pi = 3.14159265359;

Table<std::string, std::string, std::string, std::string, std::string> table;

struct
{
	//K*(892)
	//Run7AuAu
	std::vector<double> ptmin = {0.9, 1.1, 1.4, 1.7, 1.9, 2.1, 2.3, 2.6, 2.9, 3.4, 4.0, 4.5, 5.0};
	std::vector<double> ptmax = {1.1, 1.4, 1.7, 1.9, 2.1, 2.3, 2.6, 2.9, 3.4, 4.0, 4.5, 5.0, 6.5};
	//full
	//std::vector<double> ptmin = {0.9, 1.1, 1.4, 1.7, 1.9, 2.1, 2.3, 2.6, 2.9, 3.4, 4.0, 4.5, 5.0, 6.5, 8.5};
	//std::vector<double> ptmax = {1.1, 1.4, 1.7, 1.9, 2.1, 2.3, 2.6, 2.9, 3.4, 4.0, 4.5, 5.0, 6.5, 8.5, 12.}; Lambda(1520)
	//std::vector<double> ptmin = {0.6, 0.9, 1.4, 1.9, 2.3, 2.6, 3.4, 4.5, 6.5};
	//std::vector<double> ptmax = {0.9, 1.4, 1.9, 2.3, 2.6, 3.4, 4.5, 6.5, 12.};
} DefaultPt;

struct
{
	//std::string taxi_dir = "/home/sergey/mount/taxi";
	std::string taxi_dir = "../data/taxi";
	std::string run = "Run7AuAu200";
	std::string runnum = "19079";
	std::string file_name = "sum.root";
	
	KStar892 particle;
	std::string channel = particle.dname1 + particle.dname2;
	
	AuAu200CTypeMB4 CType;

	//const double int_low = particle.int_low;
	
	std::string method_name = "TOF2PID"; const double fg_corr_part_integral_threshold = 0.8;
	//std::string method_name = "EMC2PID"; const double fg_corr_part_integral_threshold = 0.7;
	//std::string method_name = "2PID"; const double fg_corr_part_integral_threshold = 0.8;
	//std::string method_name = "1PID"; const double fg_corr_part_integral_threshold = 0.98;
	//std::string method_name = "noPID"; const double fg_corr_part_integral_threshold = 0.99;
	
	std::string addit_name = "";
	
	TF1 gdf_sigma = TF1("gdf_sigma", "pol1");
	
	std::vector<double> ptmin, ptmax;
	std::vector<std::vector<double>> eff, eff_err, eff_pt_sc_sys_err, eff_acc_sys_err, eff_m2_eff_sys_err;
	
	const double fg_norm_part_integral_threshold = 0.1;	
	
	bool draw_fit = 1;
	
	//if true pt range of efficiency will be applied
	//if false spectra, rcp, and rab will not be drawn, pt is default
	bool use_eff = 1;
	
	bool do_draw_spectra = 1;
	bool do_draw_rcp = 1;
	bool do_draw_rab = 1;
	
	bool user_defined_bg_fit = 1;
	//write_bg_par is ignored and bg parameters are not written if use_defined_bg_fit = 1
	bool write_bg_par = 0;

	bool write_yields = true;
	
	bool draw_sys = 1;
	
	bool use_gdf = 1;
	bool gui_fit_mode = 0;
	
	const double ptmax_nomb = 6.5;
	
	const int rebin_num = 5;
	const double srange = 2.;
	const double wide_srange = 2.5;
	
	const int zmin = 0;
	const int zmax = 2;

	const int rmin = 0;
	const int rmax = 0;
	
	//number of digits after the point
	int nf = 1;
	
	//sets limits on sigma normalized in sigma in breit-wigner fit
	const double min_ssigma = 0.93;
	const double max_ssigma = 1.03;
	
	const double pp_sigma = 42.2;
	
	std::string pp_spectra = particle.name_nl;
} Par;
