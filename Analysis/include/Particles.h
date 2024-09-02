#pragma once

struct KStar892
{
	std::string name = "K*(892)";
	std::string name_nl = "KStar892";
	
	const double mass = 0.892;
	const double gamma = 50.8e-3;

	const double zoom_min = 0.68;
	const double zoom_max = 1.1;
	
	//noPID
	//const double int_low = 1.8;
	//2PID
	const double int_low = 1.1;
	
	const double int_upp = 2.6;

	bool has_antipart = true;

	std::string dname1 = "pi";
	std::string dname2 = "K";

	const double BR = 0.66601;
	const double BR_uncertainty = 9e-5;
};

struct phi1020
{
	std::string name = "#varphi(1020)";
	std::string name_nl = "phi1020";
	
	const double mass = 1.01945;
	const double gamma = 6.2e-3;

	const double zoom_min = 1.0;
	const double zoom_max = 1.05;
	
	const double int_low = 1.06;
	const double int_upp = 1.5;

	const double fit_min_range = 1.0;
	const double fit_max_range = 1.04;

	bool has_antipart = 0;

	std::string dname1 = "K";
	std::string dname2 = "K";

	const double BR = 0.5;
	const double BR_uncertainty=5e-3;
};

struct Lambda1520
{
	std::string name = "#Lambda(1520)";
	std::string name_nl = "Lambda1520";
	
	const double mass = 1.52;
	const double sigma = 0.0156;

	const double zoom_min = 1.46;
	const double zoom_max = 1.6;
	
	const double int_low = 1.65;
	const double int_upp = 2.1;

	const double fit_min_range = 1.48;
	const double fit_max_range = 1.6;

	bool has_antipart = 0;

	std::string dname1 = "P";
	std::string dname2 = "K";

	const double BR = 0.5;
};

struct D0
{
	std::string name = "D0";
	std::string name_nl = "D0";
	
	const double mass = 1.864;
	const double sigma = 30;

	const double zoom_min = 1.6;
	const double zoom_max = 2.;
	
	const double int_low = 1.8;
	const double int_upp = 2.;

	const double fit_min_range = 1.8;
	const double fit_max_range = 1.92;

	bool has_antipart = true;

	std::string dname1 = "pi";
	std::string dname2 = "K";

	const double BR = 0.0397;
};

struct JS1S
{
	std::string name = "#J/#psi(1S)";
	std::string name_nl = "JS1S";
	
	const double mass = 3.0969;
	const double sigma = 0.0929;

	const double zoom_min = 2.85;
	const double zoom_max = 3.25;
	
	const double int_low = 1.5;
	const double int_upp = 2.5;

	const double fit_min_range = 2.9;
	const double fit_max_range = 3.2;

	bool has_antipart = 0;

	std::string dname1 = "e";
	std::string dname2 = "e";
};

struct Lambda1670
{
	std::string name = "#Lambda(1670)";
	std::string name_nl = "Lambda1670";
	
	const double mass = 1.67;
	const double sigma = 0.03;

	const double zoom_min = 1.6;
	const double zoom_max = 1.8;
	
	const double int_low = 1.9;
	const double int_upp = 2.1;

	const double fit_min_range = 1.61;
	const double fit_max_range = 1.76;

	bool has_antipart = 1;

	std::string dname1 = "P";
	std::string dname2 = "K";
};

struct B0
{
	std::string name = "B0(5279)";
	std::string name_nl = "B05279";
	
	const double mass = 5.3;
	const double sigma = 8e-3;

	const double zoom_min = 5.2;
	const double zoom_max = 5.5;
	
	const double int_low = 1.1;
	const double int_upp = 1.5;

	const double fit_min_range = 5.25;
	const double fit_max_range = 5.35;

	bool has_antipart = 0;
};


