#pragma once

#include <string>
#include <map>
#include <vector>
#include <array>

struct 
{
	const int pion = 211;
	const int kaon = 321;
	const int proton = 2212;
	const int electron = 11;
	const int noPID = 1999;
	const int junk = -9999;
} PartId;

struct ParticleStruct
{
	std::string name;
	std::string sname;
	float mass;
	int geant_id;
	int id;
	int charge;
};

struct
{
	std::map<std::string, int> iter_map = 
	{
		{"electron", 0},
		{"positron", 1},
		{"pion", 2},
		{"apion", 3},
		{"kaon", 4},
		{"akaon", 5},
		{"proton", 6},
		{"aproton", 7}
	};

	std::vector<int> id = {11, -11, 211, -211, 321, -321, 2212, -2212};

	std::vector<std::string> sname = {"e", "e", "pi", "pi", "K", "K", "P", "P"};
	std::vector<std::string> name = {"electron", "positron", 
		"pion", "apion", "kaon", "akaon", "proton", "aproton"};
	
	//for embedding
	std::vector<std::string> abs_name = {"electron", "electron", 
		"pion", "pion", "kaon", "kaon", "proton", "proton"};
	
	std::vector<double> mass = {0.5110999e-3, 0.5110999e-3, 
		0.139570, 0.139570, 0.493677, 0.493677, 0.938272, 0.938272};
	std::vector<int> geant_id = {1, 2, 8, 9, 11, 12, 14, 15};
	std::vector<int> charge = {-1, 1, 1, -1, 1, -1, 1, -1};
} ParticleProperties;

struct Lambda1520
{
	const std::string name = "lambda1520";
	const float mass = 1.5195;
	const float gamma = 0.0165;
};

struct KStar892
{
	const std::string name = "KStar892";
	const float mass = 0.89166;
	const float gamma = 0.0508;
};

struct positron
{
	const std::string name = "positron";
	const float mass = 0.5110999e-3;
	const int geant_id = 2;
	const int id = 11;
	const int charge = 1;
};

struct electron
{
	const std::string name = "electron";
	const float mass = 0.5110999e-3;
	const int geant_id = 3;
	const int id = -11;
	const int charge = -1;
};

struct pion
{
	const std::string name = "pion";
	const float mass = 0.139570;
	const int geant_id = 8;
	const int id = 211;
	const int charge = 1;
};

struct apion
{
	const std::string name = "apion";
	const float mass = 0.139570;
	const int geant_id = 9;
	const int id = -211;
	const int charge = -1;
};

struct kaon
{
	const std::string name = "kaon";
	const float mass = 0.493677;
	const int geant_id = 11;
	const int id = 321;
	const int charge = 1;
};

struct akaon
{
	const std::string name = "akaon";
	const float mass = 0.493677;
	const int geant_id = 12;
	const int id = -321;
	const int charge = -1;
};

struct proton
{
	const std::string name = "proton";
	const float mass = 0.938272;
	const int geant_id = 14;
	const int id = 2212;
	const int charge = 1;
};

struct aproton
{
	const std::string name = "aproton";
	const float mass = 0.938272;
	const int geant_id = 15;
	const int id = -2212;
	const int charge = -1;
};
