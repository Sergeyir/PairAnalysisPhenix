#pragma once

#include <iostream>
#include <string>
#include <fstream>

#include "OutputColor.h"

void PrintError(std::string message)
{
	std::cerr << OutputColor::bold_red << "ERROR: " << OutputColor::reset << message << std::endl;
	exit(1);
}

void PrintWarning(std::string message)
{
	std::cerr << OutputColor::bold_magenta << "WARNING: " << OutputColor::reset << message << std::endl;
}

void CheckInputFile(std::string name)
{
	std::ifstream file(name);
	if(!file.is_open()) PrintError("File " + name + " not found");
}

void CheckOutputFile(std::string name, std::ios_base::openmode openmode = std::ios::out)
{
	std::ofstream file(name, openmode);
	if(!file.is_open()) PrintError("File " + name + " cannot be created");
}
