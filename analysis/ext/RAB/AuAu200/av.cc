#include "../../../lib/ErrorHandler.h"
#include "../../../lib/OutputTool.h"

void WriteFile(std::string part, std::string min_centr, std::string max_centr, std::string target_centr, const double w1, const double w2)
{
	CheckInputFile(part + "_" + min_centr + ".txt");
	CheckInputFile(part + "_" + max_centr + ".txt");
	
	ifstream if1((part + "_" + min_centr + ".txt").c_str());
	ifstream if2((part + "_" + max_centr + ".txt").c_str());
	ofstream of((part + "_" + target_centr + ".txt").c_str());
	
	double x1, y1, err1, x2, y2, err2, serr1, serr2;
	while (if1 >> x1 >> y1 >> err1 >> serr1 && if2 >> x2 >> y2 >> err2 >> serr2)
	{
		double val = y1*w1+y2*w2;
		of << x1 << "	" << val << "	" << sqrt(pow(err1/y1*w1, 2) + pow(err2/y2*w2, 2))*val << "	" << serr1*w1 + serr2*w2 << std::endl;
	}
	of.close();
	PrintInfo("File " + part + "_" + target_centr + ".txt was written");
}

int av()
{
	std::string part = "proton";
	WriteFile(part, "0020", "2040", "0040", 0.5, 0.5);
	WriteFile(part, "4060", "6093", "4093", 0.5, 0.5);
	WriteFile(part, "0040", "4093", "0093", 0.5, 0.5);
	
	return 0;
}
