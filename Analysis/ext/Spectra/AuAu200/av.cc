void WriteFile(std::string part, std::string min_centr, std::string max_centr, std::string target_centr, const double w1, const double w2)
{
	ifstream if1((part + "_" + min_centr + ".txt").c_str());
	ifstream if2((part + "_" + max_centr + ".txt").c_str());
	ofstream of((part + "_" + target_centr + ".txt").c_str());
	
	/*
	double x1, y1, err1, x2, y2, err2;
	while (if1 >> x1 >> y1 >> err1 && if2 >> x2 >> y2 >> err2)
	{
		of << x1 << "	" << y1*w1+y2*w2 << "	" << sqrt(pow(err1/y1*w1, 2) + pow(err2/y2*w2, 2))*(y1*w1+y2*w2) << std::endl;
	}
	*/
	
	double x1, y1, err1, x2, y2, err2, serr1, serr2;
	while (if1 >> x1 >> y1 >> err1 >> serr1 && if2 >> x2 >> y2 >> err2 >> serr2)
	{
		of << x1 << "	" << y1*w1+y2*w2 << "	" << sqrt(pow(err1, 2) + pow(err2, 2)) << "	" << sqrt(pow(serr1, 2) + pow(serr2, 2)) << std::endl;
	}
}

int av()
{
	/*
	WriteFile("pion", "0010", "1020", "0020", 0.5, 0.5);
	WriteFile("apion", "0010", "1020", "0020", 0.5, 0.5);
	WriteFile("kaon", "0010", "1020", "0020", 0.5, 0.5);
	WriteFile("akaon", "0010", "1020", "0020", 0.5, 0.5);
	WriteFile("proton", "0010", "1020", "0020", 0.5, 0.5);
	WriteFile("aproton", "0010", "1020", "0020", 0.5, 0.5);
	*/
	WriteFile("phi1020", "40-50", "50-60", "40-60", 0.5, 0.5);

	return 0;
}
