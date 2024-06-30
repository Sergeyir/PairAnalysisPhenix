int WriteFile(std::string part)
{
	ifstream if1((part + "_0010.txt").c_str());
	ifstream if2((part + "_1020.txt").c_str());
	ifstream if3((part + "_2040.txt").c_str());
	ifstream if4((part + "_4060.txt").c_str());
	ifstream if5((part + "_6093.txt").c_str());
	
	ofstream of((part + "_0093.txt").c_str());

	double x1, y1, err1;
	double x2, y2, err2;
	double x3, y3, err3;
	double x4, y4, err4;
	double x5, y5, err5;
	double w1 = 0.1, w2 = 0.1, w3 = 0.2, w4 = 0.2, w5 = 0.32;
	
	while	(
				if1 >> x1 >> y1 >> err1 && 
				if2 >> x2 >> y2 >> err2 &&
				if3 >> x3 >> y3 >> err3 &&
				if4 >> x4 >> y4 >> err4 &&
				if5 >> x5 >> y5 >> err5
			)
	{
		double val = y1*w1 + y2*w2 + y3*w3 + y4*w4 + y5*w5;
		of << x1 << "	" << val << "	" << sqrt(pow(err1/y1*w1, 2) + pow(err2/y2*w2, 2) + pow(err3/y3*w3, 2)+pow(err4/y4*w4, 2) + pow(err5/y5*w5, 2))*val << std::endl;
	}

	return 0;
}

void av_mb()
{
	WriteFile("pion");
	WriteFile("apion");
	WriteFile("kaon");
	WriteFile("akaon");
	WriteFile("proton");
	WriteFile("aproton");
}
