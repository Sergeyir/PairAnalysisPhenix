#pragma once

#include "TMath.h"
#include "TF1.h"

struct
{
	int number_of_iter = 1000;
	double sigma_conv_range = 10.;
} RBW_GDF_conv_par;
	
Double_t RBW_GDF_conv(Double_t* x, Double_t* par)
{
	double sum = 0.;
	double norm = 0.;
	
	int iter = 0;

	for (double t = - RBW_GDF_conv_par.sigma_conv_range*par[3]; 
		t < RBW_GDF_conv_par.sigma_conv_range*par[3]; 
		t += 2.*RBW_GDF_conv_par.sigma_conv_range*par[3]/
		static_cast<double>(RBW_GDF_conv_par.number_of_iter))
	{
		iter++;
		sum += TMath::Gaus(t, 0., par[3])*
			TMath::BreitWigner(x[0]-t, par[1], par[2]);
		norm += TMath::Gaus(t, 0., par[3])*
			TMath::BreitWigner(par[1]-t, par[1], par[2]);
	}

	return par[0]*sum/norm;
}

Double_t POL2(Double_t* x, Double_t* par)
{
	return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

Double_t POL3(Double_t* x, Double_t* par)
{
	return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
}

Double_t RBW_GDF_POL2(Double_t* x, Double_t* par)
{
	return RBW_GDF_conv(x, par) + POL2(x, &par[4]);
}

Double_t RBW_GDF_POL3(Double_t* x, Double_t* par)
{
	return RBW_GDF_conv(x, par) + POL3(x, &par[4]);
}
