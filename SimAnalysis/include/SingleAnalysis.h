#pragma once

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include "ThrObj.h"

#include "Tools.h"
#include "Particles.h"

bool IsHit(const double dval)
{
	if (dval < -900) return false;
	return true;
}

bool IsMatch(const double sdphi, const double sdz, const double max_sdphi = 2., const double max_sdz = 2.)
{
	if (abs(sdphi) > max_sdphi || abs(sdz) > max_sdz) return false;
	return true;
}	

double TransformProb(double prob)
{
	if (prob > 1.) prob = 1.;
	else if (prob < 0.) prob = 0.;
	return prob;
}

void CutDCDeadAreas(TH2F *hist, bool (*cut_func)(const double, const double, const double, const double), const double phi, const double zed)
{
	for (int i = 1; i <= hist->GetXaxis()->GetNbins(); i++)
	{
		for (int j = 1; j <= hist->GetYaxis()->GetNbins(); j++)
		{
			const double val_x = hist->GetXaxis()->GetBinCenter(i);
			const double val_y = hist->GetYaxis()->GetBinCenter(j);

			if (cut_func(phi, zed, val_x, val_y))
			{
				hist->SetBinContent(i, j, 0.);
			}
		}
	}
}

bool IsQualityCut(const int qual)
{
	if (qual != 63 && qual != 31) return true;
	return false;
}

double GetM2Mean(const double pt, const double *par)
{
	return par[0] + pt*par[1];
}

double GetM2Sigma(const double pt, const double m2_mean, const double *par)
{
	return 2.*sqrt(pow(par[0]/par[3]*m2_mean*pt, 2) + 
		pow(par[1]/par[3]*m2_mean, 2)*(1.+m2_mean/pt/pt) + 
		pow(par[2]*2.9972e-4/par[4]*pt, 2)*(m2_mean + pt*pt));
}

double GetTOFwsdphi(const int field, const float mom, const float tofwdphi, const int charge, const int strip)
{
	float p0 = -9999;
	float p1 = -9999;
	float p2 = -9999;
	float p3 = -9999;
	float mean = -9999;
	float sigma = -9999;
	float x = mom;

	if(field==1&&charge==-1)
	{
		p0 = -1.82527e-04;
		p1 = 2.01926e-04;
		p2 = -4.03119e-05;
		p3 = 1.78013e-04;
		if(strip<256) mean = p0 + p1/x + p2/x/x + p3/sqrt(x);
		p0 = -4.04259e-03;
		p1 = -2.76914e-04;
		p2 = 1.46522e-05;
		p3 = 7.66450e-04;
		if(strip>255) mean = p0 + p1/x + p2/x/x + p3/sqrt(x);
		p0 = 2.35259e-03; p1 = 8.17567e-04; p2 = 1.23846e-04; p3 = -6.09512e-04; if(strip<256) sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		p0 = 2.30133e-03;
		p1 = 5.64692e-04;
		p2 = 1.56190e-04;
		p3 = -4.27461e-04;
		if(strip>255) sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		return (tofwdphi - mean)/sigma;
	}
	if(field==1&&charge==1)
	{
		p0 = 3.40483e-04;
		p1 = 4.25949e-04;
		p2 = 1.15023e-05;
		p3 = -1.05502e-03;
		if(strip<256) mean = p0 + p1/x + p2/x/x + p3/sqrt(x);
		p0 = -3.42056e-03;
		p1 = 2.42465e-04;
		p2 = 2.94131e-05;
		p3 = -7.94733e-04;
		if(strip>255) mean = p0 + p1/x + p2/x/x + p3/sqrt(x);
		p0 = 2.35600e-03;
		p1 = 7.75032e-04;
		p2 = 1.37143e-04;
		p3 = -5.63775e-04;
		if(strip<256) sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		p0 = 2.39039e-03;
		p1 = 7.07798e-04;
		p2 = 1.31897e-04;
		p3 = -6.74388e-04;
		if(strip>255) sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		return (tofwdphi - mean)/sigma;
	}
	if(field==0&&charge==-1)
	{
		p0 = 1.66923e-04;
		p1 = 7.88784e-05;
		p2 = 1.01633e-05;
		p3 = -5.02453e-04;
		if(strip<256) mean = p0 + p1/x + p2/x/x + p3/sqrt(x);
		p0 = -3.62245e-03;
		p1 = -1.67680e-04;
		p2 = 2.62321e-05;
		p3 = -1.67802e-04;
		if(strip>255) mean = p0 + p1/x + p2/x/x + p3/sqrt(x);
		p0 = 2.38894e-03;
		p1 = 9.33554e-04;
		p2 = 1.07887e-04;
		p3 = -6.95579e-04;
		if(strip<256) sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		p0 = 2.19519e-03;
		p1 = 2.75664e-04;
		p2 = 1.77759e-04;
		p3 = -8.00175e-05;
		if(strip>255) sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		return (tofwdphi - mean)/sigma;
	}
	if(field==0&&charge==1)
	{
		p0 = -3.70869e-04;
		p1 = -8.25932e-05;
		p2 = -5.66926e-05;
		p3 = 6.93094e-04;
		if(strip<256) mean = p0 + p1/x + p2/x/x + p3/sqrt(x);
		p0 = -3.98785e-03;
		p1 = -4.51479e-05;
		p2 = -5.46694e-05;
		p3 = 5.45265e-04;
		if(strip>255) mean = p0 + p1/x + p2/x/x + p3/sqrt(x);
		p0 = 2.51701e-03;
		p1 = 1.38735e-03;
		p2 = 4.38010e-05;
		p3 = -1.21198e-03;
		if(strip<256) sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		p0 = 2.12479e-03;
		p1 = 2.35848e-04;
		p2 = 1.71448e-04;
		p3 = 8.64224e-05;
		if(strip>255) sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		return (tofwdphi - mean)/sigma;
	}

	return -9999;
}

double GetTOFwsdz(const int field, const double mom, const double tofwdz, const int charge, const int strip)
{
	float p0 = -9999;
	float p1 = -9999;
	float p2 = -9999;
	float p3 = -9999;
	float mean = -9999;
	float sigma = -9999;
	float x = mom;

	if(field==1&&charge==-1)
	{
		p0 = -1.54264e-01;
		p1 = -3.17810e-01;
		p2 = 1.99612e-02;
		p3 = 4.76464e-01;
		if(strip<256) mean = p0 + p1/x + p2/x/x + p3/sqrt(x);
		p0 = -1.16291e-01;
		p1 = -7.43883e-02;
		p2 = -8.07304e-04;
		p3 = 1.87733e-01;
		if(strip>255) mean = p0 + p1/x + p2/x/x + p3/sqrt(x);
		p0 = 1.59901e+00;
		p1 = 5.68304e-01;
		p2 = 4.03995e-02;
		p3 = -8.41961e-01;
		if(strip<256) sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		p0 = 1.49113e+00;
		p1 = 2.84760e-01;
		p2 = 5.76212e-02;
		p3 = -4.89066e-01;
		if(strip>255) sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		return (tofwdz - mean)/sigma;
	}
	if(field==1&&charge==1)
	{
		p0 = 3.22766e-02;
		p1 = 8.59880e-02;
		p2 = 1.38632e-02;
		p3 = -1.53894e-01;
		if(strip<256) mean = p0 + p1/x + p2/x/x + p3/sqrt(x);
		p0 = -6.30400e-02;
		p1 = 1.34487e-01;
		p2 = -3.32726e-02;
		p3 = -1.62875e-02;
		if(strip>255) mean = p0 + p1/x + p2/x/x + p3/sqrt(x);
		p0 = 1.65186e+00;
		p1 = 7.59812e-01;
		p2 = 1.05207e-02;
		p3 = -1.03914e+00;
		if(strip<256) sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		p0 = 1.41977e+00;
		p1 = 1.21789e-01;
		p2 = 7.75272e-02;
		p3 = -2.45346e-01;
		if(strip>255) sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		return (tofwdz - mean)/sigma;
	}
	if(field==0&&charge==-1)
	{
		p0 = -1.29874e-02;
		p1 = 2.20407e-02;
		p2 = 1.93202e-02;
		p3 = -3.96144e-02;
		if(strip<256) mean = p0 + p1/x + p2/x/x + p3/sqrt(x);
		p0 = 6.84030e-03;
		p1 = 1.40576e-01;
		p2 = -2.26645e-02;
		p3 = -1.28930e-01;
		if(strip>255) mean = p0 + p1/x + p2/x/x + p3/sqrt(x);
		p0 = 1.64226e+00;
		p1 = 7.20802e-01;
		p2 = 2.27525e-02;
		p3 = -1.01364e+00;
		if(strip<256) sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		p0 = 1.53972e+00;
		p1 = 3.22544e-01;
		p2 = 6.53821e-02;
		p3 = -5.84257e-01;
		if(strip>255) sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		return (tofwdz - mean)/sigma;
	}
	if(field==0&&charge==1)
	{
		p0 = 7.38252e-02;
		p1 = 1.55284e-01;
		p2 = -1.57147e-02;
		p3 = -2.05458e-01;
		if(strip<256) mean = p0 + p1/x + p2/x/x + p3/sqrt(x);
		p0 = -8.57157e-02;
		p1 = -1.47707e-01;
		p2 = 3.14808e-02;
		p3 = 1.75709e-01;
		if(strip>255) mean = p0 + p1/x + p2/x/x + p3/sqrt(x);
		p0 = 1.55151e+00;
		p1 = 5.19775e-01;
		p2 = 3.28126e-02;
		p3 = -7.37797e-01;
		if(strip<256) sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		p0 = 1.53939e+00;
		p1 = 3.39611e-01;
		p2 = 6.03544e-02;
		p3 = -6.03277e-01;
		if(strip>255) sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		return (tofwdz - mean)/sigma;
	}

	return -9999;
}

float RecalTOFwsdphi(const int field, const double mom, const double tofwsdphi, const int charge, const int strip)
{
	float p0 = -9999;
	float p1 = -9999;
	float p2 = -9999;
	float p3 = -9999;
	//float mean = -9999;
	float sigma = -9999;
	float x = mom;

	if(field==1&&charge==-1&&strip<256)
	{
		p0 = 1.0249;
		p1 = -0.0396688;
		p2 = 0.012552;
		p3 = -0.0542971;
		sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		return tofwsdphi/sigma;
	}
	if(field==1&&charge==-1&&strip>255)
	{
		p0 = 0.965167;
		p1 = -0.179558;
		p2 = 0.0260532;
		p3 = 0.137989;
		sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		return tofwsdphi/sigma;
	}
	if(field==1&&charge==1&&strip<256)
	{
		p0 = 1.02707;
		p1 = -0.0280164;
		p2 = 0.010193;
		p3 = -0.0678821;
		sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		return tofwsdphi/sigma;
	}
	if(field==1&&charge==1&&strip>255)
	{
		p0 = 0.965422;
		p1 = -0.122838;
		p2 = 0.0166137;
		p3 = 0.0958661;
		sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		return tofwsdphi/sigma;
	}
	if(field==0&&charge==-1&&strip<256)
	{
		p0 = 0.98246;
		p1 = -0.113141;
		p2 = 0.0204387;
		p3 = 0.0559133;
		sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		return tofwsdphi/sigma;
	}
	if(field==0&&charge==-1&&strip>255)
	{
		p0 = 1.03554;
		p1 = -0.0372324;
		p2 = 0.0134871;
		p3 = -0.0670396;
		sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		return tofwsdphi/sigma;
	}
	if(field==0&&charge==1&&strip<256)
	{
		p0 = 0.954387;
		p1 = -0.176233;
		p2 = 0.0256429;
		p3 = 0.137991;
		sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		return tofwsdphi/sigma;
	}
	if(field==0&&charge==1&&strip>255)
	{
		p0 = 1.0628;
		p1 = 0.0145043;
		p2 = 0.0075927;
		p3 = -0.139396;
		sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		return tofwsdphi/sigma;
	}

	return -9999;
}

float RecalTOFwsdz(const int field, const double mom, const double tofwsdz, const int charge, const int strip)
{
	float p0 = -9999;
	float p1 = -9999;
	float p2 = -9999;
	float p3 = -9999;
	//float mean = -9999;
	float sigma = -9999;
	float x = mom;

	if(field==1&&charge==-1&&strip<256)
	{
		p0 = 0.962477;
		p1 = -0.094472;
		p2 = 0.0292912;
		p3 = 0.0989633;
		sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		return tofwsdz/sigma;
	}
	if(field==1&&charge==-1&&strip>255)
	{
		p0 = 1.00028;
		p1 = -0.0790036;
		p2 = 0.0294143;
		p3 = 0.0411551;
		sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		return tofwsdz/sigma;
	}
	if(field==1&&charge==1&&strip<256)
	{
		p0 = 0.910849;
		p1 = -0.265478;
		p2 = 0.0547336;
		p3 = 0.288292;
		sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		return tofwsdz/sigma;
	}
	if(field==1&&charge==1&&strip>255)
	{
		p0 = 1.02504;
		p1 = -0.0398175;
		p2 = 0.0261403;
		p3 = -0.0239143;
		sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		return tofwsdz/sigma;
	}
	if(field==0&&charge==-1&&strip<256)
	{
		p0 = 0.882777;
		p1 = -0.294656;
		p2 = 0.0556337;
		p3 = 0.355325;
		sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		return tofwsdz/sigma;
	}
	if(field==0&&charge==-1&&strip>255)
	{
		p0 = 1.01236;
		p1 = 0.0102616;
		p2 = 0.0168994;
		p3 = -0.0400674;
		sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		return tofwsdz/sigma;
	}
	if(field==0&&charge==1&&strip<256)
	{
		p0 = 0.961488;
		p1 = -0.125693;
		p2 = 0.031693;
		p3 = 0.124763;
		sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		return tofwsdz/sigma;
	}
	if(field==0&&charge==1&&strip>255)
	{
		p0 = 0.973728;
		p1 = -0.115231;
		p2 = 0.0336785;
		p3 = 0.104168;
		sigma = p0 + p1/x + p2/x/x + p3/sqrt(x);
		return tofwsdz/sigma;
	}

	return -9999;
}
