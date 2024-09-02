#include <cmath>
#include <array>
#include <cmath>
#include <fstream>

#include "Tools.h"
#include "OutputColor.h"
#include "OutputTool.h"
#include "ParMS.h"
#include "StrTool.h"
#include "TCanvasPrinter.h"

const double pi = 3.14159265359;

struct FitResults
{
	double mean;
	double sigma;
};

struct FitPart
{
	TGraph *ms_shade;
	
	TGraph *means, *sigmas;
	TGraph *ms_3s_min, *ms_3s_max;

	TF1 *means_fit, *sigmas_fit;
	TF1 *ms_3s_fit_min, *ms_3s_fit_max;

	std::array<double, 200> ms_shade_min, ms_shade_max;

	ofstream raw_yield_outfile;
};

struct
{
	double y_min = 1e31;
	double y_max = -1e31;
} EffPar;

//gaus is the particles the yield of which is been measured
//gaus1 and gaus2 are partcles which are subtracted from m2 distribution
double GetYield(TH1F *hist, const double mean, const double sigma, TF1 *gaus1, TF1 *gaus2, double &err)
{
	double yield = 0;
	double yield_nosubtr = 0;
	
	for (int i = hist->GetXaxis()->FindBin(mean-Par.yield_extr_srange*sigma); 
			i <= hist->GetXaxis()->FindBin(mean+Par.yield_extr_srange*sigma); i++)
	{
		yield += hist->GetBinContent(i) - 
			gaus1->Eval(hist->GetXaxis()->GetBinCenter(i)) -
			gaus2->Eval(hist->GetXaxis()->GetBinCenter(i));
		yield_nosubtr += hist->GetBinContent(i);
	}
	
	const double norm = (
		erf((hist->GetXaxis()->GetBinCenter(hist->GetXaxis()->FindBin(mean)) - 
		hist->GetXaxis()->GetBinLowEdge(
		hist->GetXaxis()->FindBin(mean - Par.yield_extr_srange*sigma)))*
		Par.yield_extr_srange/sqrt(2.)/sigma) +
		erf((hist->GetXaxis()->GetBinUpEdge(
		hist->GetXaxis()->FindBin(mean + Par.yield_extr_srange*sigma)) - 
		hist->GetXaxis()->GetBinCenter(hist->GetXaxis()->FindBin(mean)))*
		Par.yield_extr_srange/sqrt(2.)/sigma))/2./erf(Par.yield_extr_srange/sqrt(2.));

	err = sqrt(yield_nosubtr)/yield/norm;
	
	return yield/norm;
}

//same as previous function but we also subtract the background
double GetYieldnoBG(TH1F *hist, const double mean, const double sigma, TF1 *gaus1, TF1 *gaus2, TF1 *bg_func, double &err)
{
	double yield = 0;
	double yield_nosubtr = 0;

	for (int i = hist->GetXaxis()->FindBin(mean-Par.yield_extr_srange*sigma); 
			i <= hist->GetXaxis()->FindBin(mean+Par.yield_extr_srange*sigma); i++)
	{
		yield += hist->GetBinContent(i) - 
			bg_func->Eval(hist->GetXaxis()->GetBinCenter(i)) - 
			gaus1->Eval(hist->GetXaxis()->GetBinCenter(i)) - 
			gaus2->Eval(hist->GetXaxis()->GetBinCenter(i));
		yield_nosubtr += hist->GetBinContent(i);
	}

	const double norm = (
		erf((hist->GetXaxis()->GetBinCenter(hist->GetXaxis()->FindBin(mean)) - 
		hist->GetXaxis()->GetBinLowEdge(
		hist->GetXaxis()->FindBin(mean - Par.yield_extr_srange*sigma)))*
		Par.yield_extr_srange/sqrt(2.)/sigma) +
		erf((hist->GetXaxis()->GetBinUpEdge(
		hist->GetXaxis()->FindBin(mean + Par.yield_extr_srange*sigma)) - 
		hist->GetXaxis()->GetBinCenter(hist->GetXaxis()->FindBin(mean)))*
		Par.yield_extr_srange/sqrt(2.)/sigma))/2./erf(Par.yield_extr_srange/sqrt(2.));

	err = sqrt(yield_nosubtr)/yield/norm;
	
	return yield;
}

//for calculating yields with veto applied
double GetYieldVeto(TH1F *hist, const double mean, const double sigma, TF1 *gaus1, TF1 *gaus2, double veto_low, double veto_up, double &err)
{
	double yield = 0;
	double yield_nosubtr = 0;
	
	for (int i = hist->GetXaxis()->FindBin(Maximum(mean-Par.yield_extr_srange*sigma, veto_low)); 
			i <= hist->GetXaxis()->FindBin(Minimum(mean+Par.yield_extr_srange*sigma, veto_up)); i++)
	{
		yield += hist->GetBinContent(i) - 
			gaus1->Eval(hist->GetXaxis()->GetBinCenter(i)) -
			gaus2->Eval(hist->GetXaxis()->GetBinCenter(i));
		yield_nosubtr += hist->GetBinContent(i);
	}
	
	const double norm = (
		erf((hist->GetXaxis()->GetBinCenter(hist->GetXaxis()->FindBin(mean)) - 
		hist->GetXaxis()->GetBinLowEdge(
		hist->GetXaxis()->FindBin(Maximum(mean - Par.yield_extr_srange*sigma, veto_low))))*
		Par.yield_extr_srange/sqrt(2.)/sigma) +
		erf((hist->GetXaxis()->GetBinUpEdge(
		hist->GetXaxis()->FindBin(Minimum(mean + Par.yield_extr_srange*sigma, veto_up))) - 
		hist->GetXaxis()->GetBinCenter(hist->GetXaxis()->FindBin(mean)))*
		Par.yield_extr_srange/sqrt(2.)/sigma))/2./erf(Par.yield_extr_srange/sqrt(2.));

	err = sqrt(yield_nosubtr)/yield/norm;
	
	return yield/norm;
}

//for calculating yields with veto applied
double GetYieldnoBGVeto(TH1F *hist, const double mean, const double sigma, TF1 *gaus1, TF1 *gaus2, TF1 *bg_func, double veto_low, double veto_up, double &err)
{
	double yield = 0;
	double yield_nosubtr = 0;
	
	for (int i = hist->GetXaxis()->FindBin(Maximum(mean-Par.yield_extr_srange*sigma, veto_low)); 
			i <= hist->GetXaxis()->FindBin(Minimum(mean+Par.yield_extr_srange*sigma, veto_up)); i++)
	{
		yield += hist->GetBinContent(i) - 
			bg_func->Eval(hist->GetXaxis()->GetBinCenter(i)) -
			gaus1->Eval(hist->GetXaxis()->GetBinCenter(i)) -
			gaus2->Eval(hist->GetXaxis()->GetBinCenter(i));
		yield_nosubtr += hist->GetBinContent(i);
	}

	const double norm = (
		erf((hist->GetXaxis()->GetBinCenter(hist->GetXaxis()->FindBin(mean)) - 
		hist->GetXaxis()->GetBinLowEdge(
		hist->GetXaxis()->FindBin(Maximum(mean - Par.yield_extr_srange*sigma, veto_low))))*
		Par.yield_extr_srange/sqrt(2.)/sigma) +
		erf((hist->GetXaxis()->GetBinUpEdge(
		hist->GetXaxis()->FindBin(Minimum(mean + Par.yield_extr_srange*sigma, veto_up))) - 
		hist->GetXaxis()->GetBinCenter(hist->GetXaxis()->FindBin(mean)))*
		Par.yield_extr_srange/sqrt(2.)/sigma))/2./erf(Par.yield_extr_srange/sqrt(2.));

	err = sqrt(yield_nosubtr)/yield/norm/2.;
	
	return yield/norm;
}

class MSC
{
	private:

	double part_min_x, part_max_x;
	double apart_min_x, apart_max_x;

	TF1 *pi_fit_func;
	TF1 *k_fit_func;
	TF1 *pi_bg_func;
	TF1 *k_bg_func;
	TF1 *p_fit_func;
	TF1 *p_bg_func;
	
	TF1 *pi_gaus;
	TF1 *k_gaus;
	TF1 *p_gaus;

	FitPart pion, kaon, proton, apion, akaon, aproton;

	FitResults fit_pion, fit_kaon, fit_proton, fit_deutron;

	std::array<std::string, 6> part_name = {"pion", "kaon", "proton", "apion", "akaon", "aproton"};

	void PerformSingleM2Fit(const double pt, TH1F *mass1d, Part pi, Part k, Part p, const double k_ptmax, int negc)
	{
		fit_pion.mean = Par.pion.mean;
		fit_kaon.mean = Par.kaon.mean;
		fit_proton.mean = Par.proton.mean;
		
		fit_pion.sigma = Par.def_pion_sigma[negc]->Eval(pt);
		fit_kaon.sigma = Par.def_kaon_sigma[negc]->Eval(pt);
		fit_proton.sigma = Par.def_proton_sigma[negc]->Eval(pt);
		
		//pions
		if (pt >= pi.ptmin && pt <= pi.ptmax)
		{
			for (int i = 0; i < pi_fit_func->GetNpar(); i++) pi_fit_func->ReleaseParameter(i);
			for (int i = 0; i < pi_fit_func->GetNpar(); i++) pi_fit_func->SetParameter(i, 1.);
				
			pi_fit_func->SetParameter(3, mass1d->GetBinContent(mass1d->FindBin(fit_pion.mean))/1.1);
			pi_fit_func->SetParameter(4, fit_pion.mean);
			pi_fit_func->SetParameter(5, fit_pion.sigma); 
			
			//max height of gaus
			pi_fit_func->SetParLimits(3, 
				mass1d->GetBinContent(mass1d->FindBin(fit_pion.mean))/2., 
				mass1d->GetBinContent(mass1d->FindBin(fit_pion.mean))*2.);
			
			if (Par.do_fix_mean_par)
			{
				pi_fit_func->FixParameter(4, Par.pion.mean);
			}
			else
			{
				pi_fit_func->SetParLimits(4, 
					fit_pion.mean-fit_pion.sigma, 
					fit_pion.mean+fit_pion.sigma);
			}
			//means and sigmas
			if (Par.do_fix_sigma_par)
			{	
				pi_fit_func->FixParameter(5, fit_pion.sigma);
			}
			else
			{
				pi_fit_func->SetParLimits(5, 
					fit_pion.sigma/1.5, fit_pion.sigma*1.5);
			}
			
			pi_fit_func->SetRange(fit_pion.mean - 1.*fit_pion.sigma, 
				fit_pion.mean + 1.*fit_pion.sigma);
			
			mass1d->Fit(pi_fit_func, "RQMBN");
			
			fit_pion.mean = pi_fit_func->GetParameter(4);
			fit_pion.sigma = abs(pi_fit_func->GetParameter(5));
			
			pi_fit_func->SetRange(fit_pion.mean - fit_pion.sigma, 
				fit_pion.mean + fit_pion.sigma);

			for (int i = 0; i < Par.iter; i++)
			{
				double mult = 1. - static_cast<double>(i)/Par.iter*1.;

				pi_fit_func->SetParameter(4, 
					(fit_pion.mean - fit_pion.sigma*(mult - 1.) + 
					fit_pion.mean + fit_pion.sigma*(mult - 1.))/2.);
				
				if (!Par.do_fix_mean_par)
				{
					pi_fit_func->SetParLimits(4, 
						Maximum(fit_pion.mean - fit_pion.sigma*mult, Par.pi_m2_min), 
						Minimum(fit_pion.mean + fit_pion.sigma*mult, Par.pi_m2_max));
				}
				
				if (!Par.do_fix_sigma_par)
				{
					pi_fit_func->SetParLimits(5, 
						fit_pion.sigma*(1. - mult/10.), fit_pion.sigma*(1. + mult/10.));
				}

				mass1d->Fit(pi_fit_func, "RQMBN");

				fit_pion.mean = abs(pi_fit_func->GetParameter(4));
				fit_pion.sigma = abs(pi_fit_func->GetParameter(5));
				
				pi_fit_func->SetRange(fit_pion.mean-fit_pion.sigma,
					fit_pion.mean + fit_pion.sigma);
			}

			pi_fit_func->SetRange(fit_pion.mean-fit_pion.sigma*Par.draw_srange,
				fit_pion.mean+fit_pion.sigma*Par.draw_srange);
			pi_bg_func->SetRange(fit_pion.mean-fit_pion.sigma*Par.draw_srange,
				fit_pion.mean+fit_pion.sigma*Par.draw_srange);
			
			pi_gaus->SetRange(fit_pion.mean-fit_pion.sigma*5., 
				fit_pion.mean+fit_pion.sigma*5.);

			for (int i = 0; i < pi_fit_func->GetNpar(); i++)
			{
				pi_fit_func->FixParameter(i, pi_fit_func->GetParameter(i));
				
				if (i < 3) pi_bg_func->FixParameter(i, pi_fit_func->GetParameter(i));
				else pi_gaus->FixParameter(i-3, pi_fit_func->GetParameter(i));
			}

			mass1d->Fit(pi_fit_func, "RQB+");
			mass1d->Fit(pi_bg_func, "RQB+");
			mass1d->Fit(pi_gaus, "RQB+");
		}
		
		//kaons
		if (pt >= k.ptmin && pt < k.ptmax)
		{
			for (int i = 0; i < k_fit_func->GetNpar(); i++) k_fit_func->ReleaseParameter(i);

			for (int i = 0; i < k_fit_func->GetNpar(); i++) k_fit_func->SetParameter(i, 1.);
				
			k_fit_func->SetParameter(6, mass1d->GetBinContent(mass1d->FindBin(fit_kaon.mean))/1.3);
			k_fit_func->SetParameter(7, fit_kaon.mean);
			k_fit_func->SetParameter(8, fit_kaon.sigma);

			//background
			k_fit_func->SetParameter(0, 0.);
			k_fit_func->SetParameter(1, 0.);
			k_fit_func->SetParameter(2, -1.);
			
			//correcting the first bins of background approximation since kaon signal is very small
			if ((Par.Det.name == "EMCale2" || Par.Det.name == "EMCale3") && pt < 0.4) 
				k_fit_func->SetParLimits(0, 10., 11.);
			else if (Par.Det.name == "EMCalw0" && pt < 0.4) 
				k_fit_func->SetParLimits(0, 10.5, 12.);
			else if (Par.Det.name == "EMCalw1" && pt < 0.5) 
				k_fit_func->SetParLimits(0, 10., 12.);
			else if ((Par.Det.name == "EMCalw2" || Par.Det.name == "EMCalw3") && pt < 0.4) 
				k_fit_func->SetParLimits(0, 10.7, 12.);
			else k_fit_func->SetParLimits(0, 0., 20);
			k_fit_func->SetParLimits(1, -1e2, 1e2);
			k_fit_func->SetParLimits(2, -20., 0.);

			//max height of gaus
			k_fit_func->SetParLimits(6, 
				mass1d->GetBinContent(mass1d->FindBin(fit_kaon.mean))/1.4, 
				mass1d->GetBinContent(mass1d->FindBin(fit_kaon.mean))*2.);

			if (Par.do_fix_mean_par)
			{
				k_fit_func->FixParameter(7, Par.kaon.mean);
			}
			else
			{
				k_fit_func->SetParLimits(7, 
					fit_kaon.mean - fit_kaon.sigma, 
					fit_kaon.mean + fit_kaon.sigma);
			}
			//means and sigmas
			if (Par.do_fix_sigma_par)
			{	
				k_fit_func->FixParameter(8, fit_kaon.sigma);
			}
			else
			{
				k_fit_func->SetParLimits(8, 
					fit_kaon.sigma/1.5, fit_kaon.sigma*1.5);
			}

			double upper_pt_range = 0.4;
			if (pt > 0.5) upper_pt_range = 0.55;
			else if (pt > 0.9) upper_pt_range = 0.55;
			
			k_fit_func->SetRange(
				Maximum(fit_pion.mean + fit_pion.sigma*1.5,
				fit_kaon.mean - fit_kaon.sigma), upper_pt_range);

			mass1d->Fit(k_fit_func, "RQMBN");
			
			fit_kaon.mean = k_fit_func->GetParameter(7);
			fit_kaon.sigma = abs(k_fit_func->GetParameter(8));

			k_fit_func->SetRange(
				Maximum(fit_pion.mean + fit_pion.sigma*1.5, 
				fit_kaon.mean - fit_kaon.sigma*5.),
				Minimum(fit_proton.mean - fit_proton.sigma*3., upper_pt_range));

			for (int i = 0; i < Par.iter; i++)
			{
				double mult = 1. - static_cast<double>(i)/Par.iter*1.;

				if (!Par.do_fix_mean_par)
				{
					k_fit_func->SetParLimits(7, 
						Maximum(fit_kaon.mean - fit_kaon.sigma*mult, Par.k_m2_min), 
						Minimum(fit_kaon.mean + fit_kaon.sigma*mult, Par.k_m2_max));	
				}
				
				if (!Par.do_fix_sigma_par)
				{
					k_fit_func->SetParLimits(8, 
						fit_kaon.sigma*(1. - mult/10.), fit_kaon.sigma*(1. + mult/10.));
				}

				mass1d->Fit(k_fit_func, "RQMBN");

				fit_kaon.mean = abs(k_fit_func->GetParameter(7));
				fit_kaon.sigma = abs(k_fit_func->GetParameter(8));

				k_fit_func->SetRange(
					Maximum(fit_pion.mean + fit_pion.sigma*1.5, 
					fit_kaon.mean - fit_kaon.sigma*8.),
					Minimum(fit_proton.mean - fit_proton.sigma*Par.veto_srange, upper_pt_range));
			}

			k_fit_func->SetRange(
				Maximum(fit_pion.mean + fit_pion.sigma*1.5, 
				fit_kaon.mean - fit_kaon.sigma*Par.draw_srange),
				fit_kaon.mean + fit_kaon.sigma*Par.draw_srange);
			k_bg_func->SetRange(
				Maximum(fit_pion.mean + fit_pion.sigma*1.5, 
				fit_kaon.mean - fit_kaon.sigma*Par.draw_srange), 
				fit_kaon.mean + fit_kaon.sigma*Par.draw_srange);
		
			k_gaus->SetRange(
				fit_kaon.mean - fit_kaon.sigma*5., 
				fit_kaon.mean + fit_kaon.sigma*5.);

			for (int i = 0; i < k_fit_func->GetNpar(); i++)
			{
				k_fit_func->FixParameter(i, k_fit_func->GetParameter(i));
				
				if (i < 6) k_bg_func->FixParameter(i, k_fit_func->GetParameter(i));
				else k_gaus->FixParameter(i-6, k_fit_func->GetParameter(i));
			}
			
			mass1d->Fit(k_fit_func, "RQB+");
			mass1d->Fit(k_bg_func, "RQB+");
			mass1d->Fit(k_gaus, "RQB+");
		}

		if (pt >= p.ptmin && pt <= p.ptmax)
		{
			for (int i = 0; i < p_fit_func->GetNpar(); i++) p_fit_func->ReleaseParameter(i);
				
			//protons
			p_fit_func->SetParameters(1., -0.1);
			
			p_fit_func->SetParameter(4, mass1d->GetBinContent(mass1d->FindBin(fit_proton.mean))/1.1);
			p_fit_func->SetParameter(5, fit_proton.mean);
			p_fit_func->SetParameter(6, fit_proton.sigma);

			p_fit_func->SetParLimits(4, 
				mass1d->GetBinContent(mass1d->FindBin(fit_proton.mean))/10., 
				mass1d->GetBinContent(mass1d->FindBin(fit_proton.mean)));
			
			p_fit_func->SetParLimits(5, 
				fit_proton.mean - fit_proton.sigma*3., 
				fit_proton.mean + fit_proton.sigma*3.);
			
			if (Par.do_fix_sigma_par)
			{
				p_fit_func->FixParameter(6, fit_proton.sigma);
			}
			else
			{
				p_fit_func->SetParLimits(6, 
					fit_proton.sigma/1.5, fit_proton.sigma*1.5);
			}

			p_fit_func->SetRange(fit_proton.mean - 1.*fit_proton.sigma,
				fit_proton.mean + 1.*fit_proton.sigma);
							
			mass1d->Fit(p_fit_func, "RQMBN");

			fit_proton.mean = abs(p_fit_func->GetParameter(5));
			fit_proton.sigma = abs(p_fit_func->GetParameter(6));

			p_fit_func->SetRange(fit_proton.mean-fit_proton.sigma*Par.srange*1.5,
				fit_proton.mean+fit_proton.sigma*Par.srange*2.);
			
			for (int i = 0; i < Par.iter; i++)
			{
				double mult = 1.5 - static_cast<double>(i)/Par.iter*0.5;
				
				p_fit_func->SetParameter(5, 
					(fit_proton.mean - fit_proton.sigma*(mult - 1.) + 
					fit_proton.mean + fit_proton.sigma*(mult - 1.))/2.);

				p_fit_func->SetParLimits(5, 
					fit_proton.mean - fit_proton.sigma*(mult - 1.), 
					fit_proton.mean + fit_proton.sigma*(mult - 1.));

				if (!Par.do_fix_sigma_par)
				{
					p_fit_func->SetParLimits(6, fit_proton.sigma/mult, fit_proton.sigma*mult);
				}
				
				mass1d->Fit(p_fit_func, "RQMBN");

				fit_proton.mean = abs(p_fit_func->GetParameter(5));
				fit_proton.sigma = abs(p_fit_func->GetParameter(6));

				p_fit_func->SetRange(fit_proton.mean-fit_proton.sigma*Par.srange*1.5,
					fit_proton.mean+fit_proton.sigma*Par.srange*2.);
			}
			
			p_bg_func->SetRange(fit_proton.mean-fit_proton.sigma*Par.draw_srange,
				fit_proton.mean+fit_proton.sigma*Par.draw_srange);
		
			p_gaus->SetRange(fit_proton.mean-fit_proton.sigma*5., 
				fit_proton.mean+fit_proton.sigma*5.);
			
			for (int i = 0; i < p_fit_func->GetNpar(); i++)
			{
				if (i < 4) p_bg_func->FixParameter(i, p_fit_func->GetParameter(i));
				p_gaus->FixParameter(i-4, p_fit_func->GetParameter(i));
				p_fit_func->FixParameter(i, p_fit_func->GetParameter(i));
			}

			mass1d->Fit(p_fit_func, "RQB+");
			mass1d->Fit(p_bg_func, "RQB+");
			mass1d->Fit(p_gaus, "RQB+");
		}
		
		if (pt > Par.write_yield_ptmin)
		{
			double pion_yield_err, kaon_yield_err, proton_yield_err;
			double pion_yield_err_veto, kaon_yield_err_veto, proton_yield_err_veto;
			
			const double pion_yield = GetYield(mass1d, 
				fit_pion.mean, fit_pion.sigma, k_gaus, p_gaus, pion_yield_err);
			
			const double pion_yield_veto = GetYieldVeto(mass1d, 
				fit_pion.mean, fit_pion.sigma, k_gaus, p_gaus, 
				-999, fit_kaon.mean-fit_kaon.sigma*Par.veto_srange, 
				pion_yield_err_veto);
			
			const double kaon_yield = GetYieldnoBG(mass1d, 
				fit_kaon.mean, fit_kaon.sigma, pi_gaus, p_gaus, k_bg_func, kaon_yield_err);
			
			const double kaon_yield_veto = GetYieldnoBGVeto(mass1d, 
				fit_kaon.mean, fit_kaon.sigma, pi_gaus, p_gaus, k_bg_func, 
				fit_pion.mean + fit_pion.sigma*Par.veto_srange,
				fit_proton.mean - fit_proton.sigma*Par.veto_srange,
				kaon_yield_err_veto);
			
			const double proton_yield = GetYield(mass1d, fit_proton.mean, fit_proton.sigma, 
				pi_gaus, k_gaus, proton_yield_err);
				
			if (negc == 0)
			{
				if (pt >= pi.ptmin && pt <= pi.ptmax)
				{
					pion.raw_yield_outfile << pt << " " << 
						pion_yield << " " << 
						pion_yield_err << " " << 
						pion_yield_veto << " " << 
						pion_yield_err_veto << std::endl;
				}
				if (pt >= k.ptmin && pt <= k.ptmax)
				{
					kaon.raw_yield_outfile << pt << " " << 
						kaon_yield << " " <<
						kaon_yield_err << " " << 
						kaon_yield_veto << " " << 
						kaon_yield_err_veto << std::endl;
				}
				if (pt >= p.ptmin && pt <= p.ptmax)
				{
					proton.raw_yield_outfile << pt << " " << 
						proton_yield << " " << 
						proton_yield_err << std::endl;
				}
			}
			else
			{
				if (pt >= pi.ptmin && pt <= pi.ptmax)
				{
					apion.raw_yield_outfile << pt << " " << 
						pion_yield << " " <<
						pion_yield_err << " " << 
						pion_yield_veto << " " << 
						pion_yield_err_veto << std::endl;
				}
				if (pt >= k.ptmin && pt <= k.ptmax)
				{
					akaon.raw_yield_outfile << pt << " " << 
						kaon_yield << " " << 
						kaon_yield_err << " " << 
						kaon_yield_veto << " " << 
						kaon_yield_err_veto << std::endl;
				}
				if (pt >= p.ptmin && pt <= p.ptmax)
				{
					aproton.raw_yield_outfile << pt << " " << 
						proton_yield << " " << 
						proton_yield_err << std::endl;
				}
			}
			Print(pt, pion_yield_err, kaon_yield_err);
		}
	}
	
	void AddMSPoint(FitPart &part, double pt, double mean, double sigma)
	{
		part.means->AddPoint(pt, mean);
		part.sigmas->AddPoint(pt, sigma);
		
		part.ms_3s_min->AddPoint(pt, mean-3.*sigma);
		part.ms_3s_max->AddPoint(pt, mean+3.*sigma);
	}

	void PerformSingleMSFit(FitPart &part, int parc)
	{
		part.means->Fit(part.means_fit, "RQMB+");
		
		//fixing m2 means
		for (int i = 0; i < part.means_fit->GetNpar(); i++)
		{
			part.sigmas_fit->FixParameter(i, part.means_fit->GetParameter(i));
		}

		part.sigmas_fit->SetParameter(2, Par.Det.sigma_alpha);
		part.sigmas_fit->SetParameter(3, Par.Det.sigma_ms);
		part.sigmas_fit->SetParameter(4, Par.Det.sigma_t);
		part.sigmas_fit->FixParameter(5, Par.K1);
		part.sigmas_fit->FixParameter(6, Par.Det.L);

		if (Par.do_fix_sigma_par)
		{
			part.sigmas_fit->SetParLimits(2, Par.Det.sigma_alpha/1.01, Par.Det.sigma_alpha*1.01);
			part.sigmas_fit->SetParLimits(3, Par.Det.sigma_ms/1.01, Par.Det.sigma_ms*1.01);
			part.sigmas_fit->SetParLimits(4, Par.Det.sigma_t/1.01, Par.Det.sigma_t*1.01);
		}
		else
		{
			part.sigmas_fit->SetParLimits(2, Par.Det.sigma_alpha/2., Par.Det.sigma_alpha*2.);
			part.sigmas_fit->SetParLimits(3, Par.Det.sigma_ms/2., Par.Det.sigma_ms*2.);
			part.sigmas_fit->SetParLimits(4, Par.Det.sigma_t/2., Par.Det.sigma_t*2.);
		}
		
		part.sigmas->Fit(part.sigmas_fit, "RQMB+");
		
		for (int i = 0; i < part.sigmas_fit->GetNpar(); i++)
		{
			part.ms_3s_fit_min->FixParameter(i, part.sigmas_fit->GetParameter(i));
			part.ms_3s_fit_max->FixParameter(i, part.sigmas_fit->GetParameter(i));
		}

		part.ms_3s_min->Fit(part.ms_3s_fit_min, "RQB+");
		part.ms_3s_max->Fit(part.ms_3s_fit_max, "RQB+");
	}

	void SetSingleShade(FitPart &part, double min_x, double max_x)
	{
		int points = 0;
		for (int i = 0; i < part.ms_shade_min.size(); i++)
		{
			if (abs(part.ms_shade_max[i]) < 10 && abs(part.ms_shade_min[i]) < 10 
				&& part.ms_shade_max[i] > part.ms_shade_min[i])
			{
				double x = min_x + (max_x - min_x)*i/(part.ms_shade_min.size());
				part.ms_shade->AddPoint(x, part.ms_shade_max[i]);
				points++;
			}
		}

		for (int i = part.ms_shade_min.size() - 1; i >= 0; i--)
		{
			if (abs(part.ms_shade_max[i]) < 10 && abs(part.ms_shade_min[i]) < 10 && 
				part.ms_shade_max[i] > part.ms_shade_min[i])
			{
				double x = min_x + (max_x - min_x)*i/(part.ms_shade_min.size());
				part.ms_shade->AddPoint(x, part.ms_shade_min[i]);
			}
		}

	}

	void SetShades()
	{
		for (double i = 0; i < pion.ms_shade_min.size(); i++)
		{
			double x = part_min_x + (part_max_x - part_min_x)*i/(pion.ms_shade_min.size());

			pion.ms_shade_max[i] = Minimum(
				pion.means_fit->Eval(x)+pion.sigmas_fit->Eval(x)*Par.yield_extr_srange, 
				kaon.ms_3s_fit_min->Eval(x), proton.ms_3s_fit_min->Eval(x));
			pion.ms_shade_min[i] = 
				pion.means_fit->Eval(x)-pion.sigmas_fit->Eval(x)*Par.yield_extr_srange;
			
			if (proton.ms_3s_fit_min->Eval(x) >= pion.ms_3s_fit_max->Eval(x))
			{
				kaon.ms_shade_max[i] = Minimum(
					kaon.means_fit->Eval(x)+kaon.sigmas_fit->Eval(x)*Par.yield_extr_srange, 
					proton.ms_3s_fit_min->Eval(x));
				kaon.ms_shade_min[i] = Maximum(
					kaon.means_fit->Eval(x)-kaon.sigmas_fit->Eval(x)*Par.yield_extr_srange, 
					pion.ms_3s_fit_max->Eval(x));
			}
			
			proton.ms_shade_max[i] = 
				proton.means_fit->Eval(x)+proton.sigmas_fit->Eval(x)*Par.yield_extr_srange;
			proton.ms_shade_min[i] = Maximum(
				proton.means_fit->Eval(x)-proton.sigmas_fit->Eval(x)*Par.yield_extr_srange, 
				pion.ms_3s_fit_max->Eval(x), kaon.ms_3s_fit_max->Eval(x));
		}

		for (double i = 0; i < apion.ms_shade_min.size(); i++)
		{
			double x = apart_min_x + (apart_max_x - apart_min_x)*i/(pion.ms_shade_min.size());

			apion.ms_shade_max[i] = Minimum(apion.means_fit->Eval(x)+apion.sigmas_fit->Eval(x), 
				akaon.ms_3s_fit_min->Eval(x), aproton.ms_3s_fit_min->Eval(x));
			apion.ms_shade_min[i] = apion.means_fit->Eval(x)-apion.sigmas_fit->Eval(x);
			
			if (aproton.ms_3s_fit_min->Eval(x) >= apion.ms_3s_fit_max->Eval(x))
			{
				akaon.ms_shade_max[i] = Minimum(akaon.means_fit->Eval(x)+akaon.sigmas_fit->Eval(x), 
					aproton.ms_3s_fit_min->Eval(x));
				akaon.ms_shade_min[i] = Maximum(akaon.means_fit->Eval(x)-akaon.sigmas_fit->Eval(x), 
					apion.ms_3s_fit_max->Eval(x));
			}
			
			aproton.ms_shade_max[i] = aproton.means_fit->Eval(x)+aproton.sigmas_fit->Eval(x);
			aproton.ms_shade_min[i] = Maximum(aproton.means_fit->Eval(x)-aproton.sigmas_fit->Eval(x), 
				apion.ms_3s_fit_max->Eval(x), akaon.ms_3s_fit_max->Eval(x));
		}

		SetSingleShade(pion, part_min_x, part_max_x);
		SetSingleShade(kaon, part_min_x, part_max_x);
		SetSingleShade(proton, part_min_x, part_max_x);

		SetSingleShade(apion, apart_min_x, apart_max_x);
		SetSingleShade(akaon, apart_min_x, apart_max_x);
		SetSingleShade(aproton, apart_min_x, apart_max_x);
	}

	void DrawSingleMS(FitPart &part)
	{
		part.means->Draw("P");
		part.ms_3s_min->Draw("P");
		part.ms_3s_max->Draw("P");
		part.ms_3s_fit_min->Draw("SAME L");
		part.ms_3s_fit_max->Draw("SAME L");
	}

	public:
	
	void SetPartStyle(FitPart &part, Color_t color, double ptmin, double ptmax, const int parc)
	{	
		part.means = new TGraph();
		part.sigmas = new TGraph();
		
		part.ms_3s_min = new TGraph();
		part.ms_3s_max = new TGraph();
		
		part.means_fit = new TF1("m", Par.mean_fit.c_str());
		part.sigmas_fit = new TF1("s", Par.sigma_fit.c_str());

		part.sigmas_fit->SetParameters(0, 0, Par.Det.sigma_alpha, 
			Par.Det.sigma_ms, Par.Det.sigma_t, Par.Det.L, Par.K1);
		
		//part.sigmas_fit->FixParameter(part.means_fit->GetNpar(), Par.sigma_alpha);
		//part.sigmas_fit->FixParameter(part.means_fit->GetNpar()+1, Par.sigma_ms);
		//part.sigmas_fit->FixParameter(part.means_fit->GetNpar()+2, Par.sigma_t);
		
		part.sigmas_fit->SetParLimits(part.means_fit->GetNpar(), 
			Par.Det.sigma_alpha/Par.Det.max_ms_dev, Par.Det.sigma_alpha*Par.Det.max_ms_dev);
		part.sigmas_fit->SetParLimits(part.means_fit->GetNpar()+1, 
			Par.Det.sigma_ms/Par.Det.max_ms_dev, Par.Det.sigma_ms*Par.Det.max_ms_dev);
		part.sigmas_fit->SetParLimits(part.means_fit->GetNpar()+2, 
			Par.Det.sigma_t/Par.Det.max_ms_dev, Par.Det.sigma_t*Par.Det.max_ms_dev);
		
		part.ms_3s_fit_min = new TF1("m-3s", (Par.mean_fit + 
			" - 3.*(" + Par.sigma_fit + ")").c_str());
		part.ms_3s_fit_max = new TF1("m+3s", (Par.mean_fit +
			" + 3.*(" + Par.sigma_fit + ")").c_str());

		part.ms_shade = new TGraph();

		part.ms_shade->SetFillStyle(1001);
		part.ms_shade->SetFillColorAlpha(color, 0.3);
		part.ms_shade->SetLineColor(color+1);

		part.means->SetMarkerSize(2);
		part.sigmas->SetMarkerSize(2);

		part.ms_3s_min->SetMarkerSize(2);
		part.ms_3s_max->SetMarkerSize(2);

		part.means->SetMarkerStyle(52);
		part.sigmas->SetMarkerStyle(52);

		part.ms_3s_min->SetMarkerStyle(65);
		part.ms_3s_max->SetMarkerStyle(65);
		
		part.means->SetMarkerColor(color+3);
		part.sigmas->SetMarkerColor(color+3);

		part.ms_3s_min->SetMarkerColor(color+2);
		part.ms_3s_max->SetMarkerColor(color+2);
		
		part.means_fit->SetLineColor(color+3);
		part.sigmas_fit->SetLineColor(color+3);

		part.ms_3s_fit_min->SetLineColor(color+1);
		part.ms_3s_fit_max->SetLineColor(color+1);

		part.means_fit->SetRange(ptmin, ptmax);
		part.sigmas_fit->SetRange(ptmin, ptmax);

		part.ms_3s_fit_min->SetRange(ptmin, ptmax);
		part.ms_3s_fit_max->SetRange(ptmin, ptmax);
	}
	
	MSC(double p_ptmin, double p_ptmax, double ap_ptmin, double ap_ptmax)
	{	
		system(("mkdir -p ../data/Efficiency/" + Par.run_name + "/" + Par.Det.name).c_str());

		pi_gaus = new TF1("pi_gaus", "gaus");
		k_gaus = new TF1("k_gaus", "gaus");
		p_gaus = new TF1("p_gaus", "gaus");
		
		//pions
		pi_fit_func = new TF1("pi_fit_func", 
			(Par.Det.pi_bg_func + "+gaus(3)").c_str());
		pi_bg_func = new TF1("pi_bg_func", Par.Det.pi_bg_func.c_str());
		
		//kaons
		k_fit_func = new TF1("k_fit_func", 
			(Par.Det.k_bg_func + "+gaus(6)").c_str());
		k_bg_func = new TF1("k_bg_func", Par.Det.k_bg_func.c_str());
		
		//protons
		p_fit_func = new TF1("p_fit_func", (Par.Det.p_bg_func + " + gaus(4)").c_str());
		p_bg_func = new TF1("p_bg_func", Par.Det.p_bg_func.c_str());

		pi_gaus->SetLineColor(Par.pion.color-3);
		k_gaus->SetLineColor(Par.kaon.color-3);
		p_gaus->SetLineColor(Par.proton.color-3);

		pi_fit_func->SetLineColor(Par.pion.color-5);
		k_fit_func->SetLineColor(Par.kaon.color-5);
		p_fit_func->SetLineColor(Par.proton.color-5);
		
		pi_bg_func->SetLineColor(kGray+1);
		k_bg_func->SetLineColor(kGray+1);
		p_bg_func->SetLineColor(kGray+1);

		pi_bg_func->SetLineStyle(2);
		k_bg_func->SetLineStyle(2);
		p_bg_func->SetLineStyle(2);

		pi_fit_func->SetLineWidth(3);
		k_fit_func->SetLineWidth(3);
		p_fit_func->SetLineWidth(3);
		pi_bg_func->SetLineWidth(3);
		k_bg_func->SetLineWidth(3);
		p_bg_func->SetLineWidth(3);

		pi_gaus->SetLineWidth(3);
		k_gaus->SetLineWidth(3);
		p_gaus->SetLineWidth(3);

		part_min_x = p_ptmin;
		part_max_x = p_ptmax;
		apart_min_x = ap_ptmin;
		apart_max_x = ap_ptmax;
		
		pi_gaus->SetLineWidth(2);
		k_gaus->SetLineWidth(2);
		p_gaus->SetLineWidth(2);

		SetPartStyle(pion, Par.pion.color, part_min_x, part_max_x, 0);
		SetPartStyle(kaon, Par.kaon.color, part_min_x, part_max_x, 1);
		SetPartStyle(proton, Par.proton.color, part_min_x, part_max_x, 2);

		SetPartStyle(apion, Par.apion.color, apart_min_x, apart_max_x, 3);
		SetPartStyle(akaon, Par.akaon.color, apart_min_x, apart_max_x, 4);
		SetPartStyle(aproton, Par.aproton.color, apart_min_x, apart_max_x, 5);

		for (double &val : pion.ms_shade_min) val = -9999;
		for (double &val : kaon.ms_shade_min) val = -9999;
		for (double &val : proton.ms_shade_min) val = -9999;
		
		for (double &val : apion.ms_shade_min) val = -9999;
		for (double &val : akaon.ms_shade_min) val = -9999;
		for (double &val : aproton.ms_shade_min) val = -9999;

		for (double &val : pion.ms_shade_max) val = -9999;
		for (double &val : kaon.ms_shade_max) val = -9999;
		for (double &val : proton.ms_shade_max) val = -9999;
		
		for (double &val : apion.ms_shade_max) val = -9999;
		for (double &val : akaon.ms_shade_max) val = -9999;
		for (double &val : aproton.ms_shade_max) val = -9999;

		system(("mkdir -p ../data/RawYields/" + Par.run_name + "/" + Par.Det.name).c_str());
		
		pion.raw_yield_outfile.open("../data/RawYields/" + 
			Par.run_name + "/" + Par.Det.name + "/" + Par.CType.cmin_name[Par.cnum] + 
			"_" + Par.CType.cmax_name[Par.cnum] + "_pion" + Par.magf + ".txt");
		kaon.raw_yield_outfile.open("../data/RawYields/" + 
			Par.run_name + "/" + Par.Det.name + "/" + Par.CType.cmin_name[Par.cnum] + 
			"_" + Par.CType.cmax_name[Par.cnum] + "_kaon" + Par.magf + ".txt");
		proton.raw_yield_outfile.open("../data/RawYields/" + 
			Par.run_name + "/" + Par.Det.name + "/" + Par.CType.cmin_name[Par.cnum] + 
			"_" + Par.CType.cmax_name[Par.cnum] + "_proton" + Par.magf + ".txt");
		apion.raw_yield_outfile.open("../data/RawYields/" + 
			Par.run_name + "/" + Par.Det.name + "/" + Par.CType.cmin_name[Par.cnum] + 
			"_" + Par.CType.cmax_name[Par.cnum] + "_apion" + Par.magf + ".txt");
		akaon.raw_yield_outfile.open("../data/RawYields/" + 
			Par.run_name + "/" + Par.Det.name + "/" + Par.CType.cmin_name[Par.cnum] + 
			"_" + Par.CType.cmax_name[Par.cnum] + "_akaon" + Par.magf + ".txt");
		aproton.raw_yield_outfile.open("../data/RawYields/" + 
			Par.run_name + "/" + Par.Det.name + "/" + Par.CType.cmin_name[Par.cnum] + 
			"_" + Par.CType.cmax_name[Par.cnum] + "_aproton" + Par.magf + ".txt");
	}	

	void SetHistStyle(TH1F *hist)
	{
		hist->SetLineColor(kBlack);
		hist->SetMarkerColor(kBlack);
		hist->SetMarkerStyle(8);
		hist->SetMarkerSize(0.5);
		
		hist->SetTitleSize(0.08, "X");
		hist->GetXaxis()->SetLabelSize(0.08);
		hist->GetYaxis()->SetLabelSize(0.08);

		hist->GetXaxis()->SetTitleOffset(0.85);
	}
	
	void PerformM2Fit(TH1F *mass1d_pos, TH1F *mass1d_neg, const double ptmin, const double ptmax, const double nevents, bool do_write_results = true)
	{
		const double pt = (ptmax + ptmin)/2.;
		
		mass1d_pos->GetXaxis()->SetTitle("m^{2}, GeV^{2}/c^{4}");
		mass1d_neg->GetXaxis()->SetTitle("m^{2}, GeV^{2}/c^{4}");

		SetHistStyle(mass1d_pos);
		SetHistStyle(mass1d_neg);

		mass1d_pos->SetTitle("");
		mass1d_neg->SetTitle("");

		TLatex tltext;
		TLine line_pi, line_k, line_p;
		
		line_pi.SetLineColor(Par.pion.color-3);
		line_k.SetLineColor(Par.kaon.color-3);
		line_p.SetLineColor(Par.proton.color-3);
		
		line_pi.SetLineWidth(2);
		line_k.SetLineWidth(2);
		line_p.SetLineWidth(2);

		line_pi.SetLineStyle(2);
		line_k.SetLineStyle(2);
		line_p.SetLineStyle(2);

		tltext.SetTextFont(52);
		tltext.SetTextSize(0.05);

		TText t;
		t.SetTextFont(43);
		t.SetTextSize(40);
		t.SetTextAngle(90);
		
		TCanvas canv("canv", "canv", 900, 900);
		
		canv.Divide(1, 2);
	
		t.DrawText(0.04, 0.65, "charge = +1");
		t.DrawText(0.04, 0.15, "charge = -1");

		tltext.DrawLatex(0.5, 0.9, (DtoStr(ptmin, 1) + " < p_{T} < " + DtoStr(ptmax, 1) + ", GeV/c").c_str());
		tltext.DrawLatex(0.5, 0.4, (DtoStr(ptmin, 1) + " < p_{T} < " + DtoStr(ptmax, 1) + ", GeV/c").c_str());
		
		canv.cd(1);	
		gPad->SetPad(0.02, 0.51, 1., 1.);

		PerformSingleM2Fit(pt, mass1d_pos, 
			Par.pion, Par.kaon, Par.proton, Par.kaon.ptmax, 0);
		
		if (Par.draw_log) gPad->SetLogy();
		gPad->SetBottomMargin(0.17);

		mass1d_pos->GetXaxis()->SetRange(
			mass1d_pos->GetXaxis()->FindBin(Minimum(Par.pion.mean - 
			Par.draw_srange*Par.def_pion_sigma[0]->Eval(pt), 0.005)), 
			mass1d_pos->GetXaxis()->FindBin(Minimum(Par.proton.mean + 
			Par.draw_srange*Par.def_proton_sigma[0]->Eval(pt), 1.5)));
		
		mass1d_pos->GetYaxis()->SetRange(0., mass1d_pos->GetMaximum()*1.1);

		((TH1F *) mass1d_pos->Clone())->Draw("E");
		((TH1F *) mass1d_pos->Clone())->Draw("SAME AXIS X+ Y+");
		
		if (ptmin >= Par.pion.ptmin && ptmax <= Par.pion.ptmax)
		{
			line_pi.DrawLine(fit_pion.mean - fit_pion.sigma*Par.yield_extr_srange, 0, 
				fit_pion.mean - fit_pion.sigma*Par.yield_extr_srange, mass1d_pos->GetMaximum());
			line_pi.DrawLine(fit_pion.mean + fit_pion.sigma*Par.yield_extr_srange, 0, 
				fit_pion.mean + fit_pion.sigma*Par.yield_extr_srange, mass1d_pos->GetMaximum());
		}
		if (ptmin >= Par.kaon.ptmin && ptmax <= Par.kaon.ptmax)
		{
			line_k.DrawLine(fit_kaon.mean - fit_kaon.sigma*Par.yield_extr_srange, 0, 
				fit_kaon.mean - fit_kaon.sigma*Par.yield_extr_srange, mass1d_pos->GetMaximum());
			line_k.DrawLine(fit_kaon.mean + fit_kaon.sigma*Par.yield_extr_srange, 0, 
				fit_kaon.mean + fit_kaon.sigma*Par.yield_extr_srange, mass1d_pos->GetMaximum());
		}
		if (ptmin >= Par.proton.ptmin && ptmax <= Par.proton.ptmax)
		{
			line_p.DrawLine(fit_proton.mean - fit_proton.sigma*Par.yield_extr_srange, 0, 
				fit_proton.mean - fit_proton.sigma*Par.yield_extr_srange, mass1d_pos->GetMaximum());
			line_p.DrawLine(fit_proton.mean + fit_proton.sigma*Par.yield_extr_srange, 0, 
				fit_proton.mean + fit_proton.sigma*Par.yield_extr_srange, mass1d_pos->GetMaximum());
		}
		
		if (do_write_results)
		{
			if (pt >= Par.pion.ptmin && pt <= Par.pion.ptmax) 
			{
				AddMSPoint(pion, pt, pi_gaus->GetParameter(1), pi_gaus->GetParameter(2));
			}
			
			if (pt >= Par.kaon.ptmin && pt <= Par.kaon.ptmax) 
			{
				AddMSPoint(kaon, pt, k_gaus->GetParameter(1), k_gaus->GetParameter(2));
			}
			
			if (pt >= Par.proton.ptmin && pt <= Par.proton.ptmax) 
			{
				AddMSPoint(proton, pt, p_gaus->GetParameter(1), p_gaus->GetParameter(2));
			}
		}
		
		pi_gaus->SetLineColor(Par.apion.color-3);
		k_gaus->SetLineColor(Par.akaon.color-3);
		p_gaus->SetLineColor(Par.aproton.color-3);
	
		canv.cd(2);
		gPad->SetPad(0.02, 0., 1., 0.5);

		PerformSingleM2Fit(pt, mass1d_neg, 
			Par.apion, Par.akaon, Par.aproton, Par.kaon.ptmax, 1);
		if (Par.draw_log) gPad->SetLogy();
		gPad->SetBottomMargin(0.17);
		
		mass1d_neg->GetXaxis()->SetRange(
			mass1d_neg->GetXaxis()->FindBin(Minimum(Par.pion.mean - 
			Par.draw_srange*Par.def_pion_sigma[1]->Eval(pt), 0.005)), 
			mass1d_neg->GetXaxis()->FindBin(Minimum(Par.proton.mean + 
			Par.draw_srange*Par.def_proton_sigma[1]->Eval(pt), 1.5)));

		((TH1F *) mass1d_neg->Clone())->Draw("E");
		((TH1F *) mass1d_neg->Clone())->Draw("SAME AXIS X+ Y+");

		if (ptmin >= Par.apion.ptmin && ptmax <= Par.apion.ptmax)
		{
			line_pi.DrawLine(fit_pion.mean - fit_pion.sigma*Par.yield_extr_srange, 0, 
				fit_pion.mean - fit_pion.sigma*Par.yield_extr_srange, mass1d_neg->GetMaximum());
			line_pi.DrawLine(fit_pion.mean + fit_pion.sigma*Par.yield_extr_srange, 0, 
				fit_pion.mean + fit_pion.sigma*Par.yield_extr_srange, mass1d_neg->GetMaximum());
		}
		if (ptmin >= Par.akaon.ptmin && ptmax <= Par.akaon.ptmax)
		{
			line_k.DrawLine(fit_kaon.mean - fit_kaon.sigma*Par.yield_extr_srange, 0, 
				fit_kaon.mean - fit_kaon.sigma*Par.yield_extr_srange, mass1d_neg->GetMaximum());
			line_k.DrawLine(fit_kaon.mean + fit_kaon.sigma*Par.yield_extr_srange, 0, 
				fit_kaon.mean + fit_kaon.sigma*Par.yield_extr_srange, mass1d_neg->GetMaximum());
		}
		if (ptmin >= Par.aproton.ptmin && ptmax <= Par.aproton.ptmax)
		{
			line_p.DrawLine(fit_proton.mean - fit_proton.sigma*Par.yield_extr_srange, 0, 
				fit_proton.mean - fit_proton.sigma*Par.yield_extr_srange, mass1d_neg->GetMaximum());
			line_p.DrawLine(fit_proton.mean + fit_proton.sigma*Par.yield_extr_srange, 0, 
				fit_proton.mean + fit_proton.sigma*Par.yield_extr_srange, mass1d_neg->GetMaximum());
		}
		
		PrintCanvas(&canv, 
			"../output/m2/" + Par.run_name + "/" + Par.Det.name + "_" +
			Par.CType.cmin_name[Par.cnum] + "_" + Par.CType.cmax_name[Par.cnum] + "/" +
			DtoStr(pt, 2) + Par.magf + ".png");

		tltext.SetTextSize(0.1);
		
		TCanvas pos_canv("pos", "canv", 900, 350);
		if (Par.draw_log) gPad->SetLogy();
		gPad->SetBottomMargin(0.15);
			
		((TH1F *) mass1d_pos->Clone())->Draw("E");
		((TH1F *) mass1d_pos->Clone())->Draw("SAME AXIS X+ Y+");
		tltext.DrawLatexNDC(0.5, 0.8, (DtoStr(ptmin, 1) + " < p_{T} < " + DtoStr(ptmax, 1) + ", GeV/c").c_str());
		
		PrintCanvas(&pos_canv, 
			"../output/m2/" + Par.run_name + "/" + Par.Det.name + "_" +
			Par.CType.cmin_name[Par.cnum] + "_" + Par.CType.cmax_name[Par.cnum] + "/pos_" +
			DtoStr(pt, 2) + Par.magf);

		TCanvas neg_canv("neg", "canv", 900, 350);
		if (Par.draw_log) gPad->SetLogy();
		gPad->SetBottomMargin(0.15);
		
		((TH1F *) mass1d_neg->Clone())->Draw("E");
		((TH1F *) mass1d_neg->Clone())->Draw("SAME AXIS X+ Y+");
		tltext.DrawLatexNDC(0.5, 0.8, (DtoStr(ptmin, 1) + " < p_{T} < " + DtoStr(ptmax, 1) + ", GeV/c").c_str());
		
		PrintCanvas(&neg_canv, 
			"../output/m2/" + Par.run_name + "/" + Par.Det.name + "_" +
			Par.CType.cmin_name[Par.cnum] + "_" + Par.CType.cmax_name[Par.cnum] + "/neg_" +
			DtoStr(pt, 2) + Par.magf);
		
		if (do_write_results)
		{
			if (pt >= Par.apion.ptmin && pt <= Par.apion.ptmax) 
			{
				AddMSPoint(apion, -pt, pi_gaus->GetParameter(1), pi_gaus->GetParameter(2));
			}

			if (pt >= Par.akaon.ptmin && pt <= Par.akaon.ptmax) 
			{
				AddMSPoint(akaon, -pt, k_gaus->GetParameter(1), k_gaus->GetParameter(2));
			}

			if (pt >= Par.aproton.ptmin && pt <= Par.aproton.ptmax) 
			{
				AddMSPoint(aproton, -pt, p_gaus->GetParameter(1), p_gaus->GetParameter(2));
			}
		}
	}
	
	void PerformMSFit()
	{
		PerformSingleMSFit(pion, 0);
		PerformSingleMSFit(kaon, 1);
		PerformSingleMSFit(proton, 2);

		PerformSingleMSFit(apion, 3);
		PerformSingleMSFit(akaon, 4);
		PerformSingleMSFit(aproton, 5);
	}	

	void DrawMS(std::string output_name)
	{
		TCanvas canv(Par.Det.name.c_str(), Par.Det.name.c_str(), 800, 700);

		gPad->SetLeftMargin(0.11);
		gPad->SetBottomMargin(0.11);

		double min_y = -0.2;
		double max_y = 1.5;
		
		TH2F *range_hist = new TH2F(Par.Det.name.c_str(), Par.Det.name.c_str(), 
			100, apart_min_x, part_max_x, 100, min_y, max_y);

		TLatex tltext;

		double tsize = 0.04;
		tltext.SetTextFont(52);
		tltext.SetTextSize(0.05);
		
		range_hist->GetXaxis()->SetTitle("p_{T} #times charge, GeV/c");
		range_hist->GetYaxis()->SetTitle("M^{2}, (GeV/c^{2})^{2}");
		
		range_hist->GetYaxis()->SetTitleOffset(1.);
		range_hist->GetYaxis()->SetTitleOffset(1.);

		range_hist->GetXaxis()->SetLabelSize(0.05);
		range_hist->GetYaxis()->SetLabelSize(0.05);
		
		range_hist->SetTitleSize(0.05, "X");
		range_hist->SetTitleSize(0.05, "Y");
		
		range_hist->Draw("AXIS");
		range_hist->Draw("SAME AXIS X+ Y+");

		SetShades();

		pion.ms_shade->Draw("FL");
		kaon.ms_shade->Draw("FL");
		proton.ms_shade->Draw("FL");

		apion.ms_shade->Draw("FL");
		akaon.ms_shade->Draw("FL");
		aproton.ms_shade->Draw("FL");

		DrawSingleMS(pion);
		DrawSingleMS(kaon);
		DrawSingleMS(proton);

		DrawSingleMS(apion);
		DrawSingleMS(akaon);
		DrawSingleMS(aproton);

		const double tltext_x_pos = 0.05;
		const double vert_shift = 0.5;

		tltext.DrawLatex(tltext_x_pos, 
			pion.means_fit->Eval(Par.pion.ptmin) + 
			pion.sigmas_fit->Eval(Par.pion.ptmin)*vert_shift - tsize, "#pi^{+}");
		tltext.DrawLatex(tltext_x_pos, 
			kaon.means_fit->Eval(Par.kaon.ptmin) + 
			kaon.sigmas_fit->Eval(Par.kaon.ptmin)*vert_shift - tsize, "K^{+}");
		tltext.DrawLatex(tltext_x_pos, 
			proton.means_fit->Eval(Par.kaon.ptmin) + 
			proton.sigmas_fit->Eval(Par.proton.ptmin)*vert_shift - tsize, "p^{+}");

		tltext.DrawLatex(-tltext_x_pos-tsize*5., 
			apion.means_fit->Eval(-Par.apion.ptmin) - 
			apion.sigmas_fit->Eval(-Par.apion.ptmin)*vert_shift - tsize, "#pi^{-}");
		tltext.DrawLatex(-tltext_x_pos-tsize*5., 
			akaon.means_fit->Eval(-Par.akaon.ptmin) - 
			akaon.sigmas_fit->Eval(-Par.akaon.ptmin)*vert_shift - tsize, "K^{-}");
		tltext.DrawLatex(-tltext_x_pos-tsize*5., 
			aproton.means_fit->Eval(-Par.aproton.ptmin) - 
			aproton.sigmas_fit->Eval(-Par.aproton.ptmin)*vert_shift - tsize, "p^{-}");

		PrintCanvas(&canv, output_name);

		delete range_hist;
	}

	void PrintSinglePartPar(FitPart &part, std::string part_name)
	{	
		/*
		PrintFuncPar(part.sigmas_fit, "static const double", 
			Par.Det.hist_name + "_mean_par_" + part_name);
		*/
		PrintFuncPar(part.means_fit, "static const double", 
			Par.Det.hist_name + "_mean_par_" + part_name);

		system(("mkdir -p ../par/MS/" + Par.run_name + "/" + Par.Det.name).c_str());
		std::string outfile_name = "../par/MS/" + Par.run_name + "/" + 
			Par.Det.name + "/" + part_name + ".txt";
		ofstream par_output(outfile_name.c_str());

		for (int i = 0; i < part.sigmas_fit->GetNpar()-1; i++) 
		{
			par_output << part.sigmas_fit->GetParameter(i) << " ";
		}
		par_output << part.sigmas_fit->GetParameter(part.sigmas_fit->GetNpar() - 1);
		par_output.close();
	}

	void PrintMSFunc()
	{
		std::cout << std::endl;
		
		PrintSeparator("Means and sigmas approximation functions");
		
		Print("static std::string MeanFuncName = \"" + Par.mean_fit + "\";");
		Print("static std::string SigmaFuncName = \"" + Par.sigma_fit + "\";");
		
		PrintSeparator("End of means and sigmas approximation functions");

		PrintSeparator("Means and sigmas approximation parameters in " + Par.Det.name + 
			" for centr " + Par.CType.cmin_name[Par.cnum] + "-" + Par.CType.cmax_name[Par.cnum]);
		
		PrintSinglePartPar(pion, "pion");
		PrintSinglePartPar(kaon, "kaon");
		PrintSinglePartPar(proton, "proton");

		PrintSinglePartPar(apion, "apion");
		PrintSinglePartPar(akaon, "akaon");
		PrintSinglePartPar(aproton, "aproton");

		PrintSeparator("End of means and sigmas approximation parameters in " + 
			Par.Det.name + " for centr " + 
			Par.CType.cmin_name[Par.cnum] + "-" + Par.CType.cmax_name[Par.cnum]);

		Print("pion parameters:", pion.sigmas_fit->GetParameter(2), 
				pion.sigmas_fit->GetParameter(3), pion.sigmas_fit->GetParameter(4));
		Print("apion parameters:", apion.sigmas_fit->GetParameter(2), 
				apion.sigmas_fit->GetParameter(3), apion.sigmas_fit->GetParameter(4));
		Print("kaon parameters:", kaon.sigmas_fit->GetParameter(2), 
				kaon.sigmas_fit->GetParameter(3), kaon.sigmas_fit->GetParameter(4));
		Print("akaon parameters:", akaon.sigmas_fit->GetParameter(2), 
				akaon.sigmas_fit->GetParameter(3), akaon.sigmas_fit->GetParameter(4));
		Print("proton parameters:", proton.sigmas_fit->GetParameter(2), 
				proton.sigmas_fit->GetParameter(3), proton.sigmas_fit->GetParameter(4));
		Print("aproton parameters:", aproton.sigmas_fit->GetParameter(2), 
				aproton.sigmas_fit->GetParameter(3), aproton.sigmas_fit->GetParameter(4));
		
		Print(Par.Det.name + "_sigma_alpha_" + 
			Par.CType.cmin_name[Par.cnum] + "_" + Par.CType.cmax_name[Par.cnum] + 
			" = " + DtoStr(Average(
			pion.sigmas_fit->GetParameter(2), 
			kaon.sigmas_fit->GetParameter(2), 
			proton.sigmas_fit->GetParameter(2), 
			apion.sigmas_fit->GetParameter(2), 
			akaon.sigmas_fit->GetParameter(2), 
			aproton.sigmas_fit->GetParameter(2)), 3) + ";");

		Print(Par.Det.name + "_sigma_ms_" + 
			Par.CType.cmin_name[Par.cnum] + "_" + Par.CType.cmax_name[Par.cnum] + 
			" = " + DtoStr(Average(
			pion.sigmas_fit->GetParameter(3), 
			kaon.sigmas_fit->GetParameter(3), 
			proton.sigmas_fit->GetParameter(3), 
			apion.sigmas_fit->GetParameter(3), 
			akaon.sigmas_fit->GetParameter(3), 
			aproton.sigmas_fit->GetParameter(3)), 3) + ";");

		Print(Par.Det.name + "_sigma_t_" + 
			Par.CType.cmin_name[Par.cnum] + "_" + Par.CType.cmax_name[Par.cnum] + 
			" = " + DtoStr(Average(
			pion.sigmas_fit->GetParameter(4), 
			kaon.sigmas_fit->GetParameter(4), 
			proton.sigmas_fit->GetParameter(4), 
			apion.sigmas_fit->GetParameter(4), 
			akaon.sigmas_fit->GetParameter(4), 
			aproton.sigmas_fit->GetParameter(4)), 3) + ";");			
	}
};
