#include "ParPair.h"
#include "Tools.h"
#include "OutputTool.h"

double GetEMCalEid(const double pt, const int part_id, const int charge, const double phi, const int sector)
{
	if (pt > 1.2) return 0.;
	
	double mean_pion, mean_kaon, mean_proton;
	double sigma_pion, sigma_kaon, sigma_proton;
	
	if (phi >= 1.5)
	{
		if (charge > 0)
		{
			mean_pion = GetM2Mean(pt, m2_emcale_mean_par_pion[sector-2]);
			mean_kaon = GetM2Mean(pt, m2_emcale_mean_par_kaon[sector-2]);
			mean_proton = GetM2Mean(pt, m2_emcale_mean_par_proton[sector-2]);
			
			sigma_pion = GetM2Sigma(pt, mean_pion, m2_emcale_sigma_par[sector-2]);
			sigma_kaon = GetM2Sigma(pt, mean_kaon, m2_emcale_sigma_par[sector-2]);
			sigma_proton = GetM2Sigma(pt, mean_proton, m2_emcale_sigma_par[sector-2]);
		}
		else 
		{
			mean_pion = GetM2Mean(pt, m2_emcale_mean_par_apion[sector]);
			mean_kaon = GetM2Mean(pt, m2_emcale_mean_par_akaon[sector]);
			mean_proton = GetM2Mean(pt, m2_emcale_mean_par_aproton[sector]);
			
			sigma_pion = GetM2Sigma(pt, mean_pion, m2_emcale_sigma_par[sector]);
			sigma_kaon = GetM2Sigma(pt, mean_kaon, m2_emcale_sigma_par[sector]);
			sigma_proton = GetM2Sigma(pt, mean_proton, m2_emcale_sigma_par[sector]);
		}
	}
	else
	{
		if (charge > 0)
		{
			mean_pion = GetM2Mean(pt, m2_emcalw_mean_par_pion[sector]);
			mean_kaon = GetM2Mean(pt, m2_emcalw_mean_par_kaon[sector]);
			mean_proton = GetM2Mean(pt, m2_emcalw_mean_par_proton[sector]);
			
			sigma_pion = GetM2Sigma(pt, mean_pion, m2_emcalw_sigma_par[sector]);
			sigma_kaon = GetM2Sigma(pt, mean_kaon, m2_emcalw_sigma_par[sector]);
			sigma_proton = GetM2Sigma(pt, mean_proton, m2_emcalw_sigma_par[sector]);
		}
		else 
		{
			mean_pion = GetM2Mean(pt, m2_emcalw_mean_par_apion[sector]);
			mean_kaon = GetM2Mean(pt, m2_emcalw_mean_par_akaon[sector]);
			mean_proton = GetM2Mean(pt, m2_emcalw_mean_par_aproton[sector]);
			
			sigma_pion = GetM2Sigma(pt, mean_pion, m2_emcalw_sigma_par[sector]);
			sigma_kaon = GetM2Sigma(pt, mean_kaon, m2_emcalw_sigma_par[sector]);
			sigma_proton = GetM2Sigma(pt, mean_proton, m2_emcalw_sigma_par[sector]);
		}
	}
	
	const double sigma_veto = 3.;
	const double sigma_range = 1.;
		
	double weight = 0;

	switch (part_id)
	{
		case 211:
		{
			double low_range = mean_pion - sigma_range*sigma_pion;
			double upp_range = mean_pion + sigma_range*sigma_pion;
			
			double upp_veto = mean_kaon - sigma_veto*sigma_kaon;
			if (low_range >= upp_veto) return 0.;

			weight = erf((mean_pion - low_range)/sigma_pion/TMath::Sqrt2())/2. +
				erf((Minimum(upp_range, upp_veto) - mean_pion)/sigma_pion/TMath::Sqrt2())/2.;
		}
		case 321:
		{
			if (charge == 1 && pt < 0.4) return 0.;
			
			double low_range = mean_kaon - sigma_range*sigma_kaon;
			double upp_range = mean_kaon + sigma_range*sigma_kaon;
			
			double low_veto = mean_pion + sigma_veto*sigma_pion; 
			double upp_veto = mean_proton - sigma_veto*sigma_proton; 
			
			if (low_range >= upp_veto || upp_range <= low_veto || upp_veto <= low_veto) return 0.;
			
			weight = erf((mean_kaon - Maximum(low_range, low_veto))/sigma_kaon/TMath::Sqrt2())/2. +
				erf((Minimum(upp_range, upp_veto) - mean_kaon)/sigma_kaon/TMath::Sqrt2())/2.;
		}	
		case 2212:
		{
			if (charge == 1 && pt < 0.6) return 0.;
			else if (charge == -1 && pt < 0.4) return 0.;
			
			double low_range = mean_proton - sigma_range*sigma_proton;
			double upp_range = mean_proton + sigma_range*sigma_proton;
			
			double low_veto = mean_kaon + sigma_veto*sigma_kaon;
			
			if (upp_range <= low_veto) return 0.;

			weight = erf((mean_proton - Maximum(low_range, low_veto))/sigma_proton/TMath::Sqrt2())/2. +
				erf((upp_range - mean_proton)/sigma_proton/TMath::Sqrt2())/2.;
		}
	}	
	return weight/erf(sigma_range/TMath::Sqrt2());
}

int GetTOFePID(const double pt, const int charge, const double m2)
{
	if (pt > 4.5 || m2 < -998.) return PartId.noPID;
	
	double mean_pion, mean_kaon, mean_proton;
	double sigma_pion, sigma_kaon, sigma_proton;
	
	if (charge > 0)
	{
		mean_pion = GetM2Mean(pt, m2_tofe_mean_par_pion);
		mean_kaon = GetM2Mean(pt, m2_tofe_mean_par_kaon);
		mean_proton = GetM2Mean(pt, m2_tofe_mean_par_proton);
		
		sigma_pion = GetM2Sigma(pt, mean_pion, m2_tofe_sigma_par);
		sigma_kaon = GetM2Sigma(pt, mean_kaon, m2_tofe_sigma_par);
		sigma_proton = GetM2Sigma(pt, mean_proton, m2_tofe_sigma_par);
	}
	else 
	{
		mean_pion = GetM2Mean(pt, m2_tofe_mean_par_apion);
		mean_kaon = GetM2Mean(pt, m2_tofe_mean_par_akaon);
		mean_proton = GetM2Mean(pt, m2_tofe_mean_par_aproton);
		
		sigma_pion = GetM2Sigma(pt, mean_pion, m2_tofe_sigma_par);
		sigma_kaon = GetM2Sigma(pt, mean_kaon, m2_tofe_sigma_par);
		sigma_proton = GetM2Sigma(pt, mean_proton, m2_tofe_sigma_par);
	}
	
	const double sigma_veto = 3.;
	const double sigma_range = 2.;
		
	const double pion_hi_veto = mean_pion + sigma_veto*sigma_pion;
	const double kaon_lo_veto = mean_kaon - sigma_veto*sigma_kaon;
	const double kaon_hi_veto = mean_kaon + sigma_veto*sigma_kaon;
	const double prot_lo_veto = mean_proton - sigma_veto*sigma_proton;
		
	if ((m2 > mean_pion - sigma_range*sigma_pion && 
		m2 < mean_pion + sigma_range*sigma_pion) && 
		m2 < kaon_lo_veto && m2 < prot_lo_veto) return PartId.pion;
		
	if (pt < 2. && (m2 > mean_kaon - sigma_range*sigma_kaon && 
		m2 < mean_kaon + sigma_range*sigma_kaon) && 
		m2 > pion_hi_veto && m2 < prot_lo_veto) return PartId.kaon;

	if ((m2 > mean_proton - sigma_range*sigma_proton && 
		m2 < mean_proton + sigma_range*sigma_proton) && 
		m2 > kaon_hi_veto && m2 > pion_hi_veto) return PartId.proton;
			
	return PartId.noPID;
}

int GetTOFwPID(const double pt, const int charge, const double m2)
{
	if (pt > 6. || m2 < -998.) return PartId.noPID;
	
	double mean_pion, mean_kaon, mean_proton;
	double sigma_pion, sigma_kaon, sigma_proton;
	
	if (charge > 0)
	{
		mean_pion = GetM2Mean(pt, m2_tofw_mean_par_pion);
		mean_kaon = GetM2Mean(pt, m2_tofw_mean_par_kaon);
		mean_proton = GetM2Mean(pt, m2_tofw_mean_par_proton);
		
		sigma_pion = GetM2Sigma(pt, mean_pion, m2_tofw_sigma_par);
		sigma_kaon = GetM2Sigma(pt, mean_kaon, m2_tofw_sigma_par);
		sigma_proton = GetM2Sigma(pt, mean_proton, m2_tofw_sigma_par);
	}
	else 
	{
		mean_pion = GetM2Mean(pt, m2_tofw_mean_par_apion);
		mean_kaon = GetM2Mean(pt, m2_tofw_mean_par_akaon);
		mean_proton = GetM2Mean(pt, m2_tofw_mean_par_aproton);
		
		sigma_pion = GetM2Sigma(pt, mean_pion, m2_tofw_sigma_par);
		sigma_kaon = GetM2Sigma(pt, mean_kaon, m2_tofw_sigma_par);
		sigma_proton = GetM2Sigma(pt, mean_proton, m2_tofw_sigma_par);
	}
	
	const double sigma_veto = 3.;
	const double sigma_range = 2.;
		
	const double pion_hi_veto = mean_pion + sigma_veto*sigma_pion;
	const double kaon_lo_veto = mean_kaon - sigma_veto*sigma_kaon;
	const double kaon_hi_veto = mean_kaon + sigma_veto*sigma_kaon;
	const double prot_lo_veto = mean_proton - sigma_veto*sigma_proton;
		
	if ((m2 > mean_pion - sigma_range*sigma_pion && 
		m2 < mean_pion + sigma_range*sigma_pion) && 
		m2 < kaon_lo_veto && m2 < prot_lo_veto) return PartId.pion;
		
	if (pt < 4. && (m2 > mean_kaon - sigma_range*sigma_kaon && 
		m2 < mean_kaon + sigma_range*sigma_kaon) && 
		m2 > pion_hi_veto && m2 < prot_lo_veto) return PartId.kaon;

	if ((m2 > mean_proton - sigma_range*sigma_proton && 
		m2 < mean_proton + sigma_range*sigma_proton) && 
		m2 > kaon_hi_veto && m2 > pion_hi_veto) return PartId.proton;
			
	return PartId.noPID;
}
