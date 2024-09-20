#include "ParAnalyzeResonance.hpp"
#include "IOTools.hpp"

double GetEMCalId(const double pt, const int id, const int charge, const double phi, 
                  const int sector)
{
   const double SIGMA_VETO = 3.;
   const double SIGMA_RANGE = 2.;
      
   if (pt > 1.2) return 0.;
   
   double m2MeanPion, m2MeanKaon, m2MeanProton;
   double m2SigmaPion, m2SigmaKaon, m2SigmaProton;
   
   if (phi >= 1.5)
   {
      if (charge > 0)
      {
         m2MeanPion = GetM2Mean(pt, m2_emcale_mean_par_pion[sector-2]);
         m2MeanKaon = GetM2Mean(pt, m2_emcale_mean_par_kaon[sector-2]);
         m2MeanProton = GetM2Mean(pt, m2_emcale_mean_par_proton[sector-2]);
         
         m2SigmaPion = GetM2Sigma(pt, m2MeanPion, m2_emcale_sigma_par[sector-2]);
         m2SigmaKaon = GetM2Sigma(pt, m2MeanKaon, m2_emcale_sigma_par[sector-2]);
         m2SigmaProton = GetM2Sigma(pt, m2MeanProton, m2_emcale_sigma_par[sector-2]);
      }
      else 
      {
         m2MeanPion = GetM2Mean(pt, m2_emcale_mean_par_apion[sector]);
         m2MeanKaon = GetM2Mean(pt, m2_emcale_mean_par_akaon[sector]);
         m2MeanProton = GetM2Mean(pt, m2_emcale_mean_par_aproton[sector]);
         
         m2SigmaPion = GetM2Sigma(pt, m2MeanPion, m2_emcale_sigma_par[sector]);
         m2SigmaKaon = GetM2Sigma(pt, m2MeanKaon, m2_emcale_sigma_par[sector]);
         m2SigmaProton = GetM2Sigma(pt, m2MeanProton, m2_emcale_sigma_par[sector]);
      }
   }
   else
   {
      if (charge > 0)
      {
         m2MeanPion = GetM2Mean(pt, m2_emcalw_mean_par_pion[sector]);
         m2MeanKaon = GetM2Mean(pt, m2_emcalw_mean_par_kaon[sector]);
         m2MeanProton = GetM2Mean(pt, m2_emcalw_mean_par_proton[sector]);
         
         m2SigmaPion = GetM2Sigma(pt, m2MeanPion, m2_emcalw_sigma_par[sector]);
         m2SigmaKaon = GetM2Sigma(pt, m2MeanKaon, m2_emcalw_sigma_par[sector]);
         m2SigmaProton = GetM2Sigma(pt, m2MeanProton, m2_emcalw_sigma_par[sector]);
      }
      else 
      {
         m2MeanPion = GetM2Mean(pt, m2_emcalw_mean_par_apion[sector]);
         m2MeanKaon = GetM2Mean(pt, m2_emcalw_mean_par_akaon[sector]);
         m2MeanProton = GetM2Mean(pt, m2_emcalw_mean_par_aproton[sector]);
         
         m2SigmaPion = GetM2Sigma(pt, m2MeanPion, m2_emcalw_sigma_par[sector]);
         m2SigmaKaon = GetM2Sigma(pt, m2MeanKaon, m2_emcalw_sigma_par[sector]);
         m2SigmaProton = GetM2Sigma(pt, m2MeanProton, m2_emcalw_sigma_par[sector]);
      }
   }
   
   const double SIGMA_VETO = 3.;
   const double SIGMA_RANGE = 1.;
      
   double weight = 0;

   switch (id)
   {
      case 211:
      {
         double low_range = m2MeanPion - SIGMA_RANGE*m2SigmaPion;
         double upp_range = m2MeanPion + SIGMA_RANGE*m2SigmaPion;
         
         double upp_veto = m2MeanKaon - SIGMA_VETO*m2SigmaKaon;
         if (low_range >= upp_veto) return 0.;

         weight = erf((m2MeanPion - low_range)/m2SigmaPion/TMath::Sqrt2())/2. +
            erf((Minimum(upp_range, upp_veto) - m2MeanPion)/m2SigmaPion/TMath::Sqrt2())/2.;
      }
      case 321:
      {
         if (charge == 1 && pt < 0.4) return 0.;
         
         double low_range = m2MeanKaon - SIGMA_RANGE*m2SigmaKaon;
         double upp_range = m2MeanKaon + SIGMA_RANGE*m2SigmaKaon;
         
         double low_veto = m2MeanPion + SIGMA_VETO*m2SigmaPion; 
         double upp_veto = m2MeanProton - SIGMA_VETO*m2SigmaProton; 
         
         if (low_range >= upp_veto || upp_range <= low_veto || upp_veto <= low_veto) return 0.;
         
         weight = erf((m2MeanKaon - Maximum(low_range, low_veto))/m2SigmaKaon/TMath::Sqrt2())/2. +
            erf((Minimum(upp_range, upp_veto) - m2MeanKaon)/m2SigmaKaon/TMath::Sqrt2())/2.;
      }   
      case 2212:
      {
         if (charge == 1 && pt < 0.6) return 0.;
         else if (charge == -1 && pt < 0.4) return 0.;
         
         double low_range = m2MeanProton - SIGMA_RANGE*m2SigmaProton;
         double upp_range = m2MeanProton + SIGMA_RANGE*m2SigmaProton;
         
         double low_veto = m2MeanKaon + SIGMA_VETO*m2SigmaKaon;
         
         if (upp_range <= low_veto) return 0.;

         weight = erf((m2MeanProton - Maximum(low_range, low_veto))/m2SigmaProton/TMath::Sqrt2())/2. +
            erf((upp_range - m2MeanProton)/m2SigmaProton/TMath::Sqrt2())/2.;
      }
   }   
   return weight/erf(SIGMA_RANGE/TMath::Sqrt2());
}

int GetTOFePID(const double pt, const int charge, const double m2)
{
   if (pt > 4.5 || m2 < -998.) return PartId.noPID;
   
   double m2MeanPion, m2MeanKaon, m2MeanProton;
   double m2SigmaPion, m2SigmaKaon, m2SigmaProton;
   
   if (charge > 0)
   {
      m2MeanPion = GetM2Mean(pt, M2_MEAN_PAR_PION_TOFE);
      m2MeanKaon = GetM2Mean(pt, M2_MEAN_PAR_KAON_TOFE);
      m2MeanProton = GetM2Mean(pt, M2_MEAN_PAR_PROTON_TOFE);
      
      m2SigmaPion = GetM2Sigma(pt, m2MeanPion, M2_SIGMA_PAR_TOFE);
      m2SigmaKaon = GetM2Sigma(pt, m2MeanKaon, M2_SIGMA_PAR_TOFE);
      m2SigmaProton = GetM2Sigma(pt, m2MeanProton, M2_SIGMA_PAR_TOFE);
   }
   else 
   {
      m2MeanPion = GetM2Mean(pt, M2_MEAN_PAR_TOFE_APION);
      m2MeanKaon = GetM2Mean(pt, M2_MEAN_PAR_TOFE_AKAON);
      m2MeanProton = GetM2Mean(pt, M2_MEAN_PAR_TOFE_APROTON);
      
      m2SigmaPion = GetM2Sigma(pt, m2MeanPion, M2_SIGMA_PAR_TOFE);
      m2SigmaKaon = GetM2Sigma(pt, m2MeanKaon, M2_SIGMA_PAR_TOFE);
      m2SigmaProton = GetM2Sigma(pt, m2MeanProton, M2_SIGMA_PAR_TOFE);
   }
   
   const double SIGMA_VETO = 3.;
   const double SIGMA_RANGE = 2.;
      
   const double m2HighVetoPion = m2MeanPion + SIGMA_VETO*m2SigmaPion;
   const double m2LowVetoKaon = m2MeanKaon - SIGMA_VETO*m2SigmaKaon;
   const double m2HighVetoKaon = m2MeanKaon + SIGMA_VETO*m2SigmaKaon;
   const double m2LowVetoProtong = m2MeanProton - SIGMA_VETO*m2SigmaProton;
      
   if ((m2 > m2MeanPion - SIGMA_RANGE*m2SigmaPion && 
      m2 < m2MeanPion + SIGMA_RANGE*m2SigmaPion) && 
      m2 < m2LowVetoKaon && m2 < m2LowVetoProtong) return PartId.pion;
      
   if (pt < 2. && (m2 > m2MeanKaon - SIGMA_RANGE*m2SigmaKaon && 
      m2 < m2MeanKaon + SIGMA_RANGE*m2SigmaKaon) && 
      m2 > m2HighVetoPion && m2 < m2LowVetoProtong) return PartId.kaon;

   if ((m2 > m2MeanProton - SIGMA_RANGE*m2SigmaProton && 
      m2 < m2MeanProton + SIGMA_RANGE*m2SigmaProton) && 
      m2 > m2HighVetoKaon && m2 > m2HighVetoPion) return PartId.proton;
         
   return PartId.noPID;
}

int GetTOFwPID(const double pt, const int charge, const double m2)
{
   if (pt > 6. || m2 < -998.) return PartId.noPID;
   
   double m2MeanPion, m2MeanKaon, m2MeanProton;
   double m2SigmaPion, m2SigmaKaon, m2SigmaProton;
   
   if (charge > 0)
   {
      m2MeanPion = GetM2Mean(pt, M2_MEAN_PAR_TOFW_PION);
      m2MeanKaon = GetM2Mean(pt, M2_MEAN_PAR_TOFW_KAON);
      m2MeanProton = GetM2Mean(pt, M2_MEAN_PAR_TOFW_PROTON);
      
      m2SigmaPion = GetM2Sigma(pt, m2MeanPion, M2_SIGMA_PAR_TOFW);
      m2SigmaKaon = GetM2Sigma(pt, m2MeanKaon, M2_SIGMA_PAR_TOFW);
      m2SigmaProton = GetM2Sigma(pt, m2MeanProton, M2_SIGMA_PAR_TOFW);
   }
   else 
   {
      m2MeanPion = GetM2Mean(pt, M2_MEAN_PAR_TOFW_APION);
      m2MeanKaon = GetM2Mean(pt, M2_MEAN_PAR_TOFW_AKAON);
      m2MeanProton = GetM2Mean(pt, M2_MEAN_PAR_TOFW_APROTON);
      
      m2SigmaPion = GetM2Sigma(pt, m2MeanPion, M2_SIGMA_PAR_TOFW);
      m2SigmaKaon = GetM2Sigma(pt, m2MeanKaon, M2_SIGMA_PAR_TOFW);
      m2SigmaProton = GetM2Sigma(pt, m2MeanProton, M2_SIGMA_PAR_TOFW);
   }
   
   const double m2HighVetoPion = m2MeanPion + SIGMA_VETO*m2SigmaPion;
   const double m2LowVetoKaon = m2MeanKaon - SIGMA_VETO*m2SigmaKaon;
   const double m2HighVetoKaon = m2MeanKaon + SIGMA_VETO*m2SigmaKaon;
   const double m2LowVetoProtong = m2MeanProton - SIGMA_VETO*m2SigmaProton;
      
   if ((m2 > m2MeanPion - SIGMA_RANGE*m2SigmaPion && 
      m2 < m2MeanPion + SIGMA_RANGE*m2SigmaPion) && 
      m2 < m2LowVetoKaon && m2 < m2LowVetoProtong) return PartId.pion;
      
   if (pt < 4. && (m2 > m2MeanKaon - SIGMA_RANGE*m2SigmaKaon && 
      m2 < m2MeanKaon + SIGMA_RANGE*m2SigmaKaon) && 
      m2 > m2HighVetoPion && m2 < m2LowVetoProtong) return PartId.kaon;

   if ((m2 > m2MeanProton - SIGMA_RANGE*m2SigmaProton && 
      m2 < m2MeanProton + SIGMA_RANGE*m2SigmaProton) && 
      m2 > m2HighVetoKaon && m2 > m2HighVetoPion) return PartId.proton;
         
   return PartId.noPID;
}
