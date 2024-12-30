// $SOURCE$
//------------------------------------------------------------------------------------------------
//                SigmRes functions declarations and realisations
//------------------------------------------------------------------------------------------------
// SigmRes
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic set of functions for calculating sigmalized residuals dz and dphi of tracks
 **/
//------------------------------------------------------------------------------------------------

#ifndef SDPHI_SDZ_H
#define SDPHI_SDZ_H

#include <cmath>

namespace DummySigmRes
{
   double SDPhi(const double dphi, const double pT)
   {
      return fabs(dphi)/0.03;
   }
   
   double SDZ(const double dz, const double pT)
   {
      return fabs(dz)/5.;
   }
}

namespace SimRun7AuAu200SigmRes
{
   double GetPC2SDPhi(const double dphi, const double pT, const int charge) 
   {
      return DummySigmRes::SDPhi(dphi, pT);
   }
   
   double GetPC2SDZ(const double dz, const double pT, const int charge) 
   {
      return DummySigmRes::SDZ(dz, pT);
   }
   
   double GetPC3SDPhi(const double phi, const double dphi, const double pT, const int charge) 
   {
      return DummySigmRes::SDPhi(dphi, pT);
   }
   
   double GetPC3SDZ(const double phi, const double dz, const double pT, const int charge) 
   {
      return DummySigmRes::SDZ(dz, pT);
   }
   
   double GetTOFeSDPhi(const double dphi, const double pT, const int charge) 
   {
      return DummySigmRes::SDPhi(dphi, pT);
   }
   
   double GetTOFeSDZ(const double dz, const double pT, const int charge) 
   {
      return DummySigmRes::SDZ(dz, pT);
   }
   
   double GetTOFwSDPhi(const double dphi, const double pT, const int charge) 
   {
      return DummySigmRes::SDPhi(dphi, pT);
   }
   
   double GetTOFwSDZ(const double dz, const double pT, const int charge) 
   {
      return DummySigmRes::SDZ(dz, pT);
   }
   
   double GetEMCSDPhi(const double phi, const double dphi, const double pT, 
                      const int charge, const int sect) 
   {
      return DummySigmRes::SDPhi(dphi, pT);
   }
   
   double GetEMCSDZ(const double phi, const double dz, const double pT, 
                    const int charge, const int sect) 
   {
      return DummySigmRes::SDZ(dz, pT);
   }
}

namespace SimRun14HeAu200SigmRes
{
   double GetPC2SDPhi(const double dphi, const double pT, const int charge) 
   {
      return DummySigmRes::SDPhi(dphi, pT);
   }
   
   double GetPC2SDZ(const double dz, const double pT, const int charge) 
   {
      return DummySigmRes::SDZ(dz, pT);
   }
   
   double GetPC3SDPhi(const double phi, const double dphi, const double pT, const int charge) 
   {
      return DummySigmRes::SDPhi(dphi, pT);
   }
   
   double GetPC3SDZ(const double phi, const double dz, const double pT, const int charge) 
   {
      return DummySigmRes::SDZ(dz, pT);
   }
   
   double GetTOFeSDPhi(const double dphi, const double pT, const int charge) 
   {
      return DummySigmRes::SDPhi(dphi, pT);
   }
   
   double GetTOFeSDZ(const double dz, const double pT, const int charge) 
   {
      return DummySigmRes::SDZ(dz, pT);
   }
   
   double GetTOFwSDPhi(const double dphi, const double pT, const int charge) 
   {
      return DummySigmRes::SDPhi(dphi, pT);
   }
   
   double GetTOFwSDZ(const double dz, const double pT, const int charge) 
   {
      return DummySigmRes::SDZ(dz, pT);
   }
   
   double GetEMCSDPhi(const double phi, const double dphi, const double pT, 
                      const int charge, const int sect) 
   {
      return DummySigmRes::SDPhi(dphi, pT);
   }
   
   double GetEMCSDZ(const double phi, const double dz, const double pT, 
                    const int charge, const int sect) 
   {
      return DummySigmRes::SDZ(dz, pT);
   }
}

#endif /* SDPHI_SDZ_H */
