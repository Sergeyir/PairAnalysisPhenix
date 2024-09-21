// $HEADER$
//------------------------------------------------------------------------------------------------
//                                 PTrackFun declaration
//------------------------------------------------------------------------------------------------
// PTrackFun - pair track functions
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic functions used in analysis for particle tracks pairs
 **/
//------------------------------------------------------------------------------------------------

#ifndef P_TRACK_FUN_CPP
#define P_TRACK_FUN_CPP

#include "../include/PTrackFun.hpp"

bool IsTOF2PID(idContainer *id1, idContainer *id2)
{
   if (id1->tof == id1->origId && id2->tof == id2->origId) return true;
   
   return false;
}

bool IsEMC2PID(idContainer *id1, idContainer *id2)
{
   if (id1->emc == id1->origId && id2->emc == id2->origId) return true;
   
   return false;
}

bool IsEMCnoPID(idContainer *id1, idContainer *id2)
{
   if (id1->emc != PartId.junk && id2->emc != PartId.junk) return true;
   
   return false;
}

bool Is1TOFandEMC2PID(idContainer *id1, idContainer *id2)
{
   if 
   (
      IsTOF2PID(id1, id2) ||
      (id1->emc == id1->origId && id2->tof == id2->origId) || 
      (id1->tof == id1->origId && id2->emc == id2->origId)
   ) return true;

   return false;
}

bool Is2PID(idContainer *id1, idContainer *id2)
{
   if (IsEMC2PID(id1, id2) || Is1TOFandEMC2PID(id1, id2)) return true;
   
   return false;
}

bool Is1PID(idContainer *id1, idContainer *id2)
{
   if 
   (
      (id1->tof == id1->origId && 
      (
         id2->tof != PartId.junk || 
         id2->emc != PartId.junk || 
         id2->pc2 != PartId.junk ||
         id2->pc3 != PartId.junk)
      ) || 
      (id2->tof == id2->origId && 
      (
         id1->tof != PartId.junk || 
         id1->emc != PartId.junk || 
         id1->pc2 != PartId.junk ||
         id1->pc3 != PartId.junk)
      ) 
   ) return true;
   
   return false;
}

bool IsnoPID(idContainer *id1, idContainer *id2)
{
   return true;
}

//mom = momentum
double GetPairPT(const double *mom1, const double *mom2)
{
   const double pT = sqrt((mom1[0]+mom2[0])*(mom1[0]+mom2[0]) + 
                          (mom1[1]+mom2[1])*(mom1[1]+mom2[1]));
   return pT;
}

//mom = momentum
double GetMass(const double *mom1, const double *mom2, const double m1, const double m2)
{
   const double mom1Squared = mom1[0]*mom1[0] + mom1[1]*mom1[1] + mom1[2]*mom1[2];
   const double mom2Squared = mom2[0]*mom2[0] + mom2[1]*mom2[1] + mom2[2]*mom2[2];
   
   const double energy1 = sqrt(pow(m1, 2) + mom1Squared);
   const double energy2 = sqrt(pow(m2, 2) + mom2Squared);
   const double pairEnergy = energy1 + energy2;

   double pairMom[3];
   
   pairMom[0] = mom1[0] + mom2[0];
   pairMom[1] = mom1[1] + mom2[1];
   pairMom[2] = mom1[2] + mom2[2];
   
   const double pairMass = sqrt(pairEnergy*pairEnergy - pairMom[0]*pairMom[0] - 
                                pairMom[1]*pairMom[1] - pairMom[2]*pairMom[2]);
   
   return pairMass;
}

bool IsOneArmCut(const double phi1, const double phi2)
{
   if ( !((phi1 < 1.5 && phi2 < 1.5) || (phi1 > 1.5 && phi2 > 1.5)) ) return true;
   return false;
}

bool IsGhostCut(double dZed, double dPhi, double dAlpha)
{
   //pc1
   if (fabs(dZed) < 6.0 && fabs(dPhi - (0.13*dAlpha)) < 0.015) return true;
   //x1x2_1
   if (fabs(dPhi - (0.04*dAlpha)) < 0.015) return true;
   //x1x2_2
   if (fabs(dPhi - (-0.065*dAlpha)) < 0.015) return true;
   return false;
}

bool noCut()
{
   return true;
}

#endif /* P_TRACK_FUN_CPP */
