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

bool IsTOF2PID(idContainer *Id1, idContainer *Id2)
{
   if (Id1->tof == Id1->orig_id && Id2->tof == Id2->orig_id) return true;
   
   return false;
}

bool IsEMC2PID(idContainer *Id1, idContainer *Id2)
{
   if (Id1->emc == Id1->orig_id && Id2->emc == Id2->orig_id) return true;
   
   return false;
}

bool IsEMCnoPID(idContainer *Id1, idContainer *Id2)
{
   if (Id1->emc != PartId.junk && Id2->emc != PartId.junk) return true;
   
   return false;
}

bool Is1TOFandEMC2PID(idContainer *Id1, idContainer *Id2)
{
   if 
   (
      IsTOF2PID(Id1, Id2) ||
      (Id1->emc == Id1->orig_id && Id2->tof == Id2->orig_id) || 
      (Id1->tof == Id1->orig_id && Id2->emc == Id2->orig_id)
   ) return true;

   return false;
}

bool Is2PID(idContainer *Id1, idContainer *Id2)
{
   if (IsEMC2PID(Id1, Id2) || Is1TOFandEMC2PID(Id1, Id2)) return true;
   
   return false;
}

bool Is1PID(idContainer *Id1, idContainer *Id2)
{
   if 
   (
      (Id1->tof == Id1->orig_id && 
      (
         Id2->tof != PartId.junk || 
         Id2->emc != PartId.junk || 
         Id2->pc2 != PartId.junk ||
         Id2->pc3 != PartId.junk)
      ) || 
      (Id2->tof == Id2->orig_id && 
      (
         Id1->tof != PartId.junk || 
         Id1->emc != PartId.junk || 
         Id1->pc2 != PartId.junk ||
         Id1->pc3 != PartId.junk)
      ) 
   ) return true;
   
   return false;
}

bool IsnoPID(idContainer *Id1, idContainer *Id2)
{
   return true;
}

double GetPairPt(const double *pp1, const double *pp2)
{
   const double pt = sqrt((pp1[0]+pp2[0])*(pp1[0]+pp2[0]) + (pp1[1]+pp2[1])*(pp1[1]+pp2[1]));
   return pt;
}

double GetMass(const double *pp1, const double *pp2, const double m1, const double m2)
{
   const double pm1 = pp1[0]*pp1[0] + pp1[1]*pp1[1] + pp1[2]*pp1[2];
   const double pm2 = pp2[0]*pp2[0] + pp2[1]*pp2[1] + pp2[2]*pp2[2];

   if( pm1 <= 0. || pm2 <= 0. ) return -999;

   const double e1 = sqrt(pow(m1, 2) + pm1);
   const double e2 = sqrt(pow(m2, 2) + pm2);
   const double es = e1 + e2;

   double ps[3];
   
   ps[0] = pp1[0] + pp2[0];
   ps[1] = pp1[1] + pp2[1];
   ps[2] = pp1[2] + pp2[2];
   
   const double mass = sqrt(es*es - ps[0]*ps[0] - ps[1]*ps[1] - ps[2]*ps[2]);
   
   return mass;
}

bool IsOneArmCut(const double phi1, const double phi2)
{
   if ( !((phi1 < 1.5 && phi2 < 1.5) || (phi1 > 1.5 && phi2 > 1.5)) ) return true;
   return false;
}

bool IsGhostCut(double dzed, double dphi, double dalpha)
{
   //pc1
   if (fabs(dzed) < 6.0 && fabs(dphi - (0.13*dalpha)) < 0.015) return true;
   //x1x2_1
   if (fabs(dphi - (0.04*dalpha)) < 0.015) return true;
   //x1x2_2
   if (fabs(dphi - (-0.065*dalpha)) < 0.015) return true;
   return false;
}

bool noCut()
{
   return true;
}

#endif /* P_TRACK_FUN_CPP */
