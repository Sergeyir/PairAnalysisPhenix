/** 
 *  @file   SingleTrackFunc.cpp
 *  @brief  Contains realisations of functions that are used for single tracks in PHENIX simulation analysis
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysisPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/

#ifndef SINGLE_TRACK_FUNC_CPP
#define SINGLE_TRACK_FUNC_CPP

#include "../include/SingleTrackFunc.hpp"

ChargedTrack::ChargedTrack(const double m, const SimTreeReader& simCNT, const int i)
{
   this->m = m;
   index = i;
   pX = simCNT.mom(i)*sin(simCNT.the0(i))*cos(simCNT.phi0(i));
   pY = simCNT.mom(i)*sin(simCNT.the0(i))*sin(simCNT.phi0(i));
   pZ = simCNT.mom(i)*cos(simCNT.the0(i));
   phi = simCNT.phi(i);
   alpha = simCNT.alpha(i);
   zed = simCNT.zed(i);
   p = sqrt(pX*pX + pY*pY + pZ*pZ);
   e = sqrt(p*p + m*m);
}

bool IsHit(const double dVal)
{
   if (dVal < -9998) return false;
   return true;
}

bool IsMatch(const double pT, const double sdphi, const double sdz, const double sDevSquareMax)
{
   if (sdphi*sdphi + sdz*sdz < sDevSquareMax) return true;
   return false;
}   

bool IsQualityCut(const int qual)
{
   if (qual != 63 && qual != 31) return true;
   return false;
}

double TransformProb(double prob)
{
   if (prob > 1.) prob = 1.;
   else if (prob < 0.) prob = 0.;
   return prob;
}

#endif /* SINGLE_TRACK_FUNC_CPP */
