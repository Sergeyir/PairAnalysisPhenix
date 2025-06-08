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

bool IsHit(const double dVal)
{
   if (dVal < -9998) return false;
   return true;
}

bool IsMatch(const double pT, const double sdphi, const double sdz, const double sDevSquareMax)
{
   if (pT > 3. || sdphi*sdphi + sdz*sdz < sDevSquareMax) return true;
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
