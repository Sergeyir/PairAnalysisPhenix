/** 
 *  @file   PairTrackFunc.cpp
 *  @brief  Contains realisations of functions and structs that are used for analysis of pair of tracks of PHENIX simulated data
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysisPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef PAIR_TRACK_FUNC_CPP
#define PAIR_TRACK_FUNC_CPP

#include "../include/PairTrackFunc.hpp"

ChargedTrack::ChargedTrack(const double m, const double pX, const double pY, const double pZ, 
                           const double phi, const double alpha, const double zed)
{
   this->m = m;
   this->pX = pX;
   this->pY = pY;
   this->pZ = pZ;
   this->phi = phi;
   this->alpha = alpha;
   this->zed = zed;
   e = sqrt(pX*pX + pY*pY + pZ*pZ + m*m);
}

bool IsTOF2PID(const ChargedTrack& track1, const ChargedTrack& track2, 
               const int id1, const int id2)
{
   return ((track1.idTOFe == id1 && track2.idTOFe == id2) ||
           (track1.idTOFw == id1 && track2.idTOFw == id2));
}
bool IsEMCal2PID(const ChargedTrack& track1, const ChargedTrack& track2, 
                 const int id1, const int id2)
{
   return (track1.idEMCal == id1 && track2.idEMCal == id2);
}
bool Is2PID(const ChargedTrack& track1, const ChargedTrack& track2, 
            const int id1, const int id2)
{
   return (IsTOF2PID(track1, track2, id1, id2) || 
           IsEMCal2PID(track1, track2, id1, id2) ||
           (track1.idEMCal == id1 && track2.idTOFe == id2) ||
           (track1.idEMCal == id1 && track2.idTOFw == id2) ||
           (track1.idTOFe == id1 && track2.idEMCal == id2) ||
           (track1.idEMCal == id1 && track2.idTOFw == id2));
           
}
bool Is1PID(const ChargedTrack& track1, const ChargedTrack& track2, const int id1, const int id2)
{
   return (IsNoPID(track1, track2) && 
           (track1.idTOFe == id1 || track1.idTOFw == id1 ||
            track2.idTOFe == id2 || track2.idTOFw == id2));
}
bool IsNoPID(const ChargedTrack& track1, const ChargedTrack& track2)
{
   // pairs are only selected in a way to satisfy IsNoPID check for now
   return true;
}
double GetPairPT(const ChargedTrack& track1, const ChargedTrack& track2)
{
   return sqrt((track1.pX + track2.pX)*(track1.pX + track2.pX) + 
               (track1.pY + track2.pY)*(track1.pY + track2.pY));
}
double GetPairMass(const ChargedTrack& track1, const ChargedTrack& track2)
{
   return sqrt((track1.e + track2.e)*(track1.e + track2.e) - 
               (track1.pX + track2.pX)*(track1.pX + track2.pX) - 
               (track1.pY + track2.pY)*(track1.pY + track2.pY) -
               (track1.pZ + track2.pZ)*(track1.pZ + track2.pZ));
}
bool IsOneArmCut(const ChargedTrack& track1, const ChargedTrack& track2)
{
   return ((track1.phi > M_PI/2. && track2.phi > M_PI/2.) ||
           (track1.phi < M_PI/2. && track2.phi < M_PI/2.));
}
bool IsGhostCut(const ChargedTrack& track1, const ChargedTrack& track2)
{
   const double dPhi = track1.phi - track2.phi;
   const double dAlpha = track1.alpha - track2.alpha;
   const double dZed = track1.zed - track2.zed;

   return (fabs(dZed) < 6.0 || 
           fabs(dPhi - (0.13*dAlpha)) < 0.015 ||
           fabs(dPhi - (0.04*dAlpha)) < 0.015 ||
           fabs(dPhi - (-0.065*dAlpha)) < 0.015);
}

#endif /* PAIR_TRACK_FUNC_CPP */
