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

bool IsPC2NoPID(const ChargedTrack& track1, const ChargedTrack& track2)
{
   return (track1.idPC2 != PART_ID::JUNK && track2.idPC2 != PART_ID::JUNK && 
           !IsPC2Ghost(track1, track2));
}

bool IsPC3NoPID(const ChargedTrack& track1, const ChargedTrack& track2)
{
   return (track1.idPC3 != PART_ID::JUNK && track2.idPC3 != PART_ID::JUNK && 
           !IsPC3Ghost(track1, track2));
}

bool IsEMCalNoPID(const ChargedTrack& track1, const ChargedTrack& track2)
{
   return (track1.idEMCal != PART_ID::JUNK && track2.idEMCal != PART_ID::JUNK && 
           !IsEMCalGhost(track1, track2));
}

bool IsTOFeNoPID(const ChargedTrack& track1, const ChargedTrack& track2)
{
   return (track1.idTOFe != PART_ID::JUNK && track2.idTOFe != PART_ID::JUNK && 
           !IsTOFeGhost(track1, track2));
}

bool IsTOFwNoPID(const ChargedTrack& track1, const ChargedTrack& track2)
{
   return (track1.idTOFw != PART_ID::JUNK && track2.idTOFw != PART_ID::JUNK && 
           !IsTOFwGhost(track1, track2));
}

bool IsNoPID(const ChargedTrack& track1, const ChargedTrack& track2)
{
   return (IsPC2NoPID(track1, track2) ||
           IsPC3NoPID(track1, track2) ||
           IsEMCalNoPID(track1, track2),
           IsTOFeNoPID(track1, track2) ||
           IsTOFwNoPID(track1, track2));
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
   return ((track1.phi > M_PI/2. && track2.phi < M_PI/2.) ||
           (track1.phi < M_PI/2. && track2.phi > M_PI/2.));
}

bool IsPC2Ghost(const ChargedTrack& track1, const ChargedTrack& track2)
{
   return (fabs(track1.pc2phi - track2.pc2phi) < 0.02 && fabs(track1.pc2z - track2.pc2z) < 4.);
}

bool IsPC3Ghost(const ChargedTrack& track1, const ChargedTrack& track2)
{
   return (fabs(track1.pc3phi - track2.pc3phi) < 0.02 && fabs(track1.pc3z - track2.pc3z) < 4.);
}

bool IsEMCalGhost(const ChargedTrack& track1, const ChargedTrack& track2)
{
   return (track1.sector == track2.sector && 
           abs(track1.yTower - track1.yTower) < 2 && 
           abs(track1.zTower - track2.zTower) < 2);
}

bool IsTOFeGhost(const ChargedTrack& track1, const ChargedTrack& track2)
{
   return (track1.slat == track2.slat);
}

bool IsTOFwGhost(const ChargedTrack& track1, const ChargedTrack& track2)
{
   return (track1.strip == track2.strip);
}

bool IsGhostCut(const ChargedTrack& track1, const ChargedTrack& track2)
{
   const double dPhi = track1.phi - track2.phi;
   const double dAlpha = track1.alpha - track2.alpha;
   const double dZed = track1.zed - track2.zed;

   return ((fabs(dZed) < 6.0 && fabs(dPhi - (0.13*dAlpha)) < 0.015) ||
           fabs(dPhi - (0.04*dAlpha)) < 0.015 ||
           fabs(dPhi - (-0.065*dAlpha)) < 0.015);
}

bool IsSailorCut(const ChargedTrack& track1, const ChargedTrack& track2)
{
   return (track1.phi < track2.phi);
}

bool IsCowboyCut(const ChargedTrack& track1, const ChargedTrack& track2)
{
   return !IsSailorCut(track1, track2);
}

#endif /* PAIR_TRACK_FUNC_CPP */
