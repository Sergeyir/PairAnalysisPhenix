/** 
 *  @file   PairTrackFunc.hpp
 *  @brief  Contains declarations of functions that are used for analysis of pair of tracks of PHENIX simulated data
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysisPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef PAIR_TRACK_FUNC_HPP
#define PAIR_TRACK_FUNC_HPP

#include <cmath>

#include "Particles.hpp"

/*! @namespace PART_ID
 * @brief Contains constants that store absolute values of particle ids
 */
namespace PART_ID
{
   const int PION = 211;
   const int KAON = 321;
   const int PROTON = 2212;
   const int ELECTRON = 11;
   const int NONE = 999;
   const int JUNK = 1999;
}

/*! @struct ChargedTrack
 * @brief Convenient data container for analysing simulated charged tracks
 */
struct ChargedTrack
{
   /*! @brief Constructor
    * @param[in] m mass of a particle [GeV/c^2]
    * @param[in] pX X component of momentum [GeV/c]
    * @param[in] pY Y component of momentum [GeV/c]
    * @param[in] pz Z component of momentum [GeV/c]
    */
   ChargedTrack(const double m, const double pX, const double pY, const double pZ);
   /// mass of a particle [GeV/c^2]
   double m;
   /// reconstructed X compoment of momentum of a track [GeV/c]
   double pX;
   /// reconstructed Y compoment of momentum of a track [GeV/c]
   double pY;
   /// reconstructed Z compoment of momentum of a track [GeV/c]
   double pZ;
   /// energy of a track [GeV]
   double e;
   /// reconstructed azimuthal angle of a track
   double phi;
   /// reconstructed polar angle of a track
   double theta;
   /// id of a particle obtained in PC2
   int idPC2;
   /// id of a particle obtained in PC3
   int idPC3;
   /// id of a particle obtained in EMCal
   int idEMCal;
   /// id of a particle obtained in TOFe
   int idTOFe;
   /// id of a particle obtained in TOFw
   int idTOFw
}

bool IsTOF2PID(const ChargedTrack& track1, const ChargedTrack& track2);
bool IsEMC2PID(const ChargedTrack& track1, const ChargedTrack& track2);
bool Is2PID(const ChargedTrack& track1, const ChargedTrack& track2);
bool Is1PID(const ChargedTrack& track1, const ChargedTrack& track2);
bool IsNoPID(const ChargedTrack& track1, const ChargedTrack& track2);
double GetPairPT(const ChargedTrack& track1, const ChargedTrack& track2);
double GetMass(const ChargedTrack& track1, const ChargedTrack& track2);
bool IsOneArmCut(const ChargedTrack& track1, const ChargedTrack& track2);
bool IsGhostCut(const ChargedTrack& track1, const ChargedTrack& track2);

#endif /* PAIR_TRACK_FUNC_HPP */
