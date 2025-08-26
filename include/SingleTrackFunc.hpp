/** 
 *  @file   SingleTrackFunc.hpp
 *  @brief  Contains declarations of functions that are used for single tracks in PHENIX simulation analysis
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysisPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef SINGLE_TRACK_FUNC_HPP
#define SINGLE_TRACK_FUNC_HPP

#include <cmath>

#include "SimTreeReader.hpp"

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
   /// Default constructor (deleted)
   ChargedTrack() = delete;
   /*! @brief Constructor
    * @param[in] m mass of a particle [GeV/c^2]
    * @param[in] pX X component of momentum [GeV/c]
    * @param[in] pY Y component of momentum [GeV/c]
    * @param[in] pz Z component of momentum [GeV/c]
    */
   ChargedTrack(const double m, const SimTreeReader& simCNT, const int i);
   /// mass of a particle [GeV/c^2]
   double m;
   /// index of SimTreeReader in the current event
   int index;
   /// reconstructed X compoment of momentum of a track [GeV/c]
   double pX;
   /// reconstructed Y compoment of momentum of a track [GeV/c]
   double pY;
   /// reconstructed Z compoment of momentum of a track [GeV/c]
   double pZ;
   /// reconstructed momentum of a track [GeV/c^2]
   double p;
   /// energy of a track [GeV]
   double e;
   /// reconstructed azimuthal angle of a track [rad]
   double phi;
   /// reconstructed polar angle of a track [rad]
   double alpha;
   /// reconstructed zDC coordinate [cm]
   double zed;
   /// reconstructed phi obtained from PC2
   double pc2phi;
   /// reconstructed z obtained from PC2
   double pc2z;
   /// reconstructed phi obtained from PC3
   double pc3phi;
   /// reconstructed z obtained from PC3
   double pc3z;
   /// EMCal sector
   int sector;
   /// y tower number in EMCal
   int yTower;
   /// z tower number in EMCal
   int zTower;
   /// slat number in TOFe
   int slat;
   /// strip number in TOFw
   int strip;
   /// id of a particle obtained in PC2
   int idPC2 = PART_ID::JUNK;
   /// id of a particle obtained in PC3
   int idPC3 = PART_ID::JUNK;
   /// id of a particle obtained in EMCal
   int idEMCal = PART_ID::JUNK;
   /// id of a particle obtained in TOFe
   int idTOFe = PART_ID::JUNK;
   /// id of a particle obtained in TOFw
   int idTOFw = PART_ID::JUNK;
};

/*! @brief Checks if the hit was detected in the given detector
 *
 * @param[in] dVal deviation of track (dphi or dz) of the detector the hit is needed to be checked in
 */
bool IsHit(const double dVal);
/*! @brief Checks if the charged track deviation is matched within required sigma of charged tracks deviation distribution
 *
 * @param[in] pT transverse momentum [GeV]
 * @param[in] sdphi sigmalized deviation of phi 
 * @param[in] sdz sigmalized deviation of z
 * @param[in] sDevSquareMax maximum allowed deviation of (sdphi*sdphi + sdz*sdz)
 */
bool IsMatch(const double pT, const double sdphi, const double sdz, 
             const double sDevSquareMax = 9.);
/*! @brief Checks if the quality cut for the charged track is applied
 *
 * @param[in] qual - quality of the charged track
 */
bool IsQualityCut(const int qual);
/*! @brief Transforms probability if it bigger than 1
 *
 * Used if the probability is calculated as a ratio of 2 variables in which case probability can be bigger than 1 if the ratio is taken with the wrong denominator.
 *
 * @param[in] prob probability
 */
double TransformProb(double prob);

#endif /* SINGLE_TRACK_FUNC_HPP */
