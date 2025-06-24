/** 
 *  @file   PairTrackFunc.hpp
 *  @brief  Contains declarations of functions and structs that are used for analysis of pair of tracks of PHENIX simulated data
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysisPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef PAIR_TRACK_FUNC_HPP
#define PAIR_TRACK_FUNC_HPP

#include <cmath>

#include "SingleTrackFunc.hpp"

/*! @brief Returns true if the pair of charged tracks passes TOF2PID check
 * @param[in] track1 1st charged track in a pair
 * @param[in] track1 2nd charged track in a pair
 * @param[in] id1 id that 1st track must have
 * @param[in] id2 id that 2nd track must have
 */
bool IsTOF2PID(const ChargedTrack& track1, const ChargedTrack& track2, 
               const int id1, const int id2);
/*! @brief Returns true if the pair of charged tracks passes EMC2PID check
 * @param[in] track1 1st charged track in a pair
 * @param[in] track1 2nd charged track in a pair
 * @param[in] id1 id that 1st track must have
 * @param[in] id2 id that 2nd track must have
 */
bool IsEMCal2PID(const ChargedTrack& track1, const ChargedTrack& track2, 
                 const int id1, const int id2);
/*! @brief Returns true if the pair of charged tracks passes TOF2PID check
 * @param[in] track1 1st charged track in a pair
 * @param[in] track1 2nd charged track in a pair
 * @param[in] id1 id that 1st track must have
 * @param[in] id2 id that 2nd track must have
 */
bool Is2PID(const ChargedTrack& track1, const ChargedTrack& track2, 
            const int id1, const int id2);
/*! @brief Returns true if the pair of charged tracks passes 2PID check
 * @param[in] track1 1st charged track in a p1ir
 * @param[in] track1 2nd charged track in a pair
 * @param[in] id1 id that 1st might have 
 * @param[in] id2 id that 2nd might have
 */
bool Is1PID(const ChargedTrack& track1, const ChargedTrack& track2, 
            const int id1, const int id2);
/*! @brief Returns true if the pair of charged tracks passes NoPID check
 * @param[in] track1 1st charged track in a pair
 * @param[in] track1 2nd charged track in a pair
 */
bool IsNoPID(const ChargedTrack& track1, const ChargedTrack& track2);
/*! @brief Returns pT of a pair of charged tracks
 * @param[in] track1 1st charged track in a pair
 * @param[in] track1 2nd charged track in a pair
 */
double GetPairPT(const ChargedTrack& track1, const ChargedTrack& track2);
/*! @brief Returns invariant mass of a pair of charged tracks
 * @param[in] track1 1st charged track in a pair
 * @param[in] track1 2nd charged track in a pair
 */
double GetPairMass(const ChargedTrack& track1, const ChargedTrack& track2);
/*! @brief Returns true if the pair of charged tracks does not pass 1 arm check
 * (i.e. both particles must be in one arm)
 * @param[in] track1 1st charged track in a pair
 * @param[in] track1 2nd charged track in a pair
 */
bool IsOneArmCut(const ChargedTrack& track1, const ChargedTrack& track2);
/*! @brief Returns true if the pair of charged tracks does not pass pair ghost cut
 * @param[in] track1 1st charged track in a pair
 * @param[in] track1 2nd charged track in a pair
 */
bool IsGhostCut(const ChargedTrack& track1, const ChargedTrack& track2);
/*! @brief Returns true if the pair of charged tracks does not pass sailor cut
 * @param[in] track1 1st charged track in a pair
 * @param[in] track1 2nd charged track in a pair
 */
bool IsSailorCut(const ChargedTrack& track1, const ChargedTrack& track2);
/*! @brief Returns true if the pair of charged tracks does not pass cowboy cut
 * @param[in] track1 1st charged track in a pair
 * @param[in] track1 2nd charged track in a pair
 */
bool IsCowboyCut(const ChargedTrack& track1, const ChargedTrack& track2);

#endif /* PAIR_TRACK_FUNC_HPP */
