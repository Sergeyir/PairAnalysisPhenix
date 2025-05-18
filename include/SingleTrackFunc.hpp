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

/*! @brief Checks if the hit was detected in the detector
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
bool IsMatch(const double pT, const double sdphi, const double sdz, const double sDevSquareMax = 4.);
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
