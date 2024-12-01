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

#ifndef P_TRACK_FUN_HPP
#define P_TRACK_FUN_HPP

#include <cmath>

#include "Particles.hpp"

//container for storing id of particles
struct idContainer
{
   int tof;
   int emc;
   int pc2;
   int pc3;
   int origId;
};

bool IsTOF2PID(idContainer *id1, idContainer *id2);
bool IsEMC2PID(idContainer *id1, idContainer *id2);
bool IsEMCnoPID(idContainer *id1, idContainer *id2);
bool Is1TOFandEMC2PID(idContainer *id1, idContainer *id2);
bool Is2PID(idContainer *id1, idContainer *id2);
bool Is1PID(idContainer *id1, idContainer *id2);
bool IsnoPID(idContainer *id1, idContainer *id2);
double GetPairPT(const double *mom1, const double *mom2);
double GetMass(const double *mom1, const double *mom2, const double m1, const double m2);
bool IsOneArmCut(const double phi1, const double phi2);
bool IsGhostCut(double dZed, double dPhi, double dAlpha);
bool noCut();

#endif /* P_TRACK_FUN_HPP */
