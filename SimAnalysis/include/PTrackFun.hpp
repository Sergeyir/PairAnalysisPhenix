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

//container for storing id of particles
struct idContainer
{
   int tof;
   int emc;
   int pc2;
   int pc3;
   int origId;
};

bool IsTOF2PID(idContainer *Id1, idContainer *Id2);
bool IsEMC2PID(idContainer *Id1, idContainer *Id2);
bool IsEMCnoPID(idContainer *Id1, idContainer *Id2);
bool Is1TOFandEMC2PID(idContainer *Id1, idContainer *Id2);
bool Is2PID(idContainer *Id1, idContainer *Id2);
bool Is1PID(idContainer *Id1, idContainer *Id2);
bool IsnoPID(idContainer *Id1, idContainer *Id2);
double GetPairPt(const double *pp1, const double *pp2);
double GetMass(const double *pp1, const double *pp2, const double m1, const double m2);
bool IsOneArmCut(const double phi1, const double phi2);
bool IsGhostCut(double dzed, double dphi, double dalpha);
bool noCut();

#endif /* P_TRACK_FUN_HPP */
