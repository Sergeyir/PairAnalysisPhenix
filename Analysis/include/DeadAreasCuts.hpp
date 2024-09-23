// $HEADER$
//------------------------------------------------------------------------------------------------
//                       DeadAreasCuts functions declaration
//------------------------------------------------------------------------------------------------
// DeadAreasCuts
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic set of functions for rejecting data in bad/dead areas of different detectors
 **/
//------------------------------------------------------------------------------------------------

#ifndef DEAD_AREAS_CUTS_HPP
#define DEAD_AREAS_CUTS_HPP

#include "GlobalConfiguration.h"

bool IsDeadDC(const double phi, const double zed, const double board, const double alpha);
bool IsDeadPC1(const double phi, const double pc1z, const double pc1phi);
bool IsDeadPC2(const double pc2z, const double pc2phi);
bool IsDeadPC3(const double phi, const double pc3z, const double pc3phi);
bool IsDeadEMCal(const double phi, const double zed, const int sect, 
                 const double pemcy, const double pemcz);
bool IsBadSlat(const int slat);
bool IsBadStripTOFw(const int strip);
bool IsDeadTOFe(const double zed, const double tofy, const double tofz);
bool IsDeadTOFw(const double zed, const double board, const double alpha);

#endif /* DEAD_AREAS_CUTS_HPP */
