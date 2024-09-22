// $HEADER$
//------------------------------------------------------------------------------------------------
//                              IdentFun functions declaration
//------------------------------------------------------------------------------------------------
// IdentFun
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic functions used for charged particle tracks identification
 **/
//------------------------------------------------------------------------------------------------

#ifndef IDENT_FUN_HPP
#define IDENT_FUN_HPP

#include <cmath>

#include "MathTools.hpp"

#include "Particles.hpp"
#include "STrackFun.hpp"

#include "M2IdentPar.hpp"

double GetEMCalId(const double pt, const int id, const int charge, 
                  const double phi, const int sector);
int GetTOFePID(const double pt, const int charge, const double m2);
int GetTOFwPID(const double pt, const int charge, const double m2);

#endif /* IDENT_FUN */
