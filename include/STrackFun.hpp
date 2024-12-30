// $HEADER$
//------------------------------------------------------------------------------------------------
//                                 STrackFun declaration
//------------------------------------------------------------------------------------------------
// STrackFun - single track functions
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic functions used in analysis for single particle tracks
 **/
//------------------------------------------------------------------------------------------------

#ifndef S_TRACK_FUN_HPP
#define S_TRACK_FUN_HPP

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

bool IsHit(const double dval);
bool IsMatch(const double sdphi, const double sdz, 
             const double sdphiMax = 2., const double sdzMax = 2.);
double TransformProb(double prob);
bool IsQualityCut(const int qual);
double GetM2Mean(const double pt, const double *par);
double GetM2Sigma(const double pt, const double m2Mean, const double *par);

#endif /*S_TRACK_FUN_HPP*/
