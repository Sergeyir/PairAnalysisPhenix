// $HEADER$
//------------------------------------------------------------------------------------------------
//                            RunConfiguration flags
//------------------------------------------------------------------------------------------------
// RunConfiguration
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Flags and constant variables that determine run configuration
 **/
//------------------------------------------------------------------------------------------------

#ifndef RUN_CONFIGURATION_HPP
#define RUN_CONFIGURATION_HPP

#include "GlobalConfiguration.h"

#ifdef RUN7AUAU200

static const char *NUMBERED_RUN_NAME = "Run7";
static const char *COLLISION_SYSTEM_NAME = "AuAu200";
static const char *FULL_RUN_NAME = "Au+Au@#sqrt{s_{NN}} = 200 GeV";

#endif /* RUN7AUAU200 */

#endif /* RUN_CONFIGURATION_HPP */
