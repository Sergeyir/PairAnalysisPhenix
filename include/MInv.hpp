/** 
 *  @file   MInv.hpp 
 *  @brief  Contains declarations of functions that are used for subtracting background and for merging invariant mass histograms in a given centrality, z, and r regions
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysis).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef M_INV_HPP
#define M_INV_HPP

#include <thread>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TROOT.h"
#include "TMath.h"

#include "StrTools.hpp"
#include "IOTools.hpp"
#include "MathTools.hpp"

/*! @namespace MInv
 * @brief Contains all functions and variables for MInv.cpp
 */
namespace MInv
{
   /*! Merges invariant mass distributions with subtracted background for all centrality (c in CabanaBoy), z_{vtx} (z in CabanaBoy) and r_{vtx} (r in CabanaBoy). Low resolution (LR) histograms are used for background scaling.
    *
    * @param[in] methodName name of the method that was used to extract pairs of charged tracks
    * @param[in] decayMode decay mode of which all histograms will be merged
    * @param[in] centralityBin centrality bin which will be used for invariant mass histogram merging
    * @param[in] pTBin pT bin which will be used for invariant mass histogram merging
    * @param[in] distrMInvMergedFG histogram to pass that will be filled with contents of all scaled foreground histograms
    * @param[in] distrMInvMergedBG histogram to pass that will be filled with contents of all scaled background histograms
    * @param[in] distrMInvMergedFGLR histogram to pass that will be filled with contents of all scaled foreground histograms with low resolution (for background scaling)
    * @param[in] distrMInvMergedBGLR histogram to pass that will be filled with contents of all scaled background histograms with low resolution (for background scaling)
    * @param[out] merged invariant mass distribution with background extracted
    */
   TH1D *Merge(TFile *inputFile, const std::string& methodName, const std::string& decayMode,
               const int cMin, const int cMax, const int zMin, const int zMax, 
               const int rMin, const int rMax, const double pTMin, const double pTMax,
               TH1D*& distrMInvMergedFG, TH1D*& distrMInvMergedBG,
               TH1D*& distrMInvMergedFGLR, TH1D*& distrMInvMergedBGLR, double& numberOfEvents);
   /*! Subtracts background for the specified histogram 
    *
    * @param[in] distrMInvFG foreground M_{inv} distribution from which background will be extracted
    * @param[in] distrMInvFG background M_{inv} distribution which will be extracted from foreground; in the process scaling will be applied
    * @param[out] invariant mass distribution with background subtracted
    */
   TH1D *SubtractBG(TH1D*& distrMInvFG, TH1D*& distrMInvBG, 
                    TH1D*& distrMInvFGLR, TH1D*& distrMInvBGLR);
};

#endif /* M_INV_HPP */
