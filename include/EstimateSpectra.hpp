/** 
 *  @file   EstimateSpectra.hpp 
 *  @brief  Contains declarations of functions and variables that are used for estimation of invariant pT spectra
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysis).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef ESTIMATE_SPECTRA_HPP
#define ESTIMATE_SPECTRA_HPP

#include <thread>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"

#include "StrTools.hpp"
#include "IOTools.hpp"
#include "MathTools.hpp"

#include "TCanvasTools.hpp"

#include "PBar.hpp"

#include "InputYAMLReader.hpp"
#include "FitFunc.hpp"

/*! @namespace EstimateResonanceEff
 * @brief Contains all functions and variables for EstimateResonanceEff.cpp
 */
namespace EstimateResonanceEff
{
   /*! Returns estimated invariant pT spectra for the given method
    *
    * @param[in] methodName name of the method that was used to extract pairs of charged tracsk
    * @param[in] spectraWithSysErrors histogram that will be filled with spectra values with systematic uncertainties
    *
    * @param[out] spectraWithStatErrors histogram that will be filled with spectra values with statistical uncertainties
    */
   TH1D *GetSpectraForMethod(const std::string& methodName, TH1D&* spectraWithSysErrors);
   /*! Applies bin shift correction along y axis to the passed histograms
    *
    * @param[in] spectraWithSysErrors histogram containing spectra values with systematic uncertainties
    * @param[in] spectraWithStatErrors histogram containing spectra values with statistical uncertainties
    */
   void ApplyBinShiftCorrection(TH1D&* spectraWithStatErrors, TH1D&* spectraWithSysErrors);
   /// Contents of input .yaml file for run configuration
   InputYAMLReader inputYAMLMain;
   /// Contents of input .yaml file for the information about resonance
   InputYAMLReader inputYAMLResonance;
   /// Name of run (e.g. Run14HeAu200 or Run7AuAu200)
   std::string runName;
   /// Name of the input file
   std::string inputFileName;
   /// Input file
   TFile *inputFile;
   /// Output file for writing invariant pT spectra
   TFile *outputFile;
   /// name of the resonance
   std::string resonanceName;
   /// number of pT bins
   unsigned int pTNBins;
   /// pT bins ranges [GeV/c]
   std::vector<double> pTBinRanges;
   /// TText object template for quick text insertions
   TText text;
   /// TLatex object template for quick text insertions
   TLatex texText;
};

#endif /* ESTIMATE_SPECTRA_HPP */
