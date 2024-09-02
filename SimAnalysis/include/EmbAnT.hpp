// $HEADER$
//------------------------------------------------------------------------------------------------
//                            EmbAnT functinos declaration
//------------------------------------------------------------------------------------------------
// EmbAnT - embedding analysis tool
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic class functions for embedding efficiency evaluation
 **/
//------------------------------------------------------------------------------------------------

#ifndef EMB_AN_T_HPP
#define EMB_AN_T_HPP

#include "ErrorHandler.hpp"
#include "IOTools.hpp"
#include "Box.hpp"

#include "PBar.hpp"

#include "ThrObj.hpp"

#include "Particles.hpp"
#include "ParEmb.hpp"
#include "EmbTreeReader.hpp"
#include "STrackFun.hpp"

#include "ROOT/TTreeProcessorMT.hxx"

//container for storing ThrObj histograms
struct ThrHistStruct
{
	std::unique_ptr<ThrObj<TH1F>> reg_dc_pc1, reg_pc2, reg_pc3, reg_tofe, reg_tofw;
	std::array<std::unique_ptr<ThrObj<TH1F>>, 4> reg_emcale, reg_emcalw;
	ThrHistStruct(std::string run_type);
};

int GetEmcSector(const double phi, const double pemcy);
void Analyze(std::string part, const int queue_num);
void EmbAnalyze();
int main();

#endif
