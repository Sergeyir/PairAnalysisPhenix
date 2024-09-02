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
	std::unique_ptr<ThrObj<TH1F>> regDCPC1, regPC2, regPC3, regTOFe, regTOFw;
	std::array<std::unique_ptr<ThrObj<TH1F>>, 4> regEMCale, regEMCalw;
	ThrHistStruct(std::string runType);
};

int GetEmcSector(const double phi, const double pemcy);
void Analyze(std::string part, const int queueNum);
void EmbAnalyze();
int main();

#endif
