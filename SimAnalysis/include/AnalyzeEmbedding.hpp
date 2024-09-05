// $HEADER$
//------------------------------------------------------------------------------------------------
//                          AnalyzeEmbedding function declaration
//------------------------------------------------------------------------------------------------
// AnalyzeEmbedding
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic macro for embedding study 
 * from simulation output of event-like TTrees to processed histograms
 **/
//------------------------------------------------------------------------------------------------

#ifndef ANALYZE_EMBEDDING_HPP
#define ANALYZE_EMBEDDING_HPP

#include "ErrorHandler.hpp"
#include "Box.hpp"

#include "PBar.hpp"

#include "ThrObj.hpp"

#include "Particles.hpp"
#include "ParAnalyzeEmbedding.hpp"
#include "EmbTreeReader.hpp"
#include "STrackFun.hpp"

#include "ROOT/TTreeProcessorMT.hxx"

//container for storing ThrObj histograms
struct ThrContainer
{
	std::unique_ptr<ThrObj<TH1F>> regDCPC1, regPC2, regPC3, regTOFe, regTOFw;
	std::array<std::unique_ptr<ThrObj<TH1F>>, 4> regEMCale, regEMCalw;
	ThrContainer(std::string runType);
};

int GetEmcSector(const double phi, const double pemcy);
void AnalyzeParticleEmbedding(std::string part, const int queueNum);
void AnalyzeEmbedding();
int main();

#endif /*ANALYZE_EMBEDDING_HPP*/
