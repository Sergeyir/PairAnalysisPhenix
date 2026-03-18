/** 
 *  @file   PainterHelper.cpp
 *  @brief  Contains realisation of class PainterHelper
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/CalPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef PAINTER_HELPER_CPP
#define PAINTER_HELPER_CPP

#include "PainterHelper.hpp"


PainterHelper::PainterHelper() {}

PainterHelper::PainterHelper(const double markerSize, const int lineWidth, const double sysWidth)
{
   this->markerSize = markerSize;
   this->lineWidht = lineWidth;
   this->sysWiddht = sysWidth;
}

void PainterHelper::DrawHist(TH1D *histogramWithStatErr, TH1D *histogramWithSysErr, 
                             const Color_t color, const double alpha, 
                             const Style_t markerStyle, const std::string& legendEntry)
{
   bool disableSysErr = false;

   if (!histogramWithStatErr) CppTools::PrintError("Histogram containing values and "\
                                                   "statistical uncertanities is nullptr");
   if (!histogramWithSysErr) disableSysErr = true;

   if (!TH1::CheckEqualAxes(histogramWithStatErr->GetXaxis(), 
                            histogramWithStatErr->GetXaxis())) disableSysErr = true;

   histogramWithStatErr->SetMarkerStyle(markerStyle);
   histogramWithStatErr->SetMarkerSize(markerSize);

   histogramWithStatErr->SetMarkerColorAlpha(color, alpha);
   histogramWithStatErr->SetLineColorAlpha(color, alpha);
   histogramWithStatErr->LineWidth(lineWidth);

   histogramWithStatErr->Draw("SAME P E1 X0");

   if (!disableSysErr)
   {
      TGraphErrors
   }

   histogramWithSysErr->SetMarkerStyle(9);
   histogramWithSysErr->SetMarkerSize(0);
   histogramWithSysErr->SetMarkerColorAlpha(0, 0.);
   histogramWithSysErr->SetLineColorAlpha(color, alpha);
   histogramWithSysErr->LineWidth(lineWidth);
}

void PainterHelper::DrawGraph(TGraphErrors *graphWithStatErr, TGraphErrors *graphWithSysErr, 
                              const Color_t color, const double alpha, 
                              const Style_t markerStyle, const std::string& legendEntry)
{
}

void PainterHelper::DrawHistFromROOTFile(const std::string& fileName, const std::string histogramName, 
                                         const Color_t color, const double alpha, 
                                         const Style_t markerStyle, const std::string& legendEntry)
{
}

void PainterHelper::DrawGraphFromROOTFile(const std::string& fileName, const std::string graphName, 
                                          const Color_t color, const double alpha, 
                                          const Style_t markerStyle, const std::string& legendEntry)
{
}

void PainterHelper::DrawFromYAMLFile(const std::string& fileName, const std::string& qualifier, 
                                     const Color_t color, const double alpha, 
                                     const Style_t markerStyle, const std::string& legendEntry = "")
{
}

void PainterHelper::DrawGraphFromTXTFile(const std::string& fileName, 
                                         const Color_t color, const double alpha,
                                         const Style_t markerStyle, const std::string& legendEntry)
{
}

TGraphErrors *PainterHelper::GetGraphFromYAMLFile(const std::string& fileName, const std::string& qualifier, 
                                                  TGraphErrors*& graphWithSysErrors, bool readSysErrors = true)
{
}

TGraphErrors *PainterHelper::GetGraphFromTXTFile(const std::string& fileName, TGraphErrors*& graphWithSysErrors, 
                                                 bool readSysErrors = true)
{
}

void PainterHelper::DrawTypeCUncertainty(const double value, const double xPos, const double yPos)
{
}


PainterHelper::~PainterHelper() {}


#endif /* PAINTER_HELPER_CPP */
