#pragma once

#include "TCanvasTools.hpp"

#include "PainterHelper.hpp"

void RAB_central()
{
   gROOT->SetBatch(kTRUE);
   gStyle->SetOptStat(kFALSE);

   TCanvas canv("canv", "canv", 1200, 800);

   canv.SetFillStyle(4000);
   canv.SetFrameFillColor(0);
   canv.SetFrameFillStyle(0);
   canv.SetFrameBorderMode(0);

   gPad->SetRightMargin(0.01); gPad->SetTopMargin(0.016); 
   gPad->SetLeftMargin(0.08); gPad->SetBottomMargin(0.112);

   const double xMin = 0.8;
   const double xMax = 7.;
   const double yMin = 0.;
   const double yMax = 2.;

   ROOTTools::DrawFrame(xMin, yMin, xMax, yMax, "", "p_{T} [GeV/c]", "R_{AB}", 1., 0.8);

   TLegend legend(0.4, 0.8, 0.95, 0.95);
   legend.SetLineColorAlpha(0, 0.);
   legend.SetFillColorAlpha(0, 0.);

   TLine line(xMin, 1., xMax, 1.);
   line.SetLineColor(kGray + 1);
   line.SetLineWidth(3);
   line.SetLineStyle(2);

   line.Draw();

   PainterHelper rab(&legend);
   rab.SetMarkerSize(1.6);
   rab.SetLineWidth(2);
   rab.SetSysWidth(0.08);

   rab.DrawGraphFromTXTFile("data/RAB/HeAu200/KStar892_0-20.txt", kAzure - 3, 0.9, 72,
                            "K^{*0}(892) ^{3}HeAu@200 0-20\%");
   rab.DrawGraphFromTXTFile("data/RAB/AuAu200/KStar892_0-20.txt", kRed - 3, 0.9, 71,
                            "K^{*0}(892) AuAu@200 0-20\%");

   rab.DrawLegend();

   rab.DrawTypeCUncertainty(0.1, 6.8, 1., kAzure - 3, 0.5);
   rab.DrawTypeCUncertainty(0.1, 6.6, 1., kRed - 3, 0.5);

   TLatex tlText;

   tlText.SetTextFont(52);
   tlText.SetTextSize(0.05);

   tlText.DrawLatexNDC(0.75, 0.15, "#cbar#eta#cbar < 0.5");

   ROOTTools::PrintCanvas(&canv, "tmp/RAB_central");
}
