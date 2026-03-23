#pragma once

#include "TCanvasTools.hpp"

#include "PainterHelper.hpp"

int main()
{
   gROOT->SetBatch(kTRUE);
   gStyle->SetOptStat(kFALSE);

   TCanvas canv("canv", "canv", 800, 800);

   canv.SetFillStyle(4000);
   canv.SetFrameFillColor(0);
   canv.SetFrameFillStyle(0);
   canv.SetFrameBorderMode(0);

   gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05); 
   gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.112);

   ROOTTools::DrawFrame(0.8, 0., 7., 2., "", "p_{T}", "R_{AB}");

   return 0;
}
