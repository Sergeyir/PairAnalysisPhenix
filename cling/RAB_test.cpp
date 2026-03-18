#pragma once

#include "PainterHelper.hpp"
#include "TCanvasTools.hpp"

void RAB_test()
{
   gROOT->SetBatch(kTRUE);
   gStyle->SetOptStat(kFALSE);

   TCanvas centralCanv("canv", "", 800, 800);
   centralCanv.SetFillStyle(4000);
   centralCanv.SetFrameFillColor(0);
   centralCanv.SetFrameFillStyle(0);
   centralCanv.SetFrameBorderMode(0);

   TLegend legend(0.1, 0.3, 0.5, 0.8);

   TCanvasTools::DrawFrame(0.8, 0., 9., 2., "p_{T}", "R_{AB}");

   PainterHelper rab(&legend);
}
