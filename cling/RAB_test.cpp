#pragma once

#include "TCanvasTools.hpp"

#include "PainterHelper.hpp"

void RAB_test()
{
   gROOT->SetBatch(kTRUE);
   gStyle->SetOptStat(kFALSE);

   TCanvas centralCanv("canv", "", 800, 800);
   centralCanv.SetFillStyle(4000);
   centralCanv.SetFrameFillColor(0);
   centralCanv.SetFrameFillStyle(0);
   centralCanv.SetFrameBorderMode(0);

   TLegend legend(0.65, 0.7, 0.95, 0.95);
   legend.SetLineColorAlpha(0, 0.);
   legend.SetFillColorAlpha(0, 0.);
   legend.SetNColumns(3);

   ROOTTools::DrawFrame(0.8, 0., 9., 2., "", "p_{T}", "R_{AB}");

   PainterHelper rab(&legend);

   rab.DrawGraphFromTXTFile("data/RAB/AuAu200/KStar892_0-20.txt", kRed, 0.8, 20, "K^{*0}(892) 0-20\%");
   rab.DrawGraphFromTXTFile("data/RAB/HeAu200/KStar892_0-20.txt", kCyan, 0.8, 20, "K^{*0}(892) 0-20\%");

   rab.DrawLegend();

   centralCanv.SaveAs("tmp/RAB.png");
   //ROOTTools::PrintCanvas(&centralCanv, "tmp/RAB");
}
