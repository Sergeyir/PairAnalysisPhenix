#include "TCanvasTools.hpp"

#include "PainterHelper.hpp"

void RAB_central()
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

   TLegend legend(0.4, 0.75, 0.95, 0.9);
   legend.SetLineColorAlpha(0, 0.);
   legend.SetFillColorAlpha(0, 0.);

   PainterHelper rab(&legend);
   rab.SetMarkerSize(1.6);
   rab.SetLineWidth(1);

   rab.DrawGraphFromTXTFile("data/RAB/HeAu200/KStar892_0-20.txt", kAzure-3, 0.9, 72,
                            "K^{*0}(892) ^{3}HeAu@200 0-20\%");
   rab.DrawGraphFromTXTFile("data/RAB/AuAu200/KStar892_0-20.txt", kRed-3, 0.9, 71,
                            "K^{*0}(892) AuAu@200 0-20\%");

   rab.DrawLegend();

   TLatex tlText;

   tlText.SetTextFont(52);
   tlText.SetTextSize(0.05);

   tlText.DrawLatexNDC(0.75, 0.15, "#cbar#eta#cbar < 0.5");

   ROOTTools::PrintCanvas(&canv, "tmp/RAB_central");
}
