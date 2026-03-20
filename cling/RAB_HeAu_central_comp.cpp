#include "TCanvasTools.hpp"

#include "PainterHelper.hpp"

void RAB_HeAu_central_comp()
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

   ROOTTools::DrawFrame(0.8, 0., 9., 2., "", "p_{T}", "R_{AB}");

   TLegend legend(0.5, 0.75, 0.95, 0.9);
   legend.SetLineColorAlpha(0, 0.);
   legend.SetFillColorAlpha(0, 0.);

   PainterHelper rab(&legend);
   rab.SetMarkerSize(1.6);
   rab.SetLineWidth(1);

   TLatex tlText;

   tlText.SetTextFont(52);
   tlText.SetTextSize(0.05);

   tlText.DrawLatexNDC(0.2, 0.15, "^{3}HeAu@200");
   tlText.DrawLatexNDC(0.75, 0.15, "#cbar#eta#cbar < 0.5");

   rab.DrawGraphFromTXTFile("data/RAB/HeAu200/ppbar0-20PHENIX.txt", kBlack, 0.4, 75,
                            "(p+#bar{p})/2 PRC109, 054910");
   rab.DrawLegend();
   ROOTTools::PrintCanvas(&canv, "tmp/RAB_HeAu_central1");

   rab.DrawGraphFromYAMLFile("data/RAB/HeAu200/pi0PHENIX.yaml", "0-20", kGreen - 3, 0.9, 77,
                             "#pi^{0} PRC105 064902");
   rab.DrawLegend();
   ROOTTools::PrintCanvas(&canv, "tmp/RAB_HeAu_central2");

   rab.DrawGraphFromYAMLFile("data/RAB/HeAu200/phi1020PHENIX.yaml", "0-20", kAzure - 3, 0.9, 74,
                             "#varphi(1020) PRC106, 014982");
   rab.DrawLegend();
   ROOTTools::PrintCanvas(&canv, "tmp/RAB_HeAu_central3");

   rab.DrawGraphFromTXTFile("data/RAB/HeAu200/KStar892_0-20.txt", kRed - 3, 0.9, 72,
                            "K^{*0}(892)");
   rab.DrawLegend();
   ROOTTools::PrintCanvas(&canv, "tmp/RAB_HeAu_central4");
}
