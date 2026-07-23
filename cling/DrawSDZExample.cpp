#include "TCanvasTools.hpp"

void DrawSDZExample()
{
   gROOT->SetBatch(kTRUE);
   gStyle->SetOptStat(0);

   std::filesystem::create_directories("output/other");

   TFile *inputFile = TFile::Open("data/Real/Run15pAl200/SingleTrack/sum.root");

   TH2F *distrSDZ = static_cast<TH2F *>(inputFile->Get("sdz vs pT, EMCale0, charge>0"));
   TH2F *distrSDPhi = static_cast<TH2F *>(inputFile->Get("sdphi vs pT, EMCale0, charge>0"));

   TCanvas canv("canv", "", 600, 600);

   gPad->SetRightMargin(0.11); gPad->SetTopMargin(0.07); 
   gPad->SetLeftMargin(0.105); gPad->SetBottomMargin(0.1);

   gPad->SetLogz();

   ROOTTools::DrawFrame(distrSDZ, "", "#it{sdz}", "#it{p}_{T}", 
                        0.9, 1., 0.05, 0.05, true, true, "COLZ");
   canv.SaveAs("output/other/sdz_example.pdf");

   canv.Clear();

   gPad->SetRightMargin(0.11); gPad->SetTopMargin(0.07); 
   gPad->SetLeftMargin(0.105); gPad->SetBottomMargin(0.1);

   gPad->SetLogz();

   ROOTTools::DrawFrame(distrSDPhi, "", "#it{sd#varphi}", "#it{p}_{T}", 
                        0.9, 1., 0.05, 0.05, true, true, "COLZ");
   canv.SaveAs("output/other/sdphi_example.pdf");
}
