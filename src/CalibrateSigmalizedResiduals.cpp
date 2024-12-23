// $SOURCE$
//------------------------------------------------------------------------------------------------
//                      CalibrateSigmailzedResiduals functions realisations
//------------------------------------------------------------------------------------------------
// CalibrateSigmalizedResiduals
//
// ** Free and open code for anyone to use **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic code for calibrating sigmalized residuals (sdphi, sdz)
 **/
//------------------------------------------------------------------------------------------------

#ifndef CALIBRATE_SIGMALIZED_RESIDUALS_CPP
#define CALIBRATE_SIGMALIZED_RESIDUALS_CPP

#include "../include/CalibrateSigmalizedResiduals.hpp"

int main(int argc, char **argv)
{
   if (argc != 2)
   {
      PrintError("Expected 1 parameters while " + std::to_string(argc - 1) + " were provided");
   }

   const std::string inputJSONName = argv[1];
   CheckInputFile(inputJSONName);
   std::ifstream inputJSON(inputJSONName.c_str(), std::ifstream::binary);
 
   Json::Value inputJSONContents;
   inputJSON >> inputJSONContents;

   ROOT::EnableImplicitMT();
   gErrorIgnoreLevel = kWarning;
   
   const std::string runName = inputJSONContents["RUN"]["name"].asString();
   const std::string inputFileName = "data/Real/" + runName + "/sum.root";

   TFile inputFile(inputFileName.c_str());

   TH3F *distrDPhi = static_cast<TH3F *>(inputFile.Get("dphi vs pT vs zed: EMCale0, charge < 0"));
   distrDPhi->RebinZ(3);
   distrDPhi->RebinY(2);

   TF1 fitFuncGaus("gaus", "gaus");
   fitFuncGaus.SetParameters(1., 0., 0.02);
   fitFuncGaus.SetParLimits(1, -0.02, 0.02);
   fitFuncGaus.SetParLimits(2, 0.001, 0.02);
   fitFuncGaus.SetParLimits(4, 0.02, 0.5);
   TF1 fitFuncDPhi("fitFunc", "gaus(0) + gaus(3)");
   fitFuncGaus.SetParameters(1., 0., 0.02);
   fitFuncGaus.SetParLimits(1, -0.015, 0.015);
   fitFuncDPhi.SetParLimits(2, 0.001, 0.02);
   TF1 fitFuncBG("bg", "gaus");
   fitFuncBG.SetParLimits(1, -0.15, 0.15);
   fitFuncBG.SetParLimits(2, 0.05, 0.5);

   fitFuncDPhi.SetLineColorAlpha(kRed+1, 0.7);
   fitFuncDPhi.SetLineWidth(5);
   fitFuncBG.SetLineColorAlpha(kGreen+1, 0.9);
   fitFuncBG.SetLineWidth(3);
   fitFuncBG.SetLineStyle(2);
   fitFuncGaus.SetLineColorAlpha(kAzure-3, 0.9);
   fitFuncGaus.SetLineWidth(3);
   fitFuncGaus.SetLineStyle(2);
   
   for (int i = 1; i <= distrDPhi->GetZaxis()->GetNbins(); i++)
   {
      TGraphErrors means, sigmas;

      means.SetMarkerStyle(20);
      sigmas.SetMarkerStyle(20);

      means.SetMarkerSize(0.5);
      sigmas.SetMarkerSize(0.5);

      Print(distrDPhi->GetZaxis()->GetBinLowEdge(i), distrDPhi->GetZaxis()->GetBinUpEdge(i));
      for (int j = 1; j <= distrDPhi->GetYaxis()->GetNbins(); j++)
      {
         if (distrDPhi->GetYaxis()->GetBinCenter(j) < 0.3) continue;
         
         TH1D *distrDPhiProj;
         
         if (distrDPhi->GetYaxis()->GetBinCenter(j) < 3.)
         {
            distrDPhiProj = distrDPhi->ProjectionX("", j, j, i, i);
         }
         else if (distrDPhi->GetYaxis()->GetBinCenter(j) < 6. && 
                  j < distrDPhi->GetYaxis()->GetNbins())
         {
            distrDPhiProj = distrDPhi->ProjectionX("", j, j + 1, i, i);
            j++;
         }
         else if (j < distrDPhi->GetYaxis()->GetNbins() - 4)
         {
            distrDPhiProj = distrDPhi->ProjectionX("", j, j + 4, i, i);
            j += 4;
         }
         else continue;

         double minX = 0., maxX = -1.;
         for (int k = 1; k <= distrDPhi->GetXaxis()->GetNbins(); k++)
         {
            if (distrDPhiProj->GetBinContent(k) > 1e-7)
            {
               minX = distrDPhiProj->GetXaxis()->GetBinLowEdge(k);
               break;
            }
         }

         for (int k = distrDPhi->GetXaxis()->GetNbins(); k > 1; k--)
         {
            if (distrDPhiProj->GetBinContent(k) > 1e-7)
            {
               maxX = distrDPhiProj->GetXaxis()->GetBinUpEdge(k);
               break;
            }
         }
         
         if (minX > maxX) 
         {
            PrintWarning("Hist is empty");
            continue;
         }

         distrDPhiProj->GetXaxis()->SetRange(distrDPhiProj->GetXaxis()->FindBin(minX+0.01),
                                             distrDPhiProj->GetXaxis()->FindBin(maxX-0.01));

         fitFuncDPhi.SetRange(minX, maxX);
         fitFuncBG.SetRange(minX, maxX);
         
         fitFuncGaus.SetParameter(0, distrDPhiProj->GetMaximum());
         fitFuncGaus.SetRange(-0.01, 0.01);
         distrDPhiProj->Fit(&fitFuncGaus, "RQMBN");
         distrDPhiProj->Fit(&fitFuncBG, "RQMBN");

         fitFuncGaus.SetRange(minX, maxX);
         
         for (int k = 0; k < 3; k++)
         {
            fitFuncDPhi.SetParameter(k, fitFuncGaus.GetParameter(k));
            fitFuncDPhi.SetParameter(k + 3, fitFuncBG.GetParameter(k));
         }
         
         TCanvas canv("", "", 400, 400);
         distrDPhiProj->Fit(&fitFuncDPhi, "RQMBN");

         for (int k = 0; k < 3; k++)
         {
            fitFuncGaus.SetParameter(k, fitFuncDPhi.GetParameter(k));
            fitFuncBG.SetParameter(k, fitFuncDPhi.GetParameter(k + 3));
         }
         
         distrDPhiProj->SetMarkerStyle(20);
         distrDPhiProj->SetMarkerSize(0.4);
         
         distrDPhiProj->Draw("P");
         fitFuncDPhi.Draw("SAME");
         fitFuncBG.Draw("SAME");
         fitFuncGaus.Draw("SAME");
         PrintCanvas(&canv, "tmp/" + std::to_string(i) + std::to_string(j));
         
         means.AddPoint(distrDPhi->GetYaxis()->GetBinCenter(j), fitFuncDPhi.GetParameter(1));
         means.SetPointError(means.GetN() - 1, fitFuncDPhi.GetParError(1));
         sigmas.AddPoint(distrDPhi->GetYaxis()->GetBinCenter(j), fitFuncDPhi.GetParameter(2));
         sigmas.SetPointError(sigmas.GetN() - 1, fitFuncDPhi.GetParError(2));
      }
      
      TCanvas meanCanv("means", "means", 600, 600);
      gPad->DrawFrame(0., -0.012, 10., 0.012);
      means.Draw("P");
      PrintCanvas(&meanCanv, "tmp/means" + std::to_string(i));
      
      TCanvas sigmaCanv("sigmas", "sigmas", 600, 600);
      gPad->DrawFrame(0., 0., 10., 0.015);
      sigmas.Draw("P");
      PrintCanvas(&sigmaCanv, "tmp/sigmas" + std::to_string(i));
   }
}

#endif /* CALIBRATE_SIGMALIZED_RESIDUALS_CPP */
