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

   PrintInfo("Clearing output directory");
   system(("rm -r output/TrackResiduals/" + runName + "/*").c_str());

   for (auto detector : inputJSONContents["RUN"]["uncalibrated_sigmalized_residuals_detectors"])
   {
      for (int zDC = -60; zDC < 45; zDC += 15)
      {
         PerformFits(runName, detector["name"].asString(), "dphi", zDC, zDC + 15, true);
         PerformFits(runName, detector["name"].asString(), "dphi", zDC, zDC + 15, false);
         PerformFits(runName, detector["name"].asString(), "dz", zDC, zDC + 15, true);
         PerformFits(runName, detector["name"].asString(), "dz", zDC, zDC + 15, false);
      }
   }
}
   

void PerformFits(const std::string& runName, const std::string& detectorName, 
                 const std::string& dValName, const int zDCMin, const int zDCMax,
                 const bool isPositive)
{ 
   Print(runName, detectorName, dValName, zDCMin, zDCMax, isPositive);
   const std::string inputFileName = "data/Real/" + runName + "/SingleTrack/sum.root";
   CheckInputFile(inputFileName);
   TFile inputFile(inputFileName.c_str());

   const std::string zDCRangeName = std::to_string(zDCMin) + "<zDC<" + std::to_string(zDCMax);
   const std::string chargeName = (isPositive) ? "charge>0" : "charge<0";

   const std::string distrDValName = dValName + " vs pT vs centrality: " + 
                     detectorName + ", " + chargeName + ", " + zDCRangeName;
   TH3F *distrDVal = static_cast<TH3F *>(inputFile.Get(distrDValName.c_str()));
   if (!distrDVal) PrintError("Histogram named \"" + distrDValName + 
                              "\" does not exist in file " + inputFileName);
 
   const double nBins = distrDVal->GetXaxis()->GetNbins();
   const double minBinX = distrDVal->GetXaxis()->GetBinLowEdge(1);
   const double maxBinX = distrDVal->GetXaxis()->GetBinUpEdge(nBins);
   const double binWidth = distrDVal->GetXaxis()->GetBinWidth(1);
   
   TF1 fitFuncGaus("gaus", "gaus");
   fitFuncGaus.SetParameters(1., 0., binWidth*2.);
   fitFuncGaus.SetParLimits(1, Average(minBinX, minBinX, maxBinX), 
                               Average(maxBinX, maxBinX, minBinX));
   fitFuncGaus.SetParLimits(2, binWidth, Average(maxBinX, maxBinX, minBinX));
   
   TF1 fitFuncDVal("fitFunc", "gaus(0) + gaus(3)");
   fitFuncDVal.SetParameters(1., 0., binWidth*2.);
   fitFuncDVal.SetParLimits(1, Average(minBinX, minBinX, maxBinX), 
                               Average(maxBinX, maxBinX, minBinX));
   fitFuncDVal.SetParLimits(2, binWidth, Average(maxBinX, maxBinX, minBinX));
   fitFuncDVal.SetParLimits(4, minBinX, maxBinX);
   fitFuncDVal.SetParLimits(5, Average(maxBinX, maxBinX, minBinX), maxBinX);
   
   TF1 fitFuncBG("bg", "gaus");
   fitFuncBG.SetParameters(1., 0., Average(maxBinX, maxBinX, maxBinX, minBinX));
   fitFuncBG.SetParLimits(1, minBinX, maxBinX);
   fitFuncBG.SetParLimits(2, Average(maxBinX, maxBinX, minBinX), maxBinX);

   fitFuncDVal.SetLineColorAlpha(kRed+1, 0.7);
   fitFuncDVal.SetLineWidth(5);
   fitFuncBG.SetLineColorAlpha(kGreen+1, 0.9);
   fitFuncBG.SetLineWidth(3);
   fitFuncBG.SetLineStyle(2);
   fitFuncGaus.SetLineColorAlpha(kAzure-3, 0.9);
   fitFuncGaus.SetLineWidth(3);
   fitFuncGaus.SetLineStyle(2);
 
   for (int i = 1; i <= distrDVal->GetZaxis()->GetNbins(); i++)
   {
      std::string outputDir = "output/TrackResiduals/" + runName + "/" + detectorName;
      if (isPositive) outputDir += "/Pos/";
      else outputDir += "/Neg/";
      outputDir += "/zDC" + std::to_string(zDCMin) + "_" + std::to_string(zDCMax) + 
                   "/centr" + std::to_string((i - 1)*5) + std::to_string(i*5) + "/";
      
      system(("mkdir -p " + outputDir).c_str());
      
      TGraphErrors means, sigmas;

      means.SetMarkerStyle(20);
      sigmas.SetMarkerStyle(20);

      means.SetMarkerSize(0.5);
      sigmas.SetMarkerSize(0.5);

      double meansYMin = 1e31;
      double meansYMax = -1e31;
      double sigmasYMin = 1e31;
      double sigmasYMax = -1e31;

      for (int j = 1; j <= distrDVal->GetYaxis()->GetNbins(); j++)
      {
         double pTMin, pTMax;
         if (distrDVal->GetYaxis()->GetBinCenter(j) < 0.3) continue;
         
         TH1D *distrDValProj;
         
         if (distrDVal->GetYaxis()->GetBinCenter(j) < 2.)
         {
            distrDValProj = distrDVal->ProjectionX("", j, j, i, i);
            pTMin = distrDVal->GetYaxis()->GetBinLowEdge(j);
            pTMax = distrDVal->GetYaxis()->GetBinUpEdge(j);
         }
         else if (distrDVal->GetYaxis()->GetBinCenter(j) < 3. && 
                  j < distrDVal->GetYaxis()->GetNbins())
         {
            distrDValProj = distrDVal->ProjectionX("", j, j + 1, i, i);
            distrDValProj->Rebin(2);
            pTMin = distrDVal->GetYaxis()->GetBinLowEdge(j);
            pTMax = distrDVal->GetYaxis()->GetBinUpEdge(j + 1);
            j++;
         }
         else if (distrDVal->GetYaxis()->GetBinCenter(j) < 5. && 
                  j < distrDVal->GetYaxis()->GetNbins())
         {
            distrDValProj = distrDVal->ProjectionX("", j, j + 2, i, i);
            distrDValProj->Rebin(2);
            pTMin = distrDVal->GetYaxis()->GetBinLowEdge(j);
            pTMax = distrDVal->GetYaxis()->GetBinUpEdge(j + 2);
            j += 2;
         }
         else if (distrDVal->GetYaxis()->GetBinCenter(j) < 6. && 
                  j < distrDVal->GetYaxis()->GetNbins())
         {
            distrDValProj = distrDVal->ProjectionX("", j, j + 4, i, i);
            distrDValProj->Rebin(2);
            pTMin = distrDVal->GetYaxis()->GetBinLowEdge(j);
            pTMax = distrDVal->GetYaxis()->GetBinUpEdge(j + 4);
            j += 4;
         }
         else if (j < distrDVal->GetYaxis()->GetNbins() - 9)
         {
            distrDValProj = distrDVal->ProjectionX("", j, j + 9, i, i);
            distrDValProj->Rebin(4);
            pTMin = distrDVal->GetYaxis()->GetBinLowEdge(j);
            pTMax = distrDVal->GetYaxis()->GetBinUpEdge(j + 9);
            j += 9;
         }
         else continue;

         if (distrDValProj->Integral(1, distrDValProj->GetXaxis()->GetNbins()) < 1) continue;

         double minX = 0., maxX = -1.;
         for (int k = 1; k <= distrDVal->GetXaxis()->GetNbins(); k++)
         {
            if (distrDValProj->GetBinContent(k) > 1e-7)
            {
               minX = distrDValProj->GetXaxis()->GetBinLowEdge(k);
               break;
            }
         }

         for (int k = distrDVal->GetXaxis()->GetNbins(); k > 1; k--)
         {
            if (distrDValProj->GetBinContent(k) > 1e-7)
            {
               maxX = distrDValProj->GetXaxis()->GetBinUpEdge(k);
               break;
            }
         }
         
         if (minX > maxX) 
         {
            continue;
         }

         distrDValProj->GetXaxis()->SetRange(distrDValProj->GetXaxis()->FindBin(minX+0.01),
                                             distrDValProj->GetXaxis()->FindBin(maxX-0.01));

         fitFuncDVal.SetRange(minX, maxX);
         fitFuncBG.SetRange(minX, maxX);

         fitFuncGaus.SetParameter(0, distrDValProj->GetMaximum());
         fitFuncGaus.SetRange(Average(minBinX, minBinX, minBinX, maxBinX), 
                              Average(maxBinX, maxBinX, maxBinX, minBinX));
         distrDValProj->Fit(&fitFuncGaus, "RQMBN");
         distrDValProj->Fit(&fitFuncBG, "RQMBN");

         fitFuncGaus.SetRange(minX, maxX);
         
         for (int k = 0; k < 3; k++)
         {
            fitFuncDVal.SetParameter(k, fitFuncGaus.GetParameter(k));
            fitFuncDVal.SetParameter(k + 3, fitFuncBG.GetParameter(k));
         }
         
         TCanvas canv("", "", 400, 400);
         distrDValProj->Fit(&fitFuncDVal, "RQMBN");

         for (int k = 0; k < 3; k++)
         {
            fitFuncGaus.SetParameter(k, fitFuncDVal.GetParameter(k));
            fitFuncBG.SetParameter(k, fitFuncDVal.GetParameter(k + 3));
         }
         
         distrDValProj->SetMarkerStyle(20);
         distrDValProj->SetMarkerSize(0.4);
         
         distrDValProj->Draw("P");
         fitFuncDVal.Draw("SAME");
         fitFuncBG.Draw("SAME");
         fitFuncGaus.Draw("SAME");
         PrintCanvas(&canv, outputDir + dValName + "_pT" + 
                     DtoStr(pTMin, 1) + "-" + DtoStr(pTMax, 1));
         
         means.AddPoint(distrDVal->GetYaxis()->GetBinCenter(j), fitFuncDVal.GetParameter(1));
         means.SetPointError(means.GetN() - 1, fitFuncDVal.GetParError(1));
         sigmas.AddPoint(distrDVal->GetYaxis()->GetBinCenter(j), fitFuncDVal.GetParameter(2));
         sigmas.SetPointError(sigmas.GetN() - 1, fitFuncDVal.GetParError(2));

         meansYMin = Minimum(means.GetPointY(means.GetN() - 1), meansYMin);
         meansYMax = Maximum(means.GetPointY(means.GetN() - 1), meansYMax);
         sigmasYMin = Minimum(sigmas.GetPointY(sigmas.GetN() - 1), sigmasYMin);
         sigmasYMax = Maximum(sigmas.GetPointY(sigmas.GetN() - 1), sigmasYMax);
      }

      if (means.GetN() == 0) continue;
      
      TCanvas meanCanv("means", "means", 600, 600);
      gPad->DrawFrame(0., meansYMin/1.2, 10., meansYMax*1.2);
      means.Draw("P");
      PrintCanvas(&meanCanv, outputDir + dValName + "_means");
      
      TCanvas sigmaCanv("sigmas", "sigmas", 600, 600);
      gPad->DrawFrame(0., 0., 10., sigmasYMax);
      sigmas.Draw("P");
      PrintCanvas(&sigmaCanv, outputDir + dValName + "_sigmas");
   }
}

#endif /* CALIBRATE_SIGMALIZED_RESIDUALS_CPP */
