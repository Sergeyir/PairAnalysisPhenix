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
   
   InputJSONReader inputJSONContents{argv[1], "main"};
   inputJSONContents.CheckStatus("main");
   
   const std::string runName = inputJSONContents["run_name"].asString();

   if (inputJSONContents["uncalibrated_sigmalized_residuals_detectors"].size() == 0)
   {
      PrintInfo("No detectors are specified for calibrations");
      PrintInfo("Exiting the program");
      exit(1);
   }

   if (!inputJSONContents["is_pp"].asBool())
   {
      Par.centralityMin = inputJSONContents["centrality_min"].asDouble();
      Par.centralityMax = inputJSONContents["centrality_max"].asDouble();
      Par.centralityNBins = ceil((Par.centralityMax - Par.centralityMin)/5.);
   }
   else
   {
      Par.centralityMin = 0.;
      Par.centralityMax = 100.;
      Par.centralityNBins = 1.;
   }

   ROOT::EnableImplicitMT();
   //ROOT::EnableThreadSafety();
   
   gErrorIgnoreLevel = kWarning;
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);

   ProgressBar pBar("WAVE");
   
   const unsigned long numberOfIterations = 
      inputJSONContents["uncalibrated_sigmalized_residuals_detectors"].size()*
      Par.zDCMin.size()*4.;
   
   unsigned long nCalls = 0.;

   const std::string inputFileName = "data/Real/" + runName + "/SingleTrack/sum.root";

   for (const auto& detector : 
        inputJSONContents["uncalibrated_sigmalized_residuals_detectors"])
   {
      const std::string outputDir = "output/ResidualCal/" + runName + 
                                    "/" + detector["name"].asString() + "/";
      system(("mkdir -p " + outputDir).c_str());
   }

   TFile inputFile(inputFileName.c_str(), "READ");

   Par.texText.SetTextFont(52);
   Par.texText.SetTextSize(0.06);

   for (const auto& detector : 
        inputJSONContents["uncalibrated_sigmalized_residuals_detectors"])
   {
      const std::string outputDir = "output/ResidualCal/" + runName + 
                                    "/" + detector["name"].asString() + "/";
      
      for (const std::string& dValName : std::vector<std::string>{"dphi", "dz"})
      {
         for (const bool isPositive : std::vector<bool>{true, false})
         {
            std::vector<TGraphErrors> grVMeans, grVSigmas;
            std::vector<TF1> fVMeans, fVSigmas;
            
            for (unsigned long i = 0; i < Par.zDCMin.size(); i++)
            { 
               pBar.Print(static_cast<double>(nCalls)/static_cast<double>(numberOfIterations));
               nCalls++;
                
               const std::string zDCRangeName = DtoStr(Par.zDCMin[i], 0) + "<zDC<" + 
                                                DtoStr(Par.zDCMax[i], 0);
               const std::string centralityRangeName = DtoStr(Par.centralityMin, 0) + "-" + 
                                                       DtoStr(Par.centralityMax, 0) + "%";
               const std::string detectorName = detector["name"].asString();
               const std::string chargeName = (isPositive) ? "charge>0" : "charge<0";
               
               // name of histogram
               const std::string distrDValName = dValName + " vs pT vs centrality: " + 
                                                 detectorName + ", " + chargeName + ", " + 
                                                 zDCRangeName;
               TH3F *distrDVal = static_cast<TH3F *>(inputFile.Get(distrDValName.c_str()));
               
               std::string fitsOutputFileName = outputDir + dValName + 
                                                ((isPositive) ? "_pos" : "_neg") +
                                                "_c" + DtoStr(Par.centralityMin, 0) + 
                                                "-" + DtoStr(Par.centralityMax, 0) + 
                                                "_zDC" + DtoStr(Par.zDCMin[i], 0) + 
                                                "-" + DtoStr(Par.zDCMax[i], 0);

               grVMeans.emplace_back();
               grVSigmas.emplace_back();
               
               fVMeans.emplace_back((zDCRangeName + centralityRangeName + detectorName + 
                                     chargeName + dValName).c_str(), Par.meansFitFunc.c_str());
               fVSigmas.emplace_back((zDCRangeName + centralityRangeName + detectorName + 
                                      chargeName + dValName).c_str(), Par.sigmasFitFunc.c_str());
               
               fVMeans.back().SetRange(Par.pTMin.front()/1.05, Par.pTMax.back()*1.05);
               fVSigmas.back().SetRange(Par.pTMin.front()/1.05, Par.pTMax.back()*1.05);
               
               PerformFits(distrDVal, grVMeans.back(), grVSigmas.back(), 
                           fitsOutputFileName, dValName, detector["name"].asString(),
                           zDCRangeName + " [cm]", Par.centralityMin, Par.centralityMax, 
                           chargeName);

               grVMeans.back().Fit(&fVMeans.back(), "RQMN");
               grVSigmas.back().Fit(&fVSigmas.back(), "RQMN");
            }

            double meanYMin = 1e31, meanYMax = -1e31;
            double sigmaYMax = -1e31;
            
            for (unsigned long i = 0; i < grVMeans.size(); i++)
            {
               meanYMin =
                  Minimum(meanYMin, TMath::MinElement(grVMeans[i].GetN(), grVMeans[i].GetY()));
               meanYMax = 
                  Maximum(meanYMax, TMath::MaxElement(grVMeans[i].GetN(), grVMeans[i].GetY()));
               sigmaYMax = 
                  Maximum(sigmaYMax, TMath::MaxElement(grVSigmas[i].GetN(), grVSigmas[i].GetY()));
            }
            
            TCanvas canv("", "", 800, 800);
            
            TLegend legend{0.1, 0.7, 0.88, 0.89};
            legend.SetNColumns(3);
            legend.SetLineColorAlpha(0, 0.);
            legend.SetFillColorAlpha(0, 0.);
            
            TH1 *meansFrame = 
               gPad->DrawFrame(0., meanYMin - (meanYMax - meanYMin)*0.05, 
                               Par.pTMax.back()*1.05, meanYMax + (meanYMax - meanYMin)*0.35);
            
            meansFrame->Draw("SAME AXIS X+ Y+");
            
            for (unsigned long i = 0; i < grVMeans.size(); i++)
            {
               grVMeans[i].SetMarkerStyle(Par.markerStyle[i]);
               grVMeans[i].SetMarkerSize(1.4);
               grVMeans[i].SetMarkerColorAlpha(Par.markerColor[i], 0.8);
               grVMeans[i].SetLineColorAlpha(Par.markerColor[i], 0.8);
               fVMeans[i].SetLineColorAlpha(Par.markerColor[i], 0.9);
               fVMeans[i].SetLineStyle(2);

               const std::string zDCRangeName = DtoStr(Par.zDCMin[i], 0) + "<z_{DC}<" + 
                                                DtoStr(Par.zDCMax[i], 0);
               
               legend.AddEntry(&grVMeans[i], zDCRangeName.c_str(), "P");
               grVMeans[i].Clone()->Draw("SAME P");
               fVMeans[i].Clone()->Draw("SAME");
            }

            legend.DrawClone();
            PrintCanvas(&canv, outputDir + "means_" + dValName + ((isPositive) ? "_pos" : "_neg") +
                        "_c" + DtoStr(Par.centralityMin, 0) + "-" + DtoStr(Par.centralityMax, 0));
            
            legend.Clear();
            canv.Clear();
            
            TH1 *sigmasFrame = 
               gPad->DrawFrame(0., 0., Par.pTMax.back()*1.05, sigmaYMax*1.6);
            sigmasFrame->Draw("SAME AXIS X+ Y+");
            
            for (unsigned long i = 0; i < grVSigmas.size(); i++)
            {
               grVSigmas[i].SetMarkerStyle(Par.markerStyle[i]);
               grVSigmas[i].SetMarkerSize(1.4);
               grVSigmas[i].SetMarkerColorAlpha(Par.markerColor[i], 0.8);
               grVSigmas[i].SetLineColorAlpha(Par.markerColor[i], 0.8);
               fVSigmas[i].SetLineColorAlpha(Par.markerColor[i], 0.9);
               fVSigmas[i].SetLineStyle(2);
               
               const std::string zDCRangeName = DtoStr(Par.zDCMin[i], 0) + "<z_{DC}<" + 
                                                DtoStr(Par.zDCMax[i], 0);
               
               legend.AddEntry(&grVSigmas[i], zDCRangeName.c_str(), "P");
               grVSigmas[i].Clone()->Draw("SAME P");
               fVSigmas[i].Clone()->Draw("SAME");
            }

            legend.DrawClone();
            PrintCanvas(&canv, outputDir + "sigmas_" + dValName + ((isPositive) ? "_pos" : "_neg") +
                        "_c" + DtoStr(Par.centralityMin, 0) + "-" + DtoStr(Par.centralityMax, 0));
         }
      }
   }
}
   

void PerformFits(TH3F *hist, TGraphErrors &grMeans, TGraphErrors &grSigmas,
                 const std::string& outputFileNameNoExt, const std::string& dValName,
                 const std::string& detectorName, const std::string& zDCRangeName,
                 const double centralityMin, const double centralityMax,
                 const std::string& chargeName)
{ 
   const double minBinX = hist->GetXaxis()->GetBinLowEdge(1);
   const double maxBinX = hist->GetXaxis()->GetBinUpEdge(hist->GetXaxis()->GetNbins());
   const double binWidth = hist->GetXaxis()->GetBinWidth(1);

   const std::string centralityRangeName = DtoStr(centralityMin, 0) + "-" + 
                                          DtoStr(centralityMax, 0) + "%";

   TCanvas canv("dval vs pT", "", Par.pTXNBins*400, Par.pTYNBins*400);
   canv.Divide(Par.pTXNBins, Par.pTYNBins);

   int iCanv = 1;
    
   for (unsigned long i = 0; i < Par.pTMin.size(); i++)
   {
      canv.cd(iCanv);
      iCanv++;
   
      TH1D *distrDValProj = hist->
         ProjectionX(((std::string) hist->GetName() + "_projX").c_str(), 
                     hist->GetYaxis()->FindBin(Par.pTMin[i] + 1e-6), 
                     hist->GetYaxis()->FindBin(Par.pTMax[i] - 1e-6),
                     hist->GetZaxis()->FindBin(centralityMin + 1e-6),
                     hist->GetZaxis()->FindBin(centralityMax - 1e-6));


      const std::string pTRangeName = DtoStr(Par.pTMin[i], 1) + "<pT<" + 
                                      DtoStr(Par.pTMax[i], 1) + " [GeV/c]";
      
      if (distrDValProj->Integral(1, distrDValProj->GetXaxis()->GetNbins()) < 
          Par.minIntegralValue) 
      {
         PrintInfo("Integral is insufficient for projection of " + dValName + ", " + 
                   detectorName + ", " + chargeName + " at " + zDCRangeName + ", " + 
                   centralityRangeName + ", " + pTRangeName);
         continue;
      }
      
      double minX = 0., maxX = -1.;
      for (int k = 1; k <= distrDValProj->GetXaxis()->GetNbins(); k++)
      {
         if (distrDValProj->GetBinContent(k) > 1e-7)
         {
            minX = distrDValProj->GetXaxis()->GetBinLowEdge(k);
            break;
         }
      }
      
      for (int k = distrDValProj->GetXaxis()->GetNbins(); k > 1; k--)
      {
         if (distrDValProj->GetBinContent(k) > 1e-7)
         {
            maxX = distrDValProj->GetXaxis()->GetBinUpEdge(k);
            break;
         }
      }
      
      if (minX > maxX) 
      {
         PrintWarning("Something wrong for projection of " + dValName + ", " + 
                      detectorName + ", " + chargeName + " at " + zDCRangeName + ", " + 
                      centralityRangeName + ", " + pTRangeName);
         continue;
      }

      TF1 fitFuncGaus("gaus", "gaus");
      fitFuncGaus.SetParameters(1., 0., binWidth*2.);
      fitFuncGaus.SetParLimits(1, minBinX/5., maxBinX/5.);
      fitFuncGaus.SetParLimits(2, binWidth, maxBinX/5.);
      
      TF1 fitFuncDVal("fitFunc", "gaus(0) + gaus(3)");
      fitFuncDVal.SetParameters(1., 0., binWidth*2.);
      fitFuncDVal.SetParLimits(1, minX/10., maxX/10.);
      fitFuncDVal.SetParLimits(2, binWidth, Average(maxX, maxX, minX));
      fitFuncDVal.SetParLimits(4, minX*2., maxX*2.);
      fitFuncDVal.SetParLimits(5, maxX/3., maxX*2.);
      
      TF1 fitFuncBG("bg", "gaus");
      fitFuncBG.SetParameters(1., 0., Average(maxBinX, maxBinX, maxBinX, minBinX));
      fitFuncBG.SetParLimits(1, minX*2., maxX*2.);
      fitFuncBG.SetParLimits(2, maxX/3., maxX*2.);
      
      fitFuncDVal.SetLineColorAlpha(kRed+1, 0.6);
      //fitFuncDVal.SetLineWidth(2);
      fitFuncBG.SetLineColorAlpha(kGreen+1, 0.9);
      //fitFuncBG.SetLineWidth(3);
      fitFuncBG.SetLineStyle(2);
      fitFuncGaus.SetLineColorAlpha(kAzure-3, 0.9);
      //fitFuncGaus.SetLineWidth(3);
      fitFuncGaus.SetLineStyle(2);

      distrDValProj->GetXaxis()->SetTitle(dValName.c_str());
      distrDValProj->SetTitle("");
      distrDValProj->SetTitleSize(0.06, "X");
      distrDValProj->SetTitleSize(0.06, "Y");
      distrDValProj->SetLabelSize(0.06, "X");
      distrDValProj->SetLabelSize(0.06, "Y");
 
      distrDValProj->GetXaxis()->SetRange(distrDValProj->GetXaxis()->FindBin(minX+0.01),
                                          distrDValProj->GetXaxis()->FindBin(maxX-0.01));
      
      const double maxBinVal = distrDValProj->GetBinContent(distrDValProj->GetMaximumBin());
      
      //const int expMean = distrDValProj->GetXaxis()->GetBinCenter(distrDValProj->GetMaximumBin());
      fitFuncGaus.SetParLimits(0, maxBinVal/5., maxBinVal);
      fitFuncDVal.SetParLimits(0, maxBinVal/5., maxBinVal);
      fitFuncDVal.SetParLimits(4, 0., maxBinVal);
      fitFuncBG.SetParLimits(0, 0., maxBinVal);
      
      fitFuncGaus.SetRange(minBinX/5., maxBinX/5.);
      fitFuncBG.SetRange(minBinX, maxBinX);
      
      distrDValProj->Fit(&fitFuncGaus, "RQMBN");
      distrDValProj->Fit(&fitFuncBG, "RQMBN");

      fitFuncGaus.SetRange(minBinX, maxBinX);
      
      for (int k = 0; k < 3; k++)
      {
         fitFuncDVal.SetParameter(k, fitFuncGaus.GetParameter(k));
         fitFuncDVal.SetParameter(k + 3, fitFuncBG.GetParameter(k));
      }
      
      fitFuncDVal.SetRange(minBinX, maxBinX);
      distrDValProj->Fit(&fitFuncDVal, "RQMBN");

      for (unsigned short j = 1; j <= Par.fitNTries; j++)
      {
         for (int k = 0; k < fitFuncDVal.GetNpar(); k++)
         {
            fitFuncDVal.SetParLimits(k, fitFuncDVal.GetParameter(k)/(1. + 1./static_cast<double>(j)),
                                     fitFuncDVal.GetParameter(k)*(1. + 1./static_cast<double>(j)));
         }
         
         distrDValProj->Fit(&fitFuncDVal, "RQMBN");
      }

      for (int j = 0; j < 3; j++)
      {
         fitFuncGaus.SetParameter(j, fitFuncDVal.GetParameter(j));
         fitFuncBG.SetParameter(j, fitFuncDVal.GetParameter(j + 3));
      }
      
      distrDValProj->SetMarkerStyle(20);
      distrDValProj->SetMarkerSize(0.7);
      distrDValProj->SetMarkerColorAlpha(kBlack, 0.8);
      distrDValProj->SetLineColorAlpha(kBlack, 0.8);
      distrDValProj->SetMaximum(maxBinVal*1.2);
      
      gPad->SetLeftMargin(0.155);
      gPad->SetBottomMargin(0.115);
      
      distrDValProj->Clone()->Draw("P");
      fitFuncDVal.Clone()->Draw("SAME");
      fitFuncBG.Clone()->Draw("SAME");
      fitFuncGaus.Clone()->Draw("SAME");
      
      Par.texText.DrawLatexNDC(0.17, 0.85, pTRangeName.c_str());
      Par.texText.DrawLatexNDC(0.17, 0.79, zDCRangeName.c_str());
      Par.texText.DrawLatexNDC(0.17, 0.73, chargeName.c_str());
      Par.texText.DrawLatexNDC(0.17, 0.66, centralityRangeName.c_str());
      
      grMeans.AddPoint(Average(Par.pTMin[i], Par.pTMax[i]), fitFuncDVal.GetParameter(1));
      grMeans.SetPointError(grMeans.GetN() - 1, fitFuncDVal.GetParError(1));
      grSigmas.AddPoint(Average(Par.pTMin[i], Par.pTMax[i]), fitFuncDVal.GetParameter(2));
      grSigmas.SetPointError(grSigmas.GetN() - 1, fitFuncDVal.GetParError(2));
   }
   
   if (grMeans.GetN() == 0) 
   {
      PrintError("Graph is empty for " + dValName + ", " + detectorName + ", " + 
                 chargeName + " at " + zDCRangeName + ", " + centralityRangeName);
   }
   
   PrintCanvas(&canv, outputFileNameNoExt, false, true);
}

#endif /* CALIBRATE_SIGMALIZED_RESIDUALS_CPP */
