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

   const Json::Value calibrationsInput = inputJSONContents["CALIBRATIONS"]["sigmalized_residuals"];
   
   if (calibrationsInput["detectors_to_calibrate"].size() == 0)
   {
      PrintInfo("No detectors are specified for calibrations");
      PrintInfo("Exiting the program");
      exit(1);
   }
   
   if (!inputJSONContents["is_pp"].asBool())
   {
      Par.centralityMin = 80.;//inputJSONContents["centrality_min"].asDouble();
      Par.centralityMax = 88.;//inputJSONContents["centrality_max"].asDouble();
      Par.centralityNBins = ceil((Par.centralityMax - Par.centralityMin - 1e-6)/5.);
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

   ProgressBar pBar("FANCY1");
   pBar.SetColor(PBarColor::BOLD_GREEN);
   
   const unsigned long numberOfIterations = 
      calibrationsInput["detectors_to_calibrate"].size()*Par.zDCMin.size()*4.;
   
   unsigned long nCalls = 0.;

   const std::string inputFileName = "data/Real/" + runName + "/SingleTrack/sum.root";

   for (const auto& detector : 
        calibrationsInput["detectors_to_calibrate"])
   {
      const std::string outputDir = "output/ResidualCal/" + runName + 
                                    "/" + detector["name"].asString() + "/";
      system(("mkdir -p " + outputDir).c_str());
   }

   TFile inputFile(inputFileName.c_str(), "READ");

   Par.texText.SetTextFont(52);
   Par.texText.SetTextSize(0.06);

   const std::string centralityRangeName = DtoStr(Par.centralityMin, 0) + "-" + 
                                           DtoStr(Par.centralityMax, 0) + "%"; 

   for (const auto& detector : 
        calibrationsInput["detectors_to_calibrate"])
   {
      const std::string outputDir = "output/ResidualCal/" + runName + "/";
      const std::string detectorName = detector["name"].asString();
      
      TFile outputFile((outputDir + detectorName + "/fits.root").c_str(), "RECREATE");
      
      for (const std::string& dValName : std::vector<std::string>{"dphi", "dz"})
      { 
         for (const bool isPositive : std::vector<bool>{true, false})
         {
            const std::string chargeName = (isPositive) ? "charge>0" : "charge<0";
            
            // histograms with weights representing means and sigmas; needed for quick inspection
            TH2D meansDistr("means", "#mu", 
                            Par.zDCMin.size(), Par.zDCMin.front(), Par.zDCMax.back(),
                            Par.pTMin.size(), Par.pTMin.front(), Par.pTMax.back());
            TH2D sigmasDistr("sigmas", "#sigma", 
                             Par.zDCMin.size(), Par.zDCMin.front(), Par.zDCMax.back(),
                             Par.pTMin.size(), Par.pTMin.front(), Par.pTMax.back());
            // histograms with weights representing the difference between means and sigmas and
            // a fit parameter of means and sigmas
            TH2D meansDiffDistr("means diff", "#cbar#mu - #mu_{fit}#cbar/#mu", 
                                Par.zDCMin.size(), Par.zDCMin.front(), Par.zDCMax.back(),
                                Par.pTMin.size(), Par.pTMin.front(), Par.pTMax.back());
            TH2D sigmasDiffDistr("sigmas diff", "#cbar#sigma - #sigma_{fit}#cbar/#sigma", 
                                 Par.zDCMin.size(), Par.zDCMin.front(), Par.zDCMax.back(),
                                 Par.pTMin.size(), Par.pTMin.front(), Par.pTMax.back());

            std::vector<double> pTRanges = Par.pTMin;
            pTRanges.push_back(Par.pTMax.back());

            meansDistr.GetYaxis()->Set(pTRanges.size() - 1, &pTRanges[0]);
            sigmasDistr.GetYaxis()->Set(pTRanges.size() - 1, &pTRanges[0]);
            meansDiffDistr.GetYaxis()->Set(pTRanges.size() - 1, &pTRanges[0]);
            sigmasDiffDistr.GetYaxis()->Set(pTRanges.size() - 1, &pTRanges[0]);
            
            std::vector<TGraphErrors> grVMeans, grVSigmas;
            std::vector<TF1> fVMeans, fVSigmas;
            
            for (unsigned long i = 0; i < Par.zDCMin.size(); i++)
            { 
               pBar.Print(static_cast<double>(nCalls)/static_cast<double>(numberOfIterations));
               nCalls++;
                
               const std::string zDCRangeName = DtoStr(Par.zDCMin[i], 0) + "<zDC<" + 
                                                DtoStr(Par.zDCMax[i], 0);
               // name of histogram
               const std::string distrDValName = dValName + " vs pT vs centrality: " + 
                                                 detectorName + ", " + chargeName + ", " + 
                                                 zDCRangeName;
               TH3F *distrDVal = static_cast<TH3F *>(inputFile.Get(distrDValName.c_str()));
               
               std::string fitsOutputFileName = outputDir + detectorName + "/" + dValName + 
                                                ((isPositive) ? "_pos" : "_neg") +
                                                "_c" + DtoStr(Par.centralityMin, 0) + 
                                                "-" + DtoStr(Par.centralityMax, 0) + 
                                                "_zDC" + DtoStr(Par.zDCMin[i], 0) + 
                                                "-" + DtoStr(Par.zDCMax[i], 0);

               grVMeans.emplace_back();
               grVSigmas.emplace_back();
               
               fVMeans.emplace_back((zDCRangeName + centralityRangeName + detectorName + 
                                     chargeName + dValName).c_str(), 
                                    Par.meansFitFunc.back().c_str());
               fVSigmas.emplace_back((zDCRangeName + centralityRangeName + detectorName + 
                                      chargeName + dValName).c_str(), 
                                     Par.sigmasFitFunc.back().c_str());
                              
               PerformFits(distrDVal, grVMeans.back(), grVSigmas.back(), 
                           fitsOutputFileName, dValName, detector["name"].asString(),
                           zDCRangeName + " [cm]", Par.centralityMin, Par.centralityMax, 
                           chargeName);

               fVMeans.back().SetRange(Par.pTMin.front()/1.05, Par.pTMax.back()*1.05);
               fVSigmas.back().SetRange(Par.pTMin.front()/1.05, Par.pTMax.back()*1.05);
               for (short j = 1; j <= Par.fitNTries; j++)
               {
                  grVMeans.back().Fit(&fVMeans.back(), "RQMN");
                  grVSigmas.back().Fit(&fVSigmas.back(), "RQMN");
                  for (int k = 0; k < fVMeans.back().GetNpar(); k++)
                  {
                     fVMeans.back().SetParLimits(k, fVMeans.back().GetParameter(k)/
                                                  (1. + 5./static_cast<double>(j*j)),
                                                  fVMeans.back().GetParameter(k)*
                                                  (1. + 5./static_cast<double>(j*j)));
                  }
                  for (int k = 0; k < fVSigmas.back().GetNpar(); k++)
                  {
                     fVSigmas.back().SetParLimits(k, fVSigmas.back().GetParameter(k)/
                                                   (1. + 5./static_cast<double>(j*j)),
                                                   fVSigmas.back().GetParameter(k)*
                                                   (1. + 5./static_cast<double>(j*j)));
                  }
               }

               for (int j = 0; j < grVMeans.back().GetN(); j++)
               {
                  const double x = grVMeans.back().GetPointX(j);
                  const int xBin = 
                     meansDistr.GetXaxis()->FindBin(Average(Par.zDCMin[i], Par.zDCMax[i]));
                  const int yBin = meansDistr.GetYaxis()->FindBin(x);
                  
                  meansDistr.SetBinContent(xBin, yBin, grVMeans.back().GetPointY(j));
                  sigmasDistr.SetBinContent(xBin, yBin, grVSigmas.back().GetPointY(j));
                  meansDiffDistr.SetBinContent(xBin, yBin, fabs(grVMeans.back().GetPointY(j) - 
                                                                fVMeans.back().Eval(x))/
                                               grVMeans.back().GetPointY(j));
                  sigmasDiffDistr.SetBinContent(xBin, yBin, fabs(grVSigmas.back().GetPointY(j) - 
                                                                 fVSigmas.back().Eval(x))/
                                                grVSigmas.back().GetPointY(j));
               }

               grVMeans.back().Clone()->Write(("Data means: " + dValName + 
                                               chargeName + ", " + centralityRangeName + 
                                               "%, " + zDCRangeName).c_str());
               
               grVSigmas.back().Clone()->Write(("Data sigmas: " + dValName + 
                                                chargeName + ", " + centralityRangeName +
                                                ", " + zDCRangeName).c_str());
               
               fVMeans.back().Clone()->Write(("Fit means: " + dValName + 
                                              chargeName + ", " + centralityRangeName +
                                              ", " + zDCRangeName).c_str());
               
               fVSigmas.back().Clone()->Write(("Fit sigmas: " + dValName + 
                                               chargeName + ", " + centralityRangeName +
                                               ", " + zDCRangeName).c_str());
            }

            double meanYMin = 1e31, meanYMax = -1e31;
            double sigmaYMax = -1e31;
            
            for (unsigned long i = 0; i < grVMeans.size(); i++)
            {
               double bVal = (dValName == "dphi") ? 0.04 : 10.; // bound value of dVal for graphs
               
               meanYMin = Maximum(-bVal, Minimum(meanYMin, TMath::MinElement(grVMeans[i].GetN(),
                                                                             grVMeans[i].GetY())));
               meanYMax = Minimum(bVal, Maximum(meanYMax, TMath::MaxElement(grVMeans[i].GetN(),
                                                                            grVMeans[i].GetY())));
               sigmaYMax = Minimum(bVal, Maximum(sigmaYMax, TMath::MaxElement(grVSigmas[i].GetN(), 
                                                                              grVSigmas[i].GetY())));

               grVMeans[i].SetMarkerStyle(Par.markerStyle[i]);
               grVMeans[i].SetMarkerSize(1.4);
               grVMeans[i].SetMarkerColorAlpha(Par.markerColor[i], 0.8);
               grVMeans[i].SetLineColorAlpha(Par.markerColor[i], 0.8);
               fVMeans[i].SetLineColorAlpha(Par.markerColor[i], 0.9);
               fVMeans[i].SetLineStyle(3);
               fVMeans[i].SetRange(Par.pTMin.front()/1.05, Par.pTMax.back()*1.05);

               grVSigmas[i].SetMarkerStyle(Par.markerStyle[i]);
               grVSigmas[i].SetMarkerSize(1.4);
               grVSigmas[i].SetMarkerColorAlpha(Par.markerColor[i], 0.8);
               grVSigmas[i].SetLineColorAlpha(Par.markerColor[i], 0.8);
               fVSigmas[i].SetLineColorAlpha(Par.markerColor[i], 0.9);
               fVSigmas[i].SetLineStyle(2);
               fVSigmas[i].SetRange(Par.pTMin.front()/1.05, Par.pTMax.back()*1.05);
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
               const std::string zDCRangeName = DtoStr(Par.zDCMin[i], 0) + "<z_{DC}<" + 
                                                DtoStr(Par.zDCMax[i], 0);
               
               legend.AddEntry(&grVMeans[i], zDCRangeName.c_str(), "P");
               grVMeans[i].Clone()->Draw("SAME P");
               fVMeans[i].Clone()->Draw("SAME");
            }

            legend.DrawClone();
            PrintCanvas(&canv, outputDir + detectorName + "_means_" + 
                         dValName + ((isPositive) ? "_pos" : "_neg") +
                        "_c" + DtoStr(Par.centralityMin, 0) + "-" + DtoStr(Par.centralityMax, 0));

            canv.Clone()->Write((dValName + " means: " + 
                                 chargeName + ", " + 
                                 DtoStr(Par.centralityMin, 0) + "-" + 
                                 DtoStr(Par.centralityMax, 0) + "%").c_str());
            
            legend.Clear();
            canv.Clear();
            
            TH1 *sigmasFrame = 
               gPad->DrawFrame(0., 0., Par.pTMax.back()*1.05, sigmaYMax*1.6);
            sigmasFrame->Draw("SAME AXIS X+ Y+");
            
            for (unsigned long i = 0; i < grVSigmas.size(); i++)
            { 
               const std::string zDCRangeName = DtoStr(Par.zDCMin[i], 0) + "<z_{DC}<" + 
                                                DtoStr(Par.zDCMax[i], 0);
               
               legend.AddEntry(&grVSigmas[i], zDCRangeName.c_str(), "P");
               grVSigmas[i].Clone()->Draw("SAME P");
               fVSigmas[i].Clone()->Draw("SAME");
            }

            legend.DrawClone();
            PrintCanvas(&canv, outputDir + detectorName + "_sigmas_" + 
                        dValName + ((isPositive) ? "_pos" : "_neg") +
                        "_c" + DtoStr(Par.centralityMin, 0) + "-" + DtoStr(Par.centralityMax, 0));

            canv.Clone()->Write((dValName + " sigmas: " + 
                                 ((isPositive) ? "pos, " : "neg, ") + 
                                 DtoStr(Par.centralityMin, 0) + "-" + 
                                 DtoStr(Par.centralityMax, 0) + "%").c_str());

            TCanvas parCanv("", "", 800, 800);
            parCanv.Divide(2, 2);
            
            parCanv.cd(1);
            gPad->SetLogz();
            gPad->SetRightMargin(0.13);
            meansDistr.GetXaxis()->SetTitle("z_{DC}");
            meansDistr.GetYaxis()->SetTitle("p_{T}");
            meansDistr.Draw("COLZ");

            parCanv.cd(2);
            gPad->SetLogz();
            gPad->SetRightMargin(0.13);
            sigmasDistr.GetXaxis()->SetTitle("z_{DC}");
            sigmasDistr.GetYaxis()->SetTitle("p_{T}");
            sigmasDistr.Draw("COLZ");

            parCanv.cd(3);
            gPad->SetLogz();
            gPad->SetRightMargin(0.13);
            meansDiffDistr.GetXaxis()->SetTitle("z_{DC}");
            meansDiffDistr.GetYaxis()->SetTitle("p_{T}");
            meansDiffDistr.Draw("COLZ");

            parCanv.cd(4);
            gPad->SetLogz();
            gPad->SetRightMargin(0.13);
            sigmasDiffDistr.GetXaxis()->SetTitle("z_{DC}");
            sigmasDiffDistr.GetYaxis()->SetTitle("p_{T}");
            sigmasDiffDistr.Draw("COLZ");
            
            PrintCanvas(&parCanv, outputDir + detectorName + "/fitPar_" + dValName +  + 
                        ((isPositive) ? "_pos" : "_neg") + 
                        "_c" + DtoStr(Par.centralityMin, 0) + "-" + 
                        DtoStr(Par.centralityMax, 0));
         }
      }
      outputFile.Close();
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
 
   auto PerformFitsInRange = [&](const double pTMin, const double pTMax, 
                                 TF1 *meansFit, TF1 *sigmasFit,
                                 TF1 *limitMeansFitFunc = NULL, TF1 *limitSigmasFitFunc = NULL,
                                 bool limitParByFit = false)
   {
      for (int i = grMeans.GetN() - 1; i >= 0; i--)
      {
         grMeans.RemovePoint(i);
         grSigmas.RemovePoint(i);
      }

      int iCanv = 1;
      
      for (unsigned long i = 0; i < Par.pTMin.size(); i++)
      {
         const double pT = Average(Par.pTMin[i], Par.pTMax[i]);
         if (Par.pTMin[i] < pTMin || Par.pTMax[i] > pTMax) continue;
         
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
         TF1 fitFuncDVal("fitFunc", "gaus(0) + gaus(3)");
         TF1 fitFuncBG("bg", "gaus");
         fitFuncGaus.SetParameters(1., 0., binWidth*2.);
         fitFuncDVal.SetParameters(1., 0., binWidth*2., 1., 0., maxX/2.);
         
         if (!limitParByFit)
         {
            fitFuncGaus.SetParLimits(1, minBinX/5., maxBinX/5.);
            fitFuncGaus.SetParLimits(2, binWidth, maxBinX/5.);
            fitFuncDVal.SetParLimits(1, minX/10., maxX/10.);
            fitFuncDVal.SetParLimits(2, binWidth, Average(maxX, maxX, minX));
         }
         else
         {
            fitFuncGaus.SetParLimits(1, limitMeansFitFunc->Eval(pT)/1.5, 
                                     limitMeansFitFunc->Eval(pT)*1.5);
            fitFuncGaus.SetParLimits(2, limitSigmasFitFunc->Eval(pT)/1.5, 
                                     limitSigmasFitFunc->Eval(pT)*1.5);
            fitFuncDVal.SetParLimits(1, limitMeansFitFunc->Eval(pT)/1.5, 
                                     limitMeansFitFunc->Eval(pT)*1.5);
            fitFuncDVal.SetParLimits(2, limitSigmasFitFunc->Eval(pT)/1.5, 
                                     limitSigmasFitFunc->Eval(pT)*1.5);
         }
         
         fitFuncDVal.SetParLimits(4, minX*2., maxX*2.);
         fitFuncDVal.SetParLimits(5, maxX/3., maxX*3.);
         
         fitFuncDVal.SetLineColorAlpha(kRed+1, 0.6);
         fitFuncBG.SetLineColorAlpha(kGreen+1, 0.9);
         fitFuncBG.SetLineStyle(2);
         fitFuncGaus.SetLineColorAlpha(kAzure-3, 0.9);
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
         
         // scale limits
         fitFuncGaus.SetParLimits(0, maxBinVal/5., maxBinVal);
         fitFuncDVal.SetParLimits(0, maxBinVal/5., maxBinVal);
         fitFuncDVal.SetParLimits(3, maxBinVal/20., maxBinVal);
         
         fitFuncGaus.SetRange(minBinX/5., maxBinX/5.);
         
         distrDValProj->Fit(&fitFuncGaus, "RQMBN");

         fitFuncGaus.SetRange(minBinX, maxBinX);
         fitFuncBG.SetRange(minBinX, maxBinX);
         
         for (int k = 0; k < 3; k++)
         {
            fitFuncDVal.SetParameter(k, fitFuncGaus.GetParameter(k));
         }
         
         fitFuncDVal.SetRange(minBinX, maxBinX);
         distrDValProj->Fit(&fitFuncDVal, "RQMBN");

         for (unsigned short j = 1; j <= Par.fitNTries; j++)
         {
            for (int k = 0; k < fitFuncDVal.GetNpar(); k++)
            {
               fitFuncDVal.SetParLimits(k, fitFuncDVal.GetParameter(k)/
                                        (1. + 1./static_cast<double>(j)),
                                        fitFuncDVal.GetParameter(k)*
                                        (1. + 1./static_cast<double>(j)));
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
         
         grMeans.AddPoint(pT, fitFuncDVal.GetParameter(1));
         grMeans.SetPointError(grMeans.GetN() - 1, fitFuncDVal.GetParError(1));
         grSigmas.AddPoint(Average(Par.pTMin[i], Par.pTMax[i]), fitFuncDVal.GetParameter(2));
         grSigmas.SetPointError(grSigmas.GetN() - 1, fitFuncDVal.GetParError(2));
      }
      if (grMeans.GetN() == 0) 
      {
         PrintError("Graph is empty for " + dValName + ", " + detectorName + ", " + 
                    chargeName + " at " + zDCRangeName + ", " + centralityRangeName);
      }
      
      meansFit->SetRange(pTMin/1.05, pTMax*1.05);
      
      grMeans.Fit(meansFit, "RQMN");
      grSigmas.Fit(sigmasFit, "RQMN");
   };
   
   TF1 meansFit("meansFit", Par.meansFitFunc.front().c_str());
   TF1 sigmasFit("sigmasFit", Par.sigmasFitFunc.front().c_str());
   
   PerformFitsInRange(Par.pTMinFit.front(), Par.pTMaxFit.front(), &meansFit, &sigmasFit);

   for (unsigned long i = 1; i < Par.pTMinFit.size(); i++)
   {
      PerformFitsInRange(Par.pTMinFit[i], Par.pTMaxFit[i], &meansFit, &sigmasFit,
                         (TF1 *) meansFit.Clone(), (TF1 *) sigmasFit.Clone(), true);
   }
   for (unsigned long i = 1; i < Par.meansFitFunc.size(); i++)
   {
      TF1 *oldMeansFit = (TF1 *) meansFit.Clone();
      TF1 *oldSigmasFit = (TF1 *) sigmasFit.Clone();

      meansFit.SetCurrent((TF1 *) TF1("meansFit", Par.meansFitFunc[i].c_str()).Clone());
      sigmasFit.SetCurrent((TF1 *) TF1("sigmasFit", Par.sigmasFitFunc[i].c_str()).Clone());

      PerformFitsInRange(Par.pTMinFit.back(), Par.pTMaxFit.back(), &meansFit, &sigmasFit,
                         oldMeansFit, oldSigmasFit, true);
   }
   
   PrintCanvas(&canv, outputFileNameNoExt, false, true);
   /*
   canv.Clone()->Write((dValName + ": " + chargeName + ", " +
                       centralityRangeName + ", " + zDCRangeName).c_str());
   */
}

#endif /* CALIBRATE_SIGMALIZED_RESIDUALS_CPP */
