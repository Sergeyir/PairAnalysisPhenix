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
   
   // opening input file with parameters for calibrations
   Par.inputJSONCal.OpenFile(argv[1], "calibration");
   Par.inputJSONCal.CheckStatus("calibration");
   
   Par.runName = Par.inputJSONCal["run_name"].asString();
   
   // opening input file with parameters of a run
   Par.inputJSONMain.OpenFile("input/" + Par.runName + "/main.json");
   Par.inputJSONMain.CheckStatus("main");

   const Json::Value calibrationInput = Par.inputJSONCal["SIGMALIZED_RESIDUALS"];
   
   if (calibrationInput["detectors_to_calibrate"].size() == 0)
   {
      PrintInfo("No detectors are specified for calibrations");
      PrintInfo("Exiting the program");
      exit(1);
   }

   ROOT::EnableImplicitMT();
   gErrorIgnoreLevel = kWarning;
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   
   const double pTMin = calibrationInput["pt_bins"].front()["min"].asDouble();
   const double pTMax = calibrationInput["pt_bins"].back()["max"].asDouble();
   
   ProgressBar pBar("FANCY1");
   pBar.SetColor(PBarColor::BOLD_GREEN);
   
   const unsigned long numberOfIterations = calibrationInput["detectors_to_calibrate"].size()*
                                            calibrationInput["centrality_bins"].size()*
                                            calibrationInput["variables_to_calibrate"].size()*
                                            calibrationInput["particles_to_calibrate"].size()*
                                            calibrationInput["zdc_bins"].size();
   
   unsigned long numberOfCalls = 0.;

   const std::string inputFileName = "data/Real/" + Par.runName + "/SingleTrack/sum.root";

   for (const auto& detector : calibrationInput["detectors_to_calibrate"])
   {
      const std::string outputDir = "output/ResidualCal/" + Par.runName + 
                                    "/" + detector["name"].asString() + "/";
      system(("mkdir -p " + outputDir).c_str());
   }

   TFile inputFile(inputFileName.c_str(), "READ");

   Par.texText.SetTextFont(52);
   Par.texText.SetTextSize(0.06);

   for (const Json::Value& detector : calibrationInput["detectors_to_calibrate"])
   {
      const std::string outputDir = "output/ResidualCal/" + Par.runName + "/";
      const std::string detectorName = detector["name"].asString();
      
      TFile outputFile((outputDir + detectorName + "/fits.root").c_str(), "RECREATE");

      for (const Json::Value& centralityBin : calibrationInput["centrality_bins"])
      {
         const std::string centralityRangeName = centralityBin["min"].asString() + "-" + 
                                                 centralityBin["max"].asString() + "%"; 
         // same as before but without percent; used for output file names
         const std::string centralityRangePathName = "_c" + centralityBin["min"].asString() + 
                                                       "-" + centralityBin["max"].asString(); 
         
         for (const Json::Value& variable : calibrationInput["variables_to_calibrate"])
         { 
            for (const Json::Value& particleType : calibrationInput["particles_to_calibrate"])
            {
               // histograms with weights representing means and sigmas
               TH2D meansDistr("means", "#mu", 
                               calibrationInput["zdc_bins"].size(), 
                               calibrationInput["zdc_bins"].front()["min"].asDouble(), 
                               calibrationInput["zdc_bins"].back()["max"].asDouble(),
                               1, 0., 1.);
               TH2D sigmasDistr("sigmas", "#sigma", 
                                calibrationInput["zdc_bins"].size(), 
                                calibrationInput["zdc_bins"].front()["min"].asDouble(), 
                                calibrationInput["zdc_bins"].back()["max"].asDouble(),
                                1, 0., 1.);
               // histograms with weights representing the difference between means and sigmas and
               // a fit parameter of means and sigmas
               TH2D meansDiffDistr("means diff", "#cbar#mu - #mu_{fit}#cbar/#mu", 
                                   calibrationInput["zdc_bins"].size(), 
                                   calibrationInput["zdc_bins"].front()["min"].asDouble(), 
                                   calibrationInput["zdc_bins"].back()["max"].asDouble(),
                                   1, 0., 1.);
               TH2D sigmasDiffDistr("sigmas diff", "#cbar#sigma - #sigma_{fit}#cbar/#sigma", 
                                   calibrationInput["zdc_bins"].size(), 
                                   calibrationInput["zdc_bins"].front()["min"].asDouble(), 
                                   calibrationInput["zdc_bins"].back()["max"].asDouble(),
                                   1, 0., 1.);

               std::vector<double> pTRanges;
               for (const auto & val: calibrationInput["pt_bins"])
               {
                  pTRanges.push_back(val["min"].asDouble());
               }
               pTRanges.push_back(pTMax);

               meansDistr.GetYaxis()->Set(pTRanges.size() - 1, &pTRanges[0]);
               sigmasDistr.GetYaxis()->Set(pTRanges.size() - 1, &pTRanges[0]);
               meansDiffDistr.GetYaxis()->Set(pTRanges.size() - 1, &pTRanges[0]);
               sigmasDiffDistr.GetYaxis()->Set(pTRanges.size() - 1, &pTRanges[0]);
               
               std::vector<TGraphErrors> grVMeans, grVSigmas;
               std::vector<TF1> fVMeans, fVSigmas;
               
               for (const Json::Value& zDCBin : calibrationInput["zdc_bins"])
               { 
                  pBar.Print(static_cast<double>(numberOfCalls)/
                             static_cast<double>(numberOfIterations));
                  numberOfCalls++;
                   
                  const std::string zDCRangeName = zDCBin["min"].asString() + "<zDC<" + 
                                                   zDCBin["max"].asString();
                  const std::string zDCRangePathName = "_zDC" + zDCBin["min"].asString() + 
                                                       "-" + zDCBin["max"].asString();
                  
                  const double zDCMin = zDCBin["min"].asDouble();
                  const double zDCMax = zDCBin["max"].asDouble();
                  
                  // name of histogram
                  const std::string distrVariableName =  
                     variable["name"].asString() + " vs pT vs centrality: " + 
                     detector["name"].asString() + ", " + particleType["name"].asString() + ", " + 
                     zDCRangeName;
                  
                  TH3F *distrVariable = 
                     static_cast<TH3F *>(inputFile.Get(distrVariableName.c_str()));
                  
                  if (!distrVariable) PrintError("Histogram named \"" + distrVariableName + 
                                                 "\" does not exist in file " + inputFileName);
                  
                  std::string fitsOutputFileName = outputDir + detector["name"].asString() + "/" + 
                                                   variable["name"].asString() + "_" +
                                                   particleType["name_short"].asString() +
                                                   centralityRangePathName + zDCRangePathName;

                  grVMeans.emplace_back();
                  grVSigmas.emplace_back();
                  
                  fVMeans.emplace_back((zDCRangeName + centralityRangeName + 
                                        detector["name"].asString() + 
                                        particleType["name"].asString() + 
                                        variable["name"].asString()).c_str(), 
                                       detector["means_fit"].asCString());
                  fVSigmas.emplace_back((zDCRangeName + centralityRangeName + 
                                         detector["name"].asString() + 
                                         particleType["name"].asString() +
                                         variable["name"].asString()).c_str(), 
                                        detector["sigmas_fit"].asCString());
                                 
                  PerformFits(distrVariable, grVMeans.back(), grVSigmas.back(), calibrationInput, 
                              detector, variable, zDCBin, particleType, centralityBin);
                  
                  fVMeans.back().SetRange(pTMin/1.05, pTMax*1.05);
                  fVSigmas.back().SetRange(pTMin/1.05, pTMax*1.05);
                  
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

                  // filling 2D histograms with weights as fit parameters means and sigmas
                  for (int j = 0; j < grVMeans.back().GetN(); j++)
                  {
                     const double x = grVMeans.back().GetPointX(j);
                     const int xBin = meansDistr.GetXaxis()->FindBin(Average(zDCMin, zDCMax));
                     const int yBin = meansDistr.GetYaxis()->FindBin(x);
                     
                     meansDistr.SetBinContent(xBin, yBin, grVMeans.back().GetPointY(j));
                     sigmasDistr.SetBinContent(xBin, yBin, grVSigmas.back().GetPointY(j));
                     meansDiffDistr.SetBinContent(xBin, yBin, 
                                                  fabs((grVMeans.back().GetPointY(j) - 
                                                        fVMeans.back().Eval(x))/
                                                       grVMeans.back().GetPointY(j)));
                     sigmasDiffDistr.SetBinContent(xBin, yBin, fabs(grVSigmas.back().GetPointY(j) - 
                                                                    fVSigmas.back().Eval(x))/
                                                   grVSigmas.back().GetPointY(j));
                  }

                  grVMeans.back().Clone()->
                     Write(("Data means: " + variable["name"].asString() + 
                            ", " + particleType["name"].asString() + ", " + 
                            centralityRangeName + ", " + zDCRangeName).c_str());
                  
                  grVSigmas.back().Clone()->
                     Write(("Data sigmas: " + variable["name"].asString() + 
                            ", " + particleType["name"].asString() + ", " + 
                            centralityRangeName + ", " + zDCRangeName).c_str());
                  
                  fVMeans.back().Clone()->
                     Write(("Fit means: " + variable["name"].asString() + 
                            ", " + particleType["name"].asString() + ", " + 
                            centralityRangeName + ", " + zDCRangeName).c_str());
                  
                  fVSigmas.back().Clone()->
                     Write(("Fit sigmas: " + variable["name"].asString() + 
                            ", " + particleType["name"].asString() + ", " + 
                            centralityRangeName + ", " + zDCRangeName).c_str());
               }

               double meanYMin = 1e31, meanYMax = -1e31;
               double sigmaYMin = 1e31, sigmaYMax = -1e31;
               
               for (unsigned long i = 0; i < grVMeans.size(); i++)
               {
                  meanYMin = 
                     Minimum(meanYMin, TMath::MinElement(grVMeans[i].GetN(), grVMeans[i].GetY()));
                  meanYMax = 
                     Maximum(meanYMax, TMath::MaxElement(grVMeans[i].GetN(), grVMeans[i].GetY()));
                  sigmaYMin = 
                     Minimum(sigmaYMin, TMath::MinElement(grVSigmas[i].GetN(), grVSigmas[i].GetY()));
                  sigmaYMax = 
                     Maximum(sigmaYMax, TMath::MaxElement(grVSigmas[i].GetN(), grVSigmas[i].GetY()));

                  const Json::Value zDCBin = calibrationInput["zdc_bins"][static_cast<int>(i)];
                  
                  grVMeans[i].SetMarkerStyle(zDCBin["marker_style"].asInt());
                  grVMeans[i].SetMarkerSize(1.4);
                  grVMeans[i].SetMarkerColorAlpha(zDCBin["color"].asInt(), 0.8);
                  grVMeans[i].SetLineColorAlpha(zDCBin["color"].asInt(), 0.8);
                  fVMeans[i].SetLineColorAlpha(zDCBin["color"].asInt(), 0.9);
                  fVMeans[i].SetLineStyle(3);
                  
                  grVSigmas[i].SetMarkerStyle(zDCBin["marker_style"].asInt());
                  grVSigmas[i].SetMarkerSize(1.4);
                  grVSigmas[i].SetMarkerColorAlpha(zDCBin["color"].asInt(), 0.8);
                  grVSigmas[i].SetLineColorAlpha(zDCBin["color"].asInt(), 0.8);
                  fVSigmas[i].SetLineColorAlpha(zDCBin["color"].asInt(), 0.9);
                  fVSigmas[i].SetLineStyle(2);
               }
               
               TCanvas canv("", "", 800, 800);
               
               TLegend legend{0.15, 0.7, 0.88, 0.89};
               legend.SetNColumns(3);
               legend.SetLineColorAlpha(0, 0.);
               legend.SetFillColorAlpha(0, 0.);
               

               gPad->SetLeftMargin(0.135);
               
               TH1 *meansFrame = 
                  gPad->DrawFrame(pTMin - 0.1, meanYMin - (meanYMax - meanYMin)*0.05, 
                                  pTMax*1.05, meanYMax + (meanYMax - meanYMin)*0.35);
               
               meansFrame->GetXaxis()->SetTitle("p_{T} [GeV/c]");
               meansFrame->GetYaxis()->
                  SetTitle(("#mu_{" + variable["tex_name"].asString() + "}").c_str());
               meansFrame->GetXaxis()->SetTitleOffset(1.1);
               meansFrame->GetYaxis()->SetTitleOffset(2.0);
               meansFrame->Draw("SAME AXIS X+ Y+");
               
               for (unsigned long i = 0; i < calibrationInput["zdc_bins"].size(); i++)
               {
                  const std::string zDCRangeName = 
                     calibrationInput["zdc_bins"][static_cast<int>(i)]["min"].asString() + 
                     "<z_{DC}<" + 
                     calibrationInput["zdc_bins"][static_cast<int>(i)]["max"].asString();
                  
                  legend.AddEntry(&grVMeans[i], zDCRangeName.c_str(), "P");
                  grVMeans[i].Clone()->Draw("SAME P");
                  fVMeans[i].Clone()->Draw("SAME");
               }

               legend.DrawClone();
               PrintCanvas(&canv, outputDir + detector["name"].asString() + "_means_" + 
                           variable["name"].asString() + "_" + 
                           particleType["name_short"].asString() +
                           centralityRangePathName);

               canv.Clone()->Write((variable["name"].asString() + " means: " + 
                                    particleType["name"].asString() + ", " + 
                                    centralityRangeName).c_str());
               
               legend.Clear();
               canv.Clear();

               TH1 *sigmasFrame = 
                  gPad->DrawFrame(pTMin - 0.1, sigmaYMin/1.1, 
                                  pTMax*1.05, sigmaYMax*1.4);
               sigmasFrame->GetXaxis()->SetTitle("p_{T} [GeV/c]");
               sigmasFrame->GetYaxis()->
                  SetTitle(("#sigma_{" + variable["name"].asString() + "}").c_str());
               sigmasFrame->Draw("SAME AXIS X+ Y+");
               
               for (unsigned long i = 0; i < calibrationInput["zdc_bins"].size(); i++)
               { 
                  const std::string zDCRangeName = 
                     calibrationInput["zdc_bins"][static_cast<int>(i)]["min"].asString() + 
                     "<z_{DC}<" + 
                     calibrationInput["zdc_bins"][static_cast<int>(i)]["max"].asString();
                  
                  legend.AddEntry(&grVSigmas[i], zDCRangeName.c_str(), "P");
                  grVSigmas[i].Clone()->Draw("SAME P");
                  fVSigmas[i].Clone()->Draw("SAME");
               }

               legend.DrawClone();
               PrintCanvas(&canv, outputDir + detector["name"].asString() + "_sigmas_" + 
                           variable["name"].asString() + "_" + 
                           particleType["name_short"].asString() +
                           centralityRangePathName);

               canv.Clone()->Write((variable["name"].asString() + " sigmas: " + 
                                    particleType["name"].asString() + ", " + 
                                    centralityRangeName).c_str());

               TCanvas parCanv("", "", 800, 800);
               parCanv.Divide(2, 2);
               
               parCanv.cd(1);
               gPad->SetRightMargin(0.13);
               meansDistr.GetXaxis()->SetTitle(variable["tex_name"].asCString());
               meansDistr.GetYaxis()->SetTitle("p_{T}");
               meansDistr.Draw("COLZ");

               parCanv.cd(2);
               gPad->SetRightMargin(0.13);
               sigmasDistr.GetXaxis()->SetTitle(variable["tex_name"].asCString());
               sigmasDistr.GetYaxis()->SetTitle("p_{T}");
               sigmasDistr.Draw("COLZ");

               parCanv.cd(3);
               gPad->SetLogz();
               gPad->SetRightMargin(0.13);
               meansDiffDistr.GetXaxis()->SetTitle(variable["tex_name"].asCString());
               meansDiffDistr.GetYaxis()->SetTitle("p_{T}");
               meansDiffDistr.Draw("COLZ");

               parCanv.cd(4);
               gPad->SetLogz();
               gPad->SetRightMargin(0.13);
               sigmasDiffDistr.GetXaxis()->SetTitle(variable["tex_name"].asCString());
               sigmasDiffDistr.GetYaxis()->SetTitle("p_{T}");
               sigmasDiffDistr.Draw("COLZ");
               
               PrintCanvas(&parCanv, outputDir + detectorName + 
                           "/fitPar_" + variable["name"].asString() + "_" + 
                           particleType["name_short"].asString() + 
                           centralityRangePathName);
            }
         }
      }
      outputFile.Close();
   }
   pBar.Print(1.);
}

void PerformFits(TH3F *hist, TGraphErrors &grMeans, TGraphErrors &grSigmas,
                 const Json::Value& calibrationInput,
                 const Json::Value& detector, const Json::Value& variable,
                 const Json::Value& zDCBin, const Json::Value& particleType,
                 const Json::Value& centralityBin)
{ 
   const double minBinX = hist->GetXaxis()->GetBinLowEdge(1);
   const double maxBinX = hist->GetXaxis()->GetBinUpEdge(hist->GetXaxis()->GetNbins());
   const double binWidth = hist->GetXaxis()->GetBinWidth(1);

   const std::string centralityRangeName = centralityBin["min"].asString() + "-" + 
                                           centralityBin["max"].asString() + "%";
   const std::string centralityRangePathName = "_c" + centralityBin["min"].asString() + 
                                                 "-" + centralityBin["max"].asString();

   const std::string zDCRangeName = zDCBin["min"].asString() + "<zDC<" + zDCBin["max"].asString();
   const std::string zDCRangePathName = "_zDC" + zDCBin["min"].asString() + 
                                        "-" + zDCBin["max"].asString();

   TCanvas canv("dval vs pT", "", 
                calibrationInput["pt_nbinsx"].asDouble()*400.,
                calibrationInput["pt_nbinsy"].asDouble()*400.);
   
   canv.Divide(calibrationInput["pt_nbinsx"].asInt(), calibrationInput["pt_nbinsy"].asInt());
 
   auto PerformFitsInRange = [&](const double pTMin, const double pTMax, 
                                 TF1 *meansFit, TF1 *sigmasFit,
                                 TF1 *limitMeansFitFunc = NULL, TF1 *limitSigmasFitFunc = NULL,
                                 bool limitParByFit = false)
   {
      // clearing the graphs after previous fits
      for (int i = grMeans.GetN() - 1; i >= 0; i--)
      {
         grMeans.RemovePoint(i);
         grSigmas.RemovePoint(i);
      }

      int iCanv = 1;
      
      for (const Json::Value& pTBin : calibrationInput["pt_bins"])
      {
         const double pT = Average(pTBin["min"].asDouble(), pTBin["max"].asDouble());
         if (pT < pTMin || pT > pTMax) continue;
         
         canv.cd(iCanv);
         iCanv++;
      
         TH1D *distrVariableProj = hist->
            ProjectionX(((std::string) hist->GetName() + "_projX").c_str(), 
                        hist->GetYaxis()->FindBin(pTBin["min"].asDouble() + 1e-6), 
                        hist->GetYaxis()->FindBin(pTBin["max"].asDouble() - 1e-6),
                        hist->GetZaxis()->FindBin(centralityBin["min"].asDouble() + 1e-6),
                        hist->GetZaxis()->FindBin(centralityBin["max"].asDouble() - 1e-6));

         const std::string pTRangeName = DtoStr(pTBin["min"].asDouble(), 1) + "<pT<" + 
                                         DtoStr(pTBin["max"].asDouble(), 1);
         
         if (distrVariableProj->Integral(1, distrVariableProj->GetXaxis()->GetNbins()) < 
             Par.minIntegralValue) 
         {
            PrintInfo("Integral is insufficient for projection of " + variable["name"].asString() + 
                      ", " + detector["name"].asString() + ", " + particleType["name"].asString() + 
                      " at " + zDCRangeName + ", " + centralityRangeName + ", " + pTRangeName);
            continue;
         }
         
         double minX = 0., maxX = -1.;
         for (int k = 1; k <= distrVariableProj->GetXaxis()->GetNbins(); k++)
         {
            if (distrVariableProj->GetBinContent(k) > 1e-7)
            {
               minX = distrVariableProj->GetXaxis()->GetBinLowEdge(k);
               break;
            }
         }
         
         for (int k = distrVariableProj->GetXaxis()->GetNbins(); k > 1; k--)
         {
            if (distrVariableProj->GetBinContent(k) > 1e-7)
            {
               maxX = distrVariableProj->GetXaxis()->GetBinUpEdge(k);
               break;
            }
         }
         
         if (minX > maxX) 
         {
            PrintWarning("Something wrong for projection of " + variable["name"].asString() + ", " + 
                         detector["name"].asString() + ", " + particleType["name"].asString() + 
                         " at " + zDCRangeName + ", " + centralityRangeName + ", " + pTRangeName);
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

         distrVariableProj->GetXaxis()->SetTitle(variable["tex_name"].asCString());
         distrVariableProj->SetTitle("");
         distrVariableProj->SetTitleSize(0.06, "X");
         distrVariableProj->SetTitleSize(0.06, "Y");
         distrVariableProj->SetLabelSize(0.06, "X");
         distrVariableProj->SetLabelSize(0.06, "Y");
    
         distrVariableProj->GetXaxis()->SetRange(distrVariableProj->GetXaxis()->FindBin(minX+0.01),
                                                 distrVariableProj->GetXaxis()->FindBin(maxX-0.01));
         
         const double maxBinVal = 
            distrVariableProj->GetBinContent(distrVariableProj->GetMaximumBin());
         
         // scale limits
         fitFuncGaus.SetParLimits(0, maxBinVal/5., maxBinVal);
         fitFuncDVal.SetParLimits(0, maxBinVal/5., maxBinVal);
         fitFuncDVal.SetParLimits(3, maxBinVal/20., maxBinVal);
         
         fitFuncGaus.SetRange(minBinX/5., maxBinX/5.);
         
         distrVariableProj->Fit(&fitFuncGaus, "RQMBN");

         fitFuncGaus.SetRange(minBinX, maxBinX);
         fitFuncBG.SetRange(minBinX, maxBinX);
         
         for (int k = 0; k < 3; k++)
         {
            fitFuncDVal.SetParameter(k, fitFuncGaus.GetParameter(k));
         }
         
         fitFuncDVal.SetRange(minBinX, maxBinX);
         distrVariableProj->Fit(&fitFuncDVal, "RQMBN");

         for (unsigned short j = 1; j <= Par.fitNTries; j++)
         {
            for (int k = 0; k < fitFuncDVal.GetNpar(); k++)
            {
               fitFuncDVal.SetParLimits(k, fitFuncDVal.GetParameter(k)/
                                        (1. + 1./static_cast<double>(j)),
                                        fitFuncDVal.GetParameter(k)*
                                        (1. + 1./static_cast<double>(j)));
            }
            
            distrVariableProj->Fit(&fitFuncDVal, "RQMBN");
         }

         for (int j = 0; j < 3; j++)
         {
            fitFuncGaus.SetParameter(j, fitFuncDVal.GetParameter(j));
            fitFuncBG.SetParameter(j, fitFuncDVal.GetParameter(j + 3));
         }
         
         distrVariableProj->SetMarkerStyle(20);
         distrVariableProj->SetMarkerSize(0.7);
         distrVariableProj->SetMarkerColorAlpha(kBlack, 0.8);
         distrVariableProj->SetLineColorAlpha(kBlack, 0.8);
         distrVariableProj->SetMaximum(maxBinVal*1.2);
         
         gPad->SetLeftMargin(0.155);
         gPad->SetBottomMargin(0.115);
         
         distrVariableProj->Clone()->Draw("P");
         fitFuncDVal.Clone()->Draw("SAME");
         fitFuncBG.Clone()->Draw("SAME");
         fitFuncGaus.Clone()->Draw("SAME");
         
         Par.texText.DrawLatexNDC(0.17, 0.85, pTRangeName.c_str());
         Par.texText.DrawLatexNDC(0.17, 0.79, zDCRangeName.c_str());
         Par.texText.DrawLatexNDC(0.17, 0.73, particleType["name"].asCString());
         Par.texText.DrawLatexNDC(0.17, 0.66, centralityRangeName.c_str());
         
         if (fabs(fitFuncDVal.GetParameter(1)) < variable["abs_max_fit_val"].asDouble() && 
             fabs(fitFuncDVal.GetParameter(2)) < variable["abs_max_fit_val"].asDouble())
         {
            grMeans.AddPoint(pT, fitFuncDVal.GetParameter(1));
            grSigmas.AddPoint(pT, fitFuncDVal.GetParameter(2));
            //grMeans.SetPointError(grMeans.GetN() - 1, 0, fitFuncDVal.GetParError(1));
            //grSigmas.SetPointError(grSigmas.GetN() - 1, 0, fitFuncDVal.GetParError(2));
         }
      }
      if (grMeans.GetN() == 0) 
      {
         PrintError("Graph is empty for " + variable["name"].asString() + ", " + 
                    detector["name"].asString() + ", " + particleType["name"].asString() + 
                    " at " + zDCRangeName + ", " + centralityRangeName);
      }
      
      meansFit->SetRange(pTMin, pTMax);
      sigmasFit->SetRange(pTMin, pTMax);
      
      grMeans.Fit(meansFit, "RQMNW");
      grSigmas.Fit(sigmasFit, "RQMNW");
   };
   
   TF1 meansFit("meansFit", Par.meansFitPrelimFunc.c_str());
   TF1 sigmasFit("sigmasFit", Par.sigmasFitPrelimFunc.c_str());
   
   for (const Json::Value& pT : calibrationInput["consequetive_pt_fit_ranges"])
   {
      PerformFitsInRange(pT["min"].asDouble(), pT["max"].asDouble(), &meansFit, &sigmasFit);
   }

   TF1 *oldMeansFit = (TF1 *) meansFit.Clone();
   TF1 *oldSigmasFit = (TF1 *) sigmasFit.Clone();
   
   meansFit.SetCurrent((TF1 *) TF1("meansFit", detector["means_fit"].asCString()).Clone());
   sigmasFit.SetCurrent((TF1 *) TF1("sigmasFit", detector["sigmas_fit"].asCString()).Clone());

   PerformFitsInRange(calibrationInput["pt_bins"].front()["min"].asDouble()/1.05, 
                      calibrationInput["pt_bins"].back()["max"].asDouble()*1.05, 
                      &meansFit, &sigmasFit, oldMeansFit, oldSigmasFit, true);
   
   const std::string outputFileNameNoExt = 
      "output/ResidualCal/" + Par.runName + "/" + detector["name"].asString() + "/" + 
      variable["name"].asString() + "_" + particleType["name_short"].asString() + 
      centralityRangePathName + zDCRangePathName;
   PrintCanvas(&canv, outputFileNameNoExt, false, true);
}

#endif /* CALIBRATE_SIGMALIZED_RESIDUALS_CPP */
