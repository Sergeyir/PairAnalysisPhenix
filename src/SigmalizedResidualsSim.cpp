/** 
 *  @file   SigmalizedResidualsSimCailbration.cpp 
 *  @brief  Contains realisation of functions that are used for estimation of values for calibration of sigmalized residuals dphi and dz from the PHENIX simulation
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysisPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef SIGMALIZED_RESIDUALS_SIM_CPP
#define SIGMALIZED_RESIDUALS_SIM_CPP

#include "../include/SigmalizedResidualsSim.hpp"

// this namespace is only used so that documentation will not become a mess
// so there is no need to enforce the contents inside of it 
// being accessed only via the scope resolution operator in this file
using namespace SigmalizedResidualsSim;

int main(int argc, char **argv)
{
   using namespace SigmalizedResidualsSim;

   TH1::AddDirectory(false);
   TH2::AddDirectory(false);
   TH3::AddDirectory(false);

   if (argc < 2 || argc > 3) 
   {
      std::string errMsg = "Expected 1-2 parameters while " + std::to_string(argc - 1) + 
                           " parameter(s) were provided \n";
      errMsg += "Usage: bin/SigmalizedResidualsSim inputFile numberOfThreads=" + 
                std::to_string(std::thread::hardware_concurrency()) + "*\n";
      errMsg += "*: default argument is the number of threads on the current machine \n";
      CppTools::PrintError(errMsg);
   }

   // initializing this program parameters
   inputYAMLMain.OpenFile(argv[1], "main");
   inputYAMLMain.CheckStatus("main");

   runName = inputYAMLMain["run_name"].as<std::string>();

   unsigned int numberOfThreads;
   if (argc > 2) numberOfThreads = std::stoi(argv[2]);
   else numberOfThreads = std::thread::hardware_concurrency();
   if (numberOfThreads == 0) CppTools::PrintError("Number of threads must be bigger than 0");

   // initializing ROOT parameters
   ROOT::EnableThreadSafety();
   gErrorIgnoreLevel = kWarning;
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);

   TDirectory::AddDirectory(kFALSE);

   ROOT::EnableImplicitMT(numberOfThreads);

   outputDir = "output/SigmalizedResidualsSim/" + runName + "/";
   system(("mkdir -p " + outputDir).c_str());

   const std::string inputFileName = "data/PostSim/" + runName + "/SingleTrack/all.root";
   CppTools::CheckInputFile(inputFileName);

   inputFile = TFile::Open(inputFileName.c_str(), "READ");
   outputFile = TFile::Open((outputDir + "all_fits.root").c_str(), "RECREATE");

   const std::string detectorsConfiguration = 
      inputYAMLMain["detectors_configuration"].as<std::string>();
   
   if (detectorsConfiguration.size() != 7)
   {
      CppTools::PrintError("Detector configuration size is " + 
                           std::to_string(detectorsConfiguration.size()) + 
                           " while 7 has been expected");
   }

   if (detectorsConfiguration[2] == '1') numberOfIterations += 1;
   if (detectorsConfiguration[3] == '1') numberOfIterations += 2;
   if (detectorsConfiguration[4] == '1') numberOfIterations += 1;
   if (detectorsConfiguration[5] == '1') numberOfIterations += 1;
   if (detectorsConfiguration[6] == '1') numberOfIterations += 8;

   numberOfIterations *= 4;

   system(("mkdir -p data/Parameters/SigmalizedResidualsSim/" + runName).c_str());

   PerformCalibrationsForDetector("PC2", (detectorsConfiguration[2] == '1'));
   PerformCalibrationsForDetector("PC3e", (detectorsConfiguration[3] == '1'));
   PerformCalibrationsForDetector("PC3w", (detectorsConfiguration[3] == '1'));
   PerformCalibrationsForDetector("TOFe", (detectorsConfiguration[4] == '1'));
   PerformCalibrationsForDetector("TOFw", (detectorsConfiguration[5] == '1'));
   for (int i = 0; i < 4; i++)
   {
      PerformCalibrationsForDetector("EMCale" + std::to_string(i), detectorsConfiguration[6] == '1');
      PerformCalibrationsForDetector("EMCalw" + std::to_string(i), detectorsConfiguration[6] == '1');
   }

   pBar.Finish();

   inputFile->Close();
   outputFile->Close();

   return 0;
}

void SigmalizedResidualsSim::PerformCalibrationsForDetector(const std::string& detectorName,
                                                            const bool performCalibration)
{
   parametersOutput.open("data/Parameters/SigmalizedResidualsSimSim/" + 
                         runName + "/" + detectorName + ".txt");

   if (!performCalibration)
   {
      parametersOutput << 0;
      parametersOutput.close();
      return;
   }
   
   // writing the calibration status and the number of parameters of approximation
   parametersOutput << 1 << " " << 4 << " " << 3 << std::endl;

   for (unsigned long i = 0; i < variableNames.size(); i++)
   {
      for (unsigned long j = 0; j < particleCharges.size(); j++)
      {
         PerformCalibrationsVsPT(detectorName, variableNames[i], particleCharges[j]);

         if (i < variableNames.size() - 1 || j < particleCharges.size() - 1)
         {
            parametersOutput << std::endl;
         }
      }
   }

   parametersOutput.close();
}

void SigmalizedResidualsSim::PerformCalibrationsVsPT(const std::string& detectorName, 
                                                     const std::string& variableName,
                                                     const int charge)
{
   pBar.Print(static_cast<double>(numberOfCalls)/static_cast<double>(numberOfIterations));

   const std::string chargeName = ((charge > 0) ? "charge>0" : "charge<0");
   const std::string chargeNameShort = ((charge > 0) ? "pos" : "neg");

   const std::string distrDValVsPTName = variableName + "/" + variableName + " vs pT: " + 
                                         detectorName + ", " + chargeName;

   TH2F *distrDValVsPT = static_cast<TH2F *>(inputFile->Get(distrDValVsPTName.c_str()));

   if (!distrDValVsPT)
   {
      CppTools::PrintError("No histogram named \"" + distrDValVsPTName + 
                           "\" in file data/PostSim/" + runName + "/SingleTrack/all.root");
   }

   TGraphErrors grMeansDVal;
   TGraphErrors grSigmasDVal;

   const double xMin = 
      distrDValVsPT->GetXaxis()->GetBinLowEdge(1);
   const double xMax = 
      distrDValVsPT->GetXaxis()->GetBinUpEdge(distrDValVsPT->GetXaxis()->GetNbins());

   double pTMin = 1e31;

   for (int i = 1; i < distrDValVsPT->GetYaxis()->FindBin(3.); i++) // pT < 3 GeV/c
   {
      if (distrDValVsPT->Integral(1, distrDValVsPT->GetXaxis()->GetNbins(), i, i) < 1e-15) continue;

      pTMin = CppTools::Minimum(pTMin, distrDValVsPT->GetYaxis()->GetBinCenter(i));

      TH1D *distrDValProj = distrDValVsPT->
         ProjectionX((variableName + ": " + detectorName + ", " + chargeName + ", pT" +
                      CppTools::DtoStr(distrDValVsPT->GetYaxis()->GetBinLowEdge(i), 1) + "-" + 
                      CppTools::DtoStr(distrDValVsPT->GetYaxis()->
                             GetBinLowEdge(distrDValVsPT->GetYaxis()->GetNbins()), 
                      1)).c_str(), i, i);

      TF1 fitDistrDVal((detectorName + variableName + chargeName + std::to_string(i)).c_str(), 
                       "gaus(0) + gaus(3)");

      // set of alternative fit functions used for uncertainty estimation 
      // by varying ranges of approximation around mean by n*sigma of the main fit
      // for first vector approximation ranges are varied symmetrically within mean
      // for second vector approximation ranges are varied within right of mean; left range is 1sigma
      // for third vector approximation ranges are varied within left of mean; right range is 1sigma
      std::vector<TF1> fitDistrDValAlt, fitDistrDValAltRight, fitDistrDValAltLeft;

      const double maxBinVal = distrDValProj->GetBinContent(distrDValProj->GetMaximumBin());

      for (unsigned long i = 0; i < 4; i++)
      {
         fitDistrDValAlt.emplace_back(("fitDistrDValAlt_" + std::to_string(i) + 
                                      "_" + std::to_string(i)).c_str(), "gaus(0) + gaus(3)");
         fitDistrDValAlt.back().SetParLimits(0, maxBinVal/2., maxBinVal);
         fitDistrDValAlt.back().SetParLimits(3, maxBinVal/20., maxBinVal);
         fitDistrDValAltRight.emplace_back(("fitDistrDValAltRight_" + std::to_string(i) + 
                                         "_" + std::to_string(i)).c_str(), "gaus(0) + gaus(3)");
         fitDistrDValAltRight.back().SetParLimits(0, maxBinVal/2., maxBinVal);
         fitDistrDValAltRight.back().SetParLimits(3, maxBinVal/20., maxBinVal);
         fitDistrDValAltLeft.emplace_back(("fitDistrDValAltLeft_" + std::to_string(i) + 
                                         "_" + std::to_string(i)).c_str(), "gaus(0) + gaus(3)");
         fitDistrDValAltLeft.back().SetParLimits(0, maxBinVal/2., maxBinVal);
         fitDistrDValAltLeft.back().SetParLimits(3, maxBinVal/20., maxBinVal);
      }

      fitDistrDVal.SetParameters(maxBinVal, 0., distrDValProj->GetXaxis()->GetBinWidth(1)*10.);

      fitDistrDVal.SetParLimits(0, maxBinVal/3., maxBinVal);
      fitDistrDVal.SetParLimits(1, xMin/5., xMax/5.);
      fitDistrDVal.SetParLimits(2, distrDValProj->GetXaxis()->GetBinWidth(1)*10., xMax/2.);
      fitDistrDVal.SetParLimits(3, maxBinVal/20., maxBinVal/3.);
      fitDistrDVal.SetParLimits(4, xMin/5., xMax/5.);
      fitDistrDVal.SetParLimits(5, xMax/5., xMax*5.);

      fitDistrDVal.SetRange(xMin, xMax);

      for (unsigned int i = 1; i <= fitNTries; i++)
      {
         distrDValProj->Fit(&fitDistrDVal, "RQMBN");

         fitDistrDVal.SetParLimits(0, fitDistrDVal.GetParameter(0)/
                                   (1. + 2./static_cast<double>(i*i*i)),
                                   fitDistrDVal.GetParameter(0)*
                                   (1. + 2./static_cast<double>(i*i*i)));
         fitDistrDVal.SetParLimits(1, fitDistrDVal.GetParameter(1) - 
                                   fitDistrDVal.GetParameter(2)*5./static_cast<double>(i*i*i),
                                   fitDistrDVal.GetParameter(1) + 
                                   fitDistrDVal.GetParameter(2)*5./static_cast<double>(i*i*i));
         fitDistrDVal.SetParLimits(2, fitDistrDVal.GetParameter(2)/
                                   (1. + 5./static_cast<double>(i*i*i)),
                                   fitDistrDVal.GetParameter(2)*
                                   (1. + 5./static_cast<double>(i*i*i)));
         fitDistrDVal.SetParLimits(3, fitDistrDVal.GetParameter(3)/
                                   (1. + 2./static_cast<double>(i*i*i)),
                                   fitDistrDVal.GetParameter(3)*
                                   (1. + 2./static_cast<double>(i*i*i)));
         fitDistrDVal.SetParLimits(4, fitDistrDVal.GetParameter(4) - 
                                   fitDistrDVal.GetParameter(5)/static_cast<double>(i*i*i),
                                   fitDistrDVal.GetParameter(4) + 
                                   fitDistrDVal.GetParameter(5)/static_cast<double>(i*i*i));
         fitDistrDVal.SetParLimits(5, fitDistrDVal.GetParameter(5)/
                                   (1. + 5./static_cast<double>(i*i*i)),
                                   fitDistrDVal.GetParameter(5)*
                                   (1. + 5./static_cast<double>(i*i*i)));

         fitDistrDVal.SetRange(fitDistrDVal.GetParameter(1) - 5.*fitDistrDVal.GetParameter(2), 
                               fitDistrDVal.GetParameter(1) + 5.*fitDistrDVal.GetParameter(2));
      }

      distrDValProj->Fit(&fitDistrDVal, "RQMB");

      for (unsigned long i = 0; i < fitDistrDValAlt.size(); i++)
      {
         fitDistrDValAlt[i].
            SetRange(fitDistrDVal.GetParameter(1) - fitDistrDVal.GetParameter(2)*
                     static_cast<double>(i + 1)*2., fitDistrDVal.GetParameter(1) + 
                     fitDistrDVal.GetParameter(2)*static_cast<double>(i + 1)*2.);
         fitDistrDValAltRight[i].
            SetRange(fitDistrDVal.GetParameter(1) - fitDistrDVal.GetParameter(2), 
                     fitDistrDVal.GetParameter(1) + 
                     fitDistrDVal.GetParameter(2)*static_cast<double>(i + 1)*2.);
         fitDistrDValAltLeft[i].
            SetRange(fitDistrDVal.GetParameter(1) - fitDistrDVal.GetParameter(2)*
                     static_cast<double>(i + 1)*2., fitDistrDVal.GetParameter(1) + 
                     fitDistrDVal.GetParameter(2));

         for (int j = 0; j < fitDistrDVal.GetNpar(); j++)
         {
            fitDistrDValAlt[i].SetParameter(j, fitDistrDVal.GetParameter(j)); 
            fitDistrDValAltRight[i].SetParameter(j, fitDistrDVal.GetParameter(j)); 
            fitDistrDValAltLeft[i].SetParameter(j, fitDistrDVal.GetParameter(j)); 

            if (j == 0 || j == 3)
            {
               fitDistrDValAlt[i].SetParLimits(j, fitDistrDVal.GetParameter(j)/1.2, 
                                              fitDistrDVal.GetParameter(j)*1.2); 
               fitDistrDValAltRight[i].SetParLimits(j, fitDistrDVal.GetParameter(j)/1.2, 
                                                   fitDistrDVal.GetParameter(j)*1.2); 
               fitDistrDValAltLeft[i].SetParLimits(j, fitDistrDVal.GetParameter(j)/1.2, 
                                                  fitDistrDVal.GetParameter(j)*1.2); 
            }
            else if (j == 2 || j == 4)
            {
               fitDistrDValAlt[i].SetParLimits(j, fitDistrDVal.GetParameter(j)/1.5, 
                                              fitDistrDVal.GetParameter(j)*1.5); 
               fitDistrDValAltRight[i].SetParLimits(j, fitDistrDVal.GetParameter(j)/1.5, 
                                                   fitDistrDVal.GetParameter(j)*1.5); 
               fitDistrDValAltLeft[i].SetParLimits(j, fitDistrDVal.GetParameter(j)/1.5, 
                                                   fitDistrDVal.GetParameter(j)*1.5); 
            }
         }
         distrDValProj->Fit(&fitDistrDValAlt[i], "RQMBNL");
         distrDValProj->Fit(&fitDistrDValAltRight[i], "RQMBNL");
         distrDValProj->Fit(&fitDistrDValAltLeft[i], "RQMBNL");
      }

      distrDValProj->Write();
      
      grMeansDVal.AddPoint(distrDValVsPT->GetYaxis()->GetBinCenter(i), 
                           fitDistrDVal.GetParameter(1));
      grSigmasDVal.AddPoint(distrDValVsPT->GetYaxis()->GetBinCenter(i), 
                            fitDistrDVal.GetParameter(2));

      const double meanError = 
         CppTools::StandardError(fitDistrDValAlt[0].GetParameter(1),
                                 fitDistrDValAlt[1].GetParameter(1),
                                 fitDistrDValAlt[2].GetParameter(1),
                                 fitDistrDValAlt[3].GetParameter(1), 
                                 fitDistrDValAltRight[0].GetParameter(1),
                                 fitDistrDValAltRight[1].GetParameter(1),
                                 fitDistrDValAltRight[2].GetParameter(1),
                                 fitDistrDValAltRight[3].GetParameter(1),
                                 fitDistrDValAltLeft[0].GetParameter(1),
                                 fitDistrDValAltLeft[1].GetParameter(1),
                                 fitDistrDValAltLeft[2].GetParameter(1),
                                 fitDistrDValAltLeft[3].GetParameter(1),
                                 fitDistrDVal.GetParameter(1));

      const double sigmaError = 
         CppTools::StandardError(fitDistrDValAlt[0].GetParameter(2),
                                 fitDistrDValAlt[1].GetParameter(2),
                                 fitDistrDValAlt[2].GetParameter(2),
                                 fitDistrDValAlt[3].GetParameter(2), 
                                 fitDistrDValAltRight[0].GetParameter(2),
                                 fitDistrDValAltRight[1].GetParameter(2),
                                 fitDistrDValAltRight[2].GetParameter(2),
                                 fitDistrDValAltRight[3].GetParameter(2),
                                 fitDistrDValAltLeft[0].GetParameter(2),
                                 fitDistrDValAltLeft[1].GetParameter(2),
                                 fitDistrDValAltLeft[2].GetParameter(2),
                                 fitDistrDValAltLeft[3].GetParameter(2),
                                 fitDistrDVal.GetParameter(2));

      grMeansDVal.SetPointError(grMeansDVal.GetN() - 1, 0., meanError);
      grSigmasDVal.SetPointError(grSigmasDVal.GetN() - 1, 0., sigmaError);
   }

   if (grMeansDVal.GetN() == 0) 
   {
      CppTools::PrintError("Number of points in a graph for " + detectorName + 
                           ", " + variableName + ", " + chargeName + " equals 0");
   }

   TF1 fitMeans("means fit", "[0] + [1]/x + [2]/x^2 + [3]*x", pTMin - 0.1, 3.1);
   TF1 fitSigmas("means fit", "[0] + [1]/x + [2]*x", pTMin - 0.1, 3.1);

   if (grMeansDVal.GetN() == 1) 
   {
      fitMeans.SetParameter(0, grMeansDVal.GetY()[0]);
      fitSigmas.SetParameter(0, grSigmasDVal.GetY()[0]);
   }
   else
   {
      grMeansDVal.Fit(&fitMeans, "RQMBN");
      grSigmasDVal.Fit(&fitSigmas, "RQMBN");
   }

   TCanvas meansVsPTCanv("means vs pT canvas", "", 800, 800);

   TH1F fitParVsPTFrame("means and sigmas frame", "", 10, pTMin, 3.1);

   fitParVsPTFrame.SetMinimum(CppTools::Minimum(TMath::MinElement(grMeansDVal.GetN(), 
                                                                  grMeansDVal.GetY())));
   if (fitParVsPTFrame.GetMinimum() < 0.) 
   {
      fitParVsPTFrame.SetMinimum(fitParVsPTFrame.GetMinimum()*1.2);
   }
   else
   {
      fitParVsPTFrame.SetMinimum(fitParVsPTFrame.GetMinimum()/1.2);
   }

   fitParVsPTFrame.SetMaximum(CppTools::Maximum(TMath::MaxElement(grMeansDVal.GetN(), 
                                                                  grMeansDVal.GetY())));
   if (fitParVsPTFrame.GetMaximum() < 0.) 
   {
      fitParVsPTFrame.SetMaximum(fitParVsPTFrame.GetMaximum()/1.2);
   }
   else
   {
      fitParVsPTFrame.SetMaximum(fitParVsPTFrame.GetMaximum()*1.2);
   }

   fitParVsPTFrame.Draw("AXIS");
   fitParVsPTFrame.Draw("SAME AXIS X+ Y+");

   grMeansDVal.SetMarkerStyle(20);
   grMeansDVal.SetLineColor(kBlack);
   grMeansDVal.SetLineWidth(2);
   grMeansDVal.SetMarkerColor(kBlack);
   grMeansDVal.SetMarkerSize(1.);

   fitMeans.SetLineColor(kRed-3);

   fitMeans.Clone()->Draw("SAME");
   grMeansDVal.Clone()->Draw("P");

   for (int i = 0; i < fitMeans.GetNpar(); i++)
   {
      parametersOutput << fitMeans.GetParameter(i);
      parametersOutput << " ";
   }

   ROOTTools::PrintCanvas(&meansVsPTCanv, outputDir + "means_" + variableName + "_" + 
                          detectorName + "_" + chargeNameShort);

   TCanvas sigmasVsPTCanv("means vs pT canvas", "", 800, 800);

   fitParVsPTFrame.SetMinimum(CppTools::Minimum(TMath::MinElement(grSigmasDVal.GetN(), 
                                                                  grSigmasDVal.GetY()))/1.5);
   fitParVsPTFrame.SetMaximum(CppTools::Maximum(TMath::MaxElement(grSigmasDVal.GetN(), 
                                                                  grSigmasDVal.GetY()))*1.5);

   grSigmasDVal.SetMarkerStyle(20);
   grSigmasDVal.SetLineColor(kBlack);
   grSigmasDVal.SetLineWidth(2);
   grSigmasDVal.SetMarkerColor(kBlack);
   grSigmasDVal.SetMarkerSize(1);

   fitParVsPTFrame.Draw("AXIS");
   fitParVsPTFrame.Draw("SAME AXIS X+ Y+");

   fitSigmas.SetLineColor(kRed-3);

   fitSigmas.Clone()->Draw("SAME");
   grSigmasDVal.Clone()->Draw("P");

   ROOTTools::PrintCanvas(&sigmasVsPTCanv, outputDir + "sigmas_" +variableName + "_" + 
                          detectorName + "_" + chargeNameShort);

   for (int i = 0; i < fitSigmas.GetNpar(); i++)
   {
      parametersOutput << fitSigmas.GetParameter(i);
      if (i < fitSigmas.GetNpar() - 1) parametersOutput << " ";
   }

   numberOfCalls++;
}

#endif /* SIGMALIZED_RESIDUALS_SIM_CPP */
