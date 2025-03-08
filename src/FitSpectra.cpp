#pragma once

#include "ErrorHandler.hpp"
#include "IOTools.hpp"
#include "MathTools.hpp"

#include "json/json.h"

#include "yaml-cpp/yaml.h"

const double MASS_PI = 0.139570;
const double MASS_K = 0.493677;
const double MASS_P = 0.938272;

//tsallis distribution
//https://iopscience.iop.org/article/10.1088/1742-6596/878/1/012016/pdf
const std::string TSALLIS_FUNC = 
   "[0]*x^([5]-1)*(1 + ([1] - 1)*(sqrt([4] + x^2)^([5]) - [2])/[3])^(-1./([1]-1.))";

const std::string EXP_FUNC = 
   "[0]*exp([1] + [2]*x*x + [3]*sqrt(x))";

struct
{
   const std::string collisionSystem = "HeAu200";
   const std::string particleName = "k-";
   const std::string experimentName = "PHENIX";
   const std::string inputFileName = "./data/Spectra/" + collisionSystem + 
                                     "/" + particleName + experimentName + ".yaml";
   const std::string fitFunc = TSALLIS_FUNC;
   //const std::string fitFunc = EXP_FUNC;
   
   const double mass = MASS_P;
   const std::string centralityClass = "0-88\\%";
   const double kappa = 2.;
} Par;

void FitSpectra()
{
   gROOT->SetBatch(kTRUE);
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   
   CheckInputFile(Par.inputFileName);
   TGraphErrors spectraGr = TGraphErrors();
   
   YAML::Node yamlFileContents = YAML::LoadFile(Par.inputFileName);

   int cIndex = -1;
   for (unsigned long i = 0; i < yamlFileContents["dependent_variables"].size(); i++)
   {
      Print(yamlFileContents["dependent_variables"][i]["qualifiers"][0]["value"].as<std::string>());
      if (yamlFileContents["dependent_variables"][i]["qualifiers"][0]["value"].as<std::string>() ==
          Par.centralityClass)
      {
         cIndex = i;
         break;
      }
   }
   if (cIndex == -1) PrintError("Centrality class " + Par.centralityClass + 
                                " was not found in file " + Par.inputFileName);

   auto pT = yamlFileContents["independent_variables"][0]["values"];
   auto y = yamlFileContents["dependent_variables"][cIndex]["values"];
   
   double pTMin = 1e+31;
   double pTMax = 1e-31;
   double yMin = 1e31;
   double yMax = 1e-31;
   
   for (int i = 0; i < yamlFileContents["independent_variables"][0]["values"].size(); i++)
   {
      double error = 0.;
      // summing all uncertainties
      for (auto errType : y[i]["errors"])
      {
         error += pow(errType["symerror"].as<double>(), 2);
      }
      error = sqrt(error);
      
      // some data files only have low and high edge pT bin values
      // and some only have the middle of the bin
      // therefore they need to be differentiated
      double pTValue = 0.;
      if (pT[i].size() == 1 || pT[i].size() == 3)
      {
         pTValue = pT[i]["value"].as<double>();
      }
      else if (pT[i].size() == 2)
      {
         pTValue = Average(pT[i]["low"].as<double>(), pT[i]["high"].as<double>());
      }
      else (PrintError("Unexpected number of values of pT columns in file " + Par.inputFileName));
      
      double yValue = y[i]["value"].as<double>();
      
      spectraGr.AddPoint(pTValue, yValue);
      spectraGr.SetPointError(i, 0., error);
      
      pTMin = Minimum(pTValue, pTMin);
      pTMax = Maximum(pTValue, pTMax);
      yMin = Minimum(yValue, yMin);
      yMax = Maximum(yValue, yMax);
   }
   
   TF1 fitFunc = TF1("fitFunc", Par.fitFunc.c_str());
   
   if (Par.fitFunc == TSALLIS_FUNC)
   {
      fitFunc.SetParameter(0, 10); //scale
      fitFunc.SetParameter(1, 1.01); //q
      fitFunc.SetParameter(2, 0.2); //mu
      fitFunc.SetParameter(3, 0.1); //T
      
      fitFunc.SetParLimits(1, 1.00001, 2.);
      fitFunc.SetParLimits(3, 0.01, 0.4);
      fitFunc.FixParameter(4, Par.mass*Par.mass); //m2
      fitFunc.FixParameter(5, Par.kappa); //kappa
   }
   else
   {
      fitFunc.SetParameters(1., 0., 0., 0., 0., 0., 0., 0., 0.);
   }

   spectraGr.Fit(&fitFunc, "QBN");


   const std::string particleNameTex = 
      yamlFileContents["dependent_variables"][cIndex]["qualifiers"][0]["name"].as<std::string>();
   
   TCanvas *canv = new TCanvas("canv", "canv", 600, 800);
   gPad->SetLeftMargin(0.12);
   TH1F *pad = gPad->DrawFrame(0.2, yMin/5., pTMax*1.3, yMax*5., 
                               (particleNameTex + " spectra").c_str());

   pad->Draw("SAME AXIS X+ Y+");

   pad->GetXaxis()->SetTitle("p_{T}");
   pad->GetYaxis()->SetTitle("1/(2 #pi p_{T} dp_{T}) d^{2} N / dp_{T} dy");
   
   gPad->SetLogy();
   
   spectraGr.SetMarkerStyle(20);
   spectraGr.SetMarkerColor(kRed-3);
   spectraGr.SetLineColor(kRed-3);
   
   fitFunc.SetLineColor(kRed+1);
   fitFunc.SetRange(0., 20.);
   
   spectraGr.Draw("P");
   fitFunc.Draw("SAME");
   
   system(("mkdir -p output/Fits/" + Par.collisionSystem).c_str());
   canv->SaveAs(("output/Fits/" + Par.collisionSystem + "/" + Par.particleName + ".png").c_str());

   Json::Value fitInfo, fitPar(Json::arrayValue);
   fitInfo["fit_function"] = Par.fitFunc;
   
   PrintInfo("Fit function: " + Par.fitFunc);
   PrintInfo("Fit parameters: ");
   for (int i = 0; i < fitFunc.GetNpar(); i++)
   {
      Print(fitFunc.GetParameter(i));
      fitPar.append(fitFunc.GetParameter(i));
   }
   fitInfo["fit_parameters"] = fitPar;

   const std::string outputFileName = "data/Spectra/" + Par.collisionSystem + 
                                      "/" + Par.particleName + "Fit.json";
   std::ofstream outputFile(outputFileName, std::ofstream::binary);
   outputFile << fitInfo;
}
