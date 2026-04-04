/** 
 *  @file   PainterHelper.cpp
 *  @brief  Contains realisation of class PainterHelper
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/CalPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef PAINTER_HELPER_CPP
#define PAINTER_HELPER_CPP

#include "PainterHelper.hpp"


PainterHelper::PainterHelper(TLegend *legend, const double markerSize, const int lineWidth, const double sysWidth)
{
   if (!legend) CppTools::PrintError("PainterHelper: legend pointer is nullptr");
   this->legend = legend;
   this->markerSize = markerSize;
   this->lineWidth = lineWidth;
   this->sysWidth = sysWidth;
}

void PainterHelper::DrawHist(TH1D *histogramWithStatErrors, TH1D *histogramWithSysErrors, 
                             const Color_t color, const double alpha, 
                             const Style_t markerStyle, const std::string& legendEntry)
{
   bool disableSysErrors = false;

   if (!histogramWithStatErrors) CppTools::PrintError("Histogram containing values and "\
                                                   "statistical uncertanities is nullptr");
   if (histogramWithSysErrors || 
       (histogramWithStatErrors->GetXaxis()->GetNbins() != 
        histogramWithStatErrors->GetXaxis()->GetNbins()))
   {
      disableSysErrors = true;
   }

   histogramWithStatErrors->SetMarkerStyle(markerStyle);
   histogramWithStatErrors->SetMarkerSize(markerSize);

   histogramWithStatErrors->SetMarkerColorAlpha(color, alpha);
   histogramWithStatErrors->SetLineColorAlpha(color, alpha);
   histogramWithStatErrors->SetLineWidth(lineWidth);

   histogramWithStatErrors->Clone()->Draw("SAME P E1 X0");

   if (!disableSysErrors)
   {
      TGraphErrors sysGr;

      for (int i = 1; i <= histogramWithSysErrors->GetXaxis()->GetNbins(); i++)
      {
         sysGr.AddPoint(histogramWithSysErrors->GetXaxis()->GetBinCenter(i), histogramWithSysErrors->GetBinContent(i));
         sysGr.SetPointError(i - 1, sysWidth, histogramWithSysErrors->GetBinError(i));
      }
      sysGr.SetLineWidth(lineWidth);
      sysGr.SetLineColorAlpha(color, alpha);
      sysGr.SetFillColorAlpha(0, 0.);

      sysGr.Clone()->Draw("5");
   }

   legend->AddEntry(histogramWithStatErrors->Clone(), legendEntry.c_str(), "P");
}

void PainterHelper::DrawGraph(TGraphErrors *graphWithStatErrors, TGraphErrors *graphWithSysErrors, 
                              const Color_t color, const double alpha, 
                              const Style_t markerStyle, const std::string& legendEntry)
{
   bool disableSysErrors = false;

   if (!graphWithStatErrors) CppTools::PrintError("Graph containing values and "\
                                                  "statistical uncertanities is nullptr");
   
   if (!graphWithSysErrors || graphWithStatErrors->GetN() != graphWithSysErrors->GetN()) 
   {
      disableSysErrors = true;
   }
   else
   {
      for (int i = 0; i < graphWithStatErrors->GetN(); i++)
      {
         if (fabs(graphWithStatErrors->GetPointX(i) - graphWithSysErrors->GetPointX(i)) > 
             fabs(graphWithStatErrors->GetPointX(i)/1e6))
         {
            disableSysErrors = true;
            break;
         }
         graphWithSysErrors->SetPointY(i, graphWithStatErrors->GetPointY(i));
      }
   }

   graphWithStatErrors->SetMarkerStyle(markerStyle);
   graphWithStatErrors->SetMarkerSize(markerSize);

   graphWithStatErrors->SetMarkerColorAlpha(color, alpha);
   graphWithStatErrors->SetLineColorAlpha(color, alpha);
   graphWithStatErrors->SetLineWidth(lineWidth);

   graphWithStatErrors->Clone()->Draw("P");

   if (!disableSysErrors)
   {
      graphWithSysErrors->SetLineWidth(lineWidth);
      graphWithSysErrors->SetLineColorAlpha(color, alpha);
      graphWithSysErrors->SetFillColorAlpha(0, 0.);

      graphWithSysErrors->Clone()->Draw("5");
   }

   legend->AddEntry(graphWithStatErrors->Clone(), legendEntry.c_str(), "P");
}

void PainterHelper::DrawGraphFromYAMLFile(const std::string& fileName, const std::string& qualifier, 
                                          const Color_t color, const double alpha, 
                                          const Style_t markerStyle, const std::string& legendEntry,
                                          const bool relativeErrors, const bool readSysErrors)
{
   TGraphErrors *graphWithSysErrors = nullptr;
   TGraphErrors *graphWithStatErrors = GetGraphFromYAMLFile(fileName, qualifier, graphWithSysErrors, 
                                                            relativeErrors, readSysErrors);
   DrawGraph(graphWithStatErrors, graphWithSysErrors, color, 
             alpha, markerStyle, legendEntry);
}

void PainterHelper::DrawGraphFromTXTFile(const std::string& fileName, 
                                         const Color_t color, const double alpha,
                                         const Style_t markerStyle, const std::string& legendEntry,
                                         const bool relativeErrors, const bool readSysErrors)
{
   TGraphErrors *graphWithSysErrors = nullptr;
   TGraphErrors *graphWithStatErrors = GetGraphFromTXTFile(fileName, graphWithSysErrors, 
                                                           relativeErrors, readSysErrors);
   DrawGraph(graphWithStatErrors, graphWithSysErrors, color, 
             alpha, markerStyle, legendEntry);
}

TGraphErrors *PainterHelper::GetGraphFromYAMLFile(const std::string& fileName, const std::string& qualifier, 
                                                  TGraphErrors*& graphWithSysErrors, 
                                                  const bool relativeErrors, const bool readSysErrors)
{
   InputYAMLReader yamlFileContents(fileName);

   for (const auto& dependentVariableEntry: yamlFileContents["dependent_variables"])
   {
      if (dependentVariableEntry["qualifiers"][0]["value"].as<std::string>() ==
          qualifier)
      {
         if (yamlFileContents["independent_variables"][0]["values"].size() != 
             dependentVariableEntry["values"].size()) 
         {
            CppTools::PrintError("Different amount of X and Y values in file " + 
                                 fileName + " for the field " + qualifier);
         }

         TGraphErrors *graphWithStatErrors = new TGraphErrors();
         if (!graphWithSysErrors) graphWithSysErrors = new TGraphErrors();
         else graphWithSysErrors->Clear();

         auto x = yamlFileContents["independent_variables"][0]["values"];
         auto y = dependentVariableEntry["values"];
         
         for (unsigned int i = 0; i < yamlFileContents["independent_variables"][0]["values"].size(); i++)
         {
            double xVal = 0.;
            if (x[i].size() == 1 || x[i].size() == 3)
            {
               xVal = x[i]["value"].as<double>();
            }
            else if (x[i].size() == 2)
            {
               xVal = (x[i]["low"].as<double>() + x[i]["high"].as<double>())/2.;
            }
            else CppTools::PrintError("Unexpected number of values of independed_variables columns in file " + fileName);

            graphWithStatErrors->AddPoint(xVal, y[i]["value"].as<double>());
            graphWithStatErrors->SetPointError(i, 0., y[i]["errors"][0]["symerror"].as<double>());

            if (readSysErrors)
            {
               graphWithSysErrors->AddPoint(xVal, y[i]["value"].as<double>());
               graphWithSysErrors->SetPointError(i, sysWidth, y[i]["errors"][1]["symerror"].as<double>());
            }
         }
         
         if (relativeErrors)
         {
            for (int i = 0; i < graphWithStatErrors->GetN(); i++)
            {
               graphWithStatErrors->SetPointError(i, 0., graphWithStatErrors->GetErrorY(i)*
                                                  graphWithStatErrors->GetPointY(i));

               if (readSysErrors)
               {
                  graphWithSysErrors->SetPointError(i, sysWidth, graphWithSysErrors->GetErrorY(i)*
                                                    graphWithSysErrors->GetPointY(i));
               }
            }
         }
         return graphWithStatErrors;
      }
   }
   CppTools::PrintError("No field qualifier \"" + qualifier + "\" found in file " + fileName);
   return nullptr;
}

TGraphErrors *PainterHelper::GetGraphFromTXTFile(const std::string& fileName, TGraphErrors*& graphWithSysErrors, 
                                                 const bool relativeErrors, const bool readSysErrors)
{
   CppTools::CheckInputFile(fileName);
   std::ifstream inputFile(fileName);

   TGraphErrors *graphWithStatErrors = new TGraphErrors();
   if (!graphWithSysErrors) graphWithSysErrors = new TGraphErrors();
   else graphWithSysErrors->Clear();

   double x, y, statError;
   int i = 0;
   while (inputFile >> x >> y >> statError)
   {
      if (relativeErrors) statError *= y;

      graphWithStatErrors->AddPoint(x, y);
      graphWithStatErrors->SetPointError(i, 0., statError);

      if (readSysErrors)
      {
         double sysError;
         if (!(inputFile >> sysError)) 
         {
            CppTools::PrintError("Unexpected end of file " + fileName + 
                                 ": expected systematic uncertainty value instead");
         }

         if (relativeErrors) sysError *= y;

         graphWithSysErrors->AddPoint(x, y);
         graphWithSysErrors->SetPointError(i, sysWidth, sysError);
      }

      i++;
   }

   return graphWithStatErrors;
}

void PainterHelper::DrawTypeCUncertainty(const double value, const double xPos, 
                                         const double yPos, const Color_t color, 
                                         const double alpha, const std::string& text)
{
   TGraphErrors gr;
   gr.AddPoint(xPos - sysWidth, yPos);
   gr.SetPointError(0, sysWidth, value);

   gr.SetLineColorAlpha(0, 0.);
   gr.SetFillColorAlpha(color, alpha);
   gr.SetFillStyle(1001);

   gr.Clone()->Draw("5");

   if (text != "")
   {
      TText ttext;

      ttext.SetTextFont(52);
      ttext.SetTextSize(0.05);
      ttext.SetTextAngle(90);

      ttext.DrawText(xPos, yPos + value, text.c_str());
   }
}

void PainterHelper::DrawLegend()
{
   legend->Draw();
}

void PainterHelper::SetMarkerSize(const double markerSize)
{

   this->markerSize = markerSize;
}

void PainterHelper::SetLineWidth(const int lineWidth)
{
   this->lineWidth = lineWidth;
}

void PainterHelper::SetSysWidth(const double sysWidth)
{
   this->sysWidth = sysWidth;
}

PainterHelper::~PainterHelper() {}


#endif /* PAINTER_HELPER_CPP */
