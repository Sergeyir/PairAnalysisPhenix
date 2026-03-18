/** 
 *  @file   PainterHelper.hpp 
 *  @brief  Contains declaration of class PainterHelper
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/CalPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef PAINTER_HELPER_HPP
#define PAINTER_HELPER_HPP

#include <string>
#include <filesystem>

#include "TColor.h"
#include "TStyle.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"

#include "ErrorHandler.hpp"

#include "InputYAMLReader.hpp"


/*! @class PainterHelper
 * @brief Class PainterHelper is a helper class to work with histograms, graphs, .yaml, .txt files for drawing invariant p_{T| spectra and nuclear modification factors
 */
class PainterHelper
{
   public:

   ///@brief Default constructor
   PainterHelper();
   /*! @brief Constructor with parameters
    * @param[in] sysWidth with of systematic error box of each point. By default it is set to 0.1
    */
   PainterHelper(const double sysWidth);

   /*! @brief Adds the histogram to the object list to draw
    * @param[in] histogramWithStatErr histogram containing values with statistical uncertainties
    * @param[in] histogramWithSysErr histogram containing values with systematic uncertainties
    * @param[in] color color of markers and error boxes to be set for the current histogram
    * @param[in] markerStyle marker style to be set for the current histogram
    * @param[in] legendEntry legend entry for the current histogram
    */
   void AddHist(TH1D *histogramWithStatErr, TH1D *histogramWithSysErr, 
                const Color_t color, const Style_t markerStyle, 
                const std::string& legendEntry);
   /*! @brief Adds the graph to the object list to draw
    * @param[in] graphWithStatErr graph containing values with statistical uncertainties
    * @param[in] graphWithSysErr graph containing values with systematic uncertainties
    * @param[in] color color of markers and error boxes to be set for the current graph
    * @param[in] markerStyle marker style to be set for the current graph
    * @param[in] legendEntry legend entry for the current graph
    */
   void AddGraph(TGraphErrors *graphWithStatErr, TGraphErrors *graphWithSysErr, 
                 const Color_t color, const Style_t markerStyle, 
                 const std::string& legendEntry);
   /*! @brief Adds the histogram read from the specified .root file to the object list to draw
    * @param[in] fileName name of the .root file from which histogram will be read
    * @param[in] histogramName name of the histogram that will be read
    * @param[in] color color of markers and error boxes to be set for the current histogram
    * @param[in] markerStyle marker style to be set for the current histogram
    * @param[in] legendEntry legend entry for the current histogram
    */
   void AddHistFromROOTFile(const std::string& fileName, const std::string histogramName, 
                            const Color_t color, const Style_t markerStyle, 
                            const std::string& legendEntry);
   /*! @brief Adds the graph read from the specified .root file to the object list to draw. The graph contained in ROOT file must be TGraphPainter
    * @param[in] fileName name of the .root file from which graph will be read
    * @param[in] graphName name of the graph that will be read. The graph contained in ROOT file must be TGraphPainter
    * @param[in] color color of markers and error boxes to be set for the current graph
    * @param[in] markerStyle marker style to be set for the current graph
    * @param[in] legendEntry legend entry for the current graph
    */
   void AddGraphFromROOTFile(const std::string& fileName, const std::string graphName, 
                             const Color_t color, const Style_t markerStyle, 
                             const std::string& legendEntry);
   /*! @brief Adds the graph using the data specified in .yaml file
    * @param[in] fileName name of the file from which the data will be read
    * @param[in] qualifier string qualifier to find the needed data. This qualifier must be the field value of "[dependent_variables][i][qualifiers][0][value]", where i - arbitrary integer.
    * @param[in] color color of markers and error boxes to be set for the current graph
    * @param[in] markerStyle marker style to be set for the current graph
    * @param[in] legendEntry legend entry for the current graph. If emtpy legendEntry is specified (by default) the legend entry will be read from "[source]" field from the specified .yaml file
    */
   void AddFromYAMLFile(const std::string& fileName, const std::string& qualifier, 
                        const Color_t color, const Style_t markerStyle, 
                        const std::string& legendEntry = "");

   private:

   /// width of a systematic uncertainty box of each point. By default equals to 0.1
   sysWidth = 0.1;

   /// container for storing histograms containing values with statistical uncertainties
   std::vector histsWithStatErr<TH1D *>;
   /// container for storing histograms containing values with systematic uncertainties
   std::vector histsWithSysErr<TH1D *>;
   /// container for storing graphs containing values with statistical uncertainties
   std::vector graphsWithStatErr<TH1D *>;
   /// container for storing graphs containing values with systematic uncertainties
   std::vector graphsWithSysErr<TH1D *>;
}

#endif /* PAINTER_HELPER_HPP */
