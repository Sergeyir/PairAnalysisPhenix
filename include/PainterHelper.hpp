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
#include "TLegend.h"
#include "TText.h"

#include "ErrorHandler.hpp"
#include "IOTools.hpp"

#include "InputYAMLReader.hpp"


/*! @class PainterHelper
 * @brief Class PainterHelper is a helper class to work with histograms, graphs, .yaml, .txt files for drawing invariant p_{T| spectra and nuclear modification factors
 */
class PainterHelper
{
   public:

   ///@brief Default constructor (forbidden)
   PainterHelper() = delete;
   /*! @brief Constructor with parameters
    * @param[in] legend pointer to the legend to which all painted histograms and graphs will be added
    * @param[in] markerSize size of the marker that will be set for all points of graphs and hisotgrams. By default this value is set to 1.
    * @param[in] lineWidth line width of error boxes and errr bars. By default this value is set to 1
    * @param[in] sysWidth width of systematic error box of each point. By default it is set to 0.1
    */
   PainterHelper(TLegend *legend, const double markerSize = 1., 
                 const int lineWidth = 1., const double sysWidth = 0.1);
   /*! @brief Draws the specified histogram. The histogram will be painted with the option "SAME P E0 X0"
    * @param[in] histogramWithStatErrors histogram containing values with statistical uncertainties
    * @param[in] histogramWithSysErrors histogram containing values with systematic uncertainties. If nullptr is specified systematic uncertainties will not be drawn
    * @param[in] color color of markers and error boxes to be set for the current histogram
    * @param[in] alpha alpha (opacity) for the color of markers and lines of each point
    * @param[in] markerStyle marker style to be set for the current histogram
    * @param[in] legendEntry legend entry for the current histogram
    */
   void DrawHist(TH1D *histogramWithStatErrors, TH1D *histogramWithSysErrors, 
                 const Color_t color, const double alpha, 
                 const Style_t markerStyle, const std::string& legendEntry);
   /*! @brief Draws the specified graph
    * @param[in] graphWithStatErrors graph containing values with statistical uncertainties
    * @param[in] graphWithSysErrors graph containing values with systematic uncertainties
    * @param[in] color color of markers and error boxes to be set for the current graph
    * @param[in] alpha alpha (opacity) for the color of markers and lines of each point
    * @param[in] markerStyle marker style to be set for the current graph
    * @param[in] legendEntry legend entry for the current graph
    */
   void DrawGraph(TGraphErrors *graphWithStatErrors, TGraphErrors *graphWithSysErrors, 
                  const Color_t color, const double alpha, 
                  const Style_t markerStyle, const std::string& legendEntry);
   /*! @brief Draws the graph obtained from the data read from the spcified .yaml file
    * @param[in] fileName name of the file from which the data will be read
    * @param[in] qualifier string qualifier to find the needed data. This qualifier must be the field value of "[dependent_variables][i][qualifiers][0][value]", where i - arbitrary integer.
    * @param[in] color color of markers and error boxes to be set for the current graph
    * @param[in] alpha alpha (opacity) for the color of markers and lines of each point
    * @param[in] markerStyle marker style to be set for the current graph
    * @param[in] legendEntry legend entry for the current graph. If emtpy legendEntry is specified (by default) the legend entry will be read from "[source]" field from the specified .yaml file
    * @param[in] relativeErrors shows whether uncertainties will be read as absolute or relative values. By default absolute values will be read
    * @param[in] readSysErrors shows whether systematica uncertainties are written in the file 
    */
   void DrawGraphFromYAMLFile(const std::string& fileName, const std::string& qualifier, 
                              const Color_t color, const double alpha, const Style_t markerStyle, 
                              const std::string& legendEntry, 
                              const bool relativeErrors = false, const bool readSysErrors = true);
   /*! @brief Draws the graph obtained from the date read from the specified .txt file
    * @param[in] fileName name of the .txt file from which data will be read
    * @param[in] color color of markers and error boxes to be set for the current graph
    * @param[in] alpha alpha (opacity) for the color of markers and lines of each point
    * @param[in] markerStyle marker style to be set for the current graph
    * @param[in] legendEntry legend entry for the current graph
    * @param[in] relativeErrors shows whether uncertainties will be read as absolute or relative values. By default absolute values will be read
    * @param[in] readSysErrors shows whether systematica uncertainties are written in the file 
    */
   void DrawGraphFromTXTFile(const std::string& fileName, 
                             const Color_t color, const double alpha, const Style_t markerStyle, 
                             const std::string& legendEntry, const bool relativeErrors = false,
                             const bool readSysErrors = true);
   /*! @brief Draws the legend of already painted graphs
    * @param[in] fileName name of the .root file from which graph will be read
    * @param[in] graphWithSysErrors graph containing values with systematic uncertainties
    * @param[in] relativeErrors shows whether uncertainties will be read as absolute or relative values. By default absolute values will be read
    * @param[in] readSysErrors shows whether systematica uncertainties are written in the file 
    */
   TGraphErrors *GetGraphFromYAMLFile(const std::string& fileName, const std::string& qualifier, 
                                      TGraphErrors*& graphWithSysErrors, const bool relativeErrors = false,
                                      const bool readSysErrors = true);
   /*! @brief Draws the graph read from the specified .txt file to the object list to draw
    * @param[in] fileName name of the .root file from which graph will be read
    * @param[in] graphWithSysErrors graph containing values with systematic uncertainties
    * @param[in] relativeErrors shows whether uncertainties will be read as absolute or relative values. By default absolute values will be read
    * @param[in] readSysErrors shows whether systematica uncertainties are written in the file 
    */
   TGraphErrors *GetGraphFromTXTFile(const std::string& fileName, TGraphErrors*& graphWithSysErrors, 
                                     const bool relativeErrors = false, const bool readSysErrors = true);
   /*! @brief Draws the graph containing type c systematic uncertainty box
    * @param[in] value absolute value of the uncertainty
    * @param[in] xPos x position of the center of uncertainty box
    * @param[in] yPos y position of the center of uncertainty box
    * @param[in] color color of the filled area of the systematic uncertainty box
    * @param[in] alpha alpha for the color of the filled area of the systematic uncertainty box
    * @param[in] text text that will be displayed above the uncertainty box (text will be rotated 90 degrees to avoid overlapping)
    */
   void DrawTypeCUncertainty(const double value, const double xPos, const double yPos,
                             const Color_t color, const double alpha, 
                             const std::string& text = "");
   /// Draws the legend containing all previously drawn histograms and graphs
   void DrawLegend();
   /// Sets the size of the marker that will be set for all points of graphs and hisotgrams. By default this value is set to 1.
   void SetMarkerSize(const double markerSize);
   /// Sets the line width of error boxes and errr bars. By default this value is set to 1
   void SetLineWidth(const int lineWidth);
   /// sets the width of systematic error box of each point. By default it is set to 0.1
   void SetSysWidth(const double sysWidth);
   /// Default destructor
   ~PainterHelper();

   private:

   /// legend containing entries of all histograms and graphs drawn
   TLegend *legend;
   /// marker size of each point of histograms and graphs, By Default equals to 1.
   double markerSize;
   /// line wisth of each point error bars and boxes, By Default equals to 1
   double lineWidth;
   /// width of a systematic uncertainty box of each point. By default equals to 0.1
   double sysWidth;
};

#endif /* PAINTER_HELPER_HPP */
