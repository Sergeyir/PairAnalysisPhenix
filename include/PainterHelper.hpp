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
    * @param[in] ySysWidth with of systematic error box. By default it is set to 0.1
    */
   PainterHelper(const double ySysWidth);

   /*! @brief Adds the histogram to the object list to draw
    * @param[in] histogram histogram to be added to the object list
    */
   void AddHist(TH1D *histogram);
   /*! @brief Adds the histogram read from the specified .root file to the object list to draw
    * @param[in] fileName name of the file from which histogram will be read
    * @param[in] histogramName name of the histogram that will be read
    */
   void AddHistFromFile(const std::string& fileName, const std::string histogramName);

   private:

   ySysWidth = 0.1;
}

#endif /* PAINTER_HELPER_HPP */
