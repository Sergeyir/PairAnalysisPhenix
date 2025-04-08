/** 
 *  @file   DeadMapCutter.hpp 
 *  @brief  Contains declaration of class DeadMapCutter
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/CalPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef DEAD_MAP_CUTTER_HPP
#define DEAD_MAP_CUTTER_HPP

#include <string>
#include <fstream>

#include "ErrorHandler.hpp"

/*! @class DeadMapCutter
 * @brief Class DeadMapCutter provides simple means to implement and to use dead areas cuts (deadmaps) onto heatmaps
 */
class DeadMapCutter
{
   public:
   /*! @brief Constructor
    *
    * @param[in] runName name of the run
    * @param[in] options options that show which detectors deadmaps will be read and utilized
    *
    * Detectors in options go in the following order:
    *  -# DC
    *  -# PC1
    *  -# PC2
    *  -# PC3
    *  -# TOFe
    *  -# TOFw
    *  -# EMCal
    *
    * Options example: "110111" - this one uses all detectors apart from PC2
    */
   DeadMapCutter(const std::string& runName, const std::string& options = "1111111");

   /// Returns true if data in DC is in bad/dead area
   bool IsDeadDC(const double phi, const double alpha, const double zed);
   /// Returns true if data in PC1 is in bad/dead area
   bool IsDeadPC1(const double phi, const double ppc1z, const double ppc1phi);
   /// Returns true if data in PC2 is in bad/dead area
   bool IsDeadPC2(const double ppc2z, const double ppc2phi);
   /// Returns true if data in PC3 is in bad/dead area
   bool IsDeadPC3(const double phi, const double ppc3z, const double ppc3phi);
   /// Returns true if data in TOFe is in bad/dead area
   bool IsDeadTOFe(const double ptofy, const double ptofz, const int slat);
   /// Returns true if data in TOFw is in bad/dead area
   bool IsDeadTOFw(const double ptofwy, const double ptofwz, const int striptofw);
   /// Returns true if data in EMCal is in bad/dead area
   bool IsDeadEMCal(const double phi, const int sector, const int ytower, const int ztower);

   private:
   void SetDeadAreas(const std::string& inputFileName, std::vector<std::vector>>& deadAreas);
   /// cut areas for DCe, zDC<0
   std::vector<std::vector<bool>> cutAreasDCe0;
   /// cut areas for DCe, zDC>0
   std::vector<std::vector<bool>> cutAreasDCe1;
   /// cut areas for DCw, zDC<0
   std::vector<std::vector<bool>> cutAreasDCw0;
   /// cut areas for DCw, zDC>0
   std::vector<std::vector<bool>> cutAreasDCw1;
   /// cut areas for PC1e
   std::vector<std::vector<bool>> cutAreasPC1e;
   /// cut areas for PC1w
   std::vector<std::vector<bool>> cutAreasPC1w;
   /// cut areas for PC2
   std::vector<std::vector<bool>> cutAreasPC2;
   /// cut areas for PC3e
   std::vector<std::vector<bool>> cutAreasPC3e;
   /// cut areas for PC3w
   std::vector<std::vector<bool>> cutAreasPC3w;
   /// cut areas for TOFe
   std::vector<std::vector<bool>> cutAreasTOFe;
   /// cut slats in TOFe
   std::vector<bool> cutSlatsTOFe;
   /// cut areas for TOFw
   std::vector<std::vector<bool>> cutAreasTOFw;
   /// cut strips in TOFw
   std::vector<bool> cutStripsTOFw;
   /// cut areas for EMCale(0-3)
   std::vector<std::vector<bool>> EMCale[4];
   /// cut areas for EMCalw(0-3)
   std::vector<std::vector<bool>> EMCalw[4];
};

#endif /* DEAD_MAP_CUTTER_HPP */
