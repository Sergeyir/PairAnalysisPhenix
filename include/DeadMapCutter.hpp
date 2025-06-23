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
#include <vector>
#include <fstream>

#include "ErrorHandler.hpp"
#include "IOTools.hpp"

/*! @class DeadMapCutter
 * @brief Class DeadMapCutter provides simple means to implement and to use bad/dead areas cuts (deadmaps) onto 2D heatmaps
 */
class DeadMapCutter
{
   public:
   ///@brief default constructor
   DeadMapCutter();
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
    * Options example: "1101111" - this one uses all detectors apart from PC2
    */
   DeadMapCutter(const std::string& runName, const std::string& options = "1111111");
   /*! @brief Initializes the object DeadMapCutter
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
    * Options example: "1101111" - this one uses all detectors apart from PC2
    */
   void Initialize(const std::string& runName, const std::string& options = "1111111");

   /// Returns true if data in DC is in bad/dead area
   bool IsDeadDC(const int dcarm, const double zDC, const double board, const double alpha);
   /// Returns true if data in PC1 is in bad/dead area
   bool IsDeadPC1(const int dcarm, const double ppc1z, const double ppc1phi);
   /// Returns true if data in PC2 is in bad/dead area
   bool IsDeadPC2(const double ppc2z, const double ppc2phi);
   /// Returns true if data in PC3 is in bad/dead area
   bool IsDeadPC3(const int dcarm, const double ppc3z, const double ppc3phi);
   /// Returns true if data in TOFe is in bad/dead area
   bool IsDeadTOFe(const int chamber, const int slat);
   /// Returns true if data in TOFw is in bad/dead area
   bool IsDeadTOFw(const int chamber, const int strip);
   /// Returns true if data in EMCal is in bad/dead area
   bool IsDeadEMCal(const int dcarm, const int sector, const int ytower, const int ztower);
   /// Returns true if timing distribution in TOFe is bad
   bool IsDeadTimingTOFe(const int chamber, const int slat);
   /// Returns true if timing distribution in TOFw is bad
   bool IsDeadTimingTOFw(const int chamber, const int strip);
   /// Returns true if timing distribution in EMCal is bad
   bool IsDeadTimingEMCal(const int dcarm, const int sector, const int ytower, const int ztower);

   private:
   /// read 2D arrays from the file into class attributes
   void SetDeadAreas(const std::string& inputFileName, std::vector<std::vector<bool>>& cutAreas,
                     double *ranges);
   /// shows whether option for DC was specified
   bool useDC;
   /// shows whether option for PC1 was specified
   bool usePC1;
   /// shows whether option for PC2 was specified
   bool usePC2;
   /// shows whether option for PC3 was specified
   bool usePC3;
   /// shows whether option for TOFe was specified
   bool useTOFe;
   /// shows whether option for TOFw was specified
   bool useTOFw;
   /// shows whether option for EMCal was specified
   bool useEMCal;
   /// cut areas for DCeX1, zDC>=0
   std::vector<std::vector<bool>> cutAreasDCe0X1;
   /// cut areas for DCeX1, zDC<0
   std::vector<std::vector<bool>> cutAreasDCe1X1;
   /// cut areas for DCwX1, zDC>=0
   std::vector<std::vector<bool>> cutAreasDCw0X1;
   /// cut areas for DCwX1, zDC<0
   std::vector<std::vector<bool>> cutAreasDCw1X1;
   /// cut areas for DCeX2, zDC>=0
   std::vector<std::vector<bool>> cutAreasDCe0X2;
   /// cut areas for DCeX2, zDC<0
   std::vector<std::vector<bool>> cutAreasDCe1X2;
   /// cut areas for DCwX2, zDC>=0
   std::vector<std::vector<bool>> cutAreasDCw0X2;
   /// cut areas for DCwX2, zDC<0
   std::vector<std::vector<bool>> cutAreasDCw1X2;
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
   /// cut areas for TOFw
   std::vector<std::vector<bool>> cutAreasTOFw;
   /// cut areas for EMCale(0-3)
   std::vector<std::vector<bool>> cutAreasEMCale[4];
   /// cut areas for EMCalw(0-3)
   std::vector<std::vector<bool>> cutAreasEMCalw[4];
   /// cut areas for timing for TOFe
   std::vector<std::vector<bool>> cutAreasTimingTOFe;
   /// cut areas for timing for TOFw
   std::vector<std::vector<bool>> cutAreasTimingTOFw;
   /// cut areas for timing for EMCale(0-3)
   std::vector<std::vector<bool>> cutAreasTimingEMCale[4];
   /// cut areas for timing for EMCalw(0-3)
   std::vector<std::vector<bool>> cutAreasTimingEMCalw[4];
   /// ranges of X [0], [1] and Y, [2], [3] axis for DCeX1, zDC>=0 heatmaps
   double cutAreasDCe0X1Range[4];
   /// ranges of X [0], [1] and Y, [2], [3] axis for DCeX1, zDC<0 heatmaps
   double cutAreasDCe1X1Range[4];
   /// ranges of X [0], [1] and Y, [2], [3] axis for DCwX1, zDC>=0 heatmaps
   double cutAreasDCw0X1Range[4];
   /// ranges of X [0], [1] and Y, [2], [3] axis for DCwX1, zDC<0 heatmaps
   double cutAreasDCw1X1Range[4];
   /// ranges of X [0], [1] and Y, [2], [3] axis for DCeX2, zDC>=0 heatmaps
   double cutAreasDCe0X2Range[4];
   /// ranges of X [0], [1] and Y, [2], [3] axis for DCeX2, zDC<0 heatmaps
   double cutAreasDCe1X2Range[4];
   /// ranges of X [0], [1] and Y, [2], [3] axis for DCwX2, zDC>=0 heatmaps
   double cutAreasDCw0X2Range[4];
   /// ranges of X [0], [1] and Y, [2], [3] axis for DCwX2, zDC<0 heatmaps
   double cutAreasDCw1X2Range[4];
   /// ranges of X [0], [1] and Y, [2], [3] axis for PC1e heatmaps
   double cutAreasPC1eRange[4];
   /// ranges of X [0], [1] and Y, [2], [3] axis for PC1w heatmaps
   double cutAreasPC1wRange[4];
   /// ranges of X [0], [1] and Y, [2], [3] axis for PC2 heatmaps
   double cutAreasPC2Range[4];
   /// ranges of X [0], [1] and Y, [2], [3] axis for PC3e heatmaps
   double cutAreasPC3eRange[4];
   /// ranges of X [0], [1] and Y, [2], [3] axis for PC3w heatmaps
   double cutAreasPC3wRange[4];
   /// ranges of X [0], [1] and Y, [2], [3] axis for TOFe heatmaps
   double cutAreasTOFeRange[4];
   /// ranges of X [0], [1] and Y, [2], [3] axis for TOFw heatmaps
   double cutAreasTOFwRange[4];
   /// ranges of X [0], [1] and Y, [2], [3] axis for EMCale(0-3) heatmaps
   double cutAreasEMCaleRange[4][4];
   /// ranges of X [0], [1] and Y, [2], [3] axis for EMCalw(0-3) heatmaps
   double cutAreasEMCalwRange[4][4];
   /// ranges of X [0], [1] and Y, [2], [3] axis for TOFe timing heatmaps
   double cutAreasTimingTOFeRange[4];
   /// ranges of X [0], [1] and Y, [2], [3] axis for TOFw timing heatmaps
   double cutAreasTimingTOFwRange[4];
   /// ranges of X [0], [1] and Y, [2], [3] axis for EMCale(0-3) timing heatmaps
   double cutAreasTimingEMCaleRange[4][4];
   /// ranges of X [0], [1] and Y, [2], [3] axis for EMCalw(0-3) timing heatmaps
   double cutAreasTimingEMCalwRange[4][4];
};

#endif /* DEAD_MAP_CUTTER_HPP */
