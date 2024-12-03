// $HEADER$
//------------------------------------------------------------------------------------------------
//                             DataCutsSelector class declaration
//------------------------------------------------------------------------------------------------
// DataCutsSelector
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic class used for selection of data cuts based on dataset
 **/
//------------------------------------------------------------------------------------------------

#ifndef DATA_CUT_SELECTOR_HPP
#define DATA_CUT_SELECTOR_HPP

#include "ErrorHandler.hpp"
#include "IOTools.hpp"

#include "DeadAreasCuts.h"

class DataCutsSelector
{
   public:

   DataCutsSelector(const std::string &datasetName, bool printInfo = false);

   bool IsDeadDC(const double phi, const double zed, const double board, const double alpha);
   bool IsDeadPC1(const double phi, const double pc1z, const double pc1phi);
   bool IsDeadPC2(const double pc2z, const double pc2phi);
   bool IsDeadPC3(const double phi, const double pc3z, const double pc3phi);
   bool IsDeadEMCal(const double phi, const double zed, const double sect,
                    const double pemcy, const double pemcz);
   bool IsBadSlat(const int slat);
   bool IsBadStripTOFw(const int strip);
   bool IsDeadTOFe(const double zed, const double tofy, const double tofz);
   bool IsDeadTOFw(const double zed, const double board, const double alpha);

   virtual ~DataCutsSelector();

   bool (*fIsDeadDC)(const double, const double, const double, const double);
   bool (*fIsDeadPC1)(const double, const double, const double);
   bool (*fIsDeadPC2)(const double, const double);
   bool (*fIsDeadPC3)(const double, const double, const double);
   bool (*fIsDeadEMCal)(const double, const double, const int, const double, const double);
   bool (*fIsBadSlat)(const int);
   bool (*fIsBadStripTOFw)(const int);
   bool (*fIsDeadTOFe)(const double, const double, const double);
   bool (*fIsDeadTOFw)(const double, const double, const double);
   
};

#endif /* DATA_CUT_SELECTOR_CPP */
