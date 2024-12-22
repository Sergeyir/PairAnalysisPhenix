// $SOURCE$
//------------------------------------------------------------------------------------------------
//                             DataCutsSelector class realisation
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

#ifndef DATA_CUT_SELECTOR_CPP
#define DATA_CUT_SELECTOR_CPP

#include "../include/DataCutsSelector.hpp"

DataCutsSelector::DataCutsSelector(const std::string &datasetName, bool printInfo)
{
   if (datasetName == "Run7AuAu200")
   {
      fIsDeadDC = &Run7AuAu200Cuts::IsDeadDC;
      fIsDeadPC1 = &Run7AuAu200Cuts::IsDeadPC1;
      fIsDeadPC2 = &Run7AuAu200Cuts::IsDeadPC2;
      fIsDeadPC3 = &Run7AuAu200Cuts::IsDeadPC3;
      fIsDeadEMCal = &Run7AuAu200Cuts::IsDeadEMCal;
      fIsBadSlat = &Run7AuAu200Cuts::IsBadSlat;
      fIsBadStripTOFw = &Run7AuAu200Cuts::IsBadStripTOFw;
      fIsDeadTOFe = &Run7AuAu200Cuts::IsDeadTOFe;
      fIsDeadTOFw = &Run7AuAu200Cuts::IsDeadTOFw;
   }
   else if (datasetName == "Run14HeAu200MB")
   {
      fIsDeadDC = &Run14HeAu200Cuts::IsDeadDC;
      fIsDeadPC1 = &Run14HeAu200Cuts::IsDeadPC1;
      fIsDeadPC2 = &Run14HeAu200Cuts::IsDeadPC2;
      fIsDeadPC3 = &Run14HeAu200Cuts::IsDeadPC3;
      fIsDeadEMCal = &Run14HeAu200Cuts::IsDeadEMCal;
      fIsBadSlat = &Run14HeAu200Cuts::IsBadSlat;
      fIsBadStripTOFw = &Run14HeAu200Cuts::IsBadStripTOFw;
      fIsDeadTOFe = &Run14HeAu200Cuts::IsDeadTOFe;
      fIsDeadTOFw = &Run14HeAu200Cuts::IsDeadTOFw;
   }
   else PrintError("Dataset name " + datasetName + " is not defined in DataCutsSelector");
}

bool DataCutsSelector::IsDeadDC(const double phi, const double zed, 
                                const double board, const double alpha)
{
   return fIsDeadDC(phi, zed, board, alpha);
}

bool DataCutsSelector::IsDeadPC1(const double phi, const double pc1z, const double pc1phi)
{
   return fIsDeadPC1(phi, pc1z, pc1phi);
}
bool DataCutsSelector::IsDeadPC2(const double pc2z, const double pc2phi)
{
   return fIsDeadPC2(pc2z, pc2phi);
}
bool DataCutsSelector::IsDeadPC3(const double phi, const double pc3z, const double pc3phi)
{
   return fIsDeadPC3(phi, pc3z, pc3phi);
}
bool DataCutsSelector::IsDeadEMCal(const double phi, const double zed, const double sect,
                                   const double pemcy, const double pemcz)
{
   return fIsDeadEMCal(phi, zed, sect, pemcy, pemcz);
}
bool DataCutsSelector::IsBadSlat(const int slat)
{
   return fIsBadSlat(slat);
}
bool DataCutsSelector::IsBadStripTOFw(const int strip)
{
   return fIsBadStripTOFw(strip);
}
bool DataCutsSelector::IsDeadTOFe(const double zed, const double tofy, const double tofz)
{
   return fIsDeadTOFe(zed, tofy, tofz);
}
bool DataCutsSelector::IsDeadTOFw(const double zed, const double board, const double alpha)
{
   return fIsDeadTOFw(zed, board, alpha);
}

DataCutsSelector::~DataCutsSelector() {}

#endif /* DATA_CUT_SELECTOR_CPP */
