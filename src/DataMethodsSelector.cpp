// $SOURCE$
//------------------------------------------------------------------------------------------------
//                             DataMethodsSelector class realisation
//------------------------------------------------------------------------------------------------
// DataMethodsSelector
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic class used for selection of methods based on dataset
 **/
//------------------------------------------------------------------------------------------------

#ifndef DATA_METHODS_SELECTOR_CPP
#define DATA_METHODS_SELECTOR_CPP

#include "../include/DataMethodsSelector.hpp"

DataMethodsSelector::DataMethodsSelector() {}

DataMethodsSelector::DataMethodsSelector(const std::string &datasetName, bool printInfo)
{
   Init(datasetName, printInfo);
}

void DataMethodsSelector::Init(const std::string& datasetName, bool printInfo)
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

      fPC2SDPhi = &SimRun7AuAu200SigmRes::GetPC2SDPhi;
      fPC2SDZ = &SimRun7AuAu200SigmRes::GetPC2SDZ;
      fPC3SDPhi = &SimRun7AuAu200SigmRes::GetPC3SDPhi;
      fPC3SDZ = &SimRun7AuAu200SigmRes::GetPC3SDZ;
      fTOFeSDPhi = &SimRun7AuAu200SigmRes::GetTOFeSDPhi;
      fTOFeSDZ = &SimRun7AuAu200SigmRes::GetTOFeSDZ;
      fTOFwSDPhi = &SimRun7AuAu200SigmRes::GetTOFwSDPhi;
      fTOFwSDZ = &SimRun7AuAu200SigmRes::GetTOFwSDZ;
      fEMCSDPhi = &SimRun7AuAu200SigmRes::GetEMCSDPhi;
      fEMCSDZ = &SimRun7AuAu200SigmRes::GetEMCSDZ;
   }
   else if (datasetName == "Run14HeAu200")
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

      fPC2SDPhi = &SimRun14HeAu200SigmRes::GetPC2SDPhi;
      fPC2SDZ = &SimRun14HeAu200SigmRes::GetPC2SDZ;
      fPC3SDPhi = &SimRun14HeAu200SigmRes::GetPC3SDPhi;
      fPC3SDZ = &SimRun14HeAu200SigmRes::GetPC3SDZ;
      fTOFeSDPhi = &SimRun14HeAu200SigmRes::GetTOFeSDPhi;
      fTOFeSDZ = &SimRun14HeAu200SigmRes::GetTOFeSDZ;
      fTOFwSDPhi = &SimRun14HeAu200SigmRes::GetTOFwSDPhi;
      fTOFwSDZ = &SimRun14HeAu200SigmRes::GetTOFwSDZ;
      fEMCSDPhi = &SimRun14HeAu200SigmRes::GetEMCSDPhi;
      fEMCSDZ = &SimRun14HeAu200SigmRes::GetEMCSDZ;
   }
   else 
   {
      CppTools::PrintError("Dataset name " + datasetName + " is not defined in DataMethodsSelector");
   }
}

bool DataMethodsSelector::IsDeadDC(const double phi, const double zed, 
                                   const double board, const double alpha)
{
   return fIsDeadDC(phi, zed, board, alpha);
}

bool DataMethodsSelector::IsDeadPC1(const double phi, const double pc1z, const double pc1phi)
{
   return fIsDeadPC1(phi, pc1z, pc1phi);
}

bool DataMethodsSelector::IsDeadPC2(const double pc2z, const double pc2phi)
{
   return fIsDeadPC2(pc2z, pc2phi);
}

bool DataMethodsSelector::IsDeadPC3(const double phi, const double pc3z, const double pc3phi)
{
   return fIsDeadPC3(phi, pc3z, pc3phi);
}

bool DataMethodsSelector::IsDeadEMCal(const double phi, const double zed, const double sect,
                                      const double pemcy, const double pemcz)
{
   return fIsDeadEMCal(phi, zed, sect, pemcy, pemcz);
}

bool DataMethodsSelector::IsBadSlat(const int slat)
{
   return fIsBadSlat(slat);
}

bool DataMethodsSelector::IsBadStripTOFw(const int strip)
{
   return fIsBadStripTOFw(strip);
}

bool DataMethodsSelector::IsDeadTOFe(const double zed, const double tofy, const double tofz)
{
   return fIsDeadTOFe(zed, tofy, tofz);
}

bool DataMethodsSelector::IsDeadTOFw(const double zed, const double board, const double alpha)
{
   return fIsDeadTOFw(zed, board, alpha);
}

double DataMethodsSelector::GetPC2SDPhi(const double dphi, const double pT, const int charge)
{
   return fPC2SDPhi(dphi, pT, charge);
}

double DataMethodsSelector::GetPC2SDZ(const double dz, const double pT, const int charge)
{
   return fPC2SDZ(dz, pT, charge);
}

double DataMethodsSelector::GetPC3SDPhi(const double phi, const double dphi, const double pT, 
                                        const int charge)
{
   return fPC3SDPhi(phi, dphi, pT, charge);
}

double DataMethodsSelector::GetPC3SDZ(const double phi, const double dz, const double pT, 
                                      const int charge)
{
   return fPC3SDZ(phi, dz, pT, charge);
}

double DataMethodsSelector::GetTOFeSDPhi(const double dphi, const double pT, const int charge)
{
   return fTOFeSDPhi(dphi, pT, charge);
}

double DataMethodsSelector::GetTOFeSDZ(const double dz, const double pT, const int charge)
{
   return fTOFeSDZ(dz, pT, charge);
}

double DataMethodsSelector::GetTOFwSDPhi(const double dphi, const double pT, const int charge)
{
   return fTOFwSDPhi(dphi, pT, charge);
}

double DataMethodsSelector::GetTOFwSDZ(const double dz, const double pT, const int charge)
{
   return fTOFwSDZ(dz, pT, charge);
}

double DataMethodsSelector::GetEMCSDPhi(const double phi, const double dphi, const double pT, 
                                        const int charge, const int sect)
{
   return fEMCSDPhi(phi, dphi, pT, charge, sect);
}

double DataMethodsSelector::GetEMCSDZ(const double phi, const double dz, const double pT, 
                                      const int charge, const int sect)
{
   return fEMCSDZ(phi, dz, pT, charge, sect);
}

DataMethodsSelector::~DataMethodsSelector() {}

#endif /* DATA_METHODS_SELECTOR_CPP */
