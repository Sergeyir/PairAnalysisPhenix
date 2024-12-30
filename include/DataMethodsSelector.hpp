// $HEADER$
//------------------------------------------------------------------------------------------------
//                             DataMethodsSelector class declaration
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

#ifndef DATA_METHODS_SELECTOR_HPP
#define DATA_METHODS_SELECTOR_HPP

#include "ErrorHandler.hpp"
#include "IOTools.hpp"

#include "DeadAreasCuts.h"
#include "SigmRes.h"

class DataMethodsSelector
{
   public:
   
   DataMethodsSelector();
   DataMethodsSelector(const std::string& datasetName, bool printInfo = false);

   void Init(const std::string& datasetName, bool printInfo = false);
   
   // cuts
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
   // sigmalized residuals
   double GetPC2SDPhi(const double dphi, const double pT, const int charge);
   double GetPC2SDZ(const double dz, const double pT, const int charge);
   double GetPC3SDPhi(const double phi, const double dphi, const double pT, const int charge);
   double GetPC3SDZ(const double phi, const double dz, const double pT, const int charge);
   double GetTOFwSDPhi(const double dphi, const double pT, const int charge);
   double GetTOFwSDZ(const double dz, const double pT, const int charge);
   double GetTOFeSDPhi(const double dphi, const double pT, const int charge);
   double GetTOFeSDZ(const double dz, const double pT, const int charge);
   double GetEMCSDPhi(const double phi, const double dphi, const double pT, 
                      const int charge, const int sect);
   double GetEMCSDZ(const double phi, const double dz, const double pT, 
                    const int charge, const int sect);

   virtual ~DataMethodsSelector();

   protected:

   bool (*fIsDeadDC)(const double, const double, const double, const double);
   bool (*fIsDeadPC1)(const double, const double, const double);
   bool (*fIsDeadPC2)(const double, const double);
   bool (*fIsDeadPC3)(const double, const double, const double);
   bool (*fIsDeadEMCal)(const double, const double, const int, const double, const double);
   bool (*fIsBadSlat)(const int);
   bool (*fIsBadStripTOFw)(const int);
   bool (*fIsDeadTOFe)(const double, const double, const double);
   bool (*fIsDeadTOFw)(const double, const double, const double);

   double (*fPC2SDPhi)(const double, const double, const int);
   double (*fPC2SDZ)(const double, const double, const int);
   double (*fPC3SDPhi)(const double, const double, const double, const int charge);
   double (*fPC3SDZ)(const double phi, const double dz, const double pT, const int charge);
   double (*fTOFwSDPhi)(const double dphi, const double pT, const int charge);
   double (*fTOFwSDZ)(const double dz, const double pT, const int charge);
   double (*fTOFeSDPhi)(const double dphi, const double pT, const int charge);
   double (*fTOFeSDZ)(const double dz, const double pT, const int charge);
   double (*fEMCSDPhi)(const double phi, const double dphi, const double pT, 
                   const int charge, const int sect);
   double (*fEMCSDZ)(const double phi, const double dz, const double pT, 
                    const int charge, const int sect);
};

#endif /* DATA_METHODS_SELECTOR_CPP */
