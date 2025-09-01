/** 
 *  @file   SimSigmalizedResiduals.hpp 
 *  @brief  Contains declaration of class SimSigmalizedResiduals
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/CalPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef SIM_SIGMALIZED_RESIDUALS_HPP
#define SIM_SIGMALIZED_RESIDUALS_HPP

#include <array>
#include <vector>
#include <string>
#include <fstream>

#include "ErrorHandler.hpp"
#include "IOTools.hpp"

/*! @class SimSigmalizedResiduals
 * @brief Class SimSigmalizedResiduals provides simple means to implement and to use calibrations for analysis of Trees from PHENIX simulation (such as sigmalized residuals calibrations)
 */
class SimSigmalizedResiduals
{
   public:
   ///@brief default constructor
   SimSigmalizedResiduals();
   /*! @brief Constructor
    *
    * @param[in] runName name of the run
    * @param[in] options options that show which detectors calibrations will be read and utilized
    *
    * If the detector calibration is not utilized approximate value of sigmalized residuals will be returned when called the corresponding method (mean = 0 for both dphi and dz, sigma=0.002 for dphi and sigma=2 for dz)
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
   SimSigmalizedResiduals(const std::string& runName, const std::string& options = "1111111");
   /*! @brief Initializes the object SimSigmalizedResiduals
    *
    * @param[in] runName name of the run
    * @param[in] options options that show which detectors cailbrations will be read and utilized
    *
    * If the detector calibration is not utilized approximate value of sigmalized residuals will be returned when called the corresponding method (mean = 0 for both dphi and dz, sigma=0.002 for dphi and sigma=2 for dz)
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
   /// @brief Returns sigmalized dphi from PC2
   double PC2SDPhi(const double dphi, const double pT, const int charge);
   /// @brief Returns sigmalized dz from PC2
   double PC2SDZ(const double dz, const double pT, const int charge);
   /// @brief Returns sigmalized dphi from PC3
   double PC3SDPhi(const double dphi, const double pT, const int charge, const int dcarm);
   /// @brief Returns sigmalized dz from PC3
   double PC3SDZ(const double dz, const double pT, const int charge, const int dcarm);
   /// @brief Returns sigmalized dphi from TOFe
   double TOFeSDPhi(const double dphi, const double pT, const int charge);
   /// @brief Returns sigmalized dz from TOFe
   double TOFeSDZ(const double dz, const double pT, const int charge);
   /// @brief Returns sigmalized dphi from TOFw
   double TOFwSDPhi(const double dphi, const double pT, const int charge);
   /// @brief Returns sigmalized dz from TOFw
   double TOFwSDZ(const double dz, const double pT, const int charge);
   /// @brief Returns sigmalized dphi from EMCal
   double EMCalSDPhi(const double dphi, const double pT, const int charge, 
                     const int dcarm, const int sector);
   /// @brief Returns sigmalized dz from EMCal
   double EMCalSDZ(const double dz, const double pT, const int charge, 
                   const int dcarm, const int sector);
   private:
   /// @brief Returns mean of track deviation (sdphi or sdz)
   inline double GetDValMean(double pT, const double *par);
   /// @brief Returns sigma of track deviation (sdphi or sdz)
   inline double GetDValSigma(double pT, const double *par);
   /// read 2D arrays from the file into class attributes
   void SetParameters(const std::string& inputFileName, 
                      std::array<std::vector<double>, 4>& parMeans,
                      std::array<std::vector<double>, 4>& parSigmas);
   /// shows whether option for DC was specified
   bool doCalDC;
   /// shows whether option for PC1 was specified
   bool doCalPC1;
   /// shows whether option for PC2 was specified
   bool doCalPC2;
   /// shows whether option for PC3 was specified
   bool doCalPC3;
   /// shows whether option for TOFe was specified
   bool doCalTOFe;
   /// shows whether option for TOFw was specified
   bool doCalTOFw;
   /// shows whether option for EMCal was specified
   bool doCalEMCal;
   /// fit parameters of dphi and dz means in PC2 for different charges 
   /// (indices: 0-dphi, charge>0; 1-dphi,charge<0, 2-dz,charge>0, 3-dz,charge<0)
   std::array<std::vector<double>, 4> parMeansPC2;
   /// fit parameters of dphi and dz means in PC3e for different charges
   /// (indices: 0-dphi, charge>0; 1-dphi,charge<0, 2-dz,charge>0, 3-dz,charge<0)
   std::array<std::vector<double>, 4> parMeansPC3e;
   /// fit parameters of dphi and dz means in PC3w for different charges
   /// (indices: 0-dphi, charge>0; 1-dphi,charge<0, 2-dz,charge>0, 3-dz,charge<0)
   std::array<std::vector<double>, 4> parMeansPC3w;
   /// fit parameters of dphi and dz means in TOFe for different charges
   /// (indices: 0-dphi, charge>0; 1-dphi,charge<0, 2-dz,charge>0, 3-dz,charge<0)
   std::array<std::vector<double>, 4> parMeansTOFe;
   /// fit parameters of dphi and dz means in TOFw for different charges
   /// (indices: 0-dphi, charge>0; 1-dphi,charge<0, 2-dz,charge>0, 3-dz,charge<0)
   std::array<std::vector<double>, 4> parMeansTOFw;
   /// fit parameters of dphi and dz means in EMCale for different charges
   /// (indices: 0-dphi, charge>0; 1-dphi,charge<0, 2-dz,charge>0, 3-dz,charge<0)
   std::array<std::array<std::vector<double>, 4>, 4> parMeansEMCale;
   /// fit parameters of dphi and dz means in EMCalw for different charges
   /// (indices: 0-dphi, charge>0; 1-dphi,charge<0, 2-dz,charge>0, 3-dz,charge<0)
   std::array<std::array<std::vector<double>, 4>, 4> parMeansEMCalw;
   /// fit parameters of dphi and dz sigmas in PC2 for different charges 
   /// (indices: 0-dphi, charge>0; 1-dphi,charge<0, 2-dz,charge>0, 3-dz,charge<0)
   std::array<std::vector<double>, 4> parSigmasPC2;
   /// fit parameters of dphi and dz sigmas in PC3e for different charges
   /// (indices: 0-dphi, charge>0; 1-dphi,charge<0, 2-dz,charge>0, 3-dz,charge<0)
   std::array<std::vector<double>, 4> parSigmasPC3e;
   /// fit parameters of dphi and dz sigmas in PC3w for different charges
   /// (indices: 0-dphi, charge>0; 1-dphi,charge<0, 2-dz,charge>0, 3-dz,charge<0)
   std::array<std::vector<double>, 4> parSigmasPC3w;
   /// fit parameters of dphi and dz sigmas in TOFe for different charges
   /// (indices: 0-dphi, charge>0; 1-dphi,charge<0, 2-dz,charge>0, 3-dz,charge<0)
   std::array<std::vector<double>, 4> parSigmasTOFe;
   /// fit parameters of dphi and dz sigmas in TOFw for different charges
   /// (indices: 0-dphi, charge>0; 1-dphi,charge<0, 2-dz,charge>0, 3-dz,charge<0)
   std::array<std::vector<double>, 4> parSigmasTOFw;
   /// fit parameters of dphi and dz sigmas in EMCale for different charges
   /// (indices: 0-dphi, charge>0; 1-dphi,charge<0, 2-dz,charge>0, 3-dz,charge<0)
   std::array<std::array<std::vector<double>, 4>, 4> parSigmasEMCale;
   /// fit parameters of dphi and dz sigmas in EMCalw for different charges
   /// (indices: 0-dphi, charge>0; 1-dphi,charge<0, 2-dz,charge>0, 3-dz,charge<0)
   std::array<std::array<std::vector<double>, 4>, 4> parSigmasEMCalw;
};

#endif /* SIM_SIGMALIZED_RESIDUALS_HPP */
