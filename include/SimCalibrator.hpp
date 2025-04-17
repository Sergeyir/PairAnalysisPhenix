/** 
 *  @file   SimCalibrator.hpp 
 *  @brief  Contains declaration of class SimCalibrator
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/CalPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef SIM_CALIBRATOR_HPP
#define SIM_CALIBRATOR_HPP

#include <string>
#include <vector>
#include <fstream>

#include "ErrorHandler.hpp"

/*! @class SimCalibrator
 * @brief Class SimCalibrator provides simple means to implement and to use calibrations for analysis of Trees from PHENIX simulation (such as sigmalized residuals calibrations)
 */
class SimCalibrator
{
   public:
   ///@brief default constructor
   SimCalibrator();
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
   SimCalibrator(const std::string& runName, const std::string& options = "1111111");
   /*! @brief Initializes the object SimCalibrator
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

   /// @brief Returns sigmalized dphi from PC2
   double PC2SDPhi(const double dphi, const double pT, const int charge);
   /// @brief Returns sigmalized dz from PC2
   double PC2SDZ(const double dz, const double pT, const int charge);
   /// @brief Returns sigmalized dphi from PC3
   double PC3SDPhi(const double dphi, const double pT, const int charge, const int arm);
   /// @brief Returns sigmalized dz from PC3
   double PC3SDZ(const double dz, const double pT, const int charge, const int arm);
   /// @brief Returns sigmalized dphi from TOFe
   double TOFeSDPhi(const double dphi, const double pT, const int charge);
   /// @brief Returns sigmalized dz from TOFe
   double TOFeSDZ(const double dz, const double pT, const int charge);
   /// @brief Returns sigmalized dphi from TOFw
   double TOFwSDPhi(const double dphi, const double pT, const int charge);
   /// @brief Returns sigmalized dz from TOFw
   double TOFwSDZ(const double dz, const double pT, const int charge);
   /// @brief Returns sigmalized dphi from EMCal
   double EMCalSDPhi(const double dphi, const double pT, const int charge, const int arm, const int sector);
   /// @brief Returns sigmalized dz from EMCal
   double EMCalSDZ(const double dz, const double pT, const int charge, const int arm, const int sector);

   private:
   /// read 2D arrays from the file into class attributes
   void SetParameters(const std::string& inputFileName, 
                      std::array<std::vector<double>, 4>& parMeans,
                      std::array<std::vector<double>, 4>& parSigmas);
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
   /// fit parameters of dphi and dz means in PC2 for different charges 
   /// (indices: 0-dphi, charge>0; 1-dphi,charge<0, 2-dz,charge>0, 3-dz,charge<0)
   std::array<std::vector<double>, 4> parMeansPC2
   /// fit parameters of dphi and dz means in PC3e for different charges
   /// (indices: 0-dphi, charge>0; 1-dphi,charge<0, 2-dz,charge>0, 3-dz,charge<0)
   std::array<std::vector<double>, 4> parMeansPC3e
   /// fit parameters of dphi and dz means in PC3w for different charges
   /// (indices: 0-dphi, charge>0; 1-dphi,charge<0, 2-dz,charge>0, 3-dz,charge<0)
   std::array<std::vector<double>, 4> parMeansPC3w
   /// fit parameters of dphi and dz means in TOFe for different charges
   /// (indices: 0-dphi, charge>0; 1-dphi,charge<0, 2-dz,charge>0, 3-dz,charge<0)
   std::array<std::vector<double>, 4> parMeansTOFe
   /// fit parameters of dphi and dz means in TOFw for different charges
   /// (indices: 0-dphi, charge>0; 1-dphi,charge<0, 2-dz,charge>0, 3-dz,charge<0)
   std::array<std::vector<double>, 4> parMeansTOFw
   /// fit parameters of dphi and dz means in EMCale for different charges
   /// (indices: 0-dphi, charge>0; 1-dphi,charge<0, 2-dz,charge>0, 3-dz,charge<0)
   std::array<std::array<std::vector<double>, 4>, 4> parMeansEMCale
   /// fit parameters of dphi and dz means in EMCalw for different charges
   /// (indices: 0-dphi, charge>0; 1-dphi,charge<0, 2-dz,charge>0, 3-dz,charge<0)
   std::array<std::array<std::vector<double>, 4>, 4> parMeansEMCalw
   /// fit parameters of dphi and dz sigmas in PC2 for different charges 
   /// (indices: 0-dphi, charge>0; 1-dphi,charge<0, 2-dz,charge>0, 3-dz,charge<0)
   std::array<std::vector<double>, 4> parSigmasPC2
   /// fit parameters of dphi and dz sigmas in PC3e for different charges
   /// (indices: 0-dphi, charge>0; 1-dphi,charge<0, 2-dz,charge>0, 3-dz,charge<0)
   std::array<std::vector<double>, 4> parSigmasPC3e
   /// fit parameters of dphi and dz sigmas in PC3w for different charges
   /// (indices: 0-dphi, charge>0; 1-dphi,charge<0, 2-dz,charge>0, 3-dz,charge<0)
   std::array<std::vector<double>, 4> parSigmasPC3w
   /// fit parameters of dphi and dz sigmas in TOFe for different charges
   /// (indices: 0-dphi, charge>0; 1-dphi,charge<0, 2-dz,charge>0, 3-dz,charge<0)
   std::array<std::vector<double>, 4> parSigmasTOFe
   /// fit parameters of dphi and dz sigmas in TOFw for different charges
   /// (indices: 0-dphi, charge>0; 1-dphi,charge<0, 2-dz,charge>0, 3-dz,charge<0)
   std::array<std::vector<double>, 4> parSigmasTOFw
   /// fit parameters of dphi and dz sigmas in EMCale for different charges
   /// (indices: 0-dphi, charge>0; 1-dphi,charge<0, 2-dz,charge>0, 3-dz,charge<0)
   std::array<std::array<std::vector<double>, 4>, 4> parSigmasEMCale
   /// fit parameters of dphi and dz sigmas in EMCalw for different charges
   /// (indices: 0-dphi, charge>0; 1-dphi,charge<0, 2-dz,charge>0, 3-dz,charge<0)
   std::array<std::array<std::vector<double>, 4>, 4> parSigmasEMCalw
};

#endif /* SIM_CALIBRATOR_HPP */
