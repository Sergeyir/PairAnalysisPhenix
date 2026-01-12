/** 
 *  @file   SimM2Identificator.hpp
 *  @brief  Contains declaration of class SimM2Identificator
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef SIM_M2_IDENTIFICATOR_HPP
#define SIM_M2_IDENTIFICATOR_HPP

#include <fstream>
#include <string>
#include <iostream>
#include <cmath>

#include "TMath.h"

#include "ErrorHandler.hpp"
#include "MathTools.hpp"

#include "SingleTrackFunc.hpp"

/*! @class SimM2Identificator
 * @brief Class SimM2Identificator provides simple means to identifify charged hadrons (pions, kaons, protons, and antiprotons) via m2 distribution in TOFe, TOFw, and PbSc part of EMCal
 */
class SimM2Identificator
{
   public:
   ///@brief default constructor
   SimM2Identificator();
   /*! @brief Constructor
    *
    * @param[in] moduleName name of the CVS module specified in configure.in
    * @param[in] useEMCal specifies whether the data for EMCal m2 identification will be read an used
    *
    */
   SimM2Identificator(const std::string& runName, const bool useEMCal = false);
   /*! @brief Initializes the object SimM2Identificator
    *
    * @param[in] moduleName name of the CVS module specified in configure.in
    * @param[in] useEMCal specifies whether the data for EMCal m2 identification will be read an used
    *
    */
   void Initialize(const std::string& runName, const bool useEMCal = false);
   /*! @brief Returns the probability of a particle registered in TOFe being identified 
    * with the use of approximations of signals of charged hadrons from m2 distribution
    *
    * @param[in] id of a particle for which the probability will be calculated
    * @param[in] pT transverse momentum [GeV/c]
    * @param[in] sigmalizedExtrRange sigmalized extraction range in which identification is performed
    * @param[in] sigmalizedVetoRange sigmalized veto range from other particles in which identification is not allowed
    */
   double GetTOFeIdProb(const int id, const double pT,
                        const double sigmalizedExtrRange, const double sigmalizedVetoRange);
   /*! @brief Returns the probability of a particle registered in TOFw being identified 
    * with the use of approximations of signals of charged hadrons from m2 distribution
    *
    * @param[in] id of a particle for which the probability will be calculated
    * @param[in] pT transverse momentum [GeV/c]
    * @param[in] sigmalizedExtrRange sigmalized extraction range in which identification is performed
    * @param[in] sigmalizedVetoRange sigmalized veto range from other particles in which identification is not allowed
    */
   double GetTOFwIdProb(const int id, const double pT,
                        const double sigmalizedExtrRange, const double sigmalizedVetoRange);
   /*! @brief Returns the probability of a particle registered in EMCal being identified 
    * with the use of approximations of signals of charged hadrons from m2 distribution
    *
    * @param[in] dcarm arm of a detector (0 for east, 1 for west)
    * @param[in] sector sector of EMCal (2-3 for EMCale, 0-3 for EMCalw)
    * @param[in] id of a particle for which the probability will be calculated
    * @param[in] pT transverse momentum [GeV/c]
    * @param[in] sigmalizedExtrRange sigmalized extraction range in which identification is performed
    * @param[in] sigmalizedVetoRange sigmalized veto range from other particles in which identification is not allowed
    */
   double GetEMCalIdProb(const int dcarm, const int sector, const int id, const double pT, 
                         const double sigmalizedExtrRange, const double sigmalizedVetoRange);
   /// Destructor
   ~SimM2Identificator();

   private:
   /// reads parameters from files and stores them into arrays
   void SetParameters(const std::string& inputFileName, 
                      double (& parMean)[6][2], double (& parSigma)[5]);
   /// function that was used for approximation of m2 distirbutions means
   double GetM2Mean(const double pT, double *par);
   /// function that was used for approximation of m2 distirbutions sigmas
   double GetM2Sigma(const double pT, const double m2Mean, double *par);
   /// shows whether TOFe identification will be employed
   bool useTOFe = true;
   /// shows whether TOFw identification will be employed
   bool useTOFw = true;
   /// shows whether EMCal identification will be employed
   bool useEMCal;
   /// TOFe parameters for means
   double parMeanTOFe[6][2];
   /// TOFw parameters for means
   double parMeanTOFw[6][2];
   /// EMCale(2-3) parameters for means
   double parMeanEMCale[2][6][2];
   /// EMCalw(0-3) parameters for means
   double parMeanEMCalw[4][6][2];
   /// TOFe parameters for sigmas
   double parSigmaTOFe[5];
   /// TOFw parameters for sigmas
   double parSigmaTOFw[5];
   /// EMCale(2-3) parameters for sigmas
   double parSigmaEMCale[2][5];
   /// EMCalw(0-3) parameters for sigmas
   double parSigmaEMCalw[4][5];
};

#endif /* SIM_M2_IDENTIFICATOR_HPP */
