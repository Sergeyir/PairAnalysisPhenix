/** 
 *  @file   SimM2Identificator.cpp
 *  @brief  Contains realisation of class SimM2Identificator
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/

#ifndef SIM_M2_IDENTIFICATOR_CPP
#define SIM_M2_IDENTIFICATOR_CPP

#include "SimM2Identificator.hpp"

SimM2Identificator::SimM2Identificator() {}

SimM2Identificator::SimM2Identificator(const std::string& runName, const bool useEMCal)
{
   Initialize(runName, useEMCal);
}

void SimM2Identificator::Initialize(const std::string& runName, const bool useEMCal)
{
   SetParameters("data/Parameters/M2Id/" + runName +"/M2ParTOFe.txt", parMeanTOFe, parSigmaTOFe);
   SetParameters("data/Parameters/M2Id/" + runName +"/M2ParTOFw.txt", parMeanTOFw, parSigmaTOFw);

   this->useEMCal = useEMCal;

   if (useEMCal)
   {
      for (int i = 0; i < 2; i++)
      {
         SetParameters("data/Parameters/M2Id/" + runName + "/M2ParEMCale" + 
                       std::to_string(i + 2) + ".txt", parMeanEMCale[i], parSigmaEMCale[i]);
      }
      for (int i = 0; i < 4; i++)
      {
         SetParameters("data/Parameters/M2Id/" + runName + "/M2ParEMCalw" + 
                       std::to_string(i) + ".txt", parMeanEMCalw[i], parSigmaEMCalw[i]);
      }
   }
   else 
   {
      std::cout << "\033[1m\033[32mInfo:\033[0m SimM2Identificator: "\
                   "m2 identification with EMCal was "\
                   "specified to be not initialized" << std::endl;
   }
}

double SimM2Identificator::GetTOFeIdProb(const int id, const double pT,
                                         const double sigmalizedExtrRange, 
                                         const double sigmalizedVetoRange)
{
   double meanPi, meanK, meanP;

   if (id > 0)
   {
      meanPi = GetM2Mean(pT, &parMeanTOFe[0][0]);
      meanK = GetM2Mean(pT, &parMeanTOFe[2][0]);
      meanP = GetM2Mean(pT, &parMeanTOFe[4][0]);
   }
   else
   {
      meanPi = GetM2Mean(pT, &parMeanTOFe[1][0]);
      meanK = GetM2Mean(pT, &parMeanTOFe[3][0]);
      meanP = GetM2Mean(pT, &parMeanTOFe[5][0]);
   }

   const double sigmaPi = GetM2Sigma(pT, meanPi, &parSigmaTOFe[0]);
   const double sigmaK = GetM2Sigma(pT, meanK, &parSigmaTOFe[0]);
   const double sigmaP = GetM2Sigma(pT, meanP, &parSigmaTOFe[0]);

   switch (abs(id))
   {
      case 211:
         if (meanK - sigmaK*sigmalizedVetoRange < meanPi - sigmalizedExtrRange*sigmaPi ||
             meanP - sigmaP*sigmalizedVetoRange < meanPi - sigmalizedExtrRange*sigmaPi) return 0.;
         return erf(sigmalizedExtrRange/TMath::Sqrt2())/2. + 
                erf((CppTools::Minimum(meanPi + sigmalizedExtrRange*sigmaPi, 
                                       meanK - sigmalizedVetoRange*sigmaK,
                                       meanP - sigmalizedVetoRange*sigmaP) - meanPi)/
                     sigmaPi/TMath::Sqrt2())/2.;
         break;
      case 321:
         if (meanPi + sigmaPi*sigmalizedVetoRange > meanK + sigmalizedExtrRange*sigmaK ||
             meanP - sigmaP*sigmalizedVetoRange < meanK - sigmalizedExtrRange*sigmaK ||
             meanPi + sigmaPi*sigmalizedVetoRange > meanP - sigmaP*sigmalizedVetoRange) return 0.;
         return erf((meanK - CppTools::Maximum(meanPi + sigmalizedVetoRange*sigmaPi,
                                               meanK - sigmalizedExtrRange*sigmaK))/
                    sigmaK/TMath::Sqrt2())/2. +
                erf((CppTools::Minimum(meanK + sigmalizedExtrRange*sigmaK,
                                       meanP - sigmalizedVetoRange*sigmaP) - meanK)/
                     sigmaK/TMath::Sqrt2())/2.;
         break;
      case 2212:
         if (meanPi + sigmaPi*sigmalizedVetoRange > meanP + sigmalizedExtrRange*sigmaP ||
             meanK - sigmaK*sigmalizedVetoRange > meanP - sigmalizedExtrRange*sigmaP) return 0.;
         return erf((meanP - CppTools::Maximum(meanPi + sigmalizedVetoRange*sigmaPi,
                                               meanK + sigmalizedVetoRange*sigmaK,
                                               meanP + sigmalizedExtrRange*sigmaP))/
                     sigmaP/TMath::Sqrt2())/2. + erf(sigmalizedExtrRange/TMath::Sqrt2())/2.;
         break;
   }
   return 0.;
}

double SimM2Identificator::GetTOFwIdProb(const int id, const double pT, 
                                         const double sigmalizedExtrRange, 
                                         const double sigmalizedVetoRange)
{
   double meanPi, meanK, meanP;
   if (id > 0)
   {
      meanPi = GetM2Mean(pT, &parMeanTOFw[0][0]);
      meanK = GetM2Mean(pT, &parMeanTOFw[2][0]);
      meanP = GetM2Mean(pT, &parMeanTOFw[4][0]);
   }
   else
   {
      meanPi = GetM2Mean(pT, &parMeanTOFw[1][0]);
      meanK = GetM2Mean(pT, &parMeanTOFw[3][0]);
      meanP = GetM2Mean(pT, &parMeanTOFw[5][0]);
   }

   const double sigmaPi = GetM2Sigma(pT, meanPi, &parSigmaTOFw[0]);
   const double sigmaK = GetM2Sigma(pT, meanK, &parSigmaTOFw[0]);
   const double sigmaP = GetM2Sigma(pT, meanP, &parSigmaTOFw[0]);

   switch (abs(id))
   {
      case 211:
         if (meanK - sigmaK*sigmalizedVetoRange < meanPi - sigmalizedExtrRange*sigmaPi ||
             meanP - sigmaP*sigmalizedVetoRange < meanPi - sigmalizedExtrRange*sigmaPi) return 0.;
         return erf(sigmalizedExtrRange/TMath::Sqrt2())/2. + 
                erf((CppTools::Minimum(meanPi + sigmalizedExtrRange*sigmaPi, 
                                       meanK - sigmalizedVetoRange*sigmaK,
                                       meanP - sigmalizedVetoRange*sigmaP) - meanPi)/
                     sigmaPi/TMath::Sqrt2())/2.;
         break;
      case 321:
         if (meanPi + sigmaPi*sigmalizedVetoRange > meanK + sigmalizedExtrRange*sigmaK ||
             meanP - sigmaP*sigmalizedVetoRange < meanK - sigmalizedExtrRange*sigmaK ||
             meanPi + sigmaPi*sigmalizedVetoRange > meanP - sigmaP*sigmalizedVetoRange) return 0.;
         return erf((meanK - CppTools::Maximum(meanPi + sigmalizedVetoRange*sigmaPi,
                                               meanK - sigmalizedExtrRange*sigmaK))/
                    sigmaK/TMath::Sqrt2())/2. +
                erf((CppTools::Minimum(meanK + sigmalizedExtrRange*sigmaK,
                                       meanP - sigmalizedVetoRange*sigmaP) - meanK)/
                     sigmaK/TMath::Sqrt2())/2.;
         break;
      case 2212:
         if (meanPi + sigmaPi*sigmalizedVetoRange > meanP + sigmalizedExtrRange*sigmaP ||
             meanK - sigmaK*sigmalizedVetoRange > meanP - sigmalizedExtrRange*sigmaP) return 0.;
         return erf((meanP - CppTools::Maximum(meanPi + sigmalizedVetoRange*sigmaPi,
                                               meanK + sigmalizedVetoRange*sigmaK,
                                               meanP + sigmalizedExtrRange*sigmaP))/
                     sigmaP/TMath::Sqrt2())/2. + erf(sigmalizedExtrRange/TMath::Sqrt2())/2.;
         break;
   }
   return 0.;
}

double SimM2Identificator::GetEMCalIdProb(const int dcarm, const int sector, 
                                          const int id, const double pT, 
                                          const double sigmalizedExtrRange, 
                                          const double sigmalizedVetoRange)
{
   if (!useEMCal) return 0.;

   double meanPi, meanK, meanP;
   double sigmaPi, sigmaK, sigmaP;
   if (dcarm == 0) // EMCale
   {
      if (sector < 2) return 0.; // no PbGl

      if (id > 0)
      {
         meanPi = GetM2Mean(pT, &parMeanEMCale[sector - 2][0][0]);
         meanK = GetM2Mean(pT, &parMeanEMCale[sector - 2][2][0]);
         meanP = GetM2Mean(pT, &parMeanEMCale[sector - 2][4][0]);
      }
      else
      {
         meanPi = GetM2Mean(pT, &parMeanEMCale[sector - 2][1][0]);
         meanK = GetM2Mean(pT, &parMeanEMCale[sector - 2][3][0]);
         meanP = GetM2Mean(pT, &parMeanEMCale[sector - 2][5][0]);
      }

      sigmaPi = GetM2Sigma(pT, meanPi, &parSigmaEMCale[sector - 2][0]);
      sigmaK = GetM2Sigma(pT, meanK, &parSigmaEMCale[sector - 2][0]);
      sigmaP = GetM2Sigma(pT, meanP, &parSigmaEMCale[sector - 2][0]);
   }
   else // EMCalw
   {
      if (id > 0)
      {
         meanPi = GetM2Mean(pT, &parMeanEMCalw[sector][0][0]);
         meanK = GetM2Mean(pT, &parMeanEMCalw[sector][2][0]);
         meanP = GetM2Mean(pT, &parMeanEMCalw[sector][4][0]);
      }
      else
      {
         meanPi = GetM2Mean(pT, &parMeanEMCalw[sector][1][0]);
         meanK = GetM2Mean(pT, &parMeanEMCalw[sector][3][0]);
         meanP = GetM2Mean(pT, &parMeanEMCalw[sector][5][0]);
      }

      sigmaPi = GetM2Sigma(pT, meanPi, &parSigmaEMCale[sector][0]);
      sigmaK = GetM2Sigma(pT, meanK, &parSigmaEMCale[sector][0]);
      sigmaP = GetM2Sigma(pT, meanP, &parSigmaEMCale[sector][0]);
   }

   switch (abs(id))
   {
      case 211:
         if (meanK - sigmaK*sigmalizedVetoRange < meanPi - sigmalizedExtrRange*sigmaPi ||
             meanP - sigmaP*sigmalizedVetoRange < meanPi - sigmalizedExtrRange*sigmaPi) return 0.;
         return erf(sigmalizedExtrRange/TMath::Sqrt2())/2. + 
                erf((CppTools::Minimum(meanPi + sigmalizedExtrRange*sigmaPi, 
                                       meanK - sigmalizedVetoRange*sigmaK,
                                       meanP - sigmalizedVetoRange*sigmaP) - meanPi)/
                     sigmaPi/TMath::Sqrt2())/2.;
         break;
      case 321:
         if (meanPi + sigmaPi*sigmalizedVetoRange > meanK + sigmalizedExtrRange*sigmaK ||
             meanP - sigmaP*sigmalizedVetoRange < meanK - sigmalizedExtrRange*sigmaK ||
             meanPi + sigmaPi*sigmalizedVetoRange > meanP - sigmaP*sigmalizedVetoRange) return 0.;
         return erf((meanK - CppTools::Maximum(meanPi + sigmalizedVetoRange*sigmaPi,
                                               meanK - sigmalizedExtrRange*sigmaK))/
                    sigmaK/TMath::Sqrt2())/2. +
                erf((CppTools::Minimum(meanK + sigmalizedExtrRange*sigmaK,
                                       meanP - sigmalizedVetoRange*sigmaP) - meanK)/
                     sigmaK/TMath::Sqrt2())/2.;
         break;
      case 2212:
         if (meanPi + sigmaPi*sigmalizedVetoRange > meanP + sigmalizedExtrRange*sigmaP ||
             meanK - sigmaK*sigmalizedVetoRange > meanP - sigmalizedExtrRange*sigmaP) return 0.;
         return erf((meanP - CppTools::Maximum(meanPi + sigmalizedVetoRange*sigmaPi,
                                               meanK + sigmalizedVetoRange*sigmaK,
                                               meanP + sigmalizedExtrRange*sigmaP))/
                     sigmaP/TMath::Sqrt2())/2. + erf(sigmalizedExtrRange/TMath::Sqrt2())/2.;
         break;
   }
   return 0.;
}

void SimM2Identificator::SetParameters(const std::string& inputFileName, 
                                       double (& parMean)[6][2], double (& parSigma)[5])
{
   CppTools::CheckInputFile(inputFileName);
   std::ifstream inputFile(inputFileName);

   bool isUnexpectedEndOfFile = false;

   for (int i = 0; i < 6; i++)
   {
      if (!(inputFile >> parMean[i][0] >> parMean[i][1]))
      {
         isUnexpectedEndOfFile = true;
      }
   }

   for (int i = 0; i < 5; i++)
   {
      if (!(inputFile >> parSigma[i])) isUnexpectedEndOfFile = true;
   }

   if (isUnexpectedEndOfFile) 
   {
      CppTools::PrintError("SimM2Identificator: Unexpected end of file: " + inputFileName);
   }
}

double SimM2Identificator::GetM2Mean(const double pT, double *par)
{
   return par[0] + pT*par[1];
}

double SimM2Identificator::GetM2Sigma(const double pT, const double m2Mean, double *par)
{
   return sqrt(par[0]*par[0]/(par[3]*par[3])*4.*m2Mean*m2Mean*pT*pT + 
               par[1]*par[1]/(par[3]*par[3])*4.*m2Mean*m2Mean*(1. + m2Mean/(pT*pT)) + 
               par[2]*par[2]*2.9972e-4*2.9972e-4/(par[4]*par[4])*4.*pT*pT*(m2Mean + pT*pT));
}

SimM2Identificator::~SimM2Identificator() {}

#endif /* SIM_M2_IDENTIFICATOR_CPP */
