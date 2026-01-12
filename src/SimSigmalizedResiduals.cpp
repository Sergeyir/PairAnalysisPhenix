/** 
 *  @file   SimSigmalizedResiduals.cpp 
 *  @brief  Contains realisation of class SimSigmalizedResiduals
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/CalPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/

#ifndef SIM_SIGMALIZED_RESIDUALS_CPP
#define SIM_SIGMALIZED_RESIDUALS_CPP

#include "SimSigmalizedResiduals.hpp"

SimSigmalizedResiduals::SimSigmalizedResiduals() {}

SimSigmalizedResiduals::SimSigmalizedResiduals(const std::string& runName, 
                                               const std::string& options)
{
   Initialize(runName, options);
}

void SimSigmalizedResiduals::Initialize(const std::string& runName, const std::string& options)
{
   if (options.size() != 7)
   {
      CppTools::PrintError("SimSigmalizedResiduals: options size is " + 
                           std::to_string(options.size()) + " while 7 has been expected");
   }

   this->runName = runName;

   if (options[2] == '1')
   {
      doCalPC2 = SetParameters("PC2", parMeansPC2, parSigmasPC2);
   }
   else 
   {
      CppTools::PrintInfo("SimSigmalizedResiduals: calibrations for PC2 were "\
                          "specified to not be initialized; setting rough "\
                          "estimation of sdphi and sdz for the PC2: "\
                          "sdz = dz/2., sdphi = dphi/0.002");
      doCalPC2 = false;
   }

   if (options[3] == '1')
   {
      doCalPC3 = (SetParameters("PC3e", parMeansPC3e, parSigmasPC3e) &&
                  SetParameters("PC3w", parMeansPC3w, parSigmasPC3w));
   }
   else 
   {
      CppTools::PrintInfo("SimSigmalizedResiduals: calibrations for PC3 were "\
                          "specified to not be initialized; setting rough "\
                          "estimation of sdphi and sdz for the PC3: "\
                          "sdz = dz/2., sdphi = dphi/0.002");
      doCalPC3 = false;
   }

   if (options[4] == '1')
   {
      doCalTOFe = SetParameters("TOFe", parMeansTOFe, parSigmasTOFe);
   }
   else 
   {
      CppTools::PrintInfo("SimSigmalizedResiduals: calibrations for TOFe were "\
                          "specified to not be initialized; setting rough "\
                          "estimation of sdphi and sdz for the TOFe: "\
                          "sdz = dz/2., sdphi = dphi/0.002");
      doCalTOFe = false;
   }

   if (options[5] == '1')
   {
      doCalTOFw = SetParameters("TOFw", parMeansTOFw, parSigmasTOFw);
   }
   else 
   {
      CppTools::PrintInfo("SimSigmalizedResiduals: calibrations for TOFw were "\
                          "specified to not be initialized; setting rough "\
                          "estimation of sdphi and sdz for the TOFw: "\
                          "sdz = dz/2., sdphi = dphi/0.002");
      doCalTOFw = false;
   }

   if (options[6] == '1')
   {
      doCalEMCal = true;
      for (int i = 0; i < 4 && doCalEMCal; i++)
      {
         doCalEMCal = (doCalEMCal && SetParameters("EMCale" + std::to_string(i), 
                                                   parMeansEMCale[i], parSigmasEMCale[i]));
      }
      for (int i = 0; i < 4 && doCalEMCal; i++)
      {
         doCalEMCal = (doCalEMCal && SetParameters("EMCalw" + std::to_string(i), 
                                                   parMeansEMCalw[i], parSigmasEMCalw[i]));
      }
   }
   else 
   {
      CppTools::PrintInfo("SimSigmalizedResiduals: calibrations for EMCal were "\
                          "specified to not be initialized; setting rough "\
                          "estimation of sdphi and sdz for the EMCal: "\
                          "sdz = dz/2., sdphi = dphi/0.002");
      doCalEMCal = false;
   }
}

double SimSigmalizedResiduals::PC2SDPhi(const double dphi, const double pT, const int charge)
{
   if (!doCalPC2) return dphi/0.002;

   double mean;
   double sigma;
   if (charge > 0)
   {
      mean = GetDValMean(pT, &parMeansPC2[0][0]);
      sigma = GetDValSigma(pT, &parSigmasPC2[0][0]);
   }
   else
   {
      mean = GetDValMean(pT, &parMeansPC2[1][0]);
      sigma = GetDValSigma(pT, &parSigmasPC2[1][0]);
   }
   return (mean - dphi)/sigma;
}

double SimSigmalizedResiduals::PC2SDZ(const double dz, const double pT, const int charge)
{
   if (!doCalPC2) return dz/2.;

   double mean;
   double sigma;
   if (charge > 0)
   {
      mean = GetDValMean(pT, &parMeansPC2[2][0]);
      sigma = GetDValSigma(pT, &parSigmasPC2[2][0]);
   }
   else
   {
      mean = GetDValMean(pT, &parMeansPC2[3][0]);
      sigma = GetDValSigma(pT, &parSigmasPC2[3][0]);
   }
   return (mean - dz)/sigma;
}

double SimSigmalizedResiduals::PC3SDPhi(const double dphi, const double pT, 
                                        const int charge, const int dcarm)
{
   if (!doCalPC3) return dphi/0.002;

   double mean;
   double sigma;

   if (dcarm == 0)
   {
      if (charge > 0)
      {
         mean = GetDValMean(pT, &parMeansPC3e[0][0]);
         sigma = GetDValSigma(pT, &parSigmasPC3e[0][0]);
      }
      else
      {
         mean = GetDValMean(pT, &parMeansPC3e[1][0]);
         sigma = GetDValSigma(pT, &parSigmasPC3e[1][0]);
      }
   }
   else
   {
      if (charge > 0)
      {
         mean = GetDValMean(pT, &parMeansPC3w[0][0]);
         sigma = GetDValSigma(pT, &parSigmasPC3w[0][0]);
      }
      else
      {
         mean = GetDValMean(pT, &parMeansPC3w[1][0]);
         sigma = GetDValSigma(pT, &parSigmasPC3w[1][0]);
      }
   }
   return (mean - dphi)/sigma;
}

double SimSigmalizedResiduals::PC3SDZ(const double dz, const double pT, 
                                      const int charge, const int dcarm)
{
   if (!doCalPC3) return dz/2.;

   double mean;
   double sigma;

   if (dcarm == 0)
   {
      if (charge > 0)
      {
         mean = GetDValMean(pT, &parMeansPC3e[2][0]);
         sigma = GetDValSigma(pT, &parSigmasPC3e[2][0]);
      }
      else
      {
         mean = GetDValMean(pT, &parMeansPC3e[3][0]);
         sigma = GetDValSigma(pT, &parSigmasPC3e[3][0]);
      }
   }
   else
   {
      if (charge > 0)
      {
         mean = GetDValMean(pT, &parMeansPC3w[2][0]);
         sigma = GetDValSigma(pT, &parSigmasPC3w[2][0]);
      }
      else
      {
         mean = GetDValMean(pT, &parMeansPC3w[3][0]);
         sigma = GetDValSigma(pT, &parSigmasPC3w[3][0]);
      }
   }
   return (mean - dz)/sigma;
}

double SimSigmalizedResiduals::TOFeSDPhi(const double dphi, const double pT, const int charge)
{
   if (!doCalTOFe) return dphi/0.002;

   double mean;
   double sigma;
   if (charge > 0)
   {
      mean = GetDValMean(pT, &parMeansTOFe[0][0]);
      sigma = GetDValSigma(pT, &parSigmasTOFe[0][0]);
   }
   else
   {
      mean = GetDValMean(pT, &parMeansTOFe[1][0]);
      sigma = GetDValSigma(pT, &parSigmasTOFe[1][0]);
   }
   return (mean - dphi)/sigma;
}

double SimSigmalizedResiduals::TOFeSDZ(const double dz, const double pT, const int charge)
{
   if (!doCalTOFe) return dz/2.;

   double mean;
   double sigma;
   if (charge > 0)
   {
      mean = GetDValMean(pT, &parMeansTOFe[2][0]);
      sigma = GetDValSigma(pT, &parSigmasTOFe[2][0]);
   }
   else
   {
      mean = GetDValMean(pT, &parMeansTOFe[3][0]);
      sigma = GetDValSigma(pT, &parSigmasTOFe[3][0]);
   }
   return (mean - dz)/sigma;
}

double SimSigmalizedResiduals::TOFwSDPhi(const double dphi, const double pT, const int charge)
{
   if (!doCalTOFw) return dphi/0.002;

   double mean;
   double sigma;
   if (charge > 0)
   {
      mean = GetDValMean(pT, &parMeansTOFw[0][0]);
      sigma = GetDValSigma(pT, &parSigmasTOFw[0][0]);
   }
   else
   {
      mean = GetDValMean(pT, &parMeansTOFw[1][0]);
      sigma = GetDValSigma(pT, &parSigmasTOFw[1][0]);
   }
   return (mean - dphi)/sigma;
}

double SimSigmalizedResiduals::TOFwSDZ(const double dz, const double pT, const int charge)
{
   if (!doCalTOFw) return dz/2.;

   double mean;
   double sigma;
   if (charge > 0)
   {
      mean = GetDValMean(pT, &parMeansTOFw[2][0]);
      sigma = GetDValSigma(pT, &parSigmasTOFw[2][0]);
   }
   else
   {
      mean = GetDValMean(pT, &parMeansTOFw[3][0]);
      sigma = GetDValSigma(pT, &parSigmasTOFw[3][0]);
   }
   return (mean - dz)/sigma;
}

double SimSigmalizedResiduals::EMCalSDPhi(const double dphi, double pT, const int charge, 
                                          const int dcarm, const int sector)
{
   if (!doCalEMCal) return dphi/0.002;

   double mean;
   double sigma;

   if (dcarm == 0)
   {
      if (charge > 0)
      {
         mean = GetDValMean(pT, &parMeansEMCale[sector][0][0]);
         sigma = GetDValSigma(pT, &parSigmasEMCale[sector][0][0]);
      }
      else
      {
         mean = GetDValMean(pT, &parMeansEMCale[sector][1][0]);
         sigma = GetDValSigma(pT, &parSigmasEMCale[sector][1][0]);
      }
   }
   else
   {
      if (charge > 0)
      {
         mean = GetDValMean(pT, &parMeansEMCalw[sector][0][0]);
         sigma = GetDValSigma(pT, &parSigmasEMCalw[sector][0][0]);
      }
      else
      {
         mean = GetDValMean(pT, &parMeansEMCalw[sector][1][0]);
         sigma = GetDValSigma(pT, &parSigmasEMCalw[sector][1][0]);
      }
   }
   return (mean - dphi)/sigma;
}

double SimSigmalizedResiduals::EMCalSDZ(const double dz, const double pT, const int charge, 
                                        const int dcarm, const int sector)
{
   if (!doCalEMCal) return dz/2.;

   double mean;
   double sigma;

   if (dcarm == 0)
   {
      if (charge > 0)
      {
         mean = GetDValMean(pT, &parMeansEMCale[sector][2][0]);
         sigma = GetDValSigma(pT, &parSigmasEMCale[sector][2][0]);
      }
      else
      {
         mean = GetDValMean(pT, &parMeansEMCale[sector][3][0]);
         sigma = GetDValSigma(pT, &parSigmasEMCale[sector][3][0]);
      }
   }
   else
   {
      if (charge > 0)
      {
         mean = GetDValMean(pT, &parMeansEMCalw[sector][2][0]);
         sigma = GetDValSigma(pT, &parSigmasEMCalw[sector][2][0]);
      }
      else
      {
         mean = GetDValMean(pT, &parMeansEMCalw[sector][3][0]);
         sigma = GetDValSigma(pT, &parSigmasEMCalw[sector][3][0]);
      }
   }
   return (mean - dz)/sigma;
}

double SimSigmalizedResiduals::GetDValMean(double pT, const double *par)
{
   // on pT > ~2.5-3 mean becomes constant
   if (pT > 3.) pT = 3.;
   return par[0] + par[1]/pT + par[2]/(pT*pT) + par[3]*pT;
}

double SimSigmalizedResiduals::GetDValSigma(double pT, const double *par)
{
   // on pT > ~2.5-3 sigma becomes constant
   if (pT > 3.) pT = 3.;
   return par[0] + par[1]/pT + par[2]*pT;
}

bool SimSigmalizedResiduals::SetParameters(const std::string& detectorName, 
                                           std::array<std::vector<double>, 4>& parMeans,
                                           std::array<std::vector<double>, 4>& parSigmas)
{
   const std::string inputFileName = "data/Parameters/CalibrateSimSigmalizedResiduals/" + 
                                     runName + "/" + detectorName + ".txt";

   if (!CppTools::FileExists(inputFileName))
   {
      CppTools::PrintWarning("SimSigmalizedResiduals: File " + inputFileName + 
                             " does not exists; setting rough estimation of "\
                             "sdphi and sdz for the " + detectorName + ": "\
                             "sdz = dz/2., sdphi = dphi/0.002");
      return false;
   }
   std::ifstream inputFile(inputFileName);

   bool isUsed;
   inputFile >> isUsed;

   if (!isUsed) 
   {
      CppTools::PrintError("SimSigmalizedResiduals: Usage was specified as false in file " + 
                           inputFileName);
   }

   int numberOfParametersMeans;
   int numberOfParametersSigmas;

   bool isUnexpectedEndOfFile = false;

   if (!(inputFile >> numberOfParametersMeans >> numberOfParametersSigmas))
   {
      isUnexpectedEndOfFile = true;
   }

   if (numberOfParametersMeans != 4)
   {
      CppTools::PrintError("SimSigmalizedResiduals: Mismatching number of parameters "\
                           "for means approximation in file " + inputFileName +
                           ": expected 4 while " + std::to_string(numberOfParametersMeans) + 
                           "were provided");
   }
   if (numberOfParametersSigmas != 3)
   {
      CppTools::PrintError("SimSigmalizedResiduals: Mismatching number of parameters "\
                           "for sigmas approximation in file " + inputFileName +
                           ": expected 3 while " + std::to_string(numberOfParametersSigmas) + 
                           "were provided");
   }

   for (int i = 0; i < 4; i++)
   {
      parMeans[i].resize(numberOfParametersMeans);
      parSigmas[i].resize(numberOfParametersSigmas);

      for (int j = 0; j < numberOfParametersMeans; j++)
      {
         double tmp;
         if (!(inputFile >> tmp)) 
         {
            isUnexpectedEndOfFile = true;
            break;
         }
         parMeans[i][j] = tmp;
      }
      for (int j = 0; j < numberOfParametersSigmas; j++)
      {
         double tmp;
         if (!(inputFile >> tmp)) 
         {
            isUnexpectedEndOfFile = true;
            break;
         }
         parSigmas[i][j] = tmp;
      }
   }

   if (isUnexpectedEndOfFile) 
   {
      CppTools::PrintError("SimSigmalizedResiduals: Unexpected end of file: " + inputFileName);
   }

   return true;
}

#endif /* SIM_SIGMALIZED_RESIDUALS_CPP */
