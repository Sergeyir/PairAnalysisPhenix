/** 
 *  @file   SimCalibrator.cpp 
 *  @brief  Contains realisation of class SimCalibrator
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/CalPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/

#ifndef SIM_CALIBRATOR_CPP
#define SIM_CALIBRATOR_CPP

#include "../include/SimCalibrator.hpp"

SimCalibrator::SimCalibrator() {}

SimCalibrator::SimCalibrator(const std::string& runName, const std::string& options)
{
   Initialize(runName, options);
}

void SimCalibrator::Initialize(const std::string& runName, const std::string& options)
{
   if (options.size() != 6)
   {
      CppTools::PrintError("SimCalibrator: options size is " + std::to_string(options.size()) + 
                           " while 6 has been expected");
   }

   if (options[0] == '1') useDC = true;
   else useDC = false;

   if (options[1] == '1') usePC1 = true;
   else usePC1 = false;

   if (options[2] == '1')
   {
      usePC2 = true;
      SetParameters("data/Parameters/SigmalizedResidualsSim/" + runName + "/PC2.txt", 
                    parMeansPC2, parSigmasPC2);
   }
   else usePC2 = false;

   if (options[3] == '1')
   {
      usePC3 = true;
      SetParameters("data/Parameters/SigmalizedResidualsSim/" + runName + "/PC3e.txt", 
                    parMeansPC3e, parSigmasPC3e);
      SetParameters("data/Parameters/SigmalizedResidualsSim/" + runName + "/PC3w.txt", 
                    parMeansPC3w, parSigmasPC3w);
   }
   else usePC3 = false;

   if (options[4] == '1')
   {
      useTOFe = true;
      SetParameters("data/Parameters/SigmalizedResidualsSim/" + runName + "/TOFe.txt", 
                    parMeansTOFe, parSigmasTOFe);
   }
   else useTOFe = false;

   if (options[5] == '1')
   {
      useTOFw = true;
      SetParameters("data/Parameters/SigmalizedResidualsSim/" + runName + "/TOFw.txt", 
                    parMeansTOFw, parSigmasTOFw);
   }
   else useTOFw = false;

   if (options[6] == '1')
   {
      useEMCal = true;
      for (int i = 0; i < 4; i++)
      {
         SetParameters("data/Parameters/SigmalizedResidualsSim/" + runName + 
                       "/EMCale" + std::to_string(i) + ".txt", 
                       parMeansEMCale[i], parSigmasEMCale[i]);
         SetParameters("data/Parameters/SigmalizedResidualsSim/" + runName + 
                       "/EMCalw" + std::to_string(i) + ".txt", 
                       parMeansEMCalw[i], parSigmasEMCalw[i]);
      }
   }
   else useEMCal = false;
}

void SimCalibrator::SetParameters(const std::string& inputFileName, 
                                  std::array<std::vector<bool>>& parMeans,
                                  std::array<std::vector<bool>>& parSigmas)
{
   CppTools::CheckInputFile(inputFileName);
   std::ifstream inputFile(inputFileName);

   bool isUsed;
   inputFile >> isUsed;

   if (!isUsed) 
   {
      CppTools::PrintError("SimCalibrator: Usage was specified as false in file " + inputFileName);
   }

   int numberOfParametersMeans;
   int numberOfParametersSigmas

   bool isUnexpectedEndOfFile = false;

   if (!(inputFile >> numberOfParametersMeans >> numberOfParametersSigmas))
   {
      isUnexpectedEndOfFile = true;
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
      CppTools::PrintError("SimCalibrator: Unexpected end of file: " + inputFileName);
   }
}

#endif /* SIM_CALIBRATOR_CPP */
