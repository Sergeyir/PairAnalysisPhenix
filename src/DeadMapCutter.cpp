/** 
 *  @file   DeadMapCutter.cpp 
 *  @brief  Contains realisation of class DeadMapCutter
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/CalPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/

#ifndef DEAD_MAP_CUTTER_CPP
#define DEAD_MAP_CUTTER_CPP

#include "../include/DeadMapCutter.hpp"

DeadMapCutter::DeadMapCutter() {}

DeadMapCutter::DeadMapCutter(const std::string& runName, const std::string& options)
{
   Initialize(runName, options);
}

void DeadMapCutter::Initialize(const std::string& runName, const std::string& options)
{
   if (options.size() != 7)
   {
      CppTools::PrintError("DeadMapCutter: options size is " + std::to_string(options.size()) + 
                           " while 7 has been expected");
   }

   const std::string deadmapsDir = "data/Parameters/Deadmaps/" + runName + "/";
   const std::string timingDeadmapsDir = "data/Parameters/TimingDeadmaps/" + runName + "/";

   if (options[0] == '1')
   {
      useDC = (SetDeadAreas(deadmapsDir + "DCeX1_0.txt", cutAreasDCe0X1, cutAreasDCe0X1Range) &&
               SetDeadAreas(deadmapsDir + "DCeX1_1.txt", cutAreasDCe1X1, cutAreasDCe1X1Range) &&
               SetDeadAreas(deadmapsDir + "DCwX1_0.txt", cutAreasDCw0X1, cutAreasDCw0X1Range) &&
               SetDeadAreas(deadmapsDir + "DCwX1_1.txt", cutAreasDCw1X1, cutAreasDCw1X1Range) &&
               SetDeadAreas(deadmapsDir + "DCeX2_0.txt", cutAreasDCe0X2, cutAreasDCe0X2Range) &&
               SetDeadAreas(deadmapsDir + "DCeX2_1.txt", cutAreasDCe1X2, cutAreasDCe1X2Range) &&
               SetDeadAreas(deadmapsDir + "DCwX2_0.txt", cutAreasDCw0X2, cutAreasDCw0X2Range) &&
               SetDeadAreas(deadmapsDir + "DCwX2_1.txt", cutAreasDCw1X2, cutAreasDCw1X2Range));
   }
   else CppTools::PrintInfo("DeadMapCutter: Cuts for DC were specified to be not initialized");
   if (!useDC) CppTools::PrintInfo("DeadMapCutter: No cuts for DC will be applied");

   if (options[1] == '1')
   {
      usePC1 = (SetDeadAreas(deadmapsDir + "PC1e.txt", cutAreasPC1e, cutAreasPC1eRange) &&
                SetDeadAreas(deadmapsDir + "PC1w.txt", cutAreasPC1w, cutAreasPC1wRange));
   }
   else CppTools::PrintInfo("DeadMapCutter: Cuts for PC1 were specified to be not initialized");
   if (!usePC1) CppTools::PrintInfo("DeadMapCutter: No cuts for PC1 will be applied");

   if (options[2] == '1')
   {
      usePC2 = SetDeadAreas(deadmapsDir + "PC2.txt", cutAreasPC2, cutAreasPC2Range);
   }
   else CppTools::PrintInfo("DeadMapCutter: Cuts for PC2 were specified to be not initialized");
   if (!usePC2) CppTools::PrintInfo("DeadMapCutter: No cuts for PC2 will be applied");

   if (options[3] == '1')
   {
      usePC3 = (SetDeadAreas(deadmapsDir + "PC3e.txt", cutAreasPC3e, cutAreasPC3eRange) && 
                SetDeadAreas(deadmapsDir + "PC3w.txt", cutAreasPC3w, cutAreasPC3wRange));
   }
   else CppTools::PrintInfo("DeadMapCutter: Cuts for PC3 were specified to be not initialized");
   if (!usePC3) CppTools::PrintInfo("DeadMapCutter: No cuts for PC3 will be applied");

   if (options[4] == '1')
   {
      useTOFe = SetDeadAreas(deadmapsDir + "TOFe.txt", cutAreasTOFe, cutAreasTOFeRange);
      useTOFeTiming = SetDeadAreas(timingDeadmapsDir + "TimingDeadmapTOFe.txt", 
                                   cutAreasTimingTOFe, cutAreasTimingTOFeRange);
   }
   else CppTools::PrintInfo("DeadMapCutter: Cuts for TOFe were specified to be not initialized");
   if (!useTOFe) CppTools::PrintInfo("DeadMapCutter: No cuts for TOFe will be applied");
   if (!useTOFeTiming) CppTools::PrintInfo("DeadMapCutter: No cuts for TOFe will be "\
                                           "applied for improved timing");

   if (options[5] == '1')
   {
      useTOFw = SetDeadAreas(deadmapsDir + "TOFw.txt", cutAreasTOFw, cutAreasTOFwRange);
      useTOFwTiming = SetDeadAreas(timingDeadmapsDir + "TimingDeadmapTOFw.txt", 
                                   cutAreasTimingTOFw, cutAreasTimingTOFwRange);
   }
   else CppTools::PrintInfo("DeadMapCutter: Cuts for TOFe were specified to be not initialized");
   if (!useTOFw) CppTools::PrintInfo("DeadMapCutter: No cuts for TOFe will be applied");
   if (!useTOFwTiming) CppTools::PrintInfo("DeadMapCutter: No cuts for TOFe will be "\
                                           "applied for improved timing");

   if (options[6] == '1')
   {
      useEMCal = true;
      for (int i = 0; i < 4 && useEMCal; i++)
      {
         useEMCal = (SetDeadAreas(deadmapsDir + "EMCale" + std::to_string(i) + ".txt", 
                                  cutAreasEMCale[i], cutAreasEMCaleRange[i]) &&
                     SetDeadAreas(deadmapsDir + "EMCalw" + std::to_string(i) + ".txt", 
                                  cutAreasEMCalw[i], cutAreasEMCalwRange[i]));
      }
      useEMCalTiming = true;
      for (int i = 0; i < 4 && useEMCalTiming; i++)
      {
         useEMCalTiming = (SetDeadAreas(timingDeadmapsDir + "TimingDeadmapEMCale" + 
                                        std::to_string(i) + ".txt", 
                                        cutAreasTimingEMCale[i], cutAreasTimingEMCaleRange[i]) &&
                           SetDeadAreas(timingDeadmapsDir + "/TimingDeadmapEMCalw" + 
                                        std::to_string(i) + ".txt", 
                                        cutAreasTimingEMCalw[i], cutAreasTimingEMCalwRange[i]));
      }
   }
   else CppTools::PrintInfo("DeadMapCutter: Cuts for EMCal were specified to be not initialized");
   if (!useEMCal) CppTools::PrintInfo("DeadMapCutter: No cuts for EMCal will be applied");
   if (!useEMCalTiming) CppTools::PrintInfo("DeadMapCutter: No cuts for EMCal will be "\
                                           "applied for improved timing");
}

bool DeadMapCutter::IsDeadDC(const int dcarm, const double zDC, 
                             const double board, const double alpha)
{
   if (!useDC) return false;
   if (dcarm == 0) // DCe
   {
      if (zDC >= 0)
      {
         if (board <= cutAreasDCe0X1Range[0] || board >= cutAreasDCe0X1Range[1] ||
             alpha <= cutAreasDCe0X1Range[2] || alpha >= cutAreasDCe0X1Range[3] ||
             board <= cutAreasDCe0X2Range[0] || board >= cutAreasDCe0X2Range[1] ||
             alpha <= cutAreasDCe0X2Range[2] || alpha >= cutAreasDCe0X2Range[3]) return true;

         const short yBinX1 = static_cast<short>((alpha - cutAreasDCe0X1Range[2])/
                                                (cutAreasDCe0X1Range[3] - cutAreasDCe0X1Range[2])*
                                                static_cast<double>(cutAreasDCe0X1.size()));
         const short xBinX1 = static_cast<short>((board - cutAreasDCe0X1Range[0])/
                                                (cutAreasDCe0X1Range[1] - cutAreasDCe0X1Range[0])*
                                                static_cast<double>(cutAreasDCe0X1[0].size()));
         const short yBinX2 = static_cast<short>((alpha - cutAreasDCe0X2Range[2])/
                                                (cutAreasDCe0X2Range[3] - cutAreasDCe0X2Range[2])*
                                                static_cast<double>(cutAreasDCe0X2.size()));
         const short xBinX2 = static_cast<short>((board - cutAreasDCe0X2Range[0])/
                                                (cutAreasDCe0X2Range[1] - cutAreasDCe0X2Range[0])*
                                                static_cast<double>(cutAreasDCe0X2[0].size()));

         return (cutAreasDCe0X1[yBinX1][xBinX1] || cutAreasDCe0X2[yBinX2][xBinX2]);
      }
      else
      {
         if (board <= cutAreasDCe1X1Range[0] || board >= cutAreasDCe1X1Range[1] ||
             alpha <= cutAreasDCe1X1Range[2] || alpha >= cutAreasDCe1X1Range[3] ||
             board <= cutAreasDCe1X2Range[0] || board >= cutAreasDCe1X2Range[1] ||
             alpha <= cutAreasDCe1X2Range[2] || alpha >= cutAreasDCe1X2Range[3]) return true;

         const short yBinX1 = static_cast<short>((alpha - cutAreasDCe1X1Range[2])/
                                                (cutAreasDCe1X1Range[3] - cutAreasDCe1X1Range[2])*
                                                static_cast<double>(cutAreasDCe1X1.size()));
         const short xBinX1 = static_cast<short>((board - cutAreasDCe1X1Range[0])/
                                                (cutAreasDCe1X1Range[1] - cutAreasDCe1X1Range[0])*
                                                static_cast<double>(cutAreasDCe1X1[0].size()));
         const short yBinX2 = static_cast<short>((alpha - cutAreasDCe1X2Range[2])/
                                                (cutAreasDCe1X2Range[3] - cutAreasDCe1X2Range[2])*
                                                static_cast<double>(cutAreasDCe1X2.size()));
         const short xBinX2 = static_cast<short>((board - cutAreasDCe1X2Range[0])/
                                                (cutAreasDCe1X2Range[1] - cutAreasDCe1X2Range[0])*
                                                static_cast<double>(cutAreasDCe1X2[0].size()));

         return (cutAreasDCe1X1[yBinX1][xBinX1] || cutAreasDCe1X2[yBinX2][xBinX2]);
      }
   }
   // DCw
   if (zDC >= 0)
   {
      if (board <= cutAreasDCw0X1Range[0] || board >= cutAreasDCw0X1Range[1] ||
          alpha <= cutAreasDCw0X1Range[2] || alpha >= cutAreasDCw0X1Range[3] ||
          board <= cutAreasDCw0X2Range[0] || board >= cutAreasDCw0X2Range[1] ||
          alpha <= cutAreasDCw0X2Range[2] || alpha >= cutAreasDCw0X2Range[3]) return true;

      const short yBinX1 = static_cast<short>((alpha - cutAreasDCw0X1Range[2])/
                                             (cutAreasDCw0X1Range[3] - cutAreasDCw0X1Range[2])*
                                             static_cast<double>(cutAreasDCw0X1.size()));
      const short xBinX1 = static_cast<short>((board - cutAreasDCw0X1Range[0])/
                                             (cutAreasDCw0X1Range[1] - cutAreasDCw0X1Range[0])*
                                             static_cast<double>(cutAreasDCw0X1[0].size()));
      const short yBinX2 = static_cast<short>((alpha - cutAreasDCw0X2Range[2])/
                                             (cutAreasDCw0X2Range[3] - cutAreasDCw0X2Range[2])*
                                             static_cast<double>(cutAreasDCw0X2.size()));
      const short xBinX2 = static_cast<short>((board - cutAreasDCw0X2Range[0])/
                                             (cutAreasDCw0X2Range[1] - cutAreasDCw0X2Range[0])*
                                             static_cast<double>(cutAreasDCw0X2[0].size()));

      return (cutAreasDCw0X1[yBinX1][xBinX1] || cutAreasDCw0X2[yBinX2][xBinX2]);
   }
   else
   {
      if (board <= cutAreasDCw1X1Range[0] || board >= cutAreasDCw1X1Range[1] ||
          alpha <= cutAreasDCw1X1Range[2] || alpha >= cutAreasDCw1X1Range[3] ||
          board <= cutAreasDCw1X2Range[0] || board >= cutAreasDCw1X2Range[1] ||
          alpha <= cutAreasDCw1X2Range[2] || alpha >= cutAreasDCw1X2Range[3]) return true;

      const short yBinX1 = static_cast<short>((alpha - cutAreasDCw1X1Range[2])/
                                             (cutAreasDCw1X1Range[3] - cutAreasDCw1X1Range[2])*
                                             static_cast<double>(cutAreasDCw1X1.size()));
      const short xBinX1 = static_cast<short>((board - cutAreasDCw1X1Range[0])/
                                             (cutAreasDCw1X1Range[1] - cutAreasDCw1X1Range[0])*
                                             static_cast<double>(cutAreasDCw1X1[0].size()));
      const short yBinX2 = static_cast<short>((alpha - cutAreasDCw1X2Range[2])/
                                             (cutAreasDCw1X2Range[3] - cutAreasDCw1X2Range[2])*
                                             static_cast<double>(cutAreasDCw1X2.size()));
      const short xBinX2 = static_cast<short>((board - cutAreasDCw1X2Range[0])/
                                             (cutAreasDCw1X2Range[1] - cutAreasDCw1X2Range[0])*
                                             static_cast<double>(cutAreasDCw1X2[0].size()));

      return (cutAreasDCw1X1[yBinX1][xBinX1] || cutAreasDCw1X2[yBinX2][xBinX2]);
   }
}

bool DeadMapCutter::IsDeadPC1(const int dcarm, const double ppc1z, const double ppc1phi)
{
   if (!usePC1) return false;
   if (dcarm == 0) // PC1e
   {
      if (ppc1z <= cutAreasPC1eRange[0] || ppc1z >= cutAreasPC1eRange[1] ||
          ppc1phi <= cutAreasPC1eRange[2] || ppc1phi >= cutAreasPC1eRange[3]) return true;

      const short yBin = static_cast<short>((ppc1phi - cutAreasPC1eRange[2])/
                                            (cutAreasPC1eRange[3] - cutAreasPC1eRange[2])*
                                            static_cast<double>(cutAreasPC1e.size()));
      const short xBin = static_cast<short>((ppc1z - cutAreasPC1eRange[0])/
                                            (cutAreasPC1eRange[1] - cutAreasPC1eRange[0])*
                                            static_cast<double>(cutAreasPC1e[0].size()));
      return cutAreasPC1e[yBin][xBin];
   }
   // PC1w
   if (ppc1z <= cutAreasPC1wRange[0] || ppc1z >= cutAreasPC1wRange[1] ||
       ppc1phi <= cutAreasPC1wRange[2] || ppc1phi >= cutAreasPC1wRange[3]) return true;

   const short yBin = static_cast<short>((ppc1phi - cutAreasPC1wRange[2])/
                                         (cutAreasPC1wRange[3] - cutAreasPC1wRange[2])*
                                         static_cast<double>(cutAreasPC1w.size()));
   const short xBin = static_cast<short>((ppc1z - cutAreasPC1wRange[0])/
                                         (cutAreasPC1wRange[1] - cutAreasPC1wRange[0])*
                                         static_cast<double>(cutAreasPC1w[0].size()));
   return cutAreasPC1w[yBin][xBin];
}

bool DeadMapCutter::IsDeadPC2(const double ppc2z, const double ppc2phi)
{
   if (!usePC2) return false;
   if (ppc2z <= cutAreasPC2Range[0] || ppc2z >= cutAreasPC2Range[1] ||
       ppc2phi <= cutAreasPC2Range[2] || ppc2phi >= cutAreasPC2Range[3]) return true;

   const short yBin = static_cast<short>((ppc2phi - cutAreasPC2Range[2])/
                                         (cutAreasPC2Range[3] - cutAreasPC2Range[2])*
                                         static_cast<double>(cutAreasPC2.size()));
   const short xBin = static_cast<short>((ppc2z - cutAreasPC2Range[0])/
                                         (cutAreasPC2Range[1] - cutAreasPC2Range[0])*
                                         static_cast<double>(cutAreasPC2[0].size()));
   return cutAreasPC2[yBin][xBin];
}

bool DeadMapCutter::IsDeadPC3(const int dcarm, const double ppc3z, const double ppc3phi)
{
   if (!usePC3) return false;
   if (dcarm == 0) // PC3e
   {
      if (ppc3z <= cutAreasPC3eRange[0] || ppc3z >= cutAreasPC3eRange[1] ||
          ppc3phi <= cutAreasPC3eRange[2] || ppc3phi >= cutAreasPC3eRange[3]) return true;

      const short yBin = static_cast<short>((ppc3phi - cutAreasPC3eRange[2])/
                                            (cutAreasPC3eRange[3] - cutAreasPC3eRange[2])*
                                            static_cast<double>(cutAreasPC3e.size()));
      const short xBin = static_cast<short>((ppc3z - cutAreasPC3eRange[0])/
                                            (cutAreasPC3eRange[1] - cutAreasPC3eRange[0])*
                                            static_cast<double>(cutAreasPC3e[0].size()));
      return cutAreasPC3e[yBin][xBin];
   }
   // PC3w
   if (ppc3z <= cutAreasPC3wRange[0] || ppc3z >= cutAreasPC3wRange[1] ||
       ppc3phi <= cutAreasPC3wRange[2] || ppc3phi >= cutAreasPC3wRange[3]) return true;

   const short yBin = static_cast<short>((ppc3phi - cutAreasPC3wRange[2])/
                                         (cutAreasPC3wRange[3] - cutAreasPC3wRange[2])*
                                         static_cast<double>(cutAreasPC3w.size()));
   const short xBin = static_cast<short>((ppc3z - cutAreasPC3wRange[0])/
                                         (cutAreasPC3wRange[1] - cutAreasPC3wRange[0])*
                                         static_cast<double>(cutAreasPC3w[0].size()));
   return cutAreasPC3w[yBin][xBin];
}

bool DeadMapCutter::IsDeadTOFe(const int chamber, const int slat)
{
   if (!useTOFe) return false;
   if (chamber < cutAreasTOFeRange[0] || chamber >= cutAreasTOFeRange[1] ||
       slat < cutAreasTOFeRange[2] || slat >= cutAreasTOFeRange[3]) return true;

   return cutAreasTOFe[slat][chamber];
}

bool DeadMapCutter::IsDeadTOFw(const int chamber, const int strip)
{
   if (!useTOFw) return false;

   if (chamber < cutAreasTOFwRange[0] || chamber >= cutAreasTOFwRange[1] ||
       strip < cutAreasTOFwRange[2] || strip >= cutAreasTOFwRange[3]) return true;

   return cutAreasTOFw[strip][chamber];
}

bool DeadMapCutter::IsDeadEMCal(const int dcarm, const int sector, 
                                const int yTower, const int zTower)
{
   if (!useEMCal) return false;
   if (dcarm == 0) // EMCale
   {
      if (yTower < cutAreasEMCaleRange[sector][0] || 
          yTower >= cutAreasEMCaleRange[sector][1] ||
          zTower < cutAreasEMCaleRange[sector][2] || 
          zTower >= cutAreasEMCaleRange[sector][3]) return true;

      return cutAreasEMCale[sector][zTower][yTower];
   }
   // EMCalw
   if (yTower < cutAreasEMCalwRange[sector][0] || 
       yTower >= cutAreasEMCalwRange[sector][1] ||
       zTower < cutAreasEMCalwRange[sector][2] || 
       zTower >= cutAreasEMCalwRange[sector][3]) return true;

   return cutAreasEMCalw[sector][zTower][yTower];
}

bool DeadMapCutter::IsDeadTimingTOFe(const int chamber, const int slat)
{
   if (!useTOFe) return false;
   if (chamber <= cutAreasTimingTOFeRange[0] || chamber >= cutAreasTimingTOFeRange[1] ||
       slat <= cutAreasTimingTOFeRange[2] || slat >= cutAreasTimingTOFeRange[3]) return true;

   const short yBin = static_cast<short>((slat - cutAreasTimingTOFeRange[2])/
                                         (cutAreasTimingTOFeRange[3] - cutAreasTimingTOFeRange[2])*
                                         static_cast<double>(cutAreasTimingTOFe.size()));
   const short xBin = static_cast<short>((chamber - cutAreasTimingTOFeRange[0])/
                                         (cutAreasTimingTOFeRange[1] - cutAreasTimingTOFeRange[0])*
                                         static_cast<double>(cutAreasTimingTOFe[0].size()));
   return cutAreasTimingTOFe[yBin][xBin];
}

bool DeadMapCutter::IsDeadTimingTOFw(const int chamber, const int strip)
{
   if (!useTOFw) return false;

   if (chamber <= cutAreasTimingTOFwRange[0] || chamber >= cutAreasTimingTOFwRange[1] ||
       strip <= cutAreasTimingTOFwRange[2] || strip >= cutAreasTimingTOFwRange[3]) return true;

   const short yBin = static_cast<short>((strip - cutAreasTimingTOFwRange[2])/
                                         (cutAreasTimingTOFwRange[3] - cutAreasTimingTOFwRange[2])*
                                         static_cast<double>(cutAreasTimingTOFw.size()));
   const short xBin = static_cast<short>((chamber - cutAreasTimingTOFwRange[0])/
                                         (cutAreasTimingTOFwRange[1] - cutAreasTimingTOFwRange[0])*
                                         static_cast<double>(cutAreasTimingTOFw[0].size()));
   return cutAreasTimingTOFw[yBin][xBin];
}

bool DeadMapCutter::IsDeadTimingEMCal(const int dcarm, const int sector, 
                                      const int yTower, const int zTower)
{
   if (!useEMCal) return false;
   if (dcarm == 0) // EMCale
   {
      if (yTower <= cutAreasTimingEMCaleRange[sector][0] || 
          yTower >= cutAreasTimingEMCaleRange[sector][1] ||
          zTower <= cutAreasTimingEMCaleRange[sector][2] || 
          zTower >= cutAreasTimingEMCaleRange[sector][3]) return true;

      return cutAreasTimingEMCale[sector][zTower][yTower];
   }
   // EMCalw
   if (yTower <= cutAreasTimingEMCalwRange[sector][0] || 
       yTower >= cutAreasTimingEMCalwRange[sector][1] ||
       zTower <= cutAreasTimingEMCalwRange[sector][2] || 
       zTower >= cutAreasTimingEMCalwRange[sector][3]) return true;

   return cutAreasTimingEMCalw[sector][zTower][yTower];
}

bool DeadMapCutter::SetDeadAreas(const std::string& inputFileName,
                                 std::vector<std::vector<bool>>& cutAreas, double *ranges)
{
   if (!std::filesystem::exists(inputFileName))
   {
      CppTools::PrintWarning("DeadMapCutter: File " + inputFileName + " does not exists");
      return false;
   }

   std::ifstream inputFile(inputFileName);

   int xNBins, yNBins;
   bool isUnexpectedEndOfFile = false;

   if (!(inputFile >> xNBins >> ranges[0] >> ranges[1] >> yNBins >> ranges[2] >> ranges[3]))
   {
      isUnexpectedEndOfFile = true;
   }

   cutAreas.resize(yNBins);

   for (int i = 0; i < yNBins; i++)
   {
      cutAreas[i].resize(xNBins);
      for (int j = 0; j < xNBins; j++)
      {
         bool tmp;
         if (!(inputFile >> tmp)) 
         {
            isUnexpectedEndOfFile = true;
            break;
         }
         cutAreas[i][j] = tmp;
      }
   }

   if (isUnexpectedEndOfFile) 
   {
      CppTools::PrintError("DeadMapCutter: Unexpected end of file: " + inputFileName);
   }

   return true;
}

bool DeadMapCutter::IsDeadDCLoose(const int dcarm, const double zDC, 
                                  const double board, const double alpha,
                                  const double boardVar, const double alphaVar)
{
   if (!useDC) return false;
   if (dcarm == 0) // DCe
   {
      if (zDC >= 0)
      {
         if (board <= cutAreasDCe0X1Range[0] || board >= cutAreasDCe0X1Range[1] ||
             alpha <= cutAreasDCe0X1Range[2] || alpha >= cutAreasDCe0X1Range[3] ||
             board <= cutAreasDCe0X2Range[0] || board >= cutAreasDCe0X2Range[1] ||
             alpha <= cutAreasDCe0X2Range[2] || alpha >= cutAreasDCe0X2Range[3]) return true;

         short yBinX1Min = static_cast<short>((alpha - alphaVar - cutAreasDCe0X1Range[2])/
                                              (cutAreasDCe0X1Range[3] - cutAreasDCe0X1Range[2])*
                                              static_cast<double>(cutAreasDCe0X1.size()));
         short xBinX1Min = static_cast<short>((board - boardVar - cutAreasDCe0X1Range[0])/
                                              (cutAreasDCe0X1Range[1] - cutAreasDCe0X1Range[0])*
                                              static_cast<double>(cutAreasDCe0X1[0].size()));
         short yBinX2Min = static_cast<short>((alpha - alphaVar - cutAreasDCe0X2Range[2])/
                                              (cutAreasDCe0X2Range[3] - cutAreasDCe0X2Range[2])*
                                              static_cast<double>(cutAreasDCe0X2.size()));
         short xBinX2Min = static_cast<short>((board - boardVar - cutAreasDCe0X2Range[0])/
                                              (cutAreasDCe0X2Range[1] - cutAreasDCe0X2Range[0])*
                                              static_cast<double>(cutAreasDCe0X2[0].size()));

         short yBinX1Max = static_cast<short>((alpha + alphaVar - cutAreasDCe0X1Range[2])/
                                              (cutAreasDCe0X1Range[3] - cutAreasDCe0X1Range[2])*
                                              static_cast<double>(cutAreasDCe0X1.size()));
         short xBinX1Max = static_cast<short>((board + boardVar - cutAreasDCe0X1Range[0])/
                                              (cutAreasDCe0X1Range[1] - cutAreasDCe0X1Range[0])*
                                              static_cast<double>(cutAreasDCe0X1[0].size()));
         short yBinX2Max = static_cast<short>((alpha + alphaVar - cutAreasDCe0X2Range[2])/
                                              (cutAreasDCe0X2Range[3] - cutAreasDCe0X2Range[2])*
                                              static_cast<double>(cutAreasDCe0X2.size()));
         short xBinX2Max = static_cast<short>((board + boardVar - cutAreasDCe0X2Range[0])/
                                              (cutAreasDCe0X2Range[1] - cutAreasDCe0X2Range[0])*
                                              static_cast<double>(cutAreasDCe0X2[0].size()));

         if (yBinX1Min < 0) yBinX1Min = 0;
         if (xBinX1Min < 0) xBinX1Min = 0;
         if (yBinX2Min < 0) yBinX2Min = 0;
         if (xBinX2Min < 0) xBinX2Min = 0;

         if (yBinX1Max >= static_cast<short>(cutAreasDCe0X1.size())) 
            yBinX1Max = static_cast<short>(cutAreasDCe0X1.size()) - 1;
         if (xBinX1Max >= static_cast<short>(cutAreasDCe0X1[0].size())) 
            xBinX1Max = static_cast<short>(cutAreasDCe0X1[0].size()) - 1;
         if (yBinX2Max >= static_cast<short>(cutAreasDCe0X2.size())) 
            yBinX2Max = static_cast<short>(cutAreasDCe0X2.size()) - 1;
         if (xBinX2Max >= static_cast<short>(cutAreasDCe0X2[0].size())) 
            xBinX2Max = static_cast<short>(cutAreasDCe0X2[0].size()) - 1;

         bool isDeadDCX1 = true;
         for (short xBin = xBinX1Min; xBin <= xBinX1Max; xBin++)
         {
            for (short yBin = yBinX1Min; yBin <= yBinX1Max; yBin++)
            {
               if (!cutAreasDCe0X1[yBin][xBin]) 
               {
                  isDeadDCX1 = false;
                  break;
               }
            }
         }

         bool isDeadDCX2 = true;
         for (short xBin = xBinX2Min; xBin <= xBinX2Max; xBin++)
         {
            for (short yBin = yBinX2Min; yBin <= yBinX2Max; yBin++)
            {
               if (!cutAreasDCe0X2[yBin][xBin])
               {
                  isDeadDCX2 = false;
                  break;
               }
            }
         }
         return (isDeadDCX1 || isDeadDCX2);
      }
      else
      {
         if (board <= cutAreasDCe1X1Range[0] || board >= cutAreasDCe1X1Range[1] ||
             alpha <= cutAreasDCe1X1Range[2] || alpha >= cutAreasDCe1X1Range[3] ||
             board <= cutAreasDCe1X2Range[0] || board >= cutAreasDCe1X2Range[1] ||
             alpha <= cutAreasDCe1X2Range[2] || alpha >= cutAreasDCe1X2Range[3]) return true;

         short yBinX1Min = static_cast<short>((alpha - alphaVar - cutAreasDCe1X1Range[2])/
                                              (cutAreasDCe1X1Range[3] - cutAreasDCe1X1Range[2])*
                                              static_cast<double>(cutAreasDCe1X1.size()));
         short xBinX1Min = static_cast<short>((board - boardVar - cutAreasDCe1X1Range[0])/
                                              (cutAreasDCe1X1Range[1] - cutAreasDCe1X1Range[0])*
                                              static_cast<double>(cutAreasDCe1X1[0].size()));
         short yBinX2Min = static_cast<short>((alpha - alphaVar - cutAreasDCe1X2Range[2])/
                                              (cutAreasDCe1X2Range[3] - cutAreasDCe1X2Range[2])*
                                              static_cast<double>(cutAreasDCe1X2.size()));
         short xBinX2Min = static_cast<short>((board - boardVar - cutAreasDCe1X2Range[0])/
                                              (cutAreasDCe1X2Range[1] - cutAreasDCe1X2Range[0])*
                                              static_cast<double>(cutAreasDCe1X2[0].size()));

         short yBinX1Max = static_cast<short>((alpha + alphaVar - cutAreasDCe1X1Range[2])/
                                              (cutAreasDCe1X1Range[3] - cutAreasDCe1X1Range[2])*
                                              static_cast<double>(cutAreasDCe1X1.size()));
         short xBinX1Max = static_cast<short>((board + boardVar - cutAreasDCe1X1Range[0])/
                                              (cutAreasDCe1X1Range[1] - cutAreasDCe1X1Range[0])*
                                              static_cast<double>(cutAreasDCe1X1[0].size()));
         short yBinX2Max = static_cast<short>((alpha + alphaVar - cutAreasDCe1X2Range[2])/
                                              (cutAreasDCe1X2Range[3] - cutAreasDCe1X2Range[2])*
                                              static_cast<double>(cutAreasDCe1X2.size()));
         short xBinX2Max = static_cast<short>((board + boardVar - cutAreasDCe1X2Range[0])/
                                              (cutAreasDCe1X2Range[1] - cutAreasDCe1X2Range[0])*
                                              static_cast<double>(cutAreasDCe1X2[0].size()));

         if (yBinX1Min < 0) yBinX1Min = 0;
         if (xBinX1Min < 0) xBinX1Min = 0;
         if (yBinX2Min < 0) yBinX2Min = 0;
         if (xBinX2Min < 0) xBinX2Min = 0;

         if (yBinX1Max >= static_cast<short>(cutAreasDCe1X1.size())) 
            yBinX1Max = static_cast<short>(cutAreasDCe1X1.size()) - 1;
         if (xBinX1Max >= static_cast<short>(cutAreasDCe1X1[0].size())) 
            xBinX1Max = static_cast<short>(cutAreasDCe1X1[0].size()) - 1;
         if (yBinX2Max >= static_cast<short>(cutAreasDCe1X2.size())) 
            yBinX2Max = static_cast<short>(cutAreasDCe1X2.size()) - 1;
         if (xBinX2Max >= static_cast<short>(cutAreasDCe1X2[0].size())) 
            xBinX2Max = static_cast<short>(cutAreasDCe1X2[0].size()) - 1;

         bool isDeadDCX1 = true;
         for (short xBin = xBinX1Min; xBin <= xBinX1Max; xBin++)
         {
            for (short yBin = yBinX1Min; yBin <= yBinX1Max; yBin++)
            {
               if(!cutAreasDCe1X1[yBin][xBin])
               {
                  isDeadDCX1 = false;
                  break;
               }
            }
         }

         bool isDeadDCX2 = true;
         for (short xBin = xBinX2Min; xBin <= xBinX2Max; xBin++)
         {
            for (short yBin = yBinX2Min; yBin <= yBinX2Max; yBin++)
            {
               if (!cutAreasDCe1X2[yBin][xBin])
               {
                  isDeadDCX2 = false;
                  break;
               }
            }
         }
         return (isDeadDCX1 || isDeadDCX2);
      }
   }
   // DCw
   if (zDC >= 0)
   {
      if (board <= cutAreasDCw0X1Range[0] || board >= cutAreasDCw0X1Range[1] ||
          alpha <= cutAreasDCw0X1Range[2] || alpha >= cutAreasDCw0X1Range[3] ||
          board <= cutAreasDCw0X2Range[0] || board >= cutAreasDCw0X2Range[1] ||
          alpha <= cutAreasDCw0X2Range[2] || alpha >= cutAreasDCw0X2Range[3]) return true;

      short yBinX1Min = static_cast<short>((alpha - alphaVar - cutAreasDCw0X1Range[2])/
                                           (cutAreasDCw0X1Range[3] - cutAreasDCw0X1Range[2])*
                                           static_cast<double>(cutAreasDCw0X1.size()));
      short xBinX1Min = static_cast<short>((board - boardVar - cutAreasDCw0X1Range[0])/
                                           (cutAreasDCw0X1Range[1] - cutAreasDCw0X1Range[0])*
                                           static_cast<double>(cutAreasDCw0X1[0].size()));
      short yBinX2Min = static_cast<short>((alpha - alphaVar - cutAreasDCw0X2Range[2])/
                                           (cutAreasDCw0X2Range[3] - cutAreasDCw0X2Range[2])*
                                           static_cast<double>(cutAreasDCw0X2.size()));
      short xBinX2Min = static_cast<short>((board - boardVar - cutAreasDCw0X2Range[0])/
                                           (cutAreasDCw0X2Range[1] - cutAreasDCw0X2Range[0])*
                                           static_cast<double>(cutAreasDCw0X2[0].size()));

      short yBinX1Max = static_cast<short>((alpha + alphaVar - cutAreasDCw0X1Range[2])/
                                           (cutAreasDCw0X1Range[3] - cutAreasDCw0X1Range[2])*
                                           static_cast<double>(cutAreasDCw0X1.size()));
      short xBinX1Max = static_cast<short>((board + boardVar - cutAreasDCw0X1Range[0])/
                                           (cutAreasDCw0X1Range[1] - cutAreasDCw0X1Range[0])*
                                           static_cast<double>(cutAreasDCw0X1[0].size()));
      short yBinX2Max = static_cast<short>((alpha + alphaVar - cutAreasDCw0X2Range[2])/
                                           (cutAreasDCw0X2Range[3] - cutAreasDCw0X2Range[2])*
                                           static_cast<double>(cutAreasDCw0X2.size()));
      short xBinX2Max = static_cast<short>((board + boardVar - cutAreasDCw0X2Range[0])/
                                           (cutAreasDCw0X2Range[1] - cutAreasDCw0X2Range[0])*
                                           static_cast<double>(cutAreasDCw0X2[0].size()));

      if (yBinX1Min < 0) yBinX1Min = 0;
      if (xBinX1Min < 0) xBinX1Min = 0;
      if (yBinX2Min < 0) yBinX2Min = 0;
      if (xBinX2Min < 0) xBinX2Min = 0;

      if (yBinX1Max >= static_cast<short>(cutAreasDCw0X1.size())) 
         yBinX1Max = static_cast<short>(cutAreasDCw0X1.size()) - 1;
      if (xBinX1Max >= static_cast<short>(cutAreasDCw0X1[0].size())) 
         xBinX1Max = static_cast<short>(cutAreasDCw0X1[0].size()) - 1;
      if (yBinX2Max >= static_cast<short>(cutAreasDCw0X2.size())) 
         yBinX2Max = static_cast<short>(cutAreasDCw0X2.size()) - 1;
      if (xBinX2Max >= static_cast<short>(cutAreasDCw0X2[0].size())) 
         xBinX2Max = static_cast<short>(cutAreasDCw0X2[0].size()) - 1;

      bool isDeadDCX1 = true;
      for (short xBin = xBinX1Min; xBin <= xBinX1Max; xBin++)
      {
         for (short yBin = yBinX1Min; yBin <= yBinX1Max; yBin++)
         {
            if (!cutAreasDCw0X1[yBin][xBin])
            {
               isDeadDCX1 = false;
               break;
            }
         }
      }

      bool isDeadDCX2 = true;
      for (short xBin = xBinX2Min; xBin <= xBinX2Max; xBin++)
      {
         for (short yBin = yBinX2Min; yBin <= yBinX2Max; yBin++)
         {
            if (!cutAreasDCw0X2[yBin][xBin])
            {
               isDeadDCX2 = false;
               break;
            }
         }
      }
      return (isDeadDCX1 || isDeadDCX2);
   }
   else
   {
      if (board <= cutAreasDCw1X1Range[0] || board >= cutAreasDCw1X1Range[1] ||
          alpha <= cutAreasDCw1X1Range[2] || alpha >= cutAreasDCw1X1Range[3] ||
          board <= cutAreasDCw1X2Range[0] || board >= cutAreasDCw1X2Range[1] ||
          alpha <= cutAreasDCw1X2Range[2] || alpha >= cutAreasDCw1X2Range[3]) return true;

      short yBinX1Min = static_cast<short>((alpha - alphaVar - cutAreasDCw1X1Range[2])/
                                           (cutAreasDCw1X1Range[3] - cutAreasDCw1X1Range[2])*
                                           static_cast<double>(cutAreasDCw1X1.size()));
      short xBinX1Min = static_cast<short>((board - boardVar - cutAreasDCw1X1Range[0])/
                                           (cutAreasDCw1X1Range[1] - cutAreasDCw1X1Range[0])*
                                           static_cast<double>(cutAreasDCw1X1[0].size()));
      short yBinX2Min = static_cast<short>((alpha - alphaVar - cutAreasDCw1X2Range[2])/
                                           (cutAreasDCw1X2Range[3] - cutAreasDCw1X2Range[2])*
                                           static_cast<double>(cutAreasDCw1X2.size()));
      short xBinX2Min = static_cast<short>((board - boardVar - cutAreasDCw1X2Range[0])/
                                           (cutAreasDCw1X2Range[1] - cutAreasDCw1X2Range[0])*
                                           static_cast<double>(cutAreasDCw1X2[0].size()));

      short yBinX1Max = static_cast<short>((alpha + alphaVar - cutAreasDCw1X1Range[2])/
                                           (cutAreasDCw1X1Range[3] - cutAreasDCw1X1Range[2])*
                                           static_cast<double>(cutAreasDCw1X1.size()));
      short xBinX1Max = static_cast<short>((board + boardVar - cutAreasDCw1X1Range[0])/
                                           (cutAreasDCw1X1Range[1] - cutAreasDCw1X1Range[0])*
                                           static_cast<double>(cutAreasDCw1X1[0].size()));
      short yBinX2Max = static_cast<short>((alpha + alphaVar - cutAreasDCw1X2Range[2])/
                                           (cutAreasDCw1X2Range[3] - cutAreasDCw1X2Range[2])*
                                           static_cast<double>(cutAreasDCw1X2.size()));
      short xBinX2Max = static_cast<short>((board + boardVar - cutAreasDCw1X2Range[0])/
                                           (cutAreasDCw1X2Range[1] - cutAreasDCw1X2Range[0])*
                                           static_cast<double>(cutAreasDCw1X2[0].size()));

      if (yBinX1Min < 0) yBinX1Min = 0;
      if (xBinX1Min < 0) xBinX1Min = 0;
      if (yBinX2Min < 0) yBinX2Min = 0;
      if (xBinX2Min < 0) xBinX2Min = 0;

      if (yBinX1Max >= static_cast<short>(cutAreasDCw1X1.size())) 
         yBinX1Max = static_cast<short>(cutAreasDCw1X1.size()) - 1;
      if (xBinX1Max >= static_cast<short>(cutAreasDCw1X1[0].size())) 
         xBinX1Max = static_cast<short>(cutAreasDCw1X1[0].size()) - 1;
      if (yBinX2Max >= static_cast<short>(cutAreasDCw1X2.size())) 
         yBinX2Max = static_cast<short>(cutAreasDCw1X2.size()) - 1;
      if (xBinX2Max >= static_cast<short>(cutAreasDCw1X2[0].size())) 
         xBinX2Max = static_cast<short>(cutAreasDCw1X2[0].size()) - 1;

      bool isDeadDCX1 = true;
      for (short xBin = xBinX1Min; xBin <= xBinX1Max; xBin++)
      {
         for (short yBin = yBinX1Min; yBin <= yBinX1Max; yBin++)
         {
            if (!cutAreasDCw1X1[yBin][xBin])
            {
               isDeadDCX1 = false;
               break;
            }
         }
      }

      bool isDeadDCX2 = true;
      for (short xBin = xBinX2Min; xBin <= xBinX2Max; xBin++)
      {
         for (short yBin = yBinX2Min; yBin <= yBinX2Max; yBin++)
         {
            if (!cutAreasDCw1X2[yBin][xBin])
            {
               isDeadDCX2 = false;
               break;
            }
         }
      }
      return (isDeadDCX1 || isDeadDCX2);
   }
}

bool DeadMapCutter::IsDeadPC1Loose(const int dcarm, const double ppc1z, const double ppc1phi,
                                   const double ppc1zVar, const double ppc1phiVar)
{
   if (!usePC1) return false;
   if (dcarm == 0) // PC1e
   {
      if (ppc1z <= cutAreasPC1eRange[0] || ppc1z >= cutAreasPC1eRange[1] ||
          ppc1phi <= cutAreasPC1eRange[2] || ppc1phi >= cutAreasPC1eRange[3]) return true;

      short yBinMin = static_cast<short>((ppc1phi - ppc1phiVar - cutAreasPC1eRange[2])/
                                         (cutAreasPC1eRange[3] - cutAreasPC1eRange[2])*
                                         static_cast<double>(cutAreasPC1e.size()));
      short xBinMin = static_cast<short>((ppc1z - ppc1zVar - cutAreasPC1eRange[0])/
                                         (cutAreasPC1eRange[1] - cutAreasPC1eRange[0])*
                                         static_cast<double>(cutAreasPC1e[0].size()));
      short yBinMax = static_cast<short>((ppc1phi + ppc1phiVar - cutAreasPC1eRange[2])/
                                         (cutAreasPC1eRange[3] - cutAreasPC1eRange[2])*
                                         static_cast<double>(cutAreasPC1e.size()));
      short xBinMax = static_cast<short>((ppc1z + ppc1zVar - cutAreasPC1eRange[0])/
                                         (cutAreasPC1eRange[1] - cutAreasPC1eRange[0])*
                                         static_cast<double>(cutAreasPC1e[0].size()));

      if (yBinMin < 0) yBinMin = 0;
      if (xBinMin < 0) xBinMin = 0;
      if (yBinMax >= static_cast<short>(cutAreasPC1e.size())) 
         yBinMax = static_cast<short>(cutAreasPC1e.size()) - 1;
      if (xBinMax >= static_cast<short>(cutAreasPC1e[0].size())) 
         xBinMax = static_cast<short>(cutAreasPC1e[0].size()) - 1;

      for (short xBin = xBinMin; xBin <= xBinMax; xBin++)
      {
         for (short yBin = yBinMin; yBin <= yBinMax; yBin++)
         {
            if (!cutAreasPC1e[yBin][xBin]) return false;
         }
      }
      return true;
   }
   // PC1w
   if (ppc1z <= cutAreasPC1wRange[0] || ppc1z >= cutAreasPC1wRange[1] ||
       ppc1phi <= cutAreasPC1wRange[2] || ppc1phi >= cutAreasPC1wRange[3]) return true;

   short yBinMin = static_cast<short>((ppc1phi - ppc1phiVar - cutAreasPC1wRange[2])/
                                      (cutAreasPC1wRange[3] - cutAreasPC1wRange[2])*
                                      static_cast<double>(cutAreasPC1w.size()));
   short xBinMin = static_cast<short>((ppc1z - ppc1zVar - cutAreasPC1wRange[0])/
                                      (cutAreasPC1wRange[1] - cutAreasPC1wRange[0])*
                                      static_cast<double>(cutAreasPC1w[0].size()));
   short yBinMax = static_cast<short>((ppc1phi + ppc1phiVar - cutAreasPC1wRange[2])/
                                      (cutAreasPC1wRange[3] - cutAreasPC1wRange[2])*
                                      static_cast<double>(cutAreasPC1w.size()));
   short xBinMax = static_cast<short>((ppc1z + ppc1zVar - cutAreasPC1wRange[0])/
                                      (cutAreasPC1wRange[1] - cutAreasPC1wRange[0])*
                                      static_cast<double>(cutAreasPC1w[0].size()));

   if (yBinMin < 0) yBinMin = 0;
   if (xBinMin < 0) xBinMin = 0;
   if (yBinMax >= static_cast<short>(cutAreasPC1w.size())) 
      yBinMax = static_cast<short>(cutAreasPC1w.size()) - 1;
   if (xBinMax >= static_cast<short>(cutAreasPC1w[0].size())) 
      xBinMax = static_cast<short>(cutAreasPC1w[0].size()) - 1;

   for (short xBin = xBinMin; xBin <= xBinMax; xBin++)
   {
      for (short yBin = yBinMin; yBin <= yBinMax; yBin++)
      {
         if (!cutAreasPC1w[yBin][xBin]) return false;
      }
   }
   return true;
}

bool DeadMapCutter::IsDeadPC2Loose(const double ppc2z, const double ppc2phi,
                                   const double ppc2zVar, const double ppc2phiVar)
{
   if (!usePC2) return false;
   if (ppc2z <= cutAreasPC2Range[0] || ppc2z >= cutAreasPC2Range[1] ||
       ppc2phi <= cutAreasPC2Range[2] || ppc2phi >= cutAreasPC2Range[3]) return true;

   short yBinMin = static_cast<short>((ppc2phi - ppc2phiVar - cutAreasPC2Range[2])/
                                      (cutAreasPC2Range[3] - cutAreasPC2Range[2])*
                                      static_cast<double>(cutAreasPC2.size()));
   short xBinMin = static_cast<short>((ppc2z - ppc2zVar - cutAreasPC2Range[0])/
                                      (cutAreasPC2Range[1] - cutAreasPC2Range[0])*
                                      static_cast<double>(cutAreasPC2[0].size()));
   short yBinMax = static_cast<short>((ppc2phi + ppc2phiVar - cutAreasPC2Range[2])/
                                      (cutAreasPC2Range[3] - cutAreasPC2Range[2])*
                                      static_cast<double>(cutAreasPC2.size()));
   short xBinMax = static_cast<short>((ppc2z + ppc2zVar - cutAreasPC2Range[0])/
                                      (cutAreasPC2Range[1] - cutAreasPC2Range[0])*
                                      static_cast<double>(cutAreasPC2[0].size()));

   if (yBinMin < 0) yBinMin = 0;
   if (xBinMin < 0) xBinMin = 0;
   if (yBinMax >= static_cast<short>(cutAreasPC2.size())) 
      yBinMax = static_cast<short>(cutAreasPC2.size()) - 1;
   if (xBinMax >= static_cast<short>(cutAreasPC2[0].size())) 
      xBinMax = static_cast<short>(cutAreasPC2[0].size()) - 1;

   for (short xBin = xBinMin; xBin <= xBinMax; xBin++)
   {
      for (short yBin = yBinMin; yBin <= yBinMax; yBin++)
      {
         if (!cutAreasPC2[yBin][xBin]) return false;
      }
   }
   return true;
}

bool DeadMapCutter::IsDeadPC3Loose(const int dcarm, const double ppc3z, const double ppc3phi,
                                   const double ppc3zVar, const double ppc3phiVar)
{
   if (!usePC3) return false;
   if (dcarm == 0) // PC3e
   {
      if (ppc3z <= cutAreasPC3eRange[0] || ppc3z >= cutAreasPC3eRange[1] ||
          ppc3phi <= cutAreasPC3eRange[2] || ppc3phi >= cutAreasPC3eRange[3]) return true;

      short yBinMin = static_cast<short>((ppc3phi - ppc3phiVar - cutAreasPC3eRange[2])/
                                         (cutAreasPC3eRange[3] - cutAreasPC3eRange[2])*
                                         static_cast<double>(cutAreasPC3e.size()));
      short xBinMin = static_cast<short>((ppc3z - ppc3zVar - cutAreasPC3eRange[0])/
                                         (cutAreasPC3eRange[1] - cutAreasPC3eRange[0])*
                                         static_cast<double>(cutAreasPC3e[0].size()));
      short yBinMax = static_cast<short>((ppc3phi + ppc3phiVar - cutAreasPC3eRange[2])/
                                         (cutAreasPC3eRange[3] - cutAreasPC3eRange[2])*
                                         static_cast<double>(cutAreasPC3e.size()));
      short xBinMax = static_cast<short>((ppc3z + ppc3zVar - cutAreasPC3eRange[0])/
                                         (cutAreasPC3eRange[1] - cutAreasPC3eRange[0])*
                                         static_cast<double>(cutAreasPC3e[0].size()));

      if (yBinMin < 0) yBinMin = 0;
      if (xBinMin < 0) xBinMin = 0;
      if (yBinMax >= static_cast<short>(cutAreasPC3e.size())) 
         yBinMax = static_cast<short>(cutAreasPC3e.size()) - 1;
      if (xBinMax >= static_cast<short>(cutAreasPC3e[0].size())) 
         xBinMax = static_cast<short>(cutAreasPC3e[0].size()) - 1;

      for (short xBin = xBinMin; xBin <= xBinMax; xBin++)
      {
         for (short yBin = yBinMin; yBin <= yBinMax; yBin++)
         {
            if (!cutAreasPC3e[yBin][xBin]) return false;
         }
      }
      return true;
   }
   // PC3w
   if (ppc3z <= cutAreasPC3wRange[0] || ppc3z >= cutAreasPC3wRange[1] ||
       ppc3phi <= cutAreasPC3wRange[2] || ppc3phi >= cutAreasPC3wRange[3]) return true;

   short yBinMin = static_cast<short>((ppc3phi - ppc3phiVar - cutAreasPC3wRange[2])/
                                      (cutAreasPC3wRange[3] - cutAreasPC3wRange[2])*
                                      static_cast<double>(cutAreasPC3w.size()));
   short xBinMin = static_cast<short>((ppc3z - ppc3zVar - cutAreasPC3wRange[0])/
                                      (cutAreasPC3wRange[1] - cutAreasPC3wRange[0])*
                                      static_cast<double>(cutAreasPC3w[0].size()));
   short yBinMax = static_cast<short>((ppc3phi + ppc3phiVar - cutAreasPC3wRange[2])/
                                      (cutAreasPC3wRange[3] - cutAreasPC3wRange[2])*
                                      static_cast<double>(cutAreasPC3w.size()));
   short xBinMax = static_cast<short>((ppc3z + ppc3zVar - cutAreasPC3wRange[0])/
                                      (cutAreasPC3wRange[1] - cutAreasPC3wRange[0])*
                                      static_cast<double>(cutAreasPC3w[0].size()));

   if (yBinMin < 0) yBinMin = 0;
   if (xBinMin < 0) xBinMin = 0;
   if (yBinMax >= static_cast<short>(cutAreasPC3w.size())) 
      yBinMax = static_cast<short>(cutAreasPC3w.size()) - 1;
   if (xBinMax >= static_cast<short>(cutAreasPC3w[0].size())) 
      xBinMax = static_cast<short>(cutAreasPC3w[0].size()) - 1;

   for (short xBin = xBinMin; xBin <= xBinMax; xBin++)
   {
      for (short yBin = yBinMin; yBin <= yBinMax; yBin++)
      {
         if (!cutAreasPC3w[yBin][xBin]) return false;
      }
   }
   return true;
}

bool DeadMapCutter::IsDeadTOFeLoose(const int chamber, const int slat, 
                                    const int chamberVar, const int slatVar)
{
   if (!useTOFe) return false;
   if (chamber < cutAreasTOFeRange[0] || chamber >= cutAreasTOFeRange[1] ||
       slat < cutAreasTOFeRange[2] || slat >= cutAreasTOFeRange[3]) return true;

   short yBinMin = slat - slatVar;
   short xBinMin = chamber - chamberVar;
   short yBinMax = slat + slatVar;
   short xBinMax = chamber + chamberVar;

   if (yBinMin < cutAreasTOFeRange[2]) yBinMin = cutAreasTOFeRange[2];
   if (xBinMin < cutAreasTOFeRange[0]) xBinMin = cutAreasTOFeRange[0];
   if (yBinMax >= cutAreasTOFeRange[1]) yBinMax = cutAreasTOFeRange[3] - 1;
   if (xBinMax >= cutAreasTOFeRange[3]) xBinMax = cutAreasTOFeRange[1] - 1;

   for (short xBin = xBinMin; xBin <= xBinMax; xBin++)
   {
      for (short yBin = yBinMin; yBin <= yBinMax; yBin++)
      {
         if (!cutAreasTOFe[yBin][xBin]) return false;
      }
   }
   return true;
}

bool DeadMapCutter::IsDeadTOFwLoose(const int chamber, const int strip,
                                    const int chamberVar, const int stripVar)
{
   if (!useTOFw) return false;

   if (chamber < cutAreasTOFwRange[0] || chamber >= cutAreasTOFwRange[1] ||
       strip < cutAreasTOFwRange[2] || strip >= cutAreasTOFwRange[3]) return true;

   short yBinMin = strip - stripVar;
   short xBinMin = chamber - chamberVar;
   short yBinMax = strip + stripVar;
   short xBinMax = chamber + chamberVar;

   if (yBinMin < cutAreasTOFwRange[2]) yBinMin = cutAreasTOFwRange[2];
   if (xBinMin < cutAreasTOFwRange[0]) xBinMin = cutAreasTOFwRange[0];
   if (yBinMax >= cutAreasTOFwRange[1]) yBinMax = cutAreasTOFwRange[3] - 1;
   if (xBinMax >= cutAreasTOFwRange[3]) xBinMax = cutAreasTOFwRange[1] - 1;

   for (short xBin = xBinMin; xBin <= xBinMax; xBin++)
   {
      for (short yBin = yBinMin; yBin <= yBinMax; yBin++)
      {
         if (!cutAreasTOFw[yBin][xBin]) return false;
      }
   }
   return true;
}

bool DeadMapCutter::IsDeadEMCalLoose(const int dcarm, const int sector, 
                                     const int yTower, const int zTower,
                                     const int yTowerVar, const int zTowerVar)
{
   if (!useEMCal) return false;

   if (dcarm == 0) // EMCale
   {
      if (yTower < cutAreasEMCaleRange[sector][0] || 
          yTower >= cutAreasEMCaleRange[sector][1] ||
          zTower < cutAreasEMCaleRange[sector][2] || 
          zTower >= cutAreasEMCaleRange[sector][3]) return true;

      short yBinMin = zTower - zTowerVar;
      short xBinMin = yTower - yTowerVar;
      short yBinMax = zTower + zTowerVar;
      short xBinMax = yTower + yTowerVar;

      if (yBinMin < cutAreasEMCaleRange[sector][2]) yBinMin = cutAreasEMCaleRange[sector][2];
      if (xBinMin < cutAreasEMCaleRange[sector][0]) xBinMin = cutAreasEMCaleRange[sector][0];
      if (yBinMax >= cutAreasEMCaleRange[sector][1]) yBinMax = cutAreasEMCaleRange[sector][3] - 1;
      if (xBinMax >= cutAreasEMCaleRange[sector][3]) xBinMax = cutAreasEMCaleRange[sector][1] - 1;

      for (short xBin = xBinMin; xBin <= xBinMax; xBin++)
      {
         for (short yBin = yBinMin; yBin <= yBinMax; yBin++)
         {
            if (!cutAreasEMCale[sector][yBin][xBin]) return false;
         }
      }
      return true;
   }
   // EMCalw
   if (yTower < cutAreasEMCalwRange[sector][0] || 
       yTower >= cutAreasEMCalwRange[sector][1] ||
       zTower < cutAreasEMCalwRange[sector][2] || 
       zTower >= cutAreasEMCalwRange[sector][3]) return true;

   short yBinMin = zTower - zTowerVar;
   short xBinMin = yTower - yTowerVar;
   short yBinMax = zTower + zTowerVar;
   short xBinMax = yTower + yTowerVar;

   if (yBinMin < cutAreasEMCalwRange[sector][2]) yBinMin = cutAreasEMCalwRange[sector][2];
   if (xBinMin < cutAreasEMCalwRange[sector][0]) xBinMin = cutAreasEMCalwRange[sector][0];
   if (yBinMax >= cutAreasEMCalwRange[sector][1]) yBinMax = cutAreasEMCalwRange[sector][3] - 1;
   if (xBinMax >= cutAreasEMCalwRange[sector][3]) xBinMax = cutAreasEMCalwRange[sector][1] - 1;

   for (short xBin = xBinMin; xBin <= xBinMax; xBin++)
   {
      for (short yBin = yBinMin; yBin <= yBinMax; yBin++)
      {
         if (!cutAreasEMCalw[sector][yBin][xBin]) return false;
      }
   }
   return true;
}

bool DeadMapCutter::IsDeadTimingTOFeLoose(const int chamber, const int slat, 
                                          const int chamberVar, const int slatVar)
{
   if (!useTOFe) return false;
   if (chamber < cutAreasTimingTOFeRange[0] || chamber >= cutAreasTimingTOFeRange[1] ||
       slat < cutAreasTimingTOFeRange[2] || slat >= cutAreasTimingTOFeRange[3]) return true;

   short yBinMin = slat - slatVar;
   short xBinMin = chamber - chamberVar;
   short yBinMax = slat + slatVar;
   short xBinMax = chamber + chamberVar;

   if (yBinMin < cutAreasTimingTOFeRange[2]) yBinMin = cutAreasTimingTOFeRange[2];
   if (xBinMin < cutAreasTimingTOFeRange[0]) xBinMin = cutAreasTimingTOFeRange[0];
   if (yBinMax >= cutAreasTimingTOFeRange[1]) yBinMax = cutAreasTimingTOFeRange[3] - 1;
   if (xBinMax >= cutAreasTimingTOFeRange[3]) xBinMax = cutAreasTimingTOFeRange[1] - 1;

   for (short xBin = xBinMin; xBin <= xBinMax; xBin++)
   {
      for (short yBin = yBinMin; yBin <= yBinMax; yBin++)
      {
         if (!cutAreasTimingTOFe[yBin][xBin]) return false;
      }
   }
   return true;
}

bool DeadMapCutter::IsDeadTimingTOFwLoose(const int chamber, const int strip,
                                          const int chamberVar, const int stripVar)
{
   if (!useTOFw) return false;

   if (chamber < cutAreasTimingTOFwRange[0] || chamber >= cutAreasTimingTOFwRange[1] ||
       strip < cutAreasTimingTOFwRange[2] || strip >= cutAreasTimingTOFwRange[3]) return true;

   short yBinMin = strip - stripVar;
   short xBinMin = chamber - chamberVar;
   short yBinMax = strip + stripVar;
   short xBinMax = chamber + chamberVar;

   if (yBinMin < cutAreasTimingTOFwRange[2]) yBinMin = cutAreasTimingTOFwRange[2];
   if (xBinMin < cutAreasTimingTOFwRange[0]) xBinMin = cutAreasTimingTOFwRange[0];
   if (yBinMax >= cutAreasTimingTOFwRange[1]) yBinMax = cutAreasTimingTOFwRange[3] - 1;
   if (xBinMax >= cutAreasTimingTOFwRange[3]) xBinMax = cutAreasTimingTOFwRange[1] - 1;

   for (short xBin = xBinMin; xBin <= xBinMax; xBin++)
   {
      for (short yBin = yBinMin; yBin <= yBinMax; yBin++)
      {
         if (!cutAreasTimingTOFw[yBin][xBin]) return false;
      }
   }
   return true;
}

bool DeadMapCutter::IsDeadTimingEMCalLoose(const int dcarm, const int sector, 
                                           const int yTower, const int zTower,
                                           const int yTowerVar, const int zTowerVar)
{
   if (!useEMCal) return false;

   if (dcarm == 0) // EMCale
   {
      if (yTower < cutAreasTimingEMCaleRange[sector][0] || 
          yTower >= cutAreasTimingEMCaleRange[sector][1] ||
          zTower < cutAreasTimingEMCaleRange[sector][2] || 
          zTower >= cutAreasTimingEMCaleRange[sector][3]) return true;

      short yBinMin = zTower - zTowerVar;
      short xBinMin = yTower - yTowerVar;
      short yBinMax = zTower + zTowerVar;
      short xBinMax = yTower + yTowerVar;

      if (yBinMin < cutAreasTimingEMCaleRange[sector][2]) 
         yBinMin = cutAreasTimingEMCaleRange[sector][2];
      if (xBinMin < cutAreasTimingEMCaleRange[sector][0]) 
         xBinMin = cutAreasTimingEMCaleRange[sector][0];
      if (yBinMax >= cutAreasTimingEMCaleRange[sector][1]) 
         yBinMax = cutAreasTimingEMCaleRange[sector][3] - 1;
      if (xBinMax >= cutAreasTimingEMCaleRange[sector][3]) 
         xBinMax = cutAreasTimingEMCaleRange[sector][1] - 1;

      for (short xBin = xBinMin; xBin <= xBinMax; xBin++)
      {
         for (short yBin = yBinMin; yBin <= yBinMax; yBin++)
         {
            if (!cutAreasTimingEMCale[sector][yBin][xBin]) return false;
         }
      }
      return true;
   }
   // EMCalw
   if (yTower < cutAreasTimingEMCalwRange[sector][0] || 
       yTower >= cutAreasTimingEMCalwRange[sector][1] ||
       zTower < cutAreasTimingEMCalwRange[sector][2] || 
       zTower >= cutAreasTimingEMCalwRange[sector][3]) return true;

   short yBinMin = zTower - zTowerVar;
   short xBinMin = yTower - yTowerVar;
   short yBinMax = zTower + zTowerVar;
   short xBinMax = yTower + yTowerVar;

   if (yBinMin < cutAreasTimingEMCalwRange[sector][2]) 
      yBinMin = cutAreasTimingEMCalwRange[sector][2];
   if (xBinMin < cutAreasTimingEMCalwRange[sector][0]) 
      xBinMin = cutAreasTimingEMCalwRange[sector][0];
   if (yBinMax >= cutAreasTimingEMCalwRange[sector][1]) 
      yBinMax = cutAreasTimingEMCalwRange[sector][3] - 1;
   if (xBinMax >= cutAreasTimingEMCalwRange[sector][3]) 
      xBinMax = cutAreasTimingEMCalwRange[sector][1] - 1;

   for (short xBin = xBinMin; xBin <= xBinMax; xBin++)
   {
      for (short yBin = yBinMin; yBin <= yBinMax; yBin++)
      {
         if (!cutAreasTimingEMCalw[sector][yBin][xBin]) return false;
      }
   }
   return true;
}

bool DeadMapCutter::IsDeadDCTight(const int dcarm, const double zDC, 
                                  const double board, const double alpha,
                                  const double boardVar, const double alphaVar)
{
   if (!useDC) return false;
   if (dcarm == 0) // DCe
   {
      if (zDC >= 0)
      {
         if (board <= cutAreasDCe0X1Range[0] || board >= cutAreasDCe0X1Range[1] ||
             alpha <= cutAreasDCe0X1Range[2] || alpha >= cutAreasDCe0X1Range[3] ||
             board <= cutAreasDCe0X2Range[0] || board >= cutAreasDCe0X2Range[1] ||
             alpha <= cutAreasDCe0X2Range[2] || alpha >= cutAreasDCe0X2Range[3]) return true;

         short yBinX1Min = static_cast<short>((alpha - alphaVar - cutAreasDCe0X1Range[2])/
                                              (cutAreasDCe0X1Range[3] - cutAreasDCe0X1Range[2])*
                                              static_cast<double>(cutAreasDCe0X1.size()));
         short xBinX1Min = static_cast<short>((board - boardVar - cutAreasDCe0X1Range[0])/
                                              (cutAreasDCe0X1Range[1] - cutAreasDCe0X1Range[0])*
                                              static_cast<double>(cutAreasDCe0X1[0].size()));
         short yBinX2Min = static_cast<short>((alpha - alphaVar - cutAreasDCe0X2Range[2])/
                                              (cutAreasDCe0X2Range[3] - cutAreasDCe0X2Range[2])*
                                              static_cast<double>(cutAreasDCe0X2.size()));
         short xBinX2Min = static_cast<short>((board - boardVar - cutAreasDCe0X2Range[0])/
                                              (cutAreasDCe0X2Range[1] - cutAreasDCe0X2Range[0])*
                                              static_cast<double>(cutAreasDCe0X2[0].size()));

         short yBinX1Max = static_cast<short>((alpha + alphaVar - cutAreasDCe0X1Range[2])/
                                              (cutAreasDCe0X1Range[3] - cutAreasDCe0X1Range[2])*
                                              static_cast<double>(cutAreasDCe0X1.size()));
         short xBinX1Max = static_cast<short>((board + boardVar - cutAreasDCe0X1Range[0])/
                                              (cutAreasDCe0X1Range[1] - cutAreasDCe0X1Range[0])*
                                              static_cast<double>(cutAreasDCe0X1[0].size()));
         short yBinX2Max = static_cast<short>((alpha + alphaVar - cutAreasDCe0X2Range[2])/
                                              (cutAreasDCe0X2Range[3] - cutAreasDCe0X2Range[2])*
                                              static_cast<double>(cutAreasDCe0X2.size()));
         short xBinX2Max = static_cast<short>((board + boardVar - cutAreasDCe0X2Range[0])/
                                              (cutAreasDCe0X2Range[1] - cutAreasDCe0X2Range[0])*
                                              static_cast<double>(cutAreasDCe0X2[0].size()));

         if (yBinX1Min < 0) yBinX1Min = 0;
         if (xBinX1Min < 0) xBinX1Min = 0;
         if (yBinX2Min < 0) yBinX2Min = 0;
         if (xBinX2Min < 0) xBinX2Min = 0;

         if (yBinX1Max >= static_cast<short>(cutAreasDCe0X1.size())) 
            yBinX1Max = static_cast<short>(cutAreasDCe0X1.size()) - 1;
         if (xBinX1Max >= static_cast<short>(cutAreasDCe0X1[0].size())) 
            xBinX1Max = static_cast<short>(cutAreasDCe0X1[0].size()) - 1;
         if (yBinX2Max >= static_cast<short>(cutAreasDCe0X2.size())) 
            yBinX2Max = static_cast<short>(cutAreasDCe0X2.size()) - 1;
         if (xBinX2Max >= static_cast<short>(cutAreasDCe0X2[0].size())) 
            xBinX2Max = static_cast<short>(cutAreasDCe0X2[0].size()) - 1;

         for (short xBin = xBinX1Min; xBin <= xBinX1Max; xBin++)
         {
            for (short yBin = yBinX1Min; yBin <= yBinX1Max; yBin++)
            {
               if (cutAreasDCe0X1[yBin][xBin]) return true;
            }
         }

         for (short xBin = xBinX2Min; xBin <= xBinX2Max; xBin++)
         {
            for (short yBin = yBinX2Min; yBin <= yBinX2Max; yBin++)
            {
               if (cutAreasDCe0X2[yBin][xBin]) return true;
            }
         }
         return false;
      }
      else
      {
         if (board <= cutAreasDCe1X1Range[0] || board >= cutAreasDCe1X1Range[1] ||
             alpha <= cutAreasDCe1X1Range[2] || alpha >= cutAreasDCe1X1Range[3] ||
             board <= cutAreasDCe1X2Range[0] || board >= cutAreasDCe1X2Range[1] ||
             alpha <= cutAreasDCe1X2Range[2] || alpha >= cutAreasDCe1X2Range[3]) return true;

         short yBinX1Min = static_cast<short>((alpha - alphaVar - cutAreasDCe1X1Range[2])/
                                              (cutAreasDCe1X1Range[3] - cutAreasDCe1X1Range[2])*
                                              static_cast<double>(cutAreasDCe1X1.size()));
         short xBinX1Min = static_cast<short>((board - boardVar - cutAreasDCe1X1Range[0])/
                                              (cutAreasDCe1X1Range[1] - cutAreasDCe1X1Range[0])*
                                              static_cast<double>(cutAreasDCe1X1[0].size()));
         short yBinX2Min = static_cast<short>((alpha - alphaVar - cutAreasDCe1X2Range[2])/
                                              (cutAreasDCe1X2Range[3] - cutAreasDCe1X2Range[2])*
                                              static_cast<double>(cutAreasDCe1X2.size()));
         short xBinX2Min = static_cast<short>((board - boardVar - cutAreasDCe1X2Range[0])/
                                              (cutAreasDCe1X2Range[1] - cutAreasDCe1X2Range[0])*
                                              static_cast<double>(cutAreasDCe1X2[0].size()));

         short yBinX1Max = static_cast<short>((alpha + alphaVar - cutAreasDCe1X1Range[2])/
                                              (cutAreasDCe1X1Range[3] - cutAreasDCe1X1Range[2])*
                                              static_cast<double>(cutAreasDCe1X1.size()));
         short xBinX1Max = static_cast<short>((board + boardVar - cutAreasDCe1X1Range[0])/
                                              (cutAreasDCe1X1Range[1] - cutAreasDCe1X1Range[0])*
                                              static_cast<double>(cutAreasDCe1X1[0].size()));
         short yBinX2Max = static_cast<short>((alpha + alphaVar - cutAreasDCe1X2Range[2])/
                                              (cutAreasDCe1X2Range[3] - cutAreasDCe1X2Range[2])*
                                              static_cast<double>(cutAreasDCe1X2.size()));
         short xBinX2Max = static_cast<short>((board + boardVar - cutAreasDCe1X2Range[0])/
                                              (cutAreasDCe1X2Range[1] - cutAreasDCe1X2Range[0])*
                                              static_cast<double>(cutAreasDCe1X2[0].size()));

         if (yBinX1Min < 0) yBinX1Min = 0;
         if (xBinX1Min < 0) xBinX1Min = 0;
         if (yBinX2Min < 0) yBinX2Min = 0;
         if (xBinX2Min < 0) xBinX2Min = 0;

         if (yBinX1Max >= static_cast<short>(cutAreasDCe1X1.size())) 
            yBinX1Max = static_cast<short>(cutAreasDCe1X1.size()) - 1;
         if (xBinX1Max >= static_cast<short>(cutAreasDCe1X1[0].size())) 
            xBinX1Max = static_cast<short>(cutAreasDCe1X1[0].size()) - 1;
         if (yBinX2Max >= static_cast<short>(cutAreasDCe1X2.size())) 
            yBinX2Max = static_cast<short>(cutAreasDCe1X2.size()) - 1;
         if (xBinX2Max >= static_cast<short>(cutAreasDCe1X2[0].size())) 
            xBinX2Max = static_cast<short>(cutAreasDCe1X2[0].size()) - 1;

         for (short xBin = xBinX1Min; xBin <= xBinX1Max; xBin++)
         {
            for (short yBin = yBinX1Min; yBin <= yBinX1Max; yBin++)
            {
               if (cutAreasDCe1X1[yBin][xBin]) return true;
            }
         }

         for (short xBin = xBinX2Min; xBin <= xBinX2Max; xBin++)
         {
            for (short yBin = yBinX2Min; yBin <= yBinX2Max; yBin++)
            {
               if (cutAreasDCe1X2[yBin][xBin]) return true;
            }
         }
         return false;
      }
   }
   // DCw
   if (zDC >= 0)
   {
      if (board <= cutAreasDCw0X1Range[0] || board >= cutAreasDCw0X1Range[1] ||
          alpha <= cutAreasDCw0X1Range[2] || alpha >= cutAreasDCw0X1Range[3] ||
          board <= cutAreasDCw0X2Range[0] || board >= cutAreasDCw0X2Range[1] ||
          alpha <= cutAreasDCw0X2Range[2] || alpha >= cutAreasDCw0X2Range[3]) return true;

      short yBinX1Min = static_cast<short>((alpha - alphaVar - cutAreasDCw0X1Range[2])/
                                           (cutAreasDCw0X1Range[3] - cutAreasDCw0X1Range[2])*
                                           static_cast<double>(cutAreasDCw0X1.size()));
      short xBinX1Min = static_cast<short>((board - boardVar - cutAreasDCw0X1Range[0])/
                                           (cutAreasDCw0X1Range[1] - cutAreasDCw0X1Range[0])*
                                           static_cast<double>(cutAreasDCw0X1[0].size()));
      short yBinX2Min = static_cast<short>((alpha - alphaVar - cutAreasDCw0X2Range[2])/
                                           (cutAreasDCw0X2Range[3] - cutAreasDCw0X2Range[2])*
                                           static_cast<double>(cutAreasDCw0X2.size()));
      short xBinX2Min = static_cast<short>((board - boardVar - cutAreasDCw0X2Range[0])/
                                           (cutAreasDCw0X2Range[1] - cutAreasDCw0X2Range[0])*
                                           static_cast<double>(cutAreasDCw0X2[0].size()));

      short yBinX1Max = static_cast<short>((alpha + alphaVar - cutAreasDCw0X1Range[2])/
                                           (cutAreasDCw0X1Range[3] - cutAreasDCw0X1Range[2])*
                                           static_cast<double>(cutAreasDCw0X1.size()));
      short xBinX1Max = static_cast<short>((board + boardVar - cutAreasDCw0X1Range[0])/
                                           (cutAreasDCw0X1Range[1] - cutAreasDCw0X1Range[0])*
                                           static_cast<double>(cutAreasDCw0X1[0].size()));
      short yBinX2Max = static_cast<short>((alpha + alphaVar - cutAreasDCw0X2Range[2])/
                                           (cutAreasDCw0X2Range[3] - cutAreasDCw0X2Range[2])*
                                           static_cast<double>(cutAreasDCw0X2.size()));
      short xBinX2Max = static_cast<short>((board + boardVar - cutAreasDCw0X2Range[0])/
                                           (cutAreasDCw0X2Range[1] - cutAreasDCw0X2Range[0])*
                                           static_cast<double>(cutAreasDCw0X2[0].size()));

      if (yBinX1Min < 0) yBinX1Min = 0;
      if (xBinX1Min < 0) xBinX1Min = 0;
      if (yBinX2Min < 0) yBinX2Min = 0;
      if (xBinX2Min < 0) xBinX2Min = 0;

      if (yBinX1Max >= static_cast<short>(cutAreasDCw0X1.size())) 
         yBinX1Max = static_cast<short>(cutAreasDCw0X1.size()) - 1;
      if (xBinX1Max >= static_cast<short>(cutAreasDCw0X1[0].size())) 
         xBinX1Max = static_cast<short>(cutAreasDCw0X1[0].size()) - 1;
      if (yBinX2Max >= static_cast<short>(cutAreasDCw0X2.size())) 
         yBinX2Max = static_cast<short>(cutAreasDCw0X2.size()) - 1;
      if (xBinX2Max >= static_cast<short>(cutAreasDCw0X2[0].size())) 
         xBinX2Max = static_cast<short>(cutAreasDCw0X2[0].size()) - 1;

      for (short xBin = xBinX1Min; xBin <= xBinX1Max; xBin++)
      {
         for (short yBin = yBinX1Min; yBin <= yBinX1Max; yBin++)
         {
            if (cutAreasDCw0X1[yBin][xBin]) return true;
         }
      }

      for (short xBin = xBinX2Min; xBin <= xBinX2Max; xBin++)
      {
         for (short yBin = yBinX2Min; yBin <= yBinX2Max; yBin++)
         {
            if (cutAreasDCw0X2[yBin][xBin]) return true;
         }
      }
      return false;
   }
   else
   {
      if (board <= cutAreasDCw1X1Range[0] || board >= cutAreasDCw1X1Range[1] ||
          alpha <= cutAreasDCw1X1Range[2] || alpha >= cutAreasDCw1X1Range[3] ||
          board <= cutAreasDCw1X2Range[0] || board >= cutAreasDCw1X2Range[1] ||
          alpha <= cutAreasDCw1X2Range[2] || alpha >= cutAreasDCw1X2Range[3]) return true;

      short yBinX1Min = static_cast<short>((alpha - alphaVar - cutAreasDCw1X1Range[2])/
                                           (cutAreasDCw1X1Range[3] - cutAreasDCw1X1Range[2])*
                                           static_cast<double>(cutAreasDCw1X1.size()));
      short xBinX1Min = static_cast<short>((board - boardVar - cutAreasDCw1X1Range[0])/
                                           (cutAreasDCw1X1Range[1] - cutAreasDCw1X1Range[0])*
                                           static_cast<double>(cutAreasDCw1X1[0].size()));
      short yBinX2Min = static_cast<short>((alpha - alphaVar - cutAreasDCw1X2Range[2])/
                                           (cutAreasDCw1X2Range[3] - cutAreasDCw1X2Range[2])*
                                           static_cast<double>(cutAreasDCw1X2.size()));
      short xBinX2Min = static_cast<short>((board - boardVar - cutAreasDCw1X2Range[0])/
                                           (cutAreasDCw1X2Range[1] - cutAreasDCw1X2Range[0])*
                                           static_cast<double>(cutAreasDCw1X2[0].size()));

      short yBinX1Max = static_cast<short>((alpha + alphaVar - cutAreasDCw1X1Range[2])/
                                           (cutAreasDCw1X1Range[3] - cutAreasDCw1X1Range[2])*
                                           static_cast<double>(cutAreasDCw1X1.size()));
      short xBinX1Max = static_cast<short>((board + boardVar - cutAreasDCw1X1Range[0])/
                                           (cutAreasDCw1X1Range[1] - cutAreasDCw1X1Range[0])*
                                           static_cast<double>(cutAreasDCw1X1[0].size()));
      short yBinX2Max = static_cast<short>((alpha + alphaVar - cutAreasDCw1X2Range[2])/
                                           (cutAreasDCw1X2Range[3] - cutAreasDCw1X2Range[2])*
                                           static_cast<double>(cutAreasDCw1X2.size()));
      short xBinX2Max = static_cast<short>((board + boardVar - cutAreasDCw1X2Range[0])/
                                           (cutAreasDCw1X2Range[1] - cutAreasDCw1X2Range[0])*
                                           static_cast<double>(cutAreasDCw1X2[0].size()));

      if (yBinX1Min < 0) yBinX1Min = 0;
      if (xBinX1Min < 0) xBinX1Min = 0;
      if (yBinX2Min < 0) yBinX2Min = 0;
      if (xBinX2Min < 0) xBinX2Min = 0;

      if (yBinX1Max >= static_cast<short>(cutAreasDCw1X1.size())) 
         yBinX1Max = static_cast<short>(cutAreasDCw1X1.size()) - 1;
      if (xBinX1Max >= static_cast<short>(cutAreasDCw1X1[0].size())) 
         xBinX1Max = static_cast<short>(cutAreasDCw1X1[0].size()) - 1;
      if (yBinX2Max >= static_cast<short>(cutAreasDCw1X2.size())) 
         yBinX2Max = static_cast<short>(cutAreasDCw1X2.size()) - 1;
      if (xBinX2Max >= static_cast<short>(cutAreasDCw1X2[0].size())) 
         xBinX2Max = static_cast<short>(cutAreasDCw1X2[0].size()) - 1;

      for (short xBin = xBinX1Min; xBin <= xBinX1Max; xBin++)
      {
         for (short yBin = yBinX1Min; yBin <= yBinX1Max; yBin++)
         {
            if (cutAreasDCw1X1[yBin][xBin]) return true;
         }
      }

      for (short xBin = xBinX2Min; xBin <= xBinX2Max; xBin++)
      {
         for (short yBin = yBinX2Min; yBin <= yBinX2Max; yBin++)
         {
            if (cutAreasDCw1X2[yBin][xBin]) return true;
         }
      }
      return false;
   }
}

bool DeadMapCutter::IsDeadPC1Tight(const int dcarm, const double ppc1z, const double ppc1phi,
                                   const double ppc1zVar, const double ppc1phiVar)
{
   if (!usePC1) return false;
   if (dcarm == 0) // PC1e
   {
      if (ppc1z <= cutAreasPC1eRange[0] || ppc1z >= cutAreasPC1eRange[1] ||
          ppc1phi <= cutAreasPC1eRange[2] || ppc1phi >= cutAreasPC1eRange[3]) return true;

      short yBinMin = static_cast<short>((ppc1phi - ppc1phiVar - cutAreasPC1eRange[2])/
                                         (cutAreasPC1eRange[3] - cutAreasPC1eRange[2])*
                                         static_cast<double>(cutAreasPC1e.size()));
      short xBinMin = static_cast<short>((ppc1z - ppc1zVar - cutAreasPC1eRange[0])/
                                         (cutAreasPC1eRange[1] - cutAreasPC1eRange[0])*
                                         static_cast<double>(cutAreasPC1e[0].size()));
      short yBinMax = static_cast<short>((ppc1phi + ppc1phiVar - cutAreasPC1eRange[2])/
                                         (cutAreasPC1eRange[3] - cutAreasPC1eRange[2])*
                                         static_cast<double>(cutAreasPC1e.size()));
      short xBinMax = static_cast<short>((ppc1z + ppc1zVar - cutAreasPC1eRange[0])/
                                         (cutAreasPC1eRange[1] - cutAreasPC1eRange[0])*
                                         static_cast<double>(cutAreasPC1e[0].size()));

      if (yBinMin < 0) yBinMin = 0;
      if (xBinMin < 0) xBinMin = 0;
      if (yBinMax >= static_cast<short>(cutAreasPC1e.size())) 
         yBinMax = static_cast<short>(cutAreasPC1e.size()) - 1;
      if (xBinMax >= static_cast<short>(cutAreasPC1e[0].size())) 
         xBinMax = static_cast<short>(cutAreasPC1e[0].size()) - 1;

      for (short xBin = xBinMin; xBin <= xBinMax; xBin++)
      {
         for (short yBin = yBinMin; yBin <= yBinMax; yBin++)
         {
            if (cutAreasPC1e[yBin][xBin]) return true;
         }
      }
      return false;
   }
   // PC1w
   if (ppc1z <= cutAreasPC1wRange[0] || ppc1z >= cutAreasPC1wRange[1] ||
       ppc1phi <= cutAreasPC1wRange[2] || ppc1phi >= cutAreasPC1wRange[3]) return true;

   short yBinMin = static_cast<short>((ppc1phi - ppc1phiVar - cutAreasPC1wRange[2])/
                                      (cutAreasPC1wRange[3] - cutAreasPC1wRange[2])*
                                      static_cast<double>(cutAreasPC1w.size()));
   short xBinMin = static_cast<short>((ppc1z - ppc1zVar - cutAreasPC1wRange[0])/
                                      (cutAreasPC1wRange[1] - cutAreasPC1wRange[0])*
                                      static_cast<double>(cutAreasPC1w[0].size()));
   short yBinMax = static_cast<short>((ppc1phi + ppc1phiVar - cutAreasPC1wRange[2])/
                                      (cutAreasPC1wRange[3] - cutAreasPC1wRange[2])*
                                      static_cast<double>(cutAreasPC1w.size()));
   short xBinMax = static_cast<short>((ppc1z + ppc1zVar - cutAreasPC1wRange[0])/
                                      (cutAreasPC1wRange[1] - cutAreasPC1wRange[0])*
                                      static_cast<double>(cutAreasPC1w[0].size()));

   if (yBinMin < 0) yBinMin = 0;
   if (xBinMin < 0) xBinMin = 0;
   if (yBinMax >= static_cast<short>(cutAreasPC1w.size())) 
      yBinMax = static_cast<short>(cutAreasPC1w.size()) - 1;
   if (xBinMax >= static_cast<short>(cutAreasPC1w[0].size())) 
      xBinMax = static_cast<short>(cutAreasPC1w[0].size()) - 1;

   for (short xBin = xBinMin; xBin <= xBinMax; xBin++)
   {
      for (short yBin = yBinMin; yBin <= yBinMax; yBin++)
      {
         if (cutAreasPC1w[yBin][xBin]) return true;
      }
   }
   return false;
}

bool DeadMapCutter::IsDeadPC2Tight(const double ppc2z, const double ppc2phi,
                                   const double ppc2zVar, const double ppc2phiVar)
{
   if (!usePC2) return false;
   if (ppc2z <= cutAreasPC2Range[0] || ppc2z >= cutAreasPC2Range[1] ||
       ppc2phi <= cutAreasPC2Range[2] || ppc2phi >= cutAreasPC2Range[3]) return true;

   short yBinMin = static_cast<short>((ppc2phi - ppc2phiVar - cutAreasPC2Range[2])/
                                      (cutAreasPC2Range[3] - cutAreasPC2Range[2])*
                                      static_cast<double>(cutAreasPC2.size()));
   short xBinMin = static_cast<short>((ppc2z - ppc2zVar - cutAreasPC2Range[0])/
                                      (cutAreasPC2Range[1] - cutAreasPC2Range[0])*
                                      static_cast<double>(cutAreasPC2[0].size()));
   short yBinMax = static_cast<short>((ppc2phi + ppc2phiVar - cutAreasPC2Range[2])/
                                      (cutAreasPC2Range[3] - cutAreasPC2Range[2])*
                                      static_cast<double>(cutAreasPC2.size()));
   short xBinMax = static_cast<short>((ppc2z + ppc2zVar - cutAreasPC2Range[0])/
                                      (cutAreasPC2Range[1] - cutAreasPC2Range[0])*
                                      static_cast<double>(cutAreasPC2[0].size()));

   if (yBinMin < 0) yBinMin = 0;
   if (xBinMin < 0) xBinMin = 0;
   if (yBinMax >= static_cast<short>(cutAreasPC2.size())) 
      yBinMax = static_cast<short>(cutAreasPC2.size()) - 1;
   if (xBinMax >= static_cast<short>(cutAreasPC2[0].size())) 
      xBinMax = static_cast<short>(cutAreasPC2[0].size()) - 1;

   for (short xBin = xBinMin; xBin <= xBinMax; xBin++)
   {
      for (short yBin = yBinMin; yBin <= yBinMax; yBin++)
      {
         if (cutAreasPC2[yBin][xBin]) return true;
      }
   }
   return false;
}

bool DeadMapCutter::IsDeadPC3Tight(const int dcarm, const double ppc3z, const double ppc3phi,
                                   const double ppc3zVar, const double ppc3phiVar)
{
   if (!usePC3) return false;
   if (dcarm == 0) // PC3e
   {
      if (ppc3z <= cutAreasPC3eRange[0] || ppc3z >= cutAreasPC3eRange[1] ||
          ppc3phi <= cutAreasPC3eRange[2] || ppc3phi >= cutAreasPC3eRange[3]) return true;

      short yBinMin = static_cast<short>((ppc3phi - ppc3phiVar - cutAreasPC3eRange[2])/
                                         (cutAreasPC3eRange[3] - cutAreasPC3eRange[2])*
                                         static_cast<double>(cutAreasPC3e.size()));
      short xBinMin = static_cast<short>((ppc3z - ppc3zVar - cutAreasPC3eRange[0])/
                                         (cutAreasPC3eRange[1] - cutAreasPC3eRange[0])*
                                         static_cast<double>(cutAreasPC3e[0].size()));
      short yBinMax = static_cast<short>((ppc3phi + ppc3phiVar - cutAreasPC3eRange[2])/
                                         (cutAreasPC3eRange[3] - cutAreasPC3eRange[2])*
                                         static_cast<double>(cutAreasPC3e.size()));
      short xBinMax = static_cast<short>((ppc3z + ppc3zVar - cutAreasPC3eRange[0])/
                                         (cutAreasPC3eRange[1] - cutAreasPC3eRange[0])*
                                         static_cast<double>(cutAreasPC3e[0].size()));

      if (yBinMin < 0) yBinMin = 0;
      if (xBinMin < 0) xBinMin = 0;
      if (yBinMax >= static_cast<short>(cutAreasPC3e.size())) 
         yBinMax = static_cast<short>(cutAreasPC3e.size()) - 1;
      if (xBinMax >= static_cast<short>(cutAreasPC3e[0].size())) 
         xBinMax = static_cast<short>(cutAreasPC3e[0].size()) - 1;

      for (short xBin = xBinMin; xBin <= xBinMax; xBin++)
      {
         for (short yBin = yBinMin; yBin <= yBinMax; yBin++)
         {
            if (cutAreasPC3e[yBin][xBin]) return true;
         }
      }
      return false;
   }
   // PC3w
   if (ppc3z <= cutAreasPC3wRange[0] || ppc3z >= cutAreasPC3wRange[1] ||
       ppc3phi <= cutAreasPC3wRange[2] || ppc3phi >= cutAreasPC3wRange[3]) return true;

   short yBinMin = static_cast<short>((ppc3phi - ppc3phiVar - cutAreasPC3wRange[2])/
                                      (cutAreasPC3wRange[3] - cutAreasPC3wRange[2])*
                                      static_cast<double>(cutAreasPC3w.size()));
   short xBinMin = static_cast<short>((ppc3z - ppc3zVar - cutAreasPC3wRange[0])/
                                      (cutAreasPC3wRange[1] - cutAreasPC3wRange[0])*
                                      static_cast<double>(cutAreasPC3w[0].size()));
   short yBinMax = static_cast<short>((ppc3phi + ppc3phiVar - cutAreasPC3wRange[2])/
                                      (cutAreasPC3wRange[3] - cutAreasPC3wRange[2])*
                                      static_cast<double>(cutAreasPC3w.size()));
   short xBinMax = static_cast<short>((ppc3z + ppc3zVar - cutAreasPC3wRange[0])/
                                      (cutAreasPC3wRange[1] - cutAreasPC3wRange[0])*
                                      static_cast<double>(cutAreasPC3w[0].size()));

   if (yBinMin < 0) yBinMin = 0;
   if (xBinMin < 0) xBinMin = 0;
   if (yBinMax >= static_cast<short>(cutAreasPC3w.size())) 
      yBinMax = static_cast<short>(cutAreasPC3w.size()) - 1;
   if (xBinMax >= static_cast<short>(cutAreasPC3w[0].size())) 
      xBinMax = static_cast<short>(cutAreasPC3w[0].size()) - 1;

   for (short xBin = xBinMin; xBin <= xBinMax; xBin++)
   {
      for (short yBin = yBinMin; yBin <= yBinMax; yBin++)
      {
         if (cutAreasPC3w[yBin][xBin]) return true;
      }
   }
   return false;
}

bool DeadMapCutter::IsDeadTOFeTight(const int chamber, const int slat, 
                                    const int chamberVar, const int slatVar)
{
   if (!useTOFe) return false;
   if (chamber < cutAreasTOFeRange[0] || chamber >= cutAreasTOFeRange[1] ||
       slat < cutAreasTOFeRange[2] || slat >= cutAreasTOFeRange[3]) return true;

   short yBinMin = slat - slatVar;
   short xBinMin = chamber - chamberVar;
   short yBinMax = slat + slatVar;
   short xBinMax = chamber + chamberVar;

   if (yBinMin < cutAreasTOFeRange[2]) yBinMin = cutAreasTOFeRange[2];
   if (xBinMin < cutAreasTOFeRange[0]) xBinMin = cutAreasTOFeRange[0];
   if (yBinMax >= cutAreasTOFeRange[1]) yBinMax = cutAreasTOFeRange[3] - 1;
   if (xBinMax >= cutAreasTOFeRange[3]) xBinMax = cutAreasTOFeRange[1] - 1;

   for (short xBin = xBinMin; xBin <= xBinMax; xBin++)
   {
      for (short yBin = yBinMin; yBin <= yBinMax; yBin++)
      {
         if (cutAreasTOFe[yBin][xBin]) return true;
      }
   }
   return false;
}

bool DeadMapCutter::IsDeadTOFwTight(const int chamber, const int strip,
                                    const int chamberVar, const int stripVar)
{
   if (!useTOFw) return false;

   if (chamber < cutAreasTOFwRange[0] || chamber >= cutAreasTOFwRange[1] ||
       strip < cutAreasTOFwRange[2] || strip >= cutAreasTOFwRange[3]) return true;

   short yBinMin = strip - stripVar;
   short xBinMin = chamber - chamberVar;
   short yBinMax = strip + stripVar;
   short xBinMax = chamber + chamberVar;

   if (yBinMin < cutAreasTOFwRange[2]) yBinMin = cutAreasTOFwRange[2];
   if (xBinMin < cutAreasTOFwRange[0]) xBinMin = cutAreasTOFwRange[0];
   if (yBinMax >= cutAreasTOFwRange[1]) yBinMax = cutAreasTOFwRange[3] - 1;
   if (xBinMax >= cutAreasTOFwRange[3]) xBinMax = cutAreasTOFwRange[1] - 1;

   for (short xBin = xBinMin; xBin <= xBinMax; xBin++)
   {
      for (short yBin = yBinMin; yBin <= yBinMax; yBin++)
      {
         if (cutAreasTOFw[yBin][xBin]) return true;
      }
   }
   return false;
}

bool DeadMapCutter::IsDeadEMCalTight(const int dcarm, const int sector, 
                                     const int yTower, const int zTower,
                                     const int yTowerVar, const int zTowerVar)
{
   if (!useEMCal) return false;

   if (dcarm == 0) // EMCale
   {
      if (yTower < cutAreasEMCaleRange[sector][0] || 
          yTower >= cutAreasEMCaleRange[sector][1] ||
          zTower < cutAreasEMCaleRange[sector][2] || 
          zTower >= cutAreasEMCaleRange[sector][3]) return true;

      short yBinMin = zTower - zTowerVar;
      short xBinMin = yTower - yTowerVar;
      short yBinMax = zTower + zTowerVar;
      short xBinMax = yTower + yTowerVar;

      if (yBinMin < cutAreasEMCaleRange[sector][2]) yBinMin = cutAreasEMCaleRange[sector][2];
      if (xBinMin < cutAreasEMCaleRange[sector][0]) xBinMin = cutAreasEMCaleRange[sector][0];
      if (yBinMax >= cutAreasEMCaleRange[sector][1]) yBinMax = cutAreasEMCaleRange[sector][3] - 1;
      if (xBinMax >= cutAreasEMCaleRange[sector][3]) xBinMax = cutAreasEMCaleRange[sector][1] - 1;

      for (short xBin = xBinMin; xBin <= xBinMax; xBin++)
      {
         for (short yBin = yBinMin; yBin <= yBinMax; yBin++)
         {
            if (cutAreasEMCale[sector][yBin][xBin]) return true;
         }
      }
      return false;
   }
   // EMCalw
   if (yTower < cutAreasEMCalwRange[sector][0] || 
       yTower >= cutAreasEMCalwRange[sector][1] ||
       zTower < cutAreasEMCalwRange[sector][2] || 
       zTower >= cutAreasEMCalwRange[sector][3]) return true;

   short yBinMin = zTower - zTowerVar;
   short xBinMin = yTower - yTowerVar;
   short yBinMax = zTower + zTowerVar;
   short xBinMax = yTower + yTowerVar;

   if (yBinMin < cutAreasEMCalwRange[sector][2]) yBinMin = cutAreasEMCalwRange[sector][2];
   if (xBinMin < cutAreasEMCalwRange[sector][0]) xBinMin = cutAreasEMCalwRange[sector][0];
   if (yBinMax >= cutAreasEMCalwRange[sector][1]) yBinMax = cutAreasEMCalwRange[sector][3] - 1;
   if (xBinMax >= cutAreasEMCalwRange[sector][3]) xBinMax = cutAreasEMCalwRange[sector][1] - 1;

   for (short xBin = xBinMin; xBin <= xBinMax; xBin++)
   {
      for (short yBin = yBinMin; yBin <= yBinMax; yBin++)
      {
         if (cutAreasEMCalw[sector][yBin][xBin]) return true;
      }
   }
   return false;
}

bool DeadMapCutter::IsDeadTimingTOFeTight(const int chamber, const int slat, 
                                          const int chamberVar, const int slatVar)
{
   if (!useTOFe) return false;
   if (chamber < cutAreasTimingTOFeRange[0] || chamber >= cutAreasTimingTOFeRange[1] ||
       slat < cutAreasTimingTOFeRange[2] || slat >= cutAreasTimingTOFeRange[3]) return true;

   short yBinMin = slat - slatVar;
   short xBinMin = chamber - chamberVar;
   short yBinMax = slat + slatVar;
   short xBinMax = chamber + chamberVar;

   if (yBinMin < cutAreasTimingTOFeRange[2]) yBinMin = cutAreasTimingTOFeRange[2];
   if (xBinMin < cutAreasTimingTOFeRange[0]) xBinMin = cutAreasTimingTOFeRange[0];
   if (yBinMax >= cutAreasTimingTOFeRange[1]) yBinMax = cutAreasTimingTOFeRange[3] - 1;
   if (xBinMax >= cutAreasTimingTOFeRange[3]) xBinMax = cutAreasTimingTOFeRange[1] - 1;

   for (short xBin = xBinMin; xBin <= xBinMax; xBin++)
   {
      for (short yBin = yBinMin; yBin <= yBinMax; yBin++)
      {
         if (cutAreasTimingTOFe[yBin][xBin]) return true;
      }
   }
   return false;
}

bool DeadMapCutter::IsDeadTimingTOFwTight(const int chamber, const int strip,
                                          const int chamberVar, const int stripVar)
{
   if (!useTOFw) return false;

   if (chamber < cutAreasTimingTOFwRange[0] || chamber >= cutAreasTimingTOFwRange[1] ||
       strip < cutAreasTimingTOFwRange[2] || strip >= cutAreasTimingTOFwRange[3]) return true;

   short yBinMin = strip - stripVar;
   short xBinMin = chamber - chamberVar;
   short yBinMax = strip + stripVar;
   short xBinMax = chamber + chamberVar;

   if (yBinMin < cutAreasTimingTOFwRange[2]) yBinMin = cutAreasTimingTOFwRange[2];
   if (xBinMin < cutAreasTimingTOFwRange[0]) xBinMin = cutAreasTimingTOFwRange[0];
   if (yBinMax >= cutAreasTimingTOFwRange[1]) yBinMax = cutAreasTimingTOFwRange[3] - 1;
   if (xBinMax >= cutAreasTimingTOFwRange[3]) xBinMax = cutAreasTimingTOFwRange[1] - 1;

   for (short xBin = xBinMin; xBin <= xBinMax; xBin++)
   {
      for (short yBin = yBinMin; yBin <= yBinMax; yBin++)
      {
         if (cutAreasTimingTOFw[yBin][xBin]) return true;
      }
   }
   return false;
}

bool DeadMapCutter::IsDeadTimingEMCalTight(const int dcarm, const int sector, 
                                           const int yTower, const int zTower,
                                           const int yTowerVar, const int zTowerVar)
{
   if (!useEMCal) return false;

   if (dcarm == 0) // EMCale
   {
      if (yTower < cutAreasTimingEMCaleRange[sector][0] || 
          yTower >= cutAreasTimingEMCaleRange[sector][1] ||
          zTower < cutAreasTimingEMCaleRange[sector][2] || 
          zTower >= cutAreasTimingEMCaleRange[sector][3]) return true;

      short yBinMin = zTower - zTowerVar;
      short xBinMin = yTower - yTowerVar;
      short yBinMax = zTower + zTowerVar;
      short xBinMax = yTower + yTowerVar;

      if (yBinMin < cutAreasTimingEMCaleRange[sector][2]) 
         yBinMin = cutAreasTimingEMCaleRange[sector][2];
      if (xBinMin < cutAreasTimingEMCaleRange[sector][0]) 
         xBinMin = cutAreasTimingEMCaleRange[sector][0];
      if (yBinMax >= cutAreasTimingEMCaleRange[sector][1]) 
         yBinMax = cutAreasTimingEMCaleRange[sector][3] - 1;
      if (xBinMax >= cutAreasTimingEMCaleRange[sector][3]) 
         xBinMax = cutAreasTimingEMCaleRange[sector][1] - 1;

      for (short xBin = xBinMin; xBin <= xBinMax; xBin++)
      {
         for (short yBin = yBinMin; yBin <= yBinMax; yBin++)
         {
            if (cutAreasTimingEMCale[sector][yBin][xBin]) return true;
         }
      }
      return false;
   }
   // EMCalw
   if (yTower < cutAreasTimingEMCalwRange[sector][0] || 
       yTower >= cutAreasTimingEMCalwRange[sector][1] ||
       zTower < cutAreasTimingEMCalwRange[sector][2] || 
       zTower >= cutAreasTimingEMCalwRange[sector][3]) return true;

   short yBinMin = zTower - zTowerVar;
   short xBinMin = yTower - yTowerVar;
   short yBinMax = zTower + zTowerVar;
   short xBinMax = yTower + yTowerVar;

   if (yBinMin < cutAreasTimingEMCalwRange[sector][2]) 
      yBinMin = cutAreasTimingEMCalwRange[sector][2];
   if (xBinMin < cutAreasTimingEMCalwRange[sector][0]) 
      xBinMin = cutAreasTimingEMCalwRange[sector][0];
   if (yBinMax >= cutAreasTimingEMCalwRange[sector][1]) 
      yBinMax = cutAreasTimingEMCalwRange[sector][3] - 1;
   if (xBinMax >= cutAreasTimingEMCalwRange[sector][3]) 
      xBinMax = cutAreasTimingEMCalwRange[sector][1] - 1;

   for (short xBin = xBinMin; xBin <= xBinMax; xBin++)
   {
      for (short yBin = yBinMin; yBin <= yBinMax; yBin++)
      {
         if (cutAreasTimingEMCalw[sector][yBin][xBin]) return true;
      }
   }
   return false;
}

#endif /* DEAD_MAP_CUTTER_CPP */
