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
   if (options[0] == '1')
   {
      doCutDC = true;
      SetDeadAreas("data/Parameters/Deadmaps/" + runName + "/DCeX1,zDC>=0.txt", 
                   cutAreasDCe0X1, cutAreasDCe0X1Range);
      SetDeadAreas("data/Parameters/Deadmaps/" + runName + "/DCeX1,zDC<0.txt", 
                   cutAreasDCe1X1, cutAreasDCe1X1Range);
      SetDeadAreas("data/Parameters/Deadmaps/" + runName + "/DCwX1,zDC>=0.txt", 
                   cutAreasDCw0X1, cutAreasDCw0X1Range);
      SetDeadAreas("data/Parameters/Deadmaps/" + runName + "/DCwX1,zDC<0.txt", 
                   cutAreasDCw1X1, cutAreasDCw1X1Range);
      SetDeadAreas("data/Parameters/Deadmaps/" + runName + "/DCeX2,zDC>=0.txt", 
                   cutAreasDCe0X2, cutAreasDCe0X2Range);
      SetDeadAreas("data/Parameters/Deadmaps/" + runName + "/DCeX2,zDC<0.txt", 
                   cutAreasDCe1X2, cutAreasDCe1X2Range);
      SetDeadAreas("data/Parameters/Deadmaps/" + runName + "/DCwX2,zDC>=0.txt", 
                   cutAreasDCw0X2, cutAreasDCw0X2Range);
      SetDeadAreas("data/Parameters/Deadmaps/" + runName + "/DCwX2,zDC<0.txt", 
                   cutAreasDCw1X2, cutAreasDCw1X2Range);
   }
   else 
   {
      CppTools::PrintInfo("DeadMapCutter: cuts for DC were specified to be not initialized");
      doCutDC = false;
   }

   if (options[1] == '1')
   {
      doCutPC1 = true;
      SetDeadAreas("data/Parameters/Deadmaps/" + runName + "/PC1e.txt", 
                   cutAreasPC1e, cutAreasPC1eRange);
      SetDeadAreas("data/Parameters/Deadmaps/" + runName + "/PC1w.txt", 
                   cutAreasPC1w, cutAreasPC1wRange);
   }
   else 
   {
      CppTools::PrintInfo("DeadMapCutter: cuts for PC1 were specified to be not initialized");
      doCutPC1 = false;
   }

   if (options[2] == '1')
   {
      doCutPC2 = true;
      SetDeadAreas("data/Parameters/Deadmaps/" + runName + "/PC2.txt", 
                   cutAreasPC2, cutAreasPC2Range);
   }
   else 
   {
      CppTools::PrintInfo("DeadMapCutter: cuts for PC2 were specified to be not initialized");
      doCutPC2 = false;
   }

   if (options[3] == '1')
   {
      doCutPC3 = true;
      SetDeadAreas("data/Parameters/Deadmaps/" + runName + "/PC3e.txt", 
                   cutAreasPC3e, cutAreasPC3eRange);
      SetDeadAreas("data/Parameters/Deadmaps/" + runName + "/PC3w.txt", 
                   cutAreasPC3w, cutAreasPC3wRange);
   }
   else 
   {
      CppTools::PrintInfo("DeadMapCutter: cuts for PC3 were specified to be not initialized");
      doCutPC3 = false;
   }

   if (options[4] == '1')
   {
      doCutTOFe = true;
      SetDeadAreas("data/Parameters/Deadmaps/" + runName + "/TOFe.txt", 
                   cutAreasTOFe, cutAreasTOFeRange);
   }
   else 
   {
      CppTools::PrintInfo("DeadMapCutter: cuts for TOFe were specified to be not initialized");
      doCutTOFe = false;
   }

   if (options[5] == '1')
   {
      doCutTOFw = true;
      SetDeadAreas("data/Parameters/Deadmaps/" + runName + "/TOFw.txt", 
                   cutAreasTOFw, cutAreasTOFwRange);
   }
   else 
   {
      CppTools::PrintInfo("DeadMapCutter: cuts for TOFw were specified to be not initialized");
      doCutTOFw = false;
   }

   if (options[6] == '1')
   {
      doCutEMCal = true;
      for (int i = 0; i < 4; i++)
      {
         SetDeadAreas("data/Parameters/Deadmaps/" + runName + "/EMCale" + 
                      std::to_string(i) + ".txt", cutAreasEMCale[i], cutAreasEMCaleRange[i]);
         SetDeadAreas("data/Parameters/Deadmaps/" + runName + "/EMCalw" + 
                      std::to_string(i) + ".txt", cutAreasEMCalw[i], cutAreasEMCalwRange[i]);
      }
   }
   else 
   {
      CppTools::PrintInfo("DeadMapCutter: cuts for EMCal were specified to be not initialized");
      doCutEMCal = false;
   }
}

bool DeadMapCutter::IsDeadDC(const int dcarm, const double zDC, 
                             const double board, const double alpha)
{
   if (!doCutDC) return false;
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
   if (!doCutPC1) return false;
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
   if (!doCutPC2) return false;
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
   if (!doCutPC3) return false;
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
   if (!doCutTOFe) return false;
   if (chamber <= cutAreasTOFeRange[0] || chamber >= cutAreasTOFeRange[1] ||
       slat <= cutAreasTOFeRange[2] || slat >= cutAreasTOFeRange[3]) return true;

   const short yBin = static_cast<short>((slat - cutAreasTOFeRange[2])/
                                         (cutAreasTOFeRange[3] - cutAreasTOFeRange[2])*
                                         static_cast<double>(cutAreasTOFe.size()));
   const short xBin = static_cast<short>((chamber - cutAreasTOFeRange[0])/
                                         (cutAreasTOFeRange[1] - cutAreasTOFeRange[0])*
                                         static_cast<double>(cutAreasTOFe[0].size()));
   return cutAreasTOFe[yBin][xBin];
}

bool DeadMapCutter::IsDeadTOFw(const int chamber, const int strip)
{
   if (!doCutTOFw) return false;

   if (chamber <= cutAreasTOFwRange[0] || chamber >= cutAreasTOFwRange[1] ||
       strip <= cutAreasTOFwRange[2] || strip >= cutAreasTOFwRange[3]) return true;

   const short yBin = static_cast<short>((strip - cutAreasTOFwRange[2])/
                                         (cutAreasTOFwRange[3] - cutAreasTOFwRange[2])*
                                         static_cast<double>(cutAreasTOFw.size()));
   const short xBin = static_cast<short>((chamber - cutAreasTOFwRange[0])/
                                         (cutAreasTOFwRange[1] - cutAreasTOFwRange[0])*
                                         static_cast<double>(cutAreasTOFw[0].size()));
   return cutAreasTOFw[yBin][xBin];
}

bool DeadMapCutter::IsDeadEMCal(const int dcarm, const int sector, 
                                const int ytower, const int ztower)
{
   if (!doCutEMCal) return false;
   if (dcarm == 0) // EMCale
   {
      if (ytower <= cutAreasEMCaleRange[sector][0] || 
          ytower >= cutAreasEMCaleRange[sector][1] ||
          ztower <= cutAreasEMCaleRange[sector][2] || 
          ztower >= cutAreasEMCaleRange[sector][3]) return true;

      const short yBin = static_cast<short>((ztower - cutAreasEMCaleRange[sector][2])/
                                            (cutAreasEMCaleRange[sector][3] - 
                                             cutAreasEMCaleRange[sector][2])*
                                            static_cast<double>(cutAreasEMCale[sector].size()));
      const short xBin = static_cast<short>((ytower - cutAreasEMCaleRange[sector][0])/
                                            (cutAreasEMCaleRange[sector][1] - 
                                             cutAreasEMCaleRange[sector][0])*
                                            static_cast<double>(cutAreasEMCale[sector][0].size()));
      return cutAreasEMCale[sector][yBin][xBin];
   }
   // EMCalw
   if (ytower <= cutAreasEMCalwRange[sector][0] || 
       ytower >= cutAreasEMCalwRange[sector][1] ||
       ztower <= cutAreasEMCalwRange[sector][2] || 
       ztower >= cutAreasEMCalwRange[sector][3]) return true;

   const short yBin = static_cast<short>((ztower - cutAreasEMCalwRange[sector][2])/
                                         (cutAreasEMCalwRange[sector][3] - 
                                          cutAreasEMCalwRange[sector][2])*
                                         static_cast<double>(cutAreasEMCalw[sector].size()));
   const short xBin = static_cast<short>((ytower - cutAreasEMCalwRange[sector][0])/
                                         (cutAreasEMCalwRange[sector][1] - 
                                          cutAreasEMCalwRange[sector][0])*
                                         static_cast<double>(cutAreasEMCalw[sector][0].size()));
   return cutAreasEMCalw[sector][yBin][xBin];
}

void DeadMapCutter::SetDeadAreas(const std::string& inputFileName, 
                                 std::vector<std::vector<bool>>& cutAreas, double *ranges)
{
   CppTools::CheckInputFile(inputFileName);

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
}

#endif /* DEAD_MAP_CUTTER_CPP */
