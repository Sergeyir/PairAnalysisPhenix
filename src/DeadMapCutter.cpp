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
      SetDeadAreas("data/Parameters/Deadmaps/" + runName + "/DCe, zDC>=0.txt", 
                   cutAreasDCe0, cutAreasDCe0Range);
      SetDeadAreas("data/Parameters/Deadmaps/" + runName + "/DCe, zDC<0.txt", 
                   cutAreasDCe1, cutAreasDCe1Range);
      SetDeadAreas("data/Parameters/Deadmaps/" + runName + "/DCw, zDC>=0.txt", 
                   cutAreasDCw0, cutAreasDCw0Range);
      SetDeadAreas("data/Parameters/Deadmaps/" + runName + "/DCw, zDC<0.txt", 
                   cutAreasDCw1, cutAreasDCw1Range);
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
                   cutAreasPC1e, cutAreasPC3eRange);
      SetDeadAreas("data/Parameters/Deadmaps/" + runName + "/PC3w.txt", 
                   cutAreasPC1w, cutAreasPC3wRange);
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
      SetDeadAreas("data/Parameters/Deadmaps/" + runName + "/Slat.txt", 
                   cutSlatsTOFe, cutSlatsRange);
   }
   else 
   {
      CppTools::PrintInfo("DeadMapCutter: cuts for TOFe were specified to be not initialized");
      doCutTOFe = false;
   }

   if (options[5] == '1')
   {
      doCutTOFw = true;
      SetDeadAreas("data/Parameters/Deadmaps/" + runName + "/TOFw, ptofy<100.txt", 
                   cutAreasTOFw0, cutAreasTOFw0Range);
      SetDeadAreas("data/Parameters/Deadmaps/" + runName + "/TOFw, ptofy>100.txt", 
                   cutAreasTOFw1, cutAreasTOFw1Range);
      SetDeadAreas("data/Parameters/Deadmaps/" + runName + "/Strip.txt", 
                   cutStripsTOFw, cutStripsRange);
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

bool DeadMapCutter::IsDeadDC(const int dcarm, const double zed, 
                             const double board, const double alpha)
{
   if (!doCutDC) return false;
   if (dcarm == 1) // DCe
   {
      if (zed >= 0)
      {
         if (board < cutAreasDCe0Range[0] || board > cutAreasDCe0Range[1] ||
             alpha < cutAreasDCe0Range[2] || alpha > cutAreasDCe0Range[3]) return true;

         const short yBin = static_cast<short>((alpha - cutAreasDCe0Range[2])/
                                               (cutAreasDCe0Range[3] - cutAreasDCe0Range[2])*
                                               static_cast<double>(cutAreasDCe0.size()));
         const short xBin = static_cast<short>((board - cutAreasDCe0Range[0])/
                                               (cutAreasDCe0Range[1] - cutAreasDCe0Range[0])*
                                               static_cast<double>(cutAreasDCe0[0].size()));
         return cutAreasDCe0[yBin][xBin];
      }
      else
      {
         if (board < cutAreasDCe1Range[0] || board > cutAreasDCe1Range[1] ||
             alpha < cutAreasDCe1Range[2] || alpha > cutAreasDCe1Range[3]) return true;

         const short yBin = static_cast<short>((alpha - cutAreasDCe1Range[2])/
                                               (cutAreasDCe1Range[3] - cutAreasDCe1Range[2])*
                                               static_cast<double>(cutAreasDCe1.size()));
         const short xBin = static_cast<short>((board - cutAreasDCe1Range[0])/
                                               (cutAreasDCe1Range[1] - cutAreasDCe1Range[0])*
                                               static_cast<double>(cutAreasDCe1[0].size()));
         return cutAreasDCe1[yBin][xBin];
      }
   }
   // DCw
   if (zed < 0)
   {
      if (board < cutAreasDCw0Range[0] || board > cutAreasDCw0Range[1] ||
          alpha < cutAreasDCw0Range[2] || alpha > cutAreasDCw0Range[3]) return true;

      const short yBin = static_cast<short>((alpha - cutAreasDCw0Range[2])/
                                            (cutAreasDCw0Range[3] - cutAreasDCw0Range[2])*
                                            static_cast<double>(cutAreasDCw0.size()));
      const short xBin = static_cast<short>((board - cutAreasDCw0Range[0])/
                                            (cutAreasDCw0Range[1] - cutAreasDCw0Range[0])*
                                            static_cast<double>(cutAreasDCw0[0].size()));
      return cutAreasDCw0[yBin][xBin];
   }
   else
   {
      if (board < cutAreasDCw1Range[0] || board > cutAreasDCw1Range[1] ||
          alpha < cutAreasDCw1Range[2] || alpha > cutAreasDCw1Range[3]) return true;

      const short yBin = static_cast<short>((alpha - cutAreasDCw1Range[2])/
                                            (cutAreasDCw1Range[3] - cutAreasDCw1Range[2])*
                                            static_cast<double>(cutAreasDCw1.size()));
      const short xBin = static_cast<short>((board - cutAreasDCw1Range[0])/
                                            (cutAreasDCw1Range[1] - cutAreasDCw1Range[0])*
                                            static_cast<double>(cutAreasDCw1[0].size()));
      return cutAreasDCw1[yBin][xBin];
   }
}

bool DeadMapCutter::IsDeadPC1(const int dcarm, const double ppc1z, const double ppc1phi)
{
   if (!doCutPC1) return false;
   if (dcarm == 1) // PC1e
   {
      if (ppc1z < cutAreasPC1eRange[0] || ppc1z > cutAreasPC1eRange[1] ||
          ppc1phi < cutAreasPC1eRange[2] || ppc1phi > cutAreasPC1eRange[3]) return true;

      const short yBin = static_cast<short>((ppc1phi - cutAreasPC1eRange[2])/
                                            (cutAreasPC1eRange[3] - cutAreasPC1eRange[2])*
                                            static_cast<double>(cutAreasPC1e.size()));
      const short xBin = static_cast<short>((ppc1z - cutAreasPC1eRange[0])/
                                            (cutAreasPC1eRange[1] - cutAreasPC1eRange[0])*
                                            static_cast<double>(cutAreasPC1e[0].size()));
      return cutAreasPC1e[yBin][xBin];
   }
   // PC1w
   if (ppc1z < cutAreasPC1wRange[0] || ppc1z > cutAreasPC1wRange[1] ||
       ppc1phi < cutAreasPC1wRange[2] || ppc1phi > cutAreasPC1wRange[3]) return true;

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
   if (ppc2z < cutAreasPC2Range[0] || ppc2z > cutAreasPC2Range[1] ||
       ppc2phi < cutAreasPC2Range[2] || ppc2phi > cutAreasPC2Range[3]) return true;

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
   if (dcarm == 1) // PC3e
   {
      if (ppc3z < cutAreasPC3eRange[0] || ppc3z > cutAreasPC3eRange[1] ||
          ppc3phi < cutAreasPC3eRange[2] || ppc3phi > cutAreasPC3eRange[3]) return true;

      const short yBin = static_cast<short>((ppc3phi - cutAreasPC3eRange[2])/
                                            (cutAreasPC3eRange[3] - cutAreasPC3eRange[2])*
                                            static_cast<double>(cutAreasPC3e.size()));
      const short xBin = static_cast<short>((ppc3z - cutAreasPC3eRange[0])/
                                            (cutAreasPC3eRange[1] - cutAreasPC3eRange[0])*
                                            static_cast<double>(cutAreasPC3e[0].size()));
      return cutAreasPC3e[yBin][xBin];
   }
   // PC3w
   if (ppc3z < cutAreasPC3wRange[0] || ppc3z > cutAreasPC3wRange[1] ||
       ppc3phi < cutAreasPC3wRange[2] || ppc3phi > cutAreasPC3wRange[3]) return true;

   const short yBin = static_cast<short>((ppc3phi - cutAreasPC3wRange[2])/
                                         (cutAreasPC3wRange[3] - cutAreasPC3wRange[2])*
                                         static_cast<double>(cutAreasPC3w.size()));
   const short xBin = static_cast<short>((ppc3z - cutAreasPC3wRange[0])/
                                         (cutAreasPC3wRange[1] - cutAreasPC3wRange[0])*
                                         static_cast<double>(cutAreasPC3w[0].size()));
   return cutAreasPC3w[yBin][xBin];
}

bool DeadMapCutter::IsDeadTOFe(const double ptofy, const double ptofz, const int slat)
{
   if (!doCutTOFe) return false;
   if (ptofy < cutAreasTOFeRange[0] || ptofy > cutAreasTOFeRange[1] ||
       ptofz < cutAreasTOFeRange[2] || ptofz > cutAreasTOFeRange[3]) return true;

   const short yBin = static_cast<short>((ptofz - cutAreasTOFeRange[2])/
                                         (cutAreasTOFeRange[3] - cutAreasTOFeRange[2])*
                                         static_cast<double>(cutAreasTOFe.size()));
   const short xBin = static_cast<short>((ptofy - cutAreasTOFeRange[0])/
                                         (cutAreasTOFeRange[1] - cutAreasTOFeRange[0])*
                                         static_cast<double>(cutAreasTOFe[0].size()));
   return cutAreasTOFe[yBin][xBin];
}

bool DeadMapCutter::IsDeadTOFw(const double ptofwy, const double ptofwz, const int striptofw)
{
   if (!doCutTOFw) return false;
   if (ptofwy < 100)
   {
      if (ptofwy < cutAreasTOFw0Range[0] || ptofwy > cutAreasTOFw0Range[1] ||
          ptofwz < cutAreasTOFw0Range[2] || ptofwz > cutAreasTOFw0Range[3]) return true;

      const short yBin = static_cast<short>((ptofwz - cutAreasTOFw0Range[2])/
                                            (cutAreasTOFw0Range[3] - cutAreasTOFw0Range[2])*
                                            static_cast<double>(cutAreasTOFw0.size()));
      const short xBin = static_cast<short>((ptofwy - cutAreasTOFw0Range[0])/
                                            (cutAreasTOFw0Range[1] - cutAreasTOFw0Range[0])*
                                            static_cast<double>(cutAreasTOFw0[0].size()));
      return cutAreasTOFw0[yBin][xBin];
   }
   if (ptofwy < cutAreasTOFw1Range[0] || ptofwy > cutAreasTOFw1Range[1] ||
       ptofwz < cutAreasTOFw1Range[2] || ptofwz > cutAreasTOFw1Range[3]) return true;

   const short yBin = static_cast<short>((ptofwz - cutAreasTOFw1Range[2])/
                                         (cutAreasTOFw1Range[3] - cutAreasTOFw1Range[2])*
                                         static_cast<double>(cutAreasTOFw1.size()));
   const short xBin = static_cast<short>((ptofwy - cutAreasTOFw1Range[0])/
                                         (cutAreasTOFw1Range[1] - cutAreasTOFw1Range[0])*
                                         static_cast<double>(cutAreasTOFw1[0].size()));
   return cutAreasTOFw1[yBin][xBin];
}

bool DeadMapCutter::IsDeadEMCal(const int dcarm, const int sector, const int ytower, const int ztower)
{
   if (!doCutEMCal) return false;
   if (dcarm == 1) // EMCale
   {
      if (ytower < cutAreasEMCaleRange[sector][0] || 
          ytower > cutAreasEMCaleRange[sector][1] ||
          ztower < cutAreasEMCaleRange[sector][2] || 
          ztower > cutAreasEMCaleRange[sector][3]) return true;

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
   if (ytower < cutAreasEMCalwRange[sector][0] || 
       ytower > cutAreasEMCalwRange[sector][1] ||
       ztower < cutAreasEMCalwRange[sector][2] || 
       ztower > cutAreasEMCalwRange[sector][3]) return true;

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

void DeadMapCutter::SetDeadAreas(const std::string& inputFileName, 
                                 std::vector<bool>& cutAreas, double *ranges)
{
   CppTools::CheckInputFile(inputFileName);

   std::ifstream inputFile(inputFileName);

   int xNBins;
   bool isUnexpectedEndOfFile = false;

   if (!(inputFile >> xNBins >> ranges[0] >> ranges[1]))
   {
      isUnexpectedEndOfFile = true;
   }

   cutAreas.resize(xNBins);

   for (int i = 0; i < xNBins; i++)
   {
      bool tmp;
      if (!(inputFile >> tmp)) 
      {
         isUnexpectedEndOfFile = true;
         break;
      }
      cutAreas[i] = tmp;
   }

   if (isUnexpectedEndOfFile) 
   {
      CppTools::PrintError("DeadMapCutter: Unexpected end of file: " + inputFileName);
   }
}

#endif /* DEAD_MAP_CUTTER_CPP */
