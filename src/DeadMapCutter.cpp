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

DeadMapCutter::DeadMapCutter(const std::string& runName, const std::string& options)
{
   if (options[0] == '1')
   {
      SetDeadAreas("DCe0, zDC<0.txt", cutAreasDCe0);
      SetDeadAreas("DCe1, zDC>0.txt", cutAreasDCe0);
   }
}

#endif /* DEAD_MAP_CUTTER_CPP */
