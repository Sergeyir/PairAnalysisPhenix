// $SOURCE$
//------------------------------------------------------------------------------------------------
//                                Time functions declarations
//------------------------------------------------------------------------------------------------
// Time
//
// ** Free and open code for anyone to use **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic set of functions for measuring time
 **/
//------------------------------------------------------------------------------------------------

#ifndef TIME_CPP
#define TIME_CPP

#include "../include/Time.hpp"

chrono_t GetCurrentTime()
{
   return std::chrono::high_resolution_clock::now();
}

void PrintCurrentTime()
{
   const std::time_t currentTime = std::time(nullptr);
   std::cout << std::asctime(std::localtime(&currentTime));
}

void PrintTimeDuration(const std::chrono::_V2::system_clock::time_point start, 
                       const std::chrono::_V2::system_clock::time_point stop)
{
   const unsigned int duration = (unsigned int) 
      std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
   
   std::cout << duration/3600000000000 << "h " << 
                (duration % 3600000000)/60000000000 << "m " << 
                (duration % 60000000)/1000000 << "s " << 
                (duration % 1000000)/1000 << "ms " << 
                duration % 1000 << "μs" << std::endl;
}

#endif /*TIME_CPP*/
