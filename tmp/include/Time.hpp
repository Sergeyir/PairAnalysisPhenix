// $HEADER$
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

#ifndef TIME_HPP
#define TIME_HPP

#include <ctime>
#include <chrono>
#include <iostream>

typedef std::chrono::_V2::system_clock::time_point chrono_t;

chrono_t GetCurrentTime();
void PrintCurrentTime();
void PrintTimeDuration(const std::chrono::_V2::system_clock::time_point start, 
                       const std::chrono::_V2::system_clock::time_point stop);

#endif /*TIME_HPP*/
