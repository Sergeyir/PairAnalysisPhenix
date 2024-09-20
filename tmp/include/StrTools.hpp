// $HEADER$
//------------------------------------------------------------------------------------------------
//                               StrTools functions declarations
//------------------------------------------------------------------------------------------------
// StrTools : string tools
//
// ** Free and open code for anyone to use **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic set of functions for handling cpp strings
 **/
//------------------------------------------------------------------------------------------------

#ifndef STR_TOOLS_HPP
#define STR_TOOLS_HPP

#include <string>
#include <sstream>
#include <iomanip>

int utf8_strlen(const std::string& str);
std::string DtoStr(const double val, const int precision = 2);
std::string BtoStr(const bool val);

#endif /*STR_TOOLS_HPP*/
