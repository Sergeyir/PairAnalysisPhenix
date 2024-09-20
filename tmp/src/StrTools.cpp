// $SOURCE$
//------------------------------------------------------------------------------------------------
//                               StrTools functions realisations
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

#ifndef STR_TOOLS_CPP
#define STR_TOOLS_CPP

#include "../include/StrTools.hpp"

int utf8_strlen(const std::string& str)
{
   int length = 0;
   for (long unsigned int i=0; i < str.length(); i++, length++)
   {
      const unsigned char c = (unsigned char) str[i];
      if      (c>=0 && c<=127)     i+=0;
      else if ((c & 0xE0) == 0xC0) i+=1;
      else if ((c & 0xF0) == 0xE0) i+=2;
      else if ((c & 0xF8) == 0xF0) i+=3;
      else return 0;
   }
   return length;
}

std::string DtoStr(const double val, const int precision)
{
   std::stringstream ssval;
   ssval << std::fixed << std::setprecision(precision) << val;
   return ssval.str();
}

std::string BtoStr(const bool val)
{
   if (val) return "true";
   return "false";
}

#endif /*STR_TOOLS_CPP*/
