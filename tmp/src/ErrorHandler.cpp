// $SOURCE$
//------------------------------------------------------------------------------------------------
//                               Error handler functions declarations
//------------------------------------------------------------------------------------------------
// ErrorHandler
//
// ** Free and open code for anyone to use **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic set of functions for handling errors and warnings
 **/
//------------------------------------------------------------------------------------------------

#ifndef ERROR_HANDLER_CPP
#define ERROR_HANDLER_CPP

#include "../include/ErrorHandler.hpp"

//Print an error with message
void PrintError(const std::string& message)
{
   std::cerr << OutputColor::BOLD_RED << " ERROR: " << OutputColor::WHITE << message << std::endl;
   exit(1);
}

//Print warning with message
void PrintWarning(const std::string& message)
{
   std::cerr << OutputColor::BOLD_MAGENTA << " WARNING: " << OutputColor::WHITE << message << std::endl;
}

//Checks if the file exists
bool CheckInputFile(const std::string& name, const bool closeAfterFail)
{
   std::ifstream file(name.c_str());
   if(!file.is_open()) 
   {
      if (closeAfterFail) PrintError("File " + name + " not found");
      else PrintWarning("File " + name + " not found");
      return false;
   }
   return true;
}

//Checks if the file can be created
void CheckOutputFile(const std::string& name, const std::ios_base::openmode openmode)
{
   std::ofstream file(name.c_str(), openmode);
   if(!file.is_open()) PrintError("File " + name + " cannot be created");
}

#endif /*ERROR_HANDLER_CPP*/
