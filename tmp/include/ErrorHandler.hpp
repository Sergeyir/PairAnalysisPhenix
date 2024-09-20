// $HEADER$
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

#ifndef ERROR_HANDLER_HPP
#define ERROR_HANDLER_HPP

#include <iostream>
#include <fstream>

#include "OutputColor.hpp"

void PrintError(const std::string& message);
void PrintWarning(const std::string& message);
bool CheckInputFile(const std::string& name, const bool closeAfterFail = true);
void CheckOutputFile(const std::string& name, 
                     const std::ios_base::openmode openmode = std::ios::out);

#endif /*ERROR_HANDLER_HPP*/
