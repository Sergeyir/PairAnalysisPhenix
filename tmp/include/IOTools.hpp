// $HEADER$
//------------------------------------------------------------------------------------------------
//                               IOTools functions declarations
//------------------------------------------------------------------------------------------------
// IOTools: input-output tools
//
// ** Free and open code for anyone to use **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic set of functions for handling input and output in the terminal
 **/
//------------------------------------------------------------------------------------------------

#ifndef IOTOOLS_HPP
#define IOTOOLS_HPP

#include <stdio.h>
#include <sys/ioctl.h>
#include <unistd.h>

#include <vector>

#include "ErrorHandler.hpp"
#include "StrTools.hpp"

template <typename... T>
void ReadFile(const std::string& file_name, T&... args)
{
   CheckInputFile(file_name);
   std::ifstream file(file_name);
   
   ((file >> args), ...);
}

template <typename... T>
void Print(T... args)
{
   ((std::cout << args << " "), ...);
   std::cout << std::endl;
}

int GetTerminalWidth();
double *ReadFileIntoArray(const std::string& fileName, const int size);
void PrintInfo(const std::string& message);
void PrintSimpleSeparator(const std::string& leftEdge = "|", const std::string& body = "-", 
                          const std::string& rightEdge = "|", int length = -1);
void PrintSeparator(const std::string& text, const std::string& color = "", 
                    const std::string& leftEdge = "//",  const std::string& body = "-", 
                    const std::string& rightEdge = "//", int length = -1);
void PrintEdgedLine(const std::string& entry1, const std::string& entry2, 
                    const std::string& leftEdge = "|", const std::string& rightEdge = "|", 
                    int length = -1);
void PrintBigSeparator(const std::string& text, const std::string& ULCorner = "╓", 
                       const std::string& URCorner = "╖", const std::string& horizontalLine = "─", 
                       const std::string& verticalLine = "║", const std::string& DLCorner = "╙", 
                       const std::string& DRCorner = "╜");

#endif /*IOTOOLS_HPP*/
