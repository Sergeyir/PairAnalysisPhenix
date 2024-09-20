// $SOURCE$
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

#ifndef IOTOOLS_CPP
#define IOTOOLS_CPP

#include "../include/IOTools.hpp"

int GetTerminalWidth()
{
   struct winsize terminalWindow;
   ioctl(STDOUT_FILENO, TIOCGWINSZ, &(terminalWindow));
   return terminalWindow.ws_col;
}

double *ReadFileIntoArray(const std::string& fileName, const int size)
{
   CheckInputFile(fileName);
   double *buff = new double[size];
   std::ifstream file(fileName);
   
   for (int i = 0; i < size; i++) {file >> buff[i];};
   return buff;
}

void PrintInfo(const std::string& message)
{
   std::cout << OutputColor::GREEN << " INFO: " << OutputColor::RESET << message << std::endl;
}

void PrintSimpleSeparator(const std::string& leftEdge, const std::string& body, 
                          const std::string& rightEdge, int length)
{
   if (length < 0) length = GetTerminalWidth();
   
   std::cout << leftEdge;
   for (int i = 0; i < length - utf8_strlen(leftEdge) - utf8_strlen(rightEdge); i++) 
   {
      std::cout << body;
   }
   std::cout << rightEdge << std::endl;
}

void PrintSeparator(const std::string& text, const std::string& color, const std::string& leftEdge, 
                    const std::string& body, const std::string& rightEdge, int length)
{
   if (length < 0) length = GetTerminalWidth();
   
   std::cout << leftEdge;
   for (int i = 0; i < length/2 - utf8_strlen(text)/2 - utf8_strlen(leftEdge); i++) 
   {
      std::cout << body;
   }
   std::cout << color << text << OutputColor::RESET;
   for (int i = 0; i < length - length/2 - utf8_strlen(text) + utf8_strlen(text)/2 - utf8_strlen(rightEdge); i++) 
   {
      std::cout << body;
   }
   std::cout << rightEdge << std::endl;
}

void PrintEdgedLine(const std::string& entry1, const std::string& entry2, 
                    const std::string& leftEdge, const std::string& rightEdge, int length)
{
   if (length < 0) length = GetTerminalWidth();
   
   std::cout << leftEdge << " " << entry1;
   int space_size = length - utf8_strlen(entry1) - utf8_strlen(entry2) - 
                    utf8_strlen(leftEdge) - utf8_strlen(rightEdge) - 2;
   for (int i = 0; i < space_size; i++) std::cout << " ";
   std::cout << entry2 << " " << rightEdge << std::endl;
}

void PrintBigSeparator(const std::string& text, const std::string& ULCorner, 
                       const std::string& URCorner, const std::string& horizontalLine, 
                       const std::string& verticalLine, const std::string& DLCorner, 
                       const std::string& DRCorner)
{
   PrintSimpleSeparator(" " + ULCorner, horizontalLine, URCorner);
   PrintSeparator(text, "", " " + verticalLine, " ", verticalLine);
   PrintSimpleSeparator(" " + DLCorner, horizontalLine, DRCorner);
}

#endif /*IOTOOLS_CPP*/
