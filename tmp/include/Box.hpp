// $HEADER$
//------------------------------------------------------------------------------------------------
//                                   Box class declaration
//------------------------------------------------------------------------------------------------
// Box
//
// ** Free and open code for anyone to use **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic tool for putting text output in an ascii box
 **/
//------------------------------------------------------------------------------------------------

#ifndef BOX_HPP
#define BOX_HPP

#include "OutputColor.hpp"
#include "ErrorHandler.hpp"
#include "IOTools.hpp"
#include "StrTools.hpp"

class Box
{
   public:

   Box(const int boxWidth = -1);
   Box(const std::string& boxName, const int boxWidth = -1);
   void SetName(const std::string& boxName);
   void AddEntry(const std::string& entryName, const double entryValue, 
                 const unsigned int precision = 2);
   void AddEntry(const std::string& entryName, const int entryValue);
   void AddEntry(const std::string& entryName, const std::string& entryValue);
   void AddEntry(const std::string& entryName, const bool entryValue);
   void Print(const std::string& titleColor = OutputColor::GREEN);
   
   virtual ~Box();
   
   protected:
   
   std::string name; 
   std::vector<std::string> entryNames, entryValues;
   unsigned short width;
   void CheckEntry();
};

#endif /*BOX_HPP*/
