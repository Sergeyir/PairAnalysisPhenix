// $SOURCE$
//------------------------------------------------------------------------------------------------
//                                   Table class declaration
//------------------------------------------------------------------------------------------------
// Table
//
// ** Free and open code for anyone to use **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic tool for putting text output in an ascii table
 **/
//------------------------------------------------------------------------------------------------

#ifndef TABLE_CPP
#define TABLE_CPP

#include "../include/Table.hpp"

Table::Table(const unsigned short numberOfColumns, const int width) 
{
   nColumns = numberOfColumns;
   if (width > 0) rowLength = width;
   else rowLength = GetTerminalWidth();
};

void Table::Begin(const std::string& name)
{
   PrintSimpleSeparator(" " + ULBorderCorner, horizontalBorder, URBorderCorner);
   PrintSeparator(name, OutputColor::GREEN, " " + verticalBorder, " ", verticalBorder);
}

void Table::End()
{
   std::cout << " " << DLBorderCorner;
   for (int i = 0; i < rowLength - utf8_strlen(rightDoubleAdjVerticalBorder) - utf8_strlen(leftDoubleAdjVerticalBorder) - 1; i++)
   {
      if (i != 0 && i % static_cast<int>(cellSize) == 0 && 
          rowLength - i >= cellSize) std::cout << upAdjHorizontalBorder;
      else std::cout << horizontalBorder;
   }
   std::cout << DRBorderCorner << std::endl;
}

Table::~Table() {};

#endif /*TABLE_CPP*/
