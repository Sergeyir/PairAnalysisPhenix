// $HEADER$
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

#ifndef TABLE_HPP
#define TABLE_HPP

#include "OutputColor.hpp"
#include "IOTools.hpp"
#include "ErrorHandler.hpp"

class Table
{
   public:
   
   Table(const unsigned short numberOfColumns, const int width = -1);
   
   void Begin(const std::string& name);
   
   template<typename... Ts>
   void PrintHeader(Ts... args)
   {
      constexpr unsigned short size = sizeof...(args);
      if (size != nColumns) 
      {
         PrintError("Number of arguments is not equal to the number of columns:" + 
                    std::to_string(size) + " vs " + std::to_string(nColumns));
      }

      std::string dummy[size] = {(std::string) args...};
      
      //checks whether or not next cell should be printed
      bool check = true;
      //number of the next cell to be printed
      int ncell = 0;

      //size of the cell
      cellSize = static_cast<double>(rowLength - 2*utf8_strlen(verticalBorder))/size;

      //finishing the title
      std::cout << " " << rightDoubleAdjVerticalBorder;
      for (int i = 0; i < rowLength - utf8_strlen(rightDoubleAdjVerticalBorder) - 
           utf8_strlen(leftDoubleAdjVerticalBorder) - 1; i++)
      {
         if (i != 0 && i % static_cast<int>(cellSize) == 0 && 
             rowLength - i >= cellSize) 
         {
            std::cout << downAdjHorizontalBorder;
         }
         else std::cout << horizontalBorder;
      }
      std::cout << leftDoubleAdjVerticalBorder << std::endl;

      //printing header
      std::cout << " " << verticalBorder;
      for (int i = 0; i < rowLength - utf8_strlen(rightDoubleAdjVerticalBorder) - 
           utf8_strlen(leftDoubleAdjVerticalBorder) - 1; i++)
      {
         if (i != 0 && i % static_cast<int>(cellSize) == 0 && rowLength - i >= cellSize) 
         {
            std::cout << verticalSeparator;
            check = true;
         }
         else if (check && ncell < size)
         {
            std::cout << " " << dummy[ncell];
            i += utf8_strlen(" " + dummy[ncell]) - 1;
            ncell++;
            check = false;
         }
         else std::cout << " ";
      }
      std::cout << verticalBorder << std::endl;

      //printing header separator
      std::cout << " " << leftAdjVerticalBorder;
      for (int i = 0; i < rowLength - utf8_strlen(rightDoubleAdjVerticalBorder) - 
           utf8_strlen(leftDoubleAdjVerticalBorder) - 1; i++)
      {
         if (i != 0 && i % static_cast<int>(cellSize) == 0 && 
             rowLength - i >= cellSize) std::cout << cross;
         else std::cout << horizontalSeparator;
      }
      std::cout << rightAdjVerticalBorder << std::endl;
   }

   template<typename... Ts>
   void PrintRow(Ts... args)
   {
      constexpr unsigned short size = sizeof...(args);
      
      if (size != nColumns) 
      {
         PrintError("Number of arguments is not equal to the number of columns:" + 
                    std::to_string(size) + " vs " + std::to_string(nColumns));
      }

      std::string dummy[size] = {(std::string) args...};
      
      //checks whether or not next cell shoudl be printed
      bool check = true;
      //number of the next cell to be printed
      int ncell = 0;

      //printing the row
      std::cout << " " << verticalBorder;
      for (int i = 0; i < rowLength - utf8_strlen(rightDoubleAdjVerticalBorder) - 
           utf8_strlen(leftDoubleAdjVerticalBorder) - 1; i++)
      {
         if (i != 0 && i % static_cast<int>(cellSize) == 0 && rowLength - i >= cellSize) 
         {
            std::cout << verticalSeparator;
            check = true;
         }
         else if (check && ncell < size)
         {
            std::cout << " " << dummy[ncell];
            i += utf8_strlen(" " + dummy[ncell]) - 1;
            ncell++;
            check = false;
         }
         else std::cout << " ";
      }
      std::cout << verticalBorder << std::endl;
   }
   
   void End();

   virtual ~Table();

   protected:
   
   std::string ULBorderCorner = "╔";
   std::string URBorderCorner = "╗";
   std::string horizontalBorder = "═";
   std::string verticalBorder = "║";
   std::string verticalSeparator = "│";
   std::string leftAdjVerticalBorder = "╟";
   std::string rightAdjVerticalBorder = "╢";
   std::string rightDoubleAdjVerticalBorder = "╠";
   std::string leftDoubleAdjVerticalBorder = "╣";
   std::string horizontalSeparator = "─";
   std::string DLBorderCorner = "╚";
   std::string DRBorderCorner = "╝";
   std::string downAdjHorizontalBorder = "╤";
   std::string upAdjHorizontalBorder = "╧";
   std::string cross = "┼";
   
   int nColumns;
   int rowLength;
   double cellSize;
};

#endif /*TABLE_HPP*/
