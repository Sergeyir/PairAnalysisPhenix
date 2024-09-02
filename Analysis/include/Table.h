#pragma once

#include "OutputTool.h"
#include "ErrorHandler.h"

template<typename... Ts>
class Table
{
	private:
	
	std::string ul_corner = "╔";
	std::string ur_corner = "╗";
	std::string hline = "═";
	std::string vline = "║";
	std::string vsep = "│";
	std::string clvline = "╟";
	std::string crvline = "╢";
	std::string dclvline = "╠";
	std::string dcrvline = "╣";
	std::string hsep = "─";
	std::string dl_corner = "╚";
	std::string dr_corner = "╝";
	std::string cuhline = "╤";
	std::string cdhline = "╧";
	std::string cross = "┼";
	
	int row_length;

	double cell_size;

	public:
	
	Table(int width = 80) 
	{
		row_length = width;
	};

	void Begin(std::string name)
	{
		PrintSimpleSeparator(" " + ul_corner, hline, ur_corner);
		PrintSeparator(name, OutputColor::green, " " + vline, " ", vline);
	}

	void PrintHeader(Ts... args)
	{
		constexpr int size = sizeof...(args);

		std::string dummy[size] = {(std::string) args...};
		
		//checks whether or not next cell shoudl be printed
		bool check = true;
		//number of the next cell to be printed
		int ncell = 0;

		//size of the cell
		cell_size = static_cast<double>(row_length - 2*utf8_strlen(vline))/size;

		//finishing the title
		std::cout << " " << dclvline;
		for (int i = 0; i < row_length - utf8_strlen(dclvline) - utf8_strlen(dcrvline) - 1; i++)
		{
			if (i != 0 && i % static_cast<int>(cell_size) == 0 && row_length - i >= cell_size) std::cout << cuhline;
			else std::cout << hline;
		}
		std::cout << dcrvline << std::endl;

		//printing header
		std::cout << " " << vline;
		for (int i = 0; i < row_length - utf8_strlen(dclvline) - utf8_strlen(dcrvline) - 1; i++)
		{
			if (i != 0 && i % static_cast<int>(cell_size) == 0 && row_length - i >= cell_size) 
			{
				std::cout << vsep;
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
		std::cout << vline << std::endl;

		//printing header separator
		std::cout << " " << clvline;
		for (int i = 0; i < row_length - utf8_strlen(dclvline) - utf8_strlen(dcrvline) - 1; i++)
		{
			if (i != 0 && i % static_cast<int>(cell_size) == 0 && row_length - i >= cell_size) std::cout << cross;
			else std::cout << hsep;
		}
		std::cout << crvline << std::endl;
	}

	void PrintRow(Ts... args)
	{
		constexpr int size = sizeof...(args);

		std::string dummy[size] = {(std::string) args...};
		
		//checks whether or not next cell shoudl be printed
		bool check = true;
		//number of the next cell to be printed
		int ncell = 0;

		//printing the row
		std::cout << " " << vline;
		for (int i = 0; i < row_length - utf8_strlen(dclvline) - utf8_strlen(dcrvline) - 1; i++)
		{
			if (i != 0 && i % static_cast<int>(cell_size) == 0 && row_length - i >= cell_size) 
			{
				std::cout << vsep;
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
		std::cout << vline << std::endl;
	}

	void End()
	{
		std::cout << " " << dl_corner;
		for (int i = 0; i < row_length - utf8_strlen(dclvline) - utf8_strlen(dcrvline) - 1; i++)
		{
			if (i != 0 && i % static_cast<int>(cell_size) == 0 && row_length - i >= cell_size) std::cout << cdhline;
			else std::cout << hline;
		}
		std::cout << dr_corner << std::endl;
	}

};
