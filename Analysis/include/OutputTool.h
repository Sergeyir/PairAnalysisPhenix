#pragma once

#include "ErrorHandler.h"

#include <string>
#include <vector>

template <typename... T>
void Print(T... args)
{
	((std::cout << args << " "), ...);
	std::cout << endl;
}

void PrintInfo(std::string info_line)
{
	std::cout << OutputColor::green << "INFO: " << OutputColor::reset << info_line << endl;
}

template <typename... T>
void WriteFile(std::string file_name, T... args)
{
	CheckOutputFile(file_name);
	std::ofstream file(file_name);
	
	((file << args << " "), ...);
	PrintInfo("File " + file_name + " was written");
}

template<typename T> 
void WriteArrayIntoFile(std::string file_name, T *arr, const int size)
{
	CheckOutputFile(file_name);
	std::ofstream file(file_name);
	
	for (int i = 0; i < size; i++) {file << arr[i] << "	";};
	PrintInfo("File " + file_name + " was written");
}

int utf8_strlen(const std::string& str)
{
	int len = 0;
	for (int i=0; i < str.length(); i++, len++)
	{
		const unsigned char c = (unsigned char) str[i];
		if      (c>=0   && c<=127) i+=0;
		else if ((c & 0xE0) == 0xC0) i+=1;
		else if ((c & 0xF0) == 0xE0) i+=2;
		else if ((c & 0xF8) == 0xF0) i+=3;
		else return 0;
	}
	return len;
}

void PrintSimpleSeparator(std::string left_edge = "|", 
		std::string body = "-", 
		std::string right_edge = "|",
		const int length = 80)
{
	std::cout << left_edge;
	for (int i = 0; i < length - utf8_strlen(left_edge) - utf8_strlen(right_edge); i++) std::cout << body;
	std::cout << right_edge << std::endl;
	std::cout << OutputColor::reset;
}

void PrintSeparator(std::string text, 
	std::string color = OutputColor::bold_cyan,
	std::string left_edge = "//", 
	std::string body = "-", 
	std::string right_edge = "//",
	const int length = 80)
{
	std::cout << left_edge;
	for (int i = 0; i < length/2 - utf8_strlen(text)/2 - utf8_strlen(left_edge)*2 - 1; i++) std::cout << body;
	std::cout << left_edge << " " << color << text << OutputColor::reset << " " << right_edge;
	for (int i = 0; i < length - length/2 - utf8_strlen(text) + utf8_strlen(text)/2 - utf8_strlen(right_edge)*2 - 1; i++) std::cout << body;
	std::cout << right_edge << std::endl;
}

void PrintEdgedLine(std::string entry1, std::string entry2,
	std::string left_edge = "|",
	std::string right_edge = "|",
	const int length = 80)
{
	std::cout << left_edge << " " << entry1;
	int space_size = length - utf8_strlen(entry1) - utf8_strlen(entry2) - utf8_strlen(left_edge) - utf8_strlen(right_edge) - 2;
	for (int i = 0; i < space_size; i++) std::cout << " ";
	std::cout << entry2 << " " << right_edge << std::endl;
}

void PrintBigSeparator(std::string text,
	std::string color = OutputColor::bold_cyan,
	std::string ul_corner = "╓", 
	std::string ur_corner = "╖",
	std::string horizontal_line = "─",
	std::string vertical_line = "║",
	std::string dl_corner = "╙",
	std::string dr_corner = "╜")
{
	PrintSimpleSeparator(" " + ul_corner, horizontal_line, ur_corner);
	PrintSeparator(text, color, " " + vertical_line, " ", vertical_line);
	PrintSimpleSeparator(" " + dl_corner, horizontal_line, dr_corner);
}

void PrintFuncPar(TF1 *func, std::string type, std::string arr_name)
{
	std::cout << type << " " << arr_name << "[" << to_string(func->GetNpar()) << "] = {";
	
	for (int i = 0; i < func->GetNpar(); i++)
	{
		std::cout << func->GetParameter(i);
		if (i != func->GetNpar() - 1) std::cout << ", ";
	}
	
	std::cout << "};" << std::endl;
}

void PrintFuncPar2D(std::vector<TF1 *> func_vec, std::string type, std::string arr_name)
{
	std::cout << type << " " << arr_name << 
		"[" << func_vec.size() << "][" << to_string(func_vec[0]->GetNpar()) << "] = " << std::endl;
	
	std::cout << "{" << std::endl;
	
	for (int i = 0; i < func_vec.size(); i++)
	{
		std::cout << "	{";
		for (int j = 0; j < func_vec[i]->GetNpar(); j++)
		{
			std::cout << func_vec[i]->GetParameter(j);
			if (j != func_vec[i]->GetNpar() - 1) std::cout << ", ";
		}
		std::cout << "}";
		if (i != func_vec.size() - 1) std::cout << ", " << std::endl;
		else std::cout << std::endl;
	}
	
	std::cout << "};" << std::endl;
}

template <typename T>
void PrintArr(T *arr, const int size, std::string name = "", bool next_line = true)
{
	if (name != "") std::cout << name << "[" << size << "] = ";

	std::cout << "{";
	for (int i = 0; i < size - 1; i++)
	{
		std::cout << arr[i] << ", ";
	}
	std::cout << arr[size-1] << "}";

	if (name != "") std::cout << ";" << std::endl;
}

template <typename T>
void PrintVec(std::vector<T> vec, std::string name = "", bool next_line = true)
{
	if (name != "") std::cout << name << "[" << vec.size() << "] = ";

	std::cout << "{";
	for (int i = 0; i < vec.size() - 1; i++)
	{
		std::cout << vec[i] << ", ";
	}
	std::cout << vec[vec.size()-1] << "}";

	if (name != "") std::cout << ";";
	
	if (next_line) std::cout << std::endl;
}

template <typename T>
void PrintVec2D(std::vector<std::vector<T>> vec, std::string name = "", bool next_line = true)
{
	if (name != "") std::cout << name << "[" << vec.size() << "][" << vec[0].size() << "] = ";

	std::cout << "{" << std::endl;
	for (int i = 0; i < vec.size() - 1; i++)
	{
		PrintVec(vec[i], "", false);
		std::cout << ", " << std::endl;
	}
	PrintVec(vec[vec.size()-1], "", false);
	std::cout << std::endl << "}";

	if (name != "") std::cout << ";";
	if (next_line) std::cout << std::endl;
}
