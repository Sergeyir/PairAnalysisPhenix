#pragma once

#include "TCanvas.h"

void PrintCanvas (TCanvas *canv, const std::string outputFileName, bool printPdf = true)
{
	canv->SaveAs((outputFileName + ".png").c_str());
	if (printPdf) canv->SaveAs((outputFileName + ".tmp.pdf").c_str());

	if (printPdf) system(
		("gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.5 \
		-dNOPAUSE -dQUIET -dBATCH -dPrinted=false -sOutputFile=" + 
		outputFileName + ".pdf " + outputFileName + ".tmp.pdf && rm " + 
		outputFileName + ".tmp.pdf &").c_str());
}
