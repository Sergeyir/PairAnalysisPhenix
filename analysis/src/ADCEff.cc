#include "../lib/OutputTool.h"
#include "../lib/ErrorHandler.h"

struct
{
	std::string run_name = "Run7AuAu200";
	std::string file_name = "sum.root";
} Par;

bool IsADCCut(double qtofw)
{
	if (qtofw < 60 || qtofw > 600) return true;
	return false;
}

void ADCEff()
{
	TFile file = TFile(("../data/" + Par.run_name + "/" + Par.file_name).c_str());
	
	TH2F *map2d = (TH2F *) file.Get("tofw_adc");

	TH1F *map1d = (TH1F *) map2d->ProjectionY("proj", 1, map2d->GetXaxis()->GetNbins());
	const double average = map1d->Integral(map1d->GetXaxis()->FindBin(60), map1d->GetXaxis()->FindBin(600))/
		map1d->Integral(1, map1d->GetXaxis()->GetNbins()+1);

	for (int i = 1; i <= map2d->GetXaxis()->GetNbins(); i++)
	{
		const double orig_int = map2d->Integral(i, i, 1., map2d->GetYaxis()->GetNbins()+1);
		if (orig_int <= 0) continue;
		for (int j = 1; j <= map2d->GetYaxis()->GetNbins(); j++)
		{
			if (IsADCCut(map2d->GetYaxis()->GetBinCenter(j))) map2d->SetBinContent(i, j, 0);
		}
		const double remain_int = map2d->Integral(i, i, 1., map2d->GetYaxis()->GetNbins());
		Print(map2d->GetXaxis()->GetBinCenter(i), remain_int/orig_int);
	}

	Print(average);
}
