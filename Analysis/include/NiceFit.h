#pragma once

#include <vector>
#include <string>
#include "TF1.h"

namespace NiceFit
{
	void Perform(TGraphErrors *data, TF1* fit_func, std::string quiet_options = "RQMBN", std::string final_options = "RQMBE+", unsigned const int iter_num = 100, const double min_variation = 0.01, double *fixed_par = NULL)
	{
		data->Fit(fit_func, quiet_options.c_str());
		
		for (int iter = 1; iter <= iter_num; iter++)
		{
			const double mult = (1. + min_variation)*iter_num/iter;
			for (int i = 0; i < fit_func->GetNpar(); i++)
			{
				fit_func->SetParLimits(i, fit_func->GetParameter(i)/mult, fit_func->GetParameter(i)*mult);
			}
			data->Fit(fit_func, quiet_options.c_str());
		}

		data->Fit(fit_func, final_options.c_str());
	}
};
