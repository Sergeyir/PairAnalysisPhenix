#pragma once

template <typename... Ts> 
double Maximum(Ts... args)
{
	constexpr int size = sizeof...(args);
	double entries[size] = {static_cast<double>(args)...};
	double result = entries[0];
	for (double val : entries) if (val > result) result = val;
	return result;
}

template <typename... Ts> 
double Minimum(Ts... args)
{
	constexpr int size = sizeof...(args);
	double entries[size] = {static_cast<double>(args)...};
	double result = entries[0];
	for (double val : entries) if (val < result) result = val;
	return result;
}

template <typename... Ts> 
double Average(Ts... args)
{
	constexpr int size = sizeof...(args);
	double entries[size] = {static_cast<double>(args)...};
	double result = 0.;
	for (double val : entries) result += val/static_cast<double>(size);
	return result;
}

double AverageVec(std::vector<double> vec, const int size)
{
	double result = 0;
	for (int i = 0; i < size; i++)
	{
		result += vec[i]/static_cast<double>(size);
	}
	return result;
}

double GetNormRatio(double ratio)
{
	if (ratio >= 1.) return ratio;
	else return 1./ratio;
}

//uncertainty propagation
template<typename... Ts>
double ErrPropagation(Ts... args)
{
	constexpr int size = sizeof...(args);
	double entries[size] = {static_cast<double>(args)...};
	double prop = 0;
	for (double var : entries) prop += var*var;
	return sqrt(prop);
}

template<typename... Ts>
double RMS(Ts... args)
{
	constexpr int size = sizeof...(args);
	double entries[size] = {static_cast<double>(args)...};
	double rms = 0;
	for (double var : entries) rms += var*var;
	return sqrt(rms/static_cast<double>(size));
}

double RMSv(std::vector<double> vec)
{
	double rms = 0.;
	for (double var : vec) rms += var*var;
	return sqrt(rms/static_cast<double>(vec.size()));
}
