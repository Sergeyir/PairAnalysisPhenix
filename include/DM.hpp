#include "ErrorHandler.hpp"
#include "IOTools.hpp"
#include "StrTools.hpp"
#include "Table.hpp"
#include "MathTools.hpp"

#include "DeadAreasCuts.h"

#include "TCanvasPrinter.hpp"

void CheckHists(TH2F *hist1, TH2F *hist2);
void SetHistStyle(TH2F *hist, const std::string& title, 
                  const std::string& xTitle, const std::string &yTitle);
double GetNormRatio(const double ratio);
double CalculateUncertaintyFromProj(TH1F *dataCutDistr, TH1F *simCutDistr);
double CalculateUncertainty1Proj(TH2F *dataDistr, TH2F *simDistr, 
                                 TH2F *dataCutDistr, TH2F* simCutDistr, 
                                 const std::string& detectorName, const std::string& title, 
                                 const std::string& xTitle, const std::string& yTitle,
                                 const int rebinX);
double CalculateUncertainty2Proj(TH2F *dataDistr, TH2F *simDistr, 
                                 TH2F *dataCutDistr, TH2F* simCutDistr, 
                                 const std::string& detectorName, const std::string& title, 
                                 const std::string& xTitle, const std::string& yTitle,
                                 const int rebinX, const int rebinY);
double GetUncertainty(const std::string& histName, const std::string &detectorName, 
                      const std::string &title, const std::string xTitle, 
                      const std::string yTitle, bool sysCalc2Proj, 
                      bool (*cut_func)(const double, const double), 
                      const int rebinX = 1, const int rebinY = 1);
double  GetUncertainty(const std::string& histName, const std::string &detectorName, 
                       const std::string &title, const std::string xTitle, 
                       const std::string yTitle, bool sysCalc2Proj, 
                       bool (*cut_func)(const double, const double, const double), 
                       const double auxVal, const int rebinX = 1, const int rebinY = 1);
double GetUncertainty(const std::string& histName, const std::string &detectorName, 
                      const std::string &title, const std::string xTitle, 
                      const std::string yTitle, bool sysCalc2Proj, 
                      bool (*cut_func)(const double, const double, const double, const double), 
                      const double auxVal1, const double auxVal2,
                      const int rebinX = 1, const int rebinY = 1);
double GetUncertainty(const std::string& histName, const std::string &detectorName, 
                      const std::string &title, const std::string xTitle, 
                      const std::string yTitle, bool sysCalc2Proj, 
                      bool (*cut_func)(const double, const double, const int, 
                                       const double, const double), 
                      const double auxVal1, const double auxVal2, const int auxVal3,
                      const int rebinX = 1, const int rebinY = 1);
void DM();
