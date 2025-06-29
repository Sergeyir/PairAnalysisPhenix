/** 
 *  @file   FitFunc.cpp
 *  @brief  Contains realisations of functions that are used for various fits for approximation of signals in PairAnalysisPhenix project such as relativistic breit wigner (RBW), its convolution with gaus (RBWGDF), RBWGDF + various backgrounds, etc. Functions are defined in a special way so that they can be passed in a constructor of TF1 (https://root.cern.ch/doc/master/classTF1.html#aa8905d28455ed7be02019f20b9cc5827)
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysisPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/

#ifndef FIT_FUNC_CPP
#define FIT_FUNC_CPP

#include "../include/FitFunc.hpp"

double FitFunc::RBW(double *x, double *par)
{
   return par[0]*TMath::BreitWignerRelativistic(x[0], par[1], par[2])/
          TMath::BreitWignerRelativistic(par[1], par[1], par[2]);
}

double FitFunc::RBWConvGaus(double *x, double *par)
{
   /// if sigma is small comapared to gamma then there is no need to perform convolution
   /// ff sigma << gamma function estimation may even not be performed 
   /// correctly due to limited precision of double
   if (par[3] < par[2]/1e3)
   {
      return RBW(x, par);
   }

   // integration sum
   double sum = 0.;
   // normalization constant
   double norm = 0.;
   
   // integrating over t
   for (double t = - sigmalizedConvolutionRange*par[3]; t < sigmalizedConvolutionRange*par[3]; 
        t += 2.*sigmalizedConvolutionRange*par[3]/static_cast<double>(integralNumberOfIterations))
   {
      sum += TMath::Gaus(t, 0., par[3])*
         TMath::BreitWigner(x[0] - t, par[1], par[2]);
      // TMath::Gaus already has a maximum 1 at the mean so no 
      // need to add it to normalization constant
      norm += TMath::BreitWigner(par[1] - t, par[1], par[2]);
   }

   return par[0]*sum/norm;
}

double FitFunc::Pol2(double *x, double *par)
{
   return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

double FitFunc::Pol3(double *x, double *par)
{
   return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
}

double FitFunc::Gaus(double *x, double *par)
{
   return par[0]*TMath::Gaus(x[0], par[1], par[2]);
}

double FitFunc::RBWConvGausBGPol2(double *x, double *par)
{
   return RBWConvGaus(x, par) + Pol2(x, &par[4]);
}

double FitFunc::RBWConvGausBGPol3(double *x, double *par)
{
   return RBWConvGaus(x, par) + Pol3(x, &par[4]);
}

double FitFunc::RBWConvGausBGGaus(double *x, double *par)
{
   // TMath::Gaus already has a maximum 1 at the mean so no need to normalize it
   return RBWConvGaus(x, par) + Gaus(x, &par[4]);
}

#endif /* FIT_FUNC_CPP */
