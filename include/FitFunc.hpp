/** 
 *  @file   FitFunc.hpp
 *  @brief  Contains declarations of functions that are used for various fits for approximation of signals in PairAnalysisPhenix project such as relativistic breit wigner (RBW), its convolution with gaus (RBWGDF), RBWGDF + various backgrounds, etc. Functions are defined in a special way so that they can be passed in a constructor of TF1 (https://root.cern.ch/doc/master/classTF1.html#aa8905d28455ed7be02019f20b9cc5827)
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysisPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/

#ifndef FIT_FUNC_HPP
#define FIT_FUNC_HPP

#include <cmath>

#include "TMath.h"

/* @namespace FitFunc
 *
 * @brief Contains all functions for approximation of signals, and variables for them. Functions are defined in a special way so that they can be passed in a constructor of TF1 (https://root.cern.ch/doc/master/classTF1.html#aa8905d28455ed7be02019f20b9cc5827)
 */
namespace FitFunc
{
   /// number of iterations for the integration if there is an integral in a function (for convolution, or integral functions)
   unsigned int integralNumberOfIterations = 1000;
   /// sigmalized range from the mean of a gaus to use for convolution
   double sigmalizedConvolutionRange = 10.;
   /* @brief Calculates value of relativistic Breit-Wigner distribution at x[0]. The maximum of this distribution is equal to scale parameter at the median.
    * @param[in] x[0] transverse momentum [GeV/c]
    * @param[in] par[0] scale parameter
    * @param[in] par[1] median [GeV/c^2]
    * @param[in] par[2] gamma [GeV/c^2]
    * This function is defined in a special way so that it can be passed in a constructor of TF1 (https://root.cern.ch/doc/master/classTF1.html#aa8905d28455ed7be02019f20b9cc5827)
    */
   double RBW(double *x, double *par);
   /* @brief Calculates value of relativistic Breit-Wigner distribution convoluted with gaus at x[0]. The maximum of this distribution is equal to scale parameter at the median.
    * @param[in] x[0] transverse momentum [GeV/c]
    * @param[in] par[0] scale parameter
    * @param[in] par[1] median of RBW [GeV/c^2]
    * @param[in] par[2] gamma of RBW [GeV/c^2]
    * @param[in] par[3] sigma of gaus for convolution [GeV/c^2]
    * This function is defined in a special way so that it can be passed in a constructor of TF1 (https://root.cern.ch/doc/master/classTF1.html#aa8905d28455ed7be02019f20b9cc5827)
    */
   double RBWConvGaus(double *x, double *par);
   /* @brief Calculates value of 2nd order polynomial (a + b*x + c*x^2) at x[0]
    * @param[in] par[0] a
    * @param[in] par[1] b
    * @param[in] par[2] c 
    * This function is defined in a special way so that it can be passed in a constructor of TF1 (https://root.cern.ch/doc/master/classTF1.html#aa8905d28455ed7be02019f20b9cc5827)
    */
   double Pol2(double *x, double *par);
   /* @brief 3rd order polynomial (a + b*x + c*x^2 + d*x^3)
    * @param[in] par[0] a
    * @param[in] par[1] b
    * @param[in] par[2] c 
    * @param[in] par[3] d
    * This function is defined in a special way so that it can be passed in a constructor of TF1 (https://root.cern.ch/doc/master/classTF1.html#aa8905d28455ed7be02019f20b9cc5827)
    */
   double Pol3(double *x, double *par);
   /* @brief Calculates value of gaus at x[0]. The maximum of this distribution is equal to scale parameter at the mean.
    * @param[in] scale parameter
    * @param[in] par[1] mean
    * @param[in] par[2] sigma
    * This function is defined in a special way so that it can be passed in a constructor of TF1 (https://root.cern.ch/doc/master/classTF1.html#aa8905d28455ed7be02019f20b9cc5827)
    */
   double Gaus(double *x, double *par);
   /* @brief Relativistic Breit-Wigner distribution convoluted with gaus + 2nd order polynomial background (a + b*x + c*x^2).  Used for approximation of real data signals.
    * @param[in] x[0] transverse momentum [GeV/c]
    * @param[in] par[0] scale parameter for signal
    * @param[in] par[1] median of RBW [GeV/c^2]
    * @param[in] par[2] gamma of RBW [GeV/c^2]
    * @param[in] par[3] sigma of gaus for convolution [GeV/c^2]
    * @param[in] par[4] a
    * @param[in] par[5] b
    * @param[in] par[6] c 
    * This function is defined in a special way so that it can be passed in a constructor of TF1 (https://root.cern.ch/doc/master/classTF1.html#aa8905d28455ed7be02019f20b9cc5827)
    */
   double RBWConvGausBGPol2(double *x, double *par);
   /* @brief Relativistic Breit-Wigner distribution convoluted with gaus + 3rd order polynomial background (a + b*x + c*x^2 + d*x^3). Used for approximation of real data signals.
    * @param[in] x[0] transverse momentum [GeV/c]
    * @param[in] par[0] scale parameter for signal
    * @param[in] par[1] median of RBW [GeV/c^2]
    * @param[in] par[2] gamma of RBW [GeV/c^2]
    * @param[in] par[3] sigma of gaus for convolution [GeV/c^2]
    * @param[in] par[4] a
    * @param[in] par[5] b
    * @param[in] par[6] c 
    * @param[in] par[7] d
    * This function is defined in a special way so that it can be passed in a constructor of TF1 (https://root.cern.ch/doc/master/classTF1.html#aa8905d28455ed7be02019f20b9cc5827)
    */
   double RBWConvGausBGPol3(double *x, double *par);
   /* @brief Relativistic Breit-Wigner distribution convoluted with gaus + gaus background. Used for approximation of simulated data signals.
    * @param[in] x[0] transverse momentum [GeV/c]
    * @param[in] par[0] scale parameter for signal
    * @param[in] par[1] median of RBW [GeV/c^2]
    * @param[in] par[2] gamma of RBW [GeV/c^2]
    * @param[in] par[3] sigma of gaus for convolution [GeV/c^2]
    * @param[in] par[4] scale parameter for background gaus
    * @param[in] par[5] mean of background gaus [GeV/c^2]
    * @param[in] par[6] sigma of background gaus [GeV/c^2]
    * This function is defined in a special way so that it can be passed in a constructor of TF1 (https://root.cern.ch/doc/master/classTF1.html#aa8905d28455ed7be02019f20b9cc5827)
    */
   double RBWConvGausBGGaus(double *x, double *par);
}

#endif /* FIT_FUNC_HPP */
