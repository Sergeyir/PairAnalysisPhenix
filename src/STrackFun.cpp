// $SOURCE$
//------------------------------------------------------------------------------------------------
//                                     STrackFun realisation
//------------------------------------------------------------------------------------------------
// STrackFun - single track functions
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic functions used in analysis for single particle tracks
 **/
//------------------------------------------------------------------------------------------------

#ifndef S_TRACK_FUN_CPP
#define S_TRACK_FUN_CPP

#include "../include/STrackFun.hpp"

bool IsHit(const double dval)
{
   if (dval < -900) return false;
   return true;
}

bool IsMatch(const double sdphi, const double sdz, const double sdphiMax, const double sdzMax)
{
   if (abs(sdphi) > sdphiMax || abs(sdz) > sdzMax) return false;
   return true;
}   

double TransformProb(double prob)
{
   if (prob > 1.) prob = 1.;
   else if (prob < 0.) prob = 0.;
   return prob;
}

void CutDCDeadAreas(TH2F *hist, 
                    bool (*cutFunc)(const double, const double, const double, const double), 
                    const double phi, const double zed)
{
   for (int i = 1; i <= hist->GetXaxis()->GetNbins(); i++)
   {
      for (int j = 1; j <= hist->GetYaxis()->GetNbins(); j++)
      {
         const double xVal = hist->GetXaxis()->GetBinCenter(i);
         const double yVal = hist->GetYaxis()->GetBinCenter(j);

         if (cutFunc(phi, zed, xVal, yVal))
         {
            hist->SetBinContent(i, j, 0.);
         }
      }
   }
}

bool IsQualityCut(const int qual)
{
   if (qual != 63 && qual != 31) return true;
   return false;
}

double GetM2Mean(const double pt, const double *par)
{
   return par[0] + pt*par[1];
}

double GetM2Sigma(const double pt, const double m2Mean, const double *par)
{
   return 2.*sqrt(pow(par[0]/par[3]*m2Mean*pt, 2) + 
      pow(par[1]/par[3]*m2Mean, 2)*(1.+m2Mean/pt/pt) + 
      pow(par[2]*2.9972e-4/par[4]*pt, 2)*(m2Mean + pt*pt));
}

#endif
