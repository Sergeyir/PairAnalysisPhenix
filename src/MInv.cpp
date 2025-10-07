/** 
 *  @file   MInv.cpp 
 *  @brief  Contains declarations of functions that are used for subtracting background and for merging invariant mass histograms in a given centrality, z, and r regions
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysis).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef M_INV_CPP
#define M_INV_CPP

#include "MInv.hpp"

TH1D *MInv::Merge(TFile *inputFile, const std::string& methodName, 
                  const std::string& decayMode, const int cMin, const int cMax, 
                  const int zMin, const int zMax, const int rMin, const int rMax, 
                  const double pTMin, const double pTMax,
                  TH1D*& distrMInvMergedFG, TH1D*& distrMInvMergedBG,
                  TH1D*& distrMInvMergedFGLR, TH1D*& distrMInvMergedBGLR)
{
   TH1D *distrMInvMerged = nullptr;
   // iterating over CabanaBoy centrality bins
   for (int c = cMin; c <= cMax; c++)
   {
      const std::string cName = (c > 9) ? std::to_string(c) : "0" + std::to_string(c);

      // iterating over CabanaBoy z_{vtx} bins
      for (int z = zMin; z <= zMax; z++)
      {
         const std::string zName = (z > 9) ? std::to_string(z) : "0" + std::to_string(z);

         // iterating over CabanaBoy r_{vtx} bins
         for (int r = rMin; r <= rMax; r++)
         {
            const std::string rName = (r > 9) ? std::to_string(r) : "0" + std::to_string(r);

            const std::string distrMInvVsPTFGName = 
               "c" + cName + "_z" + zName + "_r" + rName + "/" + 
                methodName + ": " + decayMode + "_FG12";

            TH2F *distrMInvVsPTFG = 
               static_cast<TH2F *>(inputFile->Get(distrMInvVsPTFGName.c_str()));

            if (!distrMInvVsPTFG)
            {
               CppTools::PrintError("Histogram named " + distrMInvVsPTFGName + 
                                    " does not exist in file" + 
                                    static_cast<std::string>(inputFile->GetName()));
            }

            const std::string distrMInvVsPTBGName = 
               "c" + cName + "_z" + zName + "_r" + rName + "/" + 
                methodName + ": " + decayMode + "_BG12";

            TH2F *distrMInvVsPTBG = 
               static_cast<TH2F *>(inputFile->Get(distrMInvVsPTBGName.c_str()));

            if (!distrMInvVsPTBG)
            {
               CppTools::PrintError("Histogram named " + distrMInvVsPTBGName + 
                                    " does not exist in file " + 
                                    static_cast<std::string>(inputFile->GetName()));
            }

            const std::string distrMInvVsPTFGLRName = 
               "c" + cName + "_z" + zName + "_r" + rName + "/LR " + 
                methodName + ": " + decayMode + "_FG12";

            TH2F *distrMInvVsPTFGLR = 
               static_cast<TH2F *>(inputFile->Get(distrMInvVsPTFGLRName.c_str()));

            if (!distrMInvVsPTFGLR)
            {
               CppTools::PrintError("Histogram named " + distrMInvVsPTFGLRName + 
                                    " does not exist in file " +
                                    static_cast<std::string>(inputFile->GetName()));
            }

            const std::string distrMInvVsPTBGLRName = 
               "c" + cName + "_z" + zName + "_r" + rName + "/LR " + 
                methodName + ": " + decayMode + "_BG12";

            TH2F *distrMInvVsPTBGLR = 
               static_cast<TH2F *>(inputFile->Get(distrMInvVsPTBGLRName.c_str()));

            if (!distrMInvVsPTBGLR)
            {
               CppTools::PrintError("Histogram named " + distrMInvVsPTBGLRName + 
                                    " does not exist in file " + 
                                    static_cast<std::string>(inputFile->GetName()));
            }

            const int xAxisMin = distrMInvVsPTFG->GetXaxis()->FindBin(pTMin + 1e-6);
            const int xAxisMax = distrMInvVsPTFG->GetXaxis()->FindBin(pTMax - 1e-6);

            TH1D *distrMInvFG = 
               distrMInvVsPTFG->ProjectionY((distrMInvVsPTFG->GetName() + 
                                             std::to_string(c) + std::to_string(z) + 
                                             std::to_string(r)).c_str(), xAxisMin, xAxisMax);
            TH1D *distrMInvBG = 
               distrMInvVsPTBG->ProjectionY((distrMInvVsPTBG->GetName() + 
                                             std::to_string(c) + std::to_string(z) + 
                                             std::to_string(r)).c_str(), xAxisMin, xAxisMax);


            TH1D *distrMInvFGLR = 
               distrMInvVsPTFGLR->ProjectionY((distrMInvVsPTFGLR->GetName() + 
                                               std::to_string(c) + std::to_string(z) + 
                                               std::to_string(r)).c_str(), xAxisMin, xAxisMax);
            TH1D *distrMInvBGLR = 
               distrMInvVsPTBGLR->ProjectionY((distrMInvVsPTBGLR->GetName() + 
                                               std::to_string(c) + std::to_string(z) + 
                                               std::to_string(r)).c_str(), xAxisMin, xAxisMax);

            if (distrMInvFG->GetEntries() < 1e-3) continue;

            distrMInvFG->Sumw2();
            distrMInvBG->Sumw2();

            if (!distrMInvMerged) 
            {
               if (distrMInvBG->GetEntries() < 1e-3) 
               {
                  distrMInvMerged = static_cast<TH1D *>(distrMInvFG->Clone());
               }
               else distrMInvMerged = SubtractBG(distrMInvFG, distrMInvBG, 
                                                 distrMInvFGLR, distrMInvBGLR);

            }
            else 
            {
               if (distrMInvBG->GetEntries() < 1e-3) distrMInvMerged->Add(distrMInvFG);
               else distrMInvMerged->Add(SubtractBG(distrMInvFG, distrMInvBG, 
                                                    distrMInvFGLR, distrMInvBGLR));
            }
            if (!distrMInvMergedFG)
            {
               distrMInvMergedFG = static_cast<TH1D *>(distrMInvFG->Clone());
               distrMInvMergedBG = static_cast<TH1D *>(distrMInvBG->Clone());
               distrMInvMergedFGLR = static_cast<TH1D *>(distrMInvFGLR->Clone());
               distrMInvMergedBGLR = static_cast<TH1D *>(distrMInvBGLR->Clone());
            }
            else
            {
               distrMInvMergedFG->Add(distrMInvFG);
               distrMInvMergedBG->Add(distrMInvBG);
               distrMInvMergedFGLR->Add(distrMInvFGLR);
               distrMInvMergedBGLR->Add(distrMInvBGLR);
            }
         }
      }
   }

   return distrMInvMerged;
}

TH1D *MInv::SubtractBG(TH1D*& distrMInvFG, TH1D*& distrMInvBG, 
                       TH1D*& distrMInvFGLR, TH1D*& distrMInvBGLR)
{
   const double integralMInvFGLR = distrMInvFGLR->Integral(1, distrMInvFGLR->GetXaxis()->GetNbins());
   const double integralMInvBGLR = distrMInvBGLR->Integral(1, distrMInvBGLR->GetXaxis()->GetNbins());

   double partIntegralMInvFGLR = 0.;
   double partIntegralMInvBGLR = 0.;

   for (int i = distrMInvFGLR->GetXaxis()->GetNbins(); i >= 1; i--)
   {
      if (partIntegralMInvFGLR > integralMInvFGLR*0.1 && 
          partIntegralMInvBGLR > integralMInvBGLR*0.1) break;

      if (distrMInvBGLR->GetBinContent(i) < 1e-6) continue;

      partIntegralMInvFGLR += distrMInvFGLR->GetBinContent(i);
      partIntegralMInvBGLR += distrMInvBGLR->GetBinContent(i);
   }

   double scaleFactorBG = partIntegralMInvFGLR/partIntegralMInvBGLR;

   /*
   partIntegralMInvFGLR = 0.;
   partIntegralMInvBGLR = 0.;

   // if background was overestimated, some bins will be negative after BG subtraction
   // this needs to be resolved by rescaling (if needed)
   for (int i = 1; i <= distrMInvFGLR->GetXaxis()->GetNbins(); i++)
   {
      if (distrMInvBGLR->GetBinContent(i) < 1e-3) continue;
      if (distrMInvFGLR->GetBinContent(i) < 1e-3) continue;

      if (partIntegralMInvFGLR > integralMInvFGLR*0.1 && 
          partIntegralMInvBGLR > integralMInvBGLR*0.1) break;

      partIntegralMInvFGLR += distrMInvFGLR->GetBinContent(i);
      partIntegralMInvBGLR += distrMInvBGLR->GetBinContent(i);

      if ((distrMInvBGLR->GetBinContent(i) - distrMInvBGLR->GetBinError(i))*scaleFactorBG > 
          distrMInvFGLR->GetBinContent(i) + distrMInvFGLR->GetBinError(i))
      {
         scaleFactorBG *= distrMInvFGLR->GetBinContent(i)/distrMInvBGLR->GetBinContent(i);
      }
   }
   */

   // 0.95 is a rescale so that background is not overestimated
   distrMInvBG->Scale(scaleFactorBG*0.95);
   distrMInvBGLR->Scale(scaleFactorBG*0.95);

   TH1D *distrMInvSubtr = static_cast<TH1D *>(distrMInvFG->Clone());
   distrMInvSubtr->Add(distrMInvBG, -1.);

   return distrMInvSubtr;
}

#endif /* M_INV_CPP */
