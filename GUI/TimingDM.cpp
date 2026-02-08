#include "ErrorHandler.hpp"
#include "IOTools.hpp"
#include "StrTools.hpp"

#include "PBar.hpp"

#include "GUIDistrCutter2D.hpp"

void TimingDM()
{
   gStyle->SetPalette(kWaterMelon);
	gStyle->SetOptStat(0);
   TDirectory::AddDirectory(kFALSE);
   gErrorIgnoreLevel = kWarning;
   ROOT::EnableImplicitMT();

   const std::string heatmapIdentifierName = "Timing heatmap";

   CppTools::PrintInfo("List of runs in data/Real directory");
   system("ls data/Real/");

   CppTools::Print("Choose the run from the above and type it in");
   std::string runName;
   std::cout << ">> ";
   std::cin >> runName;

   system(("mkdir -p data/Parameters/TimingDeadmaps/" + runName).c_str());
   system(("mkdir -p data/Parameters/TimingOffsets/" + runName).c_str());

   const std::string realInputFileName = "data/Real/" + runName + "/SingleTrack/sum.root";
   CppTools::CheckInputFile(realInputFileName);
   TFile realInputFile(realInputFileName.c_str());

   CppTools::PrintInfo("List of heatmaps in " + realInputFileName + " file");
   realInputFile.ls((heatmapIdentifierName + ":*").c_str());

   CppTools::Print("Type in the detector name");
   std::string detectorName;
   std::cout << ">> ";
   std::cin >> std::ws;
   std::getline(std::cin, detectorName);

   TH3D *realHist = static_cast<TH3D *>
      (realInputFile.Get((heatmapIdentifierName + ": " + detectorName).c_str()));

   if (!realHist) 
   {
      CppTools::PrintError("No histogram named" + heatmapIdentifierName + ": " + detectorName + 
                           " in file " + realInputFileName);
   }

   TH2D meanTimeHist("mean T in real data", "#mu", 
                     realHist->GetYaxis()->GetNbins(), 
                     realHist->GetYaxis()->GetBinLowEdge(1), 
                     realHist->GetYaxis()->GetBinUpEdge(realHist->GetYaxis()->GetNbins()),
                     realHist->GetZaxis()->GetNbins(), 
                     realHist->GetZaxis()->GetBinLowEdge(1), 
                     realHist->GetZaxis()->GetBinUpEdge(realHist->GetZaxis()->GetNbins()));

   TH2D sigmaTimeHist("sigma T in real data", "#sigma", 
                      realHist->GetYaxis()->GetNbins(), 
                      realHist->GetYaxis()->GetBinLowEdge(1), 
                      realHist->GetYaxis()->GetBinUpEdge(realHist->GetYaxis()->GetNbins()),
                      realHist->GetZaxis()->GetNbins(), 
                      realHist->GetZaxis()->GetBinLowEdge(1), 
                      realHist->GetZaxis()->GetBinUpEdge(realHist->GetZaxis()->GetNbins()));

   TH2D ratioFGToFullIntHist("FG to full integral ratio", "FG/(full integral)", 
                             realHist->GetYaxis()->GetNbins(), 
                             realHist->GetYaxis()->GetBinLowEdge(1), 
                             realHist->GetYaxis()->GetBinUpEdge(realHist->GetYaxis()->GetNbins()),
                             realHist->GetZaxis()->GetNbins(), 
                             realHist->GetZaxis()->GetBinLowEdge(1), 
                             realHist->GetZaxis()->GetBinUpEdge(realHist->GetZaxis()->GetNbins()));

   const std::string outputDir = "output/TimingDeadmaps/" + runName + "/" + detectorName;
   system(("mkdir -p " + outputDir).c_str());

   ProgressBar pBar("BLOCK", "Preparing heatmaps", PBarColor::BOLD_MAGENTA);

   const double numberOfIterations = static_cast<double>(realHist->GetYaxis()->GetNbins()*
                                                         realHist->GetZaxis()->GetNbins());

   gROOT->SetBatch(kTRUE);

   const std::string deadmapOutputFileName = "data/Parameters/TimingDeadmaps/" + 
                                             runName + "/TimingDeadmap" + detectorName + ".txt";
   const std::string offsetOutputFileName = "data/Parameters/TimingOffsets/" + 
                                             runName + "/TimingOffset" + detectorName + ".txt";

   bool deadmapOutputFileExists = CppTools::FileExists(deadmapOutputFileName);
   bool offsetOutputFileExists = CppTools::FileExists(offsetOutputFileName);

   // bins with low statistics should be considered dead areas thus they need to be cut
   // this file will contain these bins at first; the user will also add fiducial cuts
   // by applying them via GUI
   std::ofstream deadmapOutputFile;

   // means*(-1) i.e. offsets of pions signals will be written 
   // into separate file for later correction
   std::ofstream offsetOutputFile;

   if (!deadmapOutputFileExists)
   {
      deadmapOutputFile.open(deadmapOutputFileName);

      CppTools::PrintInfo("Low statistics bins will be automatically written in file " + 
                          deadmapOutputFileName);

      deadmapOutputFile << 
         realHist->GetYaxis()->GetNbins() << " " <<
         realHist->GetYaxis()->GetBinLowEdge(1) << " " <<
         realHist->GetYaxis()->GetBinUpEdge(realHist->GetYaxis()->GetNbins()) << " " <<
         realHist->GetZaxis()->GetNbins() << " " <<
         realHist->GetZaxis()->GetBinLowEdge(1) << " " <<
         realHist->GetZaxis()->GetBinUpEdge(realHist->GetZaxis()->GetNbins()) << std::endl;
   }
   else
   {
      CppTools::PrintInfo("File " + deadmapOutputFileName + " already exists:\n"\
                          " low statistics bins will not be automatically written");
   }

   if (!offsetOutputFileExists)
   {
      offsetOutputFile.open(offsetOutputFileName);

      CppTools::PrintInfo("Timing offset will be written in file " + offsetOutputFileName);

      offsetOutputFile << 
         realHist->GetYaxis()->GetNbins() << " " <<
         realHist->GetYaxis()->GetBinLowEdge(1) << " " <<
         realHist->GetYaxis()->GetBinUpEdge(realHist->GetYaxis()->GetNbins()) << " " <<
         realHist->GetZaxis()->GetNbins() << " " <<
         realHist->GetZaxis()->GetBinLowEdge(1) << " " <<
         realHist->GetZaxis()->GetBinUpEdge(realHist->GetZaxis()->GetNbins()) << std::endl;
   }
   else
   {
      CppTools::PrintInfo("File " + offsetOutputFileName + " already exists:\n"\
                          " timing offset will not be written");
   }

   TCanvas fitCanv("fit canv", "", 400, 400);

   int numberOfBinsWithLowStat = 0;

   for (int i = 1; i <= realHist->GetZaxis()->GetNbins(); i++)
   {
      for (int j = 1; j <= realHist->GetYaxis()->GetNbins(); j++)
      {
         pBar.Print(static_cast<double>(i + j*i - 2)/numberOfIterations);

         // distribution for a single smallest unit of a detector (tower, slat, strip, etc.)
         TH1D *distrTime = realHist->
            ProjectionX((realHist->GetName() + std::to_string(i) + std::to_string(j)).c_str(), 
                        j, j, i, i);

         const double fullIntegral = distrTime->Integral(1, distrTime->GetXaxis()->GetNbins());
         // removing units with insufficient statistics
         if (fullIntegral < 100.)
         {
            numberOfBinsWithLowStat++;
            if (!deadmapOutputFileExists)
            {
               deadmapOutputFile << 1;
               if (j < realHist->GetYaxis()->GetNbins()) deadmapOutputFile << " ";
               else deadmapOutputFile << std::endl;
            }
            if (!offsetOutputFileExists)
            {
               offsetOutputFile << 0.;
               if (j < realHist->GetYaxis()->GetNbins()) offsetOutputFile << " ";
               else offsetOutputFile << std::endl;
            }
            continue;
         }

         // approximation of t-t_exp for pi+
         TF1 fitFG("fg fit", "gaus(0) + gaus(3) + [6]");
         TF1 fitBG("bg fit", "gaus(0) + [3]");

         fitBG.SetLineColor(kBlue);
         fitBG.SetLineStyle(2);

         // usually the time is distributed by gausses of pions, kaons, and protons
         // by taking maximum value we extract the value close to the center of pions signal
         const int maximumBin = distrTime->GetMaximumBin();

         fitFG.SetParameters(distrTime->GetBinContent(maximumBin), 
                             distrTime->GetXaxis()->GetBinCenter(maximumBin), 
                             distrTime->GetXaxis()->GetBinWidth(1)*2., 
                             distrTime->GetBinContent(maximumBin + 10),
                             distrTime->GetXaxis()->GetBinLowEdge(maximumBin + 6),
                             distrTime->GetXaxis()->GetBinWidth(1)*5.);

         fitFG.SetParLimits(0, distrTime->GetBinContent(maximumBin)/3., 
                            distrTime->GetBinContent(maximumBin));
         fitFG.SetParLimits(1, distrTime->GetXaxis()->GetBinLowEdge(maximumBin - 1), 
                            distrTime->GetXaxis()->GetBinUpEdge(maximumBin));
         fitFG.SetParLimits(3, distrTime->GetBinContent(maximumBin + 10)/2., 
                            distrTime->GetBinContent(maximumBin)/2.);
         fitFG.SetParLimits(4, distrTime->GetXaxis()->GetBinLowEdge(maximumBin + 2),
                            distrTime->GetXaxis()->GetBinCenter(maximumBin + 10));

         fitFG.SetRange(distrTime->GetXaxis()->GetBinLowEdge(maximumBin - 10), 
                        distrTime->GetXaxis()->GetBinUpEdge(maximumBin + 10));
         fitBG.SetRange(distrTime->GetXaxis()->GetBinLowEdge(maximumBin - 10), 
                        distrTime->GetXaxis()->GetBinUpEdge(maximumBin + 10));

         distrTime->Fit(&fitFG, "RQMBN");
         meanTimeHist.SetBinContent(j, i, fitFG.GetParameter(1));
         sigmaTimeHist.SetBinContent(j, i, fabs(fitFG.GetParameter(2)));

         for (int k = 0; k < fitBG.GetNpar(); k++)
         {
            fitBG.SetParameter(k, fitFG.GetParameter(k + 3));
         }

         double integralFG = 0.;
         for (int k = maximumBin - 5; k <= maximumBin + 5; k++)
         {
            integralFG += fitFG.Eval(distrTime->GetXaxis()->GetBinCenter(k)) - 
                          fitBG.Eval(distrTime->GetXaxis()->GetBinCenter(k)); 
         }

         ratioFGToFullIntHist.SetBinContent(j, i, integralFG/fullIntegral);

         distrTime->SetLineColor(kBlack);

         gPad->Clear();
         distrTime->Draw();
         fitFG.Draw("SAME");
         fitBG.Draw("SAME");
         
         fitCanv.SaveAs((outputDir + "/fit_" + std::to_string(i) + "_" + 
                         std::to_string(j) + ".png").c_str());

         if (!deadmapOutputFileExists)
         {
            deadmapOutputFile << 0;
            if (j < realHist->GetYaxis()->GetNbins()) deadmapOutputFile << " ";
            else deadmapOutputFile << std::endl;
         }
         if (!offsetOutputFileExists)
         {
            offsetOutputFile << -1.*fitFG.GetParameter(1);
            if (j < realHist->GetYaxis()->GetNbins()) offsetOutputFile << " ";
            else offsetOutputFile << std::endl;
         }
      }
   }

   CppTools::PrintInfo("There were " + std::to_string(numberOfBinsWithLowStat) + 
                       " bins with low statistics");

   pBar.Clear();
   CppTools::PrintInfo("Preparing heatmaps: done");

   gROOT->SetBatch(kFALSE);

	TCanvas *canv = new TCanvas("", "", 1080, 1080);
   
   GUIDistrCutter2D::AddHistogram(&meanTimeHist);
   GUIDistrCutter2D::AddHistogram(&sigmaTimeHist);
   GUIDistrCutter2D::AddHistogram(&ratioFGToFullIntHist);

   while (detectorName.find(" ") < detectorName.size())
   {
      const unsigned int spacePos = detectorName.find(" ");
      detectorName.erase(spacePos, 1);
   }

   if (CppTools::FileExists(deadmapOutputFileName))
   {
      GUIDistrCutter2D::ReadCutAreas(deadmapOutputFileName);
   }
   GUIDistrCutter2D::SetOutputFile(deadmapOutputFileName);

	gPad->AddExec("exec", "GUIDistrCutter2D::Exec()");
}
