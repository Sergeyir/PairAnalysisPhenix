#include "ErrorHandler.hpp"
#include "IOTools.hpp"
#include "StrTools.hpp"

#include "GUIDistrCutter2D.hpp"

void TimingDM()
{
	gStyle->SetOptStat(0);
   TDirectory::AddDirectory(kFALSE);

   const std::string heatmapIdentifierName = "Timing heatmap";

   CppTools::PrintInfo("List of runs in data/Real directory");
   system("ls data/Real/");

   CppTools::Print("Choose the run from the above ant type it in");
   std::string runName = "Run14HeAu200";
   //std::cout << ">> ";
   //std::cin >> runName;

   const std::string realInputFileName = "data/Real/" + runName + "/SingleTrack/sum.root";
   const std::string simInputFileName = "data/PostSim/" + runName + "/SingleTrack/all.root";
   CppTools::CheckInputFile(realInputFileName);
   CppTools::CheckInputFile(simInputFileName);
   TFile realInputFile(realInputFileName.c_str());
   TFile simInputFile(simInputFileName.c_str());

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

   TH2D averageTimeHist("average T in real data", "", 
                        realHist->GetYaxis()->GetNbins(), 
                        realHist->GetYaxis()->GetBinLowEdge(1), 
                        realHist->GetYaxis()->GetBinUpEdge(realHist->GetYaxis()->GetNbins()),
                        realHist->GetZaxis()->GetNbins(), 
                        realHist->GetZaxis()->GetBinLowEdge(1), 
                        realHist->GetZaxis()->GetBinUpEdge(realHist->GetZaxis()->GetNbins()));

   for (int i = 1; i <= realHist->GetYaxis()->GetNbins(); i++)
   {
      for (int j = 1; j <= realHist->GetZaxis()->GetNbins(); j++)
      {
         // distribution for a single smallest unit of a detector (tower, slat, strip, etc.)
         TH1D *distrTime = realHist->
            ProjectionX((realHist->GetName() + std::to_string(i) + std::to_string(j)).c_str(), 
                        i, i, j, j);

         //CppTools::Print(i, j, distrTime->Integral(1, distrTime->GetXaxis()->GetNbins());
         // removing units with insufficient statistics
         if (distrTime->Integral(1, distrTime->GetXaxis()->GetNbins()) < 10.)
         {
            //averageTimeHist.SetBinContent(i, j, -9999.);
            continue;
         }

         // usually the time is distributed by gausses of pions, kaons, and protons
         // by taking maximum value we extract the value close to the center of pions signal
         const double mean = distrTime->GetBinCenter(distrTime->GetMaximumBin());
         averageTimeHist.SetBinContent(i, j, mean);
      }
   }

   /*
   TH2D *simHist = static_cast<TH2D *>
      (simInputFile.Get((heatmapIdentifierName + ": " + detectorName).c_str()));

   if (!simHist) 
   {
      CppTools::PrintError("No histogram named" + heatmapIdentifierName + ": " + 
                           detectorName + " in file " + simInputFileName);
   }
   */

	TCanvas *canv = new TCanvas("", "", 900, 900);
   
   GUIDistrCutter2D::AddHistogram(&averageTimeHist);
   //GUIDistrCutter2D::AddHistogram(static_cast<TH2D *>(simHist->Clone("sim")));
   system(("mkdir -p data/Parameters/TimingDeadmaps/" + runName).c_str());

   while (detectorName.find(" ") < detectorName.size())
   {
      const unsigned int spacePos = detectorName.find(" ");
      detectorName.erase(spacePos, 1);
   }

   const std::string outputCutsFileName = "data/Parameters/TimingDeadmaps/" + 
                                          runName + "/TimingDeadmap" + detectorName + ".txt";

   if (CppTools::FileExists(outputCutsFileName))
   {
      GUIDistrCutter2D::ReadCutAreas(outputCutsFileName);
   }
   GUIDistrCutter2D::SetOutputFile(outputCutsFileName);

	gPad->AddExec("exec", "GUIDistrCutter2D::Exec()");
}
