#include "ErrorHandler.hpp"
#include "IOTools.hpp"
#include "StrTools.hpp"

#include "GUIDistrCutter2D.hpp"

void DM()
{
	gStyle->SetOptStat(0);

   std::string runName;
   CppTools::Print("Type in the run name");
   std::cin >> runName;

   const std::string realInputFileName = "data/Real/" + runName + "/SingleTrack/sum.root";
   CppTools::CheckInputFile(realInputFileName);
   TFile realInputFile(realInputFileName.c_str());

   realInputFile.ls("Heatmap*");

   std::string detectorName;
   CppTools::Print("Type in the detector name");
   std::cin >> std::ws;
   std::getline(std::cin, detectorName);

   TH2D *realHist = static_cast<TH2D *>(realInputFile.Get(("Heatmap: " + detectorName).c_str()));
   if (!realHist) 
   {
      CppTools::PrintError("No histogram named Heatmap: " + detectorName + " in file " + realInputFileName);
   }

	TCanvas *canv = new TCanvas("", "", 900, 900);
   
   GUIDistrCutter2D::AddHistogram(realHist);
   system(("mkdir -p data/Deadmaps/" + runName).c_str());
   GUIDistrCutter2D::SetOutputFile("data/Deadmaps/" + runName + "/" + detectorName + ".txt");

   GUIDistrCutter2D::Exec();
	
	gPad->AddExec("exec", "GUIDistrCutter2D::Exec()");
}
