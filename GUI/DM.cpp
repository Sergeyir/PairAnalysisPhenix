#include "ErrorHandler.hpp"
#include "IOTools.hpp"
#include "StrTools.hpp"

#include "GUIDistrCutter2D.hpp"

void DM()
{
	gStyle->SetOptStat(0);
   TDirectory::AddDirectory(kFALSE);

   CppTools::PrintInfo("List of runs in data/Real directory");
   system("ls data/Real/");

   CppTools::Print("Choose the run from the above ant type it in");
   std::string runName = "Run14HeAu200";
   //std::cout << ">> ";
   //std::cin >> runName;

   const std::string realInputFileName = "data/Real/" + runName + "/SingleTrack/sum.root";
   CppTools::CheckInputFile(realInputFileName);
   TFile realInputFile(realInputFileName.c_str());

   CppTools::PrintInfo("List of heatmaps in " + realInputFileName + " file");
   realInputFile.ls("Heatmap:*");

   CppTools::Print("Type in the detector name");
   std::string detectorName;
   std::cout << ">> ";
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

   const std::string outputCutsFileName = "data/Deadmaps/" + runName + "/" + detectorName + ".txt";
   if (CppTools::FileExists(outputCutsFileName))
   {
      GUIDistrCutter2D::ReadCutAreas(outputCutsFileName);
   }
   GUIDistrCutter2D::SetOutputFile(outputCutsFileName);

	gPad->AddExec("exec", "GUIDistrCutter2D::Exec()");
}
