#include "ErrorHandler.hpp"
#include "IOTools.hpp"
#include "StrTools.hpp"

#include "GUIDistrCutter2D.hpp"

void DM()
{
	gStyle->SetOptStat(0);
   TDirectory::AddDirectory(kFALSE);

   const std::string heatmapIdentifierName = "Heatmap";

   CppTools::PrintInfo("List of runs in data/Real directory");
   system("ls data/Real/");

   CppTools::Print("Choose the run from the above and type it in");
   std::string runName;
   std::cout << ">> ";
   std::cin >> runName;

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

   TH2D *realHist = static_cast<TH2D *>
      (realInputFile.Get((heatmapIdentifierName + ": " + detectorName).c_str()));

   if (!realHist) 
   {
      CppTools::PrintError("No histogram named \"" + heatmapIdentifierName + ": " + detectorName + 
                           "\" in file " + realInputFileName);
   }

   TH2D *simHist = static_cast<TH2D *>
      (simInputFile.Get((heatmapIdentifierName + ": " + detectorName).c_str()));

   if (!simHist) 
   {
      CppTools::PrintError("No histogram named" + heatmapIdentifierName + ": " + 
                           detectorName + " in file " + simInputFileName);
   }

   //simHist->Scale(realHist->Integral()/simHist->Integral());

	TCanvas *canv = new TCanvas("", "", 900, 900);
   
   GUIDistrCutter2D::AddHistogram(realHist);
   GUIDistrCutter2D::AddHistogram(static_cast<TH2D *>(simHist->Clone("sim")));
   system(("mkdir -p data/Parameters/Deadmaps/" + runName).c_str());

   while (detectorName.find(" ") < detectorName.size())
   {
      const unsigned int spacePos = detectorName.find(" ");
      detectorName.erase(spacePos, 1);
   }

   const std::string outputCutsFileName = "data/Parameters/Deadmaps/" + 
                                          runName + "/" + detectorName + ".txt";

   if (CppTools::FileExists(outputCutsFileName))
   {
      GUIDistrCutter2D::ReadCutAreas(outputCutsFileName);
   }
   GUIDistrCutter2D::SetOutputFile(outputCutsFileName);

	gPad->AddExec("exec", "GUIDistrCutter2D::Exec()");
}
