{
   gInterpreter->AddIncludePath("include");
   gInterpreter->AddIncludePath("CppTools/include");
   gInterpreter->AddIncludePath("ROOTTools/include");
   gInterpreter->AddIncludePath("ProgressBar/include");
   gInterpreter->AddIncludePath("yaml-cpp/include");
   
   gSystem->Load("lib/libErrorHandler.so");
   gSystem->Load("lib/libStrTools.so");
   gSystem->Load("lib/libIOTools.so");
   gSystem->Load("lib/libBox.so");
   gSystem->Load("lib/libTable.so");

   gSystem->Load("lib/libTCanvasTools.so");
   gSystem->Load("lib/libGUIDistrCutter2D.so");

   gSystem->Load("lib/libPBar.so");

   gSystem->Load("lib/libDeadMapCutter.so");
   gSystem->Load("lib/libFitFunc.so");

   //gSystem->Load("yaml-cpp/build/libyaml-cpp.so");

   //gSystem->Load("libjsoncpp.so");
}
