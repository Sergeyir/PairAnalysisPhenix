{
   gInterpreter->AddIncludePath("include");
   gInterpreter->AddIncludePath("CppTools/include");
   gInterpreter->AddIncludePath("ROOTTools/include");
   gInterpreter->AddIncludePath("ProgressBar/include");
   
   gSystem->Load("CppTools/lib/ErrorHandler.so");
   gSystem->Load("CppTools/lib/StrTools.so");
   gSystem->Load("CppTools/lib/IOTools.so");
   gSystem->Load("CppTools/lib/Box.so");
   gSystem->Load("CppTools/lib/Table.so");

   gSystem->Load("ROOTTools/lib/TCanvasPrinter.so");
   gSystem->Load("ROOTTools/lib/GUIFit.so");

   gSystem->Load("ProgressBar/lib/PBar.so");

   gSystem->Load("lib/DataCutsSelector.so");
}
