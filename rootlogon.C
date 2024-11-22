{
   gInterpreter->AddIncludePath("Analysis/include");
   gInterpreter->AddIncludePath("CppTools/include");
   gInterpreter->AddIncludePath("ROOTTools/include");
   gInterpreter->AddIncludePath("ProgressBar/include");
   
   gSystem->Load("CppTools/lib/ErrorHandler.so");
   gSystem->Load("CppTools/lib/StrTools.so");
   gSystem->Load("CppTools/lib/IOTools.so");
   gSystem->Load("CppTools/lib/Box.so");
   gSystem->Load("CppTools/lib/Table.so");

   gSystem->Load("Analysis/lib/DeadAreasCuts.so");
   
   gSystem->Load("ProgressBar/lib/PBar.so");
}
