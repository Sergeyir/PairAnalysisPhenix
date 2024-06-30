{
   /*
   gSystem->Setenv("LD_LIBRARY_PATH", 
                  (static_cast<std::string>(gSystem->Getenv("LD_LIBRARY_PATH")) + ":" + 
                   static_cast<std::string>(std::filesystem::current_path()) + 
                   "/../../ProgressBar/lib").c_str());
    */
   
   //system("echo $LD_LIBRARY_PATH");
   gInterpreter->AddIncludePath("../../ProgressBar/include");
   
   gSystem->Load("../../ProgressBar/lib/PBar.so");
}
