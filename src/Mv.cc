// $SOURCE$
//------------------------------------------------------------------------------------------------
//                              Mv macro
//------------------------------------------------------------------------------------------------
// Mv - move
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic macro for renaming TTree branches
 **/
//------------------------------------------------------------------------------------------------

#pragma once

#include "ErrorHandler.hpp"
#include "IOTools.hpp"

struct
{
   std::string run_name = "Run7AuAu200";
   std::string file_name = "akaon.root";
   std::string old_branch_name = "striptow";
   std::string new_branch_name = "striptofw";
} Par;

void Mv()
{
   CheckInputFile("../data/" + Par.run_name + "/" + Par.file_name);
   TFile file(("../data/" + Par.run_name + "/" + Par.file_name).c_str(), "UPDATE");
   TTree *tree = (TTree *) file.Get("Tree");
   TBranch *branch = (TBranch *) tree->GetBranch(Par.old_branch_name.c_str());
   branch->SetName(Par.new_branch_name.c_str());
   tree->Write();
   file.Close();
   PrintInfo("Name has changed and file has been written");
}
