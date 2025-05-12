#include "../lib/MS.h"
#include "../lib/ErrorHandler.h"
#include "../lib/StrTool.h"

void FitMS(std::string, std::string, std::string, TFile *, const double, const double, std::string, std::string = "");

int MS()
{
   ROOT::EnableImplicitMT();
   
   gROOT->SetBatch(kTRUE);
   gROOT->ProcessLine("gErrorIngonreLevel = 1001;");
   gErrorIgnoreLevel = kWarning;

   gStyle->SetOptStat(0);

   SetPartPar();

   std::string input_file_name = Par.run_dir + Par.run_name + "/" + Par.file_name +Par.magf + ".root";
   CheckInputFile(input_file_name);
   TFile *input = new TFile(input_file_name.c_str());

   TH1F *centr = (TH1F*) input->Get("central_bin");

   TH3F *mass3d_pos = (TH3F*) input->Get((Par.Det.hist_name + "_pos").c_str());
   TH3F *mass3d_neg = (TH3F*) input->Get((Par.Det.hist_name + "_neg").c_str());
   
   std::string output_dir = "../output/m2/" + Par.run_name;

   system(("mkdir -p " + output_dir + "/" + Par.Det.name + 
      "_" + Par.CType.cmin_name[Par.cnum] + "_" + Par.CType.cmax_name[Par.cnum]).c_str());

   const double nevents = centr->Integral(centr->FindBin(Par.CType.cmin[Par.cnum]*10.+0.01),
   centr->FindBin((Par.CType.cmax[Par.cnum]+1)*10-0.01));

   const double ptmin = Minimum(Par.pion.ptmin, Par.kaon.ptmin, Par.proton.ptmin);
   const double ptmax = Maximum(Par.pion.ptmax, Par.kaon.ptmax, Par.proton.ptmax);

   MSC *mass_fit = new MSC(ptmin, ptmax*Par.Det.ptmax_mult, -ptmax*Par.Det.ptmax_mult, -ptmin);

   for (int i = 0; i < Par.ptmin.size(); i++)
   {
      double pt = (Par.ptmin[i] + Par.ptmax[i])/2.;

      if (pt > ptmax) break;
      if (pt < ptmin) continue;
      
      std::string name_pos = "mass1d_pos" + DtoStr(pt, 2);
      std::string name_neg = "mass1d_neg" + DtoStr(pt, 2);

      TH1F *mass1d_pos = (TH1F*) mass3d_pos->ProjectionY(name_pos.c_str(), 
         mass3d_pos->GetXaxis()->FindBin(Par.ptmin[i]+0.01), 
         mass3d_pos->GetXaxis()->FindBin(Par.ptmax[i]-0.01), 
         mass3d_pos->GetZaxis()->FindBin(Par.CType.cmin[Par.cnum]*10.+0.1), 
         mass3d_pos->GetZaxis()->FindBin((Par.CType.cmax[Par.cnum]+1.)*10.-0.1));

      TH1F *mass1d_neg = (TH1F*) mass3d_neg->ProjectionY(name_neg.c_str(), 
         mass3d_neg->GetXaxis()->FindBin(Par.ptmin[i]+0.01), 
         mass3d_neg->GetXaxis()->FindBin(Par.ptmax[i]-0.01), 
         mass3d_neg->GetZaxis()->FindBin(Par.CType.cmin[Par.cnum]*10.+0.1),
         mass3d_neg->GetZaxis()->FindBin((Par.CType.cmax[Par.cnum]+1.)*10.-0.1));

      mass1d_pos->Rebin(Par.Det.rebin);
      mass1d_neg->Rebin(Par.Det.rebin);

      mass_fit->PerformM2Fit(mass1d_pos, mass1d_neg, Par.ptmin[i], Par.ptmax[i], nevents);
   }

   mass_fit->PerformMSFit();

   mass_fit->DrawMS(("../output/m2/" + Par.run_name + "/" + Par.Det.name + "_" + 
      Par.CType.cmin_name[Par.cnum] + "-" + Par.CType.cmax_name[Par.cnum]).c_str());

   mass_fit->PrintMSFunc();

   return 0;
}
