#include "../lib/ParInv.h"
#include "../lib/RABPainter.h"
#include "../lib/TCanvasPrinter.h"
#include "../lib/LogoDrawer.h"

void DrawRAB() 
{
   gROOT->SetBatch(kTRUE);
   gStyle->SetOptStat(kFALSE);

   const double ptmin = 0.9;
   const double ptmax = 6.5;

   TCanvas canv("canv", "canv", 1200, 800);
   canv.Divide(3, 2, 0., 0.);

   RABPainter single_rab = RABPainter("", 
      ptmin/1.1, ptmax*1.1, "p_{T} [GeV/c]", "R_{AA}", false, 0.05);
   
   system(("mkdir -p ../output/Tables/RAB/" + Par.run + "/" + Par.runnum).c_str());
   
   for (int i = 0; i < Par.CType.cmin.size(); i++)
   {
      if (i >= 2) canv.cd(i+2);
      else canv.cd(i+1);
      if (i < 2) gPad->SetTopMargin(0.016);
      else gPad->SetBottomMargin(0.14);

      switch (i)
      {
         case 0: 
            gPad->SetPad(0.01, 0.52, 0.35, 0.99);
            break;
         case 1: 
            gPad->SetPad(0.35, 0.52, 0.67, 0.99);
            break;
         case 2: 
            gPad->SetPad(0.01, 0., 0.35, 0.52);
            break;
         case 3: 
            gPad->SetPad(0.35, 0., 0.67, 0.52);
            break;
         case 4: gPad->SetPad(0.67, 0., 0.99, 0.52);
            break;
      }
      
      single_rab.AddGraph("../data/Spectra/" + Par.run + "/" + 
         Par.runnum + "/" + Par.particle.name_nl + "_" + 
         Par.CType.cname_nop[i] + ".txt",
         "../ext/Spectra/pp200/KStar892.txt",
         Par.particle.name, kRed+1, 53, Par.CType.ncolls[i], Par.pp_sigma,
         ptmin, ptmax);

      single_rab.PrintTable("../output/Tables/RAB/" + Par.run + "/" + 
         Par.runnum + "/KStar892_" + Par.CType.cname_nop[i] + ".tex");

      single_rab.AddScalingUncertainty(ptmax*1.05, 
         Par.CType.ncolls_uncertainty[i]/Par.CType.ncolls[i], kGray+2, 0.5, 0.1);
      
      switch (i)
      {
         case 0:
            single_rab.Draw(Par.CType.cname[i].c_str(), true, false, true, false);
            break;

         case 1:
            single_rab.Draw(Par.CType.cname[i].c_str(), true, false, false, false);
            break;

         case 2:
            single_rab.Draw(Par.CType.cname[i].c_str(), true, true, true, false);
            break;

         //case falltrhough
         case 3:
         case 4:
            single_rab.Draw(Par.CType.cname[i].c_str(), true, true, false, false);
            break;
      }
      
      DrawPHENIXLogoPreliminary(0.65, 0.68, 1.2);
   }

   canv.cd(3);
   gPad->SetPad(0.67, 0.52, 0.99, 0.99);
   gPad->SetTopMargin(0.016);
   single_rab.Draw(Par.CType.name, false, false, false, false);

   PrintCanvas(&canv, "../output/InvM/" + Par.run + "/" + 
      Par.runnum + "/SingleRAB_" + Par.particle.name_nl);

   RABPainter rab = RABPainter("", 
      ptmin/1.1, ptmax*1.1, "p_{T} [GeV/c]", "R_{AA}", false, 0.05);
   
   for (int i = 0; i < Par.CType.cmin.size(); i++)
   {
      if (i >= 2) canv.cd(i+2);
      else canv.cd(i+1);
      if (i < 2) gPad->SetTopMargin(0.016);
      else gPad->SetBottomMargin(0.14);

      switch (i)
      {
         case 0: 
            gPad->SetPad(0.01, 0.52, 0.35, 0.99);
            break;
         case 1: 
            gPad->SetPad(0.35, 0.52, 0.67, 0.99);
            break;
         case 2: 
            gPad->SetPad(0.01, 0., 0.35, 0.52);
            break;
         case 3: 
            gPad->SetPad(0.35, 0., 0.67, 0.52);
            break;
         case 4: gPad->SetPad(0.67, 0., 0.99, 0.52);
            break;
      }
      
      rab.AddGraph("../ext/RAB/AuAu200/phi1020_" + 
         Par.CType.cname_nop[i] + ".txt",   
         //"#varphi(1020), PPG096", kAzure+2, 54,
         "#varphi(1020), PRC72, 014903", kAzure+2, 54,
         ptmin, ptmax);

      rab.AddGraph("../ext/RAB/AuAu200/pi0_" + 
         Par.CType.cname_nop[i] + ".txt",   
         //"#pi^{0}, PPG014", kViolet-8, 55,
         "#pi^{0}, PRC91, 072301", kViolet-8, 55,
         ptmin, ptmax);

      rab.AddGraph("../ext/RAB/AuAu200/proton_" + 
         Par.CType.cname_nop[i] + ".txt",   
         //"p+#bar{p}, PPG146", kSpring+2, 59,
         "p+#bar{p}, PRC88, 024906", kSpring+2, 59,
         ptmin, ptmax);
      
      rab.AddGraph("../data/Spectra/" + Par.run + "/" + 
         Par.runnum + "/" + Par.particle.name_nl + "_" + 
         Par.CType.cname_nop[i] + ".txt",
         "../ext/Spectra/pp200/KStar892.txt",
         Par.particle.name, kRed+1, 53, Par.CType.ncolls[i], Par.pp_sigma,
         ptmin, ptmax);

      rab.PrintTable("../output/Tables/RAB/" + Par.run + "/" + 
         Par.runnum + "/KStar892_" + Par.CType.cname_nop[i] + ".tex");

      rab.AddScalingUncertainty(ptmax*1.05, 
         Par.CType.ncolls_uncertainty[i]/Par.CType.ncolls[i], kGray+2, 0.5, 0.1);
      
      switch (i)
      {
         case 0:
            rab.Draw(Par.CType.cname[i].c_str(), true, false, true);
            break;

         case 1:
            rab.Draw(Par.CType.cname[i].c_str(), true, false, false);
            break;

         case 2:
            rab.Draw(Par.CType.cname[i].c_str());
            break;

         //case falltrhough
         case 3:
         case 4:
            rab.Draw(Par.CType.cname[i].c_str(), true, true, false);
            break;
      }
      
      DrawPHENIXLogoPreliminary(0.65, 0.68, 1.2);
   }

   canv.cd(3);
   gPad->SetPad(0.67, 0.52, 0.99, 0.99);
   gPad->SetTopMargin(0.016);
   rab.Draw(Par.CType.name, false, false, false);

   PrintCanvas(&canv, "../output/InvM/" + Par.run + "/" + 
      Par.runnum + "/RAB_" + Par.particle.name_nl);

   RABPainter star_comp_rab = RABPainter("", 
      0.5, ptmax*1.1, "p_{T} [GeV/c]", "R_{AA}", false, 0.05);
   
   system(("mkdir -p ../output/Tables/RAB/" + Par.run + "/" + Par.runnum).c_str());

   TCanvas star_comp_canv("star_comp_canv", "canv", 600, 450);

   gPad->SetTopMargin(0.02);
   gPad->SetLeftMargin(0.1);
   gPad->SetBottomMargin(0.14);
   
   star_comp_rab.AddGraph("../data/Spectra/" + Par.run + "/" + 
      Par.runnum + "/" + Par.particle.name_nl + "_" + 
      Par.CType.cname_nop[0] + ".txt",
      "../ext/Spectra/pp200/KStar892.txt",
      Par.particle.name + " " + Par.CType.cname[0], kRed+1, 53, Par.CType.ncolls[0], Par.pp_sigma,
      0.5, ptmax);

   star_comp_rab.AddGraph("../ext/RAB/AuAu200/star_kstar892_0-10.txt",   
      "K*(892) 0-10\%, PRC71, 064902", kBlue+2, 65,
      0.5, ptmax, true, false);

   star_comp_rab.AddScalingUncertainty(ptmax*1.01, 
      Par.CType.ncolls_uncertainty[0]/Par.CType.ncolls[0], kRed+2, 0.5, 0.1);

   star_comp_rab.AddScalingUncertainty(ptmax*1.05, 
      0.097969, kBlue+2, 0.5, 0.1);

   star_comp_rab.Draw("", true, true, true);
      
   PrintCanvas(&star_comp_canv, "../output/InvM/" + Par.run + "/" + 
      Par.runnum + "/STARCompRAB_" + Par.particle.name_nl);
}   
