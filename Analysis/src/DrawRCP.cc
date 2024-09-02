#include "../lib/ParInv.h"
#include "../lib/RCPPainter.h"
#include "../lib/TCanvasPrinter.h"
#include "../lib/LogoDrawer.h"

void DrawRCP() 
{
   gROOT->SetBatch(kTRUE);
   gStyle->SetOptStat(kFALSE);

   const double ptmin = 0.5;
   const double ptmax = 6.5;

   TCanvas canv("canv", "canv", 1200, 600);

   canv.Divide(2);

   RCPPainter single_rcp = RCPPainter("", 
      ptmin/1.1, ptmax*1.1, "p_{T} [GeV/c]", "R_{CP}", false, 0.05);

   system(("mkdir -p ../output/Tables/RCP/" + Par.run + "/" + Par.runnum).c_str());
   
   single_rcp.AddGraph("../data/Spectra/" + Par.run + "/" + 
      Par.runnum + "/" + Par.particle.name_nl + "_" + 
      Par.CType.cname_nop[Par.CType.central] + ".txt",
      "../data/Spectra/" + Par.run + "/" + 
      Par.runnum + "/" + Par.particle.name_nl + "_" + 
      Par.CType.cname_nop[Par.CType.peripheral-1] + ".txt",
      Par.particle.name, kRed-3, 53, 
      Par.CType.ncolls[Par.CType.central], Par.CType.ncolls[Par.CType.peripheral-1],
      ptmin, ptmax);

   single_rcp.PrintTable("../output/Tables/RCP/" + Par.run + "/" + Par.runnum + "/KStar892_" + 
      Par.CType.cname_nop[Par.CType.central] + "_" + Par.CType.cname_nop[Par.CType.peripheral-1] + ".tex");

   single_rcp.AddScalingUncertainty(ptmax*1.01, ErrPropagation(
      Par.CType.ncolls_uncertainty[Par.CType.central]/
      Par.CType.ncolls[Par.CType.central], 
      Par.CType.ncolls_uncertainty[Par.CType.peripheral-1]/
      Par.CType.ncolls[Par.CType.peripheral-1]), 
      kGray+2, 0.5, 0.08);
   
   canv.cd(1);
   gPad->SetLeftMargin(0.16);
   gPad->SetBottomMargin(0.14);
   
   single_rcp.Draw(Par.CType.cname_nop[Par.CType.central] + "/" + 
      Par.CType.cname_nop[Par.CType.peripheral-1], false);

   DrawPHENIXLogoPreliminary(0.62, 0.57, 1.2);

   single_rcp.AddGraph("../data/Spectra/" + Par.run + "/" + 
      Par.runnum + "/" + Par.particle.name_nl + "_" + 
      Par.CType.cname_nop[Par.CType.central] + ".txt",
      "../data/Spectra/" + Par.run + "/" + 
      Par.runnum + "/" + Par.particle.name_nl + "_" + 
      Par.CType.cname_nop[Par.CType.peripheral] + ".txt",
      Par.particle.name, kRed-3, 53, 
      Par.CType.ncolls[Par.CType.central], Par.CType.ncolls[Par.CType.peripheral],
      ptmin, ptmax);

   single_rcp.PrintTable("../output/Tables/RCP/" + Par.run + "/" + Par.runnum + "/KStar892_" + 
      Par.CType.cname_nop[Par.CType.central] + "_" + Par.CType.cname_nop[Par.CType.peripheral] + ".tex");

   single_rcp.AddScalingUncertainty(ptmax*1.01, ErrPropagation(
      Par.CType.ncolls_uncertainty[Par.CType.central]/
      Par.CType.ncolls[Par.CType.central], 
      Par.CType.ncolls_uncertainty[Par.CType.peripheral]/
      Par.CType.ncolls[Par.CType.peripheral]), 
      kGray+2, 0.5, 0.08);

   canv.cd(2);
   gPad->SetLeftMargin(0.16);
   gPad->SetBottomMargin(0.14);
   
   single_rcp.Draw(Par.CType.cname_nop[Par.CType.central] + "/" + 
      Par.CType.cname_nop[Par.CType.peripheral], false);

   DrawPHENIXLogoPreliminary(0.62, 0.57, 1.2);

   PrintCanvas(&canv, "../output/InvM/" + Par.run + "/" + 
      Par.runnum + "/SingleRCP_" + Par.particle.name_nl);

   RCPPainter rcp = RCPPainter("", 
      ptmin/1.1, ptmax*1.1, "p_{T} [GeV/c]", "R_{CP}", false, 0.05);

   rcp.AddGraph("../ext/RCP/AuAu200/pion_0010_4060.txt",   
      //"0-10/40-60 #pi^{+}+#pi^{-} (PPG146)", kViolet-8, 55,
      "0-10/40-60 #pi^{+}+#pi^{-}, PRC88, 024906", kViolet-8, 55,
      ptmin, ptmax);

   /*
   rcp.AddGraph("../ext/RCP/AuAu200/kaon_0010_4060.txt",   
      "0-10/40-60 K^{+}+K^{-} (PPG146)", kGray+3, 56,
      ptmin, ptmax);
   */

   rcp.AddGraph("../ext/RCP/AuAu200/proton_0010_4060.txt",   
      //"0-10/40-60 p+#bar{p} (PPG146)", kSpring+5, 59,
      "0-10/40-60 p+#bar{p}, PRC88, 024906", kSpring+5, 59,
      ptmin, ptmax);

   rcp.AddGraph("../ext/Spectra/AuAu200/phi1020_0-10.txt",
      "../ext/Spectra/AuAu200/phi1020_40-60.txt",
      //"0-10/40-60 #varphi(1020) (PPG096)", kAzure+2, 55, 
      "0-10/40-60 #varphi(1020), PRC72, 014903", kAzure+2, 55, 
      955.4, Par.CType.ncolls[Par.CType.peripheral-1],
      ptmin, ptmax);

   rcp.AddGraph("../data/Spectra/" + Par.run + "/" + 
      Par.runnum + "/" + Par.particle.name_nl + "_" + 
      Par.CType.cname_nop[Par.CType.central] + ".txt",
      "../data/Spectra/" + Par.run + "/" + 
      Par.runnum + "/" + Par.particle.name_nl + "_" + 
      Par.CType.cname_nop[Par.CType.peripheral-1] + ".txt",
      Par.particle.name, kRed-3, 53, 
      Par.CType.ncolls[Par.CType.central], Par.CType.ncolls[Par.CType.peripheral-1],
      ptmin, ptmax);

   rcp.AddScalingUncertainty(ptmax*1.01, ErrPropagation(
      Par.CType.ncolls_uncertainty[Par.CType.central]/
      Par.CType.ncolls[Par.CType.central], 
      Par.CType.ncolls_uncertainty[Par.CType.peripheral-1]/
      Par.CType.ncolls[Par.CType.peripheral-1]), 
      kRed-6, 0.8, 0.08);

   rcp.AddScalingUncertainty(ptmax*1.01+0.16, 0.169568, 
      kAzure-7, 0.8, 0.08);
   
   canv.cd(1);
   gPad->SetLeftMargin(0.16);
   gPad->SetBottomMargin(0.14);
   
   rcp.Draw(Par.CType.cname_nop[Par.CType.central] + "/" + 
      Par.CType.cname_nop[Par.CType.peripheral-1]);

   DrawPHENIXLogoPreliminary(0.62, 0.57, 1.2);

   rcp.AddGraph("../ext/RCP/AuAu200/pion_0010_6093.txt",   
      //"0-10/60-93 #pi^{+}+#pi^{-} (PPG146)", kViolet-8, 55,
      "0-10/60-93 #pi^{+}+#pi^{-}, PRC88, 024906", kViolet-8, 55,
      ptmin, ptmax);

   /*
   rcp.AddGraph("../ext/RCP/AuAu200/kaon_0010_6093.txt",   
      "0-10/60-93 K^{+}+K^{-} (PPG146)", kGray+3, 56,
      ptmin, ptmax);
      */
   
   rcp.AddGraph("../ext/RCP/AuAu200/proton_0010_6093.txt",   
      //"0-10/60-93 p+#bar{p} (PPG146)", kSpring+5, 59,
      "0-10/60-93 p+#bar{p}, PRC88, 024906", kSpring+5, 59,
      ptmin, ptmax);

   rcp.AddGraph("../ext/Spectra/AuAu200/phi1020_0-10.txt",
      "../ext/Spectra/AuAu200/phi1020_60-93.txt",
      //"0-10/60-93 #varphi(1020) (PPG096)", kAzure+2, 55, 
      "0-10/60-93 #varphi(1020), PRC72, 014903", kAzure+2, 55, 
      955.4, Par.CType.ncolls[Par.CType.peripheral],
      ptmin, ptmax);

   rcp.AddGraph("../data/Spectra/" + Par.run + "/" + 
      Par.runnum + "/" + Par.particle.name_nl + "_" + 
      Par.CType.cname_nop[Par.CType.central] + ".txt",
      "../data/Spectra/" + Par.run + "/" + 
      Par.runnum + "/" + Par.particle.name_nl + "_" + 
      Par.CType.cname_nop[Par.CType.peripheral] + ".txt",
      Par.particle.name, kRed-3, 53, 
      Par.CType.ncolls[Par.CType.central], Par.CType.ncolls[Par.CType.peripheral],
      ptmin, ptmax);

   rcp.AddScalingUncertainty(ptmax*1.01, ErrPropagation(
      Par.CType.ncolls_uncertainty[Par.CType.central]/
      Par.CType.ncolls[Par.CType.central], 
      Par.CType.ncolls_uncertainty[Par.CType.peripheral]/
      Par.CType.ncolls[Par.CType.peripheral]), 
      kRed-6, 0.8, 0.08);

   rcp.AddScalingUncertainty(ptmax*1.01+0.16, 0.2927, 
      kAzure-7, 0.8, 0.08);
   
   canv.cd(2);
   gPad->SetLeftMargin(0.16);
   gPad->SetBottomMargin(0.14);
   
   rcp.Draw(Par.CType.cname_nop[Par.CType.central] + "/" + 
      Par.CType.cname_nop[Par.CType.peripheral]);

   DrawPHENIXLogoPreliminary(0.62, 0.57, 1.2);

   PrintCanvas(&canv, "../output/InvM/" + Par.run + "/" + 
      Par.runnum + "/RCP_" + Par.particle.name_nl);
}   
