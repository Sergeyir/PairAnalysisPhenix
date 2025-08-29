/** 
 *  @file   M2IdentFit.cpp 
 *  @brief  Contains realistation of functions that are used for approximation of m2 distributions for the estimation of timing parameters that are used for identification via m2
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/CalPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef M2_IDENT_FIT_CPP
#define M2_IDENT_FIT_CPP

#include "M2IdentFit.hpp"

// this namespace is only used so that documentation will not become a mess
// so there is no need to enforce the contents inside of it 
// being accessed only via the scope resolution operator in this file
using namespace M2IdentFit;

int main(int argc, char **argv)
{
   if (argc < 2 || argc > 3) 
   {
      std::string errMsg = "Expected 1-2 parameters while " + std::to_string(argc) + " ";
      errMsg += "parameter(s) were provided \n Usage: bin/M2IdentFit \
                 inputYAMLName detectorName=all* \n\
                 * you can specify the name of a detector you want for this program to run for or \
                 to leave it empty for the program to run over all detectors that are specified in \
                 the input file";
      CppTools::PrintError(errMsg);
   }
 
   CppTools::CheckInputFile(argv[1]);
 
   ROOT::EnableImplicitMT();

   inputYAMLM2Id.OpenFile(argv[1], "m2id");
   inputYAMLM2Id.CheckStatus("m2id");

   runName = inputYAMLM2Id["run_name"].as<std::string>();

   inputYAMLMain.OpenFile("input/" + runName + "/main.yaml");
   inputYAMLMain.CheckStatus("main");
 
   collisionSystemName = inputYAMLMain["collision_system_name"].as<std::string>();

   const std::string inputDataFileName = 
      "data/Real/" + runName + "/SingleTrack/sum.root";
   CppTools::CheckInputFile(inputDataFileName);
   inputDataFile = TFile::Open(inputDataFileName.c_str(), "READ");

   K1 = inputYAMLM2Id["K1"].as<double>();

   outputDir = "output/M2Id/" + runName;

   parametersDir = "data/Parameters/M2Id/" + runName;
   system(("mkdir -p " + parametersDir).c_str());

   rawYieldsDir = "data/RawYields/SingleTrack/" + runName;
   system(("mkdir -p " + rawYieldsDir).c_str());

   gErrorIgnoreLevel = kWarning;
   gStyle->SetOptStat(0);

   // calculating the number of iterations needed for this program
   // this step is needed by progress bar
   for (const YAML::Node& detector : inputYAMLM2Id["detectors"])
   {
      if (argc == 3 && detector["name"].as<std::string>() != 
          static_cast<std::string>(argv[2])) continue;

      TH3F* m2DistrPos = static_cast<TH3F *>
         (inputDataFile->Get(("m2, " + detector["name"].as<std::string>() + ", charge>0").c_str()));

      for (const YAML::Node& pTBin : inputYAMLM2Id["pt_bins"])
      {
         if (pTBin["min"].as<double>() + 1e-3 < m2DistrPos->GetXaxis()->GetBinLowEdge(1)) continue;
         if (pTBin["max"].as<double>() - 1e-3 > m2DistrPos->GetXaxis()->
             GetBinUpEdge(m2DistrPos->GetXaxis()->GetNbins())) break;
         numberOfIterations++;
      }
   }
   numberOfIterations++;

   for (const YAML::Node& detector : inputYAMLM2Id["detectors"])
   {
      if (argc == 3 && detector["name"].as<std::string>() != 
          static_cast<std::string>(argv[2])) continue;
      system(("mkdir -p " + outputDir + "/" + detector["name"].as<std::string>()).c_str());
      PerformFitsForDetector(detector, inputYAMLMain["centrality_min"].as<double>(), 
                             inputYAMLMain["centrality_max"].as<double>());
      pBar.HandleOutput(" " + detector["name"].as<std::string>() + " done");
   }
   pBar.Clear();
   CppTools::PrintInfo("M2IdentFit has finished running succesfully");
}

void M2IdentFit::PerformFitsForDetector(const YAML::Node& detector, 
                                        const double centralityMin, 
                                        const double centralityMax)
{
   TH3F* m2DistrPos = static_cast<TH3F *>
      (inputDataFile->Get(("m2, " + detector["name"].as<std::string>() + ", charge>0").c_str()));
   TH3F* m2DistrNeg = static_cast<TH3F *>
      (inputDataFile->Get(("m2, " + detector["name"].as<std::string>() + ", charge<0").c_str()));

   // minimum pT in whole pT range
   const double pTMin = CppTools::Minimum(detector["pi+_pt_bounds"][0].as<double>(), 
                                          detector["k+_pt_bounds"][0].as<double>(),
                                          detector["p_pt_bounds"][0].as<double>(),
                                          detector["pi-_pt_bounds"][0].as<double>(), 
                                          detector["k-_pt_bounds"][0].as<double>(),
                                          detector["pbar_pt_bounds"][0].as<double>());
   // maximum pT in whole pT range
   const double pTMax = CppTools::Maximum(detector["pi+_pt_bounds"][1].as<double>(), 
                                          detector["k+_pt_bounds"][1].as<double>(),
                                          detector["p_pt_bounds"][1].as<double>(),
                                          detector["pi-_pt_bounds"][1].as<double>(), 
                                          detector["k-_pt_bounds"][1].as<double>(),
                                          detector["pbar_pt_bounds"][1].as<double>());

   // container that holds fit parameters and yields for pi^+
   FitParameters fitPiPlus("pi+", pow(0.139570, 2), detector, true, kRed);
   // container that holds fit parameters and yields for K^+
   FitParameters fitKPlus("K+", pow(0.493677, 2), detector, true, kGreen);
   // container that holds fit parameters and yields for protons
   FitParameters fitP("p", pow(0.938272, 2), detector, true, kAzure);
   // container that holds fit parameters and yields for pi^-
   FitParameters fitPiMinus("pi-", pow(0.139570, 2), detector, false, kRed);
   // container that holds fit parameters and yields for K^-
   FitParameters fitKMinus("K-", pow(0.493677, 2), detector, false, kGreen);
   // container that holds fit parameters and yields for antiprotons
   FitParameters fitPBar("pbar", pow(0.938272, 2), detector, false, kAzure);

   for (const YAML::Node& pTBin : inputYAMLM2Id["pt_bins"])
   {
      pBar.Print(static_cast<double>(numberOfCalls)/static_cast<double>(numberOfIterations));

      // minium pT for the current pT bin
      const double binPTMin = pTBin["min"].as<double>();
      // maximum pT for the current pT bin
      const double binPTMax = pTBin["max"].as<double>();

      const double pT = (binPTMin + binPTMax)/2.;

      if (binPTMin + 1e-3 < m2DistrPos->GetXaxis()->GetBinLowEdge(1)) continue;
      if (binPTMax - 1e-3 > m2DistrPos->GetXaxis()->
          GetBinUpEdge(m2DistrPos->GetXaxis()->GetNbins())) break;

      const std::string pTRangeName = 
         CppTools::DtoStr(pTBin["min"].as<double>(), 1) + " < p_{T} < " +
         CppTools::DtoStr(pTBin["max"].as<double>(), 1) + " [GeV/c]";

      TH1D *m2DistrPosProj = m2DistrPos->
         ProjectionY((m2DistrPos->GetName() + std::to_string((binPTMin + binPTMax)/2.)).c_str(),
                     m2DistrPos->GetXaxis()->FindBin(binPTMin + 1e-15),
                     m2DistrPos->GetXaxis()->FindBin(binPTMax - 1e-15),
                     m2DistrPos->GetXaxis()->FindBin(centralityMin + 1e-15),
                     m2DistrPos->GetXaxis()->FindBin(centralityMax - 1e-15));

      if (m2DistrPosProj->Integral(1, m2DistrPosProj->GetXaxis()->GetNbins()) < 1.)
      {
         CppTools::PrintWarning("Histogram for positive tracks is empty for " + pTRangeName +
                                " in " + detector["name"].as<std::string>());
      }

      TH1D *m2DistrNegProj = m2DistrNeg->
         ProjectionY((m2DistrNeg->GetName() + std::to_string((binPTMin + binPTMax)/2.)).c_str(),
                     m2DistrNeg->GetXaxis()->FindBin(binPTMin + 1e-15),
                     m2DistrNeg->GetXaxis()->FindBin(binPTMax - 1e-15),
                     m2DistrNeg->GetXaxis()->FindBin(centralityMin + 1e-15),
                     m2DistrNeg->GetXaxis()->FindBin(centralityMax - 1e-15));

      if (m2DistrNegProj->Integral(1, m2DistrPosProj->GetXaxis()->GetNbins()) < 1.)
      {
         CppTools::PrintWarning("Histogram for negative tracks is empty for " + pTRangeName +
                                " in " + detector["name"].as<std::string>());
      }

      TCanvas m2SingleFitCanv(("m2 pos and neg fits" + pTRangeName).c_str(), "", 900, 900);
		m2SingleFitCanv.Divide(1, 2);

		TText text;

		text.SetTextFont(43);
		text.SetTextSize(40);
		text.SetTextAngle(90);
	
		text.DrawText(0.04, 0.65, "charge = +1");
		text.DrawText(0.04, 0.15, "charge = -1");

		TLatex tlText;

		tlText.SetTextFont(52);
		tlText.SetTextSize(0.05);

		tlText.DrawLatex(0.45, 0.9, pTRangeName.c_str());
		tlText.DrawLatex(0.45, 0.4, pTRangeName.c_str());

      m2DistrPosProj->GetXaxis()->SetRange(m2DistrPosProj->GetXaxis()->FindBin(-0.4), 
                                           m2DistrPosProj->GetXaxis()->FindBin(1.4));
      m2DistrPosProj->GetYaxis()->SetRange(0., m2DistrPosProj->GetMaximumBin());
      m2DistrNegProj->GetXaxis()->SetRange(m2DistrNegProj->GetXaxis()->FindBin(-0.4), 
                                           m2DistrNegProj->GetXaxis()->FindBin(1.4));
      m2DistrNegProj->GetYaxis()->SetRange(0., m2DistrNegProj->GetMaximumBin());

      m2DistrPosProj->SetTitle("");
      m2DistrPosProj->GetXaxis()->SetTitle("m^{2} [GeV/c^{2}]");
      m2DistrPosProj->SetLineColor(kBlack);
      m2DistrPosProj->SetMarkerColor(kBlack);
      m2DistrPosProj->SetMarkerStyle(8);
      m2DistrPosProj->SetMarkerSize(0.5);
      m2DistrPosProj->SetTitleSize(0.08, "X");
      m2DistrPosProj->GetXaxis()->SetLabelSize(0.08);
      m2DistrPosProj->GetYaxis()->SetLabelSize(0.08);
      m2DistrPosProj->GetXaxis()->SetTitleOffset(0.85);

      m2DistrNegProj->SetTitle("");
      m2DistrNegProj->GetXaxis()->SetTitle("m^{2} [GeV/c^{2}]");
      m2DistrNegProj->SetLineColor(kBlack);
      m2DistrNegProj->SetMarkerColor(kBlack);
      m2DistrNegProj->SetMarkerStyle(8);
      m2DistrNegProj->SetMarkerSize(0.5);
      m2DistrNegProj->SetTitleSize(0.08, "X");
      m2DistrNegProj->GetXaxis()->SetLabelSize(0.08);
      m2DistrNegProj->GetYaxis()->SetLabelSize(0.08);
      m2DistrNegProj->GetXaxis()->SetTitleOffset(0.85);
		
      // positive tracks
		m2SingleFitCanv.cd(1);	
		gPad->SetPad(0.02, 0.51, 1., 1.);
		gPad->SetBottomMargin(0.17);
      gPad->SetLogy();

      m2DistrPosProj->Draw("E");
      m2DistrPosProj->Draw("SAME AXIS X+ Y+");

      if (binPTMin > detector["pi+_pt_bounds"][0].as<double>() - 1e-3 && 
          binPTMax < detector["pi+_pt_bounds"][1].as<double>() + 1e-3)
      {
         PerformSingleM2Fit(m2DistrPosProj, detector["sigmalized_identification_range"].as<double>(),
                            pT, fitPiPlus, detector["pion_bg_func"].as<std::string>());
      }
      if (binPTMin > detector["k+_pt_bounds"][0].as<double>() - 1e-3 && 
          binPTMax < detector["k+_pt_bounds"][1].as<double>() + 1e-3)
      {
         PerformSingleM2Fit(m2DistrPosProj, detector["sigmalized_identification_range"].as<double>(),
                            pT, fitKPlus, detector["kaon_bg_func"].as<std::string>());
      }
      if (binPTMin > detector["p_pt_bounds"][0].as<double>() - 1e-3 && 
          binPTMax < detector["p_pt_bounds"][1].as<double>() + 1e-3)
      {
         PerformSingleM2Fit(m2DistrPosProj, detector["sigmalized_identification_range"].as<double>(),
                            pT, fitP, detector["proton_bg_func"].as<std::string>());
      }

      // negative tracks
		m2SingleFitCanv.cd(2);
		gPad->SetPad(0.02, 0., 1., 0.5);
		gPad->SetBottomMargin(0.17);
      gPad->SetLogy();

      m2DistrNegProj->Draw("E");
      m2DistrNegProj->Draw("SAME AXIS X+ Y+");

      if (binPTMin > detector["pi-_pt_bounds"][0].as<double>() - 1e-3 && 
          binPTMax < detector["pi-_pt_bounds"][1].as<double>() + 1e-3)
      {
         PerformSingleM2Fit(m2DistrNegProj, detector["sigmalized_identification_range"].as<double>(),
                            pT, fitPiMinus, detector["pion_bg_func"].as<std::string>());
      }
      if (binPTMin > detector["k-_pt_bounds"][0].as<double>() - 1e-3 && 
          binPTMax < detector["k-_pt_bounds"][1].as<double>() + 1e-3)
      {
         PerformSingleM2Fit(m2DistrNegProj, detector["sigmalized_identification_range"].as<double>(),
                            pT, fitKMinus, detector["kaon_bg_func"].as<std::string>());
      }
      if (binPTMin > detector["pbar_pt_bounds"][0].as<double>() - 1e-3 && 
          binPTMax < detector["pbar_pt_bounds"][1].as<double>() + 1e-3)
      {
         PerformSingleM2Fit(m2DistrNegProj, detector["sigmalized_identification_range"].as<double>(),
                            pT, fitPBar, detector["proton_bg_func"].as<std::string>());
      }

      ROOTTools::PrintCanvas(&m2SingleFitCanv, outputDir + "/" + 
                             detector["name"].as<std::string>() + "/" + 
                             CppTools::DtoStr(pTBin["min"].as<double>(), 1) + "-" +
                             CppTools::DtoStr(pTBin["max"].as<double>(), 1));
      numberOfCalls++;
   }
   
   pBar.Finish();

   fitPiPlus.meansVsPTFit->SetRange(detector["pi+_pt_bounds"][0].as<double>() - 0.01, 
                                    detector["pi+_pt_bounds"][1].as<double>() + 0.01);
   fitKPlus.meansVsPTFit->SetRange(detector["k+_pt_bounds"][0].as<double>() - 0.01, 
                                   detector["k+_pt_bounds"][1].as<double>() + 0.01);
   fitP.meansVsPTFit->SetRange(detector["p_pt_bounds"][0].as<double>() - 0.01, 
                               detector["p_pt_bounds"][1].as<double>() + 0.01);
   fitPiMinus.meansVsPTFit->SetRange(detector["pi-_pt_bounds"][0].as<double>() - 0.01, 
                                     detector["pi-_pt_bounds"][1].as<double>() + 0.01);
   fitKMinus.meansVsPTFit->SetRange(detector["k-_pt_bounds"][0].as<double>() - 0.01, 
                                    detector["k-_pt_bounds"][1].as<double>() + 0.01);
   fitPBar.meansVsPTFit->SetRange(detector["pbar_pt_bounds"][0].as<double>() - 0.01, 
                                  detector["pbar_pt_bounds"][1].as<double>() + 0.01);

   fitPiPlus.sigmasVsPTFit->SetRange(detector["pi+_pt_bounds"][0].as<double>() - 0.01, 
                                     detector["pi+_pt_bounds"][1].as<double>() + 0.01);
   fitKPlus.sigmasVsPTFit->SetRange(detector["k+_pt_bounds"][0].as<double>() - 0.01, 
                                    detector["k+_pt_bounds"][1].as<double>() + 0.01);
   fitP.sigmasVsPTFit->SetRange(detector["p_pt_bounds"][0].as<double>() - 0.01, 
                                detector["p_pt_bounds"][1].as<double>() + 0.01);
   fitPiMinus.sigmasVsPTFit->SetRange(detector["pi-_pt_bounds"][0].as<double>() - 0.01, 
                                      detector["pi-_pt_bounds"][1].as<double>() + 0.01);
   fitKMinus.sigmasVsPTFit->SetRange(detector["k-_pt_bounds"][0].as<double>() - 0.01, 
                                     detector["k-_pt_bounds"][1].as<double>() + 0.01);
   fitPBar.sigmasVsPTFit->SetRange(detector["pbar_pt_bounds"][0].as<double>() - 0.01, 
                                   detector["pbar_pt_bounds"][1].as<double>() + 0.01);

   fitPiPlus.meansVsPT.Fit(fitPiPlus.meansVsPTFit.get(), "RQMBN");
   fitKPlus.meansVsPT.Fit(fitKPlus.meansVsPTFit.get(), "RQMBN");
   fitP.meansVsPT.Fit(fitP.meansVsPTFit.get(), "RQMBN");
   fitPiMinus.meansVsPT.Fit(fitPiMinus.meansVsPTFit.get(), "RQMBN");
   fitKMinus.meansVsPT.Fit(fitKMinus.meansVsPTFit.get(), "RQMBN");
   fitPBar.meansVsPT.Fit(fitPBar.meansVsPTFit.get(), "RQMBN");

   for (int i = 0; i < fitPBar.meansVsPTFit->GetNpar(); i++)
   {
      fitPiPlus.sigmasVsPTFit->FixParameter(i, fitPiPlus.meansVsPTFit->GetParameter(i));
      fitKPlus.sigmasVsPTFit->FixParameter(i, fitKPlus.meansVsPTFit->GetParameter(i));
      fitP.sigmasVsPTFit->FixParameter(i, fitP.meansVsPTFit->GetParameter(i));
      fitPiMinus.sigmasVsPTFit->FixParameter(i, fitPiMinus.meansVsPTFit->GetParameter(i));
      fitKMinus.sigmasVsPTFit->FixParameter(i, fitKMinus.meansVsPTFit->GetParameter(i));
      fitPBar.sigmasVsPTFit->FixParameter(i, fitPBar.meansVsPTFit->GetParameter(i));
   }

   fitPiPlus.sigmasVsPT.Fit(fitPiPlus.sigmasVsPTFit.get(), "RQMBN");
   fitKPlus.sigmasVsPT.Fit(fitKPlus.sigmasVsPTFit.get(), "RQMBN");
   fitPBar.sigmasVsPT.Fit(fitPBar.sigmasVsPTFit.get(), "RQMBN");
   fitPiMinus.sigmasVsPT.Fit(fitPiMinus.sigmasVsPTFit.get(), "RQMBN");
   fitKMinus.sigmasVsPT.Fit(fitKMinus.sigmasVsPTFit.get(), "RQMBN");
   fitP.sigmasVsPT.Fit(fitP.sigmasVsPTFit.get(), "RQMBN");

   if (!detector["is_calibrated"].as<bool>())
   {
      CppTools:: PrintInfo("From apporximations for " + detector["name"].as<std::string>());
      CppTools::Print(" pions: sigma_alpha=" + 
                      std::to_string(CppTools::Average(fitPiPlus.sigmasVsPTFit->GetParameter(2),
                                                       fitPiMinus.sigmasVsPTFit->GetParameter(2)))+ 
                      ", sigma_ms=" + 
                      std::to_string(CppTools::Average(fitPiPlus.sigmasVsPTFit->GetParameter(3),
                                                       fitPiMinus.sigmasVsPTFit->GetParameter(3))) + 
                      ", sigma_t=" + 
                      std::to_string(CppTools::Average(fitPiPlus.sigmasVsPTFit->GetParameter(4),
                                                       fitPiMinus.sigmasVsPTFit->GetParameter(4))));
      CppTools::Print(" kaons: sigma_alpha=" + 
                      std::to_string(CppTools::Average(fitKPlus.sigmasVsPTFit->GetParameter(2),
                                                       fitKMinus.sigmasVsPTFit->GetParameter(2)))+ 
                      ", sigma_ms=" + 
                      std::to_string(CppTools::Average(fitKPlus.sigmasVsPTFit->GetParameter(3),
                                                       fitKMinus.sigmasVsPTFit->GetParameter(3))) + 
                      ", sigma_t=" + 
                      std::to_string(CppTools::Average(fitKPlus.sigmasVsPTFit->GetParameter(4),
                                                       fitKMinus.sigmasVsPTFit->GetParameter(4))));
      CppTools::Print(" p and pbar: sigma_alpha=" + 
                      std::to_string(CppTools::Average(fitP.sigmasVsPTFit->GetParameter(2),
                                                       fitPBar.sigmasVsPTFit->GetParameter(2)))+ 
                      ", sigma_ms=" + 
                      std::to_string(CppTools::Average(fitP.sigmasVsPTFit->GetParameter(3),
                                                       fitPBar.sigmasVsPTFit->GetParameter(3))) + 
                      ", sigma_t=" + 
                      std::to_string(CppTools::Average(fitP.sigmasVsPTFit->GetParameter(4),
                                                       fitPBar.sigmasVsPTFit->GetParameter(4))));
   }

   double sigmaAlpha = CppTools::Average(fitPiPlus.sigmasVsPTFit->GetParameter(2), 
                                         fitPiMinus.sigmasVsPTFit->GetParameter(2), 
                                         fitKPlus.sigmasVsPTFit->GetParameter(2), 
                                         fitKMinus.sigmasVsPTFit->GetParameter(2),
                                         fitP.sigmasVsPTFit->GetParameter(2),
                                         fitPBar.sigmasVsPTFit->GetParameter(2));

   double sigmaMS = CppTools::Average(fitPiPlus.sigmasVsPTFit->GetParameter(3), 
                                      fitPiMinus.sigmasVsPTFit->GetParameter(3), 
                                      fitKPlus.sigmasVsPTFit->GetParameter(3), 
                                      fitKMinus.sigmasVsPTFit->GetParameter(3), 
                                      fitP.sigmasVsPTFit->GetParameter(3),
                                      fitPBar.sigmasVsPTFit->GetParameter(3));

   double sigmaT = CppTools::Average(fitPiPlus.sigmasVsPTFit->GetParameter(4), 
                                     fitPiMinus.sigmasVsPTFit->GetParameter(4), 
                                     fitKPlus.sigmasVsPTFit->GetParameter(4), 
                                     fitKMinus.sigmasVsPTFit->GetParameter(4), 
                                     fitP.sigmasVsPTFit->GetParameter(4),
                                     fitPBar.sigmasVsPTFit->GetParameter(4));

   fitPiPlus.sigmasVsPTFit->FixParameter(2, sigmaAlpha);
   fitKPlus.sigmasVsPTFit->FixParameter(2, sigmaAlpha);
   fitPiMinus.sigmasVsPTFit->FixParameter(2, sigmaAlpha);
   fitKMinus.sigmasVsPTFit->FixParameter(2, sigmaAlpha);
   fitP.sigmasVsPTFit->FixParameter(2, sigmaAlpha);
   fitPBar.sigmasVsPTFit->FixParameter(2, sigmaAlpha);

   fitPiPlus.sigmasVsPTFit->FixParameter(3, sigmaMS);
   fitKPlus.sigmasVsPTFit->FixParameter(3, sigmaMS);
   fitPiMinus.sigmasVsPTFit->FixParameter(3, sigmaMS);
   fitKMinus.sigmasVsPTFit->FixParameter(3, sigmaMS);
   fitP.sigmasVsPTFit->FixParameter(3, sigmaMS);
   fitPBar.sigmasVsPTFit->FixParameter(3, sigmaMS);

   fitPiPlus.sigmasVsPTFit->FixParameter(4, sigmaT);
   fitKPlus.sigmasVsPTFit->FixParameter(4, sigmaT);
   fitPiMinus.sigmasVsPTFit->FixParameter(4, sigmaT);
   fitKMinus.sigmasVsPTFit->FixParameter(4, sigmaT);
   fitP.sigmasVsPTFit->FixParameter(4, sigmaT);
   fitPBar.sigmasVsPTFit->FixParameter(4, sigmaT);

   if (!detector["is_calibrated"].as<bool>())
   {
      CppTools:: PrintInfo("Suggested new parameters for " + detector["name"].as<std::string>() + 
                           ":\n sigma_alpha=" + std::to_string(sigmaAlpha) + 
                           ", sigma_ms=" + std::to_string(sigmaMS) + 
                           ", sigma_t=" + std::to_string(sigmaT));
   }

   TLine trueM2LinePi(pTMin, fitPiPlus.m2, pTMax, fitPiPlus.m2);
   TLine trueM2LineK(pTMin, fitKPlus.m2, pTMax, fitKPlus.m2);
   TLine trueM2LineP(pTMin, fitP.m2, pTMax, fitP.m2);

   trueM2LinePi.SetLineColorAlpha(kGray + 3, 0.5);
   trueM2LinePi.SetLineStyle(6);
   trueM2LinePi.SetLineWidth(3);

   trueM2LineK.SetLineColorAlpha(kGray + 3, 0.5);
   trueM2LineK.SetLineStyle(6);
   trueM2LineK.SetLineWidth(3);

   trueM2LineP.SetLineColorAlpha(kGray + 3, 0.5);
   trueM2LineP.SetLineStyle(6);
   trueM2LineP.SetLineWidth(3);

   TCanvas fitParVsPTCanv("m2 id fit parameters vs pT", "", 1000, 1000);

   fitParVsPTCanv.Divide(2, 2);

   TText text;

   text.SetTextFont(43);
   text.SetTextSize(40);
   text.SetTextAngle(90);

   text.DrawText(0.04, 0.65, "charge = +1");
   text.DrawText(0.04, 0.15, "charge = -1");

   // positive means
   fitParVsPTCanv.cd(1);
   gPad->SetPad(0.04, 0.51, 0.52, 1.);

   gPad->SetLeftMargin(0.155);
   gPad->SetBottomMargin(0.11);

   ROOTTools::DrawFrame(pTMin - 0.05, -0.1, pTMax + 0.05, 1.1,
                        "", "p_{T} [GeV/c]", "#mu_{m^{2}} [GeV/c^{2})^{2}]");

   trueM2LinePi.Draw();
   trueM2LineK.Draw();
   trueM2LineP.Draw();

   fitPiPlus.meansVsPTFit->Draw("SAME");
   fitKPlus.meansVsPTFit->Draw("SAME");
   fitP.meansVsPTFit->Draw("SAME");

   fitPiPlus.meansVsPT.Draw("P");
   fitKPlus.meansVsPT.Draw("P");
   fitP.meansVsPT.Draw("P");

   // positive sigmas
   fitParVsPTCanv.cd(2);
   gPad->SetPad(0.52, 0.51, 1., 1.);

   gPad->SetLeftMargin(0.155);
   gPad->SetBottomMargin(0.11);

   ROOTTools::DrawFrame(pTMin - 0.05, 0., pTMax + 0.05, fitP.sigmasVsPTFit->Eval(pTMax)*1.3, 
                        "", "p_{T} [GeV/c]", "#sigma_{m^{2}} [GeV/c^{2})^{2}]");

   fitPiPlus.sigmasVsPTFit->Draw("SAME");
   fitKPlus.sigmasVsPTFit->Draw("SAME");
   fitP.sigmasVsPTFit->Draw("SAME");

   fitPiPlus.sigmasVsPT.Draw("P");
   fitKPlus.sigmasVsPT.Draw("P");
   fitP.sigmasVsPT.Draw("P");

   // negative means
   fitParVsPTCanv.cd(3);
   gPad->SetPad(0.04, 0.02, 0.52, 0.51);

   gPad->SetLeftMargin(0.155);
   gPad->SetBottomMargin(0.11);

   ROOTTools::DrawFrame(pTMin - 0.05, -0.1, pTMax + 0.05, 1.1, 
                        "", "p_{T} [GeV/c]", "#mu_{m^{2}} [GeV/c^{2})^{2}]");

   trueM2LinePi.Draw();
   trueM2LineK.Draw();
   trueM2LineP.Draw();

   fitPiMinus.meansVsPTFit->Draw("SAME");
   fitKMinus.meansVsPTFit->Draw("SAME");
   fitPBar.meansVsPTFit->Draw("SAME");

   fitPiMinus.meansVsPT.Draw("P");
   fitKMinus.meansVsPT.Draw("P");
   fitPBar.meansVsPT.Draw("P");

   // negative sigmas
   fitParVsPTCanv.cd(4);
   gPad->SetPad(0.52, 0.02, 1., 0.51);

   gPad->SetLeftMargin(0.155);
   gPad->SetBottomMargin(0.11);

   ROOTTools::DrawFrame(pTMin - 0.05, 0., pTMax + 0.05, fitPBar.sigmasVsPTFit->Eval(pTMax)*1.3, 
                        "", "p_{T} [GeV/c]", "#sigma_{m^{2}} [GeV/c^{2})^{2}]");

   fitPiMinus.sigmasVsPTFit->Draw("SAME");
   fitKMinus.sigmasVsPTFit->Draw("SAME");
   fitPBar.sigmasVsPTFit->Draw("SAME");

   fitPiMinus.sigmasVsPT.Draw("P");
   fitKMinus.sigmasVsPT.Draw("P");
   fitPBar.sigmasVsPT.Draw("P");

   ROOTTools::PrintCanvas(&fitParVsPTCanv, outputDir + "/" + detector["name"].as<std::string>() + 
                          "/fitParameters");

   TCanvas m2IdVsPTCanv("m2 identification vs pT", "", 800, 800);

   gPad->SetLeftMargin(0.11);
   gPad->SetBottomMargin(0.11);

   ROOTTools::DrawFrame(-1.*pTMax - 0.05, -0.4, pTMax + 0.05, 1.4, 
                        "", "p_{T} #times charge [GeV/c]", "m^{2} [GeV/c^{2})^{2}]");

   fitPiPlus.extractionRangeLowVsPT.Clone()->Draw("P");
   fitPiPlus.extractionRangeUpVsPT.Clone()->Draw("P");
   fitKPlus.extractionRangeLowVsPT.Clone()->Draw("P");
   fitKPlus.extractionRangeUpVsPT.Clone()->Draw("P");
   fitP.extractionRangeLowVsPT.Clone()->Draw("P");
   fitP.extractionRangeUpVsPT.Clone()->Draw("P");
   fitPiMinus.extractionRangeLowVsPT.Clone()->Draw("P");
   fitPiMinus.extractionRangeUpVsPT.Clone()->Draw("P");
   fitKMinus.extractionRangeLowVsPT.Clone()->Draw("P");
   fitKMinus.extractionRangeUpVsPT.Clone()->Draw("P");
   fitPBar.extractionRangeLowVsPT.Clone()->Draw("P");
   fitPBar.extractionRangeUpVsPT.Clone()->Draw("P");

   // writing parameters in output file
   if (detector["is_calibrated"].as<bool>())
   {
      std::ofstream 
         parametersOutputFile(parametersDir + "/M2Par" + detector["name"].as<std::string>() + ".txt");

      parametersOutputFile << fitPiPlus.meansVsPTFit->GetParameter(0) << " " << 
                             fitPiPlus.meansVsPTFit->GetParameter(1) << std::endl;;
      parametersOutputFile << fitPiMinus.meansVsPTFit->GetParameter(0) << " " << 
                             fitPiMinus.meansVsPTFit->GetParameter(1) << std::endl;
      parametersOutputFile << fitKPlus.meansVsPTFit->GetParameter(0) << " " << 
                             fitKPlus.meansVsPTFit->GetParameter(1) << std::endl;;
      parametersOutputFile << fitKMinus.meansVsPTFit->GetParameter(0) << " " << 
                             fitKMinus.meansVsPTFit->GetParameter(1) << std::endl;
      parametersOutputFile << fitP.meansVsPTFit->GetParameter(0) << " " << 
                             fitP.meansVsPTFit->GetParameter(1) << std::endl;;
      parametersOutputFile << fitPBar.meansVsPTFit->GetParameter(0) << " " << 
                             fitPBar.meansVsPTFit->GetParameter(1) << std::endl;
      parametersOutputFile << sigmaAlpha << " " << sigmaMS << " " << sigmaT << " " <<
                             fitP.sigmasVsPTFit->GetParameter(5) << " " <<
                             fitP.sigmasVsPTFit->GetParameter(6);
   }
   else
   {
      pBar.HandleOutput(CppTools::INFO_PROMPT, "Detector " + detector["name"].as<std::string>() + 
                        " was specified as not being calibrated; change \"is_calibrated\" field\
                        to true in file " + inputYAMLM2Id.GetFileName() + 
                        " for the approximation parameters to be written");
   }

   ROOTTools::PrintCanvas(&m2IdVsPTCanv, outputDir + "/" + 
                          detector["name"].as<std::string>() + "/ms");
}

void M2IdentFit::PerformSingleM2Fit(TH1D *massDistr, const double sigmalizedYieldExtractionRange,
                                    const double pT, FitParameters& fitPar, 
                                    const std::string& funcBG)
{
   TF1 m2FitBG("bg fit", funcBG.c_str());
   TF1 m2FitGaus("fg fit", "gaus");
   TF1 m2Fit("bg+fg fit", (funcBG + "+ gaus(3)").c_str());

   m2Fit.SetParameter(3, massDistr->GetBinContent(massDistr->GetXaxis()->
                                                  FindBin(fitPar.meansVsPTFit->Eval(pT))));
   m2Fit.SetParameter(4, fitPar.meansVsPTFit->Eval(pT));
   m2Fit.SetParameter(5, fitPar.sigmasVsPTFit->Eval(pT));

   m2Fit.SetParLimits(3, 0., massDistr->GetBinContent(massDistr->GetMaximumBin()));
   if (!fitPar.useM2MeansPrevFit)
   {
      m2Fit.SetParLimits(4, fitPar.meansVsPTFit->Eval(pT) - fitPar.sigmasVsPTFit->Eval(pT)/2., 
                         fitPar.meansVsPTFit->Eval(pT) + fitPar.sigmasVsPTFit->Eval(pT)/2.);
   }
   else
   {
      m2Fit.SetParLimits(4, fitPar.meansVsPTFit->Eval(pT)/1.1, fitPar.meansVsPTFit->Eval(pT)*1.1);
   }

   if (!fitPar.isCalibrated)
   {
      m2Fit.SetParLimits(5, fitPar.sigmasVsPTFit->Eval(pT)/1.2, 
                         fitPar.sigmasVsPTFit->Eval(pT)*1.2);
   }
   else
   {
      m2Fit.SetParLimits(5, fitPar.sigmasVsPTFit->Eval(pT)/1.05, 
                         fitPar.sigmasVsPTFit->Eval(pT)*1.05);
   }

   massDistr->Fit(&m2Fit, "RQMBN");

   for (unsigned int i = 1; i <= nFitTries; i++)
   {
      m2Fit.SetParLimits(3, m2Fit.GetParameter(3)/(1. + 1./static_cast<double>(i*i*i)),
                         m2Fit.GetParameter(3)*(1. + 1./static_cast<double>(i*i*i)));
      m2Fit.SetParLimits(4, m2Fit.GetParameter(4) - 
                         m2Fit.GetParameter(5)*(0.25 + 0.25/static_cast<double>(i*i*i)),
                         m2Fit.GetParameter(4) + 
                         m2Fit.GetParameter(5)*(0.25 + 0.25/static_cast<double>(i*i*i)));
      if (!fitPar.isCalibrated)
      {
         m2Fit.SetParLimits(5, m2Fit.GetParameter(5)/(1. + 0.25/static_cast<double>(i*i*i)),
                            m2Fit.GetParameter(5)*(1. + 0.25/static_cast<double>(i*i*i)));
      }
      else
      {
         m2Fit.SetParLimits(5, m2Fit.GetParameter(5)/(1. + 0.025/static_cast<double>(i*i*i)),
                            m2Fit.GetParameter(5)*(1. + 0.025/static_cast<double>(i*i*i)));
      }
      m2Fit.SetRange(m2Fit.GetParameter(4) - m2Fit.GetParameter(5)*
                     (sigmalizedYieldExtractionRange + 0.25/static_cast<double>(i*i*i)),
                     m2Fit.GetParameter(4) + m2Fit.GetParameter(5)*
                     (sigmalizedYieldExtractionRange + 0.25/static_cast<double>(i*i*i)));
      massDistr->Fit(&m2Fit, "RQMBN");
   }

   for (int i = 0; i < m2FitGaus.GetNpar(); i++)
   {
      m2FitGaus.SetParameter(i, m2Fit.GetParameter(i + 3));
   }
   for (int i = 0; i < m2FitBG.GetNpar(); i++)
   {
      m2FitBG.SetParameter(i, m2Fit.GetParameter(i));
   }

   fitPar.meansVsPT.AddPoint(pT, m2Fit.GetParameter(4));
   fitPar.sigmasVsPT.AddPoint(pT, m2Fit.GetParameter(5));

   const double chargeMult = fitPar.isPositive ? 1. : -1;
   fitPar.extractionRangeLowVsPT.AddPoint(chargeMult*pT, m2Fit.GetParameter(4) -
                                          m2Fit.GetParameter(5)*sigmalizedYieldExtractionRange);
   fitPar.extractionRangeUpVsPT.AddPoint(chargeMult*pT, m2Fit.GetParameter(4) +
                                         m2Fit.GetParameter(5)*sigmalizedYieldExtractionRange);

   m2Fit.SetRange(m2Fit.GetParameter(4) - m2Fit.GetParameter(5)*sigmalizedYieldExtractionRange, 
                  m2Fit.GetParameter(4) + m2Fit.GetParameter(5)*sigmalizedYieldExtractionRange);
   m2FitGaus.SetRange(m2Fit.GetParameter(4) - m2Fit.GetParameter(5)*5., 
                      m2Fit.GetParameter(4) + m2Fit.GetParameter(5)*5.);
   m2FitBG.SetRange(m2Fit.GetParameter(4) - m2Fit.GetParameter(5)*sigmalizedYieldExtractionRange,
                    m2Fit.GetParameter(4) + m2Fit.GetParameter(5)*sigmalizedYieldExtractionRange);

   m2Fit.SetLineColor(fitPar.color - 5);
   m2FitGaus.SetLineColor(fitPar.color - 3);
   m2FitBG.SetLineColor(kGray + 1);

   m2FitBG.SetLineStyle(2);
   m2FitGaus.SetLineStyle(2);

   m2Fit.SetLineWidth(3);
   m2FitGaus.SetLineWidth(3);
   m2FitBG.SetLineWidth(3);

   m2FitBG.Clone()->Draw("SAME");
   m2FitGaus.Clone()->Draw("SAME");
   m2Fit.Clone()->Draw("SAME");
}

double M2IdentFit::GetYield(TH1F *hist, const double mean, const double sigma, 
                            const double sigmalizedYieldExtractionRange, TF1 *fitBG,
                            const double vetoLow, const double vetoUp, double& err)
{
   double yield = 0; // yield of a signal
   double yieldNoBGSubtr = 0; // yield of a signal + bg i.e. without bg subtraction

   for (int i = hist->GetXaxis()->
           FindBin(CppTools::Maximum(mean - sigmalizedYieldExtractionRange*sigma, vetoLow)); 
        i <= hist->GetXaxis()->
           FindBin(CppTools::Minimum(mean + sigmalizedYieldExtractionRange*sigma, vetoUp)); i++)
   {
      // this is needed for bin shift correction
      double binAverageFitBG = 0.;
      // taking into account the bin shift correction
      for (double m2 = hist->GetXaxis()->GetBinLowEdge(i); 
           m2 < hist->GetXaxis()->GetBinUpEdge(i) + 1e-15; 
           m2 += hist->GetXaxis()->GetBinWidth(i)/99.)
      {
         binAverageFitBG = fitBG->Eval(m2);
      }

      // 100 iterations in the previous loop so the division is by 100
      binAverageFitBG /= 100.;

      yield += hist->GetBinContent(i) - binAverageFitBG;
      yieldNoBGSubtr += hist->GetBinContent(i);
   }

   // correction to the yield (see M2IdentFit::PerformSingleM2Fit definition)
   const double yieldCorrection = 
      (erf((hist->GetXaxis()->GetBinCenter(hist->GetXaxis()->FindBin(mean)) - 
       hist->GetXaxis()->GetBinLowEdge
       (hist->GetXaxis()->FindBin(CppTools::Maximum(mean - sigmalizedYieldExtractionRange*sigma, 
                                                    vetoLow))))*
       sigmalizedYieldExtractionRange/sqrt(2.)/sigma) +
       erf((hist->GetXaxis()->GetBinUpEdge
       (hist->GetXaxis()->FindBin(CppTools::Minimum(mean + sigmalizedYieldExtractionRange*sigma, 
                                                    vetoUp))) - 
        hist->GetXaxis()->GetBinCenter(hist->GetXaxis()->FindBin(mean)))*
       sigmalizedYieldExtractionRange/sqrt(2.)/sigma))/
      2./erf(sigmalizedYieldExtractionRange/sqrt(2.));

   err = sqrt(yieldNoBGSubtr)/yield/yieldCorrection;

   return yield/yieldCorrection;
}

M2IdentFit::FitParameters::FitParameters(const std::string& particleName, const double massSquared,
                                         const YAML::Node& detector, const bool isPositive, 
                                         const Color_t color)
{
   name = particleName;

   if (name != "pi+" && name != "pi-" && name != "K+" && 
       name != "K-" && name != "p" && name != "pbar") 
   {
      CppTools::PrintError("Unknown charged hadron " + name);
   }

   m2 = massSquared;
   this->isPositive = isPositive;
   this->color = color;

   isCalibrated = detector["is_calibrated"].as<bool>();
   useM2MeansPrevFit = detector["use_m2_mean_par"].as<bool>();

   meansVsPTFit = 
      std::make_unique<TF1>(("means vs pT fit " + name).c_str(),
                            inputYAMLM2Id["means_vs_pt_fit_func"].as<std::string>().c_str());
   sigmasVsPTFit = 
      std::make_unique<TF1>(("sigmas vs pT fit " + name).c_str(),
                            inputYAMLM2Id["sigmas_vs_pt_fit_func"].as<std::string>().c_str());

   meansVsPT.SetMarkerStyle(20);
   sigmasVsPT.SetMarkerStyle(20);
   extractionRangeLowVsPT.SetMarkerStyle(20);
   extractionRangeUpVsPT.SetMarkerStyle(20);

   meansVsPT.SetMarkerColor(color - 3);
   sigmasVsPT.SetMarkerColor(color - 3);
   extractionRangeLowVsPT.SetMarkerColor(color - 3);
   extractionRangeUpVsPT.SetMarkerColor(color - 3);

   meansVsPTFit->SetLineWidth(3);
   sigmasVsPTFit->SetLineWidth(3);

   meansVsPTFit->SetLineColor(color + 3);
   sigmasVsPTFit->SetLineColor(color + 3);

   meansVsPTFit->SetLineStyle(2);
   sigmasVsPTFit->SetLineStyle(2);

   // expected parameters
   const double sigmaAlpha = detector["sigma_alpha"].as<double>();
   const double sigmaMS = detector["sigma_ms"].as<double>();
   const double sigmaT = detector["sigma_t"].as<double>();

   if (!useM2MeansPrevFit)
   {
      meansVsPTFit->SetParameter(0, m2);
      meansVsPTFit->SetParameter(1, 0.);
      sigmasVsPTFit->SetParameter(0, m2);
      sigmasVsPTFit->SetParameter(1, 0.);
   }
   else
   {
      std::ifstream parametersInputFile("data/Parameters/M2Id/" + runName + "/M2Par" + 
                                        detector["name"].as<std::string>() + ".txt");
      double tmp[2];
      parametersInputFile >> tmp[0] >> tmp[1];
      if (name != "pi+") 
      {
         parametersInputFile >> tmp[0] >> tmp[1];
         if (name != "pi-") 
         {
            parametersInputFile >> tmp[0] >> tmp[1];
            if (name != "K+") 
            {
               parametersInputFile >> tmp[0] >> tmp[1];
               if (name != "K-") 
               {
                  parametersInputFile >> tmp[0] >> tmp[1];
                  if (name != "p") parametersInputFile >> tmp[0] >> tmp[1];
               }
            }
         }
      }

      for (int i = 0; i < 2; i++)
      {
         meansVsPTFit->SetParameter(i, tmp[i]);
         sigmasVsPTFit->SetParameter(i, tmp[i]);
         meansVsPTFit->SetParLimits(i, tmp[i]/1.05, tmp[i]*1.05);
      }
   }

   sigmasVsPTFit->FixParameter(5, K1);
   sigmasVsPTFit->FixParameter(6, detector["L"].as<double>());

   if (!isCalibrated)
   {
      sigmasVsPTFit->SetParameter(2, sigmaAlpha);
      sigmasVsPTFit->SetParameter(3, sigmaMS);
      sigmasVsPTFit->SetParameter(4, sigmaT);

      sigmasVsPTFit->SetParLimits(2, sigmaAlpha/1.2, sigmaAlpha*1.2);
      sigmasVsPTFit->SetParLimits(3, sigmaMS/1.2, sigmaMS*1.2);
      sigmasVsPTFit->SetParLimits(4, sigmaT/1.2, sigmaT*1.2);
   }
   else
   {
      sigmasVsPTFit->FixParameter(2, sigmaAlpha);
      sigmasVsPTFit->FixParameter(3, sigmaMS);
      sigmasVsPTFit->FixParameter(4, sigmaT);
   }

   rawYieldsOutputFile.open(rawYieldsDir + "/" + name + ".txt");
}

#endif /* M2_IDENT_FIT_CPP */
