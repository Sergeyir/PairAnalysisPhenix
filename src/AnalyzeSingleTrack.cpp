// $SOURCE$
//------------------------------------------------------------------------------------------------
//                         AnalyzeHeatMaps functions realisation
//------------------------------------------------------------------------------------------------
// AnalyzeHeatMaps
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic macro for evaluating heat maps distributions for different detectors
 * from simulation output of event-like TTrees to processed histograms
 **/
//------------------------------------------------------------------------------------------------

#ifndef ANALYZE_HEAT_MAPS_CPP
#define ANALYZE_HEAT_MAPS_CPP

#include "../include/AnalyzeHeatMaps.hpp"

void Parameters::Init(const std::string inputFileName, const int nThr)
{
   CppTools::CheckInputFile(inputFileName);
   
   numberOfThreads = nThr;
   
   std::ifstream simInputJSON(inputFileName.c_str(), std::ifstream::binary);
   Json::Value simInputJSONContents;
   simInputJSON >> simInputJSONContents;

   runName = simInputJSONContents["run_name"].asString();

   std::ifstream mainInputJSON("input/" + runName + "/main.json", std::ifstream::binary);
   Json::Value mainInputJSONContents;
   mainInputJSON >> mainInputJSONContents;

   collisionSystemName = mainInputJSONContents["collision_system_name"].asString();

   for (auto particle : simInputJSONContents["single_track"]["particles"])
   {
      partQueue.push_back(particle["name"].asString());
   }
   for (auto magf : simInputJSONContents["single_track"]["magnetic_field_configurations"])
   {
      magfQueue.push_back(magf["name"].asString());
   }
   for (auto pTRange : simInputJSONContents["single_track"]["pt_ranges"])
   {
      pTRangeQueue.push_back(pTRange["name"].asString());
   }

   pTMin = simInputJSONContents["single_track"]["pt_min"].asDouble();
   pTMax = simInputJSONContents["single_track"]["pt_max"].asDouble();

   reweightForSpectra = 
      simInputJSONContents["single_track"]["heatmaps_options"]["reweight_for_spectra"].asBool();
   reweightForAlpha = 
      simInputJSONContents["single_track"]["heatmaps_options"]["reweight_for_alpha"].asBool();

   dms.Init(runName);
}

Parameters Par;

TH2F *GetDCHeatmap(TFile *file, const std::string& histName)
{
   TH2F *hist = (TH2F *) file->Get(histName.c_str());
   if (!hist) CppTools::PrintError("Histogram " + histName + " does not exist in file " + 
                                   (std::string) file->GetName());
   // ROOT things to make histograms to not be automaticaly deleted when the file is closed
   hist->SetDirectory(0);
   return hist;
}

void CheckHistsAxis(TH2F *hist1, TH2F *hist2)
{
   // heatmaps from real data and sim are required to have the same axis ranges and number of bins
   if (hist1->GetXaxis()->GetNbins() != hist2->GetXaxis()->GetNbins())
   {
      CppTools::PrintError("Histograms \"" + (std::string) hist1->GetName() + "\" and \"" + 
                           (std::string) hist2->GetName() + "\" have different number of bins on X axis");
   }
   if (hist1->GetYaxis()->GetNbins() != hist2->GetYaxis()->GetNbins())
   {
      CppTools::PrintError("Histograms \"" + (std::string) hist1->GetName() + "\" and \"" + 
                           (std::string) hist2->GetName() + "\" have different number of bins on Y axis");
   }
   if (fabs(hist1->GetXaxis()->GetBinLowEdge(1) - hist2->GetXaxis()->GetBinLowEdge(1)) > 1e-7 ||
       fabs(hist1->GetXaxis()->GetBinLowEdge(hist1->GetXaxis()->GetNbins()) - 
            hist2->GetXaxis()->GetBinLowEdge(hist2->GetXaxis()->GetNbins())) > 1e-7)
   {
      CppTools::PrintError("Histograms \"" + (std::string) hist1->GetName() + "\" and \"" + 
                           (std::string) hist2->GetName() + "\" have different ranges on X axis");
   }
   if (fabs(hist1->GetYaxis()->GetBinLowEdge(1) - hist2->GetYaxis()->GetBinLowEdge(1)) > 1e-7 ||
       fabs(hist1->GetYaxis()->GetBinLowEdge(hist1->GetYaxis()->GetNbins()) - 
            hist2->GetYaxis()->GetBinLowEdge(hist2->GetYaxis()->GetNbins())) > 1e-7)
   {
      CppTools::PrintError("Histograms \"" + (std::string) hist1->GetName() + "\" and \"" + 
                           (std::string) hist2->GetName() + "\" have different ranges on Y axis");
   }
}

// pTRange can be used to specify different statistics of the same dataset e.g. low pT or high pT
void AnalyzeConfiguration(ThrContainer *thrContainer, const std::string& part, 
                          const std::string& magf, const std::string &pTRange)  
{   

   std::string simInputFileName = "data/SimTrees/" + 
      Par.runName + "/SingleTrack/" + part + "_" + pTRange + magf + ".root";
   std::string realDataFileName = "data/Real/" + 
      Par.runName + "/SingleTrack/sum" + magf + ".root";

   TFile simInputFile = TFile(simInputFileName.c_str());
   TFile realDataFile = TFile(realDataFileName.c_str());
    
   TH1F *origPtHist = (TH1F *) simInputFile.Get("orig_pt");
   
   //thredhold is needed since there can be a little noise in the historgram
   const double origPtThreshold = origPtHist->Integral()/
      static_cast<double>(origPtHist->GetXaxis()->GetNbins())/2.;
   
   const double lowPtBound = origPtHist->GetXaxis()->GetBinLowEdge(
      origPtHist->FindFirstBinAbove(origPtThreshold));
   const double upPtBound = origPtHist->GetXaxis()->GetBinUpEdge(
      origPtHist->FindLastBinAbove(origPtThreshold));
   
   double eventNormWeight = 1.;
   if (Par.reweightForSpectra)
   {
      TH1F *centrHist = (TH1F *) realDataFile.Get("centrality");
      if (!centrHist) 
      {
         CppTools::PrintError("Histogram \"centrality\" does not exist in file" + 
                              static_cast<std::string>(realDataFile.GetName()));
      }
      
      eventNormWeight = origPtHist->Integral(
         origPtHist->GetXaxis()->FindBin(Par.pTMin), 
         origPtHist->GetXaxis()->FindBin(Par.pTMax))/
         centrHist->Integral(1, centrHist->GetXaxis()->GetNbins());

      //this normalization is needed to merge 2 files with flat pT distribution with different ranges
      eventNormWeight *= (Par.pTMax - Par.pTMin)/(upPtBound-lowPtBound);
   }

   /*
   Par.pBar.HandleOutput(ErrorHandlerSnippet::INFO + "Processing file " + 
                         simInputFileName + " with original pT distribution: " +
                         DtoStr(lowPtBound) + " < pT < " + DtoStr(upPtBound));
                         */
   
   // weight function for spectra
   std::unique_ptr<TF1> weightFunc;
   if (Par.reweightForSpectra)
   { 
      std::ifstream inputWeightFunc(("data/Spectra/" + Par.collisionSystemName + "/" + 
                                     part + "Fit.json").c_str(), std::ifstream::binary);
      Json::Value inputWeightFuncContents;
      inputWeightFunc >> inputWeightFuncContents;

      weightFunc = std::make_unique<TF1>(
         "weightFunc", inputWeightFuncContents["fit_function"].asString().c_str());
      
      int iPar = 0;
      for (auto par : inputWeightFuncContents["fit_parameters"])
      {
         weightFunc->SetParameter(iPar, par.asDouble());
         iPar++;
      }
   }
   
   ROOT::EnableImplicitMT();
   ROOT::TTreeProcessorMT tp(simInputFileName.c_str());
   
   auto ProcessMP = [&](TTreeReader &reader)
   {   
      std::shared_ptr<TH1F> distrOrigPT = thrContainer->distrOrigPT.Get();
      std::shared_ptr<TH2F> distrOrigPTVsRecPT = thrContainer->distrOrigPTVsRecPT.Get();
      
      std::shared_ptr<TH2F> heatmapDCe0 = thrContainer->heatmapDCe0.Get();
      std::shared_ptr<TH2F> heatmapDCe1 = thrContainer->heatmapDCe1.Get();
      std::shared_ptr<TH2F> heatmapDCw0 = thrContainer->heatmapDCw0.Get();
      std::shared_ptr<TH2F> heatmapDCw1 = thrContainer->heatmapDCw1.Get();
      
      std::shared_ptr<TH2F> heatmapUnscaledDCe0 = thrContainer->heatmapUnscaledDCe0.Get();
      std::shared_ptr<TH2F> heatmapUnscaledDCw0 = thrContainer->heatmapUnscaledDCw0.Get();
      std::shared_ptr<TH2F> heatmapUnscaledDCe1 = thrContainer->heatmapUnscaledDCe1.Get();
      std::shared_ptr<TH2F> heatmapUnscaledDCw1 = thrContainer->heatmapUnscaledDCw1.Get();
      
      std::shared_ptr<TH2F> heatmapPC1e = thrContainer->heatmapPC1e.Get();
      std::shared_ptr<TH2F> heatmapPC1w = thrContainer->heatmapPC1w.Get();
      
      std::shared_ptr<TH2F> heatmapPC2 = thrContainer->heatmapPC2.Get();
      std::shared_ptr<TH2F> heatmapPC3e = thrContainer->heatmapPC3e.Get();
      std::shared_ptr<TH2F> heatmapPC3w = thrContainer->heatmapPC3w.Get();
      
      std::array<std::shared_ptr<TH2F>, 4> heatmapEMCale, heatmapEMCalw;

      std::array<std::shared_ptr<TH2F>, 4> distrECoreVsPTEMCale, distrECoreVsPTEMCalw;
      
      std::shared_ptr<TH1F> distrStripTOFw = thrContainer->distrStripTOFw.Get();
      std::shared_ptr<TH1F> distrSlatTOFe = thrContainer->distrSlatTOFe.Get();
      std::shared_ptr<TH2F> distrELossTOFe = thrContainer->distrELossTOFe.Get();
      
      std::shared_ptr<TH2F> heatmapTOFe = thrContainer->heatmapTOFe.Get();
      
      std::shared_ptr<TH2F> heatmapTOFw0 = thrContainer->heatmapTOFw0.Get();
      std::shared_ptr<TH2F> heatmapTOFw1 = thrContainer->heatmapTOFw1.Get();
      
      for (int i = 0; i < 4; i++)
      {
         heatmapEMCale[i] = thrContainer->heatmapEMCale[i].Get();
         heatmapEMCalw[i] = thrContainer->heatmapEMCalw[i].Get();
         distrECoreVsPTEMCale[i] = thrContainer->distrECoreVsPTEMCale[i].Get();
         distrECoreVsPTEMCalw[i] = thrContainer->distrECoreVsPTEMCalw[i].Get();
      }
   
      EffTreeReader T(reader);
   
      while (reader.Next())
      {   
         Par.numberOfCalls++;
         const double origPt = sqrt(pow(T.mom_orig(0), 2) + pow(T.mom_orig(1), 2));

         double eventWeight;
         
         if (Par.reweightForSpectra) 
         {
            eventWeight = weightFunc->Eval(origPt)/eventNormWeight;
         }
         else eventWeight = 1.;
         
         distrOrigPT->Fill(origPt, eventWeight);
         
         const double bbcz = T.bbcz();
         if (fabs(bbcz) > 30) continue;

         for(int i = 0; i < T.nch(); i++)
         {
            const double the0 = T.the0(i);
            const double pT = (T.mom(i))*sin(the0);

            if (pT < Par.pTMin) continue;
            if (IsQualityCut(T.qual(i))) continue;
            
            int charge = T.charge(i);
            if (charge != -1 && charge != 1) continue;

            const double zed = T.zed(i);
            if (fabs(zed) > 75 && fabs(zed) < 3) continue;
            
            if (!(fabs(the0)<100 &&
               ((bbcz > 0 && ((bbcz - 250*tan(the0 - 3.1416/2)) > 2 ||
               (bbcz - 200*tan(the0 - 3.1416/2)) < -2)) ||
               (bbcz < 0 && ((bbcz - 250*tan(the0 - 3.1416/2))< -2 ||
               (bbcz - 200*tan(the0 - 3.1416/2)) > 2))))) continue;
   
            //end of basic cuts

            const double alpha = T.alpha(i);
            const double phi = T.phi(i);
            double board;
            
            if (phi > M_PI/2.) board = ((3.72402 - phi + 0.008047*cos(phi + 0.87851))/0.01963496);
            else board = ((0.573231 + phi - 0.0046 * cos(phi + 0.05721))/0.01963496);

            double particleWeight = eventWeight;

            if (phi > M_PI/2.)
            {
               if (zed >=0) 
               {
                  heatmapUnscaledDCe0->Fill(board, alpha, eventWeight);
                  if (Par.reweightForAlpha) particleWeight *= 
                     Par.alphaReweightDCe0->GetBinContent(
                     Par.alphaReweightDCe0->FindBin(alpha));
                  heatmapDCe0->Fill(board, alpha, particleWeight);
               }
               else 
               {
                  heatmapUnscaledDCe1->Fill(board, alpha, eventWeight);
                  if (Par.reweightForAlpha) particleWeight *= 
                     Par.alphaReweightDCe1->GetBinContent(
                     Par.alphaReweightDCe1->FindBin(alpha));
                  heatmapDCe1->Fill(board, alpha, particleWeight);
               }
            }
            else
            {
               if (zed >=0) 
               {
                  heatmapUnscaledDCw0->Fill(board, alpha, eventWeight);
                  if (Par.reweightForAlpha) particleWeight *= 
                     Par.alphaReweightDCw0->GetBinContent(
                     Par.alphaReweightDCw0->FindBin(alpha));
                  heatmapDCw0->Fill(board, alpha, particleWeight);
               }
               else 
               {
                  heatmapUnscaledDCw1->Fill(board, alpha, eventWeight);
                  if (Par.reweightForAlpha) particleWeight *= 
                     Par.alphaReweightDCw1->GetBinContent(
                     Par.alphaReweightDCw1->FindBin(alpha));
                  heatmapDCw1->Fill(board, alpha, particleWeight);
               }
            }
            
            distrOrigPTVsRecPT->Fill(origPt, pT, eventWeight);
            const double pc1phi = atan2(T.ppc1y(i), T.ppc1x(i));
            
            if (phi < M_PI/2.) 
            {
               heatmapPC1w->Fill(T.ppc1z(i), pc1phi, particleWeight);
            }
            else
            {
               if (pc1phi < 0) heatmapPC1e->Fill(T.ppc1z(i), pc1phi + 2.*M_PI, particleWeight);
               else heatmapPC1e->Fill(T.ppc1z(i), pc1phi, particleWeight);
            }
            
            if (IsMatch(Par.dms.GetPC2SDPhi(T.pc2dphi(i), pT, charge), 
                        Par.dms.GetPC2SDZ(T.pc2dz(i), pT, charge), 2., 2.))
            {
               const double pc2z = T.ppc2z(i) - T.pc2dz(i);
               const double pc2phi = atan2(T.ppc2y(i), T.ppc2x(i)) - T.pc2dphi(i);
               
               heatmapPC2->Fill(pc2z, pc2phi, particleWeight);
            }

            if (IsMatch(Par.dms.GetPC3SDPhi(phi, T.pc3dphi(i), pT, charge), 
                        Par.dms.GetPC3SDZ(phi, T.pc3dz(i), pT, charge), 2., 2.))
            {
               const double pc3z = T.ppc3z(i) - T.pc3dz(i);
               double pc3phi = atan2(T.ppc3y(i), T.ppc3x(i) - T.pc3dphi(i));
               
               if (phi > M_PI/2.) 
               {
                  if (pc3phi < 0) heatmapPC3e->Fill(pc3z, pc3phi + 2.*M_PI, particleWeight);
                  else heatmapPC3e->Fill(pc3z, pc3phi, particleWeight);
               }
               else heatmapPC3w->Fill(pc3z, pc3phi, particleWeight);
            }
            
            if (IsHit(T.tofdz(i)))
            {
               const double beta = T.pltof(i)/T.ttof(i)/29.97;
               const double eloss = 0.0014*pow(beta, -1.66);

               distrSlatTOFe->Fill(T.slat(i), particleWeight);
               distrELossTOFe->Fill(beta, T.etof(i));
               
               if (IsMatch(Par.dms.GetTOFeSDPhi(T.tofdphi(i), pT, charge), 
                           Par.dms.GetTOFeSDZ(T.tofdz(i), pT, charge), 2., 2.) && 
                  !Par.dms.IsBadSlat(T.slat(i)) &&
                  T.etof(i) > eloss) 
               {
                  heatmapTOFe->Fill(T.ptofy(i), T.ptofz(i), T.etof(i)*particleWeight);
               }
            }
            else if (IsHit(T.tofwdz(i)))
            {
               if (IsMatch(Par.dms.GetTOFwSDPhi(T.tofwdphi(i), pT, charge), 
                           Par.dms.GetTOFwSDZ(T.tofwdz(i), pT, charge), 2., 2.))
               {
                  distrStripTOFw->Fill(T.striptofw(i), particleWeight);
                  if (!Par.dms.IsBadStripTOFw(static_cast<int>(T.striptofw(i))))
                  {
                     if (T.ptofwy(i) < 100.) 
                     {
                        heatmapTOFw0->Fill(T.ptofwy(i), T.ptofwz(i), particleWeight);
                     }
                     else 
                     {
                        heatmapTOFw1->Fill(T.ptofwy(i), T.ptofwz(i), particleWeight);
                     }
                  }
               }
            }
            
            if (IsHit(T.emcdz(i)))
            {
               if (IsMatch(Par.dms.GetEMCSDPhi(phi, T.emcdphi(i), pT, charge, T.sect(i)), 
                           Par.dms.GetEMCSDZ(phi, T.emcdz(i), pT, charge, T.sect(i)), 2., 2.))
               {   

                  if (phi > M_PI/2.) distrECoreVsPTEMCale[T.sect(i)]->Fill(pT, T.ecore(i));
                  else distrECoreVsPTEMCalw[T.sect(i)]->Fill(pT, T.ecore(i));

                  if (phi > M_PI/2.) 
                  {
                     heatmapEMCale[T.sect(i)]->Fill(T.pemcy(i), T.pemcz(i), 
                                                  T.ecore(i)*particleWeight);
                  }
                  else
                  {   
                     heatmapEMCalw[T.sect(i)]->Fill(T.pemcy(i), T.pemcz(i), 
                                                  T.ecore(i)*particleWeight);
                  }
               }
            }

         }
      }
   };

   tp.Process(ProcessMP);
}

int main(int argc, char **argv)
{

   if (argc < 2 || argc > 3) 
   {
      const std::string errMsg = 
         "Expected 1-2 parameters while " + std::to_string(argc) + 
         " parameter(s) were provided \n Usage: bin/AnalyzeHeatMaps inputJSONName numberOfThreads=std::thread::hardware_concurrency()";
      CppTools::PrintError(errMsg);
   }
   
   if (argc == 2) Par.Init(argv[1], std::thread::hardware_concurrency());
   else Par.Init(argv[1], std::stoi(argv[2]));

   if (Par.reweightForSpectra)
   {
      for (std::string part : Par.partQueue)
      {
         CppTools::CheckInputFile("data/Spectra/" + Par.collisionSystemName + "/" + part + "Fit.json");
      }
   }
   
   for (std::string magf : Par.magfQueue)
   {
      CppTools::CheckInputFile("data/Real/" + Par.runName + "/SingleTrack/sum" + magf + ".root");
      for (std::string part : Par.partQueue)
      {
         for (std::string pTRange : Par.pTRangeQueue)
         {
            const std::string simInputFileName = "data/SimTrees/" + Par.runName + 
               "/SingleTrack/" + part + "_" + pTRange + magf + ".root";
            CppTools::CheckInputFile(simInputFileName);
            const unsigned long numberOfEvents = 
               static_cast<unsigned long>(((TTree *) TFile::Open(simInputFileName.c_str())->
                                                        Get("Tree"))->GetEntries());
            if (numberOfEvents <= 0)
            {
               CppTools::PrintError("Number of events is equal or less than 0 in file " + simInputFileName);
            }
            Par.numberOfEvents += numberOfEvents;
         }
      }
   }
    
   if (Par.reweightForAlpha)
   {
      std::string realDataInputFileName = "data/Real/" + Par.runName + "/SingleTrack/sum.root";
      std::string alphaReweightInputFileName = "data/PostSim/" + Par.runName + "/Heatmaps/all.root";
      
      CppTools::CheckInputFile(realDataInputFileName);
      
      if (CppTools::CheckInputFile(alphaReweightInputFileName, false)) 
      {
         TFile realDataInputFile(realDataInputFileName.c_str());
         TFile alphaReweightInputFile(alphaReweightInputFileName.c_str());

         TH2F *realDataDCe0 = GetDCHeatmap(&realDataInputFile, "Heatmap: DCe, zDC>=0");
         TH2F *realDataDCe1 = GetDCHeatmap(&realDataInputFile, "Heatmap: DCe, zDC<0");
         TH2F *realDataDCw0 = GetDCHeatmap(&realDataInputFile, "Heatmap: DCw, zDC>=0");
         TH2F *realDataDCw1 = GetDCHeatmap(&realDataInputFile, "Heatmap: DCw, zDC<0");

         TH2F *simDCe0 = GetDCHeatmap(&alphaReweightInputFile, "Unscaled heatmap: DCe, zDC>=0");
         TH2F *simDCe1 = GetDCHeatmap(&alphaReweightInputFile, "Unscaled heatmap: DCe, zDC<0");
         TH2F *simDCw0 = GetDCHeatmap(&alphaReweightInputFile, "Unscaled heatmap: DCw, zDC>=0");
         TH2F *simDCw1 = GetDCHeatmap(&alphaReweightInputFile, "Unscaled heatmap: DCw, zDC<0");
         
         CheckHistsAxis(realDataDCe0, simDCe0);
         CheckHistsAxis(realDataDCe1, simDCe1);
         CheckHistsAxis(realDataDCw0, simDCw0);
         CheckHistsAxis(realDataDCw1, simDCw1);
         
         for (int i = 1; i <= realDataDCe0->GetXaxis()->GetNbins(); i++)
         {
            for (int j = 1; j < realDataDCe0->GetYaxis()->GetNbins(); j++)
            {
               const double xVal = realDataDCe0->GetXaxis()->GetBinCenter(i);
               const double yVal = realDataDCe0->GetYaxis()->GetBinCenter(i);

               if (Par.dms.IsDeadDC(2., 1., xVal, yVal))
               {
                  realDataDCe0->SetBinContent(i, j, 0.);
                  simDCe0->SetBinContent(i, j, 0.);
               }
               if (Par.dms.IsDeadDC(2., -1., xVal, yVal))
               {
                  realDataDCe1->SetBinContent(i, j, 0.);
                  simDCe1->SetBinContent(i, j, 0.);
               }
               if (Par.dms.IsDeadDC(1., 1., xVal, yVal))
               {
                  realDataDCw0->SetBinContent(i, j, 0.);
                  simDCw0->SetBinContent(i, j, 0.);
               }
               if (Par.dms.IsDeadDC(1., -1., xVal, yVal))
               {
                  realDataDCw1->SetBinContent(i, j, 0.);
                  simDCw1->SetBinContent(i, j, 0.);
               }
            }
         }

         Par.alphaReweightDCe0 = (TH1F *)
            realDataDCe0->ProjectionY("dce0_reweight",
            1, realDataDCe0->GetXaxis()->GetNbins())->Clone();
         Par.alphaReweightDCe1 = (TH1F *) 
            realDataDCe1->ProjectionY("dce1_reweight",
            1, realDataDCe1->GetXaxis()->GetNbins())->Clone();
         Par.alphaReweightDCw0 = (TH1F *) 
            realDataDCw0->ProjectionY("dcw0_reweight",
            1, realDataDCw0->GetXaxis()->GetNbins())->Clone();
         Par.alphaReweightDCw1 = (TH1F *) 
            realDataDCw1->ProjectionY("dcw1_reweight",
            1, realDataDCw1->GetXaxis()->GetNbins())->Clone();

         Par.alphaReweightDCe0->Scale(simDCe0->Integral()/Par.alphaReweightDCe0->Integral());
         Par.alphaReweightDCe1->Scale(simDCe1->Integral()/Par.alphaReweightDCe1->Integral());
         Par.alphaReweightDCw0->Scale(simDCw0->Integral()/Par.alphaReweightDCw0->Integral());
         Par.alphaReweightDCw1->Scale(simDCw1->Integral()/Par.alphaReweightDCw1->Integral());

         Par.alphaReweightDCe0->Divide(simDCe0->ProjectionY("dce0_proj",
            1, simDCe0->GetXaxis()->GetNbins()));
         Par.alphaReweightDCe1->Divide(simDCe1->ProjectionY("dce1_proj",
            1, simDCe1->GetXaxis()->GetNbins()));
         Par.alphaReweightDCw0->Divide(simDCw0->ProjectionY("dcw0_proj",
            1, simDCw0->GetXaxis()->GetNbins()));
         Par.alphaReweightDCw1->Divide(simDCw1->ProjectionY("dcw1_proj",
            1, simDCw1->GetXaxis()->GetNbins()));

         //capping alpha reweight since experiments extends to higher pT than simulation
         //capped values do not affect the simulation since they are statisticaly insufficient
         //capping is only needed to exclude very big weights in some points in the DC map
         //that are not used for better visibility
         for (int i = 0; i < Par.alphaReweightDCe0->GetXaxis()->GetNbins(); i++)
         {
            if (Par.alphaReweightDCe0->GetBinContent(i) > 100.)
            {
               Par.alphaReweightDCe0->SetBinContent(i, 100.);
            }
         }
         for (int i = 0; i < Par.alphaReweightDCe1->GetXaxis()->GetNbins(); i++)
         {
            if (Par.alphaReweightDCe1->GetBinContent(i) > 100.)
            {
               Par.alphaReweightDCe1->SetBinContent(i, 100.);
            }
         }
         for (int i = 0; i < Par.alphaReweightDCw0->GetXaxis()->GetNbins(); i++)
         {
            if (Par.alphaReweightDCw0->GetBinContent(i) > 100.)
            {
               Par.alphaReweightDCw0->SetBinContent(i, 100.);
            }
         }
         for (int i = 0; i < Par.alphaReweightDCw1->GetXaxis()->GetNbins(); i++)
         {
            if (Par.alphaReweightDCw1->GetBinContent(i) > 100.)
            {
               Par.alphaReweightDCw1->SetBinContent(i, 100.);
            }
         }

         std::string alphaReweightOutputFileName = "data/PostSim/" + 
            Par.runName + "/alpha_reweight.root";
         TFile alphaReweightOutputFile = TFile(alphaReweightOutputFileName.c_str(), "RECREATE");

         alphaReweightOutputFile.cd();
         
         Par.alphaReweightDCe0->Write();
         Par.alphaReweightDCe1->Write();
         Par.alphaReweightDCw0->Write();
         Par.alphaReweightDCw1->Write();

         Par.alphaReweightDCe0->SetDirectory(0);
         Par.alphaReweightDCe1->SetDirectory(0);
         Par.alphaReweightDCw0->SetDirectory(0);
         Par.alphaReweightDCw1->SetDirectory(0);
         
         alphaReweightOutputFile.Close();
         CppTools::PrintInfo("File " + alphaReweightOutputFileName + " was written");
      }
      else 
      {
         CppTools::PrintInfo("alpha reweight is now disabled");
         Par.reweightForAlpha = false;
      }
   }

   CppTools::PrintInfo("Clearing output directory: data/PostSim/" + Par.runName + "/Heatmaps/");
   system(("mkdir -p data/PostSim/" + Par.runName + "/Heatmaps").c_str());
   system(("rm -r data/PostSim/" + Par.runName + "/Heatmaps/*").c_str());

   CppTools::Box box{"Parameters"};
   
   box.AddEntry("Run name", Par.runName);
   box.AddEntry("Particles list", Par.partQueue);
   if (Par.magfQueue.size() == 1 && Par.magfQueue.front() == "")
   {
      box.AddEntry("Magnetic field list", "run default");
   }
   else box.AddEntry("Magnetic field list", Par.magfQueue);
   box.AddEntry("pT ranges list", Par.pTRangeQueue);
   box.AddEntry("Minimum p_T, GeV", Par.pTMin);
   box.AddEntry("Maximum p_T, GeV", Par.pTMax);
   box.AddEntry("Reweight for pT spectra", Par.reweightForSpectra);
   box.AddEntry("Reweight for alpha", Par.reweightForAlpha);
   box.AddEntry("Number of threads", Par.numberOfThreads);
   box.AddEntry("Number of events to be analyzed, 1e6", 
                static_cast<double>(Par.numberOfEvents)/1e6, 3);
   box.Print();

   bool isProcessFinished = false;

   auto pBarCall = [&]()
   {
      while (!isProcessFinished)
      {
         Par.pBar.Print(static_cast<double>(Par.numberOfCalls)/
                        static_cast<double>(Par.numberOfEvents));
         std::this_thread::sleep_for(std::chrono::milliseconds(100));
      }
      Par.pBar.Print(1.);
   };
   
   std::thread pBarThread(pBarCall);
   
   for (std::string part : Par.partQueue)
   {
      for (std::string magf : Par.magfQueue)
      {
         ThrContainer thrContainer;
         for (std::string pTRange : Par.pTRangeQueue)
         {
            AnalyzeConfiguration(&thrContainer, part, magf, pTRange);
         }
         // wriiting the result
         std::string outputFileName = "data/PostSim/" + Par.runName + "/Heatmaps/" + part;
         if (magf != "") outputFileName += "magf";
         outputFileName += magf + ".root";
         
         TFile outfile(outputFileName.c_str(), "RECREATE");
         outfile.cd();
         ThrObjHolder.Write();
         outfile.Close();
      }
   }

   isProcessFinished = true;
   pBarThread.join();

   CppTools::PrintInfo("Merging output files into one");

   system(("hadd -f data/PostSim/" + 
           Par.runName + "/Heatmaps/all.root data/PostSim/" + 
           Par.runName + "/Heatmaps/*.root").c_str());

   return 0;
}

#endif /* ANALYZE_HEAT_MAPS_CPP */
