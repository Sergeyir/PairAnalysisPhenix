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

//auxName can be used to specify different statistics of the same dataset e.g. low pT or high pT
void Analyze(ThrContainer *thrContainer, const std::string& part, const std::string& magf, 
             const std::string &auxName, const int procNum)  
{   
   Box box = Box("Parameters of run " + std::to_string(procNum) + 
      " out of " + std::to_string(Par.partQueue.size()*
      Par.magfQueue.size()*Par.auxNameQueue.size()));

   std::string simInputFileName = Par.simDataDir + 
      Par.runName + "/SingleTrack/" + part + magf + auxName + ".root";
   std::string realDataFileName = Par.realDataDir + 
      Par.runName + "/sum" + magf + ".root";

   TFile simInputFile = TFile(simInputFileName.c_str());
   TFile realDataFile = TFile(realDataFileName.c_str());
   
   const double nevents = static_cast<double>(((TTree *) simInputFile.Get("Tree"))->GetEntries());

   if (nevents <= 0)
   {
      Print("Error: Number of events is equal or less than 0!");
      exit(1);
   }

   
   TH1F *origPtHist = (TH1F *) simInputFile.Get("orig_pt");
   
   //thredhold is needed since there can be a little noise in the historgram
   const double origPtThreshold = origPtHist->Integral()/
      static_cast<double>(origPtHist->GetXaxis()->GetNbins())/2.;
   
   const double lowPtBound = origPtHist->GetXaxis()->GetBinLowEdge(
      origPtHist->FindFirstBinAbove(origPtThreshold));
   const double upPtBound = origPtHist->GetXaxis()->GetBinUpEdge(
      origPtHist->FindLastBinAbove(origPtThreshold));
   
   double eventNormWeight = 1.;
   if (Par.doUseWeightFunc)
   {
      TH1F *centrHist = (TH1F *) realDataFile.Get("centrality");
      
      eventNormWeight = origPtHist->Integral(
         origPtHist->GetXaxis()->FindBin(Par.pTMin), 
         origPtHist->GetXaxis()->FindBin(Par.pTMax))/
         centrHist->Integral(1, centrHist->GetXaxis()->GetNbins());

      //this normalization is needed to merge 2 files with flat pT distribution with different ranges
      eventNormWeight *= (Par.pTMax - Par.pTMin)/(upPtBound-lowPtBound);
   }
   Print("bruh");
   
   box.AddEntry("Run name", Par.runName);
   box.AddEntry("Orig particle", part);
   box.AddEntry("Magnetic field", magf);
   box.AddEntry("Orig pT distribution span", 
      DtoStr(lowPtBound) + " < pT < " + DtoStr(upPtBound));
   box.AddEntry("Use weight function", Par.doUseWeightFunc);
   box.AddEntry("Reweight alpha", Par.doReweightAlpha);
   
   box.AddEntry("Minimum p_T, GeV", Par.pTMin);
   box.AddEntry("Maximum p_T, GeV", Par.pTMax);
   box.AddEntry("Number of events to be analyzed, 1e6", nevents/1e6, 3);

   box.AddEntry("Number of threads", Par.nThreads);
   
   box.Print();
   
   /*
   //tsallis weight function
   std::unique_ptr<TF1> weightFunc = std::make_unique<TF1>(TF1("weightFunc", 
      "[0]*x^([5]-1)*(1 + ([1] - 1)*(sqrt([4] + x^2)^([5]) - [2])/[3])^(-1./([1]-1.))"));
   weightFunc->SetParameters(
      ReadFileIntoArray("../input/SimAnalysis/Spectra/" + Par.system + "/" + part + ".txt", 6));
      */
   
   ROOT::EnableImplicitMT(Par.nThreads);
   ROOT::TTreeProcessorMT tp(simInputFileName.c_str());
   
   double ncalls = 0.;

   bool isProcessFinished = false;

   auto ProcessMP = [&](TTreeReader &reader)
   {   
      std::shared_ptr<TH1F> distrOrigPT = thrContainer->distrOrigPT.Get();
      std::shared_ptr<TH2F> distrOrigPTVsRecPT = thrContainer->distrOrigPTVsRecPT.Get();
      
      std::shared_ptr<TH2F> distrDCe0 = thrContainer->distrDCe0.Get();
      std::shared_ptr<TH2F> distrDCe1 = thrContainer->distrDCe1.Get();
      std::shared_ptr<TH2F> distrDCw0 = thrContainer->distrDCw0.Get();
      std::shared_ptr<TH2F> distrDCw1 = thrContainer->distrDCw1.Get();

      std::shared_ptr<TH2F> distrUnscaledDCe0 = thrContainer->distrUnscaledDCe0.Get();
      std::shared_ptr<TH2F> distrUnscaledDCw0 = thrContainer->distrUnscaledDCw0.Get();
      std::shared_ptr<TH2F> distrUnscaledDCe1 = thrContainer->distrUnscaledDCe1.Get();
      std::shared_ptr<TH2F> distrUnscaledDCw1 = thrContainer->distrUnscaledDCw1.Get();

      std::shared_ptr<TH2F> distrPC1e = thrContainer->distrPC1e.Get();
      std::shared_ptr<TH2F> distrPC1w = thrContainer->distrPC1w.Get();
      
      std::shared_ptr<TH2F> distrPC2 = thrContainer->distrPC2.Get();
      std::shared_ptr<TH2F> distrPC3e = thrContainer->distrPC3e.Get();
      std::shared_ptr<TH2F> distrPC3w = thrContainer->distrPC3w.Get();
      
      std::array<std::shared_ptr<TH2F>, 4> distrEMCalePos, distrEMCaleNeg, distrEMCalwPos, distrEMCalwNeg;

      std::shared_ptr<TH1F> distrStripTOFw = thrContainer->distrStripTOFw.Get();
      std::shared_ptr<TH1F> distrSlatTOFe = thrContainer->distrSlatTOFe.Get();
      
      std::shared_ptr<TH2F> distrTOFe0 = thrContainer->distrTOFe0.Get();
      std::shared_ptr<TH2F> distrTOFe1 = thrContainer->distrTOFe1.Get();
      
      std::shared_ptr<TH2F> distrTOFw0 = thrContainer->distrTOFw0.Get();
      std::shared_ptr<TH2F> distrTOFw1 = thrContainer->distrTOFw1.Get();

      for (int i = 0; i < 4; i++)
      {
         distrEMCalePos[i] = thrContainer->distrEMCalePos[i].Get();
         distrEMCaleNeg[i] = thrContainer->distrEMCaleNeg[i].Get();
         distrEMCalwPos[i] = thrContainer->distrEMCalwPos[i].Get();
         distrEMCalwNeg[i] = thrContainer->distrEMCalwNeg[i].Get();
      }
   
      EffTreeReader T(reader);
   
      while (reader.Next())
      {   
         ncalls += 1.;
         const double origPt = sqrt(pow(T.mom_orig(0), 2) + pow(T.mom_orig(1), 2));

         double eventWeight;
         
         if (Par.doUseWeightFunc) 
         {
            eventWeight = 1;//weightFunc->Eval(origPt)/eventNormWeight;
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
                  distrUnscaledDCe0->Fill(board, alpha, eventWeight);
                  if (Par.doReweightAlpha) particleWeight *= 
                     Par.alphaReweightDCe0->GetBinContent(
                     Par.alphaReweightDCe0->FindBin(alpha));
                  distrDCe0->Fill(board, alpha, particleWeight);
               }
               else 
               {
                  distrUnscaledDCe1->Fill(board, alpha, eventWeight);
                  if (Par.doReweightAlpha) particleWeight *= 
                     Par.alphaReweightDCe1->GetBinContent(
                     Par.alphaReweightDCe1->FindBin(alpha));
                  distrDCe1->Fill(board, alpha, particleWeight);
               }
            }
            else
            {
               if (zed >=0) 
               {
                  distrUnscaledDCw0->Fill(board, alpha, eventWeight);
                  if (Par.doReweightAlpha) particleWeight *= 
                     Par.alphaReweightDCw0->GetBinContent(
                     Par.alphaReweightDCw0->FindBin(alpha));
                  distrDCw0->Fill(board, alpha, particleWeight);
               }
               else 
               {
                  distrUnscaledDCw1->Fill(board, alpha, eventWeight);
                  if (Par.doReweightAlpha) particleWeight *= 
                     Par.alphaReweightDCw1->GetBinContent(
                     Par.alphaReweightDCw1->FindBin(alpha));
                  distrDCw1->Fill(board, alpha, particleWeight);
               }
            }
            
            //if (IsDeadDC(phi, zed, board, alpha)) continue;

            const double pc1phi = atan2(T.ppc1y(i), T.ppc1x(i));
            
            if (phi < M_PI/2.) 
            {
               distrPC1w->Fill(T.ppc1z(i), pc1phi, particleWeight);
               //if (IsDeadPC1(phi, T.ppc1z(i), pc1phi)) continue;
            }
            else
            {
               if (pc1phi < 0)
               {
                  distrPC1e->Fill(T.ppc1z(i), pc1phi + 2.*M_PI, particleWeight);
                  //if (IsDeadPC1(phi, T.ppc1z(i), pc1phi + 2.*M_PI)) continue;
               }
               else
               {
                  distrPC1e->Fill(T.ppc1z(i), pc1phi, particleWeight);
                  //if (IsDeadPC1(phi, T.ppc1z(i), pc1phi)) continue;
               }
            }

            continue;
            

            distrOrigPTVsRecPT->Fill(origPt, pT, eventWeight);
            
            if (IsMatch(T.pc2sdz(i), T.pc2sdphi(i), 2., 2.))
            {
               const double pc2z = T.ppc2z(i) - T.pc2dz(i);
               const double pc2phi = atan2(T.ppc2y(i), T.ppc2x(i)) - T.pc2dphi(i);
               
               distrPC2->Fill(pc2z, pc2phi, particleWeight);
            }

            if (IsMatch(T.pc3sdz(i), T.pc3sdphi(i), 2., 2.))
            {
               const double pc3z = T.ppc3z(i) - T.pc3dz(i);
               double pc3phi = atan2(T.ppc3y(i), T.ppc3x(i) - T.pc3dphi(i));
               
               if (phi > M_PI/2.) 
               {
                  if (pc3phi < 0) pc3phi += M_PI*2.;
                  distrPC3e->Fill(pc3z, pc3phi, particleWeight);
               }
               else distrPC3w->Fill(pc3z, pc3phi, particleWeight);
            }
            
            if (IsHit(T.tofdz(i)))
            {
               const double beta = T.pltof(i)/T.ttof(i)/29.97;
               const double eloss = 0.0014*pow(beta, -1.66);

               distrSlatTOFe->Fill(T.slat(i), particleWeight);
               
               if (IsMatch(T.tofsdz(i), T.tofsdphi(i), 2., 2.) && 
                  !IsBadSlat(T.slat(i)) &&
                  T.etof(i) > eloss) 
               {
                  if (zed >= 0) 
                  {
                     distrTOFe0->Fill(T.ptofy(i), T.ptofz(i), T.etof(i)*particleWeight);
                  }
                  else 
                  {
                     distrTOFe1->Fill(T.ptofy(i), T.ptofz(i), T.etof(i)*particleWeight);
                  }
               }
            }
            else if (IsHit(T.tofwdz(i)))
            {
               double tofwsdphi = GetTOFwsdphi(0, T.mom(i), T.tofwdphi(i), charge, T.striptofw(i));
               double tofwsdz = GetTOFwsdz(0, T.mom(i), T.tofwdz(i), charge, T.striptofw(i));

               tofwsdphi = RecalTOFwsdphi(0, T.mom(i), tofwsdphi, charge, T.striptofw(i));
               tofwsdz = RecalTOFwsdz(0, T.mom(i), tofwsdz, charge, T.striptofw(i));
   
               if (IsMatch(tofwsdphi, tofwsdz, 3., 3.))
               {
                  distrStripTOFw->Fill(T.striptofw(i), particleWeight);
                  if (!IsBadStripTOFw(static_cast<int>(T.striptofw(i))))
                  {
                     if (zed >=0) 
                     {
                        distrTOFw0->Fill(board, alpha, particleWeight);
                     }
                     else 
                     {
                        distrTOFw1->Fill(board, alpha, particleWeight);
                     }
                  }
               }
            }
            
            if (IsHit(T.emcdz(i)))
            {
               if (IsMatch(T.emcsdz(i), T.emcsdphi(i), 2., 2.))
               {   
                  if (T.ecore(i) > 0.25)
                  {
                     if (phi > M_PI) 
                     {
                        if (zed >= 0) distrEMCalePos[T.sect(i)]->Fill(
                           T.pemcy(i), T.pemcz(i), T.ecore(i)*particleWeight);
                        else distrEMCaleNeg[T.sect(i)]->Fill(
                           T.pemcy(i), T.pemcz(i), T.ecore(i)*particleWeight);
                     }
                     else
                     {   
                        if (zed >= 0) distrEMCalwPos[T.sect(i)]->Fill(
                           T.pemcy(i), T.pemcz(i), T.ecore(i)*particleWeight);
                        else distrEMCalwNeg[T.sect(i)]->Fill(
                           T.pemcy(i), T.pemcz(i), T.ecore(i)*particleWeight);
                     }
                  }
               }
            }

         }
      }
   };

   auto PbarCall = [&]()
   {
      ProgressBar pbar = ProgressBar("BLOCK");
      
      while (!isProcessFinished)
      {
         pbar.Print(ncalls/nevents);
         std::this_thread::sleep_for(std::chrono::milliseconds(20));
      }
      pbar.Print(1.);
   };
   
   std::thread pBarThread(PbarCall);
   tp.Process(ProcessMP);
   isProcessFinished = true;
   pBarThread.join();
}

//for ROOT CINT call
void HeatMapper()
{
   for (std::string magf : Par.magfQueue)
   {
      CheckInputFile(Par.realDataDir + Par.runName + "/sum" + magf + ".root");
      
      for (std::string part : Par.partQueue)
      {
         for (std::string auxName : Par.auxNameQueue)
         {
            CheckInputFile(Par.simDataDir + Par.runName + 
               "/SingleTrack/" + part + magf + auxName + ".root");
         }
      }
   }
   if (Par.doUseWeightFunc)
   {
      /*
      for (std::string part : Par.partQueue)
      {
         CheckInputFile("../input/SimAnalysis/Spectra/" + Par.system + "/" + part + ".txt");
      }
      */
   }
   
   if (Par.doReweightAlpha)
   {
      std::string realDataInputFileName = Par.realDataDir + Par.runName + "/sum.root";
      std::string alphaReweightInputFileName = Par.outputDir + 
         Par.runName + "/Heatmaps/all.root";
      
      CheckInputFile(realDataInputFileName);
      if (CheckInputFile(alphaReweightInputFileName, false)) 
      {
         TFile realDataInputFile(realDataInputFileName.c_str());
         TFile alphaReweightInputFile(alphaReweightInputFileName.c_str());

         TH2F *realDataDCe0 = (TH2F *) realDataInputFile.Get("dceast0");
         TH2F *realDataDCe1 = (TH2F *) realDataInputFile.Get("dceast1");
         TH2F *realDataDCw0 = (TH2F *) realDataInputFile.Get("dcwest0");
         TH2F *realDataDCw1 = (TH2F *) realDataInputFile.Get("dcwest1");
         
         TH2F *simDCe0 = (TH2F *) alphaReweightInputFile.Get("unscaled_dceast0");
         TH2F *simDCe1 = (TH2F *) alphaReweightInputFile.Get("unscaled_dceast1");
         TH2F *simDCw0 = (TH2F *) alphaReweightInputFile.Get("unscaled_dcwest0");
         TH2F *simDCw1 = (TH2F *) alphaReweightInputFile.Get("unscaled_dcwest1");

         CutDCDeadAreas(realDataDCe0, &IsDeadDC, 2., 1.);
         CutDCDeadAreas(realDataDCe1, &IsDeadDC, 2., -1.);
         CutDCDeadAreas(realDataDCw0, &IsDeadDC, 1., 1.);
         CutDCDeadAreas(realDataDCw1, &IsDeadDC, 1., -1.);

         CutDCDeadAreas(simDCe0, &IsDeadDC, 2., 1.);
         CutDCDeadAreas(simDCe1, &IsDeadDC, 2., -1.);
         CutDCDeadAreas(simDCw0, &IsDeadDC, 1., 1.);
         CutDCDeadAreas(simDCw1, &IsDeadDC, 1., -1.);

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

         std::string alphaReweightOutputFileName = Par.outputDir + 
            Par.runName + "/alpha_reweight.root";
         TFile alphaReweightOutputFile = TFile(alphaReweightOutputFileName.c_str(), "RECREATE");

         alphaReweightOutputFile.cd();
         
         Par.alphaReweightDCe0->Write();
         Par.alphaReweightDCe1->Write();
         Par.alphaReweightDCw0->Write();
         Par.alphaReweightDCw1->Write();

         //decoupling histograms from the open file directory
         //so they would not be deleted at the end of the scope
         Par.alphaReweightDCe0->SetDirectory(0);
         Par.alphaReweightDCe1->SetDirectory(0);
         Par.alphaReweightDCw0->SetDirectory(0);
         Par.alphaReweightDCw1->SetDirectory(0);         

         alphaReweightOutputFile.Close();
         PrintInfo("File " + alphaReweightOutputFileName + " was written");
      }
      else 
      {
         PrintInfo("alpha reweight is now disabled");
         Par.doReweightAlpha = false;
      }
   }
   
   PrintInfo("Clearing output directory: " + Par.outputDir + 
      Par.runName + "/Heatmaps/");
   system(("mkdir -p " + Par.outputDir + Par.runName + "/Heatmaps").c_str());
   system(("rm -r " + Par.outputDir + Par.runName + "/Heatmaps/*").c_str());

   int num = 1;
   for (std::string part : Par.partQueue)
   {
      for (std::string magf : Par.magfQueue)
      {
         ThrContainer thrContainer;
         for (std::string auxName : Par.auxNameQueue)
         {
            Analyze(&thrContainer, part, magf, auxName, num);
            num++;
         }
         //wiriting the result
         const std::string outputFileName = 
            Par.outputDir + Par.runName + 
            "/Heatmaps/" + part + magf + ".root";
         
         TFile outfile(outputFileName.c_str(), "RECREATE");
         outfile.cd();
         ThrObjHolder.Write();
         outfile.Close();
         PrintInfo("File " + outputFileName + " was written");
      }
   }

   PrintInfo("Merging output files into one");

   system(("hadd -f " + Par.outputDir + 
      Par.runName + "/Heatmaps/all.root " + Par.outputDir + 
      Par.runName + "/Heatmaps/*.root").c_str());
}

int main()
{
   HeatMapper();
   return 0;
}

#endif
