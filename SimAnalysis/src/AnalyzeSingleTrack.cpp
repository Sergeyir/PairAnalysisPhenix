// $SOURCE$
//------------------------------------------------------------------------------------------------
//                     AnalyzeSingleTrack functions realisation
//------------------------------------------------------------------------------------------------
// AnalyzeSingleTrack
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic macro used for evaluation of registration and/or identification of single tracks
 * from simulation output of event-like TTrees to processed histograms 
 * for further efficiency evaluation
 **/
//------------------------------------------------------------------------------------------------

#ifndef ANALYZE_SINGLE_TRACK_CPP
#define ANALYZE_SINGLE_TRACK_CPP

#include "../include/AnalyzeSingleTrack.hpp"

void AnalyzeParticle(ThrContainer *thrContainer, const std::string &particle, 
                     const std::string &magf, const std::string& auxName, const int procNum)
{
   Box box = Box("Parameters of run " + std::to_string(procNum) + 
                 " out of " + std::to_string(Par.partQueue.size()*
                 Par.magfQueue.size()*Par.auxNameQueue.size()));
   
   double nParticles = 0;
   
   const int partIndex = ParticleProperties.iterMap[particle];
   const int partGeantId = ParticleProperties.geantId[partIndex];
   const double partMass = ParticleProperties.mass[partIndex];
   const int partCharge = ParticleProperties.charge[partIndex];

   const std::string simInputFileName = 
      Par.simDataDir + Par.runName + "/Single/" + particle + magf + auxName + ".root";
   const std::string realInputFileName = 
      Par.realDataDir + Par.runName + "/sum" + magf + ".root";

   TFile simInputFile = TFile(simInputFileName.c_str());
   TFile realInputFile = TFile(realInputFileName.c_str());

   const double nevents = static_cast<double>(((TTree *) simInputFile.Get("Tree"))->GetEntries());

   if (nevents <= 0)
   {
      Print("Error: Number of events is equal or less than 0!");
      exit(1);
   }

   TH1F *origPtDistr = (TH1F *) simInputFile.Get("orig_pt");

   const double origPtThreshold = origPtDistr->Integral()/
      static_cast<double>(origPtDistr->GetXaxis()->GetNbins())/2.;

   const double lowPtBound = origPtDistr->GetXaxis()->GetBinLowEdge(
      origPtDistr->FindFirstBinAbove(origPtThreshold));
   const double upPtBound = origPtDistr->GetXaxis()->GetBinUpEdge(
      origPtDistr->FindLastBinAbove(origPtThreshold));

   double eventNorm = 1.;
   if (Par.doUseWeightFunc)
   {
      TH1F *centrDistr = (TH1F *) realInputFile.Get("central_bin");
      
      eventNorm = origPtDistr->Integral(
         origPtDistr->GetXaxis()->FindBin(Par.ptMin),
         origPtDistr->GetXaxis()->FindBin(Par.ptMax))/
         centrDistr->Integral(1, centrDistr->GetXaxis()->GetNbins());

         // this normalization is needed for merging 2 files 
         // with flat pt distribution with different ranges
         eventNorm *= (Par.ptMax - Par.ptMin)/(upPtBound - lowPtBound);
   }

   box.AddEntry("Run name", Par.runName);
   box.AddEntry("Orig particle", particle);
   box.AddEntry("Magnetic field", magf);
   box.AddEntry("Orig pT distribution span", 
                DtoStr(lowPtBound) + " < pT < " + DtoStr(upPtBound));
   box.AddEntry("Use weight function", Par.doUseWeightFunc);
   
   box.AddEntry("Minimum p_T, GeV", Par.ptMin);
   box.AddEntry("Maximum p_T, GeV", Par.ptMax);
   box.AddEntry("Number of events to be analyzed, 1e6", nevents/1e6, 3);

   box.AddEntry("Number of threads", Par.nthreads);
   
   box.Print();

   //tsallis weight function
   std::unique_ptr<TF1> weightFun(
      new TF1("weightFun", 
      "[0]*x^([5]-1)*(1 + ([1] - 1)*(sqrt([4] + x^2)^([5]) - [2])/[3])^(-1./([1]-1.))"));
   weightFun->SetParameters(
      ReadFileIntoArray("../input/SimAnalysis/Spectra/" + Par.system + "/" + particle + ".txt", 6));

   ROOT::EnableImplicitMT(Par.nthreads);
   ROOT::TTreeProcessorMT tp(simInputFileName.c_str());
   
   double ncalls = 0;

   bool isProcessFinished = false;
   auto ProcessMP = [&](TTreeReader &reader)
   {   
      std::shared_ptr<TH2> nPartDistr = thrContainer->nPartDistr.Get();
      
      std::shared_ptr<TH1F> origPtDistr = thrContainer->origPtDistr.Get();

      std::shared_ptr<TH2F> m2TOFwDistr = thrContainer->m2TOFwDistr.Get();
      std::shared_ptr<TH2F> m2TOFeDistr = thrContainer->m2TOFeDistr.Get();
      std::shared_ptr<TH2F> timeTOFeDistr = thrContainer->timeTOFeDistr.Get();
      std::shared_ptr<TH2F> timeTOFwDistr = thrContainer->timeTOFwDistr.Get();

      std::shared_ptr<TH2F> origPtVsRecPtDistr = thrContainer->origPtVsRecPtDistr.Get();

      std::shared_ptr<TH1F> regTOFeDistr = thrContainer->regTOFeDistr.Get();
      std::shared_ptr<TH1F> regTOFwDistr = thrContainer->regTOFwDistr.Get();
      
      std::array<std::shared_ptr<TH1F>, 4> regEMCaleDistr, regEMCalwDistr;
      
      for (int i = 0; i < 4; i++)
      {
         regEMCaleDistr[i] = thrContainer->regEMCaleDistr[i].Get();
         regEMCalwDistr[i] = thrContainer->regEMCalwDistr[i].Get();
      }
      
      EffTreeReader T(reader);
      
      while (reader.Next())
      {
         ncalls += 1.;
         const double origPt = sqrt(pow(T.mom_orig(0), 2) + pow(T.mom_orig(1), 2));
         
         double eventWeight;
         if (Par.doUseWeightFunc) 
         {
            eventWeight = weightFun->Eval(origPt)/eventNorm;
         }
         else eventWeight = 1.;
         
         origPtDistr->Fill(origPt, eventWeight);
         
         nPartDistr->Fill(T.nch()-0.5, origPt, eventWeight);
         
         const double bbcz = T.bbcz();
         if (fabs(bbcz) > 30) continue;

         if (T.nch() <= 0 || T.nch() > 49) continue;

         nParticles += T.nch()*eventWeight;

         for(int i = 0; i < T.nch(); i++)
         {
            const double the0 = T.the0(i);
            const double pt = (T.mom(i))*sin(the0);
            
            if (pt < Par.ptMin || pt > Par.ptMax) continue;
            
            const int charge = T.charge(i);
            if (charge != partCharge) continue;
               
            if (IsQualityCut(T.qual(i))) continue;
   
            const double zed = T.zed(i);
            
            if (abs(zed) > 75 || abs(zed) < 3) continue;
            
            if (!(fabs(the0)<100 &&
               ((bbcz > 0. && ((bbcz - 250.*tan(the0 - TMath::Pi()/2.)) > 2. ||
               (bbcz - 200.*tan(the0 - TMath::Pi()/2.)) < -2.)) ||
               (bbcz < 0. && ((bbcz - 250.*tan(the0 - TMath::Pi()/2.))< -2. ||
               (bbcz - 200.*tan(the0 - TMath::Pi()/2.)) > 2.))))) continue;

            const double alpha = T.alpha(i);
            const double phi = T.phi(i);

            double particleWeight = eventWeight;
            
            double board;
            if (phi>1.5) board = ((3.72402-phi+0.008047*cos(phi+0.87851))/0.01963496);
            else board = ((0.573231 + phi - 0.0046 * cos(phi + 0.05721)) / 0.01963496);

            if (IsDeadDC(phi, zed, board, alpha)) continue;

            double pc1phi = atan2(T.ppc1y(i), T.ppc1x(i));
            if (phi >= 1.5 && pc1phi < 0) pc1phi += M_PI*2.;

            if (IsDeadPC1(phi, T.ppc1z(i), pc1phi)) continue;

            origPtVsRecPtDistr->Fill(pt, origPt, particleWeight);

            if (IsHit(T.tofdz(i)))
            {
               const double beta = T.pltof(i)/T.ttof(i)/29.97;
               const double eLoss = 0.0014*pow(beta, -1.66);

               if (IsMatch(T.tofsdz(i), T.tofsdphi(i)) &&
                  T.etof(i) > eLoss &&
                  !IsBadSlat(T.slat(i)) &&
                  !IsDeadTOFe(zed, T.ptofy(i), T.ptofz(i)))
               {
                  if (T.particle_id(i) == partGeantId && T.primary_id(i) == -999)
                  {
                     regTOFeDistr->Fill(pt, particleWeight);
                  }
                  
                  const double expectedTime = sqrt(pow(partMass/T.mom(i), 2) + 1.)*T.pltof(i)/29.979;
                  const double m2 = pow(T.mom(i), 2)*(pow((T.ttof(i))*29.979/T.pltof(i), 2) - 1.);
                  
                  timeTOFeDistr->Fill(pt, T.ttof(i) - expectedTime, particleWeight);
                  m2TOFeDistr->Fill(pt, m2, particleWeight*0.909);
               }
            }
            else if (IsHit(T.tofwdz(i)))
            {
               double tofwsdphi = GetTOFwsdphi(0, T.mom(i), T.tofwdphi(i), charge, T.striptofw(i));
               double tofwsdz = GetTOFwsdz(0, T.mom(i), T.tofwdz(i), charge, T.striptofw(i));
               
               tofwsdphi = RecalTOFwsdphi(0, T.mom(i), tofwsdphi, charge, T.striptofw(i));
               tofwsdz = RecalTOFwsdz(0, T.mom(i), tofwsdz, charge, T.striptofw(i));
               
               if (IsMatch(tofwsdphi, tofwsdz) && 
                  !IsBadStripTOFw(static_cast<int>(T.striptofw(i))) &&
                  !IsDeadTOFw(zed, board, alpha))
               {
                  //adc cut + tofw tracking efficiency correction
                  const double weightCorrection = 0.7796;
                  
                  double pltofw = T.pltofw(i);
                  const int ichamber = int(T.striptofw(i)/4)%32;

                  if (ichamber<16&&ichamber%2==1) pltofw += 3.358;
                  else if (ichamber>=16&&ichamber%2==0) pltofw += 3.358;
                  double ttofw = T.ttofw(i)-1.12;
                  
                  const double expectedTime = sqrt(
                     pow(partMass/T.mom(i), 2) + 1.)*pltofw/29.979;   
                  timeTOFwDistr->Fill(pt, ttofw-expectedTime, particleWeight*weightCorrection);
                  
                  if (T.particle_id(i) == partGeantId && T.primary_id(i) == -999)
                  {
                     regTOFwDistr->Fill(pt, particleWeight*weightCorrection);
                  }
                  
                  const double m2 = pow(T.mom(i), 2)*
                     (pow((ttofw)*29.979/pltofw, 2) - 1.);
                  
                  m2TOFwDistr->Fill(pt, m2, particleWeight*weightCorrection);
               }
            }

            if (IsHit(T.emcdz(i)))
            {
               if (T.particle_id(i) == partGeantId && T.primary_id(i) == -999 &&
                  IsMatch(T.emcsdz(i), T.emcsdphi(i), 2., 2.) && T.ecore(i) > 0.25 &&
                  !IsDeadEMCal(phi, zed, T.sect(i), T.pemcy(i), T.pemcz(i)))
               {
                  if (phi > 1.5) 
                  {
                     regEMCaleDistr[T.sect(i)]->Fill(pt, particleWeight);
                  }
                  else 
                  {   
                     regEMCalwDistr[T.sect(i)]->Fill(pt, particleWeight);
                  }
               }
            }
         }
      }
   };
   
   auto PbarCall = [&]()
   {
      ProgressBar pbar = ProgressBar("Block");
      while (!isProcessFinished)
      {
         pbar.Print(ncalls/nevents);
         std::this_thread::sleep_for(std::chrono::milliseconds(20));
      }
      pbar.Print(1.);
   };
   
   std::thread pbarThread(PbarCall);
   tp.Process(ProcessMP);
   isProcessFinished = true;
   pbarThread.join();
}

void AnalyzeSingleTrack()
{
   if (Par.doUseWeightFunc)
   {
      for (std::string magf : Par.magfQueue)
      {
         CheckInputFile(Par.realDataDir + Par.runName + "/sum" + magf + ".root");
         
         for (std::string part : Par.partQueue)
         {
            for (std::string auxName : Par.auxNameQueue)
            {
               CheckInputFile(Par.simDataDir + Par.runName + "/Single/" + 
                              part + magf + auxName + ".root");
            }
         }
      }
      for (std::string part : Par.partQueue)
      {
         CheckInputFile("../input/SimAnalysis/Spectra/" + Par.system + "/" + part + ".txt");
      }
   }
   
   system(("mkdir -p " + Par.outputDir + Par.runName).c_str());   
   
   int num = 1;
   for (std::string part : Par.partQueue)
   {
      for (std::string magf : Par.magfQueue)
      {
         ThrContainer thrContainer;
         
         for (std::string auxName : Par.auxNameQueue)
         {
            AnalyzeParticle(&thrContainer, part, magf, auxName, num);
            num++;
         }
         const std::string outputFileName = 
            Par.outputDir + Par.runName + "/" + part + magf + ".root";

         TFile outfile(outputFileName.c_str(), "RECREATE");
         outfile.cd();
         ThrObjHolder.Write();
         outfile.Close();
         PrintInfo("File " + outputFileName + " was written");
      }
      
      system(("hadd -f " + Par.outputDir + Par.runName + "/" + part + 
              ".root " + Par.outputDir + Par.runName + "/" + part + "*").c_str());
   }
}

int main()
{
   AnalyzeSingleTrack();
   return 0;
}

#endif /*ANALYZE_SINGLE_TRACK_CPP*/
