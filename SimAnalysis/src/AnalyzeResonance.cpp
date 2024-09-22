// $SOURCE$
//------------------------------------------------------------------------------------------------
//                       AnalyzeResonance function realisation
//------------------------------------------------------------------------------------------------
// AnalyzeResonance
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic macro used for evaluation of registration of pair of tracks for different methods
 * from simulation output of event-like TTrees to processed histograms 
 * for further track pair from resonance registering correction evaluation
 **/
//------------------------------------------------------------------------------------------------

#ifndef ANALYZE_RESONANCE_CPP
#define ANALYZE_RESONANCE_CPP

#include "../include/AnalyzeResonance.hpp"

//Analyze specific track pair configuration
void AnalyzeConfiguration(ThrContainer *thrContainer, const std::string& daughter1, 
                          const std::string& daughter2, const std::string& magf, 
                          const std::string& auxName, const double pTDeviation, 
                          const int procNum)
{
   Box box = Box("Parameters of run " + std::to_string(procNum) + " out of " + std::to_string(
      Par.daughter1Queue.size()*Par.magfQueue.size()*
      Par.pTDeviationQueue.size()*Par.auxNameQueue.size()));

   
   const std::string decayChannel = 
      ParticleProperties.shortName[ParticleProperties.iterMap[daughter1]] + 
      ParticleProperties.shortName[ParticleProperties.iterMap[daughter2]];

   const std::string inputFileName = "../data/" + Par.runName + "/Resonances/" + 
      Par.origParticle.name + "_" + decayChannel + magf + auxName + ".root";
   
   TFile inputFile = TFile(inputFileName.c_str());

   const double nevents= static_cast<double>(((TTree *) inputFile.Get("Tree"))->GetEntries());
   
   if (nevents <= 0)
   {
      Print("Error: Number of events is less or equal to 0");
      exit(1);
   }
   
   TH1F *distrOrigPT = (TH1F *) inputFile.Get("orig_pt");

   const double origPTThreshold = distrOrigPT->Integral()/
      static_cast<double>(distrOrigPT->GetXaxis()->GetNbins())/2.;

   const double lowPTBound = distrOrigPT->GetXaxis()->GetBinLowEdge(
      distrOrigPT->FindFirstBinAbove(origPTThreshold));
   const double UpPTBound = distrOrigPT->GetXaxis()->GetBinUpEdge(
      distrOrigPT->FindLastBinAbove(origPTThreshold));

   double eventNormWeight = 1.;
   if (Par.doUseWeightFunc)
   {
      eventNormWeight = distrOrigPT->Integral(
         distrOrigPT->GetXaxis()->FindBin(Par.pTMinPair),
         distrOrigPT->GetXaxis()->FindBin(Par.pTMaxPair));

      //this normalization is needed to merge 2 files with flapt pt distributin with different ranges
      eventNormWeight *= (Par.pTMaxPair - Par.pTMinPair)/(UpPTBound - lowPTBound);
   }

   box.AddEntry("Run name", Par.runName);
   box.AddEntry("Orig particle", Par.origParticle.name);
   box.AddEntry("Magnetic field configuration", magf);
   box.AddEntry("Use weight function", Par.doUseWeightFunc);
   
   box.AddEntry("Daughter particle 1", daughter1);
   box.AddEntry("Daughter particle 2", daughter2);
   
   box.AddEntry("Minimum p_T, GeV", Par.pTMin);
   box.AddEntry("Maximum p_T, GeV", Par.pTMax);

   box.AddEntry("Input file orig p_T range, GeV", DtoStr(lowPTBound, 1) + "-" + DtoStr(UpPTBound, 1));
   
   box.AddEntry("Number of events to be analyzed, 1e6", nevents/1e6, 3);

   box.AddEntry("pT deviation", DtoStr(pTDeviation, 3));
   
   box.AddEntry("Number of threads", Par.nThreads);
   
   box.Print();

   TF1 weightFunc = TF1("weightFunc", 
      "[0]*x^([5]-1)*(1 + ([1] - 1)*(sqrt([4] + x^2)^([5]) - [2])/[3])^(-1./([1]-1.))");
   weightFunc.SetParameters(
      ReadFileIntoArray("../input/Spectra/" + Par.system + "/" + Par.origParticle.name + ".txt", 6));
   
   ROOT::EnableImplicitMT(Par.nThreads);
   ROOT::TTreeProcessorMT tp(inputFileName.c_str());
   
   double ncalls = 0.;

   //acceptance systematic uncertainties of detectors
   double sysDCe0, sysDCe1, sysDCw0, sysDCw1,
      sysPC1e, sysPC1w, sysPC2, sysPC3e, sysPC3w,
      sysTOFe0, sysTOFe1, sysTOFw0, sysTOFw1;
   std::array<double, 4> sysEMCalePos, sysEMCaleNeg, sysEMCalwPos, sysEMCalwNeg;

   ReadFile("../input/Systematics/" + Par.runName + "/acceptance.txt",
      sysDCe0, sysDCe1, sysDCw0, sysDCw1,
      sysPC1e, sysPC1w, sysPC2, sysPC3e, sysPC3w,
      sysTOFe0, sysTOFe1, sysTOFw0, sysTOFw1,
      sysEMCalePos[0], sysEMCaleNeg[0], sysEMCalePos[1], sysEMCaleNeg[1],
      sysEMCalePos[2], sysEMCaleNeg[2], sysEMCalePos[3], sysEMCaleNeg[3],
      sysEMCalwPos[0], sysEMCalwNeg[0], sysEMCalwPos[1], sysEMCalwNeg[1],
      sysEMCalwPos[2], sysEMCalwNeg[2], sysEMCalwPos[3], sysEMCalwNeg[3]);
   
   bool isProcessFinished = false;
   auto ProcessMP = [&](TTreeReader &reader)
   {
      //arrays and constants for positives and negative particles and their tracks
      ParticleContainer pTrack, nTrack;

      //arrays for keeping weights for varied acceptancies
      ParticleContainer pTrackDecreasedAcceptance, pTrackIncreasedAcceptance, 
                        nTrackDecreasedAcceptance, nTrackIncreasedAcceptance;
      ParticleContainer pTrackDecreasedM2Eff, pTrackIncreasedM2Eff, 
                        nTrackDecreasedM2Eff, nTrackIncreasedM2Eff;

      //arrays for keeping weights for varied /epsilon_{m^2}
      
      //setting daughter particles
      pTrack.iter = ParticleProperties.iterMap[daughter1];
      nTrack.iter = ParticleProperties.iterMap[daughter2];
      
      pTrack.origId = abs(ParticleProperties.id[pTrack.iter]);
      nTrack.origId = abs(ParticleProperties.id[nTrack.iter]);

      idContainer id1, id2;

      id1.origId = pTrack.origId;
      id2.origId = nTrack.origId;

      pTrack.geantId = ParticleProperties.geantId[pTrack.iter];
      nTrack.geantId = ParticleProperties.geantId[nTrack.iter];

      pTrack.mass = ParticleProperties.mass[pTrack.iter];
      nTrack.mass = ParticleProperties.mass[nTrack.iter];

      //setting embedding
      pTrack.embDCPC1 = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
         ParticleProperties.absName[pTrack.iter] + "_dc_pc1.txt", 
         CentralityContainer.size));
      pTrack.embPC2 = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
         ParticleProperties.absName[pTrack.iter] + "_pc2.txt", 
         CentralityContainer.size));
      pTrack.embPC3 = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
         ParticleProperties.absName[pTrack.iter] + "_pc3.txt", 
         CentralityContainer.size));
      pTrack.embTOFe = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
         ParticleProperties.absName[pTrack.iter] + "_tofe.txt", 
         CentralityContainer.size));
      pTrack.embTOFw = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
         ParticleProperties.absName[pTrack.iter] + "_tofw.txt", 
         CentralityContainer.size));

      nTrack.embDCPC1 = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
         ParticleProperties.absName[nTrack.iter] + "_dc_pc1.txt", 
         CentralityContainer.size));
      nTrack.embPC2 = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
         ParticleProperties.absName[nTrack.iter] + "_pc2.txt", 
         CentralityContainer.size));
      nTrack.embPC3 = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
         ParticleProperties.absName[nTrack.iter] + "_pc3.txt", 
         CentralityContainer.size));
      nTrack.embTOFe = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
         ParticleProperties.absName[nTrack.iter] + "_tofe.txt", 
         CentralityContainer.size));
      nTrack.embTOFw = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
         ParticleProperties.absName[nTrack.iter] + "_tofw.txt", 
         CentralityContainer.size));
      
      for (int i = 0; i < 4; i++)
      {
         pTrack.embEMCale[i] = std::unique_ptr<double>(
            ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
            ParticleProperties.absName[pTrack.iter] + "_emcale" + std::to_string(i) + ".txt", 
            CentralityContainer.size));
         pTrack.embEMCalw[i] = std::unique_ptr<double>(
            ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
            ParticleProperties.absName[pTrack.iter] + "_emcalw" + std::to_string(i) + ".txt", 
            CentralityContainer.size));
         
         nTrack.embEMCale[i] = std::unique_ptr<double>(
            ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
            ParticleProperties.absName[nTrack.iter] + "_emcale" + std::to_string(i) + ".txt", 
            CentralityContainer.size));
         nTrack.embEMCalw[i] = std::unique_ptr<double>(
            ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
            ParticleProperties.absName[nTrack.iter] + "_emcalw" + std::to_string(i) + ".txt", 
            CentralityContainer.size));
      }

      pTrack.m2EffEMCale = std::unique_ptr<double>(
         ReadFileIntoArray("../input/M2Eff/" + Par.runName + "/" +
         ParticleProperties.name[pTrack.iter] + "_EMCale" + magf + ".txt", 2));
      pTrack.m2EffEMCalw = std::unique_ptr<double>(
         ReadFileIntoArray("../input/M2Eff/" + Par.runName + "/" +
         ParticleProperties.name[pTrack.iter] + "_EMCalw" + magf + ".txt", 4));
      
      nTrack.m2EffEMCale = std::unique_ptr<double>(
         ReadFileIntoArray("../input/M2Eff/" + Par.runName + "/" +
         ParticleProperties.name[nTrack.iter] + "_EMCale" + magf + ".txt", 2));
      nTrack.m2EffEMCalw = std::unique_ptr<double>(
         ReadFileIntoArray("../input/M2Eff/" + Par.runName + "/" +
         ParticleProperties.name[nTrack.iter] + "_EMCalw" + magf + ".txt", 4));

      pTrack.m2EffSysEMCale = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Systematics/" + Par.runName + "/m2eff_" +
         ParticleProperties.name[pTrack.iter] + "_EMCale" + magf + ".txt", 2));
      pTrack.m2EffSysEMCalw = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Systematics/" + Par.runName + "/m2eff_" +
         ParticleProperties.name[pTrack.iter] + "_EMCalw" + magf + ".txt", 4));
      
      nTrack.m2EffSysEMCale = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Systematics/" + Par.runName + "/m2eff_" +
         ParticleProperties.name[nTrack.iter] + "_EMCale" + magf + ".txt", 2));
      nTrack.m2EffSysEMCalw = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Systematics/" + Par.runName + "/m2eff_" +
         ParticleProperties.name[nTrack.iter] + "_EMCalw" + magf + ".txt", 4));

      std::shared_ptr<TH1F> origPtDistr = thrContainer->origPtDistr.Get();
      
      std::shared_ptr<TH2F> distrOrigPTVsRecPT = thrContainer->distrOrigPTVsRecPT.Get();
      
      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> distrInvMNoPID;
      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> distrInvMNoPIDDecreasedAcceptance;
      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> distrInvMNoPIDIncreasedAcceptance;

      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> distrInvM1PID;
      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> distrInvM1PIDDecreasedAcceptance;
      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> distrInvM1PIDIncreasedAcceptance;

      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> distrInvM2PID;
      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> distrInvM2PIDDecreasedAcceptance;
      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> distrInvM2PIDIncreasedAcceptance;
      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> distrInvM2PIDDecreasedM2Eff;
      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> distrInvM2PIDIncreasedM2Eff;

      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> distrInvMTOF2PID;
      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> distrInvMTOF2PIDDecreasedAcceptance;
      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> distrInvMTOF2PIDIncreasedAcceptance;

      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> distrInvMEMC2PID;
      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> distrInvMEMC2PIDDecreasedAcceptance;
      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> distrInvMEMC2PIDIncreasedAcceptance;
      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> distrInvMEMC2PIDDecreasedM2Eff;
      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> distrInvMEMC2PIDIncreasedM2Eff;

      for (int i = 0; i < CentralityContainer.size; i++)
      {
         distrInvMNoPID[i] = thrContainer->distrInvMNoPID[i]->Get();
         distrInvMNoPIDDecreasedAcceptance[i] = thrContainer->distrInvMNoPIDDecreasedAcceptance[i]->Get();
         distrInvMNoPIDIncreasedAcceptance[i] = thrContainer->distrInvMNoPIDIncreasedAcceptance[i]->Get();

         distrInvM1PID[i] = thrContainer->distrInvM1PID[i]->Get();
         distrInvM1PIDDecreasedAcceptance[i] = thrContainer->distrInvM1PIDDecreasedAcceptance[i]->Get();
         distrInvM1PIDIncreasedAcceptance[i] = thrContainer->distrInvM1PIDIncreasedAcceptance[i]->Get();

         distrInvM2PID[i] = thrContainer->distrInvM2PID[i]->Get();
         distrInvM2PIDDecreasedAcceptance[i] = thrContainer->distrInvM2PIDDecreasedAcceptance[i]->Get();
         distrInvM2PIDIncreasedAcceptance[i] = thrContainer->distrInvM2PIDIncreasedAcceptance[i]->Get();
         distrInvM2PIDDecreasedM2Eff[i] = thrContainer->distrInvM2PIDDecreasedM2Eff[i]->Get();
         distrInvM2PIDIncreasedM2Eff[i] = thrContainer->distrInvM2PIDIncreasedM2Eff[i]->Get();

         distrInvMTOF2PID[i] = thrContainer->distrInvMTOF2PID[i]->Get();
         distrInvMTOF2PIDDecreasedAcceptance[i] = thrContainer->distrInvMTOF2PIDDecreasedAcceptance[i]->Get();
         distrInvMTOF2PIDIncreasedAcceptance[i] = thrContainer->distrInvMTOF2PIDIncreasedAcceptance[i]->Get();

         distrInvMEMC2PID[i] = thrContainer->distrInvMEMC2PID[i]->Get();
         distrInvMEMC2PIDDecreasedAcceptance[i] = thrContainer->distrInvMEMC2PIDDecreasedAcceptance[i]->Get();
         distrInvMEMC2PIDIncreasedAcceptance[i] = thrContainer->distrInvMEMC2PIDIncreasedAcceptance[i]->Get();
         distrInvMEMC2PIDDecreasedM2Eff[i] = thrContainer->distrInvMEMC2PIDDecreasedM2Eff[i]->Get();
         distrInvMEMC2PIDIncreasedM2Eff[i] = thrContainer->distrInvMEMC2PIDIncreasedM2Eff[i]->Get();
      }
      
      EffTreeReader T(reader);
      
      while (reader.Next())
      {   
         ncalls += 1.;
         const double origPt = sqrt(pow(T.mom_orig(0), 2) + pow(T.mom_orig(1), 2))*pTDeviation;
         
         double eventWeight = 1.;
         
         if (Par.doUseWeightFunc) eventWeight = weightFunc.Eval(origPT)/eventNormWeight*4e11;
         
         origPtDistr->Fill(origPt, eventWeight);
         
         if (T.nch() > 49 || T.nch() <= 0) continue;
         
         const double bbcz = T.bbcz();
         if (fabs(bbcz) > 30) continue;

         pTrack.size = 0;
         nTrack.size = 0;

         for(int i = 0; i < T.nch(); i++)
         {   
            const double the0 = T.the0(i);
            const double pT = (T.mom(i))*sin(the0)*pTDeviation;

            if (pT < Par.pTMin || pT > Par.pTMax) continue;
            if (IsQualityCut(T.qual(i))) continue;
            
            const int charge = T.charge(i);

            if (charge != -1 && charge != 1) continue;
   
            if (!(fabs(the0)<100 &&
               ((bbcz > 0 && ((bbcz - 250*tan(the0 - 3.1416/2)) > 2 ||
               (bbcz - 200*tan(the0 - 3.1416/2)) < -2)) ||
               (bbcz < 0 && ((bbcz - 250*tan(the0 - 3.1416/2))< -2 ||
               (bbcz - 200*tan(the0 - 3.1416/2)) > 2))))) continue;

            const double zed = T.zed(i);
            if (abs(zed) > 75 && abs(zed) < 3) continue;
            
            const double alpha = T.alpha(i);
            const double phi = T.phi(i);
            double board = -9999;
            
            if (phi>1.5) board = ((3.72402-phi+0.008047*cos(phi+0.87851))/0.01963496);
            else board = ((0.573231 + phi - 0.0046 * cos(phi + 0.05721)) / 0.01963496);

            if (IsDeadDC(phi, zed, board, alpha)) continue;

            double pc1phi = atan2(T.ppc1y(i), T.ppc1x(i));
            if (phi >= 1.5 && pc1phi < 0) pc1phi += M_PI*2.;

            if (IsDeadPC1(phi, T.ppc1z(i), pc1phi)) continue;

            pTrack.ResetTrack(pTrack.size);
            nTrack.ResetTrack(nTrack.size);
            pTrackDecreasedAcceptance.ResetTrack(pTrack.size);
            nTrackDecreasedAcceptance.ResetTrack(nTrack.size);
            pTrackIncreasedAcceptance.ResetTrack(pTrack.size);
            nTrackIncreasedAcceptance.ResetTrack(nTrack.size);
            pTrackDecreasedM2Eff.ResetTrack(pTrack.size);
            nTrackDecreasedM2Eff.ResetTrack(nTrack.size);
            pTrackIncreasedM2Eff.ResetTrack(pTrack.size);
            nTrackIncreasedM2Eff.ResetTrack(nTrack.size);

            if (IsHit(T.tofdz(i)))
            {
               const double beta = T.pltof(i)/T.ttof(i)/29.97;
               const double eloss = 0.0014*pow(beta, -1.66);
               
               if (IsMatch(T.tofsdz(i), T.tofsdphi(i)) && 
                  T.etof(i) > eloss &&
                  !IsBadSlat(T.slat(i)) &&
                  !IsDeadTOFe(zed, T.ptofy(i), T.ptofz(i)))
               {
                  const double m2 = pow(T.mom(i), 2)*(pow((T.ttof(i))*29.979/T.pltof(i), 2) - 1.);
                  
                  const double id = GetTOFePID(pT, charge, m2);
                  
                  //variation of acceptance
                  double acceptanceVariation;
                  if (zed >= 0) acceptanceVariation = ErrPropagation(sysTOFe0, sysDCe0, sysPC1e);
                  else acceptanceVariation = ErrPropagation(sysTOFe1, sysDCe1, sysPC1e);
   
                  if (charge == 1)
                  {
                     pTrack.idTOF[pTrack.size] = id;
                     
                     for (int j = 0; j < CentralityContainer.size; j++)
                     {
                        pTrack.weightTOF[j][pTrack.size] = 
                           pTrack.embTOFe.get()[j]/pTrack.embDCPC1.get()[j];
                        pTrackDecreasedAcceptance.weightTOF[j][pTrack.size] = 
                           pTrack.weightTOF[j][pTrack.size]*(1. - acceptanceVariation);
                        pTrackIncreasedAcceptance.weightTOF[j][pTrack.size] = 
                           Minimum(1., pTrack.weightTOF[j][pTrack.size]*(1.+acceptanceVariation));

                        if (id == pTrack.origId) 
                        {
                           pTrack.weightIdTOF[j][pTrack.size] = 
                              pTrack.weightTOF[j][pTrack.size];
                           pTrackDecreasedAcceptance.weightIdTOF[j][pTrack.size] = 
                              pTrack.weightIdTOF[j][pTrack.size]*(1. - acceptanceVariation);
                           pTrackIncreasedAcceptance.weightIdTOF[j][pTrack.size] = 
                              Minimum(1., pTrack.weightIdTOF[j][pTrack.size]*(1.+acceptanceVariation));
                        }   
                     }
                  }
                  else
                  {
                     nTrack.idTOF[nTrack.size] = id;
                     
                     for (int j = 0; j < CentralityContainer.size; j++)
                     {
                        nTrack.weightTOF[j][nTrack.size] = 
                           nTrack.embTOFe.get()[j]/nTrack.embDCPC1.get()[j];
                        nTrackDecreasedAcceptance.weightTOF[j][nTrack.size] = 
                           nTrack.weightTOF[j][nTrack.size]*(1. - acceptanceVariation);
                        nTrackIncreasedAcceptance.weightTOF[j][nTrack.size] = 
                           Minimum(1., nTrack.weightTOF[j][nTrack.size]*(1.+acceptanceVariation));
                        
                        if (id == nTrack.origId)
                        {
                           nTrack.weightIdTOF[j][nTrack.size] = 
                              nTrack.weightTOF[j][nTrack.size];
                           nTrackDecreasedAcceptance.weightIdTOF[j][nTrack.size] = 
                              nTrack.weightIdTOF[j][nTrack.size]*(1. - acceptanceVariation);
                           nTrackIncreasedAcceptance.weightIdTOF[j][nTrack.size] = Minimum(1.,
                              nTrack.weightIdTOF[j][nTrack.size]*(1.+acceptanceVariation));
                        }
                     }
                  }
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
                  double pltofw = T.pltofw(i);
                  const int ichamber = int(T.striptofw(i)/4)%32;

                  if (ichamber < 16 && ichamber % 2 == 1) pltofw += 3.358;
                  else if (ichamber >= 16 && ichamber % 2 == 0) pltofw += 3.358;
                  double ttofw = T.ttofw(i)-1.12;

                  const double m2 = pow(T.mom(i), 2)*
                     (pow((ttofw)*29.979/pltofw, 2) - 1.);
                  
                  const double id = GetTOFePID(pT, charge, m2);

                  double acceptanceVariation;
                  if (zed >= 0) acceptanceVariation = ErrPropagation(sysTOFw0, sysDCw0, sysPC1w);
                  else acceptanceVariation = ErrPropagation(sysTOFw1, sysDCw1, sysPC1w);
                  
                  if (charge == 1)
                  {
                     pTrack.idTOF[pTrack.size] = id;
                     
                     for (unsigned int j = 0; j < CentralityContainer.size; j++)
                     {
                        pTrack.weightTOF[j][pTrack.size] = 
                           Par.correctionTOFw*pTrack.embTOFw.get()[j]/pTrack.embDCPC1.get()[j];
                        pTrackDecreasedAcceptance.weightTOF[j][pTrack.size] = 
                           pTrack.weightTOF[j][pTrack.size]*(1. - acceptanceVariation);
                        pTrackIncreasedAcceptance.weightTOF[j][pTrack.size] = Minimum(1.,
                           pTrack.weightTOF[j][pTrack.size]*(1. + acceptanceVariation));

                        if (id == pTrack.origId) 
                        {
                           pTrack.weightIdTOF[j][pTrack.size] = 
                               pTrack.weightTOF[j][pTrack.size];
                           pTrackDecreasedAcceptance.weightIdTOF[j][pTrack.size] = 
                               pTrack.weightTOF[j][pTrack.size]*(1. - acceptanceVariation);
                           pTrackIncreasedAcceptance.weightIdTOF[j][pTrack.size] = Minimum(1.,
                              pTrack.weightTOF[j][pTrack.size]*(1.+acceptanceVariation));
                        }
                     }
                  }
                  else
                  {
                     nTrack.idTOF[nTrack.size] = id;
                     
                     for (unsigned int j = 0; j < CentralityContainer.size; j++)
                     {
                        nTrack.weightTOF[j][nTrack.size] = 
                           Par.correctionTOFw*nTrack.embTOFw.get()[j]/nTrack.embDCPC1.get()[j];
                        nTrackDecreasedAcceptance.weightTOF[j][nTrack.size] = 
                           nTrack.weightTOF[j][nTrack.size]*(1. - acceptanceVariation);
                        nTrackIncreasedAcceptance.weightTOF[j][nTrack.size] = Minimum(1.,
                           nTrack.weightTOF[j][nTrack.size]*(1.+acceptanceVariation));

                        if (id == nTrack.origId) 
                        {
                           nTrack.weightIdTOF[j][nTrack.size] = 
                              nTrack.weightTOF[j][nTrack.size];
                           nTrackDecreasedAcceptance.weightIdTOF[j][nTrack.size] = 
                              nTrack.weightIdTOF[j][nTrack.size]*(1. - acceptanceVariation);
                           nTrackIncreasedAcceptance.weightIdTOF[j][nTrack.size] = Minimum(1.,
                              nTrack.weightIdTOF[j][nTrack.size]*(1.+acceptanceVariation));
                        }
                     }
                  }
               }
            }

            if (IsHit(T.emcdz(i)))
            {
               if (T.ecore(i) > 0.25 && 
                  IsMatch(T.emcsdz(i), T.emcsdphi(i)) &&
                  !IsDeadEMCal(phi, zed, T.sect(i), T.pemcy(i), T.pemcz(i)) &&
                  T.sect(i) >= 0 && T.sect(i) < 4)
               {   
                  double acceptanceVariation;
                  double m2EffVariation = 0.;

                  if (charge == 1) pTrack.idEMC[pTrack.size] = PartId.noPID;
                  else nTrack.idEMC[nTrack.size] = PartId.noPID;

                  if (phi > 1.5)
                  {
                     if (zed >= 0) acceptanceVariation = ErrPropagation(
                        sysEMCalePos[T.sect(i)], sysDCe0, sysPC1e);
                     else acceptanceVariation = ErrPropagation(
                        sysEMCalwPos[T.sect(i)], sysDCe1, sysPC1e);
                     
                     if (charge == 1)
                     {
                        for (unsigned int j = 0; j < CentralityContainer.size; j++)
                        {
                           pTrack.weightEMC[j][pTrack.size] = 
                              pTrack.embEMCale[T.sect(i)].get()[j]/pTrack.embDCPC1.get()[j];
                        }
         
                        if (T.particle_id(i) == pTrack.geantId && T.sect(i) > 1)
                        {
                           pTrack.idEMC[pTrack.size] = pTrack.origId;
                           for (unsigned int j = 0; j < CentralityContainer.size; j++)
                           {
                              pTrack.weightIdEMC[j][pTrack.size] = 
                                 pTrack.weightEMC[j][pTrack.size]*
                                 pTrack.m2EffEMCale.get()[T.sect(i)-2]*
                                 GetEMCalId(pT, pTrack.idEMC[pTrack.size], charge, phi, T.sect(i));
                              m2EffVariation = pTrack.m2EffSysEMCale.get()[T.sect(i)-2];
                           }
                        }
                     }
                     else
                     {   
                        for (unsigned int j = 0; j < CentralityContainer.size; j++)
                        {
                           nTrack.weightEMC[j][nTrack.size] = 
                              nTrack.embEMCale[T.sect(i)].get()[j]/nTrack.embDCPC1.get()[j];
                        }
                           
                        if (T.particle_id(i) == nTrack.geantId && T.sect(i) > 1)
                        {
                           nTrack.idEMC[nTrack.size] = nTrack.origId;
                           for (unsigned int j = 0; j < CentralityContainer.size; j++)
                           {
                              nTrack.weightIdEMC[j][nTrack.size] = 
                                 nTrack.weightEMC[j][nTrack.size]*
                                 nTrack.m2EffEMCale.get()[T.sect(i)-2]*
                                 GetEMCalId(pT, nTrack.idEMC[nTrack.size], charge, phi, T.sect(i));
                              m2EffVariation = nTrack.m2EffSysEMCale.get()[T.sect(i)-2];
                           }
                        }
                     }
                     if (pT > 1.1 && (T.particle_id(i) == 11 || T.particle_id(i) == 12)) 
                     {
                        switch(T.sect(i))
                        {
                           case 0: 
                              m2EffVariation *= 1.1227;
                              break;
                           case 1: 
                              m2EffVariation *= 1.0639;
                              break;
                        }
                     }
                  }
                  else 
                  {
                     if (zed >= 0) acceptanceVariation = ErrPropagation(
                        sysEMCalwPos[T.sect(i)], sysDCw0, sysPC1w);
                     else acceptanceVariation = ErrPropagation(
                        sysEMCalwNeg[T.sect(i)], sysDCw1, sysPC1w);
                     
                     if (charge == 1)
                     {
                        for (unsigned int j = 0; j < CentralityContainer.size; j++)
                        {
                           pTrack.weightEMC[j][pTrack.size] = 
                              pTrack.embEMCalw[T.sect(i)].get()[j]/pTrack.embDCPC1.get()[j];
                        }
                        
                        if (T.particle_id(i) == pTrack.geantId)
                        {
                           pTrack.idEMC[pTrack.size] = pTrack.origId;
                           for (unsigned int j = 0; j < CentralityContainer.size; j++)
                           {
                              pTrack.weightIdEMC[j][pTrack.size] = 
                                 pTrack.weightEMC[j][pTrack.size]*
                                 pTrack.m2EffEMCalw.get()[T.sect(i)]*
                                 GetEMCalId(pT, pTrack.idEMC[pTrack.size], charge, phi, T.sect(i));
                              m2EffVariation = pTrack.m2EffSysEMCalw.get()[T.sect(i)];
                           }
                        }
                     }
                     else
                     {
                        for (unsigned int j = 0; j < CentralityContainer.size; j++)
                        {
                           nTrack.weightEMC[j][nTrack.size] = 
                              nTrack.embEMCalw[T.sect(i)].get()[j]/nTrack.embDCPC1.get()[j];
                        }
                        
                        if (T.particle_id(i) == nTrack.geantId)
                        {
                           nTrack.idEMC[nTrack.size] = nTrack.origId;
                           for (unsigned int j = 0; j < CentralityContainer.size; j++)
                           {
                              nTrack.weightIdEMC[j][nTrack.size] = 
                                 nTrack.weightEMC[j][nTrack.size]*
                                 nTrack.m2EffEMCalw.get()[T.sect(i)]*
                                 GetEMCalId(pT, nTrack.idEMC[nTrack.size], charge, phi, T.sect(i));
                              m2EffVariation = nTrack.m2EffSysEMCalw.get()[T.sect(i)];
                           }
                        }
                     }
                     if (pT > 1.1 && (T.particle_id(i) == 11 || T.particle_id(i) == 12)) 
                     {
                        switch(T.sect(i))
                        {
                           case 0: 
                              m2EffVariation *= 1.2885;
                              break;
                           case 1: 
                              m2EffVariation *= 1.396;
                              break;
                        }
                     }
                  }
                  
                  if (charge == 1)
                  {
                     for (unsigned int j = 0; j < CentralityContainer.size; j++)
                     {
                        pTrackDecreasedAcceptance.weightEMC[j][pTrack.size] = 
                           pTrack.weightEMC[j][pTrack.size]*(1. - acceptanceVariation);
                        pTrackIncreasedAcceptance.weightEMC[j][pTrack.size] = Minimum(1.,
                           pTrack.weightEMC[j][pTrack.size]*(1. + acceptanceVariation));

                        pTrackDecreasedAcceptance.weightIdEMC[j][pTrack.size] = 
                           pTrack.weightIdEMC[j][pTrack.size]*(1. - acceptanceVariation);
                        pTrackIncreasedAcceptance.weightIdEMC[j][pTrack.size] = Minimum(1.,
                           pTrack.weightIdEMC[j][pTrack.size]*(1. + acceptanceVariation));

                        pTrackDecreasedM2Eff.weightIdEMC[j][pTrack.size] = 
                           pTrack.weightIdEMC[j][pTrack.size]*(1. - m2EffVariation);
                        pTrackIncreasedM2Eff.weightIdEMC[j][pTrack.size] = Minimum(1.,
                           pTrack.weightIdEMC[j][pTrack.size]*(1. + m2EffVariation));
                     }
                  }
                  else
                  {
                     for (unsigned int j = 0; j < CentralityContainer.size; j++)
                     {
                        nTrackDecreasedAcceptance.weightEMC[j][nTrack.size] = 
                           nTrack.weightEMC[j][nTrack.size]*(1. - acceptanceVariation);
                        nTrackIncreasedAcceptance.weightEMC[j][nTrack.size] = Minimum(1.,
                           nTrack.weightEMC[j][nTrack.size]*(1. + acceptanceVariation));

                        nTrackDecreasedAcceptance.weightIdEMC[j][nTrack.size] = 
                           nTrack.weightIdEMC[j][nTrack.size]*(1. - acceptanceVariation);
                        nTrackIncreasedAcceptance.weightIdEMC[j][nTrack.size] = Minimum(1.,
                           nTrack.weightIdEMC[j][nTrack.size]*(1. + acceptanceVariation));

                        nTrackDecreasedM2Eff.weightIdEMC[j][nTrack.size] = 
                           nTrack.weightIdEMC[j][nTrack.size]*(1. - m2EffVariation);
                        nTrackIncreasedM2Eff.weightIdEMC[j][nTrack.size] = Minimum(1.,
                           nTrack.weightIdEMC[j][nTrack.size]*(1. + m2EffVariation));
                     }
                  }
               }
            }

            if (IsMatch(T.pc2sdz(i), T.pc2sdphi(i)))
            {
               if (IsMatch(T.pc2sdz(i), T.pc2sdphi(i)))
               {
                  const double pc2z = T.ppc2z(i) - T.pc2dz(i);
                  const double pc2phi = atan2(T.ppc2y(i), T.ppc2x(i)) - T.pc2dphi(i);
                  
                  if (!IsDeadPC2(pc2z, pc2phi))
                  {
                     double acceptanceVariation;

                     if (zed >= 0) acceptanceVariation = ErrPropagation(sysPC2, sysDCe0, sysPC1e);
                     else acceptanceVariation = ErrPropagation(sysPC2, sysDCe1, sysPC1e);
                     
                     if (charge == 1)
                     {
                        pTrack.idPC2[pTrack.size] = PartId.noPID;
   
                        for (int j = 0; j < CentralityContainer.size; j++)
                        {
                           pTrack.weightPC2[j][pTrack.size] = 
                              pTrack.embPC2.get()[j]/pTrack.embDCPC1.get()[j];
                           
                           pTrackDecreasedAcceptance.weightPC2[j][pTrack.size] = 
                              pTrack.weightPC2[j][pTrack.size]*(1. - acceptanceVariation);
                           pTrackIncreasedAcceptance.weightPC2[j][pTrack.size] = Minimum(1.,
                              pTrack.weightPC2[j][pTrack.size]*(1. + acceptanceVariation));
                        }
                     }
                     else
                     {
                        nTrack.idPC2[nTrack.size] = PartId.noPID;
                        
                        for (int j = 0; j < CentralityContainer.size; j++)
                        {
                           nTrack.weightPC2[j][nTrack.size] = 
                              nTrack.embPC2.get()[j]/nTrack.embDCPC1.get()[j];

                           nTrackDecreasedAcceptance.weightPC2[j][nTrack.size] = 
                              nTrack.weightPC2[j][nTrack.size]*(1. - acceptanceVariation);
                           nTrackIncreasedAcceptance.weightPC2[j][nTrack.size] = Minimum(1.,
                              nTrack.weightPC2[j][nTrack.size]*(1. + acceptanceVariation));
                        }
                     }
                  }
               }
            }

            if (IsHit(T.pc3dz(i)))
            {
               if (IsMatch(T.pc3sdz(i), T.pc3sdphi(i)))
               {
                  const double pc3z = T.ppc3z(i) - T.pc3dz(i);
                  double pc3phi = atan2(T.ppc3y(i), T.ppc3x(i) - T.pc3dphi(i));

                  if (pc3phi < 0) pc3phi += 2.*M_PI;
                  
                  if (!IsDeadPC3(phi, pc3z, pc3phi))
                  {
                     double acceptanceVariation;

                     if (pc3phi < 0) 
                     {
                        if (zed >= 0) acceptanceVariation = ErrPropagation(sysPC3w, sysDCw0, sysPC1w);
                        else acceptanceVariation = ErrPropagation(sysPC3w, sysDCw1, sysPC1w);
                     }
                     else
                     {
                        if (zed >= 0) acceptanceVariation = ErrPropagation(sysPC3e, sysDCe0, sysPC1e);
                        else acceptanceVariation = ErrPropagation(sysPC3e, sysDCe1, sysPC1e);
                     }
                     
                     if (charge == 1)
                     {
                        pTrack.idPC3[pTrack.size] = PartId.noPID;
                        
                        for (int j = 0; j < CentralityContainer.size; j++)
                        {
                           pTrack.weightPC3[j][pTrack.size] = 
                              pTrack.embPC3.get()[j]/pTrack.embDCPC1.get()[j];

                           pTrackDecreasedAcceptance.weightPC3[j][pTrack.size] = 
                              pTrack.weightPC3[j][pTrack.size]*(1. - acceptanceVariation); 
                           pTrackIncreasedAcceptance.weightPC3[j][pTrack.size] = Minimum(1.,
                              pTrack.weightPC3[j][pTrack.size]*(1. + acceptanceVariation)); 
                        }
                     }
                     else
                     {
                        nTrack.idPC3[nTrack.size] = PartId.noPID;
                        
                        for (int j = 0; j < CentralityContainer.size; j++)
                        {
                           nTrack.weightPC3[j][nTrack.size] = 
                              nTrack.embPC3.get()[j]/nTrack.embDCPC1.get()[j];
                           
                           nTrackDecreasedAcceptance.weightPC3[j][nTrack.size] = 
                              nTrack.weightPC3[j][nTrack.size]*(1. - acceptanceVariation); 
                           nTrackIncreasedAcceptance.weightPC3[j][nTrack.size] = Minimum(1.,
                              nTrack.weightPC3[j][nTrack.size]*(1. + acceptanceVariation)); 
                        }
                     }
                  }
               }
            }
   
            //At least 1 detector is required
            if (charge == 1)
            {
               if (pTrack.idPC2[pTrack.size] == PartId.junk &&
                     pTrack.idPC3[pTrack.size] == PartId.junk &&
                     pTrack.idTOF[pTrack.size] == PartId.junk &&
                     pTrack.idEMC[pTrack.size] == PartId.junk) 
               {
                  continue;
               }
               else 
               {
                  pTrack.index[pTrack.size] = i;
                  pTrack.size += 1;
               }
            }
            else if (charge == -1)
            {
               if (nTrack.idPC2[nTrack.size] == PartId.junk &&
                     nTrack.idPC3[nTrack.size] == PartId.junk &&
                     nTrack.idTOF[nTrack.size] == PartId.junk &&
                     nTrack.idEMC[nTrack.size] == PartId.junk) 
               {
                  continue;
               }
               else 
               {
                  nTrack.index[nTrack.size] = i;
                  nTrack.size += 1;
               }
            }
         }

         for (int i = 0; i < pTrack.size; i++)
         {
            for (int j = 0; j < nTrack.size; j++)
            {
               const int ppc = pTrack.index[i]; //positive particle counter
               const int npc = nTrack.index[j]; //negative particle counter
               
               pTrack.mom[0] = T.mom(ppc)*sin(T.the0(ppc))*cos(T.phi0(ppc));
               pTrack.mom[1] = T.mom(ppc)*sin(T.the0(ppc))*sin(T.phi0(ppc));
               pTrack.mom[2] = T.mom(ppc)*cos(T.the0(ppc));
               
               nTrack.mom[0]   = T.mom(npc)*sin(T.the0(npc))*cos(T.phi0(npc));
               nTrack.mom[1]   = T.mom(npc)*sin(T.the0(npc))*sin(T.phi0(npc));
               nTrack.mom[2]   = T.mom(npc)*cos(T.the0(npc));
               
               if (IsGhostCut(T.zed(ppc) - T.zed(npc), 
                     T.alpha(ppc) - T.alpha(npc), 
                     T.phi(ppc) - T.phi(npc))) continue;
                  
               if (IsOneArmCut(T.phi(ppc), T.phi(npc))) continue;
               
               const double mass = GetMass(pTrack.mom, nTrack.mom, pTrack.mass, nTrack.mass);
               const double pT = GetPairPT(pTrack.mom, nTrack.mom)
               
               id1.pc2 = pTrack.idPC2[i];
               id1.pc3 = pTrack.idPC3[i];
               id1.emc = pTrack.idEMC[i];
               id1.tof = pTrack.idTOF[i];

               id2.pc2 = nTrack.idPC2[j];
               id2.pc3 = nTrack.idPC3[j];
               id2.emc = nTrack.idEMC[j];
               id2.tof = nTrack.idTOF[j];

               if (IsnoPID(&id1, &id2))
               {
                  for (int k = 0; k < CentralityContainer.size; k++)
                  {
                     //probabilities of not registering dc-pc1 tracks
                     double probNoTrack1 = Product(
                        1. - pTrack.weightPC2[k][i], 
                        1. - pTrack.weightPC3[k][i], 
                        1. - pTrack.weightTOF[k][i], 
                        1. - pTrack.weightEMC[k][i]);
                     
                     double probNoTrack2 = Product(
                        1. - nTrack.weightPC2[k][j], 
                        1. - nTrack.weightPC3[k][j], 
                        1. - nTrack.weightTOF[k][j], 
                        1. - nTrack.weightEMC[k][j]);

                     double pairWeight = pTrack.embDCPC1.get()[k]*nTrack.embDCPC1.get()[k]*
                        (1. - probNoTrack1)*(1. - probNoTrack2);
      
                     distrInvMNoPID[k]->Fill(pT, mass, eventWeight*pairWeight);

                     distrOrigPTVsRecPT->Fill(pT, origPt, eventWeight*pairWeight);

                     probNoTrack1 = Product(
                        1. - pTrackDecreasedAcceptance.weightPC2[k][i], 
                        1. - pTrackDecreasedAcceptance.weightPC3[k][i], 
                        1. - pTrackDecreasedAcceptance.weightTOF[k][i], 
                        1. - pTrackDecreasedAcceptance.weightEMC[k][i]);
                     
                     probNoTrack2 = Product(
                        1. - nTrackDecreasedAcceptance.weightPC2[k][j], 
                        1. - nTrackDecreasedAcceptance.weightPC3[k][j], 
                        1. - nTrackDecreasedAcceptance.weightTOF[k][j], 
                        1. - nTrackDecreasedAcceptance.weightEMC[k][j]);

                     pairWeight = pTrack.embDCPC1.get()[k]*nTrack.embDCPC1.get()[k]*
                        (1. - probNoTrack1)*(1. - probNoTrack2);

                     distrInvMNoPIDDecreasedAcceptance[k]->Fill(pT, mass, eventWeight*pairWeight);

                     probNoTrack1 = Product(
                        1. - pTrackIncreasedAcceptance.weightPC2[k][i], 
                        1. - pTrackIncreasedAcceptance.weightPC3[k][i], 
                        1. - pTrackIncreasedAcceptance.weightTOF[k][i], 
                        1. - pTrackIncreasedAcceptance.weightEMC[k][i]);
                     
                     probNoTrack2 = Product(
                        1. - nTrackIncreasedAcceptance.weightPC2[k][j], 
                        1. - nTrackIncreasedAcceptance.weightPC3[k][j], 
                        1. - nTrackIncreasedAcceptance.weightTOF[k][j], 
                        1. - nTrackIncreasedAcceptance.weightEMC[k][j]);
                     
                     pairWeight = pTrack.embDCPC1.get()[k]*nTrack.embDCPC1.get()[k]*
                        (1. - probNoTrack1)*(1. - probNoTrack2);

                     distrInvMNoPIDIncreasedAcceptance[k]->Fill(pT, mass, eventWeight*pairWeight);
                  }
               }
               else continue;

               if (Is1PID(&id1, &id2))
               {
                  for (int k = 0; k < CentralityContainer.size; k++)
                  {
                     double probNoTrack1 = Product(
                        1. - pTrack.weightPC2[k][i], 
                        1. - pTrack.weightPC3[k][i], 
                        1. - pTrack.weightTOF[k][i], 
                        1. - pTrack.weightEMC[k][i]);

                     double probNoIdTrack1 = 1. - pTrack.weightIdTOF[k][i];
                     
                     double probNoTrack2 = Product(
                        1. - nTrack.weightPC2[k][j], 
                        1. - nTrack.weightPC3[k][j], 
                        1. - nTrack.weightTOF[k][j], 
                        1. - nTrack.weightEMC[k][j]);

                     double probNoIdTrack2 = 1. - nTrack.weightIdTOF[k][j];
                     
                     double pairWeight = pTrack.embDCPC1.get()[k]*nTrack.embDCPC1.get()[k]*(1. -
                        (1. - (1. - probNoIdTrack1)*(1. - probNoTrack2))*
                        (1. - (1. - probNoTrack1)*(1. - probNoIdTrack2)));

                     distrInvM1PID[k]->Fill(pT, mass, eventWeight*pairWeight);

                     probNoTrack1 = Product(
                        1. - pTrackDecreasedAcceptance.weightPC2[k][i], 
                        1. - pTrackDecreasedAcceptance.weightPC3[k][i], 
                        1. - pTrackDecreasedAcceptance.weightTOF[k][i], 
                        1. - pTrackDecreasedAcceptance.weightEMC[k][i]);

                     probNoIdTrack1 = 1. - pTrackDecreasedAcceptance.weightIdTOF[k][i];
                     
                     probNoTrack2 = Product(
                        1. - nTrackDecreasedAcceptance.weightPC2[k][j], 
                        1. - nTrackDecreasedAcceptance.weightPC3[k][j], 
                        1. - nTrackDecreasedAcceptance.weightTOF[k][j], 
                        1. - nTrackDecreasedAcceptance.weightEMC[k][j]);

                     probNoIdTrack2 = 1. - nTrackDecreasedAcceptance.weightIdTOF[k][j];
                     
                     pairWeight = pTrack.embDCPC1.get()[k]*nTrack.embDCPC1.get()[k]*(1. -
                        (1. - (1. - probNoIdTrack1)*(1. - probNoTrack2))*
                        (1. - (1. - probNoTrack1)*(1. - probNoIdTrack2)));

                     distrInvM1PIDDecreasedAcceptance[k]->Fill(pT, mass, eventWeight*pairWeight);

                     probNoTrack1 = Product(
                        1. - pTrackIncreasedAcceptance.weightPC2[k][i], 
                        1. - pTrackIncreasedAcceptance.weightPC3[k][i], 
                        1. - pTrackIncreasedAcceptance.weightTOF[k][i], 
                        1. - pTrackIncreasedAcceptance.weightEMC[k][i]);

                     probNoIdTrack1 = 1. - pTrackIncreasedAcceptance.weightIdTOF[k][i];
                     
                     probNoTrack2 = Product(
                        1. - nTrackIncreasedAcceptance.weightPC2[k][j], 
                        1. - nTrackIncreasedAcceptance.weightPC3[k][j], 
                        1. - nTrackIncreasedAcceptance.weightTOF[k][j], 
                        1. - nTrackIncreasedAcceptance.weightEMC[k][j]);
                     
                     probNoIdTrack2 = 1. - nTrackIncreasedAcceptance.weightIdTOF[k][j];
                     
                     pairWeight = pTrack.embDCPC1.get()[k]*nTrack.embDCPC1.get()[k]*(1. -
                        (1. - (1. - probNoIdTrack1)*(1. - probNoTrack2))*
                        (1. - (1. - probNoTrack1)*(1. - probNoIdTrack2)));
                     
                     distrInvM1PIDIncreasedAcceptance[k]->Fill(pT, mass, eventWeight*pairWeight);
                  }
               }

               if (Is2PID(&id1, &id2))
               {
                  for (int k = 0; k < CentralityContainer.size; k++)
                  {
                     double probNoTrack1 = Product(
                        1. - pTrack.weightIdTOF[k][i], 
                        1. - pTrack.weightIdEMC[k][i]);
                     
                     double probNoTrack2 = Product(
                        1. - nTrack.weightIdTOF[k][j], 
                        1. - nTrack.weightIdEMC[k][j]);
                     
                     double pairWeight = pTrack.embDCPC1.get()[k]*nTrack.embDCPC1.get()[k]*
                        (1. - probNoTrack1)*(1. - probNoTrack2);

                     distrInvM2PID[k]->Fill(pT, mass, eventWeight*pairWeight);

                     probNoTrack1 = Product(
                        1. - pTrackDecreasedAcceptance.weightIdTOF[k][i], 
                        1. - pTrackDecreasedAcceptance.weightIdEMC[k][i]);
                     
                     probNoTrack2 = Product(
                        1. - nTrackDecreasedAcceptance.weightIdTOF[k][j], 
                        1. - nTrackDecreasedAcceptance.weightIdEMC[k][j]);
                     
                     pairWeight = pTrack.embDCPC1.get()[k]*nTrack.embDCPC1.get()[k]*
                        (1. - probNoTrack1)*(1. - probNoTrack2);

                     distrInvM2PIDDecreasedAcceptance[k]->Fill(pT, mass, eventWeight*pairWeight);

                     probNoTrack1 = Product(
                        1. - pTrackIncreasedAcceptance.weightIdTOF[k][i], 
                        1. - pTrackIncreasedAcceptance.weightIdEMC[k][i]);
                     
                     probNoTrack2 = Product(
                        1. - nTrackIncreasedAcceptance.weightIdTOF[k][j], 
                        1. - nTrackIncreasedAcceptance.weightIdEMC[k][j]);
                     
                     pairWeight = pTrack.embDCPC1.get()[k]*nTrack.embDCPC1.get()[k]*
                        (1. - probNoTrack1)*(1. - probNoTrack2);

                     distrInvM2PIDIncreasedAcceptance[k]->Fill(pT, mass, eventWeight*pairWeight);

                     probNoTrack1 = Product(
                        1. - pTrack.weightIdTOF[k][i], 
                        1. - pTrackDecreasedM2Eff.weightIdEMC[k][i]);
                     
                     probNoTrack2 = Product(
                        1. - nTrack.weightIdTOF[k][j], 
                        1. - nTrackDecreasedM2Eff.weightIdEMC[k][j]);
                     
                     pairWeight = pTrack.embDCPC1.get()[k]*nTrack.embDCPC1.get()[k]*
                        (1. - probNoTrack1)*(1. - probNoTrack2);

                     distrInvM2PIDDecreasedM2Eff[k]->Fill(pT, mass, eventWeight*pairWeight);

                     probNoTrack1 = Product(
                        1. - pTrack.weightIdTOF[k][i], 
                        1. - pTrackIncreasedM2Eff.weightIdEMC[k][i]);
                     
                     probNoTrack2 = Product(
                        1. - nTrack.weightIdTOF[k][j], 
                        1. - nTrackIncreasedM2Eff.weightIdEMC[k][j]);
                     
                     pairWeight = pTrack.embDCPC1.get()[k]*nTrack.embDCPC1.get()[k]*
                        (1. - probNoTrack1)*(1. - probNoTrack2);

                     distrInvM2PIDIncreasedM2Eff[k]->Fill(pT, mass, eventWeight*pairWeight);
                  }
               }
               else continue;

               if (IsTOF2PID(&id1, &id2))
               {
                  for (int k = 0; k < CentralityContainer.size; k++)
                  {
                     double pairWeight = 
                        pTrack.embDCPC1.get()[k]*nTrack.embDCPC1.get()[k]*
                        pTrack.weightIdTOF[k][i]*nTrack.weightIdTOF[k][j];
                     
                     distrInvMTOF2PID[k]->Fill(pT, mass, eventWeight*pairWeight);
                     
                     pairWeight = 
                        pTrack.embDCPC1.get()[k]*nTrack.embDCPC1.get()[k]*
                        pTrackDecreasedAcceptance.weightIdTOF[k][i]*nTrackDecreasedAcceptance.weightIdTOF[k][j];
                     
                     distrInvMTOF2PIDDecreasedAcceptance[k]->Fill(pT, mass, eventWeight*pairWeight);

                     pairWeight = 
                        pTrack.embDCPC1.get()[k]*nTrack.embDCPC1.get()[k]*
                        pTrackIncreasedAcceptance.weightIdTOF[k][i]*nTrackIncreasedAcceptance.weightIdTOF[k][j];
                     
                     distrInvMTOF2PIDIncreasedAcceptance[k]->Fill(pT, mass, eventWeight*pairWeight);
                  }
               }

               if (IsEMC2PID(&id1, &id2))
               {
                  for (int k = 0; k < CentralityContainer.size; k++)
                  {
                     double pairWeight = 
                        pTrack.embDCPC1.get()[k]*nTrack.embDCPC1.get()[k]*
                        pTrack.weightIdEMC[k][i]*nTrack.weightIdEMC[k][j];

                     distrInvMEMC2PID[k]->Fill(pT, mass, eventWeight*pairWeight);

                     pairWeight = 
                        pTrack.embDCPC1.get()[k]*nTrack.embDCPC1.get()[k]*
                        pTrackDecreasedAcceptance.weightIdEMC[k][i]*nTrackDecreasedAcceptance.weightIdEMC[k][j];

                     distrInvMEMC2PIDDecreasedAcceptance[k]->Fill(pT, mass, eventWeight*pairWeight);

                     pairWeight = 
                        pTrack.embDCPC1.get()[k]*nTrack.embDCPC1.get()[k]*
                        pTrackIncreasedAcceptance.weightIdEMC[k][i]*nTrackIncreasedAcceptance.weightIdEMC[k][j];

                     distrInvMEMC2PIDIncreasedAcceptance[k]->Fill(pT, mass, eventWeight*pairWeight);
                     
                     pairWeight = 
                        pTrack.embDCPC1.get()[k]*nTrack.embDCPC1.get()[k]*
                        pTrackDecreasedM2Eff.weightIdEMC[k][i]*nTrackDecreasedM2Eff.weightIdEMC[k][j];

                     distrInvMEMC2PIDDecreasedM2Eff[k]->Fill(pT, mass, eventWeight*pairWeight);

                     pairWeight = 
                        pTrack.embDCPC1.get()[k]*nTrack.embDCPC1.get()[k]*
                        pTrackIncreasedM2Eff.weightIdEMC[k][i]*nTrackIncreasedM2Eff.weightIdEMC[k][j];

                     distrInvMEMC2PIDIncreasedM2Eff[k]->Fill(pT, mass, eventWeight*pairWeight);
                  }
               }
            }//cylce by negative tracks
         }//cycle by positive traks 
      }//cycle by events
   };

   auto PbarCall = [&]()
   {
      ProgressBar pBar = ProgressBar("Fancy");
      while (!isProcessFinished)
      {
         pBar.Print(ncalls/nevents);
         std::this_thread::sleep_for(std::chrono::milliseconds(20));
      }
      pBar.Print(1.);
   };
   
   std::thread pBarThread(PbarCall);
   tp.Process(ProcessMP);
   isProcessFinished = true;
   pBarThread.join();
}

void AnalyzeTrackPair()
{
   if (Par.doUseWeightFunc) CheckInputFile("../input/Spectra/" + 
      Par.system + "/" + Par.origParticle.name + ".txt");
   
   for (long unsigned int i = 0; i < Par.daughter1Queue.size(); i++)
   {
      const int indexDaughter2 = ParticleProperties.iterMap[Par.daughter1Queue[i]];
      const int indexDaughter2 = ParticleProperties.iterMap[Par.daughter2Queue[i]];   
      
      const std::string decayChannel = 
         ParticleProperties.shortName[indexDaughter2] + 
         ParticleProperties.shortName[indexDaughter2];
      
      for (std::string detector : Par.detectors)
      {
         CheckInputFile("../input/Embedding/" + Par.runName + "/" + 
            ParticleProperties.absName[indexDaughter2] + "_" + detector + ".txt");
         CheckInputFile("../input/Embedding/" + Par.runName + "/" + 
            ParticleProperties.absName[indexDaughter2] + "_" + detector + ".txt");   
      }
      
      for (std::string magf : Par.magfQueue)
      {
         CheckInputFile("../input/M2Eff/" + Par.runName + "/" +
               ParticleProperties.name[indexDaughter2] + "_EMCale" + magf + ".txt");
         CheckInputFile("../input/M2Eff/" + Par.runName + "/" +
               ParticleProperties.name[indexDaughter2] + "_EMCalw" + magf + ".txt");
         CheckInputFile("../input/Systematics/" + Par.runName + "/m2eff_" +
               ParticleProperties.name[indexDaughter2] + "_EMCale" + magf + ".txt");
         CheckInputFile("../input/Systematics/" + Par.runName + "/m2eff_" +
               ParticleProperties.name[indexDaughter2] + "_EMCalw" + magf + ".txt");

         CheckInputFile("../input/M2Eff/" + Par.runName + "/" +
               ParticleProperties.name[indexDaughter2] + "_EMCale" + magf + ".txt");
         CheckInputFile("../input/M2Eff/" + Par.runName + "/" +
               ParticleProperties.name[indexDaughter2] + "_EMCalw" + magf + ".txt");
         CheckInputFile("../input/Systematics/" + Par.runName + "/m2eff_" +
               ParticleProperties.name[indexDaughter2] + "_EMCale" + magf + ".txt");
         CheckInputFile("../input/Systematics/" + Par.runName + "/m2eff_" +
               ParticleProperties.name[indexDaughter2] + "_EMCalw" + magf + ".txt");
         for (std::string auxName : Par.auxNameQueue)
         {
            const std::string inputFileName = "../data/" + Par.runName + "/Resonances/" + 
               Par.origParticle.name + "_" + decayChannel + magf + auxName + ".root";
            
            CheckInputFile(inputFileName);
         }
      }
   }

   CheckInputFile("../input/Systematics/" + Par.runName + "/acceptance.txt");
   
   int procNum = 1;
   for (double pTDeviation : Par.pTDeviationQueue)
   {
      thrContainer thrContainer;

      const double invMassMin = 0.;

      for (long unsigned int j = 0; j < CentralityContainer.nameWithoutPercent.size(); j++)
      {
         std::string dir = CentralityContainer.nameWithoutPercent[j];
         
         thrContainer.distrInvMNoPID[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "noPID_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "noPID", Par.pTNBins, Par.pTMinPair, Par.pTMaxPair, 
               Par.invMNBins, invMassMin, 4., dir));
         
         thrContainer.distrInvMNoPIDDecreasedAcceptance[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "noPID_A-_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
                "noPID", Par.pTNBins, Par.pTMinPair, Par.pTMaxPair, 
                Par.invMNBins, invMassMin, 4., dir));
         
         thrContainer.distrInvMNoPIDIncreasedAcceptance[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "noPID_A+_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "noPID", Par.pTNBins, Par.pTMinPair, Par.pTMaxPair, 
               Par.invMNBins, invMassMin, 4., dir));

         thrContainer.distrInvM1PID[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "1PID_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "1PID", Par.pTNBins, Par.pTMinPair, Par.pTMaxPair, 
               Par.invMNBins, invMassMin, 4., dir));
         
         thrContainer.distrInvM1PIDDecreasedAcceptance[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "1PID_A-_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "1PID", Par.pTNBins, Par.pTMinPair, Par.pTMaxPair, 
               Par.invMNBins, invMassMin, 4., dir));
         
         thrContainer.distrInvM1PIDIncreasedAcceptance[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "1PID_A+_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "1PID", Par.pTNBins, Par.pTMinPair, Par.pTMaxPair, 
               Par.invMNBins, invMassMin, 4., dir));

         thrContainer.distrInvM2PID[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "2PID_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "2PID", Par.pTNBins, Par.pTMinPair, Par.pTMaxPair, 
               Par.invMNBins, invMassMin, 4., dir));
         
         thrContainer.distrInvM2PIDDecreasedAcceptance[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "2PID_A-_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "2PID", Par.pTNBins, Par.pTMinPair, Par.pTMaxPair, 
               Par.invMNBins, invMassMin, 4., dir));
         
         thrContainer.distrInvM2PIDIncreasedAcceptance[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "2PID_A+_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "2PID", Par.pTNBins, Par.pTMinPair, Par.pTMaxPair, 
               Par.invMNBins, invMassMin, 4., dir));
         
         thrContainer.distrInvM2PIDDecreasedM2Eff[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "2PID_m2Eff-_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "2PID", Par.pTNBins, Par.pTMinPair, Par.pTMaxPair, 
               Par.invMNBins, invMassMin, 4., dir));
         
         thrContainer.distrInvM2PIDIncreasedM2Eff[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "2PID_m2Eff+_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "2PID", Par.pTNBins, Par.pTMinPair, Par.pTMaxPair, 
               Par.invMNBins, invMassMin, 4., dir));
         
         thrContainer.distrInvMTOF2PID[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "TOF2PID_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "TOF2PID", Par.pTNBins, Par.pTMinPair, Par.pTMaxPair, 
               Par.invMNBins, invMassMin, 4., dir));
         
         thrContainer.distrInvMTOF2PIDDecreasedAcceptance[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "TOF2PID_A-_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "TOF2PID", Par.pTNBins, Par.pTMinPair, Par.pTMaxPair, 
               Par.invMNBins, invMassMin, 4., dir));
         
         thrContainer.distrInvMTOF2PIDIncreasedAcceptance[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "TOF2PID_A+_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "TOF2PID", Par.pTNBins, Par.pTMinPair, Par.pTMaxPair, 
               Par.invMNBins, invMassMin, 4., dir));

         thrContainer.distrInvMEMC2PID[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "EMC2PID_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "EMC2PID", Par.pTNBins, Par.pTMinPair, Par.pTMaxPair, 
               Par.invMNBins, invMassMin, 4., dir));
         
         thrContainer.distrInvMEMC2PIDDecreasedAcceptance[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "EMC2PID_A-_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "EMC2PID", Par.pTNBins, Par.pTMinPair, Par.pTMaxPair, 
               Par.invMNBins, invMassMin, 4., dir));
         
         thrContainer.distrInvMEMC2PIDIncreasedAcceptance[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "EMC2PID_A+_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "EMC2PID", Par.pTNBins, Par.pTMinPair, Par.pTMaxPair, 
               Par.invMNBins, invMassMin, 4., dir));
         
         thrContainer.distrInvMEMC2PIDDecreasedM2Eff[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "EMC2PID_m2Eff-_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "EMC2PID", Par.pTNBins, Par.pTMinPair, Par.pTMaxPair, 
               Par.invMNBins, invMassMin, 4., dir));
         
         thrContainer.distrInvMEMC2PIDIncreasedM2Eff[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "EMC2PID_m2Eff+_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "EMC2PID", Par.pTNBins, Par.pTMinPair, Par.pTMaxPair, 
               Par.invMNBins, invMassMin, 4., dir));
      }
      
      for (long unsigned int i = 0; i < Par.daughter1Queue.size(); i++)
      {
         for (std::string magf : Par.magfQueue)
         {   
            for (std::string auxName : Par.auxNameQueue)
            {
               AnalyzeConfiguration(&thrContainer, Par.daughter1Queue[i], Par.daughter2Queue[i], 
                  magf, auxName, pTDeviation, procNum);
               procNum++;
            }
         }   
      }
      
      system(("mkdir -p ../../analysis/data/phenix_sim/" + Par.runName).c_str());

      std::string pTDeviationName;
      if (abs(pTDeviation - 1) < 1e-7) pTDeviationName = "";
      else pTDeviationName = "_pt" + DtoStr(pTDeviation, 3);
   
      const std::string outputFileName = 
         "../../analysis/data/phenix_sim/" + 
         Par.runName + "/" + Par.origParticle.name + 
         pTDeviationName + ".root";
      
      ThrObjHolder.Write(outputFileName);
      ThrObjHolder.Clear();

      PrintInfo("File " + outputFileName + " was written");
   }
}

int main()
{
   AnalyzeTrackPair();
   return 0;
}

#endif /* ANALYZE_RESONANCE_CPP */
