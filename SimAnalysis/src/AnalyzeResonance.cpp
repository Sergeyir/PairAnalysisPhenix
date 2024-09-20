// $SOURCE$
//------------------------------------------------------------------------------------------------
//                     AnalyzeResonance function realisation
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
                          const std::string& auxName, const double ptDeviation, 
                          const int procNum)
{
   Box box = Box("Parameters of run " + std::to_string(procNum) + " out of " + std::to_string(
      Par.daughter1Queue.size()*Par.magfQueue.size()*
      Par.ptDeviationQueue.size()*Par.auxNameQueue.size()));

   
   const std::string decay_channel = 
      ParticleProperties.sname[ParticleProperties.iter_map[daughter1]] + 
      ParticleProperties.sname[ParticleProperties.iter_map[daughter2]];

   const std::string input_file_name = "../data/" + Par.runName + "/Resonances/" + 
      Par.origParticle.name + "_" + decay_channel + magf + auxName + ".root";
   
   TFile input_file = TFile(input_file_name.c_str());

   const double nevents= static_cast<double>(((TTree *) input_file.Get("Tree"))->GetEntries());
   
   if (nevents <= 0)
   {
      Print("Error: Number of events is less or equal to 0");
      exit(1);
   }
   
   TH1F *origPtDistribution = (TH1F *) input_file.Get("orig_pt");

   const double orig_pt_threshold = origPtDistribution->Integral()/
      static_cast<double>(origPtDistribution->GetXaxis()->GetNbins())/2.;

   const double low_pt_bound = origPtDistribution->GetXaxis()->GetBinLowEdge(
      origPtDistribution->FindFirstBinAbove(orig_pt_threshold));
   const double up_pt_bound = origPtDistribution->GetXaxis()->GetBinUpEdge(
      origPtDistribution->FindLastBinAbove(orig_pt_threshold));

   double event_norm = 1.;
   if (Par.doUseWeightFunc)
   {
      event_norm = origPtDistribution->Integral(
         origPtDistribution->GetXaxis()->FindBin(Par.ptMinPair),
         origPtDistribution->GetXaxis()->FindBin(Par.ptMaxPair));

      //this normalization is needed to merge 2 files with flapt pt distributin with different ranges
      event_norm *= (Par.ptMaxPair - Par.ptMinPair)/(up_pt_bound - low_pt_bound);
   }

   box.AddEntry("Run name", Par.runName);
   box.AddEntry("Orig particle", Par.origParticle.name);
   box.AddEntry("Magnetic field configuration", magf);
   box.AddEntry("Use weight function", Par.doUseWeightFunc);
   
   box.AddEntry("Daughter particle 1", daughter1);
   box.AddEntry("Daughter particle 2", daughter2);
   
   box.AddEntry("Minimum p_T, GeV", Par.ptMin);
   box.AddEntry("Maximum p_T, GeV", Par.ptMax);

   box.AddEntry("Input file orig p_T range, GeV", DtoStr(low_pt_bound, 1) + "-" + DtoStr(up_pt_bound, 1));
   
   box.AddEntry("Number of events to be analyzed, 1e6", nevents/1e6, 3);

   box.AddEntry("pT deviation", DtoStr(ptDeviation, 3));
   
   box.AddEntry("Number of threads", Par.nThreads);
   
   box.Print();

   TF1 weight_func = TF1("weight_func", 
      "[0]*x^([5]-1)*(1 + ([1] - 1)*(sqrt([4] + x^2)^([5]) - [2])/[3])^(-1./([1]-1.))");
   weight_func.SetParameters(
      ReadFileIntoArray("../input/Spectra/" + Par.system + "/" + Par.origParticle.name + ".txt", 6));
   
   ROOT::EnableImplicitMT(Par.nThreads);
   ROOT::TTreeProcessorMT tp(input_file_name.c_str());
   
   double ncalls = 0.;

   //acceptance systematic uncertainties of detectors
   double dceast0_sys, dceast1_sys, dcwest0_sys, dcwest1_sys,
      pc1e_sys, pc1w_sys, pc2_sys, pc3e_sys, pc3w_sys,
      tofe0_sys, tofe1_sys, tofw0_sys, tofw1_sys;
   std::array<double, 4> emcale_pos_sys, emcale_neg_sys, emcalw_pos_sys, emcalw_neg_sys;

   ReadFile("../input/Systematics/" + Par.runName + "/acceptance.txt",
      dceast0_sys, dceast1_sys, dcwest0_sys, dcwest1_sys,
      pc1e_sys, pc1w_sys, pc2_sys, pc3e_sys, pc3w_sys,
      tofe0_sys, tofe1_sys, tofw0_sys, tofw1_sys,
      emcale_pos_sys[0], emcale_neg_sys[0], emcale_pos_sys[1], emcale_neg_sys[1],
      emcale_pos_sys[2], emcale_neg_sys[2], emcale_pos_sys[3], emcale_neg_sys[3],
      emcalw_pos_sys[0], emcalw_neg_sys[0], emcalw_pos_sys[1], emcalw_neg_sys[1],
      emcalw_pos_sys[2], emcalw_neg_sys[2], emcalw_pos_sys[3], emcalw_neg_sys[3]);
   
   bool process_finished = false;
   auto ProcessMP = [&](TTreeReader &reader)
   {
      //arrays and constants for positives and negative particles and their tracks
      ParticleContainer ptrack, ntrack;

      //arrays for keeping weights for varied acceptancies
      ParticleContainer ptrackDecreasedAcceptance, ptrackIncreasedAcceptance, 
                        ntrackDecreasedAcceptance, ntrackIncreasedAcceptance;
      ParticleContainer ptrackDecreasedM2Eff, ptrackIncreasedM2Eff, 
                        ntrackDecreasedM2Eff, ntrackIncreasedM2Eff;

      //arrays for keeping weights for varied /epsilon_{m^2}
      
      //setting daughter particles
      ptrack.iter = ParticleProperties.iter_map[daughter1];
      ntrack.iter = ParticleProperties.iter_map[daughter2];
      
      ptrack.origId = abs(ParticleProperties.id[ptrack.iter]);
      ntrack.origId = abs(ParticleProperties.id[ntrack.iter]);

      idContainer Id1, Id2;

      Id1.origId = ptrack.origId;
      Id2.origId = ntrack.origId;

      ptrack.geantId = ParticleProperties.geantId[ptrack.iter];
      ntrack.geantId = ParticleProperties.geantId[ntrack.iter];

      ptrack.mass = ParticleProperties.mass[ptrack.iter];
      ntrack.mass = ParticleProperties.mass[ntrack.iter];

      //setting embedding
      ptrack.embDCPC1 = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
         ParticleProperties.abs_name[ptrack.iter] + "_dc_pc1.txt", 
         CentralityContainer.size));
      ptrack.embPC2 = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
         ParticleProperties.abs_name[ptrack.iter] + "_pc2.txt", 
         CentralityContainer.size));
      ptrack.embPC3 = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
         ParticleProperties.abs_name[ptrack.iter] + "_pc3.txt", 
         CentralityContainer.size));
      ptrack.embTOFe = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
         ParticleProperties.abs_name[ptrack.iter] + "_tofe.txt", 
         CentralityContainer.size));
      ptrack.embTOFw = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
         ParticleProperties.abs_name[ptrack.iter] + "_tofw.txt", 
         CentralityContainer.size));

      ntrack.embDCPC1 = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
         ParticleProperties.abs_name[ntrack.iter] + "_dc_pc1.txt", 
         CentralityContainer.size));
      ntrack.embPC2 = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
         ParticleProperties.abs_name[ntrack.iter] + "_pc2.txt", 
         CentralityContainer.size));
      ntrack.embPC3 = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
         ParticleProperties.abs_name[ntrack.iter] + "_pc3.txt", 
         CentralityContainer.size));
      ntrack.embTOFe = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
         ParticleProperties.abs_name[ntrack.iter] + "_tofe.txt", 
         CentralityContainer.size));
      ntrack.embTOFw = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
         ParticleProperties.abs_name[ntrack.iter] + "_tofw.txt", 
         CentralityContainer.size));
      
      for (int i = 0; i < 4; i++)
      {
         ptrack.embEMCale[i] = std::unique_ptr<double>(
            ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
            ParticleProperties.abs_name[ptrack.iter] + "_emcale" + std::to_string(i) + ".txt", 
            CentralityContainer.size));
         ptrack.embEMCalw[i] = std::unique_ptr<double>(
            ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
            ParticleProperties.abs_name[ptrack.iter] + "_emcalw" + std::to_string(i) + ".txt", 
            CentralityContainer.size));
         
         ntrack.embEMCale[i] = std::unique_ptr<double>(
            ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
            ParticleProperties.abs_name[ntrack.iter] + "_emcale" + std::to_string(i) + ".txt", 
            CentralityContainer.size));
         ntrack.embEMCalw[i] = std::unique_ptr<double>(
            ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
            ParticleProperties.abs_name[ntrack.iter] + "_emcalw" + std::to_string(i) + ".txt", 
            CentralityContainer.size));
      }

      ptrack.m2EffEMCale = std::unique_ptr<double>(
         ReadFileIntoArray("../input/M2Eff/" + Par.runName + "/" +
         ParticleProperties.name[ptrack.iter] + "_EMCale" + magf + ".txt", 2));
      ptrack.m2EffEMCalw = std::unique_ptr<double>(
         ReadFileIntoArray("../input/M2Eff/" + Par.runName + "/" +
         ParticleProperties.name[ptrack.iter] + "_EMCalw" + magf + ".txt", 4));
      
      ntrack.m2EffEMCale = std::unique_ptr<double>(
         ReadFileIntoArray("../input/M2Eff/" + Par.runName + "/" +
         ParticleProperties.name[ntrack.iter] + "_EMCale" + magf + ".txt", 2));
      ntrack.m2EffEMCalw = std::unique_ptr<double>(
         ReadFileIntoArray("../input/M2Eff/" + Par.runName + "/" +
         ParticleProperties.name[ntrack.iter] + "_EMCalw" + magf + ".txt", 4));

      ptrack.m2EffSysEMCale = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Systematics/" + Par.runName + "/m2eff_" +
         ParticleProperties.name[ptrack.iter] + "_EMCale" + magf + ".txt", 2));
      ptrack.m2EffSysEMCalw = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Systematics/" + Par.runName + "/m2eff_" +
         ParticleProperties.name[ptrack.iter] + "_EMCalw" + magf + ".txt", 4));
      
      ntrack.m2EffSysEMCale = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Systematics/" + Par.runName + "/m2eff_" +
         ParticleProperties.name[ntrack.iter] + "_EMCale" + magf + ".txt", 2));
      ntrack.m2EffSysEMCalw = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Systematics/" + Par.runName + "/m2eff_" +
         ParticleProperties.name[ntrack.iter] + "_EMCalw" + magf + ".txt", 4));

      std::shared_ptr<TH1F> origPtDistr = thrContainer->origPtDistr.Get();
      
      std::shared_ptr<TH2F> origPtVsPtDistr = thrContainer->origPtVsPtDistr.Get();
      
      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> invMDistrNoPID;
      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> invMDistrNoPIDDecreasedAcceptance;
      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> invMDistrNoPIDIncreasedAcceptance;

      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> invMDistr1PID;
      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> invMDistr1PIDDecreasedAcceptance;
      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> invMDistr1PIDIncreasedAcceptance;

      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> invMDistr2PID;
      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> invMDistr2PIDDecreasedAcceptance;
      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> invMDistr2PIDIncreasedAcceptance;
      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> invMDistr2PIDDecreasedM2Eff;
      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> invMDistr2PIDIncreasedM2Eff;

      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> invMDistrTOF2PID;
      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> invMDistrTOF2PIDDecreasedAcceptance;
      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> invMDistrTOF2PIDIncreasedAcceptance;

      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> invMDistrEMC2PID;
      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> invMDistrEMC2PIDDecreasedAcceptance;
      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> invMDistrEMC2PIDIncreasedAcceptance;
      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> invMDistrEMC2PIDDecreasedM2Eff;
      std::array<std::shared_ptr<TH2F>, CentralityContainer.size> invMDistrEMC2PIDIncreasedM2Eff;

      for (int i = 0; i < CentralityContainer.size; i++)
      {
         invMDistrNoPID[i] = thrContainer->invMDistrNoPID[i]->Get();
         invMDistrNoPIDDecreasedAcceptance[i] = thrContainer->invMDistrNoPIDDecreasedAcceptance[i]->Get();
         invMDistrNoPIDIncreasedAcceptance[i] = thrContainer->invMDistrNoPIDIncreasedAcceptance[i]->Get();

         invMDistr1PID[i] = thrContainer->invMDistr1PID[i]->Get();
         invMDistr1PIDDecreasedAcceptance[i] = thrContainer->invMDistr1PIDDecreasedAcceptance[i]->Get();
         invMDistr1PIDIncreasedAcceptance[i] = thrContainer->invMDistr1PIDIncreasedAcceptance[i]->Get();

         invMDistr2PID[i] = thrContainer->invMDistr2PID[i]->Get();
         invMDistr2PIDDecreasedAcceptance[i] = thrContainer->invMDistr2PIDDecreasedAcceptance[i]->Get();
         invMDistr2PIDIncreasedAcceptance[i] = thrContainer->invMDistr2PIDIncreasedAcceptance[i]->Get();
         invMDistr2PIDDecreasedM2Eff[i] = thrContainer->invMDistr2PIDDecreasedM2Eff[i]->Get();
         invMDistr2PIDIncreasedM2Eff[i] = thrContainer->invMDistr2PIDIncreasedM2Eff[i]->Get();

         invMDistrTOF2PID[i] = thrContainer->invMDistrTOF2PID[i]->Get();
         invMDistrTOF2PIDDecreasedAcceptance[i] = thrContainer->invMDistrTOF2PIDDecreasedAcceptance[i]->Get();
         invMDistrTOF2PIDIncreasedAcceptance[i] = thrContainer->invMDistrTOF2PIDIncreasedAcceptance[i]->Get();

         invMDistrEMC2PID[i] = thrContainer->invMDistrEMC2PID[i]->Get();
         invMDistrEMC2PIDDecreasedAcceptance[i] = thrContainer->invMDistrEMC2PIDDecreasedAcceptance[i]->Get();
         invMDistrEMC2PIDIncreasedAcceptance[i] = thrContainer->invMDistrEMC2PIDIncreasedAcceptance[i]->Get();
         invMDistrEMC2PIDDecreasedM2Eff[i] = thrContainer->invMDistrEMC2PIDDecreasedM2Eff[i]->Get();
         invMDistrEMC2PIDIncreasedM2Eff[i] = thrContainer->invMDistrEMC2PIDIncreasedM2Eff[i]->Get();
      }
      
      EffTreeReader T(reader);
      
      while (reader.Next())
      {   
         ncalls += 1.;
         const double origPt = sqrt(pow(T.mom_orig(0), 2) + pow(T.mom_orig(1), 2))*ptDeviation;
         
         double event_weight = 1.;
         
         if (Par.doUseWeightFunc) event_weight = weight_func.Eval(orig_pt)/event_norm*4e11;
         
         origPtDistr->Fill(origPt, event_weight);
         
         if (T.nch() > 49 || T.nch() <= 0) continue;
         
         const double bbcz = T.bbcz();
         if (fabs(bbcz) > 30) continue;

         ptrack.size = 0;
         ntrack.size = 0;

         for(int i = 0; i < T.nch(); i++)
         {   
            const double the0 = T.the0(i);
            const double pt = (T.mom(i))*sin(the0)*ptDeviation;

            if (pt < Par.ptMin || pt > Par.ptMax) continue;
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

            ptrack.ResetTrack(ptrack.size);
            ntrack.ResetTrack(ntrack.size);
            ptrackDecreasedAcceptance.ResetTrack(ptrack.size);
            ntrackDecreasedAcceptance.ResetTrack(ntrack.size);
            ptrackIncreasedAcceptance.ResetTrack(ptrack.size);
            ntrackIncreasedAcceptance.ResetTrack(ntrack.size);
            ptrackDecreasedM2Eff.ResetTrack(ptrack.size);
            ntrackDecreasedM2Eff.ResetTrack(ntrack.size);
            ptrackIncreasedM2Eff.ResetTrack(ptrack.size);
            ntrackIncreasedM2Eff.ResetTrack(ntrack.size);

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
                  
                  const double id = GetTOFePID(pt, charge, m2);
                  
                  //variation of acceptance
                  double acc_var;
                  if (zed >= 0) acc_var = ErrPropagation(tofe0_sys, dceast0_sys, pc1e_sys);
                  else acc_var = ErrPropagation(tofe1_sys, dceast1_sys, pc1e_sys);
   
                  if (charge == 1)
                  {
                     ptrack.idTOF[ptrack.size] = id;
                     
                     for (int j = 0; j < CentralityContainer.size; j++)
                     {
                        ptrack.weightTOF[j][ptrack.size] = 
                           ptrack.embTOFe.get()[j]/ptrack.embDCPC1.get()[j];
                        ptrackDecreasedAcceptance.weightTOF[j][ptrack.size] = 
                           ptrack.weightTOF[j][ptrack.size]*(1. - acc_var);
                        ptrackIncreasedAcceptance.weightTOF[j][ptrack.size] = 
                           Minimum(1., ptrack.weightTOF[j][ptrack.size]*(1.+acc_var));

                        if (id == ptrack.origId) 
                        {
                           ptrack.weightIdTOF[j][ptrack.size] = 
                              ptrack.weightTOF[j][ptrack.size];
                           ptrackDecreasedAcceptance.weightIdTOF[j][ptrack.size] = 
                              ptrack.weightIdTOF[j][ptrack.size]*(1. - acc_var);
                           ptrackIncreasedAcceptance.weightIdTOF[j][ptrack.size] = 
                              Minimum(1., ptrack.weightIdTOF[j][ptrack.size]*(1.+acc_var));
                        }   
                     }
                  }
                  else
                  {
                     ntrack.idTOF[ntrack.size] = id;
                     
                     for (int j = 0; j < CentralityContainer.size; j++)
                     {
                        ntrack.weightTOF[j][ntrack.size] = 
                           ntrack.embTOFe.get()[j]/ntrack.embDCPC1.get()[j];
                        ntrackDecreasedAcceptance.weightTOF[j][ntrack.size] = 
                           ntrack.weightTOF[j][ntrack.size]*(1. - acc_var);
                        ntrackIncreasedAcceptance.weightTOF[j][ntrack.size] = 
                           Minimum(1., ntrack.weightTOF[j][ntrack.size]*(1.+acc_var));
                        
                        if (id == ntrack.origId)
                        {
                           ntrack.weightIdTOF[j][ntrack.size] = 
                              ntrack.weightTOF[j][ntrack.size];
                           ntrackDecreasedAcceptance.weightIdTOF[j][ntrack.size] = 
                              ntrack.weightIdTOF[j][ntrack.size]*(1. - acc_var);
                           ntrackIncreasedAcceptance.weightIdTOF[j][ntrack.size] = Minimum(1.,
                              ntrack.weightIdTOF[j][ntrack.size]*(1.+acc_var));
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
                  
                  const double id = GetTOFePID(pt, charge, m2);

                  double acc_var;
                  if (zed >= 0) acc_var = ErrPropagation(tofw0_sys, dcwest0_sys, pc1w_sys);
                  else acc_var = ErrPropagation(tofw1_sys, dcwest1_sys, pc1w_sys);
                  
                  if (charge == 1)
                  {
                     ptrack.idTOF[ptrack.size] = id;
                     
                     for (unsigned int j = 0; j < CentralityContainer.size; j++)
                     {
                        ptrack.weightTOF[j][ptrack.size] = 
                           Par.correctionTOFw*ptrack.embTOFw.get()[j]/ptrack.embDCPC1.get()[j];
                        ptrackDecreasedAcceptance.weightTOF[j][ptrack.size] = 
                           ptrack.weightTOF[j][ptrack.size]*(1. - acc_var);
                        ptrackIncreasedAcceptance.weightTOF[j][ptrack.size] = Minimum(1.,
                           ptrack.weightTOF[j][ptrack.size]*(1. + acc_var));

                        if (id == ptrack.origId) 
                        {
                           ptrack.weightIdTOF[j][ptrack.size] = 
                               ptrack.weightTOF[j][ptrack.size];
                           ptrackDecreasedAcceptance.weightIdTOF[j][ptrack.size] = 
                               ptrack.weightTOF[j][ptrack.size]*(1. - acc_var);
                           ptrackIncreasedAcceptance.weightIdTOF[j][ptrack.size] = Minimum(1.,
                              ptrack.weightTOF[j][ptrack.size]*(1.+acc_var));
                        }
                     }
                  }
                  else
                  {
                     ntrack.idTOF[ntrack.size] = id;
                     
                     for (unsigned int j = 0; j < CentralityContainer.size; j++)
                     {
                        ntrack.weightTOF[j][ntrack.size] = 
                           Par.correctionTOFw*ntrack.embTOFw.get()[j]/ntrack.embDCPC1.get()[j];
                        ntrackDecreasedAcceptance.weightTOF[j][ntrack.size] = 
                           ntrack.weightTOF[j][ntrack.size]*(1. - acc_var);
                        ntrackIncreasedAcceptance.weightTOF[j][ntrack.size] = Minimum(1.,
                           ntrack.weightTOF[j][ntrack.size]*(1.+acc_var));

                        if (id == ntrack.origId) 
                        {
                           ntrack.weightIdTOF[j][ntrack.size] = 
                              ntrack.weightTOF[j][ntrack.size];
                           ntrackDecreasedAcceptance.weightIdTOF[j][ntrack.size] = 
                              ntrack.weightIdTOF[j][ntrack.size]*(1. - acc_var);
                           ntrackIncreasedAcceptance.weightIdTOF[j][ntrack.size] = Minimum(1.,
                              ntrack.weightIdTOF[j][ntrack.size]*(1.+acc_var));
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
                  double acc_var;
                  double m2_eff_var = 0.;

                  if (charge == 1) ptrack.idEMC[ptrack.size] = PartId.noPID;
                  else ntrack.idEMC[ntrack.size] = PartId.noPID;

                  if (phi > 1.5)
                  {
                     if (zed >= 0) acc_var = ErrPropagation(
                        emcale_pos_sys[T.sect(i)], dceast0_sys, pc1e_sys);
                     else acc_var = ErrPropagation(
                        emcalw_pos_sys[T.sect(i)], dceast1_sys, pc1e_sys);
                     
                     if (charge == 1)
                     {
                        for (unsigned int j = 0; j < CentralityContainer.size; j++)
                        {
                           ptrack.weightEMC[j][ptrack.size] = 
                              ptrack.embEMCale[T.sect(i)].get()[j]/ptrack.embDCPC1.get()[j];
                        }
         
                        if (T.particle_id(i) == ptrack.geantId && T.sect(i) > 1)
                        {
                           ptrack.idEMC[ptrack.size] = ptrack.origId;
                           for (unsigned int j = 0; j < CentralityContainer.size; j++)
                           {
                              ptrack.weightIdEMC[j][ptrack.size] = 
                                 ptrack.weightEMC[j][ptrack.size]*
                                 ptrack.m2EffEMCale.get()[T.sect(i)-2]*
                                 GetEMCalId(pt, ptrack.idEMC[ptrack.size], charge, phi, T.sect(i));
                              m2_eff_var = ptrack.m2EffSysEMCale.get()[T.sect(i)-2];
                           }
                        }
                     }
                     else
                     {   
                        for (unsigned int j = 0; j < CentralityContainer.size; j++)
                        {
                           ntrack.weightEMC[j][ntrack.size] = 
                              ntrack.embEMCale[T.sect(i)].get()[j]/ntrack.embDCPC1.get()[j];
                        }
                           
                        if (T.particle_id(i) == ntrack.geantId && T.sect(i) > 1)
                        {
                           ntrack.idEMC[ntrack.size] = ntrack.origId;
                           for (unsigned int j = 0; j < CentralityContainer.size; j++)
                           {
                              ntrack.weightIdEMC[j][ntrack.size] = 
                                 ntrack.weightEMC[j][ntrack.size]*
                                 ntrack.m2EffEMCale.get()[T.sect(i)-2]*
                                 GetEMCalId(pt, ntrack.idEMC[ntrack.size], charge, phi, T.sect(i));
                              m2_eff_var = ntrack.m2EffSysEMCale.get()[T.sect(i)-2];
                           }
                        }
                     }
                     if (pt > 1.1 && (T.particle_id(i) == 11 || T.particle_id(i) == 12)) 
                     {
                        switch(T.sect(i))
                        {
                           case 0: 
                              m2_eff_var *= 1.1227;
                              break;
                           case 1: 
                              m2_eff_var *= 1.0639;
                              break;
                        }
                     }
                  }
                  else 
                  {
                     if (zed >= 0) acc_var = ErrPropagation(
                        emcalw_pos_sys[T.sect(i)], dcwest0_sys, pc1w_sys);
                     else acc_var = ErrPropagation(
                        emcalw_neg_sys[T.sect(i)], dcwest1_sys, pc1w_sys);
                     
                     if (charge == 1)
                     {
                        for (unsigned int j = 0; j < CentralityContainer.size; j++)
                        {
                           ptrack.weightEMC[j][ptrack.size] = 
                              ptrack.embEMCalw[T.sect(i)].get()[j]/ptrack.embDCPC1.get()[j];
                        }
                        
                        if (T.particle_id(i) == ptrack.geantId)
                        {
                           ptrack.idEMC[ptrack.size] = ptrack.origId;
                           for (unsigned int j = 0; j < CentralityContainer.size; j++)
                           {
                              ptrack.weightIdEMC[j][ptrack.size] = 
                                 ptrack.weightEMC[j][ptrack.size]*
                                 ptrack.m2EffEMCalw.get()[T.sect(i)]*
                                 GetEMCalId(pt, ptrack.idEMC[ptrack.size], charge, phi, T.sect(i));
                              m2_eff_var = ptrack.m2EffSysEMCalw.get()[T.sect(i)];
                           }
                        }
                     }
                     else
                     {
                        for (unsigned int j = 0; j < CentralityContainer.size; j++)
                        {
                           ntrack.weightEMC[j][ntrack.size] = 
                              ntrack.embEMCalw[T.sect(i)].get()[j]/ntrack.embDCPC1.get()[j];
                        }
                        
                        if (T.particle_id(i) == ntrack.geantId)
                        {
                           ntrack.idEMC[ntrack.size] = ntrack.origId;
                           for (unsigned int j = 0; j < CentralityContainer.size; j++)
                           {
                              ntrack.weightIdEMC[j][ntrack.size] = 
                                 ntrack.weightEMC[j][ntrack.size]*
                                 ntrack.m2EffEMCalw.get()[T.sect(i)]*
                                 GetEMCalId(pt, ntrack.idEMC[ntrack.size], charge, phi, T.sect(i));
                              m2_eff_var = ntrack.m2EffSysEMCalw.get()[T.sect(i)];
                           }
                        }
                     }
                     if (pt > 1.1 && (T.particle_id(i) == 11 || T.particle_id(i) == 12)) 
                     {
                        switch(T.sect(i))
                        {
                           case 0: 
                              m2_eff_var *= 1.2885;
                              break;
                           case 1: 
                              m2_eff_var *= 1.396;
                              break;
                        }
                     }
                  }
                  
                  if (charge == 1)
                  {
                     for (unsigned int j = 0; j < CentralityContainer.size; j++)
                     {
                        ptrackDecreasedAcceptance.weightEMC[j][ptrack.size] = 
                           ptrack.weightEMC[j][ptrack.size]*(1. - acc_var);
                        ptrackIncreasedAcceptance.weightEMC[j][ptrack.size] = Minimum(1.,
                           ptrack.weightEMC[j][ptrack.size]*(1. + acc_var));

                        ptrackDecreasedAcceptance.weightIdEMC[j][ptrack.size] = 
                           ptrack.weightIdEMC[j][ptrack.size]*(1. - acc_var);
                        ptrackIncreasedAcceptance.weightIdEMC[j][ptrack.size] = Minimum(1.,
                           ptrack.weightIdEMC[j][ptrack.size]*(1. + acc_var));

                        ptrackDecreasedM2Eff.weightIdEMC[j][ptrack.size] = 
                           ptrack.weightIdEMC[j][ptrack.size]*(1. - m2_eff_var);
                        ptrackIncreasedM2Eff.weightIdEMC[j][ptrack.size] = Minimum(1.,
                           ptrack.weightIdEMC[j][ptrack.size]*(1. + m2_eff_var));
                     }
                  }
                  else
                  {
                     for (unsigned int j = 0; j < CentralityContainer.size; j++)
                     {
                        ntrackDecreasedAcceptance.weightEMC[j][ntrack.size] = 
                           ntrack.weightEMC[j][ntrack.size]*(1. - acc_var);
                        ntrackIncreasedAcceptance.weightEMC[j][ntrack.size] = Minimum(1.,
                           ntrack.weightEMC[j][ntrack.size]*(1. + acc_var));

                        ntrackDecreasedAcceptance.weightIdEMC[j][ntrack.size] = 
                           ntrack.weightIdEMC[j][ntrack.size]*(1. - acc_var);
                        ntrackIncreasedAcceptance.weightIdEMC[j][ntrack.size] = Minimum(1.,
                           ntrack.weightIdEMC[j][ntrack.size]*(1. + acc_var));

                        ntrackDecreasedM2Eff.weightIdEMC[j][ntrack.size] = 
                           ntrack.weightIdEMC[j][ntrack.size]*(1. - m2_eff_var);
                        ntrackIncreasedM2Eff.weightIdEMC[j][ntrack.size] = Minimum(1.,
                           ntrack.weightIdEMC[j][ntrack.size]*(1. + m2_eff_var));
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
                     double acc_var;

                     if (zed >= 0) acc_var = ErrPropagation(pc2_sys, dceast0_sys, pc1e_sys);
                     else acc_var = ErrPropagation(pc2_sys, dceast1_sys, pc1e_sys);
                     
                     if (charge == 1)
                     {
                        ptrack.idPC2[ptrack.size] = PartId.noPID;
   
                        for (int j = 0; j < CentralityContainer.size; j++)
                        {
                           ptrack.weightPC2[j][ptrack.size] = 
                              ptrack.embPC2.get()[j]/ptrack.embDCPC1.get()[j];
                           
                           ptrackDecreasedAcceptance.weightPC2[j][ptrack.size] = 
                              ptrack.weightPC2[j][ptrack.size]*(1. - acc_var);
                           ptrackIncreasedAcceptance.weightPC2[j][ptrack.size] = Minimum(1.,
                              ptrack.weightPC2[j][ptrack.size]*(1. + acc_var));
                        }
                     }
                     else
                     {
                        ntrack.idPC2[ntrack.size] = PartId.noPID;
                        
                        for (int j = 0; j < CentralityContainer.size; j++)
                        {
                           ntrack.weightPC2[j][ntrack.size] = 
                              ntrack.embPC2.get()[j]/ntrack.embDCPC1.get()[j];

                           ntrackDecreasedAcceptance.weightPC2[j][ntrack.size] = 
                              ntrack.weightPC2[j][ntrack.size]*(1. - acc_var);
                           ntrackIncreasedAcceptance.weightPC2[j][ntrack.size] = Minimum(1.,
                              ntrack.weightPC2[j][ntrack.size]*(1. + acc_var));
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
                     double acc_var;

                     if (pc3phi < 0) 
                     {
                        if (zed >= 0) acc_var = ErrPropagation(pc3w_sys, dcwest0_sys, pc1w_sys);
                        else acc_var = ErrPropagation(pc3w_sys, dcwest1_sys, pc1w_sys);
                     }
                     else
                     {
                        if (zed >= 0) acc_var = ErrPropagation(pc3e_sys, dceast0_sys, pc1e_sys);
                        else acc_var = ErrPropagation(pc3e_sys, dceast1_sys, pc1e_sys);
                     }
                     
                     if (charge == 1)
                     {
                        ptrack.idPC3[ptrack.size] = PartId.noPID;
                        
                        for (int j = 0; j < CentralityContainer.size; j++)
                        {
                           ptrack.weightPC3[j][ptrack.size] = 
                              ptrack.embPC3.get()[j]/ptrack.embDCPC1.get()[j];

                           ptrackDecreasedAcceptance.weightPC3[j][ptrack.size] = 
                              ptrack.weightPC3[j][ptrack.size]*(1. - acc_var); 
                           ptrackIncreasedAcceptance.weightPC3[j][ptrack.size] = Minimum(1.,
                              ptrack.weightPC3[j][ptrack.size]*(1. + acc_var)); 
                        }
                     }
                     else
                     {
                        ntrack.idPC3[ntrack.size] = PartId.noPID;
                        
                        for (int j = 0; j < CentralityContainer.size; j++)
                        {
                           ntrack.weightPC3[j][ntrack.size] = 
                              ntrack.embPC3.get()[j]/ntrack.embDCPC1.get()[j];
                           
                           ntrackDecreasedAcceptance.weightPC3[j][ntrack.size] = 
                              ntrack.weightPC3[j][ntrack.size]*(1. - acc_var); 
                           ntrackIncreasedAcceptance.weightPC3[j][ntrack.size] = Minimum(1.,
                              ntrack.weightPC3[j][ntrack.size]*(1. + acc_var)); 
                        }
                     }
                  }
               }
            }
   
            //At least 1 detector is required
            if (charge == 1)
            {
               if (ptrack.idPC2[ptrack.size] == PartId.junk &&
                     ptrack.idPC3[ptrack.size] == PartId.junk &&
                     ptrack.idTOF[ptrack.size] == PartId.junk &&
                     ptrack.idEMC[ptrack.size] == PartId.junk) 
               {
                  continue;
               }
               else 
               {
                  ptrack.index[ptrack.size] = i;
                  ptrack.size += 1;
               }
            }
            else if (charge == -1)
            {
               if (ntrack.idPC2[ntrack.size] == PartId.junk &&
                     ntrack.idPC3[ntrack.size] == PartId.junk &&
                     ntrack.idTOF[ntrack.size] == PartId.junk &&
                     ntrack.idEMC[ntrack.size] == PartId.junk) 
               {
                  continue;
               }
               else 
               {
                  ntrack.index[ntrack.size] = i;
                  ntrack.size += 1;
               }
            }
         }

         for (int i = 0; i < ptrack.size; i++)
         {
            for (int j = 0; j < ntrack.size; j++)
            {
               const int ppc = ptrack.index[i]; //positive particle counter
               const int npc = ntrack.index[j]; //negative particle counter
               
               ptrack.mom[0] = T.mom(ppc)*sin(T.the0(ppc))*cos(T.phi0(ppc));
               ptrack.mom[1] = T.mom(ppc)*sin(T.the0(ppc))*sin(T.phi0(ppc));
               ptrack.mom[2] = T.mom(ppc)*cos(T.the0(ppc));
               
               ntrack.mom[0]   = T.mom(npc)*sin(T.the0(npc))*cos(T.phi0(npc));
               ntrack.mom[1]   = T.mom(npc)*sin(T.the0(npc))*sin(T.phi0(npc));
               ntrack.mom[2]   = T.mom(npc)*cos(T.the0(npc));
               
               if (IsGhostCut(T.zed(ppc) - T.zed(npc), 
                     T.alpha(ppc) - T.alpha(npc), 
                     T.phi(ppc) - T.phi(npc))) continue;
                  
               if (IsOneArmCut(T.phi(ppc), T.phi(npc))) continue;
               
               const double mass = GetMass(ptrack.mom, ntrack.mom, ptrack.mass, ntrack.mass);
               
               const double pt = sqrt((ptrack.mom[0] + ntrack.mom[0])*(ptrack.mom[0] + ntrack.mom[0]) +
                  (ptrack.mom[1] + ptrack.mom[1])*(ntrack.mom[1] + ntrack.mom[1]));

               Id1.pc2 = ptrack.idPC2[i];
               Id1.pc3 = ptrack.idPC3[i];
               Id1.emc = ptrack.idEMC[i];
               Id1.tof = ptrack.idTOF[i];

               Id2.pc2 = ntrack.idPC2[j];
               Id2.pc3 = ntrack.idPC3[j];
               Id2.emc = ntrack.idEMC[j];
               Id2.tof = ntrack.idTOF[j];

               if (IsnoPID(&Id1, &Id2))
               {
                  for (int k = 0; k < CentralityContainer.size; k++)
                  {
                     //probabilities of not registering dc-pc1 tracks
                     double no_track1_prob = Product(
                        1. - ptrack.weightPC2[k][i], 
                        1. - ptrack.weightPC3[k][i], 
                        1. - ptrack.weightTOF[k][i], 
                        1. - ptrack.weightEMC[k][i]);
                     
                     double no_track2_prob = Product(
                        1. - ntrack.weightPC2[k][j], 
                        1. - ntrack.weightPC3[k][j], 
                        1. - ntrack.weightTOF[k][j], 
                        1. - ntrack.weightEMC[k][j]);

                     double pair_weight = ptrack.embDCPC1.get()[k]*ntrack.embDCPC1.get()[k]*
                        (1. - no_track1_prob)*(1. - no_track2_prob);
      
                     invMDistrNoPID[k]->Fill(pt, mass, event_weight*pair_weight);

                     origPtVsPtDistr->Fill(pt, origPt, event_weight*pair_weight);

                     no_track1_prob = Product(
                        1. - ptrackDecreasedAcceptance.weightPC2[k][i], 
                        1. - ptrackDecreasedAcceptance.weightPC3[k][i], 
                        1. - ptrackDecreasedAcceptance.weightTOF[k][i], 
                        1. - ptrackDecreasedAcceptance.weightEMC[k][i]);
                     
                     no_track2_prob = Product(
                        1. - ntrackDecreasedAcceptance.weightPC2[k][j], 
                        1. - ntrackDecreasedAcceptance.weightPC3[k][j], 
                        1. - ntrackDecreasedAcceptance.weightTOF[k][j], 
                        1. - ntrackDecreasedAcceptance.weightEMC[k][j]);

                     pair_weight = ptrack.embDCPC1.get()[k]*ntrack.embDCPC1.get()[k]*
                        (1. - no_track1_prob)*(1. - no_track2_prob);

                     invMDistrNoPIDDecreasedAcceptance[k]->Fill(pt, mass, event_weight*pair_weight);

                     no_track1_prob = Product(
                        1. - ptrackIncreasedAcceptance.weightPC2[k][i], 
                        1. - ptrackIncreasedAcceptance.weightPC3[k][i], 
                        1. - ptrackIncreasedAcceptance.weightTOF[k][i], 
                        1. - ptrackIncreasedAcceptance.weightEMC[k][i]);
                     
                     no_track2_prob = Product(
                        1. - ntrackIncreasedAcceptance.weightPC2[k][j], 
                        1. - ntrackIncreasedAcceptance.weightPC3[k][j], 
                        1. - ntrackIncreasedAcceptance.weightTOF[k][j], 
                        1. - ntrackIncreasedAcceptance.weightEMC[k][j]);
                     
                     pair_weight = ptrack.embDCPC1.get()[k]*ntrack.embDCPC1.get()[k]*
                        (1. - no_track1_prob)*(1. - no_track2_prob);

                     invMDistrNoPIDIncreasedAcceptance[k]->Fill(pt, mass, event_weight*pair_weight);
                  }
               }
               else continue;

               if (Is1PID(&Id1, &Id2))
               {
                  for (int k = 0; k < CentralityContainer.size; k++)
                  {
                     double no_track1_prob = Product(
                        1. - ptrack.weightPC2[k][i], 
                        1. - ptrack.weightPC3[k][i], 
                        1. - ptrack.weightTOF[k][i], 
                        1. - ptrack.weightEMC[k][i]);

                     double no_id_track1_prob = 1. - ptrack.weightIdTOF[k][i];
                     
                     double no_track2_prob = Product(
                        1. - ntrack.weightPC2[k][j], 
                        1. - ntrack.weightPC3[k][j], 
                        1. - ntrack.weightTOF[k][j], 
                        1. - ntrack.weightEMC[k][j]);

                     double no_id_track2_prob = 1. - ntrack.weightIdTOF[k][j];
                     
                     double pair_weight = ptrack.embDCPC1.get()[k]*ntrack.embDCPC1.get()[k]*(1. -
                        (1. - (1. - no_id_track1_prob)*(1. - no_track2_prob))*
                        (1. - (1. - no_track1_prob)*(1. - no_id_track2_prob)));

                     invMDistr1PID[k]->Fill(pt, mass, event_weight*pair_weight);

                     no_track1_prob = Product(
                        1. - ptrackDecreasedAcceptance.weightPC2[k][i], 
                        1. - ptrackDecreasedAcceptance.weightPC3[k][i], 
                        1. - ptrackDecreasedAcceptance.weightTOF[k][i], 
                        1. - ptrackDecreasedAcceptance.weightEMC[k][i]);

                     no_id_track1_prob = 1. - ptrackDecreasedAcceptance.weightIdTOF[k][i];
                     
                     no_track2_prob = Product(
                        1. - ntrackDecreasedAcceptance.weightPC2[k][j], 
                        1. - ntrackDecreasedAcceptance.weightPC3[k][j], 
                        1. - ntrackDecreasedAcceptance.weightTOF[k][j], 
                        1. - ntrackDecreasedAcceptance.weightEMC[k][j]);

                     no_id_track2_prob = 1. - ntrackDecreasedAcceptance.weightIdTOF[k][j];
                     
                     pair_weight = ptrack.embDCPC1.get()[k]*ntrack.embDCPC1.get()[k]*(1. -
                        (1. - (1. - no_id_track1_prob)*(1. - no_track2_prob))*
                        (1. - (1. - no_track1_prob)*(1. - no_id_track2_prob)));

                     invMDistr1PIDDecreasedAcceptance[k]->Fill(pt, mass, event_weight*pair_weight);

                     no_track1_prob = Product(
                        1. - ptrackIncreasedAcceptance.weightPC2[k][i], 
                        1. - ptrackIncreasedAcceptance.weightPC3[k][i], 
                        1. - ptrackIncreasedAcceptance.weightTOF[k][i], 
                        1. - ptrackIncreasedAcceptance.weightEMC[k][i]);

                     no_id_track1_prob = 1. - ptrackIncreasedAcceptance.weightIdTOF[k][i];
                     
                     no_track2_prob = Product(
                        1. - ntrackIncreasedAcceptance.weightPC2[k][j], 
                        1. - ntrackIncreasedAcceptance.weightPC3[k][j], 
                        1. - ntrackIncreasedAcceptance.weightTOF[k][j], 
                        1. - ntrackIncreasedAcceptance.weightEMC[k][j]);
                     
                     no_id_track2_prob = 1. - ntrackIncreasedAcceptance.weightIdTOF[k][j];
                     
                     pair_weight = ptrack.embDCPC1.get()[k]*ntrack.embDCPC1.get()[k]*(1. -
                        (1. - (1. - no_id_track1_prob)*(1. - no_track2_prob))*
                        (1. - (1. - no_track1_prob)*(1. - no_id_track2_prob)));
                     
                     invMDistr1PIDIncreasedAcceptance[k]->Fill(pt, mass, event_weight*pair_weight);
                  }
               }

               if (Is2PID(&Id1, &Id2))
               {
                  for (int k = 0; k < CentralityContainer.size; k++)
                  {
                     double no_track1_prob = Product(
                        1. - ptrack.weightIdTOF[k][i], 
                        1. - ptrack.weightIdEMC[k][i]);
                     
                     double no_track2_prob = Product(
                        1. - ntrack.weightIdTOF[k][j], 
                        1. - ntrack.weightIdEMC[k][j]);
                     
                     double pair_weight = ptrack.embDCPC1.get()[k]*ntrack.embDCPC1.get()[k]*
                        (1. - no_track1_prob)*(1. - no_track2_prob);

                     invMDistr2PID[k]->Fill(pt, mass, event_weight*pair_weight);

                     no_track1_prob = Product(
                        1. - ptrackDecreasedAcceptance.weightIdTOF[k][i], 
                        1. - ptrackDecreasedAcceptance.weightIdEMC[k][i]);
                     
                     no_track2_prob = Product(
                        1. - ntrackDecreasedAcceptance.weightIdTOF[k][j], 
                        1. - ntrackDecreasedAcceptance.weightIdEMC[k][j]);
                     
                     pair_weight = ptrack.embDCPC1.get()[k]*ntrack.embDCPC1.get()[k]*
                        (1. - no_track1_prob)*(1. - no_track2_prob);

                     invMDistr2PIDDecreasedAcceptance[k]->Fill(pt, mass, event_weight*pair_weight);

                     no_track1_prob = Product(
                        1. - ptrackIncreasedAcceptance.weightIdTOF[k][i], 
                        1. - ptrackIncreasedAcceptance.weightIdEMC[k][i]);
                     
                     no_track2_prob = Product(
                        1. - ntrackIncreasedAcceptance.weightIdTOF[k][j], 
                        1. - ntrackIncreasedAcceptance.weightIdEMC[k][j]);
                     
                     pair_weight = ptrack.embDCPC1.get()[k]*ntrack.embDCPC1.get()[k]*
                        (1. - no_track1_prob)*(1. - no_track2_prob);

                     invMDistr2PIDIncreasedAcceptance[k]->Fill(pt, mass, event_weight*pair_weight);

                     no_track1_prob = Product(
                        1. - ptrack.weightIdTOF[k][i], 
                        1. - ptrackDecreasedM2Eff.weightIdEMC[k][i]);
                     
                     no_track2_prob = Product(
                        1. - ntrack.weightIdTOF[k][j], 
                        1. - ntrackDecreasedM2Eff.weightIdEMC[k][j]);
                     
                     pair_weight = ptrack.embDCPC1.get()[k]*ntrack.embDCPC1.get()[k]*
                        (1. - no_track1_prob)*(1. - no_track2_prob);

                     invMDistr2PIDDecreasedM2Eff[k]->Fill(pt, mass, event_weight*pair_weight);

                     no_track1_prob = Product(
                        1. - ptrack.weightIdTOF[k][i], 
                        1. - ptrackIncreasedM2Eff.weightIdEMC[k][i]);
                     
                     no_track2_prob = Product(
                        1. - ntrack.weightIdTOF[k][j], 
                        1. - ntrackIncreasedM2Eff.weightIdEMC[k][j]);
                     
                     pair_weight = ptrack.embDCPC1.get()[k]*ntrack.embDCPC1.get()[k]*
                        (1. - no_track1_prob)*(1. - no_track2_prob);

                     invMDistr2PIDIncreasedM2Eff[k]->Fill(pt, mass, event_weight*pair_weight);
                  }
               }
               else continue;

               if (IsTOF2PID(&Id1, &Id2))
               {
                  for (int k = 0; k < CentralityContainer.size; k++)
                  {
                     double pair_weight = 
                        ptrack.embDCPC1.get()[k]*ntrack.embDCPC1.get()[k]*
                        ptrack.weightIdTOF[k][i]*ntrack.weightIdTOF[k][j];
                     
                     invMDistrTOF2PID[k]->Fill(pt, mass, event_weight*pair_weight);
                     
                     pair_weight = 
                        ptrack.embDCPC1.get()[k]*ntrack.embDCPC1.get()[k]*
                        ptrackDecreasedAcceptance.weightIdTOF[k][i]*ntrackDecreasedAcceptance.weightIdTOF[k][j];
                     
                     invMDistrTOF2PIDDecreasedAcceptance[k]->Fill(pt, mass, event_weight*pair_weight);

                     pair_weight = 
                        ptrack.embDCPC1.get()[k]*ntrack.embDCPC1.get()[k]*
                        ptrackIncreasedAcceptance.weightIdTOF[k][i]*ntrackIncreasedAcceptance.weightIdTOF[k][j];
                     
                     invMDistrTOF2PIDIncreasedAcceptance[k]->Fill(pt, mass, event_weight*pair_weight);
                  }
               }

               if (IsEMC2PID(&Id1, &Id2))
               {
                  for (int k = 0; k < CentralityContainer.size; k++)
                  {
                     double pair_weight = 
                        ptrack.embDCPC1.get()[k]*ntrack.embDCPC1.get()[k]*
                        ptrack.weightIdEMC[k][i]*ntrack.weightIdEMC[k][j];

                     invMDistrEMC2PID[k]->Fill(pt, mass, event_weight*pair_weight);

                     pair_weight = 
                        ptrack.embDCPC1.get()[k]*ntrack.embDCPC1.get()[k]*
                        ptrackDecreasedAcceptance.weightIdEMC[k][i]*ntrackDecreasedAcceptance.weightIdEMC[k][j];

                     invMDistrEMC2PIDDecreasedAcceptance[k]->Fill(pt, mass, event_weight*pair_weight);

                     pair_weight = 
                        ptrack.embDCPC1.get()[k]*ntrack.embDCPC1.get()[k]*
                        ptrackIncreasedAcceptance.weightIdEMC[k][i]*ntrackIncreasedAcceptance.weightIdEMC[k][j];

                     invMDistrEMC2PIDIncreasedAcceptance[k]->Fill(pt, mass, event_weight*pair_weight);
                     
                     pair_weight = 
                        ptrack.embDCPC1.get()[k]*ntrack.embDCPC1.get()[k]*
                        ptrackDecreasedM2Eff.weightIdEMC[k][i]*ntrackDecreasedM2Eff.weightIdEMC[k][j];

                     invMDistrEMC2PIDDecreasedM2Eff[k]->Fill(pt, mass, event_weight*pair_weight);

                     pair_weight = 
                        ptrack.embDCPC1.get()[k]*ntrack.embDCPC1.get()[k]*
                        ptrackIncreasedM2Eff.weightIdEMC[k][i]*ntrackIncreasedM2Eff.weightIdEMC[k][j];

                     invMDistrEMC2PIDIncreasedM2Eff[k]->Fill(pt, mass, event_weight*pair_weight);
                  }
               }
            }//cylce by negative tracks
         }//cycle by positive traks 
      }//cycle by events
   };

   auto PbarCall = [&]()
   {
      ProgressBar pbar = ProgressBar("Fancy");
      while (!process_finished)
      {
         pbar.Print(ncalls/nevents);
         std::this_thread::sleep_for(std::chrono::milliseconds(20));
      }
      pbar.Print(1.);
   };
   
   std::thread pbar_thread(PbarCall);
   tp.Process(ProcessMP);
   process_finished = true;
   pbar_thread.join();
}

void AnalyzeTrackPair()
{
   if (Par.doUseWeightFunc) CheckInputFile("../input/Spectra/" + 
      Par.system + "/" + Par.origParticle.name + ".txt");
   
   for (long unsigned int i = 0; i < Par.daughter1Queue.size(); i++)
   {
      const int daughter1_iter = ParticleProperties.iter_map[Par.daughter1Queue[i]];
      const int daughter2_iter = ParticleProperties.iter_map[Par.daughter2Queue[i]];   
      
      const std::string decay_channel = 
         ParticleProperties.sname[daughter1_iter] + 
         ParticleProperties.sname[daughter2_iter];
      
      for (std::string detector : Par.detectors)
      {
         CheckInputFile("../input/Embedding/" + Par.runName + "/" + 
            ParticleProperties.abs_name[daughter1_iter] + "_" + detector + ".txt");
         CheckInputFile("../input/Embedding/" + Par.runName + "/" + 
            ParticleProperties.abs_name[daughter2_iter] + "_" + detector + ".txt");   
      }
      
      for (std::string magf : Par.magfQueue)
      {
         CheckInputFile("../input/M2Eff/" + Par.runName + "/" +
               ParticleProperties.name[daughter1_iter] + "_EMCale" + magf + ".txt");
         CheckInputFile("../input/M2Eff/" + Par.runName + "/" +
               ParticleProperties.name[daughter1_iter] + "_EMCalw" + magf + ".txt");
         CheckInputFile("../input/Systematics/" + Par.runName + "/m2eff_" +
               ParticleProperties.name[daughter1_iter] + "_EMCale" + magf + ".txt");
         CheckInputFile("../input/Systematics/" + Par.runName + "/m2eff_" +
               ParticleProperties.name[daughter1_iter] + "_EMCalw" + magf + ".txt");

         CheckInputFile("../input/M2Eff/" + Par.runName + "/" +
               ParticleProperties.name[daughter2_iter] + "_EMCale" + magf + ".txt");
         CheckInputFile("../input/M2Eff/" + Par.runName + "/" +
               ParticleProperties.name[daughter2_iter] + "_EMCalw" + magf + ".txt");
         CheckInputFile("../input/Systematics/" + Par.runName + "/m2eff_" +
               ParticleProperties.name[daughter2_iter] + "_EMCale" + magf + ".txt");
         CheckInputFile("../input/Systematics/" + Par.runName + "/m2eff_" +
               ParticleProperties.name[daughter2_iter] + "_EMCalw" + magf + ".txt");
         for (std::string auxName : Par.auxNameQueue)
         {
            const std::string input_file_name = "../data/" + Par.runName + "/Resonances/" + 
               Par.origParticle.name + "_" + decay_channel + magf + auxName + ".root";
            
            CheckInputFile(input_file_name);
         }
      }
   }

   CheckInputFile("../input/Systematics/" + Par.runName + "/acceptance.txt");
   
   int procNum = 1;
   for (double ptDeviation : Par.ptDeviationQueue)
   {
      thrContainer thrContainer;

      const double min_inv_mass = 0.;

      for (long unsigned int j = 0; j < CentralityContainer.nameWithoutPercent.size(); j++)
      {
         std::string dir = CentralityContainer.nameWithoutPercent[j];
         
         thrContainer.invMDistrNoPID[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "noPID_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "noPID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));
         
         thrContainer.invMDistrNoPIDDecreasedAcceptance[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "noPID_A-_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
                "noPID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
                Par.invMNBins, min_inv_mass, 4., dir));
         
         thrContainer.invMDistrNoPIDIncreasedAcceptance[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "noPID_A+_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "noPID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));

         thrContainer.invMDistr1PID[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "1PID_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "1PID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));
         
         thrContainer.invMDistr1PIDDecreasedAcceptance[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "1PID_A-_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "1PID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));
         
         thrContainer.invMDistr1PIDIncreasedAcceptance[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "1PID_A+_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "1PID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));

         thrContainer.invMDistr2PID[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "2PID_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "2PID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));
         
         thrContainer.invMDistr2PIDDecreasedAcceptance[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "2PID_A-_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "2PID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));
         
         thrContainer.invMDistr2PIDIncreasedAcceptance[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "2PID_A+_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "2PID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));
         
         thrContainer.invMDistr2PIDDecreasedM2Eff[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "2PID_m2Eff-_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "2PID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));
         
         thrContainer.invMDistr2PIDIncreasedM2Eff[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "2PID_m2Eff+_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "2PID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));
         
         thrContainer.invMDistrTOF2PID[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "TOF2PID_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "TOF2PID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));
         
         thrContainer.invMDistrTOF2PIDDecreasedAcceptance[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "TOF2PID_A-_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "TOF2PID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));
         
         thrContainer.invMDistrTOF2PIDIncreasedAcceptance[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "TOF2PID_A+_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "TOF2PID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));

         thrContainer.invMDistrEMC2PID[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "EMC2PID_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "EMC2PID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));
         
         thrContainer.invMDistrEMC2PIDDecreasedAcceptance[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "EMC2PID_A-_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "EMC2PID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));
         
         thrContainer.invMDistrEMC2PIDIncreasedAcceptance[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "EMC2PID_A+_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "EMC2PID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));
         
         thrContainer.invMDistrEMC2PIDDecreasedM2Eff[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "EMC2PID_m2Eff-_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "EMC2PID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));
         
         thrContainer.invMDistrEMC2PIDIncreasedM2Eff[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "EMC2PID_m2Eff+_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "EMC2PID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));
      }
      
      for (long unsigned int i = 0; i < Par.daughter1Queue.size(); i++)
      {
         for (std::string magf : Par.magfQueue)
         {   
            for (std::string auxName : Par.auxNameQueue)
            {
               AnalyzeConfiguration(&thrContainer, Par.daughter1Queue[i], Par.daughter2Queue[i], 
                  magf, auxName, ptDeviation, procNum);
               procNum++;
            }
         }   
      }
      
      system(("mkdir -p ../../analysis/data/phenix_sim/" + Par.runName).c_str());

      std::string ptDeviation_name;
      if (abs(ptDeviation - 1) < 1e-7) ptDeviation_name = "";
      else ptDeviation_name = "_pt" + DtoStr(ptDeviation, 3);
   
      const std::string output_file_name = 
         "../../analysis/data/phenix_sim/" + 
         Par.runName + "/" + Par.origParticle.name + 
         ptDeviation_name + ".root";
      
      ThrObjHolder.Write(output_file_name);
      ThrObjHolder.Clear();

      PrintInfo("File " + output_file_name + " was written");
   }
}

int main()
{
   AnalyzeTrackPair();
   return 0;
}

#endif /* ANALYZE_RESONANCE_CPP */
