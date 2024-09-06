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

//Analyze specific track pair configuration
void AnalyzeConfiguration(ThrContainerStruct *ThrContainer, const std::string& daughter1, 
                          const std::string& daughter2, const std::string magf, 
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
      PartArray ptrack, ntrack;

      //arrays for keeping weights for varied acceptancies
      PartArray ptrackDecreasedAcceptance, ptrackIncreasedAcceptance, ntrackDecreasedAcceptance, ntrackIncreasedAcceptance;
      PartArray ptrackDecreasedM2Eff, ptrackIncreasedM2Eff, ntrackDecreasedM2Eff, ntrackIncreasedM2Eff;

      //arrays for keeping weights for varied /epsilon_{m^2}
      
      //setting daughter particles
      ptrack.iter = ParticleProperties.iter_map[daughter1];
      ntrack.iter = ParticleProperties.iter_map[daughter2];
      
      ptrack.origId = abs(ParticleProperties.id[ptrack.iter]);
      ntrack.origId = abs(ParticleProperties.id[ntrack.iter]);

      idContainer Id1, Id2;

      Id1.origId = ptrack.origId;
      Id2.origId = ntrack.origId;

      ptrack.geant_id = ParticleProperties.geant_id[ptrack.iter];
      ntrack.geant_id = ParticleProperties.geant_id[ntrack.iter];

      ptrack.mass = ParticleProperties.mass[ptrack.iter];
      ntrack.mass = ParticleProperties.mass[ntrack.iter];

      //setting embedding
      ptrack.dc_pc1_emb = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
         ParticleProperties.abs_name[ptrack.iter] + "_dc_pc1.txt", 
         CentralityContainer.size));
      ptrack.pc2_emb = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
         ParticleProperties.abs_name[ptrack.iter] + "_pc2.txt", 
         CentralityContainer.size));
      ptrack.pc3_emb = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
         ParticleProperties.abs_name[ptrack.iter] + "_pc3.txt", 
         CentralityContainer.size));
      ptrack.tofe_emb = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
         ParticleProperties.abs_name[ptrack.iter] + "_tofe.txt", 
         CentralityContainer.size));
      ptrack.tofw_emb = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
         ParticleProperties.abs_name[ptrack.iter] + "_tofw.txt", 
         CentralityContainer.size));

      ntrack.dc_pc1_emb = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
         ParticleProperties.abs_name[ntrack.iter] + "_dc_pc1.txt", 
         CentralityContainer.size));
      ntrack.pc2_emb = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
         ParticleProperties.abs_name[ntrack.iter] + "_pc2.txt", 
         CentralityContainer.size));
      ntrack.pc3_emb = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
         ParticleProperties.abs_name[ntrack.iter] + "_pc3.txt", 
         CentralityContainer.size));
      ntrack.tofe_emb = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
         ParticleProperties.abs_name[ntrack.iter] + "_tofe.txt", 
         CentralityContainer.size));
      ntrack.tofw_emb = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
         ParticleProperties.abs_name[ntrack.iter] + "_tofw.txt", 
         CentralityContainer.size));
      
      for (int i = 0; i < 4; i++)
      {
         ptrack.emcale_emb[i] = std::unique_ptr<double>(
            ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
            ParticleProperties.abs_name[ptrack.iter] + "_emcale" + std::to_string(i) + ".txt", 
            CentralityContainer.size));
         ptrack.emcalw_emb[i] = std::unique_ptr<double>(
            ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
            ParticleProperties.abs_name[ptrack.iter] + "_emcalw" + std::to_string(i) + ".txt", 
            CentralityContainer.size));
         
         ntrack.emcale_emb[i] = std::unique_ptr<double>(
            ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
            ParticleProperties.abs_name[ntrack.iter] + "_emcale" + std::to_string(i) + ".txt", 
            CentralityContainer.size));
         ntrack.emcalw_emb[i] = std::unique_ptr<double>(
            ReadFileIntoArray("../input/Embedding/" + Par.runName + "/" + 
            ParticleProperties.abs_name[ntrack.iter] + "_emcalw" + std::to_string(i) + ".txt", 
            CentralityContainer.size));
      }

      ptrack.emcale_m2_eff = std::unique_ptr<double>(
         ReadFileIntoArray("../input/M2Eff/" + Par.runName + "/" +
         ParticleProperties.name[ptrack.iter] + "_EMCale" + magf + ".txt", 2));
      ptrack.emcalw_m2_eff = std::unique_ptr<double>(
         ReadFileIntoArray("../input/M2Eff/" + Par.runName + "/" +
         ParticleProperties.name[ptrack.iter] + "_EMCalw" + magf + ".txt", 4));
      
      ntrack.emcale_m2_eff = std::unique_ptr<double>(
         ReadFileIntoArray("../input/M2Eff/" + Par.runName + "/" +
         ParticleProperties.name[ntrack.iter] + "_EMCale" + magf + ".txt", 2));
      ntrack.emcalw_m2_eff = std::unique_ptr<double>(
         ReadFileIntoArray("../input/M2Eff/" + Par.runName + "/" +
         ParticleProperties.name[ntrack.iter] + "_EMCalw" + magf + ".txt", 4));

      ptrack.emcale_m2_eff_sys = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Systematics/" + Par.runName + "/m2eff_" +
         ParticleProperties.name[ptrack.iter] + "_EMCale" + magf + ".txt", 2));
      ptrack.emcalw_m2_eff_sys = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Systematics/" + Par.runName + "/m2eff_" +
         ParticleProperties.name[ptrack.iter] + "_EMCalw" + magf + ".txt", 4));
      
      ntrack.emcale_m2_eff_sys = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Systematics/" + Par.runName + "/m2eff_" +
         ParticleProperties.name[ntrack.iter] + "_EMCale" + magf + ".txt", 2));
      ntrack.emcalw_m2_eff_sys = std::unique_ptr<double>(
         ReadFileIntoArray("../input/Systematics/" + Par.runName + "/m2eff_" +
         ParticleProperties.name[ntrack.iter] + "_EMCalw" + magf + ".txt", 4));

      std::shared_ptr<TH1F> origPtDistr = ThrContainer->origPtDistr.Get();
      
      std::shared_ptr<TH2F> origPtVsPtDistr = ThrContainer->origPtVsPtDistr.Get();
      
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
         invMDistrNoPID[i] = ThrContainer->invMDistrNoPID[i]->Get();
         invMDistrNoPIDDecreasedAcceptance[i] = ThrContainer->invMDistrNoPIDDecreasedAcceptance[i]->Get();
         invMDistrNoPIDIncreasedAcceptance[i] = ThrContainer->invMDistrNoPIDIncreasedAcceptance[i]->Get();

         invMDistr1PID[i] = ThrContainer->invMDistr1PID[i]->Get();
         invMDistr1PIDDecreasedAcceptance[i] = ThrContainer->invMDistr1PIDDecreasedAcceptance[i]->Get();
         invMDistr1PIDIncreasedAcceptance[i] = ThrContainer->invMDistr1PIDIncreasedAcceptance[i]->Get();

         invMDistr2PID[i] = ThrContainer->invMDistr2PID[i]->Get();
         invMDistr2PIDDecreasedAcceptance[i] = ThrContainer->invMDistr2PIDDecreasedAcceptance[i]->Get();
         invMDistr2PIDIncreasedAcceptance[i] = ThrContainer->invMDistr2PIDIncreasedAcceptance[i]->Get();
         invMDistr2PIDDecreasedM2Eff[i] = ThrContainer->invMDistr2PIDDecreasedM2Eff[i]->Get();
         invMDistr2PIDIncreasedM2Eff[i] = ThrContainer->invMDistr2PIDIncreasedM2Eff[i]->Get();

         invMDistrTOF2PID[i] = ThrContainer->invMDistrTOF2PID[i]->Get();
         invMDistrTOF2PIDDecreasedAcceptance[i] = ThrContainer->invMDistrTOF2PIDDecreasedAcceptance[i]->Get();
         invMDistrTOF2PIDIncreasedAcceptance[i] = ThrContainer->invMDistrTOF2PIDIncreasedAcceptance[i]->Get();

         invMDistrEMC2PID[i] = ThrContainer->invMDistrEMC2PID[i]->Get();
         invMDistrEMC2PIDDecreasedAcceptance[i] = ThrContainer->invMDistrEMC2PIDDecreasedAcceptance[i]->Get();
         invMDistrEMC2PIDIncreasedAcceptance[i] = ThrContainer->invMDistrEMC2PIDIncreasedAcceptance[i]->Get();
         invMDistrEMC2PIDDecreasedM2Eff[i] = ThrContainer->invMDistrEMC2PIDDecreasedM2Eff[i]->Get();
         invMDistrEMC2PIDIncreasedM2Eff[i] = ThrContainer->invMDistrEMC2PIDIncreasedM2Eff[i]->Get();
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
                     ptrack.tof_id[ptrack.size] = id;
                     
                     for (int j = 0; j < CentralityContainer.size; j++)
                     {
                        ptrack.tof_weight[j][ptrack.size] = 
                           ptrack.tofe_emb.get()[j]/ptrack.dc_pc1_emb.get()[j];
                        ptrackDecreasedAcceptance.tof_weight[j][ptrack.size] = 
                           ptrack.tof_weight[j][ptrack.size]*(1. - acc_var);
                        ptrackIncreasedAcceptance.tof_weight[j][ptrack.size] = 
                           Minimum(1., ptrack.tof_weight[j][ptrack.size]*(1.+acc_var));

                        if (id == ptrack.origId) 
                        {
                           ptrack.tof_id_weight[j][ptrack.size] = 
                              ptrack.tof_weight[j][ptrack.size];
                           ptrackDecreasedAcceptance.tof_id_weight[j][ptrack.size] = 
                              ptrack.tof_id_weight[j][ptrack.size]*(1. - acc_var);
                           ptrackIncreasedAcceptance.tof_id_weight[j][ptrack.size] = 
                              Minimum(1., ptrack.tof_id_weight[j][ptrack.size]*(1.+acc_var));
                        }   
                     }
                  }
                  else
                  {
                     ntrack.tof_id[ntrack.size] = id;
                     
                     for (int j = 0; j < CentralityContainer.size; j++)
                     {
                        ntrack.tof_weight[j][ntrack.size] = 
                           ntrack.tofe_emb.get()[j]/ntrack.dc_pc1_emb.get()[j];
                        ntrackDecreasedAcceptance.tof_weight[j][ntrack.size] = 
                           ntrack.tof_weight[j][ntrack.size]*(1. - acc_var);
                        ntrackIncreasedAcceptance.tof_weight[j][ntrack.size] = 
                           Minimum(1., ntrack.tof_weight[j][ntrack.size]*(1.+acc_var));
                        
                        if (id == ntrack.origId)
                        {
                           ntrack.tof_id_weight[j][ntrack.size] = 
                              ntrack.tof_weight[j][ntrack.size];
                           ntrackDecreasedAcceptance.tof_id_weight[j][ntrack.size] = 
                              ntrack.tof_id_weight[j][ntrack.size]*(1. - acc_var);
                           ntrackIncreasedAcceptance.tof_id_weight[j][ntrack.size] = Minimum(1.,
                              ntrack.tof_id_weight[j][ntrack.size]*(1.+acc_var));
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
                     ptrack.tof_id[ptrack.size] = id;
                     
                     for (unsigned int j = 0; j < CentralityContainer.size; j++)
                     {
                        ptrack.tof_weight[j][ptrack.size] = 
                           Par.correctionTOFw*ptrack.tofw_emb.get()[j]/ptrack.dc_pc1_emb.get()[j];
                        ptrackDecreasedAcceptance.tof_weight[j][ptrack.size] = 
                           ptrack.tof_weight[j][ptrack.size]*(1. - acc_var);
                        ptrackIncreasedAcceptance.tof_weight[j][ptrack.size] = Minimum(1.,
                           ptrack.tof_weight[j][ptrack.size]*(1. + acc_var));

                        if (id == ptrack.origId) 
                        {
                           ptrack.tof_id_weight[j][ptrack.size] = 
                               ptrack.tof_weight[j][ptrack.size];
                           ptrackDecreasedAcceptance.tof_id_weight[j][ptrack.size] = 
                               ptrack.tof_weight[j][ptrack.size]*(1. - acc_var);
                           ptrackIncreasedAcceptance.tof_id_weight[j][ptrack.size] = Minimum(1.,
                              ptrack.tof_weight[j][ptrack.size]*(1.+acc_var));
                        }
                     }
                  }
                  else
                  {
                     ntrack.tof_id[ntrack.size] = id;
                     
                     for (unsigned int j = 0; j < CentralityContainer.size; j++)
                     {
                        ntrack.tof_weight[j][ntrack.size] = 
                           Par.correctionTOFw*ntrack.tofw_emb.get()[j]/ntrack.dc_pc1_emb.get()[j];
                        ntrackDecreasedAcceptance.tof_weight[j][ntrack.size] = 
                           ntrack.tof_weight[j][ntrack.size]*(1. - acc_var);
                        ntrackIncreasedAcceptance.tof_weight[j][ntrack.size] = Minimum(1.,
                           ntrack.tof_weight[j][ntrack.size]*(1.+acc_var));

                        if (id == ntrack.origId) 
                        {
                           ntrack.tof_id_weight[j][ntrack.size] = 
                              ntrack.tof_weight[j][ntrack.size];
                           ntrackDecreasedAcceptance.tof_id_weight[j][ntrack.size] = 
                              ntrack.tof_id_weight[j][ntrack.size]*(1. - acc_var);
                           ntrackIncreasedAcceptance.tof_id_weight[j][ntrack.size] = Minimum(1.,
                              ntrack.tof_id_weight[j][ntrack.size]*(1.+acc_var));
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

                  if (charge == 1) ptrack.emc_id[ptrack.size] = PartId.noPID;
                  else ntrack.emc_id[ntrack.size] = PartId.noPID;

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
                           ptrack.emc_weight[j][ptrack.size] = 
                              ptrack.emcale_emb[T.sect(i)].get()[j]/ptrack.dc_pc1_emb.get()[j];
                        }
         
                        if (T.particle_id(i) == ptrack.geant_id && T.sect(i) > 1)
                        {
                           ptrack.emc_id[ptrack.size] = ptrack.origId;
                           for (unsigned int j = 0; j < CentralityContainer.size; j++)
                           {
                              ptrack.emc_id_weight[j][ptrack.size] = 
                                 ptrack.emc_weight[j][ptrack.size]*
                                 ptrack.emcale_m2_eff.get()[T.sect(i)-2]*
                                 GetEMCalEid(pt, ptrack.emc_id[ptrack.size], charge, phi, T.sect(i));
                              m2_eff_var = ptrack.emcale_m2_eff_sys.get()[T.sect(i)-2];
                           }
                        }
                     }
                     else
                     {   
                        for (unsigned int j = 0; j < CentralityContainer.size; j++)
                        {
                           ntrack.emc_weight[j][ntrack.size] = 
                              ntrack.emcale_emb[T.sect(i)].get()[j]/ntrack.dc_pc1_emb.get()[j];
                        }
                           
                        if (T.particle_id(i) == ntrack.geant_id && T.sect(i) > 1)
                        {
                           ntrack.emc_id[ntrack.size] = ntrack.origId;
                           for (unsigned int j = 0; j < CentralityContainer.size; j++)
                           {
                              ntrack.emc_id_weight[j][ntrack.size] = 
                                 ntrack.emc_weight[j][ntrack.size]*
                                 ntrack.emcale_m2_eff.get()[T.sect(i)-2]*
                                 GetEMCalEid(pt, ntrack.emc_id[ntrack.size], charge, phi, T.sect(i));
                              m2_eff_var = ntrack.emcale_m2_eff_sys.get()[T.sect(i)-2];
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
                           ptrack.emc_weight[j][ptrack.size] = 
                              ptrack.emcalw_emb[T.sect(i)].get()[j]/ptrack.dc_pc1_emb.get()[j];
                        }
                        
                        if (T.particle_id(i) == ptrack.geant_id)
                        {
                           ptrack.emc_id[ptrack.size] = ptrack.origId;
                           for (unsigned int j = 0; j < CentralityContainer.size; j++)
                           {
                              ptrack.emc_id_weight[j][ptrack.size] = 
                                 ptrack.emc_weight[j][ptrack.size]*
                                 ptrack.emcalw_m2_eff.get()[T.sect(i)]*
                                 GetEMCalEid(pt, ptrack.emc_id[ptrack.size], charge, phi, T.sect(i));
                              m2_eff_var = ptrack.emcalw_m2_eff_sys.get()[T.sect(i)];
                           }
                        }
                     }
                     else
                     {
                        for (unsigned int j = 0; j < CentralityContainer.size; j++)
                        {
                           ntrack.emc_weight[j][ntrack.size] = 
                              ntrack.emcalw_emb[T.sect(i)].get()[j]/ntrack.dc_pc1_emb.get()[j];
                        }
                        
                        if (T.particle_id(i) == ntrack.geant_id)
                        {
                           ntrack.emc_id[ntrack.size] = ntrack.origId;
                           for (unsigned int j = 0; j < CentralityContainer.size; j++)
                           {
                              ntrack.emc_id_weight[j][ntrack.size] = 
                                 ntrack.emc_weight[j][ntrack.size]*
                                 ntrack.emcalw_m2_eff.get()[T.sect(i)]*
                                 GetEMCalEid(pt, ntrack.emc_id[ntrack.size], charge, phi, T.sect(i));
                              m2_eff_var = ntrack.emcalw_m2_eff_sys.get()[T.sect(i)];
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
                        ptrackDecreasedAcceptance.emc_weight[j][ptrack.size] = 
                           ptrack.emc_weight[j][ptrack.size]*(1. - acc_var);
                        ptrackIncreasedAcceptance.emc_weight[j][ptrack.size] = Minimum(1.,
                           ptrack.emc_weight[j][ptrack.size]*(1. + acc_var));

                        ptrackDecreasedAcceptance.emc_id_weight[j][ptrack.size] = 
                           ptrack.emc_id_weight[j][ptrack.size]*(1. - acc_var);
                        ptrackIncreasedAcceptance.emc_id_weight[j][ptrack.size] = Minimum(1.,
                           ptrack.emc_id_weight[j][ptrack.size]*(1. + acc_var));

                        ptrackDecreasedM2Eff.emc_id_weight[j][ptrack.size] = 
                           ptrack.emc_id_weight[j][ptrack.size]*(1. - m2_eff_var);
                        ptrackIncreasedM2Eff.emc_id_weight[j][ptrack.size] = Minimum(1.,
                           ptrack.emc_id_weight[j][ptrack.size]*(1. + m2_eff_var));
                     }
                  }
                  else
                  {
                     for (unsigned int j = 0; j < CentralityContainer.size; j++)
                     {
                        ntrackDecreasedAcceptance.emc_weight[j][ntrack.size] = 
                           ntrack.emc_weight[j][ntrack.size]*(1. - acc_var);
                        ntrackIncreasedAcceptance.emc_weight[j][ntrack.size] = Minimum(1.,
                           ntrack.emc_weight[j][ntrack.size]*(1. + acc_var));

                        ntrackDecreasedAcceptance.emc_id_weight[j][ntrack.size] = 
                           ntrack.emc_id_weight[j][ntrack.size]*(1. - acc_var);
                        ntrackIncreasedAcceptance.emc_id_weight[j][ntrack.size] = Minimum(1.,
                           ntrack.emc_id_weight[j][ntrack.size]*(1. + acc_var));

                        ntrackDecreasedM2Eff.emc_id_weight[j][ntrack.size] = 
                           ntrack.emc_id_weight[j][ntrack.size]*(1. - m2_eff_var);
                        ntrackIncreasedM2Eff.emc_id_weight[j][ntrack.size] = Minimum(1.,
                           ntrack.emc_id_weight[j][ntrack.size]*(1. + m2_eff_var));
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
                        ptrack.pc2_id[ptrack.size] = PartId.noPID;
   
                        for (int j = 0; j < CentralityContainer.size; j++)
                        {
                           ptrack.pc2_weight[j][ptrack.size] = 
                              ptrack.pc2_emb.get()[j]/ptrack.dc_pc1_emb.get()[j];
                           
                           ptrackDecreasedAcceptance.pc2_weight[j][ptrack.size] = 
                              ptrack.pc2_weight[j][ptrack.size]*(1. - acc_var);
                           ptrackIncreasedAcceptance.pc2_weight[j][ptrack.size] = Minimum(1.,
                              ptrack.pc2_weight[j][ptrack.size]*(1. + acc_var));
                        }
                     }
                     else
                     {
                        ntrack.pc2_id[ntrack.size] = PartId.noPID;
                        
                        for (int j = 0; j < CentralityContainer.size; j++)
                        {
                           ntrack.pc2_weight[j][ntrack.size] = 
                              ntrack.pc2_emb.get()[j]/ntrack.dc_pc1_emb.get()[j];

                           ntrackDecreasedAcceptance.pc2_weight[j][ntrack.size] = 
                              ntrack.pc2_weight[j][ntrack.size]*(1. - acc_var);
                           ntrackIncreasedAcceptance.pc2_weight[j][ntrack.size] = Minimum(1.,
                              ntrack.pc2_weight[j][ntrack.size]*(1. + acc_var));
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
                        ptrack.pc3_id[ptrack.size] = PartId.noPID;
                        
                        for (int j = 0; j < CentralityContainer.size; j++)
                        {
                           ptrack.pc3_weight[j][ptrack.size] = 
                              ptrack.pc3_emb.get()[j]/ptrack.dc_pc1_emb.get()[j];

                           ptrackDecreasedAcceptance.pc3_weight[j][ptrack.size] = 
                              ptrack.pc3_weight[j][ptrack.size]*(1. - acc_var); 
                           ptrackIncreasedAcceptance.pc3_weight[j][ptrack.size] = Minimum(1.,
                              ptrack.pc3_weight[j][ptrack.size]*(1. + acc_var)); 
                        }
                     }
                     else
                     {
                        ntrack.pc3_id[ntrack.size] = PartId.noPID;
                        
                        for (int j = 0; j < CentralityContainer.size; j++)
                        {
                           ntrack.pc3_weight[j][ntrack.size] = 
                              ntrack.pc3_emb.get()[j]/ntrack.dc_pc1_emb.get()[j];
                           
                           ntrackDecreasedAcceptance.pc3_weight[j][ntrack.size] = 
                              ntrack.pc3_weight[j][ntrack.size]*(1. - acc_var); 
                           ntrackIncreasedAcceptance.pc3_weight[j][ntrack.size] = Minimum(1.,
                              ntrack.pc3_weight[j][ntrack.size]*(1. + acc_var)); 
                        }
                     }
                  }
               }
            }
   
            //At least 1 detector is required
            if (charge == 1)
            {
               if (ptrack.pc2_id[ptrack.size] == PartId.junk &&
                     ptrack.pc3_id[ptrack.size] == PartId.junk &&
                     ptrack.tof_id[ptrack.size] == PartId.junk &&
                     ptrack.emc_id[ptrack.size] == PartId.junk) 
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
               if (ntrack.pc2_id[ntrack.size] == PartId.junk &&
                     ntrack.pc3_id[ntrack.size] == PartId.junk &&
                     ntrack.tof_id[ntrack.size] == PartId.junk &&
                     ntrack.emc_id[ntrack.size] == PartId.junk) 
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

               Id1.pc2 = ptrack.pc2_id[i];
               Id1.pc3 = ptrack.pc3_id[i];
               Id1.emc = ptrack.emc_id[i];
               Id1.tof = ptrack.tof_id[i];

               Id2.pc2 = ntrack.pc2_id[j];
               Id2.pc3 = ntrack.pc3_id[j];
               Id2.emc = ntrack.emc_id[j];
               Id2.tof = ntrack.tof_id[j];

               if (IsnoPID(&Id1, &Id2))
               {
                  for (int k = 0; k < CentralityContainer.size; k++)
                  {
                     //probabilities of not registering dc-pc1 tracks
                     double no_track1_prob = Product(
                        1. - ptrack.pc2_weight[k][i], 
                        1. - ptrack.pc3_weight[k][i], 
                        1. - ptrack.tof_weight[k][i], 
                        1. - ptrack.emc_weight[k][i]);
                     
                     double no_track2_prob = Product(
                        1. - ntrack.pc2_weight[k][j], 
                        1. - ntrack.pc3_weight[k][j], 
                        1. - ntrack.tof_weight[k][j], 
                        1. - ntrack.emc_weight[k][j]);

                     double pair_weight = ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*
                        (1. - no_track1_prob)*(1. - no_track2_prob);
      
                     invMDistrNoPID[k]->Fill(pt, mass, event_weight*pair_weight);

                     origPtVsPtDistr->Fill(pt, origPt, event_weight*pair_weight);

                     no_track1_prob = Product(
                        1. - ptrackDecreasedAcceptance.pc2_weight[k][i], 
                        1. - ptrackDecreasedAcceptance.pc3_weight[k][i], 
                        1. - ptrackDecreasedAcceptance.tof_weight[k][i], 
                        1. - ptrackDecreasedAcceptance.emc_weight[k][i]);
                     
                     no_track2_prob = Product(
                        1. - ntrackDecreasedAcceptance.pc2_weight[k][j], 
                        1. - ntrackDecreasedAcceptance.pc3_weight[k][j], 
                        1. - ntrackDecreasedAcceptance.tof_weight[k][j], 
                        1. - ntrackDecreasedAcceptance.emc_weight[k][j]);

                     pair_weight = ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*
                        (1. - no_track1_prob)*(1. - no_track2_prob);

                     invMDistrNoPIDDecreasedAcceptance[k]->Fill(pt, mass, event_weight*pair_weight);

                     no_track1_prob = Product(
                        1. - ptrackIncreasedAcceptance.pc2_weight[k][i], 
                        1. - ptrackIncreasedAcceptance.pc3_weight[k][i], 
                        1. - ptrackIncreasedAcceptance.tof_weight[k][i], 
                        1. - ptrackIncreasedAcceptance.emc_weight[k][i]);
                     
                     no_track2_prob = Product(
                        1. - ntrackIncreasedAcceptance.pc2_weight[k][j], 
                        1. - ntrackIncreasedAcceptance.pc3_weight[k][j], 
                        1. - ntrackIncreasedAcceptance.tof_weight[k][j], 
                        1. - ntrackIncreasedAcceptance.emc_weight[k][j]);
                     
                     pair_weight = ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*
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
                        1. - ptrack.pc2_weight[k][i], 
                        1. - ptrack.pc3_weight[k][i], 
                        1. - ptrack.tof_weight[k][i], 
                        1. - ptrack.emc_weight[k][i]);

                     double no_id_track1_prob = 1. - ptrack.tof_id_weight[k][i];
                     
                     double no_track2_prob = Product(
                        1. - ntrack.pc2_weight[k][j], 
                        1. - ntrack.pc3_weight[k][j], 
                        1. - ntrack.tof_weight[k][j], 
                        1. - ntrack.emc_weight[k][j]);

                     double no_id_track2_prob = 1. - ntrack.tof_id_weight[k][j];
                     
                     double pair_weight = ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*(1. -
                        (1. - (1. - no_id_track1_prob)*(1. - no_track2_prob))*
                        (1. - (1. - no_track1_prob)*(1. - no_id_track2_prob)));

                     invMDistr1PID[k]->Fill(pt, mass, event_weight*pair_weight);

                     no_track1_prob = Product(
                        1. - ptrackDecreasedAcceptance.pc2_weight[k][i], 
                        1. - ptrackDecreasedAcceptance.pc3_weight[k][i], 
                        1. - ptrackDecreasedAcceptance.tof_weight[k][i], 
                        1. - ptrackDecreasedAcceptance.emc_weight[k][i]);

                     no_id_track1_prob = 1. - ptrackDecreasedAcceptance.tof_id_weight[k][i];
                     
                     no_track2_prob = Product(
                        1. - ntrackDecreasedAcceptance.pc2_weight[k][j], 
                        1. - ntrackDecreasedAcceptance.pc3_weight[k][j], 
                        1. - ntrackDecreasedAcceptance.tof_weight[k][j], 
                        1. - ntrackDecreasedAcceptance.emc_weight[k][j]);

                     no_id_track2_prob = 1. - ntrackDecreasedAcceptance.tof_id_weight[k][j];
                     
                     pair_weight = ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*(1. -
                        (1. - (1. - no_id_track1_prob)*(1. - no_track2_prob))*
                        (1. - (1. - no_track1_prob)*(1. - no_id_track2_prob)));

                     invMDistr1PIDDecreasedAcceptance[k]->Fill(pt, mass, event_weight*pair_weight);

                     no_track1_prob = Product(
                        1. - ptrackIncreasedAcceptance.pc2_weight[k][i], 
                        1. - ptrackIncreasedAcceptance.pc3_weight[k][i], 
                        1. - ptrackIncreasedAcceptance.tof_weight[k][i], 
                        1. - ptrackIncreasedAcceptance.emc_weight[k][i]);

                     no_id_track1_prob = 1. - ptrackIncreasedAcceptance.tof_id_weight[k][i];
                     
                     no_track2_prob = Product(
                        1. - ntrackIncreasedAcceptance.pc2_weight[k][j], 
                        1. - ntrackIncreasedAcceptance.pc3_weight[k][j], 
                        1. - ntrackIncreasedAcceptance.tof_weight[k][j], 
                        1. - ntrackIncreasedAcceptance.emc_weight[k][j]);
                     
                     no_id_track2_prob = 1. - ntrackIncreasedAcceptance.tof_id_weight[k][j];
                     
                     pair_weight = ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*(1. -
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
                        1. - ptrack.tof_id_weight[k][i], 
                        1. - ptrack.emc_id_weight[k][i]);
                     
                     double no_track2_prob = Product(
                        1. - ntrack.tof_id_weight[k][j], 
                        1. - ntrack.emc_id_weight[k][j]);
                     
                     double pair_weight = ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*
                        (1. - no_track1_prob)*(1. - no_track2_prob);

                     invMDistr2PID[k]->Fill(pt, mass, event_weight*pair_weight);

                     no_track1_prob = Product(
                        1. - ptrackDecreasedAcceptance.tof_id_weight[k][i], 
                        1. - ptrackDecreasedAcceptance.emc_id_weight[k][i]);
                     
                     no_track2_prob = Product(
                        1. - ntrackDecreasedAcceptance.tof_id_weight[k][j], 
                        1. - ntrackDecreasedAcceptance.emc_id_weight[k][j]);
                     
                     pair_weight = ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*
                        (1. - no_track1_prob)*(1. - no_track2_prob);

                     invMDistr2PIDDecreasedAcceptance[k]->Fill(pt, mass, event_weight*pair_weight);

                     no_track1_prob = Product(
                        1. - ptrackIncreasedAcceptance.tof_id_weight[k][i], 
                        1. - ptrackIncreasedAcceptance.emc_id_weight[k][i]);
                     
                     no_track2_prob = Product(
                        1. - ntrackIncreasedAcceptance.tof_id_weight[k][j], 
                        1. - ntrackIncreasedAcceptance.emc_id_weight[k][j]);
                     
                     pair_weight = ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*
                        (1. - no_track1_prob)*(1. - no_track2_prob);

                     invMDistr2PIDIncreasedAcceptance[k]->Fill(pt, mass, event_weight*pair_weight);

                     no_track1_prob = Product(
                        1. - ptrack.tof_id_weight[k][i], 
                        1. - ptrackDecreasedM2Eff.emc_id_weight[k][i]);
                     
                     no_track2_prob = Product(
                        1. - ntrack.tof_id_weight[k][j], 
                        1. - ntrackDecreasedM2Eff.emc_id_weight[k][j]);
                     
                     pair_weight = ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*
                        (1. - no_track1_prob)*(1. - no_track2_prob);

                     invMDistr2PIDDecreasedM2Eff[k]->Fill(pt, mass, event_weight*pair_weight);

                     no_track1_prob = Product(
                        1. - ptrack.tof_id_weight[k][i], 
                        1. - ptrackIncreasedM2Eff.emc_id_weight[k][i]);
                     
                     no_track2_prob = Product(
                        1. - ntrack.tof_id_weight[k][j], 
                        1. - ntrackIncreasedM2Eff.emc_id_weight[k][j]);
                     
                     pair_weight = ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*
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
                        ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*
                        ptrack.tof_id_weight[k][i]*ntrack.tof_id_weight[k][j];
                     
                     invMDistrTOF2PID[k]->Fill(pt, mass, event_weight*pair_weight);
                     
                     pair_weight = 
                        ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*
                        ptrackDecreasedAcceptance.tof_id_weight[k][i]*ntrackDecreasedAcceptance.tof_id_weight[k][j];
                     
                     invMDistrTOF2PIDDecreasedAcceptance[k]->Fill(pt, mass, event_weight*pair_weight);

                     pair_weight = 
                        ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*
                        ptrackIncreasedAcceptance.tof_id_weight[k][i]*ntrackIncreasedAcceptance.tof_id_weight[k][j];
                     
                     invMDistrTOF2PIDIncreasedAcceptance[k]->Fill(pt, mass, event_weight*pair_weight);
                  }
               }

               if (IsEMC2PID(&Id1, &Id2))
               {
                  for (int k = 0; k < CentralityContainer.size; k++)
                  {
                     double pair_weight = 
                        ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*
                        ptrack.emc_id_weight[k][i]*ntrack.emc_id_weight[k][j];

                     invMDistrEMC2PID[k]->Fill(pt, mass, event_weight*pair_weight);

                     pair_weight = 
                        ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*
                        ptrackDecreasedAcceptance.emc_id_weight[k][i]*ntrackDecreasedAcceptance.emc_id_weight[k][j];

                     invMDistrEMC2PIDDecreasedAcceptance[k]->Fill(pt, mass, event_weight*pair_weight);

                     pair_weight = 
                        ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*
                        ptrackIncreasedAcceptance.emc_id_weight[k][i]*ntrackIncreasedAcceptance.emc_id_weight[k][j];

                     invMDistrEMC2PIDIncreasedAcceptance[k]->Fill(pt, mass, event_weight*pair_weight);
                     
                     pair_weight = 
                        ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*
                        ptrackDecreasedM2Eff.emc_id_weight[k][i]*ntrackDecreasedM2Eff.emc_id_weight[k][j];

                     invMDistrEMC2PIDDecreasedM2Eff[k]->Fill(pt, mass, event_weight*pair_weight);

                     pair_weight = 
                        ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*
                        ptrackIncreasedM2Eff.emc_id_weight[k][i]*ntrackIncreasedM2Eff.emc_id_weight[k][j];

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
      ThrContainerStruct ThrContainer;

      const double min_inv_mass = 0.;

      for (long unsigned int j = 0; j < CentralityContainer.nameWithoutPercent.size(); j++)
      {
         std::string dir = CentralityContainer.nameWithoutPercent[j];
         
         ThrContainer.invMDistrNoPID[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "noPID_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "noPID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));
         
         ThrContainer.invMDistrNoPIDDecreasedAcceptance[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "noPID_A-_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
                "noPID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
                Par.invMNBins, min_inv_mass, 4., dir));
         
         ThrContainer.invMDistrNoPIDIncreasedAcceptance[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "noPID_A+_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "noPID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));

         ThrContainer.invMDistr1PID[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "1PID_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "1PID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));
         
         ThrContainer.invMDistr1PIDDecreasedAcceptance[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "1PID_A-_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "1PID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));
         
         ThrContainer.invMDistr1PIDIncreasedAcceptance[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "1PID_A+_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "1PID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));

         ThrContainer.invMDistr2PID[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "2PID_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "2PID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));
         
         ThrContainer.invMDistr2PIDDecreasedAcceptance[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "2PID_A-_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "2PID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));
         
         ThrContainer.invMDistr2PIDIncreasedAcceptance[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "2PID_A+_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "2PID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));
         
         ThrContainer.invMDistr2PIDDecreasedM2Eff[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "2PID_m2Eff-_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "2PID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));
         
         ThrContainer.invMDistr2PIDIncreasedM2Eff[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "2PID_m2Eff+_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "2PID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));
         
         ThrContainer.invMDistrTOF2PID[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "TOF2PID_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "TOF2PID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));
         
         ThrContainer.invMDistrTOF2PIDDecreasedAcceptance[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "TOF2PID_A-_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "TOF2PID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));
         
         ThrContainer.invMDistrTOF2PIDIncreasedAcceptance[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "TOF2PID_A+_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "TOF2PID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));

         ThrContainer.invMDistrEMC2PID[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "EMC2PID_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "EMC2PID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));
         
         ThrContainer.invMDistrEMC2PIDDecreasedAcceptance[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "EMC2PID_A-_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "EMC2PID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));
         
         ThrContainer.invMDistrEMC2PIDIncreasedAcceptance[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "EMC2PID_A+_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "EMC2PID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));
         
         ThrContainer.invMDistrEMC2PIDDecreasedM2Eff[j] = 
            std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>((
               "EMC2PID_m2Eff-_" + CentralityContainer.nameWithoutPercent[j]).c_str(), 
               "EMC2PID", Par.ptNBins, Par.ptMinPair, Par.ptMaxPair, 
               Par.invMNBins, min_inv_mass, 4., dir));
         
         ThrContainer.invMDistrEMC2PIDIncreasedM2Eff[j] = 
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
               AnalyzeConfiguration(&ThrContainer, Par.daughter1Queue[i], Par.daughter2Queue[i], 
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
