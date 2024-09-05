// $SOURCE$
//------------------------------------------------------------------------------------------------
//                     AnalyzeTrackPair function realisation
//------------------------------------------------------------------------------------------------
// AnalyzeTrackPair
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic macro used for evaluation of registration of pair of tracks for different methods
 * from simulation output of event-like TTrees to processed histograms 
 * for further track pair registering correction evaluation
 **/
//------------------------------------------------------------------------------------------------

#ifndef ANALYZE_TRACK_PAIR_CPP
#define ANALYZE_TRACK_PAIR_CPP

//Analyze specific track pair configuration
void AnalyzeConfiguration(ThrContainerStruct *ThrContainer, const std::string& daughter1, 
                          const std::string& daughter2, const std::string magf, 
                          const std::string& auxName, const double ptDeviation, 
                          const int procNum)
{
   Box box = Box("Parameters of run " + std::to_string(procNum) + " out of " + std::to_string(
      Par.daughter1_queue.size()*Par.magf_queue.size()*
      Par.ptDeviation_queue.size()*Par.auxName_queue.size()));

   
   const std::string decay_channel = 
      ParticleProperties.sname[ParticleProperties.iter_map[daughter1]] + 
      ParticleProperties.sname[ParticleProperties.iter_map[daughter2]];

   const std::string input_file_name = "../data/" + Par.run_name + "/Resonances/" + 
      Par.orig_part.name + "_" + decay_channel + magf + auxName + ".root";
   
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
   if (Par.do_use_weight_func)
   {
      event_norm = origPtDistribution->Integral(
         origPtDistribution->GetXaxis()->FindBin(Par.pair_ptmin),
         origPtDistribution->GetXaxis()->FindBin(Par.pair_ptmax));

      //this normalization is needed to merge 2 files with flapt pt distributin with different ranges
      event_norm *= (Par.pair_ptmax - Par.pair_ptmin)/(up_pt_bound - low_pt_bound);
   }

   box.AddEntry("Run name", Par.run_name);
   box.AddEntry("Orig particle", Par.orig_part.name);
   box.AddEntry("Magnetic field configuration", magf);
   box.AddEntry("Use weight function", Par.do_use_weight_func);
   
   box.AddEntry("Daughter particle 1", daughter1);
   box.AddEntry("Daughter particle 2", daughter2);
   
   box.AddEntry("Minimum p_T, GeV", Par.ptmin);
   box.AddEntry("Maximum p_T, GeV", Par.ptmax);

   box.AddEntry("Input file orig p_T range, GeV", DtoStr(low_pt_bound, 1) + "-" + DtoStr(up_pt_bound, 1));
   
   box.AddEntry("Number of events to be analyzed, 1e6", nevents/1e6, 3);

   box.AddEntry("pT deviation", DtoStr(ptDeviation, 3));
   
   box.AddEntry("Number of threads", Par.nthreads);
   
   box.Print();

   TF1 weight_func = TF1("weight_func", 
      "[0]*x^([5]-1)*(1 + ([1] - 1)*(sqrt([4] + x^2)^([5]) - [2])/[3])^(-1./([1]-1.))");
   weight_func.SetParameters(
      ReadFileIntoArray("../input/Spectra/" + Par.system + "/" + Par.orig_part.name + ".txt", 6));
   
   ROOT::EnableImplicitMT(Par.nthreads);
   ROOT::TTreeProcessorMT tp(input_file_name.c_str());
   
   double ncalls = 0.;

   //acceptance systematic uncertainties of detectors
   double dceast0_sys, dceast1_sys, dcwest0_sys, dcwest1_sys,
      pc1e_sys, pc1w_sys, pc2_sys, pc3e_sys, pc3w_sys,
      tofe0_sys, tofe1_sys, tofw0_sys, tofw1_sys;
   std::array<double, 4> emcale_pos_sys, emcale_neg_sys, emcalw_pos_sys, emcalw_neg_sys;

   ReadFile("../input/Systematics/" + Par.run_name + "/acceptance.txt",
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
      PartArray ptrack_acc_decreased, ptrack_acc_increased, ntrack_acc_decreased, ntrack_acc_increased;
      PartArray ptrack_m2_eff_decreased, ptrack_m2_eff_increased, ntrack_m2_eff_decreased, ntrack_m2_eff_increased;

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
            ReadFileIntoArray("../input/Embedding/" + Par.run_name + "/" + 
         ParticleProperties.abs_name[ptrack.iter] + "_dc_pc1.txt", Par.CType.size));
      ptrack.pc2_emb = std::unique_ptr<double>(
            ReadFileIntoArray("../input/Embedding/" + Par.run_name + "/" + 
         ParticleProperties.abs_name[ptrack.iter] + "_pc2.txt", Par.CType.size));
      ptrack.pc3_emb = std::unique_ptr<double>(
            ReadFileIntoArray("../input/Embedding/" + Par.run_name + "/" + 
         ParticleProperties.abs_name[ptrack.iter] + "_pc3.txt", Par.CType.size));
      ptrack.tofe_emb = std::unique_ptr<double>(
            ReadFileIntoArray("../input/Embedding/" + Par.run_name + "/" + 
         ParticleProperties.abs_name[ptrack.iter] + "_tofe.txt", Par.CType.size));
      ptrack.tofw_emb = std::unique_ptr<double>(
            ReadFileIntoArray("../input/Embedding/" + Par.run_name + "/" + 
         ParticleProperties.abs_name[ptrack.iter] + "_tofw.txt", Par.CType.size));

      ntrack.dc_pc1_emb = std::unique_ptr<double>(
            ReadFileIntoArray("../input/Embedding/" + Par.run_name + "/" + 
         ParticleProperties.abs_name[ntrack.iter] + "_dc_pc1.txt", Par.CType.size));
      ntrack.pc2_emb = std::unique_ptr<double>(
            ReadFileIntoArray("../input/Embedding/" + Par.run_name + "/" + 
         ParticleProperties.abs_name[ntrack.iter] + "_pc2.txt", Par.CType.size));
      ntrack.pc3_emb = std::unique_ptr<double>(
            ReadFileIntoArray("../input/Embedding/" + Par.run_name + "/" + 
         ParticleProperties.abs_name[ntrack.iter] + "_pc3.txt", Par.CType.size));
      ntrack.tofe_emb = std::unique_ptr<double>(
            ReadFileIntoArray("../input/Embedding/" + Par.run_name + "/" + 
         ParticleProperties.abs_name[ntrack.iter] + "_tofe.txt", Par.CType.size));
      ntrack.tofw_emb = std::unique_ptr<double>(
            ReadFileIntoArray("../input/Embedding/" + Par.run_name + "/" + 
         ParticleProperties.abs_name[ntrack.iter] + "_tofw.txt", Par.CType.size));
      
      for (int i = 0; i < 4; i++)
      {
         ptrack.emcale_emb[i] = std::unique_ptr<double>(
               ReadFileIntoArray("../input/Embedding/" + Par.run_name + "/" + 
            ParticleProperties.abs_name[ptrack.iter] + "_emcale" + std::to_string(i) + ".txt", Par.CType.size));
         ptrack.emcalw_emb[i] = std::unique_ptr<double>(
            ReadFileIntoArray("../input/Embedding/" + Par.run_name + "/" + 
            ParticleProperties.abs_name[ptrack.iter] + "_emcalw" + std::to_string(i) + ".txt", Par.CType.size));
         
         ntrack.emcale_emb[i] = std::unique_ptr<double>(
            ReadFileIntoArray("../input/Embedding/" + Par.run_name + "/" + 
            ParticleProperties.abs_name[ntrack.iter] + "_emcale" + std::to_string(i) + ".txt", Par.CType.size));
         ntrack.emcalw_emb[i] = std::unique_ptr<double>(
            ReadFileIntoArray("../input/Embedding/" + Par.run_name + "/" + 
            ParticleProperties.abs_name[ntrack.iter] + "_emcalw" + std::to_string(i) + ".txt", Par.CType.size));
      }

      ptrack.emcale_m2_eff = std::unique_ptr<double>(
            ReadFileIntoArray("../input/M2Eff/" + Par.run_name + "/" +
            ParticleProperties.name[ptrack.iter] + "_EMCale" + magf + ".txt", 2));
      ptrack.emcalw_m2_eff = std::unique_ptr<double>(
            ReadFileIntoArray("../input/M2Eff/" + Par.run_name + "/" +
            ParticleProperties.name[ptrack.iter] + "_EMCalw" + magf + ".txt", 4));
      
      ntrack.emcale_m2_eff = std::unique_ptr<double>(
            ReadFileIntoArray("../input/M2Eff/" + Par.run_name + "/" +
            ParticleProperties.name[ntrack.iter] + "_EMCale" + magf + ".txt", 2));
      ntrack.emcalw_m2_eff = std::unique_ptr<double>(
            ReadFileIntoArray("../input/M2Eff/" + Par.run_name + "/" +
            ParticleProperties.name[ntrack.iter] + "_EMCalw" + magf + ".txt", 4));

      ptrack.emcale_m2_eff_sys = std::unique_ptr<double>(
            ReadFileIntoArray("../input/Systematics/" + Par.run_name + "/m2eff_" +
            ParticleProperties.name[ptrack.iter] + "_EMCale" + magf + ".txt", 2));
      ptrack.emcalw_m2_eff_sys = std::unique_ptr<double>(
            ReadFileIntoArray("../input/Systematics/" + Par.run_name + "/m2eff_" +
            ParticleProperties.name[ptrack.iter] + "_EMCalw" + magf + ".txt", 4));
      
      ntrack.emcale_m2_eff_sys = std::unique_ptr<double>(
            ReadFileIntoArray("../input/Systematics/" + Par.run_name + "/m2eff_" +
            ParticleProperties.name[ntrack.iter] + "_EMCale" + magf + ".txt", 2));
      ntrack.emcalw_m2_eff_sys = std::unique_ptr<double>(
            ReadFileIntoArray("../input/Systematics/" + Par.run_name + "/m2eff_" +
            ParticleProperties.name[ntrack.iter] + "_EMCalw" + magf + ".txt", 4));

      std::shared_ptr<TH1F> orig = ThrContainer->origPtDistr.Get();
      
      std::shared_ptr<TH2F> orig_pt_vs_pt = ThrContainer->orig_pt_vs_pt.Get();
      
      std::array<std::shared_ptr<TH2F>, Par.CType.size> InvM_noPID;
      std::array<std::shared_ptr<TH2F>, Par.CType.size> InvM_noPID_acc_decreased;
      std::array<std::shared_ptr<TH2F>, Par.CType.size> InvM_noPID_acc_increased;

      std::array<std::shared_ptr<TH2F>, Par.CType.size> InvM_1PID;
      std::array<std::shared_ptr<TH2F>, Par.CType.size> InvM_1PID_acc_decreased;
      std::array<std::shared_ptr<TH2F>, Par.CType.size> InvM_1PID_acc_increased;

      std::array<std::shared_ptr<TH2F>, Par.CType.size> InvM_2PID;
      std::array<std::shared_ptr<TH2F>, Par.CType.size> InvM_2PID_acc_decreased;
      std::array<std::shared_ptr<TH2F>, Par.CType.size> InvM_2PID_acc_increased;
      std::array<std::shared_ptr<TH2F>, Par.CType.size> InvM_2PID_m2_eff_decreased;
      std::array<std::shared_ptr<TH2F>, Par.CType.size> InvM_2PID_m2_eff_increased;

      std::array<std::shared_ptr<TH2F>, Par.CType.size> InvM_TOF2PID;
      std::array<std::shared_ptr<TH2F>, Par.CType.size> InvM_TOF2PID_acc_decreased;
      std::array<std::shared_ptr<TH2F>, Par.CType.size> InvM_TOF2PID_acc_increased;

      std::array<std::shared_ptr<TH2F>, Par.CType.size> InvM_EMC2PID;
      std::array<std::shared_ptr<TH2F>, Par.CType.size> InvM_EMC2PID_acc_decreased;
      std::array<std::shared_ptr<TH2F>, Par.CType.size> InvM_EMC2PID_acc_increased;
      std::array<std::shared_ptr<TH2F>, Par.CType.size> InvM_EMC2PID_m2_eff_decreased;
      std::array<std::shared_ptr<TH2F>, Par.CType.size> InvM_EMC2PID_m2_eff_increased;

      for (int i = 0; i < Par.CType.size; i++)
      {
         InvM_noPID[i] = ThrContainer->InvM_noPID[i]->Get();
         InvM_noPID_acc_decreased[i] = ThrContainer->InvM_noPID_acc_decreased[i]->Get();
         InvM_noPID_acc_increased[i] = ThrContainer->InvM_noPID_acc_increased[i]->Get();

         InvM_1PID[i] = ThrContainer->InvM_1PID[i]->Get();
         InvM_1PID_acc_decreased[i] = ThrContainer->InvM_1PID_acc_decreased[i]->Get();
         InvM_1PID_acc_increased[i] = ThrContainer->InvM_1PID_acc_increased[i]->Get();

         InvM_2PID[i] = ThrContainer->InvM_2PID[i]->Get();
         InvM_2PID_acc_decreased[i] = ThrContainer->InvM_2PID_acc_decreased[i]->Get();
         InvM_2PID_acc_increased[i] = ThrContainer->InvM_2PID_acc_increased[i]->Get();
         InvM_2PID_m2_eff_decreased[i] = ThrContainer->InvM_2PID_m2_eff_decreased[i]->Get();
         InvM_2PID_m2_eff_increased[i] = ThrContainer->InvM_2PID_m2_eff_increased[i]->Get();

         InvM_TOF2PID[i] = ThrContainer->InvM_TOF2PID[i]->Get();
         InvM_TOF2PID_acc_decreased[i] = ThrContainer->InvM_TOF2PID_acc_decreased[i]->Get();
         InvM_TOF2PID_acc_increased[i] = ThrContainer->InvM_TOF2PID_acc_increased[i]->Get();

         InvM_EMC2PID[i] = ThrContainer->InvM_EMC2PID[i]->Get();
         InvM_EMC2PID_acc_decreased[i] = ThrContainer->InvM_EMC2PID_acc_decreased[i]->Get();
         InvM_EMC2PID_acc_increased[i] = ThrContainer->InvM_EMC2PID_acc_increased[i]->Get();
         InvM_EMC2PID_m2_eff_decreased[i] = ThrContainer->InvM_EMC2PID_m2_eff_decreased[i]->Get();
         InvM_EMC2PID_m2_eff_increased[i] = ThrContainer->InvM_EMC2PID_m2_eff_increased[i]->Get();
      }
      
      EffTreeReader T(reader);
      
      while (reader.Next())
      {   
         ncalls += 1.;
         const double orig_pt = sqrt(pow(T.mom_orig(0), 2) + pow(T.mom_orig(1), 2))*ptDeviation;
         
         double event_weight = 1.;
         
         if (Par.do_use_weight_func) event_weight = weight_func.Eval(orig_pt)/event_norm*4e11;
         
         orig->Fill(orig_pt, event_weight);
         
         if (T.nch() > 49 || T.nch() <= 0) continue;
         
         const double bbcz = T.bbcz();
         if (fabs(bbcz) > 30) continue;

         ptrack.size = 0;
         ntrack.size = 0;

         for(int i = 0; i < T.nch(); i++)
         {   
            const double the0 = T.the0(i);
            const double pt = (T.mom(i))*sin(the0)*ptDeviation;

            if (pt < Par.ptmin || pt > Par.ptmax) continue;
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
            ptrack_acc_decreased.ResetTrack(ptrack.size);
            ntrack_acc_decreased.ResetTrack(ntrack.size);
            ptrack_acc_increased.ResetTrack(ptrack.size);
            ntrack_acc_increased.ResetTrack(ntrack.size);
            ptrack_m2_eff_decreased.ResetTrack(ptrack.size);
            ntrack_m2_eff_decreased.ResetTrack(ntrack.size);
            ptrack_m2_eff_increased.ResetTrack(ptrack.size);
            ntrack_m2_eff_increased.ResetTrack(ntrack.size);

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
                     
                     for (int j = 0; j < Par.CType.size; j++)
                     {
                        ptrack.tof_weight[j][ptrack.size] = 
                           ptrack.tofe_emb.get()[j]/ptrack.dc_pc1_emb.get()[j];
                        ptrack_acc_decreased.tof_weight[j][ptrack.size] = 
                           ptrack.tof_weight[j][ptrack.size]*(1. - acc_var);
                        ptrack_acc_increased.tof_weight[j][ptrack.size] = 
                           Minimum(1., ptrack.tof_weight[j][ptrack.size]*(1.+acc_var));

                        if (id == ptrack.origId) 
                        {
                           ptrack.tof_id_weight[j][ptrack.size] = 
                              ptrack.tof_weight[j][ptrack.size];
                           ptrack_acc_decreased.tof_id_weight[j][ptrack.size] = 
                              ptrack.tof_id_weight[j][ptrack.size]*(1. - acc_var);
                           ptrack_acc_increased.tof_id_weight[j][ptrack.size] = 
                              Minimum(1., ptrack.tof_id_weight[j][ptrack.size]*(1.+acc_var));
                        }   
                     }
                  }
                  else
                  {
                     ntrack.tof_id[ntrack.size] = id;
                     
                     for (int j = 0; j < Par.CType.size; j++)
                     {
                        ntrack.tof_weight[j][ntrack.size] = 
                           ntrack.tofe_emb.get()[j]/ntrack.dc_pc1_emb.get()[j];
                        ntrack_acc_decreased.tof_weight[j][ntrack.size] = 
                           ntrack.tof_weight[j][ntrack.size]*(1. - acc_var);
                        ntrack_acc_increased.tof_weight[j][ntrack.size] = 
                           Minimum(1., ntrack.tof_weight[j][ntrack.size]*(1.+acc_var));
                        
                        if (id == ntrack.origId)
                        {
                           ntrack.tof_id_weight[j][ntrack.size] = 
                              ntrack.tof_weight[j][ntrack.size];
                           ntrack_acc_decreased.tof_id_weight[j][ntrack.size] = 
                              ntrack.tof_id_weight[j][ntrack.size]*(1. - acc_var);
                           ntrack_acc_increased.tof_id_weight[j][ntrack.size] = Minimum(1.,
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
                     
                     for (unsigned int j = 0; j < Par.CType.size; j++)
                     {
                        ptrack.tof_weight[j][ptrack.size] = 
                           Par.tofw_correction*ptrack.tofw_emb.get()[j]/ptrack.dc_pc1_emb.get()[j];
                        ptrack_acc_decreased.tof_weight[j][ptrack.size] = 
                           ptrack.tof_weight[j][ptrack.size]*(1. - acc_var);
                        ptrack_acc_increased.tof_weight[j][ptrack.size] = Minimum(1.,
                           ptrack.tof_weight[j][ptrack.size]*(1. + acc_var));

                        if (id == ptrack.origId) 
                        {
                           ptrack.tof_id_weight[j][ptrack.size] = 
                               ptrack.tof_weight[j][ptrack.size];
                           ptrack_acc_decreased.tof_id_weight[j][ptrack.size] = 
                               ptrack.tof_weight[j][ptrack.size]*(1. - acc_var);
                           ptrack_acc_increased.tof_id_weight[j][ptrack.size] = Minimum(1.,
                              ptrack.tof_weight[j][ptrack.size]*(1.+acc_var));
                        }
                     }
                  }
                  else
                  {
                     ntrack.tof_id[ntrack.size] = id;
                     
                     for (unsigned int j = 0; j < Par.CType.size; j++)
                     {
                        ntrack.tof_weight[j][ntrack.size] = 
                           Par.tofw_correction*ntrack.tofw_emb.get()[j]/ntrack.dc_pc1_emb.get()[j];
                        ntrack_acc_decreased.tof_weight[j][ntrack.size] = 
                           ntrack.tof_weight[j][ntrack.size]*(1. - acc_var);
                        ntrack_acc_increased.tof_weight[j][ntrack.size] = Minimum(1.,
                           ntrack.tof_weight[j][ntrack.size]*(1.+acc_var));

                        if (id == ntrack.origId) 
                        {
                           ntrack.tof_id_weight[j][ntrack.size] = 
                              ntrack.tof_weight[j][ntrack.size];
                           ntrack_acc_decreased.tof_id_weight[j][ntrack.size] = 
                              ntrack.tof_id_weight[j][ntrack.size]*(1. - acc_var);
                           ntrack_acc_increased.tof_id_weight[j][ntrack.size] = Minimum(1.,
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
                        for (unsigned int j = 0; j < Par.CType.size; j++)
                        {
                           ptrack.emc_weight[j][ptrack.size] = 
                              ptrack.emcale_emb[T.sect(i)].get()[j]/ptrack.dc_pc1_emb.get()[j];
                        }
         
                        if (T.particle_id(i) == ptrack.geant_id && T.sect(i) > 1)
                        {
                           ptrack.emc_id[ptrack.size] = ptrack.origId;
                           for (unsigned int j = 0; j < Par.CType.size; j++)
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
                        for (unsigned int j = 0; j < Par.CType.size; j++)
                        {
                           ntrack.emc_weight[j][ntrack.size] = 
                              ntrack.emcale_emb[T.sect(i)].get()[j]/ntrack.dc_pc1_emb.get()[j];
                        }
                           
                        if (T.particle_id(i) == ntrack.geant_id && T.sect(i) > 1)
                        {
                           ntrack.emc_id[ntrack.size] = ntrack.origId;
                           for (unsigned int j = 0; j < Par.CType.size; j++)
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
                        for (unsigned int j = 0; j < Par.CType.size; j++)
                        {
                           ptrack.emc_weight[j][ptrack.size] = 
                              ptrack.emcalw_emb[T.sect(i)].get()[j]/ptrack.dc_pc1_emb.get()[j];
                        }
                        
                        if (T.particle_id(i) == ptrack.geant_id)
                        {
                           ptrack.emc_id[ptrack.size] = ptrack.origId;
                           for (unsigned int j = 0; j < Par.CType.size; j++)
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
                        for (unsigned int j = 0; j < Par.CType.size; j++)
                        {
                           ntrack.emc_weight[j][ntrack.size] = 
                              ntrack.emcalw_emb[T.sect(i)].get()[j]/ntrack.dc_pc1_emb.get()[j];
                        }
                        
                        if (T.particle_id(i) == ntrack.geant_id)
                        {
                           ntrack.emc_id[ntrack.size] = ntrack.origId;
                           for (unsigned int j = 0; j < Par.CType.size; j++)
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
                     for (unsigned int j = 0; j < Par.CType.size; j++)
                     {
                        ptrack_acc_decreased.emc_weight[j][ptrack.size] = 
                           ptrack.emc_weight[j][ptrack.size]*(1. - acc_var);
                        ptrack_acc_increased.emc_weight[j][ptrack.size] = Minimum(1.,
                           ptrack.emc_weight[j][ptrack.size]*(1. + acc_var));

                        ptrack_acc_decreased.emc_id_weight[j][ptrack.size] = 
                           ptrack.emc_id_weight[j][ptrack.size]*(1. - acc_var);
                        ptrack_acc_increased.emc_id_weight[j][ptrack.size] = Minimum(1.,
                           ptrack.emc_id_weight[j][ptrack.size]*(1. + acc_var));

                        ptrack_m2_eff_decreased.emc_id_weight[j][ptrack.size] = 
                           ptrack.emc_id_weight[j][ptrack.size]*(1. - m2_eff_var);
                        ptrack_m2_eff_increased.emc_id_weight[j][ptrack.size] = Minimum(1.,
                           ptrack.emc_id_weight[j][ptrack.size]*(1. + m2_eff_var));
                     }
                  }
                  else
                  {
                     for (unsigned int j = 0; j < Par.CType.size; j++)
                     {
                        ntrack_acc_decreased.emc_weight[j][ntrack.size] = 
                           ntrack.emc_weight[j][ntrack.size]*(1. - acc_var);
                        ntrack_acc_increased.emc_weight[j][ntrack.size] = Minimum(1.,
                           ntrack.emc_weight[j][ntrack.size]*(1. + acc_var));

                        ntrack_acc_decreased.emc_id_weight[j][ntrack.size] = 
                           ntrack.emc_id_weight[j][ntrack.size]*(1. - acc_var);
                        ntrack_acc_increased.emc_id_weight[j][ntrack.size] = Minimum(1.,
                           ntrack.emc_id_weight[j][ntrack.size]*(1. + acc_var));

                        ntrack_m2_eff_decreased.emc_id_weight[j][ntrack.size] = 
                           ntrack.emc_id_weight[j][ntrack.size]*(1. - m2_eff_var);
                        ntrack_m2_eff_increased.emc_id_weight[j][ntrack.size] = Minimum(1.,
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
   
                        for (int j = 0; j < Par.CType.size; j++)
                        {
                           ptrack.pc2_weight[j][ptrack.size] = 
                              ptrack.pc2_emb.get()[j]/ptrack.dc_pc1_emb.get()[j];
                           
                           ptrack_acc_decreased.pc2_weight[j][ptrack.size] = 
                              ptrack.pc2_weight[j][ptrack.size]*(1. - acc_var);
                           ptrack_acc_increased.pc2_weight[j][ptrack.size] = Minimum(1.,
                              ptrack.pc2_weight[j][ptrack.size]*(1. + acc_var));
                        }
                     }
                     else
                     {
                        ntrack.pc2_id[ntrack.size] = PartId.noPID;
                        
                        for (int j = 0; j < Par.CType.size; j++)
                        {
                           ntrack.pc2_weight[j][ntrack.size] = 
                              ntrack.pc2_emb.get()[j]/ntrack.dc_pc1_emb.get()[j];

                           ntrack_acc_decreased.pc2_weight[j][ntrack.size] = 
                              ntrack.pc2_weight[j][ntrack.size]*(1. - acc_var);
                           ntrack_acc_increased.pc2_weight[j][ntrack.size] = Minimum(1.,
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
                        
                        for (int j = 0; j < Par.CType.size; j++)
                        {
                           ptrack.pc3_weight[j][ptrack.size] = 
                              ptrack.pc3_emb.get()[j]/ptrack.dc_pc1_emb.get()[j];

                           ptrack_acc_decreased.pc3_weight[j][ptrack.size] = 
                              ptrack.pc3_weight[j][ptrack.size]*(1. - acc_var); 
                           ptrack_acc_increased.pc3_weight[j][ptrack.size] = Minimum(1.,
                              ptrack.pc3_weight[j][ptrack.size]*(1. + acc_var)); 
                        }
                     }
                     else
                     {
                        ntrack.pc3_id[ntrack.size] = PartId.noPID;
                        
                        for (int j = 0; j < Par.CType.size; j++)
                        {
                           ntrack.pc3_weight[j][ntrack.size] = 
                              ntrack.pc3_emb.get()[j]/ntrack.dc_pc1_emb.get()[j];
                           
                           ntrack_acc_decreased.pc3_weight[j][ntrack.size] = 
                              ntrack.pc3_weight[j][ntrack.size]*(1. - acc_var); 
                           ntrack_acc_increased.pc3_weight[j][ntrack.size] = Minimum(1.,
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
                  for (int k = 0; k < Par.CType.size; k++)
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
      
                     InvM_noPID[k]->Fill(pt, mass, event_weight*pair_weight);

                     orig_pt_vs_pt->Fill(pt, orig_pt, event_weight*pair_weight);

                     no_track1_prob = Product(
                        1. - ptrack_acc_decreased.pc2_weight[k][i], 
                        1. - ptrack_acc_decreased.pc3_weight[k][i], 
                        1. - ptrack_acc_decreased.tof_weight[k][i], 
                        1. - ptrack_acc_decreased.emc_weight[k][i]);
                     
                     no_track2_prob = Product(
                        1. - ntrack_acc_decreased.pc2_weight[k][j], 
                        1. - ntrack_acc_decreased.pc3_weight[k][j], 
                        1. - ntrack_acc_decreased.tof_weight[k][j], 
                        1. - ntrack_acc_decreased.emc_weight[k][j]);

                     pair_weight = ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*
                        (1. - no_track1_prob)*(1. - no_track2_prob);

                     InvM_noPID_acc_decreased[k]->Fill(pt, mass, event_weight*pair_weight);

                     no_track1_prob = Product(
                        1. - ptrack_acc_increased.pc2_weight[k][i], 
                        1. - ptrack_acc_increased.pc3_weight[k][i], 
                        1. - ptrack_acc_increased.tof_weight[k][i], 
                        1. - ptrack_acc_increased.emc_weight[k][i]);
                     
                     no_track2_prob = Product(
                        1. - ntrack_acc_increased.pc2_weight[k][j], 
                        1. - ntrack_acc_increased.pc3_weight[k][j], 
                        1. - ntrack_acc_increased.tof_weight[k][j], 
                        1. - ntrack_acc_increased.emc_weight[k][j]);
                     
                     pair_weight = ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*
                        (1. - no_track1_prob)*(1. - no_track2_prob);

                     InvM_noPID_acc_increased[k]->Fill(pt, mass, event_weight*pair_weight);
                  }
               }
               else continue;

               if (Is1PID(&Id1, &Id2))
               {
                  for (int k = 0; k < Par.CType.size; k++)
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

                     InvM_1PID[k]->Fill(pt, mass, event_weight*pair_weight);

                     no_track1_prob = Product(
                        1. - ptrack_acc_decreased.pc2_weight[k][i], 
                        1. - ptrack_acc_decreased.pc3_weight[k][i], 
                        1. - ptrack_acc_decreased.tof_weight[k][i], 
                        1. - ptrack_acc_decreased.emc_weight[k][i]);

                     no_id_track1_prob = 1. - ptrack_acc_decreased.tof_id_weight[k][i];
                     
                     no_track2_prob = Product(
                        1. - ntrack_acc_decreased.pc2_weight[k][j], 
                        1. - ntrack_acc_decreased.pc3_weight[k][j], 
                        1. - ntrack_acc_decreased.tof_weight[k][j], 
                        1. - ntrack_acc_decreased.emc_weight[k][j]);

                     no_id_track2_prob = 1. - ntrack_acc_decreased.tof_id_weight[k][j];
                     
                     pair_weight = ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*(1. -
                        (1. - (1. - no_id_track1_prob)*(1. - no_track2_prob))*
                        (1. - (1. - no_track1_prob)*(1. - no_id_track2_prob)));

                     InvM_1PID_acc_decreased[k]->Fill(pt, mass, event_weight*pair_weight);

                     no_track1_prob = Product(
                        1. - ptrack_acc_increased.pc2_weight[k][i], 
                        1. - ptrack_acc_increased.pc3_weight[k][i], 
                        1. - ptrack_acc_increased.tof_weight[k][i], 
                        1. - ptrack_acc_increased.emc_weight[k][i]);

                     no_id_track1_prob = 1. - ptrack_acc_increased.tof_id_weight[k][i];
                     
                     no_track2_prob = Product(
                        1. - ntrack_acc_increased.pc2_weight[k][j], 
                        1. - ntrack_acc_increased.pc3_weight[k][j], 
                        1. - ntrack_acc_increased.tof_weight[k][j], 
                        1. - ntrack_acc_increased.emc_weight[k][j]);
                     
                     no_id_track2_prob = 1. - ntrack_acc_increased.tof_id_weight[k][j];
                     
                     pair_weight = ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*(1. -
                        (1. - (1. - no_id_track1_prob)*(1. - no_track2_prob))*
                        (1. - (1. - no_track1_prob)*(1. - no_id_track2_prob)));
                     
                     InvM_1PID_acc_increased[k]->Fill(pt, mass, event_weight*pair_weight);
                  }
               }

               if (Is2PID(&Id1, &Id2))
               {
                  for (int k = 0; k < Par.CType.size; k++)
                  {
                     double no_track1_prob = Product(
                        1. - ptrack.tof_id_weight[k][i], 
                        1. - ptrack.emc_id_weight[k][i]);
                     
                     double no_track2_prob = Product(
                        1. - ntrack.tof_id_weight[k][j], 
                        1. - ntrack.emc_id_weight[k][j]);
                     
                     double pair_weight = ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*
                        (1. - no_track1_prob)*(1. - no_track2_prob);

                     InvM_2PID[k]->Fill(pt, mass, event_weight*pair_weight);

                     no_track1_prob = Product(
                        1. - ptrack_acc_decreased.tof_id_weight[k][i], 
                        1. - ptrack_acc_decreased.emc_id_weight[k][i]);
                     
                     no_track2_prob = Product(
                        1. - ntrack_acc_decreased.tof_id_weight[k][j], 
                        1. - ntrack_acc_decreased.emc_id_weight[k][j]);
                     
                     pair_weight = ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*
                        (1. - no_track1_prob)*(1. - no_track2_prob);

                     InvM_2PID_acc_decreased[k]->Fill(pt, mass, event_weight*pair_weight);

                     no_track1_prob = Product(
                        1. - ptrack_acc_increased.tof_id_weight[k][i], 
                        1. - ptrack_acc_increased.emc_id_weight[k][i]);
                     
                     no_track2_prob = Product(
                        1. - ntrack_acc_increased.tof_id_weight[k][j], 
                        1. - ntrack_acc_increased.emc_id_weight[k][j]);
                     
                     pair_weight = ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*
                        (1. - no_track1_prob)*(1. - no_track2_prob);

                     InvM_2PID_acc_increased[k]->Fill(pt, mass, event_weight*pair_weight);

                     no_track1_prob = Product(
                        1. - ptrack.tof_id_weight[k][i], 
                        1. - ptrack_m2_eff_decreased.emc_id_weight[k][i]);
                     
                     no_track2_prob = Product(
                        1. - ntrack.tof_id_weight[k][j], 
                        1. - ntrack_m2_eff_decreased.emc_id_weight[k][j]);
                     
                     pair_weight = ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*
                        (1. - no_track1_prob)*(1. - no_track2_prob);

                     InvM_2PID_m2_eff_decreased[k]->Fill(pt, mass, event_weight*pair_weight);

                     no_track1_prob = Product(
                        1. - ptrack.tof_id_weight[k][i], 
                        1. - ptrack_m2_eff_increased.emc_id_weight[k][i]);
                     
                     no_track2_prob = Product(
                        1. - ntrack.tof_id_weight[k][j], 
                        1. - ntrack_m2_eff_increased.emc_id_weight[k][j]);
                     
                     pair_weight = ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*
                        (1. - no_track1_prob)*(1. - no_track2_prob);

                     InvM_2PID_m2_eff_increased[k]->Fill(pt, mass, event_weight*pair_weight);
                  }
               }
               else continue;

               if (IsTOF2PID(&Id1, &Id2))
               {
                  for (int k = 0; k < Par.CType.size; k++)
                  {
                     double pair_weight = 
                        ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*
                        ptrack.tof_id_weight[k][i]*ntrack.tof_id_weight[k][j];
                     
                     InvM_TOF2PID[k]->Fill(pt, mass, event_weight*pair_weight);
                     
                     pair_weight = 
                        ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*
                        ptrack_acc_decreased.tof_id_weight[k][i]*ntrack_acc_decreased.tof_id_weight[k][j];
                     
                     InvM_TOF2PID_acc_decreased[k]->Fill(pt, mass, event_weight*pair_weight);

                     pair_weight = 
                        ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*
                        ptrack_acc_increased.tof_id_weight[k][i]*ntrack_acc_increased.tof_id_weight[k][j];
                     
                     InvM_TOF2PID_acc_increased[k]->Fill(pt, mass, event_weight*pair_weight);
                  }
               }

               if (IsEMC2PID(&Id1, &Id2))
               {
                  for (int k = 0; k < Par.CType.size; k++)
                  {
                     double pair_weight = 
                        ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*
                        ptrack.emc_id_weight[k][i]*ntrack.emc_id_weight[k][j];

                     InvM_EMC2PID[k]->Fill(pt, mass, event_weight*pair_weight);

                     pair_weight = 
                        ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*
                        ptrack_acc_decreased.emc_id_weight[k][i]*ntrack_acc_decreased.emc_id_weight[k][j];

                     InvM_EMC2PID_acc_decreased[k]->Fill(pt, mass, event_weight*pair_weight);

                     pair_weight = 
                        ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*
                        ptrack_acc_increased.emc_id_weight[k][i]*ntrack_acc_increased.emc_id_weight[k][j];

                     InvM_EMC2PID_acc_increased[k]->Fill(pt, mass, event_weight*pair_weight);
                     
                     pair_weight = 
                        ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*
                        ptrack_m2_eff_decreased.emc_id_weight[k][i]*ntrack_m2_eff_decreased.emc_id_weight[k][j];

                     InvM_EMC2PID_m2_eff_decreased[k]->Fill(pt, mass, event_weight*pair_weight);

                     pair_weight = 
                        ptrack.dc_pc1_emb.get()[k]*ntrack.dc_pc1_emb.get()[k]*
                        ptrack_m2_eff_increased.emc_id_weight[k][i]*ntrack_m2_eff_increased.emc_id_weight[k][j];

                     InvM_EMC2PID_m2_eff_increased[k]->Fill(pt, mass, event_weight*pair_weight);
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
   if (Par.do_use_weight_func) CheckInputFile("../input/Spectra/" + 
      Par.system + "/" + Par.orig_part.name + ".txt");
   
   for (long unsigned int i = 0; i < Par.daughter1_queue.size(); i++)
   {
      const int daughter1_iter = ParticleProperties.iter_map[Par.daughter1_queue[i]];
      const int daughter2_iter = ParticleProperties.iter_map[Par.daughter2_queue[i]];   
      
      const std::string decay_channel = 
         ParticleProperties.sname[daughter1_iter] + 
         ParticleProperties.sname[daughter2_iter];
      
      for (std::string detector : Par.detectors)
      {
         CheckInputFile("../input/Embedding/" + Par.run_name + "/" + 
            ParticleProperties.abs_name[daughter1_iter] + "_" + detector + ".txt");
         CheckInputFile("../input/Embedding/" + Par.run_name + "/" + 
            ParticleProperties.abs_name[daughter2_iter] + "_" + detector + ".txt");   
      }
      
      for (std::string magf : Par.magf_queue)
      {
         CheckInputFile("../input/M2Eff/" + Par.run_name + "/" +
               ParticleProperties.name[daughter1_iter] + "_EMCale" + magf + ".txt");
         CheckInputFile("../input/M2Eff/" + Par.run_name + "/" +
               ParticleProperties.name[daughter1_iter] + "_EMCalw" + magf + ".txt");
         CheckInputFile("../input/Systematics/" + Par.run_name + "/m2eff_" +
               ParticleProperties.name[daughter1_iter] + "_EMCale" + magf + ".txt");
         CheckInputFile("../input/Systematics/" + Par.run_name + "/m2eff_" +
               ParticleProperties.name[daughter1_iter] + "_EMCalw" + magf + ".txt");

         CheckInputFile("../input/M2Eff/" + Par.run_name + "/" +
               ParticleProperties.name[daughter2_iter] + "_EMCale" + magf + ".txt");
         CheckInputFile("../input/M2Eff/" + Par.run_name + "/" +
               ParticleProperties.name[daughter2_iter] + "_EMCalw" + magf + ".txt");
         CheckInputFile("../input/Systematics/" + Par.run_name + "/m2eff_" +
               ParticleProperties.name[daughter2_iter] + "_EMCale" + magf + ".txt");
         CheckInputFile("../input/Systematics/" + Par.run_name + "/m2eff_" +
               ParticleProperties.name[daughter2_iter] + "_EMCalw" + magf + ".txt");
         for (std::string auxName : Par.auxName_queue)
         {
            const std::string input_file_name = "../data/" + Par.run_name + "/Resonances/" + 
               Par.orig_part.name + "_" + decay_channel + magf + auxName + ".root";
            
            CheckInputFile(input_file_name);
         }
      }
   }

   CheckInputFile("../input/Systematics/" + Par.run_name + "/acceptance.txt");
   
   int procNum = 1;
   for (double ptDeviation : Par.ptDeviation_queue)
   {
      ThrContainerStruct ThrContainer;

      const double min_inv_mass = 0.;

      for (long unsigned int j = 0; j < Par.CType.cname_nop.size(); j++)
      {
         std::string dir = Par.CType.cname_nop[j];
         
         ThrContainer.InvM_noPID[j] = std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>
            (("noPID_" + Par.CType.cname_nop[j]).c_str(), 
             "noPID", Par.pt_nbins, Par.pair_ptmin, Par.pair_ptmax, Par.invm_nbins, min_inv_mass, 4., dir));
         ThrContainer.InvM_noPID_acc_decreased[j] = std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>
            (("noPID_acc_decreased_" + Par.CType.cname_nop[j]).c_str(), 
             "noPID", Par.pt_nbins, Par.pair_ptmin, Par.pair_ptmax, Par.invm_nbins, min_inv_mass, 4., dir));
         ThrContainer.InvM_noPID_acc_increased[j] = std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>
            (("noPID_acc_increased_" + Par.CType.cname_nop[j]).c_str(), 
             "noPID", Par.pt_nbins, Par.pair_ptmin, Par.pair_ptmax, Par.invm_nbins, min_inv_mass, 4., dir));

         ThrContainer.InvM_1PID[j] = std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>
            (("1PID_" + Par.CType.cname_nop[j]).c_str(), 
             "1PID", Par.pt_nbins, Par.pair_ptmin, Par.pair_ptmax, Par.invm_nbins, min_inv_mass, 4., dir));
         ThrContainer.InvM_1PID_acc_decreased[j] = std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>
            (("1PID_acc_decreased_" + Par.CType.cname_nop[j]).c_str(), 
             "1PID", Par.pt_nbins, Par.pair_ptmin, Par.pair_ptmax, Par.invm_nbins, min_inv_mass, 4., dir));
         ThrContainer.InvM_1PID_acc_increased[j] = std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>
            (("1PID_acc_increased_" + Par.CType.cname_nop[j]).c_str(), 
             "1PID", Par.pt_nbins, Par.pair_ptmin, Par.pair_ptmax, Par.invm_nbins, min_inv_mass, 4., dir));

         ThrContainer.InvM_2PID[j] = std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>
            (("2PID_" + Par.CType.cname_nop[j]).c_str(), 
             "2PID", Par.pt_nbins, Par.pair_ptmin, Par.pair_ptmax, Par.invm_nbins, min_inv_mass, 4., dir));
         ThrContainer.InvM_2PID_acc_decreased[j] = std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>
            (("2PID_acc_decreased_" + Par.CType.cname_nop[j]).c_str(), 
             "2PID", Par.pt_nbins, Par.pair_ptmin, Par.pair_ptmax, Par.invm_nbins, min_inv_mass, 4., dir));
         ThrContainer.InvM_2PID_acc_increased[j] = std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>
            (("2PID_acc_increased_" + Par.CType.cname_nop[j]).c_str(), 
             "2PID", Par.pt_nbins, Par.pair_ptmin, Par.pair_ptmax, Par.invm_nbins, min_inv_mass, 4., dir));
         ThrContainer.InvM_2PID_m2_eff_decreased[j] = std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>
            (("2PID_m2_eff_decreased_" + Par.CType.cname_nop[j]).c_str(), 
             "2PID", Par.pt_nbins, Par.pair_ptmin, Par.pair_ptmax, Par.invm_nbins, min_inv_mass, 4., dir));
         ThrContainer.InvM_2PID_m2_eff_increased[j] = std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>
            (("2PID_m2_eff_increased_" + Par.CType.cname_nop[j]).c_str(), 
             "2PID", Par.pt_nbins, Par.pair_ptmin, Par.pair_ptmax, Par.invm_nbins, min_inv_mass, 4., dir));

         ThrContainer.InvM_TOF2PID[j] = std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>
            (("TOF2PID_" + Par.CType.cname_nop[j]).c_str(), 
             "TOF2PID", Par.pt_nbins, Par.pair_ptmin, Par.pair_ptmax, Par.invm_nbins, min_inv_mass, 4., dir));
         ThrContainer.InvM_TOF2PID_acc_decreased[j] = std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>
            (("TOF2PID_acc_decreased_" + Par.CType.cname_nop[j]).c_str(), 
             "TOF2PID", Par.pt_nbins, Par.pair_ptmin, Par.pair_ptmax, Par.invm_nbins, min_inv_mass, 4., dir));
         ThrContainer.InvM_TOF2PID_acc_increased[j] = std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>
            (("TOF2PID_acc_increased_" + Par.CType.cname_nop[j]).c_str(), 
             "TOF2PID", Par.pt_nbins, Par.pair_ptmin, Par.pair_ptmax, Par.invm_nbins, min_inv_mass, 4., dir));

         ThrContainer.InvM_EMC2PID[j] = std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>
            (("EMC2PID_" + Par.CType.cname_nop[j]).c_str(), 
             "EMC2PID", Par.pt_nbins, Par.pair_ptmin, Par.pair_ptmax, Par.invm_nbins, min_inv_mass, 4., dir));
         ThrContainer.InvM_EMC2PID_acc_decreased[j] = std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>
            (("EMC2PID_acc_decreased_" + Par.CType.cname_nop[j]).c_str(), 
             "EMC2PID", Par.pt_nbins, Par.pair_ptmin, Par.pair_ptmax, Par.invm_nbins, min_inv_mass, 4., dir));
         ThrContainer.InvM_EMC2PID_acc_increased[j] = std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>
            (("EMC2PID_acc_increased_" + Par.CType.cname_nop[j]).c_str(), 
             "EMC2PID", Par.pt_nbins, Par.pair_ptmin, Par.pair_ptmax, Par.invm_nbins, min_inv_mass, 4., dir));
         ThrContainer.InvM_EMC2PID_m2_eff_decreased[j] = std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>
            (("EMC2PID_m2_eff_decreased_" + Par.CType.cname_nop[j]).c_str(), 
             "EMC2PID", Par.pt_nbins, Par.pair_ptmin, Par.pair_ptmax, Par.invm_nbins, min_inv_mass, 4., dir));
         ThrContainer.InvM_EMC2PID_m2_eff_increased[j] = std::unique_ptr<ThrObj<TH2F>>(new ThrObj<TH2F>
            (("EMC2PID_m2_eff_increased_" + Par.CType.cname_nop[j]).c_str(), 
             "EMC2PID", Par.pt_nbins, Par.pair_ptmin, Par.pair_ptmax, Par.invm_nbins, min_inv_mass, 4., dir));
      }
      
      for (long unsigned int i = 0; i < Par.daughter1_queue.size(); i++)
      {
         for (std::string magf : Par.magf_queue)
         {   
            for (std::string auxName : Par.auxName_queue)
            {
               AnalyzeConfiguration(&ThrContainer, Par.daughter1_queue[i], Par.daughter2_queue[i], 
                  magf, auxName, ptDeviation, procNum);
               procNum++;
            }
         }   
      }
      system(("mkdir -p ../../analysis/data/phenix_sim/" + Par.run_name).c_str());

      std::string ptDeviation_name;
      if (abs(ptDeviation - 1) < 1e-7) ptDeviation_name = "";
      else ptDeviation_name = "_pt" + DtoStr(ptDeviation, 3);
   
      const std::string output_file_name = 
         "../../analysis/data/phenix_sim/" + 
         Par.run_name + "/" + Par.orig_part.name + 
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

#endif /* ANALYZE_TRACK_PAIR_CPP */
