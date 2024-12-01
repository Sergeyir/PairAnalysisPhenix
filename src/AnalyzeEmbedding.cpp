// $SOURCE$
//------------------------------------------------------------------------------------------------
//                           AnalyzeEmbedding functions realisation
//------------------------------------------------------------------------------------------------
// AnalyzeEmbedding
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic macro for embedding study
 * from simulation output of event-like TTrees to processed histograms
 **/
//------------------------------------------------------------------------------------------------

#ifndef ANALYZE_EMBEDDING_CPP
#define ANALYZE_EMBEDDING_CPP

#include "../include/AnalyzeEmbedding.hpp"

ThrContainer::ThrContainer(std::string runType)
{
   regDCPC1 = std::unique_ptr<ThrObj<TH1F>>(new ThrObj<TH1F>
      (("regDCPC1" + runType).c_str(), "dc_pc1", 
       Par.centrNBins, 0., static_cast<float>(Par.centrNBins)));
   regPC2 = std::unique_ptr<ThrObj<TH1F>>(new ThrObj<TH1F>
      (("regPC2" + runType).c_str(), "pc2", 
       Par.centrNBins, 0., static_cast<float>(Par.centrNBins)));
   regPC3 = std::unique_ptr<ThrObj<TH1F>>(new ThrObj<TH1F>
      (("regPC3" + runType).c_str(), "pc3", 
       Par.centrNBins, 0., static_cast<float>(Par.centrNBins)));
   regTOFe = std::unique_ptr<ThrObj<TH1F>>(new ThrObj<TH1F>
      (("regTOFe" + runType).c_str(), "tofe", 
       Par.centrNBins, 0., static_cast<float>(Par.centrNBins)));
   regTOFw = std::unique_ptr<ThrObj<TH1F>>(new ThrObj<TH1F>
      (("regTOFw" + runType).c_str(), "tofw", 
       Par.centrNBins, 0., static_cast<float>(Par.centrNBins)));
   
   for (int i = 0; i < 4; i++)
   {
      regEMCale[i] = std::unique_ptr<ThrObj<TH1F>>(new ThrObj<TH1F>
         (("regEMCale" + std::to_string(i) + runType).c_str(), "emcale", 
         Par.centrNBins, 0., static_cast<float>(Par.centrNBins)));
      regEMCalw[i] = std::unique_ptr<ThrObj<TH1F>>(new ThrObj<TH1F>
         (("regEMCalw" + std::to_string(i) + runType).c_str(), "emcalw", 
         Par.centrNBins, 0., static_cast<float>(Par.centrNBins)));
   }
}

int GetEmcSector(const double phi, const double pemcy)
{
   if (phi > 1.5)
   {
      if (pemcy < -98) return 0;
      else if (pemcy < 102) return 1;
      else if (pemcy < 290) return 2;
      else return 3;
   }
   else   
   {
      if (pemcy < -100) return 0;
      else if (pemcy < 102) return 1;
      else if (pemcy < 290) return 2;
      else return 3;
   }
}

void AnalyzeParticleEmbedding(std::string part, const int queueNum)
{
   Box box = Box("Parameters of run " + std::to_string(queueNum) + 
      " out of " + std::to_string(Par.partQueue.size()));
   
   const double charge = ParticleProperties.charge[ParticleProperties.iterMap[part]];

   std::vector<std::string> inputFileName;
   double nevents = 0;

   box.AddEntry("Run name", Par.runName);
   box.AddEntry("Orig particle", part);
   box.AddEntry("Number of input files", static_cast<int>(
      Par.centrNBins*Par.magfQueue.size()*Par.partChargeQueue.size()));

   for (std::string partCharge : Par.partChargeQueue)
   {
      for (std::string magf : Par.magfQueue)
      {
         for (long unsigned int i = 0; i < Par.centrNBins; i++)
         {
            inputFileName.push_back(Par.dataDir + Par.runName + "/Embedding/" + 
               partCharge + part + magf + "_" + Par.centrQueue[i] +  ".root");
            
            TFile inputFile = TFile(inputFileName.back().c_str());
            
            const double neventsCurrent = static_cast<double>(((TTree *) 
               inputFile.Get("EmbedMcRecoTrack"))->GetEntries());
            if (neventsCurrent <= 0)
            {
               Print("Error: Number of events is equal or less than 0! in a file " + inputFileName.back());
               exit(1);
            }
            nevents += neventsCurrent;
         }
      }
   }
   
   box.AddEntry("Number of events to be analyzed, 1e6", nevents/1e6, 3);
   box.AddEntry("Minimum p_T, GeV", Par.ptMin);
   box.AddEntry("Maximum p_T, GeV", Par.ptMax);

   box.AddEntry("Number of threads", Par.nthreads);
   
   box.Print();

   ThrContainer SimThrContainer("_sim");
   ThrContainer RealThrContainer("_real");

   ROOT::EnableImplicitMT(Par.nthreads);

   double ncalls = 0;
   
   bool isProcessFinished = false;

   auto PbarCall = [&]()
   {
      ProgressBar pBar = ProgressBar("FANCY1", "", PBarColor::BOLD_RED);
      while (!isProcessFinished)
      {
         pBar.Print(ncalls/nevents);
         std::this_thread::sleep_for(std::chrono::milliseconds(20));
      }
      pBar.Print(1.);
   };
   
   std::thread pBarThread(PbarCall);
   
   for (long unsigned int f = 0; f < inputFileName.size(); f++)
   {
      ROOT::TTreeProcessorMT tp(inputFileName[f].c_str());
      
      auto ProcessMP = [&](TTreeReader &reader)
      {   
         //stuff to optimize multithreading
         std::shared_ptr<TH1F> regDCPC1Sim = SimThrContainer.regDCPC1->Get();
         std::shared_ptr<TH1F> regPC2Sim = SimThrContainer.regPC2->Get();
         std::shared_ptr<TH1F> regPC3Sim = SimThrContainer.regPC3->Get();
         std::shared_ptr<TH1F> regTOFeSim = SimThrContainer.regTOFe->Get();
         std::shared_ptr<TH1F> regTOFwSim = SimThrContainer.regTOFw->Get();
         
         std::shared_ptr<TH1F> regDCPC1Real = RealThrContainer.regDCPC1->Get();
         std::shared_ptr<TH1F> regPC2Real = RealThrContainer.regPC2->Get();
         std::shared_ptr<TH1F> regPC3Real = RealThrContainer.regPC3->Get();
         std::shared_ptr<TH1F> regTOFeReal = RealThrContainer.regTOFe->Get();
         std::shared_ptr<TH1F> regTOFwReal = RealThrContainer.regTOFw->Get();

         std::array<std::shared_ptr<TH1F>, 4> regEMCaleSim, regEMCalwSim;
         std::array<std::shared_ptr<TH1F>, 4> regEMCaleReal, regEMCalwReal;
         
         for (int i = 0; i < 4; i++)
         {
            regEMCaleSim[i] = SimThrContainer.regEMCale[i]->Get();
            regEMCalwSim[i] = SimThrContainer.regEMCalw[i]->Get();
            regEMCaleReal[i] = RealThrContainer.regEMCale[i]->Get();
            regEMCalwReal[i] = RealThrContainer.regEMCalw[i]->Get();
         }
         
         EventReader ET(reader);
         EmbTreeReader ST(reader, "S");
         EmbTreeReader RT(reader, "R");
         
         while (reader.Next())
         {
            ncalls += 1.;
            
            if (ET.type() != 2) continue;
            const double weight = 1.;

            const double bbcz = ET.bbcz();

            bool isGoodReal = true;
            bool isGoodTOFwReal = false;

            if (abs(RT.zed()) < 75 && abs(RT.zed()) > 3 && !IsQualityCut(RT.qual())) 
            {
               const double the0 = RT.the0();
               const double pt = (RT.mom())*sin(the0);
               
               if (pt < Par.ptMin || pt > Par.ptMax) isGoodReal = false;

               if (!(fabs(the0)<100 &&
                  ((bbcz > 0. && ((bbcz - 250.*tan(the0 - TMath::Pi()/2.)) > 2. ||
                  (bbcz - 200.*tan(the0 - TMath::Pi()/2.)) < -2.)) ||
                  (bbcz < 0. && ((bbcz - 250.*tan(the0 - TMath::Pi()/2.))< -2. ||
                  (bbcz - 200.*tan(the0 - TMath::Pi()/2.)) > 2.))))) isGoodReal = false;
               
               const double alpha = RT.alpha();
               const double phi = RT.phi();
               
               double board;
               if (phi>1.5) board = ((3.72402-phi+0.008047*cos(phi+0.87851))/0.01963496);
               else board = ((0.573231 + phi - 0.0046 * cos(phi + 0.05721)) / 0.01963496);

               double pc1phi = atan2(RT.ppc1y(), RT.ppc1x());
               if (phi >= 1.5 && pc1phi < 0) pc1phi += M_PI*2.;
               
               if (IsDeadDC(phi, RT.zed(), board, alpha) || 
                     IsDeadPC1(phi, RT.ppc1z(), pc1phi)) isGoodReal = false;

               if (IsHit(RT.tofwdz()))
               {
                  double tofwsdphi = GetTOFwsdphi(0, RT.mom(), RT.tofwdphi(), charge, RT.striptofw());
                  double tofwsdz = GetTOFwsdz(0, RT.mom(), RT.tofwdz(), charge, RT.striptofw());
                  
                  tofwsdphi = RecalTOFwsdphi(0, RT.mom(), tofwsdphi, charge, RT.striptofw());
                  tofwsdz = RecalTOFwsdz(0, RT.mom(), tofwsdz, charge, RT.striptofw());
                  
                  if (IsMatch(tofwsdphi, tofwsdz) && 
                     !IsBadStripTOFw(static_cast<int>(RT.striptofw())) &&
                     !IsDeadTOFw(RT.zed(), board, RT.alpha()))
                  {
                     isGoodTOFwReal = true;
                  }
               }
            }
            else isGoodReal = false;

            if (abs(ST.zed()) < 75 && abs(ST.zed()) > 3 && !IsQualityCut(ST.qual())) 
            {
               const double the0 = ST.the0();
               const double pt = (ST.mom())*sin(the0);
               
               if (pt < Par.ptMin || pt > Par.ptMax) continue;

               if (!(fabs(the0)<100 &&
                  ((bbcz > 0. && ((bbcz - 250.*tan(the0 - TMath::Pi()/2.)) > 2. ||
                  (bbcz - 200.*tan(the0 - TMath::Pi()/2.)) < -2.)) ||
                  (bbcz < 0. && ((bbcz - 250.*tan(the0 - TMath::Pi()/2.))< -2. ||
                  (bbcz - 200.*tan(the0 - TMath::Pi()/2.)) > 2.))))) continue;
               
               const double alpha = ST.alpha();
               const double phi = ST.phi();
               const double phiReal = ST.phi();

               if ((phi >= 1.5 && phiReal < 1.5) || 
                  (phi < 1.5 && phiReal >= 1.5)) isGoodReal = false;

               double board;
               if (phi>1.5) board = ((3.72402-phi+0.008047*cos(phi+0.87851))/0.01963496);
               else board = ((0.573231 + phi - 0.0046 * cos(phi + 0.05721)) / 0.01963496);

               double pc1phi = atan2(ST.ppc1y(), ST.ppc1x());
               if (phi >= 1.5 && pc1phi < 0) pc1phi += M_PI*2.;
               
               if (!IsDeadDC(phi, ST.zed(), board, alpha) && !IsDeadPC1(phi, ST.ppc1z(), pc1phi))
               {
                  regDCPC1Sim->Fill(f, weight);
                  if (isGoodReal) regDCPC1Real->Fill(f, weight);
                  
                  if (IsMatch(ST.tofsdz(), ST.tofsdphi()))
                  {
                     const double beta = ST.pltof()/ST.ttof()/29.97;
                     const double eloss = 0.0014*pow(beta, -1.66);

                     if (ST.tofe() > eloss &&
                        //!IsBadSlat(ST.slat()) &&
                        !IsDeadTOFe(ST.zed(), ST.ptofy(), ST.ptofz()))
                     {
                        regTOFeSim->Fill(f, weight);

                        if (isGoodReal && IsMatch(RT.tofsdz(), RT.tofsdphi()))
                        {
                           const double betaReal = RT.pltof()/RT.ttof()/29.97;
                           const double elossReal = 0.0014*pow(betaReal, -1.66);

                           if (RT.tofe() > elossReal &&
                              //!IsBadSlat(RT.slat()) &&
                              !IsDeadTOFe(RT.zed(), RT.ptofy(), RT.ptofz()))
                           {
                              regTOFeReal->Fill(f, weight);
                           }
                        }
                        
                     }
                  }
                  if (IsHit(ST.tofwdz()))
                  {
                     double tofwsdphi = GetTOFwsdphi(0, ST.mom(), ST.tofwdphi(), charge, ST.striptofw());
                     double tofwsdz = GetTOFwsdz(0, ST.mom(), ST.tofwdz(), charge, ST.striptofw());
                     
                     tofwsdphi = RecalTOFwsdphi(0, ST.mom(), tofwsdphi, charge, ST.striptofw());
                     tofwsdz = RecalTOFwsdz(0, ST.mom(), tofwsdz, charge, ST.striptofw());
                     
                     if (IsMatch(tofwsdphi, tofwsdz) && 
                        !IsBadStripTOFw(static_cast<int>(ST.striptofw())) &&
                        !IsDeadTOFw(ST.zed(), board, ST.alpha()))
                     {
                        regTOFwSim->Fill(f, weight);
                        if (isGoodTOFwReal) regTOFwReal->Fill(f, weight);
                     }
                  }
                  if (IsMatch(ST.emcsdz(), ST.emcsdphi()) && ST.ecore() > 0.25)
                  {
                     const int sect = GetEmcSector(phi, ST.pemcy());
                     
                     if (!IsDeadEMCal(phi, ST.zed(), sect, ST.pemcy(), ST.pemcz()))
                     {
                        if (phi > 1.5) 
                        {
                           regEMCaleSim[sect]->Fill(f, weight);
                        }
                        else 
                        {   
                           regEMCalwSim[sect]->Fill(f, weight);
                        }
                        if (isGoodReal && IsMatch(RT.emcsdz(), RT.emcsdphi()) && RT.ecore() > 0.25)
                        {
                           const int sectReal = GetEmcSector(phiReal, RT.pemcy());
                           
                           if (sectReal == sect && !IsDeadEMCal(phiReal, RT.zed(), sectReal, RT.pemcy(), RT.pemcz()))
                           {
                              if (phiReal > 1.5) 
                              {
                                 regEMCaleReal[sectReal]->Fill(f, weight);
                              }
                              else 
                              {   
                                 regEMCalwReal[sectReal]->Fill(f, weight);
                              }
                           }
                        }
                     }
                  }
                  if (IsMatch(ST.pc2sdz(), ST.pc2sdphi()))
                  {
                     const double pc2z = ST.ppc2z() - ST.pc2dz();
                     double pc2phi = atan2(ST.ppc2y(), ST.ppc2x() - ST.pc2dphi());

                     if (phi > 1.5 && pc2phi < 0) pc2phi += M_PI*2.;
                     if (!IsDeadPC2(pc2z, pc2phi))
                     {
                        regPC2Sim->Fill(f, weight);
                        
                        if (isGoodReal && IsMatch(RT.pc2sdz(), RT.pc2sdphi()))
                        {
                           const double pc2zReal = RT.ppc2z() - RT.pc2dz();
                           double pc2phiReal = atan2(RT.ppc2y(), RT.ppc2x() - RT.pc2dphi());

                           if (phiReal > 1.5 && pc2phiReal < 0) pc2phiReal += M_PI*2.;
                           if (!IsDeadPC2(pc2zReal, pc2phiReal))
                           {
                              regPC2Real->Fill(f, weight);
                           }
                        }
                     }
                  }
                  if (IsMatch(ST.pc3sdz(), ST.pc3sdphi()))
                  {
                     const double pc3z = ST.ppc3z() - ST.pc3dz();
                     double pc3phi = atan2(ST.ppc3y(), ST.ppc3x() - ST.pc3dphi());

                     if (phi > 1.5 && pc3phi < 0) pc3phi += M_PI*2.;
                     if (!IsDeadPC3(phi, pc3z, pc3phi))
                     {
                        regPC3Sim->Fill(f, weight);
                        
                        if (isGoodReal && IsMatch(RT.pc3sdz(), RT.pc3sdphi()))
                        {
                           const double pc3zReal = RT.ppc3z() - RT.pc3dz();
                           double pc3phiReal = atan2(RT.ppc3y(), RT.ppc3x() - RT.pc3dphi());

                           if (phiReal > 1.5 && pc3phiReal < 0) pc3phiReal += M_PI*2.;
                           if (!IsDeadPC3(phiReal, pc3zReal, pc3phiReal))
                           {
                              regPC3Real->Fill(f, weight);
                           }
                        }
                     }
                  }
               }
            }
         }
      };
   
      tp.Process(ProcessMP);
   }
   
   isProcessFinished = true;
   pBarThread.join();
   
   std::string fileName = Par.outputDir + 
      Par.runName + "/Embedding/" + part + ".root";
   
   TFile outfile(fileName.c_str(), "RECREATE");

   ThrObjHolder.Write();

   outfile.Close();
   
   PrintInfo("File " + fileName + " was written");
}

void AnalyzeEmbedding()
{
   for (std::string magf : Par.magfQueue)
   {
      for (std::string part : Par.partQueue)
      {
         for (std::string partCharge : Par.partChargeQueue)
         {
            for (std::string centr: Par.centrQueue)
            {
               CheckInputFile(Par.dataDir + Par.runName + "/Embedding/" + 
                              partCharge + part + magf + "_" + centr + ".root");
            }
         }
      }
   }
   
   system(("mkdir -p " + Par.outputDir + Par.runName + "/Embedding").c_str());   
   
   int queueNum = 1;
   for (std::string part : Par.partQueue)
   {
      AnalyzeParticleEmbedding(part, queueNum);
      queueNum++;
   }
}

int main()
{
   AnalyzeEmbedding();
   return 0;
}

#endif /*ANALYZE_EMBEDDING_CPP*/
