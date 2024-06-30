#include "../lib/ParEmb.h"

#include "../lib/EmbTreeReader.h"

#include "PBar.hpp"

#include "../lib/ThrObj.h"

#include "../lib/SingleAnalysis.h"

#include "../lib/Tools.h"
#include "../lib/Box.h"
#include "../lib/InputTools.h"

#include "ROOT/TTreeProcessorMT.hxx"


struct ThrHistStruct
{
	std::unique_ptr<ThrObj<TH1F>> reg_dc_pc1, reg_pc2, reg_pc3, reg_tofe, reg_tofw;
	std::array<std::unique_ptr<ThrObj<TH1F>>, 4> reg_emcale, reg_emcalw;
	
	ThrHistStruct(std::string run_type)
	{
		reg_dc_pc1 = std::unique_ptr<ThrObj<TH1F>>(new ThrObj<TH1F>
			(("reg_dc_pc1" + run_type).c_str(), "dc_pc1", 
			 Par.centr_nbins, 0., static_cast<float>(Par.centr_nbins)));
		reg_pc2 = std::unique_ptr<ThrObj<TH1F>>(new ThrObj<TH1F>
			(("reg_pc2" + run_type).c_str(), "pc2", 
			 Par.centr_nbins, 0., static_cast<float>(Par.centr_nbins)));
		reg_pc3 = std::unique_ptr<ThrObj<TH1F>>(new ThrObj<TH1F>
			(("reg_pc3" + run_type).c_str(), "pc3", 
			 Par.centr_nbins, 0., static_cast<float>(Par.centr_nbins)));
		reg_tofe = std::unique_ptr<ThrObj<TH1F>>(new ThrObj<TH1F>
			(("reg_tofe" + run_type).c_str(), "tofe", 
			 Par.centr_nbins, 0., static_cast<float>(Par.centr_nbins)));
		reg_tofw = std::unique_ptr<ThrObj<TH1F>>(new ThrObj<TH1F>
			(("reg_tofw" + run_type).c_str(), "tofw", 
			 Par.centr_nbins, 0., static_cast<float>(Par.centr_nbins)));
		
		for (int i = 0; i < 4; i++)
		{
			reg_emcale[i] = std::unique_ptr<ThrObj<TH1F>>(new ThrObj<TH1F>
				(("reg_emcale" + std::to_string(i) + run_type).c_str(), "emcale", 
				Par.centr_nbins, 0., static_cast<float>(Par.centr_nbins)));
			reg_emcalw[i] = std::unique_ptr<ThrObj<TH1F>>(new ThrObj<TH1F>
				(("reg_emcalw" + std::to_string(i) + run_type).c_str(), "emcalw", 
				Par.centr_nbins, 0., static_cast<float>(Par.centr_nbins)));
		}
	}
};

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

void Analyze(std::string part, const int proc_num)
{
	Box box = Box("Parameters of run " + std::to_string(proc_num) + 
		" out of " + std::to_string(Par.part_queue.size()));
	
	const double charge = ParticleProperties.charge[ParticleProperties.iter_map[part]];

	std::vector<std::string> input_file_name;
	double nevents = 0;

	box.AddEntry("Run name", Par.run_name);
	box.AddEntry("Orig particle", part);
	box.AddEntry("Number of input files", static_cast<int>(
		Par.centr_nbins*Par.magf_queue.size()*Par.part_charge_queue.size()));

	for (std::string part_charge : Par.part_charge_queue)
	{
		for (std::string magf : Par.magf_queue)
		{
			for (long unsigned int i = 0; i < Par.centr_nbins; i++)
			{
				input_file_name.push_back("../data/" + Par.run_name + "/Embedding/" + 
					part_charge + part + magf + "_" + Par.centr_queue[i] +  ".root");
				TFile input_file = TFile(input_file_name[i].c_str());
				
				const double nevents_current = static_cast<double>(((TTree *) 
					input_file.Get("EmbedMcRecoTrack"))->GetEntries());
				if (nevents_current <= 0)
				{
					Print("Error: Number of events is equal or less than 0! in a file " + input_file_name[i]);
					exit(1);
				}
				nevents += nevents_current;
			}
		}
	}
	
	box.AddEntry("Number of events to be analyzed, 1e6", nevents/1e6, 3);
	box.AddEntry("Minimum p_T, GeV", Par.ptmin);
	box.AddEntry("Maximum p_T, GeV", Par.ptmax);

	box.AddEntry("Number of threads", Par.nthreads);
	
	box.Print();

	ThrHistStruct SimThrHist("_sim");
	ThrHistStruct RealThrHist("_real");

	ROOT::EnableImplicitMT(Par.nthreads);

	double ncalls = 0;
	
	bool process_finished = false;

	auto PbarCall = [&]()
	{
		ProgressBar pbar = ProgressBar("FANCY1", "", PBarColor::BOLD_RED);
		while (!process_finished)
		{
			pbar.Print(ncalls/nevents);
			std::this_thread::sleep_for(std::chrono::milliseconds(20));
		}
		pbar.Print(1.);
	};
	
	std::thread pbar_thread(PbarCall);
	
	for (long unsigned int f = 0; f < input_file_name.size(); f++)
	{
		ROOT::TTreeProcessorMT tp(input_file_name[f].c_str());
		
		auto ProcessMP = [&](TTreeReader &reader)
		{	
			//stuff to optimize multithreading
			std::shared_ptr<TH1F> reg_dc_pc1_sim = SimThrHist.reg_dc_pc1->Get();
			std::shared_ptr<TH1F> reg_pc2_sim = SimThrHist.reg_pc2->Get();
			std::shared_ptr<TH1F> reg_pc3_sim = SimThrHist.reg_pc3->Get();
			std::shared_ptr<TH1F> reg_tofe_sim = SimThrHist.reg_tofe->Get();
			std::shared_ptr<TH1F> reg_tofw_sim = SimThrHist.reg_tofw->Get();
			
			std::shared_ptr<TH1F> reg_dc_pc1_real = RealThrHist.reg_dc_pc1->Get();
			std::shared_ptr<TH1F> reg_pc2_real = RealThrHist.reg_pc2->Get();
			std::shared_ptr<TH1F> reg_pc3_real = RealThrHist.reg_pc3->Get();
			std::shared_ptr<TH1F> reg_tofe_real = RealThrHist.reg_tofe->Get();
			std::shared_ptr<TH1F> reg_tofw_real = RealThrHist.reg_tofw->Get();

			std::array<std::shared_ptr<TH1F>, 4> reg_emcale_sim, reg_emcalw_sim;
			std::array<std::shared_ptr<TH1F>, 4> reg_emcale_real, reg_emcalw_real;
			
			for (int i = 0; i < 4; i++)
			{
				reg_emcale_sim[i] = SimThrHist.reg_emcale[i]->Get();
				reg_emcalw_sim[i] = SimThrHist.reg_emcalw[i]->Get();
				reg_emcale_real[i] = RealThrHist.reg_emcale[i]->Get();
				reg_emcalw_real[i] = RealThrHist.reg_emcalw[i]->Get();
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

				bool good_real = true;
				bool good_tofw_real = false;

				if (abs(RT.zed()) < 75 && abs(RT.zed()) > 3 && !IsQualityCut(RT.qual())) 
				{
					const double the0 = RT.the0();
					const double pt = (RT.mom())*sin(the0);
					
					if (pt < Par.ptmin || pt > Par.ptmax) good_real = false;

					if (!(fabs(the0)<100 &&
						((bbcz > 0. && ((bbcz - 250.*tan(the0 - TMath::Pi()/2.)) > 2. ||
						(bbcz - 200.*tan(the0 - TMath::Pi()/2.)) < -2.)) ||
						(bbcz < 0. && ((bbcz - 250.*tan(the0 - TMath::Pi()/2.))< -2. ||
						(bbcz - 200.*tan(the0 - TMath::Pi()/2.)) > 2.))))) good_real = false;
					
					const double alpha = RT.alpha();
					const double phi = RT.phi();
               
					double board;
					if (phi>1.5) board = ((3.72402-phi+0.008047*cos(phi+0.87851))/0.01963496);
					else board = ((0.573231 + phi - 0.0046 * cos(phi + 0.05721)) / 0.01963496);

					double pc1phi = atan2(RT.ppc1y(), RT.ppc1x());
					if (phi >= 1.5 && pc1phi < 0) pc1phi += M_PI*2.;
					
					if (IsDeadDC(phi, RT.zed(), board, alpha) || 
							IsDeadPC1(phi, RT.ppc1z(), pc1phi)) good_real = false;

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
							good_tofw_real = true;
						}
					}
				}
				else good_real = false;

				if (abs(ST.zed()) < 75 && abs(ST.zed()) > 3 && !IsQualityCut(ST.qual())) 
				{
					const double the0 = ST.the0();
					const double pt = (ST.mom())*sin(the0);
					
					if (pt < Par.ptmin || pt > Par.ptmax) continue;

					if (!(fabs(the0)<100 &&
						((bbcz > 0. && ((bbcz - 250.*tan(the0 - TMath::Pi()/2.)) > 2. ||
						(bbcz - 200.*tan(the0 - TMath::Pi()/2.)) < -2.)) ||
						(bbcz < 0. && ((bbcz - 250.*tan(the0 - TMath::Pi()/2.))< -2. ||
						(bbcz - 200.*tan(the0 - TMath::Pi()/2.)) > 2.))))) continue;
					
					const double alpha = ST.alpha();
					const double phi = ST.phi();
					const double phi_real = ST.phi();

					if ((phi >= 1.5 && phi_real < 1.5) || 
						(phi < 1.5 && phi_real >= 1.5)) good_real = false;

					double board;
					if (phi>1.5) board = ((3.72402-phi+0.008047*cos(phi+0.87851))/0.01963496);
					else board = ((0.573231 + phi - 0.0046 * cos(phi + 0.05721)) / 0.01963496);

					double pc1phi = atan2(ST.ppc1y(), ST.ppc1x());
					if (phi >= 1.5 && pc1phi < 0) pc1phi += M_PI*2.;
					
					if (!IsDeadDC(phi, ST.zed(), board, alpha) && !IsDeadPC1(phi, ST.ppc1z(), pc1phi))
					{
						reg_dc_pc1_sim->Fill(f, weight);
						if (good_real) reg_dc_pc1_real->Fill(f, weight);
						
						if (IsMatch(ST.tofsdz(), ST.tofsdphi()))
						{
							const double beta = ST.pltof()/ST.ttof()/29.97;
							const double eloss = 0.0014*pow(beta, -1.66);

							if (ST.tofe() > eloss &&
								//!IsBadSlat(ST.slat()) &&
								!IsDeadTOFe(ST.zed(), ST.ptofy(), ST.ptofz()))
							{
								reg_tofe_sim->Fill(f, weight);

								if (good_real && IsMatch(RT.tofsdz(), RT.tofsdphi()))
								{
									const double beta_real = RT.pltof()/RT.ttof()/29.97;
									const double eloss_real = 0.0014*pow(beta_real, -1.66);

									if (RT.tofe() > eloss_real &&
										//!IsBadSlat(RT.slat()) &&
										!IsDeadTOFe(RT.zed(), RT.ptofy(), RT.ptofz()))
									{
										reg_tofe_real->Fill(f, weight);
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
								reg_tofw_sim->Fill(f, weight);
								if (good_tofw_real) reg_tofw_real->Fill(f, weight);
							}
						}
						if (IsMatch(ST.emcsdz(), ST.emcsdphi()) && ST.ecore() > 0.25)
						{
							const int sect = GetEmcSector(phi, ST.pemcy());
							
							if (!IsDeadEMCal(phi, ST.zed(), sect, ST.pemcy(), ST.pemcz()))
							{
								if (phi > 1.5) 
								{
									reg_emcale_sim[sect]->Fill(f, weight);
								}
								else 
								{	
									reg_emcalw_sim[sect]->Fill(f, weight);
								}
								if (good_real && IsMatch(RT.emcsdz(), RT.emcsdphi()) && RT.ecore() > 0.25)
								{
									const int sect_real = GetEmcSector(phi_real, RT.pemcy());
									
									if (sect_real == sect && !IsDeadEMCal(phi_real, RT.zed(), sect_real, RT.pemcy(), RT.pemcz()))
									{
										if (phi_real > 1.5) 
										{
											reg_emcale_real[sect_real]->Fill(f, weight);
										}
										else 
										{	
											reg_emcalw_real[sect_real]->Fill(f, weight);
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
								reg_pc2_sim->Fill(f, weight);
								
								if (good_real && IsMatch(RT.pc2sdz(), RT.pc2sdphi()))
								{
									const double pc2z_real = RT.ppc2z() - RT.pc2dz();
									double pc2phi_real = atan2(RT.ppc2y(), RT.ppc2x() - RT.pc2dphi());

									if (phi_real > 1.5 && pc2phi_real < 0) pc2phi_real += M_PI*2.;
									if (!IsDeadPC2(pc2z_real, pc2phi_real))
									{
										reg_pc2_real->Fill(f, weight);
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
								reg_pc3_sim->Fill(f, weight);
								
								if (good_real && IsMatch(RT.pc3sdz(), RT.pc3sdphi()))
								{
									const double pc3z_real = RT.ppc3z() - RT.pc3dz();
									double pc3phi_real = atan2(RT.ppc3y(), RT.ppc3x() - RT.pc3dphi());

									if (phi_real > 1.5 && pc3phi_real < 0) pc3phi_real += M_PI*2.;
									if (!IsDeadPC3(phi_real, pc3z_real, pc3phi_real))
									{
										reg_pc3_real->Fill(f, weight);
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
	
	process_finished = true;
	pbar_thread.join();
	
	std::string file_name = "../../analysis/data/phenix_sim/" + 
		Par.run_name + "/Emb_" + part + ".root";
	
	TFile outfile(file_name.c_str(), "RECREATE");

	ThrObjHolder.Write();

	outfile.Close();
	
	PrintInfo("File " + file_name + " was written");
}

void EmbAnalyze()
{
	for (std::string magf : Par.magf_queue)
	{
		for (std::string part : Par.part_queue)
		{
			for (std::string part_charge : Par.part_charge_queue)
			{
				for (std::string centr: Par.centr_queue)
				{
					CheckInputFile("../data/" + Par.run_name + 
						"/Embedding/" + part_charge + part + magf + "_" + centr + ".root");
				}
			}
		}
	}
	
	system(("mkdir -p ../../analysis/data/phenix_sim/" + Par.run_name).c_str());	
	
	int num = 1;
	for (std::string part : Par.part_queue)
	{
		Analyze(part, num);
		num++;
	}
}

int main()
{
   EmbAnalyze();
   return 0;
}
