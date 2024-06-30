#include "../lib/ParSingle.h"

#include "../lib/EffTreeReader.h"

#include "PBar.hpp"

#include "../lib/ThrObj.h"

#include "../lib/SingleAnalysis.h"

#include "../lib/Tools.h"
#include "../lib/Box.h"
#include "../lib/InputTools.h"

#include "ROOT/TTreeProcessorMT.hxx"

struct ThrHistStruct
{
	ThrObj<TH2F> npart_hist = ThrObj<TH2F>
		("npart", "npart", 50, 0, 50, 100, Par.ptmin, Par.ptmax);

	ThrObj<TH1F> orig = ThrObj<TH1F>
		("orig","orig", Par.pt_nbins, Par.ptmin, Par.ptmax);

	ThrObj<TH2F> m2_tofw = ThrObj<TH2F>
		("m2_tofw", "m2", Par.pt_nbins, Par.ptmin, Par.ptmax, 2000, -4., 4.);
	ThrObj<TH2F> m2_tofe = ThrObj<TH2F>
		("m2_tofe", "m2", Par.pt_nbins, Par.ptmin, Par.ptmax, 2000, -4., 4.);

	ThrObj<TH2F> t_tofe = ThrObj<TH2F>
		("t_tofe", "t", Par.pt_nbins, Par.ptmin, Par.ptmax, 2000, -40, 40.);
	ThrObj<TH2F> t_tofw = ThrObj<TH2F>
		("t_tofw", "t", Par.pt_nbins, Par.ptmin, Par.ptmax, 2000, -40, 40.);

	ThrObj<TH2F> orig_pt_vs_pt = ThrObj<TH2F>
		("orig_pt_vs_pt", "pt", Par.pt_nbins, Par.ptmin, Par.ptmax, Par.pt_nbins, Par.ptmin, Par.ptmax);

	ThrObj<TH1F> reg_tofe = ThrObj<TH1F>
		("reg_tofe", "tofe", Par.pt_nbins, Par.ptmin, Par.ptmax);
	ThrObj<TH1F> reg_tofw = ThrObj<TH1F>
		("reg_tofw", "tofw", Par.pt_nbins, Par.ptmin, Par.ptmax);
	
	std::array<ThrObj<TH1F>, 4> reg_emcale = {
		ThrObj<TH1F>("reg_emcale0", "emc", Par.pt_nbins, Par.ptmin, Par.ptmax),
		ThrObj<TH1F>("reg_emcale1", "emc", Par.pt_nbins, Par.ptmin, Par.ptmax),
		ThrObj<TH1F>("reg_emcale2", "emc", Par.pt_nbins, Par.ptmin, Par.ptmax),
		ThrObj<TH1F>("reg_emcale3", "emc", Par.pt_nbins, Par.ptmin, Par.ptmax)};

	std::array<ThrObj<TH1F>, 4> reg_emcalw = {
		ThrObj<TH1F>("reg_emcalw0", "emc", Par.pt_nbins, Par.ptmin, Par.ptmax),
		ThrObj<TH1F>("reg_emcalw1", "emc", Par.pt_nbins, Par.ptmin, Par.ptmax),
		ThrObj<TH1F>("reg_emcalw2", "emc", Par.pt_nbins, Par.ptmin, Par.ptmax),
		ThrObj<TH1F>("reg_emcalw3", "emc", Par.pt_nbins, Par.ptmin, Par.ptmax)};
};

void Analyze(ThrHistStruct *ThrHist, std::string part, std::string magf, std::string aux_name, const int proc_num)
{
	Box box = Box("Parameters of run " + std::to_string(proc_num) + 
		" out of " + std::to_string(Par.part_queue.size()*Par.magf_queue.size()*Par.aux_name_queue.size()));
	
	double nparticles = 0;
	
	const int part_iter = ParticleProperties.iter_map[part];
	const int part_geant_id = ParticleProperties.geant_id[part_iter];
	const double part_mass = ParticleProperties.mass[part_iter];
	const int part_charge = ParticleProperties.charge[part_iter];

	const std::string input_file_name = "../data/" + Par.run_name + "/Single/" + part + magf + aux_name + ".root";
	const std::string se_file_name = "../../analysis/data/" + Par.run_name + "/sum" + magf + ".root";

	TFile input_file = TFile(input_file_name.c_str());
	TFile se_file = TFile(se_file_name.c_str());

	const double nevents = static_cast<double>(((TTree *) input_file.Get("Tree"))->GetEntries());

	if (nevents <= 0)
	{
		Print("Error: Number of events is equal or less than 0!");
		exit(1);
	}

	TH1F *orig_hist = (TH1F *) input_file.Get("orig_pt");

	const double orig_pt_threshold = orig_hist->Integral()/
		static_cast<double>(orig_hist->GetXaxis()->GetNbins())/2.;

	const double low_pt_bound = orig_hist->GetXaxis()->GetBinLowEdge(
		orig_hist->FindFirstBinAbove(orig_pt_threshold));
	const double up_pt_bound = orig_hist->GetXaxis()->GetBinUpEdge(
		orig_hist->FindLastBinAbove(orig_pt_threshold));

	double event_norm = 1.;
	if (Par.do_use_weight_func)
	{
		TH1F *centr_hist = (TH1F *) se_file.Get("central_bin");
		
		event_norm = orig_hist->Integral(
			orig_hist->GetXaxis()->FindBin(Par.ptmin),
			orig_hist->GetXaxis()->FindBin(Par.ptmax))/
			centr_hist->Integral(1, centr_hist->GetXaxis()->GetNbins());

			//this normalization is needed to merge 2 files with flapt pt distributin with different ranges
			event_norm *= (Par.ptmax - Par.ptmin)/(up_pt_bound - low_pt_bound);
	}

	box.AddEntry("Run name", Par.run_name);
	box.AddEntry("Orig particle", part);
	box.AddEntry("Magnetic field", magf);
	box.AddEntry("Orig pT distribution span", 
		DtoStr(low_pt_bound) + " < pT < " + DtoStr(up_pt_bound));
	box.AddEntry("Use weight function", Par.do_use_weight_func);
	
	box.AddEntry("Minimum p_T, GeV", Par.ptmin);
	box.AddEntry("Maximum p_T, GeV", Par.ptmax);
	box.AddEntry("Number of events to be analyzed, 1e6", nevents/1e6, 3);

	box.AddEntry("Number of threads", Par.nthreads);
	
	box.Print();

	//tsallis weight function
	std::unique_ptr<TF1> weight_func(new TF1("weight_func", 
		"[0]*x^([5]-1)*(1 + ([1] - 1)*(sqrt([4] + x^2)^([5]) - [2])/[3])^(-1./([1]-1.))"));
	weight_func->SetParameters(
		ReadFileIntoArray("../input/Spectra/" + Par.system + "/" + part + ".txt", 6));

	ROOT::EnableImplicitMT(Par.nthreads);
	ROOT::TTreeProcessorMT tp(input_file_name.c_str());
	
	double ncalls = 0;

	bool process_finished = false;
	auto ProcessMP = [&](TTreeReader &reader)
	{	
		std::shared_ptr<TH2> npart_hist = ThrHist->npart_hist.Get();
		
		std::shared_ptr<TH1F> orig = ThrHist->orig.Get();

		std::shared_ptr<TH2F> m2_tofw = ThrHist->m2_tofw.Get();
		std::shared_ptr<TH2F> m2_tofe = ThrHist->m2_tofe.Get();
		std::shared_ptr<TH2F> t_tofe = ThrHist->t_tofe.Get();
		std::shared_ptr<TH2F> t_tofw = ThrHist->t_tofw.Get();

		std::shared_ptr<TH2F> orig_pt_vs_pt = ThrHist->orig_pt_vs_pt.Get();

		std::shared_ptr<TH1F> reg_tofe = ThrHist->reg_tofe.Get();
		std::shared_ptr<TH1F> reg_tofw = ThrHist->reg_tofw.Get();
		
		std::array<std::shared_ptr<TH1F>, 4> reg_emcale, reg_emcalw;
		
		for (int i = 0; i < 4; i++)
		{
			reg_emcale[i] = ThrHist->reg_emcale[i].Get();
			reg_emcalw[i] = ThrHist->reg_emcalw[i].Get();
		}
		
		EffTreeReader T(reader);
		
		while (reader.Next())
		{
			ncalls += 1.;
			const double orig_pt = sqrt(pow(T.mom_orig(0), 2) + pow(T.mom_orig(1), 2));
			
			double event_weight;
			if (Par.do_use_weight_func) 
			{
				event_weight = weight_func->Eval(orig_pt)/event_norm;
			}
			else event_weight = 1.;
			
			orig->Fill(orig_pt, event_weight);
			
			npart_hist->Fill(T.nch()-0.5, orig_pt, event_weight);
			
			const double bbcz = T.bbcz();
			if (fabs(bbcz) > 30) continue;

			if (T.nch() <= 0 || T.nch() > 49) continue;

			nparticles += T.nch()*event_weight;

			for(int i = 0; i < T.nch(); i++)
			{
				const double the0 = T.the0(i);
				const double pt = (T.mom(i))*sin(the0);
				
				if (pt < Par.ptmin || pt > Par.ptmax) continue;
				
				const int charge = T.charge(i);
				if (charge != part_charge) continue;
					
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

				double particle_weight = event_weight;
				
				double board;
				if (phi>1.5) board = ((3.72402-phi+0.008047*cos(phi+0.87851))/0.01963496);
				else board = ((0.573231 + phi - 0.0046 * cos(phi + 0.05721)) / 0.01963496);

				if (IsDeadDC(phi, zed, board, alpha)) continue;

				double pc1phi = atan2(T.ppc1y(i), T.ppc1x(i));
				if (phi >= 1.5 && pc1phi < 0) pc1phi += M_PI*2.;

				if (IsDeadPC1(phi, T.ppc1z(i), pc1phi)) continue;

				orig_pt_vs_pt->Fill(pt, orig_pt, particle_weight);

				if (IsHit(T.tofdz(i)))
				{
					const double beta = T.pltof(i)/T.ttof(i)/29.97;
					const double eloss = 0.0014*pow(beta, -1.66);

					if (IsMatch(T.tofsdz(i), T.tofsdphi(i)) &&
						T.etof(i) > eloss &&
						!IsBadSlat(T.slat(i)) &&
						!IsDeadTOFe(zed, T.ptofy(i), T.ptofz(i)))
					{
						if (T.particle_id(i) == part_geant_id && T.primary_id(i) == -999)
						{
							reg_tofe->Fill(pt, particle_weight);
						}
						
						const double t_exp = sqrt(pow(part_mass/T.mom(i), 2) + 1.)*T.pltof(i)/29.979;
						const double m2 = pow(T.mom(i), 2)*(pow((T.ttof(i))*29.979/T.pltof(i), 2) - 1.);
						
						t_tofe->Fill(pt, T.ttof(i)-t_exp, particle_weight);
						m2_tofe->Fill(pt, m2, particle_weight*0.909);
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
						const double weight_corr = 0.7796;
						
						double pltofw = T.pltofw(i);
						const int ichamber = int(T.striptofw(i)/4)%32;

						if (ichamber<16&&ichamber%2==1) pltofw += 3.358;
						else if (ichamber>=16&&ichamber%2==0) pltofw += 3.358;
						double ttofw = T.ttofw(i)-1.12;
						
						const double t_exp = sqrt(
							pow(part_mass/T.mom(i), 2) + 1.)*pltofw/29.979;	
						t_tofw->Fill(pt, ttofw-t_exp, particle_weight*weight_corr);
						
						if (T.particle_id(i) == part_geant_id && T.primary_id(i) == -999)
						{
							reg_tofw->Fill(pt, particle_weight*weight_corr);
						}
						
						const double m2 = pow(T.mom(i), 2)*
							(pow((ttofw)*29.979/pltofw, 2) - 1.);
						
						m2_tofw->Fill(pt, m2, particle_weight*weight_corr);
					}
				}

				if (IsHit(T.emcdz(i)))
				{
					if (T.particle_id(i) == part_geant_id && T.primary_id(i) == -999 &&
						IsMatch(T.emcsdz(i), T.emcsdphi(i), 2., 2.) && T.ecore(i) > 0.25 &&
						!IsDeadEMCal(phi, zed, T.sect(i), T.pemcy(i), T.pemcz(i)))
					{
						if (phi > 1.5) 
						{
							reg_emcale[T.sect(i)]->Fill(pt, particle_weight);
						}
						else 
						{	
							reg_emcalw[T.sect(i)]->Fill(pt, particle_weight);
						}
					}
				}
			}
		}
	};
	
	auto PbarCall = [&]()
	{
		ProgressBar pbar = ProgressBar("Block");
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

void SingleAnalyze()
{
	if (Par.do_use_weight_func)
	{
		for (std::string magf : Par.magf_queue)
		{
			CheckInputFile("../../analysis/data/" + Par.run_name + "/sum" + magf + ".root");
			for (std::string part : Par.part_queue)
			{
				for (std::string aux_name : Par.aux_name_queue)
				{
					CheckInputFile("../data/" + Par.run_name + "/Single/" + part + magf + aux_name + ".root");
				}
			}
		}
		for (std::string part : Par.part_queue)
		{
			CheckInputFile("../input/Spectra/" + Par.system + "/" + part + ".txt");
		}
	}
	
	system(("mkdir -p ../../analysis/data/phenix_sim/" + Par.run_name).c_str());	
	
	int num = 1;
	for (std::string part : Par.part_queue)
	{
		for (std::string magf : Par.magf_queue)
		{
			ThrHistStruct ThrHist;
			for (std::string aux_name : Par.aux_name_queue)
			{
				Analyze(&ThrHist, part, magf, aux_name, num);
				num++;
			}
			const std::string output_file_name = 
				"../../analysis/data/phenix_sim/" + Par.run_name + "/" + part + magf + ".root";

			TFile outfile(output_file_name.c_str(), "RECREATE");
			outfile.cd();
			ThrObjHolder.Write();
			outfile.Close();
			PrintInfo("File " + output_file_name + " was written");
		}
		system(("hadd -f ../../analysis/data/phenix_sim/" + Par.run_name + "/" + 
			part + ".root ../../analysis/data/phenix_sim/" + Par.run_name + "/" + part + "*").c_str());
	}
}

int main()
{
   SingleAnalyze();
   return 0;
}
