#include <signal.h>

#include "../lib/ParHeatMapper.h"

#include "../lib/SingleAnalysis.h"
#include "../lib/EffTreeReader.h"
#include "../lib/ThrObj.h"

#include "PBar.hpp"

#include "../lib/Tools.h"
#include "../lib/Box.h"
#include "../lib/InputTools.h"
#include "../lib/StrTool.h"

#include "ROOT/TTreeProcessorMT.hxx"

struct ThrHistStruct
{
	ThrObj<TH1F> orig = ThrObj<TH1F>
		("orig","orig", Par.pt_nbins, Par.ptmin, Par.ptmax);

	ThrObj<TH2F> orig_pt_vs_pt = ThrObj<TH2F>
		("orig_pt_vs_pt","orig_pt vs pt", 
		100, 0., 10., Par.pt_nbins, Par.ptmin, Par.ptmax);
		
	ThrObj<TH2F> dce0 = ThrObj<TH2F>
		("dceast0", "map", 400, -1.5, 80.5, 200, -0.3, 0.3);
	ThrObj<TH2F> dcw0 = ThrObj<TH2F>
		("dcwest0", "map", 400, -1.5, 80.5, 200, -0.3, 0.3);
	ThrObj<TH2F> dce1 = ThrObj<TH2F>
		("dceast1", "map", 400, -1.5, 80.5, 200, -0.3, 0.3);
	ThrObj<TH2F> dcw1 = ThrObj<TH2F>
		("dcwest1", "map", 400, -1.5, 80.5, 200, -0.3, 0.3);

	ThrObj<TH2F> dce0_unsc = ThrObj<TH2F>
		("unscaled_dceast0", "map", 400, -1.5, 80.5, 200, -0.3, 0.3);
	ThrObj<TH2F> dcw0_unsc = ThrObj<TH2F>
		("unscaled_dcwest0", "map", 400, -1.5, 80.5, 200, -0.3, 0.3);
	ThrObj<TH2F> dce1_unsc = ThrObj<TH2F>
		("unscaled_dceast1", "map", 400, -1.5, 80.5, 200, -0.3, 0.3);
	ThrObj<TH2F> dcw1_unsc = ThrObj<TH2F>
		("unscaled_dcwest1", "map", 400, -1.5, 80.5, 200, -0.3, 0.3);

	ThrObj<TH2F> pc1e_z_vs_phi = ThrObj<TH2F>
		("pc1e_z_vs_phi", "z vs phi", 180, -90., 90., 330, 2.1, 3.75);
	ThrObj<TH2F> pc1w_z_vs_phi = ThrObj<TH2F>
		("pc1w_z_vs_phi", "z vs phi", 180, -90., 90., 330, -0.6, 1.05);

	ThrObj<TH2F> pc2_z_vs_phi = ThrObj<TH2F>
		("pc2_z_vs_phi", "z vs phi", 160, -160., 160., 330, -0.6, 1.05);

	ThrObj<TH2F> pc3e_z_vs_phi = ThrObj<TH2F>
		("pc3e_z_vs_phi", "z vs phi", 190, -190., 190., 320, 2.15, 3.75);
	ThrObj<TH2F> pc3w_z_vs_phi = ThrObj<TH2F>
		("pc3w_z_vs_phi", "z vs phi", 190, -190., 190., 330, -0.6, 1.05);

	std::array<ThrObj<TH2F>, 4> emcale_pos = 
	{
		ThrObj<TH2F>("emcale0_pos", "y vs z", 100, -300., -90., 210, -30., 210.),
		ThrObj<TH2F>("emcale1_pos", "y vs z", 100, -110., 110., 210, -30., 210.),
		ThrObj<TH2F>("emcale2_pos", "y vs z", 100, 90., 300., 210, -30., 210.),
		ThrObj<TH2F>("emcale3_pos", "y vs z", 100, 280., 440., 210, -30., 210.)
	};

	std::array<ThrObj<TH2F>, 4> emcale_neg = 
	{
		ThrObj<TH2F>("emcale0_neg", "y vs z", 100, -300., -90., 210, -210., 30.),
		ThrObj<TH2F>("emcale1_neg", "y vs z", 100, -110., 110., 210, -210., 30.),
		ThrObj<TH2F>("emcale2_neg", "y vs z", 100, 90., 300., 210, -210., 30.),
		ThrObj<TH2F>("emcale3_neg", "y vs z", 100, 280., 440., 210, -210., 30.)
	};
	
	std::array<ThrObj<TH2F>, 4> emcalw_pos = 
	{
		ThrObj<TH2F>("emcalw0_pos", "y vs z", 100, -300., -90., 210, -30., 210.),
		ThrObj<TH2F>("emcalw1_pos", "y vs z", 100, -110., 110., 210, -30., 210.),
		ThrObj<TH2F>("emcalw2_pos", "y vs z", 100, 90., 300., 210, -30., 210.),
		ThrObj<TH2F>("emcalw3_pos", "y vs z", 100, 280., 440., 210, -30., 210.)
	};

	std::array<ThrObj<TH2F>, 4> emcalw_neg = 
	{
		ThrObj<TH2F>("emcalw0_neg", "y vs z", 100, -300., -90., 210, -210., 30.),
		ThrObj<TH2F>("emcalw1_neg", "y vs z", 100, -110., 110., 210, -210., 30.),
		ThrObj<TH2F>("emcalw2_neg", "y vs z", 100, 90., 300., 210, -210., 30.),
		ThrObj<TH2F>("emcalw3_neg", "y vs z", 100, 280., 440., 210, -210., 30.)
	};

	ThrObj<TH1F> strip_tofw_hist = ThrObj<TH1F>
		("strip_tofw", "strip", 550, 0., 550.);
	ThrObj<TH1F> slat_hist = ThrObj<TH1F>
		("slat", "slat", 1000, 0., 1000.);

	ThrObj<TH2F> tofe0 = ThrObj<TH2F>
		("tofe0", "tofe0", 200, -300., 100., 200, -40., 210.);
	ThrObj<TH2F> tofe1 = ThrObj<TH2F>
		("tofe1", "tofe1", 200, -300., 100., 200, -210., 40.);

	ThrObj<TH2F> tofw0 = ThrObj<TH2F>
		("tofw0", "tofw0", 200, 15, 68, 200, -0.3, 0.3);
	ThrObj<TH2F> tofw1 = ThrObj<TH2F>
		("tofw1", "tofw1", 200, 15, 68, 200, -0.3, 0.3);
};

void Analyze(ThrHistStruct *ThrHist, std::string part, std::string magf, std::string aux_name, const int proc_num)
{	
	Box box = Box("Parameters of run " + std::to_string(proc_num) + 
		" out of " + std::to_string(Par.part_queue.size()*
		Par.magf_queue.size()*Par.aux_name_queue.size()));

	std::string input_file_name = "../data/" + 
		Par.run_name + "/Single/" + part + magf + aux_name + ".root";
	std::string se_file_name = "../../analysis/data/" + 
		Par.run_name + "/sum" + magf + ".root";

	TFile input_file = TFile(input_file_name.c_str());
	TFile se_file = TFile(se_file_name.c_str());
	
	const double nevents = static_cast<double>(((TTree *) input_file.Get("Tree"))->GetEntries());

	if (nevents <= 0)
	{
		Print("Error: Number of events is equal or less than 0!");
		exit(1);
	}
	
	TH1F *orig_hist = (TH1F *) input_file.Get("orig_pt");
	
	//thredhold is needed since there can be a little noise in the historgram
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

		//this normalization is needed to merge 2 files with flat pt distribution with different ranges
		event_norm *= (Par.ptmax - Par.ptmin)/(up_pt_bound-low_pt_bound);
	}
	
	box.AddEntry("Run name", Par.run_name);
	box.AddEntry("Orig particle", part);
	box.AddEntry("Magnetic field", magf);
	box.AddEntry("Orig pT distribution span", 
		DtoStr(low_pt_bound) + " < pT < " + DtoStr(up_pt_bound));
	box.AddEntry("Use weight function", Par.do_use_weight_func);
	box.AddEntry("Reweight alpha", Par.do_reweight_alpha);
	
	box.AddEntry("Minimum p_T, GeV", Par.ptmin);
	box.AddEntry("Maximum p_T, GeV", Par.ptmax);
	box.AddEntry("Number of events to be analyzed, 1e6", nevents/1e6, 3);

	box.AddEntry("Number of threads", Par.nthreads);
	
	box.Print();
	
	//tsallis weight function
	std::unique_ptr<TF1> weight_func = std::make_unique<TF1>(TF1("weight_func", 
		"[0]*x^([5]-1)*(1 + ([1] - 1)*(sqrt([4] + x^2)^([5]) - [2])/[3])^(-1./([1]-1.))"));
	weight_func->SetParameters(
		ReadFileIntoArray("../input/Spectra/" + Par.system + "/" + part + ".txt", 6));
	
	ROOT::EnableImplicitMT(Par.nthreads);
	ROOT::TTreeProcessorMT tp(input_file_name.c_str());
	
	double ncalls = 0.;

	bool process_finished = false;

	auto ProcessMP = [&](TTreeReader &reader)
	{	
		std::shared_ptr<TH1F> orig = ThrHist->orig.Get();
		std::shared_ptr<TH2F> orig_pt_vs_pt = ThrHist->orig_pt_vs_pt.Get();
		
		std::shared_ptr<TH2F> dce0 = ThrHist->dce0.Get();
		std::shared_ptr<TH2F> dcw0 = ThrHist->dcw0.Get();
		std::shared_ptr<TH2F> dce1 = ThrHist->dce1.Get();
		std::shared_ptr<TH2F> dcw1 = ThrHist->dcw1.Get();

		std::shared_ptr<TH2F> dce0_unsc = ThrHist->dce0_unsc.Get();
		std::shared_ptr<TH2F> dcw0_unsc = ThrHist->dcw0_unsc.Get();
		std::shared_ptr<TH2F> dce1_unsc = ThrHist->dce1_unsc.Get();
		std::shared_ptr<TH2F> dcw1_unsc = ThrHist->dcw1_unsc.Get();

		std::shared_ptr<TH2F> pc1e_z_vs_phi = ThrHist->pc1e_z_vs_phi.Get();
		std::shared_ptr<TH2F> pc1w_z_vs_phi = ThrHist->pc1w_z_vs_phi.Get();
		
		std::shared_ptr<TH2F> pc2_z_vs_phi = ThrHist->pc2_z_vs_phi.Get();
		std::shared_ptr<TH2F> pc3e_z_vs_phi = ThrHist->pc3e_z_vs_phi.Get();
		std::shared_ptr<TH2F> pc3w_z_vs_phi = ThrHist->pc3w_z_vs_phi.Get();
		
		std::array<std::shared_ptr<TH2F>, 4> emcale_pos, emcale_neg, emcalw_pos, emcalw_neg;

		std::shared_ptr<TH1F> strip_tofw_hist = ThrHist->strip_tofw_hist.Get();
		std::shared_ptr<TH1F> slat_hist = ThrHist->slat_hist.Get();
		
		std::shared_ptr<TH2F> tofe0 = ThrHist->tofe0.Get();
		std::shared_ptr<TH2F> tofe1 = ThrHist->tofe1.Get();
		
		std::shared_ptr<TH2F> tofw0 = ThrHist->tofw0.Get();
		std::shared_ptr<TH2F> tofw1 = ThrHist->tofw1.Get();

		for (int i = 0; i < 4; i++)
		{
			emcale_pos[i] = ThrHist->emcale_pos[i].Get();
			emcale_neg[i] = ThrHist->emcale_neg[i].Get();
			emcalw_pos[i] = ThrHist->emcalw_pos[i].Get();
			emcalw_neg[i] = ThrHist->emcalw_neg[i].Get();
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
			
			const double bbcz = T.bbcz();
			if (fabs(bbcz) > 30) continue;

			for(int i = 0; i < T.nch(); i++)
			{
				const double the0 = T.the0(i);
				const double pt = (T.mom(i))*sin(the0);

				if (pt < Par.ptmin || pt > Par.ptmax) continue;
				if (IsQualityCut(T.qual(i))) continue;
				
				int charge = T.charge(i);
				if (charge != -1 && charge != 1) continue;

				const double zed = T.zed(i);
				if (abs(zed) > 75 && abs(zed) < 3) continue;
				
				if (!(fabs(the0)<100 &&
					((bbcz > 0 && ((bbcz - 250*tan(the0 - 3.1416/2)) > 2 ||
					(bbcz - 200*tan(the0 - 3.1416/2)) < -2)) ||
					(bbcz < 0 && ((bbcz - 250*tan(the0 - 3.1416/2))< -2 ||
					(bbcz - 200*tan(the0 - 3.1416/2)) > 2))))) continue;
	
				//end of basic cuts

				const double alpha = T.alpha(i);
				const double phi = T.phi(i);
				double board;
				
				if (phi>1.5) board = ((3.72402-phi+0.008047*cos(phi+0.87851))/0.01963496);
				else board = ((0.573231 + phi - 0.0046 * cos(phi + 0.05721)) / 0.01963496);

				double particle_weight = event_weight;

				if (phi > 1.5)
				{
					if (zed >=0) 
					{
						dce0_unsc->Fill(board, alpha, event_weight);
						if (Par.do_reweight_alpha) particle_weight *= 
							Par.dce0_alpha_reweight->GetBinContent(
							Par.dce0_alpha_reweight->FindBin(alpha));
						dce0->Fill(board, alpha, particle_weight);
					}
					else 
					{
						dce1_unsc->Fill(board, alpha, event_weight);
						if (Par.do_reweight_alpha) particle_weight *= 
							Par.dce1_alpha_reweight->GetBinContent(
							Par.dce1_alpha_reweight->FindBin(alpha));
						dce1->Fill(board, alpha, particle_weight);
					}
				}
				else
				{
					if (zed >=0) 
					{
						dcw0_unsc->Fill(board, alpha, event_weight);
						if (Par.do_reweight_alpha) particle_weight *= 
							Par.dcw0_alpha_reweight->GetBinContent(
							Par.dcw0_alpha_reweight->FindBin(alpha));
						dcw0->Fill(board, alpha, particle_weight);
					}
					else 
					{
						dcw1_unsc->Fill(board, alpha, event_weight);
						if (Par.do_reweight_alpha) particle_weight *= 
							Par.dcw1_alpha_reweight->GetBinContent(
							Par.dcw1_alpha_reweight->FindBin(alpha));
						dcw1->Fill(board, alpha, particle_weight);
					}
				}
				
				if (IsDeadDC(phi, zed, board, alpha)) continue;

				double pc1phi = atan2(T.ppc1y(i), T.ppc1x(i));
				if (phi >= 1.5 && pc1phi < 0) pc1phi += M_PI*2.;
				
				if (phi < 1.5) pc1w_z_vs_phi->Fill(T.ppc1z(i), pc1phi, particle_weight);
				else pc1e_z_vs_phi->Fill(T.ppc1z(i), pc1phi, particle_weight);
				
				if (IsDeadPC1(phi, T.ppc1z(i), pc1phi)) continue;

				orig_pt_vs_pt->Fill(orig_pt, pt, event_weight);
				
				if (IsMatch(T.pc2sdz(i), T.pc2sdphi(i), 2., 2.))
				{
					const double pc2z = T.ppc2z(i) - T.pc2dz(i);
					const double pc2phi = atan2(T.ppc2y(i), T.ppc2x(i)) - T.pc2dphi(i);
					
					pc2_z_vs_phi->Fill(pc2z, pc2phi, particle_weight);
				}

				if (IsMatch(T.pc3sdz(i), T.pc3sdphi(i), 2., 2.))
				{
					const double pc3z = T.ppc3z(i) - T.pc3dz(i);
					double pc3phi = atan2(T.ppc3y(i), T.ppc3x(i) - T.pc3dphi(i));
					
					if (phi > 1.5) 
					{
						if (pc3phi < 0) pc3phi += 6.2831853;
						pc3e_z_vs_phi->Fill(pc3z, pc3phi, particle_weight);
					}
					else pc3w_z_vs_phi->Fill(pc3z, pc3phi, particle_weight);
				}
				
				if (IsHit(T.tofdz(i)))
				{
					const double beta = T.pltof(i)/T.ttof(i)/29.97;
					const double eloss = 0.0014*pow(beta, -1.66);

					slat_hist->Fill(T.slat(i), particle_weight);
					
					if (IsMatch(T.tofsdz(i), T.tofsdphi(i), 2., 2.) && 
						!IsBadSlat(T.slat(i)) &&
						T.etof(i) > eloss) 
					{
						if (zed >= 0) 
						{
							tofe0->Fill(T.ptofy(i), T.ptofz(i), T.etof(i)*particle_weight);
						}
						else 
						{
							tofe1->Fill(T.ptofy(i), T.ptofz(i), T.etof(i)*particle_weight);
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
						strip_tofw_hist->Fill(T.striptofw(i), particle_weight*0.877);
						if (!IsBadStripTOFw(static_cast<int>(T.striptofw(i))))
						{
							if (zed >=0) 
							{
								tofw0->Fill(board, alpha, particle_weight*0.877);
							}
							else 
							{
								tofw1->Fill(board, alpha, particle_weight*0.877);
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
							if (phi > 1.5) 
							{
								if (zed >= 0) emcale_pos[T.sect(i)]->Fill(
									T.pemcy(i), T.pemcz(i), T.ecore(i)*particle_weight);
								else emcale_neg[T.sect(i)]->Fill(
									T.pemcy(i), T.pemcz(i), T.ecore(i)*particle_weight);
							}
							else
							{	
								if (zed >= 0) emcalw_pos[T.sect(i)]->Fill(
									T.pemcy(i), T.pemcz(i), T.ecore(i)*particle_weight);
								else emcalw_neg[T.sect(i)]->Fill(
									T.pemcy(i), T.pemcz(i), T.ecore(i)*particle_weight);
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

//for ROOT CINT call
void HeatMapper()
{
	for (std::string magf : Par.magf_queue)
	{
		CheckInputFile("../../analysis/data/" + Par.run_name + "/sum" + magf + ".root");
		
		for (std::string part : Par.part_queue)
		{
			for (std::string aux_name : Par.aux_name_queue)
			{
				CheckInputFile("../data/" + Par.run_name + 
					"/Single/" + part + magf + aux_name + ".root");
			}
		}
	}
	if (Par.do_use_weight_func)
	{
		for (std::string part : Par.part_queue)
		{
			CheckInputFile("../input/Spectra/" + Par.system + "/" + part + ".txt");
		}
	}
	
	if (Par.do_reweight_alpha)
	{
		std::string data_input_file_name = "../../analysis/data/" + Par.run_name + "/sum.root";
		std::string sim_input_file_name = "../../analysis/data/phenix_sim/" + 
			Par.run_name + "/heatmaps.root";
		
		CheckInputFile(data_input_file_name);
		if (CheckInputFile(sim_input_file_name, false)) 
		{
			Par.data_input_file = std::unique_ptr<TFile>(TFile::Open(data_input_file_name.c_str()));
			Par.sim_input_file = std::unique_ptr<TFile>(TFile::Open(sim_input_file_name.c_str()));

			TH2F *data_dceast0 = (TH2F *) Par.data_input_file->Get("dceast0");
			TH2F *data_dceast1 = (TH2F *) Par.data_input_file->Get("dceast1");
			TH2F *data_dcwest0 = (TH2F *) Par.data_input_file->Get("dcwest0");
			TH2F *data_dcwest1 = (TH2F *) Par.data_input_file->Get("dcwest1");

			TH2F *sim_dceast0 = (TH2F *) Par.sim_input_file->Get("unscaled_dceast0");
			TH2F *sim_dceast1 = (TH2F *) Par.sim_input_file->Get("unscaled_dceast1");
			TH2F *sim_dcwest0 = (TH2F *) Par.sim_input_file->Get("unscaled_dcwest0");
			TH2F *sim_dcwest1 = (TH2F *) Par.sim_input_file->Get("unscaled_dcwest1");

			CutDCDeadAreas(data_dceast0, &IsDeadDC, 2., 1.);
			CutDCDeadAreas(data_dceast1, &IsDeadDC, 2., -1.);
			CutDCDeadAreas(data_dcwest0, &IsDeadDC, 1., 1.);
			CutDCDeadAreas(data_dcwest1, &IsDeadDC, 1., -1.);

			CutDCDeadAreas(sim_dceast0, &IsDeadDC, 2., 1.);
			CutDCDeadAreas(sim_dceast1, &IsDeadDC, 2., -1.);
			CutDCDeadAreas(sim_dcwest0, &IsDeadDC, 1., 1.);
			CutDCDeadAreas(sim_dcwest1, &IsDeadDC, 1., -1.);

			Par.dce0_alpha_reweight = (TH1F *)
				data_dceast0->ProjectionY("dce0_reweight",
				1, data_dceast0->GetXaxis()->GetNbins())->Clone();
			Par.dce1_alpha_reweight = (TH1F *) 
				data_dceast1->ProjectionY("dce1_reweight",
				1, data_dceast1->GetXaxis()->GetNbins())->Clone();
			Par.dcw0_alpha_reweight = (TH1F *) 
				data_dcwest0->ProjectionY("dcw0_reweight",
				1, data_dcwest0->GetXaxis()->GetNbins())->Clone();
			Par.dcw1_alpha_reweight = (TH1F *) 
				data_dcwest1->ProjectionY("dcw1_reweight",
				1, data_dcwest1->GetXaxis()->GetNbins())->Clone();

			Par.dce0_alpha_reweight->Scale(sim_dceast0->Integral()/Par.dce0_alpha_reweight->Integral());
			Par.dce1_alpha_reweight->Scale(sim_dceast1->Integral()/Par.dce1_alpha_reweight->Integral());
			Par.dcw0_alpha_reweight->Scale(sim_dcwest0->Integral()/Par.dcw0_alpha_reweight->Integral());
			Par.dcw1_alpha_reweight->Scale(sim_dcwest1->Integral()/Par.dcw1_alpha_reweight->Integral());

			Par.dce0_alpha_reweight->Divide(sim_dceast0->ProjectionY("dce0_proj",
				1, sim_dceast0->GetXaxis()->GetNbins()));
			Par.dce1_alpha_reweight->Divide(sim_dceast1->ProjectionY("dce1_proj",
				1, sim_dceast1->GetXaxis()->GetNbins()));
			Par.dcw0_alpha_reweight->Divide(sim_dcwest0->ProjectionY("dcw0_proj",
				1, sim_dcwest0->GetXaxis()->GetNbins()));
			Par.dcw1_alpha_reweight->Divide(sim_dcwest1->ProjectionY("dcw1_proj",
				1, sim_dcwest1->GetXaxis()->GetNbins()));

			//capping alpha reweight since experiments extends to higher pt than simulation
			//capped values do not affect the simulation since they are statisticaly insufficient
			//this cap is only needed to exclude very big values in some points in the DC map
			//for better visibility
			for (int i = 0; i < Par.dce0_alpha_reweight->GetXaxis()->GetNbins(); i++)
			{
				if (Par.dce0_alpha_reweight->GetBinContent(i) > 100.)
				{
					Par.dce0_alpha_reweight->SetBinContent(i, 100.);
				}
			}
			for (int i = 0; i < Par.dce1_alpha_reweight->GetXaxis()->GetNbins(); i++)
			{
				if (Par.dce1_alpha_reweight->GetBinContent(i) > 100.)
				{
					Par.dce1_alpha_reweight->SetBinContent(i, 100.);
				}
			}
			for (int i = 0; i < Par.dcw0_alpha_reweight->GetXaxis()->GetNbins(); i++)
			{
				if (Par.dcw0_alpha_reweight->GetBinContent(i) > 100.)
				{
					Par.dcw0_alpha_reweight->SetBinContent(i, 100.);
				}
			}
			for (int i = 0; i < Par.dcw1_alpha_reweight->GetXaxis()->GetNbins(); i++)
			{
				if (Par.dcw1_alpha_reweight->GetBinContent(i) > 100.)
				{
					Par.dcw1_alpha_reweight->SetBinContent(i, 100.);
				}
			}

			std::string alpha_reweight_output_name = "../../analysis/data/phenix_sim/" + 
				Par.run_name + "/alpha_reweight.root";
			TFile alpha_reweight_output = TFile(alpha_reweight_output_name.c_str(), "RECREATE");

			alpha_reweight_output.cd();
			
			Par.dce0_alpha_reweight->Write();
			Par.dce1_alpha_reweight->Write();
			Par.dcw0_alpha_reweight->Write();
			Par.dcw1_alpha_reweight->Write();

			alpha_reweight_output.Close();
			PrintInfo("File " + alpha_reweight_output_name + " was written");
		}
	}
	
	PrintInfo("Clearing output directory: ../../analysis/data/phenix_sim/" + 
		Par.run_name + "/heatmaps/");
	system(("mkdir -p ../../analysis/data/phenix_sim/" + Par.run_name + "/heatmaps").c_str());
	system(("rm -r ../../analysis/data/phenix_sim/" + Par.run_name + "/heatmaps/*").c_str());
	
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
			//wiriting the result
			const std::string output_file_name = 
				"../../analysis/data/phenix_sim/" + Par.run_name + 
				"/heatmaps/" + part + magf + ".root";
			
			TFile outfile(output_file_name.c_str(), "RECREATE");
			outfile.cd();
			ThrObjHolder.Write();
			outfile.Close();
			PrintInfo("File " + output_file_name + " was written");
		}
	}

	PrintInfo("Merging output files into one");

	system(("hadd -f ../../analysis/data/phenix_sim/" + 
		Par.run_name + "/heatmaps.root ../../analysis/data/phenix_sim/" + 
		Par.run_name + "/heatmaps/*.root").c_str());
}

//main calls the same function CINT would call
int main()
{
	HeatMapper();
	return 0;
}
