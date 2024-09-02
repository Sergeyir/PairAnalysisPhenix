#include "../lib/ParWidth.h"

#include "../lib/EffTreeReader.h"

#include "PBar.hpp"

#include "../lib/ThrObj.h"

#include "../lib/SingleAnalysis.h"
#include "../lib/PairAnalysis.h"

#include "../lib/Tools.h"
#include "../lib/Box.h"
#include "../lib/InputTools.h"

#include "ROOT/TTreeProcessorMT.hxx"

struct ThrHistStruct
{
	ThrObj<TH2F> npart_hist = ThrObj<TH2F>("npart", "npart", 
		50, 0, 50, 100, Par.ptmin, Par.ptmax);

	ThrObj<TH1F> orig = ThrObj<TH1F>("orig","orig", 
		Par.pt_nbins, Par.ptmin, Par.ptmax);
	
	ThrObj<TH2F> orig_pt_vs_pt = ThrObj<TH2F>("orig_pt_vs_pt", "pt", 
		Par.pt_nbins, Par.ptmin, Par.ptmax, Par.pt_nbins, Par.ptmin, Par.ptmax);

	ThrObj<TH2F> inv_mass_hist = ThrObj<TH2F>("InvM", "M_{inv}", 
		100, 0., 10., 2000., 0., 2.);
};

//struct for storing charged particles (positive and negative are separated)
//to combine them to the pairs later
struct PartArray
{
	unsigned int size = 0;
	double mom[3];
	std::array<int, 50> index;
	std::array<double, 50> weight;

	double mass;
};

void Analyze(std::string magf, std::string daughter1, std::string daughter2, const int proc_num)
{
	Box box = Box("Parameters of run " + std::to_string(proc_num) + 
		" out of " + std::to_string(Par.daughter1_queue.size()*Par.magf_queue.size()));

	const std::string channel = 
		ParticleProperties.sname[ParticleProperties.iter_map[daughter1]] + 
		ParticleProperties.sname[ParticleProperties.iter_map[daughter2]];
	
	const std::string input_file_name = "../data/" + Par.run_name + 
		"/Widthless/" + Par.part + "_" + channel + magf + ".root";
	
	TFile input_file = TFile(input_file_name.c_str());
	
	const double nevents = static_cast<double>(((TTree *) input_file.Get("Tree"))->GetEntries());

	if (nevents <= 0)
	{
		Print("Error: Number of events is equal or less than 0!");
		exit(1);
	}

	double event_norm = 1.;
	if (Par.do_use_weight_func)
	{
		TH1F *orig_hist = (TH1F *) input_file.Get("orig_pt");
		
		event_norm = orig_hist->Integral(
			orig_hist->GetXaxis()->FindBin(Par.ptmin),
			orig_hist->GetXaxis()->FindBin(Par.ptmax));
	}

	box.AddEntry("Run name", Par.run_name);
	box.AddEntry("Particle", Par.part);
	box.AddEntry("Decay channel", channel);
	box.AddEntry("Magnetic field", magf);
	box.AddEntry("Use weight function", Par.do_use_weight_func);
	
	box.AddEntry("Minimum p_T, GeV", Par.ptmin);
	box.AddEntry("Maximum p_T, GeV", Par.ptmax);
	box.AddEntry("Number of events to be analyzed, 1e6", nevents/1e6, 3);

	box.AddEntry("Number of threads", Par.nthreads);
	
	box.Print();

	TF1 weight_func = TF1("weight_func", 
		"[0]*x^([5]-1)*(1 + ([1] - 1)*(sqrt([4] + x^2)^([5]) - [2])/[3])^(-1./([1]-1.))");
	weight_func.SetParameters(
		ReadFileIntoArray("../input/Spectra/" + Par.system + "/" + Par.part + ".txt", 6));

	ROOT::EnableImplicitMT(Par.nthreads);
	ROOT::TTreeProcessorMT tp(input_file_name.c_str());
	
	double ncalls = 0;

	ThrHistStruct ThrHist;
	
	bool process_finished = false;
	auto ProcessMP = [&](TTreeReader &reader)
	{	
		std::shared_ptr<TH2> npart_hist = ThrHist.npart_hist.Get();
		std::shared_ptr<TH1F> orig = ThrHist.orig.Get();
		std::shared_ptr<TH2F> orig_pt_vs_pt = ThrHist.orig_pt_vs_pt.Get();
		
		std::shared_ptr<TH2F> inv_mass_hist = ThrHist.inv_mass_hist.Get();

		EffTreeReader T(reader);

		//arrays for positives and negative tracks
		PartArray ptrack, ntrack;

		ptrack.mass = ParticleProperties.mass[ParticleProperties.iter_map[daughter1]];
		ntrack.mass = ParticleProperties.mass[ParticleProperties.iter_map[daughter2]];
		
		while (reader.Next())
		{
			ncalls += 1.;
			const double orig_pt = sqrt(pow(T.mom_orig(0), 2) + pow(T.mom_orig(1), 2));
			
			double event_weight = 1.;
			if (Par.do_use_weight_func) event_weight = weight_func.Eval(orig_pt)/event_norm;
			
			orig->Fill(orig_pt, event_weight);
			
			npart_hist->Fill(T.nch()-0.5, orig_pt, event_weight);
			
			const double bbcz = T.bbcz();
			if (fabs(bbcz) > 30) continue;

			if (T.nch() <= 0 || T.nch() > 49) continue;

			ptrack.size = 0;
			ntrack.size = 0;
	
			for(int i = 0; i < T.nch(); i++)
			{
				const double the0 = T.the0(i);
				const double pt = (T.mom(i))*sin(the0);
				
				if (pt < Par.ptmin || pt > Par.ptmax) continue;
				
				const int charge = T.charge(i);
				if (charge != 1 && charge != -1) continue;
					
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

				double particle_weight = 1.;;
				
				double board;
				if (phi>1.5) board = ((3.72402-phi+0.008047*cos(phi+0.87851))/0.01963496);
				else board = ((0.573231 + phi - 0.0046 * cos(phi + 0.05721)) / 0.01963496);

				if (IsDeadDC(phi, zed, board, alpha)) continue;

				double pc1phi = atan2(T.ppc1y(i), T.ppc1x(i));
				if (phi >= 1.5 && pc1phi < 0) pc1phi += M_PI*2.;

				if (IsDeadPC1(phi, T.ppc1z(i), pc1phi)) continue;
				
				if (charge == 1) 
				{
					ptrack.index[ptrack.size] = i;
					ptrack.weight[ptrack.size] = particle_weight;
					ptrack.size++;
				}
				else
				{
					ntrack.index[ntrack.size] = i;
					ntrack.weight[ntrack.size] = particle_weight;
					ntrack.size++;
				}
			}
			
			//combining particles into pairs
			for (unsigned int i = 0; i < ptrack.size; i++)
			{
				for (unsigned int j = 0; j < ntrack.size; j++)
				{
					const int ppc = ptrack.index[i]; //positive particle counter
					const int npc = ntrack.index[j]; //negative particle counter
					
					ptrack.mom[0] = T.mom(ppc)*sin(T.the0(ppc))*cos(T.phi0(ppc));
					ptrack.mom[1] = T.mom(ppc)*sin(T.the0(ppc))*sin(T.phi0(ppc));
					ptrack.mom[2] = T.mom(ppc)*cos(T.the0(ppc));
					
					ntrack.mom[0]	= T.mom(npc)*sin(T.the0(npc))*cos(T.phi0(npc));
					ntrack.mom[1]	= T.mom(npc)*sin(T.the0(npc))*sin(T.phi0(npc));
					ntrack.mom[2]	= T.mom(npc)*cos(T.the0(npc));

					if (IsGhostCut(T.zed(ppc) - T.zed(npc), 
						T.alpha(ppc) - T.alpha(npc), 
						T.phi(ppc) - T.phi(npc))) continue;
						
					if (IsOneArmCut(T.phi(ppc), T.phi(npc))) continue;
					
					const double mass = GetMass(ptrack.mom, ntrack.mom, ptrack.mass, ntrack.mass);

					const double pt = sqrt((ptrack.mom[0] + ntrack.mom[0])*(ptrack.mom[0] + ntrack.mom[0]) +
						(ptrack.mom[1] + ptrack.mom[1])*(ntrack.mom[1] + ntrack.mom[1]));

					orig_pt_vs_pt->Fill(pt, orig_pt, event_weight*ptrack.weight[i]*ntrack.weight[j]);
					inv_mass_hist->Fill(pt, mass, event_weight*ptrack.weight[i]*ntrack.weight[j]);
				}
			}
		}
	};

	auto PbarCall = [&]()
	{
		ProgressBar pbar = ProgressBar("Fancy2");
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

	std::string file_name = "../../analysis/data/phenix_sim/" + 
		Par.run_name + "/Widthless/" + Par.part + magf + ".root";
	
	TFile outfile(file_name.c_str(), "RECREATE");

	ThrObjHolder.Write();

	outfile.Close();
	
	PrintInfo("File " + file_name + " was written");
}

void WidthAnalyze()
{
	int nchannels = Minimum(Par.daughter1_queue.size(), Par.daughter2_queue.size());
	{
		for (std::string magf : Par.magf_queue)
		{
			for (int i = 0; i < nchannels; i++)
			{
				const std::string channel = 
					ParticleProperties.sname[ParticleProperties.iter_map[Par.daughter1_queue[i]]] + 
					ParticleProperties.sname[ParticleProperties.iter_map[Par.daughter2_queue[i]]];
				
				CheckInputFile("../data/" + Par.run_name +  
					"/Widthless/" + Par.part + "_" + channel + magf + ".root");
			}
		}
	}

	if (Par.do_use_weight_func) CheckInputFile("../input/Spectra/" + Par.system + "/" + Par.part + ".txt");
	
	system(("mkdir -p ../../analysis/data/phenix_sim/" + Par.run_name + "/Widthless").c_str());	
	system(("rm -r ../../analysis/data/phenix_sim/" + Par.run_name + "Widthless/*").c_str());
	
	int num = 1;
	for (std::string magf : Par.magf_queue)
	{
		for (int i = 0; i < nchannels; i++)
		{
			Analyze(magf, Par.daughter1_queue[i], Par.daughter2_queue[i], num);
			num++;
		}
	}

	system(("hadd -f ../../analysis/data/phenix_sim/" + Par.run_name + 
		"/Widthless_" + Par.part + ".root ../../analysis/data/phenix_sim/" +
		Par.run_name + "/Widthless/*.root").c_str());
}

int main()
{
   WidthAnalyze();  
   return 0;
}
