#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <cmath>

using namespace std;

class HepMCReader
{
	private:

		ifstream input_file;
		string input_name;

		double Nlines;
		bool is_end{0};
		unsigned long nEvent{0};

		//event properties
		int ev_num, ev_multiparton_num;
		double ev_scale;

		//particles properties vectors
		vector<int> p_id, p_status;
		vector<double> p_px, p_py, p_pz, p_e, p_pt, p_m, p_pol_theta, p_pol_phi;

		//cross section in pb
		double sigma{0}, sigma_err{0};
		double sigma_sum{0}, sigma_sum_err{0};

		void CountFileNlines(string input_file_name)
		{

			char line[200];
			ifstream input_file_Nlines(input_file_name.c_str());

			while(input_file_Nlines.getline(line, 200)) {Nlines++;}

			cout << "Input file has " << Nlines << endl;
		}

		void PrepareForNewEvent()
		{
			char junk[200];

			//if (nEvent % 1000 == 0) cout << "Events passed " << nEvent << endl;

			double cross_section, cross_section_err;

			char input_type;

			input_file >> input_type >> ev_num >> ev_multiparton_num >> ev_scale;

			input_file.getline(junk, 200);
			input_file.getline(junk, 200); //N
			input_file.getline(junk, 200); //U

			input_file >> input_type >> cross_section >> cross_section_err;

			sigma = cross_section;
			sigma_err = cross_section_err;

			sigma_sum += cross_section;
			sigma_sum_err += cross_section_err;

			input_file.getline(junk, 200); //U
		}

		void ClearPreviousEvent()
		{
			p_id.clear();
			p_status.clear();
			
			p_px.clear();
			p_py.clear();
			p_pz.clear();
			p_e.clear();

			p_m.clear();

			p_pol_theta.clear();
			p_pol_phi.clear();

		}

	public:
	
		HepMCReader(string input_file_name)
		{
			input_file.open(input_file_name.c_str());	

			input_name = input_file_name;
			
			if (!input_file.is_open())
			{
				printf("\033[31mError: ");
				printf("\033[37m");
				cout << "File " << input_file_name << " not found" << endl;
				exit(0);
			}
	
			//CountFileNlines(input_file_name);
			
			char junk[200];
			
			input_file.getline(junk, 200); //skip line
			input_file.getline(junk, 200);
			input_file.getline(junk, 200);

			PrepareForNewEvent();

			input_file.getline(junk, 200);
		}
		
		void ReadNextEvent()
		{
			ClearPreviousEvent();

			bool skip = 0;
			char input_type;
			char junk[200];
	
			input_file.getline(junk, 200);

			nEvent++;
	
			while(!skip && !input_file.eof())
			{
				input_file >> input_type;
				if (input_type == 'P')
				{
					// Reading particles input

					int in_barcode, in_id, in_status;
					double in_px, in_py, in_pz, in_e, in_m, in_pol_theta, in_pol_phi;

					input_file >> in_barcode >> in_id >> in_px >> in_py >> in_pz >> in_e >> in_m >> in_status >> in_pol_theta >> in_pol_phi;

					input_file.getline(junk, 200);
					
					p_id.push_back(in_id);
					p_status.push_back(in_status);

					p_px.push_back(in_px);
					p_py.push_back(in_py);
					p_pz.push_back(in_pz);
					p_e.push_back(in_e);

					p_m.push_back(in_m);

					p_pol_theta.push_back(in_pol_theta);
					p_pol_phi.push_back(in_pol_phi);
				}
			
				else if (input_type == 'E')
				{
					PrepareForNewEvent();

					skip = 1;
					continue;
				}
				
				else
				{
					char junk[200];
					input_file.getline(junk, 200);
				}
			}

			if (input_file.eof())
			{
				is_end = 1;
			}
		}

		bool isStable(int ipart)
		{
			if (p_status[ipart] == 1) return 1;
			return 0;
		}

		int id(int ipart)
		{
		return p_id[ipart];
		}

		unsigned long eventSize()
		{
			return p_px.size();
		}

		double px(int ipart)
		{
			return p_px[ipart];
		}
			
		double py(int ipart)
		{
			return p_py[ipart];
		}
			
		double pz(int ipart)
		{
			return p_pz[ipart];
		}
			
		double e(int ipart)
		{
			return p_e[ipart];
		}
			
		double pT(int ipart)
		{
			double pt = sqrt(pow(p_px[ipart], 2) + pow(p_py[ipart], 2));
			return pt;
		}

		double sigmaGen()
		{
			return sigma_sum/nEvent;
		}

		double sigmaGenErr()
		{
			return sqrt(sigma_sum_err/nEvent);
		}

		float weight()
		{
			return 1;
		}

		bool isEnd()
		{
			return is_end;
		}
};
