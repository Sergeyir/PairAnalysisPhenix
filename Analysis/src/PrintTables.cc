#include "../lib/ErrorHandler.h"
#include "../lib/Tools.h"
#include "../lib/InputTool.h"
#include "../lib/OutputTool.h"
#include "../lib/ParInv.h"

void PrintTables()
{
	std::string spectra_input_dir = "../data/Spectra/" + Par.run + "/" + Par.runnum + "/";
	std::vector<std::string> method_names = {"EMC2PID", "TOF2PID", "2PID", "1PID", "noPID"};

	for (int i = 0; i < Par.CType.size; i++)
	{
		std::vector<std::vector<double>> pt, val, val_stat, val_sys;
		
		pt.resize(method_names.size());
		val.resize(method_names.size());
		val_stat.resize(method_names.size());
		val_sys.resize(method_names.size());
		
		for (int j = 0; j < method_names.size(); j++)
		{
			pt[j].resize(DefaultPt.ptmin.size());
			val[j].resize(DefaultPt.ptmin.size());
			val_stat[j].resize(DefaultPt.ptmin.size());
			val_sys[j].resize(DefaultPt.ptmin.size());

			for (int ptc = 0; ptc < DefaultPt.ptmin.size(); ptc++)
			{
				pt[j][ptc] = -9999.;
				val[j][ptc] = -9999.;
				val_stat[j][ptc] = -9999.;
				val_sys[j][ptc] = -9999.;
			}
			
			std::string input_file_name = spectra_input_dir + method_names[j] + "_" + 
				Par.particle.name_nl + "_" + Par.CType.cname_nop[i] + ".txt";
			
			CheckInputFile(input_file_name);
			ifstream input_file(input_file_name);

			int ptc = 0;
			while (input_file >> pt[j][ptc] >> val[j][ptc] >> val_stat[j][ptc] >> val_sys[j][ptc])
			{
				ptc++;
			}
		}
		
		system(("mkdir -p ../output/Tables/Spectra/" + Par.run + "/" + Par.runnum).c_str());
		std::string output_file_name = "../output/Tables/Spectra/" + Par.run + "/" + Par.runnum +
			"/SYS_comp_" + Par.particle.name_nl + "_" + Par.CType.cname_nop[i] + ".tex";
		
		ofstream output_file(output_file_name);

		
		output_file << "\\begin{table}[h!]" << std::endl;
		output_file << "\\centering" << std::endl;
		
		output_file << "\\begin{tabular}{c|";
		
		for (int j = 0; j < method_names.size(); j++)
		{
			output_file << "|c";
		}
		output_file << "}" << std::endl;

		output_file << "	$p_T$, GeV/c ";
		
		for (int j = 0; j < method_names.size(); j++)
		{
			output_file << " & " + method_names[j];
		}

		output_file << " \\\\" << std::endl;

		output_file << "	\\hline" << std::endl;
		output_file << "	\\hline" << std::endl;
		
		for (int ptc = 0; ptc < DefaultPt.ptmin.size(); ptc++)
		{
			output_file << "	" << DefaultPt.ptmin[ptc] << " - " << DefaultPt.ptmax[ptc];
			for (int j = 0; j < method_names.size(); j++)
			{
				bool printed = false;
				for (int k = 0; k < pt[j].size(); k++)
				{
					if (abs(pt[j][k] - Average(DefaultPt.ptmin[ptc], DefaultPt.ptmax[ptc])) < 0.05)
					{
						output_file << " & " << ErrPropagation(val_stat[j][k], val_sys[j][k]);
						printed = true;
					}
				}

				if (!printed)
				{
					output_file << " & -";
				}
			}
			output_file << " \\\\" << std::endl;
		}

		output_file << "\\end{tabular}";
		
		output_file << "	\\caption{Propagated statistical and systematic uncertainties (relative, \\%) for different pair mixing methods in " + Par.CType.cname_latex[i] + " centrality class}" << std::endl;
		
		output_file << "\\end{table}";
		
		output_file.close();
		PrintInfo(output_file_name + " was written");
	}
}
