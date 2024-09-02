#include "ErrorHandler.h"

struct
{
	std::vector<TH1F *> hists;
	std::vector<std::string> hists_x_val;

	std::vector<std::string> output_files_name;
	std::vector<ofstream> output_files;
	
	std::vector<std::vector<TF1 *>> fits, bg_fits;
	std::vector<std::vector<int>> bg_fit_par_min;

	TLatex fits_name_text = TLatex(0.2, 0.2, "");
	TLatex chi2ndf_text = TLatex(0.8, 0.9, "");
	
	//graph that contains points for background tweaking
	TGraph *bg_points_gr = new TGraph();
	
	int current_hist = -1;
	int current_fit_mode = 0;
	
	bool fit_point_active = false;
	int current_active_point = -1;
	TGraph *current_active_point_gr = new TGraph();

	Color_t active_color = kRed+1;

	void AddFitType(std::string output_file_name)
	{
		system(("cp " + output_file_name + " " + output_file_name + ".backup").c_str());
		CheckOutputFile(output_file_name, std::ios::app);
		output_files_name.push_back(output_file_name);
		output_files.push_back(ofstream("../../tmp/a.txt"));
	}
	
	void AddHist(TH1F *data_hist, std::string x_val)
	{
		fits.resize(fits.size()+1);
		bg_fits.resize(bg_fits.size()+1);
		bg_fit_par_min.resize(bg_fit_par_min.size()+1);
		
		current_hist = 0;
		hists.push_back((TH1F *) data_hist->Clone());
		hists_x_val.push_back(x_val);
	}
	
	void AddFit(TF1 *fit, TF1 *bg_fit, const int bg_par_min_start, const unsigned int hist_num)
	{
		fits[hist_num].push_back((TF1 *) fit->Clone());
		bg_fits[hist_num].push_back((TF1 *) bg_fit->Clone());
		bg_fit_par_min[hist_num].push_back(bg_par_min_start);

		for (int i = 0; i < bg_fit->GetNpar(); i++)
		{
			bg_fits[hist_num].back()->ReleaseParameter(i);
			fits[hist_num].back()->ReleaseParameter(bg_par_min_start+i);
		}
	}

	void AddFit(TF1 *fit, TF1 *bg_fit, Double_t (*fcn)(Double_t *, Double_t *), const int bg_par_min_start, const unsigned int hist_num)
	{
		double x_min, x_max;
		fit->GetRange(x_min, x_max);
		
		fits[hist_num].push_back(new TF1(fit->GetName(), fcn, x_min, x_max, fit->GetNpar()));
		bg_fits[hist_num].push_back((TF1 *) bg_fit->Clone());
		bg_fit_par_min[hist_num].push_back(bg_par_min_start);


		for (int i = 0; i < fit->GetNpar(); i++)
		{
			fits[hist_num].back()->SetParameter(i, fit->GetParameter(i));
			
			double par_min, par_max;
			fit->GetParLimits(i, par_min, par_max);
			fits[hist_num].back()->SetParLimits(i, par_min, par_max);
		}

		for (int i = 0; i < bg_fit->GetNpar(); i++)
		{
			bg_fits[hist_num].back()->ReleaseParameter(i);
			fits[hist_num].back()->ReleaseParameter(bg_par_min_start+i);
		}
	}

	//Adding fit for the last hist added
	void AddFit(TF1 *fit, TF1 *bg_fit, const int bg_par_min_start)
	{
		AddFit(fit, bg_fit, bg_par_min_start, hists.size()-1);
	}

	//Adding fit for the last hist added
	void AddFit(TF1 *fit, TF1 *bg_fit, Double_t (*fcn)(Double_t *, Double_t *), const int bg_par_min_start)
	{
		AddFit(fit, bg_fit, fcn, bg_par_min_start, hists.size()-1);
	}

	void SetBgPoints()
	{
		double x_min, x_max;
		fits[current_hist][current_fit_mode-1]->GetRange(x_min, x_max);
		
		//resetting bg points graph
		delete bg_points_gr;
		bg_points_gr = new TGraph();
		
		bg_points_gr->SetEditable(kFALSE);
		
		bg_points_gr->SetMarkerStyle(89);
		bg_points_gr->SetMarkerSize(3.);
		bg_points_gr->SetMarkerColor(active_color);
		
		const double x_shift = (x_max - x_min)/
			(static_cast<double>(bg_fits[current_hist][current_fit_mode-1]->GetNpar())-1.);
			
		bg_points_gr->AddPoint(x_min + 1e-15, 
			bg_fits[current_hist][current_fit_mode-1]->Eval(x_min));
		for (int j = 1; j < bg_fits[current_hist][current_fit_mode-1]->GetNpar()-1; j++)
		{
			const double x_pos = x_min+x_shift*static_cast<double>(j);
			bg_points_gr->AddPoint(x_pos, bg_fits[current_hist][current_fit_mode-1]->Eval(x_pos));
		 }
		bg_points_gr->AddPoint(x_max - 1e-15, 
			bg_fits[current_hist][current_fit_mode-1]->Eval(x_max));
	}

	void SetActiveFit()
	{
		bool active_set_check = false;
		
		for (int i = 0; i < fits[current_hist].size(); i++)
		{
			if (i == current_fit_mode-1)
			{
				active_set_check = true;
				
				fits[current_hist][i]->SetLineColor(active_color);
				bg_fits[current_hist][i]->SetLineColor(active_color);
				fits[current_hist][i]->SetLineWidth(4);
				bg_fits[current_hist][i]->SetLineWidth(4);	

				SetBgPoints();
			}
			else
			{
				fits[current_hist][i]->SetLineColor(kGray+1);
				bg_fits[current_hist][i]->SetLineColor(kGray+1);
				fits[current_hist][i]->SetLineWidth(1);
				bg_fits[current_hist][i]->SetLineWidth(1);
			}
		}

		if (!active_set_check)
		{
			delete bg_points_gr;
			bg_points_gr = new TGraph();
		}
	}

	void PerformFreeFit()
	{
		if (current_fit_mode == 0 || current_fit_mode > fits[current_hist].size()) return;
		
		for (int i = 0; i < bg_fits[current_hist].back()->GetNpar(); i++)
		{
			fits[current_hist][current_fit_mode-1]->ReleaseParameter(
				bg_fit_par_min[current_hist][current_fit_mode-1]+i);
		}
		
		hists[current_hist]->Fit(fits[current_hist][current_fit_mode-1], "RQMNB");
		
		for (int i = 0; i < bg_fits[current_hist].back()->GetNpar(); i++)
		{
			bg_fits[current_hist][current_fit_mode-1]->SetParameter(i, 
				fits[current_hist][current_fit_mode-1]->GetParameter(
				bg_fit_par_min[current_hist][current_fit_mode-1]+i));
		}
		
		SetBgPoints();
	}
	
	void ResetFit()
	{
		if (current_fit_mode == 0 || current_fit_mode > fits[current_hist].size()) return;
		
		double x_min, x_max;
		fits[current_hist][current_fit_mode-1]->GetRange(x_min, x_max);
		
		//resetting bg points graph
		delete bg_points_gr;
		bg_points_gr = new TGraph();
		
		bg_points_gr->SetEditable(kFALSE);
		
		bg_points_gr->SetMarkerStyle(89);
		bg_points_gr->SetMarkerSize(3.);
		bg_points_gr->SetMarkerColor(active_color);
		
		const double x_shift = (x_max - x_min)/
			(static_cast<double>(bg_fits[current_hist][current_fit_mode-1]->GetNpar())-1.);
			
		bg_points_gr->AddPoint(x_min + 1e-15, 
			hists[current_hist]->GetBinContent(hists[current_hist]->GetXaxis()->FindBin(x_min)));
		for (int j = 1; j < bg_fits[current_hist][current_fit_mode-1]->GetNpar()-1; j++)
		{
			const double x_pos = x_min+x_shift*static_cast<double>(j);
			bg_points_gr->AddPoint(x_pos, 
				hists[current_hist]->GetBinContent(hists[current_hist]->GetXaxis()->FindBin(x_pos)));
		}
		bg_points_gr->AddPoint(x_max - 1e-15, 
			hists[current_hist]->GetBinContent(hists[current_hist]->GetXaxis()->FindBin(x_max)));
		
		bg_points_gr->Fit(bg_fits[current_hist][current_fit_mode-1], "RQMN");
		
		for (int i = 0; i < bg_fits[current_hist][current_fit_mode-1]->GetNpar(); i++)
		{
			fits[current_hist][current_fit_mode-1]->FixParameter(
				bg_fit_par_min[current_hist][current_fit_mode-1]+i, 
				bg_fits[current_hist][current_fit_mode-1]->GetParameter(i));
		}
		fits[current_hist][current_fit_mode-1]->Update();
		hists[current_hist]->Fit(
			fits[current_hist][current_fit_mode-1], "RQMN");
	}

	bool IsActivePoint() 
	{
		return fit_point_active;
	}
	
	void ActivatePoint(const int point_number, const double x, const double y) 
	{
		fit_point_active = true; 
		current_active_point = point_number;
		current_active_point_gr->SetPointX(0, x);
		current_active_point_gr->SetPointY(0, y);
	}
	
	void DeactivatePoint() 
	{
		fit_point_active = false; 
		current_active_point = -1;
		current_active_point_gr->SetPointX(0, -1.);
		current_active_point_gr->SetPointY(0, -1.);
	}

	void PrintHelp()
	{
		std::cout << std::endl;

		std::cout << " ╔══════════════════════════════════════════════════════════════╗" << std::endl;
		std::cout << " ║                         General info                         ║" << std::endl;
		std::cout << " ║                                                              ║" << std::endl;
		std::cout << " ║ Move cursor to the canvas to focus on the window. This will  ║" << std::endl; 
		std::cout << " ║ make the keyboard configuration work.                        ║" << std::endl;
		std::cout << " ║                                                              ║" << std::endl;
		std::cout << " ║ Window focus may be removed when switching histograms.       ║" << std::endl;
		std::cout << " ║ Moving a mouse over the canvas will bring it back.           ║" << std::endl;
		std::cout << " ║                                                              ║" << std::endl;
		std::cout << " ║ When you focus on a fit you can grab one of the points to    ║" << std::endl;
		std::cout << " ║ move it vertically. This will change the background fit.     ║" << std::endl;
		std::cout << " ║                                                              ║" << std::endl;
		std::cout << " ║ When changing the approximation all parameters except for    ║" << std::endl;
		std::cout << " ║ parameters for backgroud approximations are printed.         ║" << std::endl;
		std::cout << " ╟──────────────────────────────────────────────────────────────╢" << std::endl;
		std::cout << " ║                    Keyboard configuration                    ║" << std::endl;
		std::cout << " ║                                                              ║" << std::endl;
		std::cout << " ║ h:    Print help (general info and keyboar configuration)    ║" << std::endl;
		std::cout << " ║ j, k: Switch focus between histograms                        ║" << std::endl;
		std::cout << " ║ 1-9:  Switch focus between fits                              ║" << std::endl;
		std::cout << " ║ f:    Perform fit with free backround fit parameters         ║" << std::endl;
		std::cout << " ║ n:    Reset background fit parameters                        ║" << std::endl;
		std::cout << " ║ p:    Write results in specified files                       ║" << std::endl;
		std::cout << " ╚══════════════════════════════════════════════════════════════╝" << std::endl;
		std::cout << std::endl;
	}

	void Draw(bool draw_hist = true, bool first_draw = true)
	{	
		if (first_draw)
		{
			if (current_hist < 0) PrintError("No histograms were added");
			fits_name_text.SetNDC();
			chi2ndf_text.SetNDC();

			if (current_active_point_gr->GetN() == 0) 
			{
				current_active_point_gr->AddPoint(-1., -1.);
				current_active_point_gr->SetMarkerSize(2.5);
				current_active_point_gr->SetMarkerStyle(89);
				current_active_point_gr->SetMarkerColor(active_color);
				current_active_point_gr->SetEditable(kFALSE);
			}

			PrintHelp();
		}
		if (draw_hist)
		{
			DeactivatePoint();
			current_fit_mode = 0;
			
			hists[current_hist]->Draw();
				
			SetActiveFit();	
		}
		
		for (int i = 0; i < fits[current_hist].size(); i++)
		{
			fits[current_hist][i]->Draw("SAME");
			bg_fits[current_hist][i]->Draw("SAME");
		}
		
		if (IsActivePoint()) current_active_point_gr->Draw("P");
		if (bg_points_gr->GetN() > 0) bg_points_gr->Draw("P");

		gPad->GetListOfPrimitives()->Remove(&fits_name_text);
		gPad->GetListOfPrimitives()->Remove(&chi2ndf_text);
		
		if (current_fit_mode > 0 && current_fit_mode <= fits[current_hist].size())
		{
			std::stringstream chi2ndf;
			chi2ndf << std::fixed << std::setprecision(2) << 
				fits[current_hist][current_fit_mode-1]->GetChisquare()/
				fits[current_hist][current_fit_mode-1]->GetNDF();
			
			fits_name_text.SetText(0.2, 0.15, fits[current_hist][current_fit_mode-1]->GetName());
			chi2ndf_text.SetText(0.75, 0.8, ("#chi^2/NDF " + chi2ndf.str()).c_str());
			
			fits_name_text.Draw();
			chi2ndf_text.Draw();
		}
		
		gPad->Modified();
		gPad->Update();
		
		gPad->GetFrame()->SetBit(TBox::kCannotMove);
	}
} GuiFitPar;

void GuiFit()
{
	const int event = gPad->GetEvent();
	const int px = gPad->GetEventX();
	const int py = gPad->GetEventY();
	const double x = gPad->PadtoX(gPad->AbsPixeltoX(px));
	const double y = gPad->PadtoX(gPad->AbsPixeltoY(py));
	
	switch (event)
	{
		case kKeyPress:
			if (px > 47 && px < 58) 
			{
				GuiFitPar.current_fit_mode = px-48;
				GuiFitPar.DeactivatePoint();
				GuiFitPar.SetActiveFit();
			}
			
			switch (px)
			{
				case 'j':
					if (GuiFitPar.hists.size() == 1)
					{
						PrintWarning("Cannot change histogram since only one was added");
					}
					else if (GuiFitPar.current_hist < GuiFitPar.hists.size()-1)
					{
						GuiFitPar.current_hist++;
					}
					else
					{
						GuiFitPar.current_hist = 0;
					}

					GuiFitPar.Draw(true, false);
					
					break;

				case 'k':
					if (GuiFitPar.hists.size() == 1)
					{
						PrintWarning("Cannot change histogram since only one was added");
					}
					else if (GuiFitPar.current_hist > 0)
					{
						GuiFitPar.current_hist--;
					}
					else
					{
						GuiFitPar.current_hist = GuiFitPar.hists.size()-1;
					}

					GuiFitPar.Draw(true, false);
					
					break;
					
				case 'p':
					for (int j = 0; j < GuiFitPar.output_files.size(); j++)
					{
						GuiFitPar.output_files[j].close();
						GuiFitPar.output_files[j].open(GuiFitPar.output_files_name[j]);
						
						for (int i = 0; i < GuiFitPar.hists.size(); i++)
						{
							GuiFitPar.output_files[j] << GuiFitPar.hists_x_val[i] << " ";

							double x_min, x_max;
							GuiFitPar.fits[i][j]->GetRange(x_min, x_max);

							const double x_shift = (x_max - x_min)/
								(static_cast<double>(GuiFitPar.bg_fits[i][j]->GetNpar())-1.);
							
							GuiFitPar.output_files[j] << 
								GuiFitPar.bg_fits[i][j]->GetParameter(0) << " ";
							
							for (int k = 1; k < GuiFitPar.bg_fits[i][j]->GetNpar() - 1; k++)
							{
								const double x_pos = x_min + x_shift*static_cast<double>(k);
								GuiFitPar.output_files[j] << GuiFitPar.bg_fits[i][j]->GetParameter(k) << " ";
							}
							
							GuiFitPar.output_files[j] << 
								GuiFitPar.bg_fits[i][j]->GetParameter(
								GuiFitPar.bg_fits[i][j]->GetNpar() - 1) << std::endl;
						}
					}
					PrintInfo("BG approximation parameters were written");
					break;

				case 'f':
					GuiFitPar.PerformFreeFit();
					break;

				case 'n':
					GuiFitPar.ResetFit();
					break;

				case 'h':
					GuiFitPar.PrintHelp();
					break;
			}
				
			break;
			
			
		case kButton1Down:
			if (!GuiFitPar.IsActivePoint()) 
			{
				for (int i = 0; i < GuiFitPar.bg_points_gr->GetN(); i++)
				{
				const double marker_position_x = gPad->XtoPixel(GuiFitPar.bg_points_gr->GetPointX(i));
				const double marker_position_y = gPad->YtoPixel(GuiFitPar.bg_points_gr->GetPointY(i));

				if (pow(marker_position_x - px, 2) + pow(marker_position_y - py, 2) < 
						pow(4.*GuiFitPar.bg_points_gr->GetMarkerSize(), 2)) 
				{
						GuiFitPar.ActivatePoint(i, 
							GuiFitPar.bg_points_gr->GetPointX(i), 
							GuiFitPar.bg_points_gr->GetPointY(i));
					}
				}
			}
			break;

		case kButton1Motion:
			if (GuiFitPar.IsActivePoint())
			{
				GuiFitPar.current_active_point_gr->SetPointY(0, y);
				GuiFitPar.bg_points_gr->SetPointY(GuiFitPar.current_active_point, y);
				
				GuiFitPar.bg_points_gr->Fit(
					GuiFitPar.bg_fits[GuiFitPar.current_hist][GuiFitPar.current_fit_mode-1], "RQNM");
				
				for (int i = 0; i < GuiFitPar.bg_fits[GuiFitPar.current_hist][GuiFitPar.current_fit_mode-1]->GetNpar(); i++)
				{
					GuiFitPar.fits[GuiFitPar.current_hist][GuiFitPar.current_fit_mode-1]->FixParameter(
						GuiFitPar.bg_fit_par_min[GuiFitPar.current_hist][GuiFitPar.current_fit_mode-1]+i, 
						GuiFitPar.bg_fits[GuiFitPar.current_hist][GuiFitPar.current_fit_mode-1]->GetParameter(i));
				}
				GuiFitPar.fits[GuiFitPar.current_hist][GuiFitPar.current_fit_mode-1]->Update();
				
				GuiFitPar.hists[GuiFitPar.current_hist]->Fit(
					GuiFitPar.fits[GuiFitPar.current_hist][GuiFitPar.current_fit_mode-1], "RQMN");

				GuiFitPar.fits[GuiFitPar.current_hist][GuiFitPar.current_fit_mode-1]->Update();

				std::cout << "Fit parameters: ";

				for (int i = 0; i < GuiFitPar.bg_fit_par_min[GuiFitPar.current_hist][GuiFitPar.current_fit_mode-1]; i++)
				{
					std::cout << 
						GuiFitPar.fits[GuiFitPar.current_hist][GuiFitPar.current_fit_mode-1]->GetParameter(i) 
						<< " ";
				}
				std::cout << std::endl;
			}
			else return;
			break;

		case kButton1Up:
			GuiFitPar.DeactivatePoint();
			break;
	}
	
	GuiFitPar.Draw(false, false);
};
