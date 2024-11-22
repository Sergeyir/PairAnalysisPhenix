#include "IOTools.hpp"
#include "StrTools.hpp"
#include "DeadAreasCuts.hpp"

#define RUN14HEAU200

struct
{
	std::vector<int> key = {1, 2, 3, 4};
	std::vector<std::string> name = {"Rectangle", "LineX", "LineY", "Inverse rectangle"};
	
	std::vector<double> rect_xmin, rect_xmax, rect_ymin, rect_ymax;
	std::vector<double> line_xmin, line_xmax;
	std::vector<double> line_ymin, line_ymax;
	std::vector<double> inv_rect_xmin, inv_rect_xmax, inv_rect_ymin, inv_rect_ymax;
	
	int current_cut_mode = -1;
} CutMode;

struct 
{
	double orig_integral;
	std::array<bool, 4> is_min = {true, true, true, true};
	TFile file = TFile("data/Analysis/Run14HeAu200/sum.root");
	//TFile file = TFile("../data/PHENIX_sim/Run7AuAu200/hotmaps.root");
	TH2F *hist;

	std::string xval_name = "pemcy";
	std::string yval_name = "pemcz";
} Par;

void CutDeadAreas(TH2F *hist, bool (*cut_func)(const double, const double))
{
	for (int i = 1; i <= hist->GetXaxis()->GetNbins(); i++)
	{
		for (int j = 1; j <= hist->GetYaxis()->GetNbins(); j++)
		{
			const double val_x = hist->GetXaxis()->GetBinCenter(i);
			const double val_y = hist->GetYaxis()->GetBinCenter(j);
			
			if (cut_func(val_x, val_y)) 
			{
				hist->SetBinContent(i, j, 0.);
			}
		}
	}
}

void CutDeadAreas(TH2F *hist, bool (*cut_func)(const double, const double, const double), const double aux_val)
{
	for (int i = 1; i <= hist->GetXaxis()->GetNbins(); i++)
	{
		for (int j = 1; j <= hist->GetYaxis()->GetNbins(); j++)
		{
			const double val_x = hist->GetXaxis()->GetBinCenter(i);
			const double val_y = hist->GetYaxis()->GetBinCenter(j);
			
			if (cut_func(aux_val, val_x, val_y)) 
			{
				hist->SetBinContent(i, j, 0.);
			}
		}
	}
}

void CutDeadAreas(TH2F *hist, bool (*cut_func)(const double, const double, const double, const double), const double aux_val1, const double aux_val2)
{
	for (int i = 1; i <= hist->GetXaxis()->GetNbins(); i++)
	{
		for (int j = 1; j <= hist->GetYaxis()->GetNbins(); j++)
		{
			const double val_x = hist->GetXaxis()->GetBinCenter(i);
			const double val_y = hist->GetYaxis()->GetBinCenter(j);
			
			if (cut_func(aux_val1, aux_val2, val_x, val_y)) 
			{
				hist->SetBinContent(i, j, 0.);
			}
		}
	}
}

void CutDeadAreas(TH2F *hist, bool (*cut_func)(const double, const double, const int, const double, const double), const double aux_val1, const double aux_val2, const int aux_val3)
{
	for (int i = 1; i <= hist->GetXaxis()->GetNbins(); i++)
	{
		for (int j = 1; j <= hist->GetYaxis()->GetNbins(); j++)
		{
			const double val_x = hist->GetXaxis()->GetBinCenter(i);
			const double val_y = hist->GetYaxis()->GetBinCenter(j);
			
			if (cut_func(aux_val1, aux_val2, aux_val3, val_x, val_y)) 
			{
				hist->SetBinContent(i, j, 0.);
			}
		}
	}
}

void drawdm(bool is_fixed_range = false)
{
	TH2F *hist_dm = (TH2F *) Par.hist->Clone();
	
	for (int i = 0; i < hist_dm->GetXaxis()->GetNbins(); i++)
	{
		for (int j = 0; j < hist_dm->GetYaxis()->GetNbins(); j++)
		{
			double xval = hist_dm->GetXaxis()->GetBinCenter(i);
			double yval = hist_dm->GetYaxis()->GetBinCenter(j);
			
			for (int k = 0; k < CutMode.rect_xmax.size(); k++)
			{
				if (xval > CutMode.rect_xmin[k] && xval < CutMode.rect_xmax[k] && 
					yval > CutMode.rect_ymin[k] && yval < CutMode.rect_ymax[k]) 
					hist_dm->SetBinContent(hist_dm->GetBin(i, j), 0.);
			}
			for (int k = 0; k < CutMode.line_xmax.size(); k++)
			{
				if (xval > CutMode.line_xmin[k] && xval < CutMode.line_xmax[k]) 
					hist_dm->SetBinContent(hist_dm->GetBin(i, j), 0.);
			}
			for (int k = 0; k < CutMode.line_ymax.size(); k++)
			{
				if (yval > CutMode.line_ymin[k] && yval < CutMode.line_ymax[k]) 
					hist_dm->SetBinContent(hist_dm->GetBin(i, j), 0.);
			}
			for (int k = 0; k < CutMode.inv_rect_xmax.size(); k++)
			{
				if (xval < CutMode.inv_rect_xmin[k] || xval > CutMode.inv_rect_xmax[k] ||
					yval < CutMode.inv_rect_ymin[k] || yval > CutMode.inv_rect_ymax[k]) 
					hist_dm->SetBinContent(hist_dm->GetBin(i, j), 0.);
			}
		}
	}

	if (is_fixed_range)
	{
		const int xmin = hist_dm->GetXaxis()->FindBin(gPad->GetUxmin());
		const int xmax = hist_dm->GetXaxis()->FindBin(gPad->GetUxmax())-1;
		const int ymin = hist_dm->GetYaxis()->FindBin(gPad->GetUymin())+1;
		const int ymax = hist_dm->GetYaxis()->FindBin(gPad->GetUymax())-1;
		
		hist_dm->GetXaxis()->SetRange(xmin, xmax);
		hist_dm->GetYaxis()->SetRange(ymin, ymax);
	}
	
	hist_dm->Draw("COLZ");
	gPad->Modified();
	gPad->Update();

	PrintInfo("Data lost " + to_string((1.-hist_dm->Integral()/Par.orig_integral)*100.) + "\%");
}

void SetTLineStyle(TLine *line)
{
	line->SetLineColor(kRed);
	line->SetLineStyle(2);
	line->SetLineWidth(2);
}

void DrawCutLines(double xmax, double ymax, const int current_cut_mode)
{
	switch (CutMode.current_cut_mode)
	{
		case 1:
			{
				TLine line1 = TLine(CutMode.rect_xmin.back(), CutMode.rect_ymin.back(), 
					xmax, CutMode.rect_ymin.back());
				TLine line2 = TLine(CutMode.rect_xmin.back(), CutMode.rect_ymin.back(), 
					CutMode.rect_xmin.back(), ymax);
				TLine line3 = TLine(xmax, CutMode.rect_ymin.back(), xmax, ymax);
				TLine line4 = TLine(CutMode.rect_xmin.back(), ymax, xmax, ymax);

				SetTLineStyle(&line1);
				SetTLineStyle(&line2);
				SetTLineStyle(&line3);
				SetTLineStyle(&line4);
				
				line1.Draw();
				line2.Draw();
				line3.Draw();
				line4.Draw();

				gPad->Modified();
				gPad->Update();
			}
			break;
		case 2:
			{
				const double ymin = Par.hist->GetYaxis()->GetBinLowEdge(1);
				ymax = Par.hist->GetYaxis()->GetBinUpEdge(Par.hist->GetYaxis()->GetNbins());
				
				TLine line1 = TLine(CutMode.line_xmin.back(), ymin, CutMode.line_xmin.back(), ymax);
				TLine line2 = TLine(xmax, ymin, xmax, ymax);

				SetTLineStyle(&line1);
				SetTLineStyle(&line2);

				line1.Draw();
				line2.Draw();

				gPad->Modified();
				gPad->Update();
			}
			break;
		case 3:
			{
				const double xmin = Par.hist->GetXaxis()->GetBinLowEdge(1);
				xmax = Par.hist->GetXaxis()->GetBinUpEdge(Par.hist->GetXaxis()->GetNbins());
				
				TLine line1 = TLine(xmin, CutMode.line_ymin.back(), xmax, CutMode.line_ymin.back());
				TLine line2 = TLine(xmin, ymax, xmax, ymax);

				SetTLineStyle(&line1);
				SetTLineStyle(&line2);

				line1.Draw();
				line2.Draw();

				gPad->Modified();
				gPad->Update();
			}
			break;
		case 4:
			{
				TLine line1 = TLine(CutMode.inv_rect_xmin.back(), CutMode.inv_rect_ymin.back(), 
					xmax, CutMode.inv_rect_ymin.back());
				TLine line2 = TLine(CutMode.inv_rect_xmin.back(), CutMode.inv_rect_ymin.back(), 
					CutMode.inv_rect_xmin.back(), ymax);
				TLine line3 = TLine(xmax, CutMode.inv_rect_ymin.back(), xmax, ymax);
				TLine line4 = TLine(CutMode.inv_rect_xmin.back(), ymax, xmax, ymax);

				SetTLineStyle(&line1);
				SetTLineStyle(&line2);
				SetTLineStyle(&line3);
				SetTLineStyle(&line4);
				
				line1.Draw();
				line2.Draw();
				line3.Draw();
				line4.Draw();

				gPad->Modified();
				gPad->Update();
			}
			break;
	}
}

void exec()
{
	int event = gPad->GetEvent();
	int px = gPad->GetEventX();
	int py = gPad->GetEventY();

	const double x = gPad->PadtoX(gPad->AbsPixeltoX(px));
	const double y = gPad->PadtoY(gPad->AbsPixeltoY(py));

	if (event != kKeyPress && (!Par.is_min[0] && CutMode.current_cut_mode == CutMode.key[0]))
		DrawCutLines(x, y, CutMode.current_cut_mode);
	else if (event != kKeyPress && (!Par.is_min[1] && CutMode.current_cut_mode == CutMode.key[1]))
		DrawCutLines(x, 0., CutMode.current_cut_mode);
	else if (event != kKeyPress && (!Par.is_min[2] && CutMode.current_cut_mode == CutMode.key[2]))
		DrawCutLines(0., y, CutMode.current_cut_mode);
	else if (event != kKeyPress && (!Par.is_min[3] && CutMode.current_cut_mode == CutMode.key[3]))
		DrawCutLines(x, y, CutMode.current_cut_mode);
	//Print(event);

	if (event == kButton1Down)
	{
		switch (CutMode.current_cut_mode)
		{
			case 1:
				if (Par.is_min[0])
				{
					CutMode.rect_xmin.push_back(x);
					CutMode.rect_ymin.push_back(y);
					
					Par.is_min[0] = false;

					PrintInfo("Setting the first point");	
				}
				else
				{
					if (x >= CutMode.rect_xmin.back()) 
					{
						CutMode.rect_xmax.push_back(Par.hist->GetXaxis()->GetBinUpEdge(
							Par.hist->GetXaxis()->FindBin(x)));
						CutMode.rect_xmin.back() = Par.hist->GetXaxis()->GetBinLowEdge(
							Par.hist->GetXaxis()->FindBin(CutMode.rect_xmin.back()));
					}
					else 
					{
						CutMode.rect_xmax.push_back(Par.hist->GetXaxis()->GetBinUpEdge(
							Par.hist->GetXaxis()->FindBin(CutMode.rect_xmin.back())));
						CutMode.rect_xmin.back() = Par.hist->GetXaxis()->GetBinLowEdge(
							Par.hist->GetXaxis()->FindBin(x));
					}
					
					if (y >= CutMode.rect_ymin.back()) 
					{
						CutMode.rect_ymax.push_back(Par.hist->GetYaxis()->GetBinUpEdge(
							Par.hist->GetYaxis()->FindBin(y)));
						CutMode.rect_ymin.back() = Par.hist->GetYaxis()->GetBinLowEdge(
							Par.hist->GetYaxis()->FindBin(CutMode.rect_ymin.back()));
					}
					else 
					{
						CutMode.rect_ymax.push_back(Par.hist->GetYaxis()->GetBinUpEdge(
							Par.hist->GetYaxis()->FindBin(CutMode.rect_ymin.back())));
						CutMode.rect_ymin.back() = Par.hist->GetYaxis()->GetBinLowEdge(
							Par.hist->GetYaxis()->FindBin(y));
					}

					Par.is_min[0] = true;
					PrintInfo("Setting the second point");
					if (Par.is_min[0]) drawdm(true);
				}
				break;

			case 2:
				if (Par.is_min[1])
				{
					CutMode.line_xmin.push_back(x);
					
					Par.is_min[1] = false;

					PrintInfo("Setting the first line");	
				}
				else
				{
					if (x >= CutMode.line_xmin.back()) 
					{
						CutMode.line_xmax.push_back(Par.hist->GetXaxis()->GetBinUpEdge(
							Par.hist->GetXaxis()->FindBin(x)));
						CutMode.line_xmin.back() = Par.hist->GetXaxis()->GetBinLowEdge(
							Par.hist->GetXaxis()->FindBin(CutMode.line_xmin.back()));
					}
					else 
					{
						CutMode.line_xmax.push_back(Par.hist->GetXaxis()->GetBinUpEdge(
							Par.hist->GetXaxis()->FindBin(CutMode.line_xmin.back())));
						CutMode.line_xmin.back() = Par.hist->GetXaxis()->GetBinLowEdge(
							Par.hist->GetXaxis()->FindBin(x));
					}

					Par.is_min[1] = true;
					PrintInfo("Setting the second point");
					if (Par.is_min[1]) drawdm(true);
				}
				break;

			case 3:
				if (Par.is_min[2])
				{
					CutMode.line_ymin.push_back(y);
					
					Par.is_min[2] = false;

					PrintInfo("Setting the first line");	
				}
				else
				{
					if (y >= CutMode.line_ymin.back()) 
					{
						CutMode.line_ymax.push_back(Par.hist->GetYaxis()->GetBinUpEdge(
							Par.hist->GetYaxis()->FindBin(y)));
						CutMode.line_ymin.back() = Par.hist->GetYaxis()->GetBinLowEdge(
							Par.hist->GetYaxis()->FindBin(CutMode.line_ymin.back()));
					}
					else 
					{
						CutMode.line_ymax.push_back(Par.hist->GetYaxis()->GetBinUpEdge(
							Par.hist->GetYaxis()->FindBin(CutMode.line_ymin.back())));
						CutMode.line_ymin.back() = Par.hist->GetYaxis()->GetBinLowEdge(
							Par.hist->GetYaxis()->FindBin(y));
					}

					Par.is_min[2] = true;
					PrintInfo("Setting the second point");
					if (Par.is_min[2]) drawdm(true);
				}
				break;

			case 4:
				if (Par.is_min[3])
				{
					CutMode.inv_rect_xmin.push_back(x);
					CutMode.inv_rect_ymin.push_back(y);
					
					Par.is_min[3] = false;

					PrintInfo("Setting the first point");	
				}
				else
				{
					if (x >= CutMode.inv_rect_xmin.back()) 
					{
						CutMode.inv_rect_xmax.push_back(Par.hist->GetXaxis()->GetBinLowEdge(
							Par.hist->GetXaxis()->FindBin(x)));
						CutMode.inv_rect_xmin.back() = Par.hist->GetXaxis()->GetBinUpEdge(
							Par.hist->GetXaxis()->FindBin(CutMode.inv_rect_xmin.back()));
					}
					else 
					{
						CutMode.inv_rect_xmax.push_back(Par.hist->GetXaxis()->GetBinLowEdge(
							Par.hist->GetXaxis()->FindBin(CutMode.inv_rect_xmin.back())));
						CutMode.inv_rect_xmin.back() = Par.hist->GetXaxis()->GetBinUpEdge(
							Par.hist->GetXaxis()->FindBin(x));
					}
					
					if (y >= CutMode.inv_rect_ymin.back()) 
					{
						CutMode.inv_rect_ymax.push_back(Par.hist->GetYaxis()->GetBinLowEdge(
							Par.hist->GetYaxis()->FindBin(y)));
						CutMode.inv_rect_ymin.back() = Par.hist->GetYaxis()->GetBinUpEdge(
							Par.hist->GetYaxis()->FindBin(CutMode.inv_rect_ymin.back()));
					}
					else 
					{
						CutMode.inv_rect_ymax.push_back(Par.hist->GetYaxis()->GetBinLowEdge(
							Par.hist->GetYaxis()->FindBin(CutMode.inv_rect_ymin.back())));
						CutMode.inv_rect_ymin.back() = Par.hist->GetYaxis()->GetBinUpEdge(
							Par.hist->GetYaxis()->FindBin(y));
					}

					Par.is_min[3] = true;
					PrintInfo("Setting the second point");
					if (Par.is_min[3]) drawdm(true);
				}
				break;
		}
	}

	switch(px)
	{
		case 'u':
			if (CutMode.current_cut_mode == 1)
			{
				if (CutMode.rect_xmin.size() != 0)
				{
					if (Par.is_min[0])
					{
						CutMode.rect_xmin.pop_back();
						CutMode.rect_ymin.pop_back();
						CutMode.rect_xmax.pop_back();
						CutMode.rect_ymax.pop_back();

						PrintInfo("Deleting last minimum and maximum points");
					}
					else
					{
						CutMode.rect_xmin.pop_back();
						CutMode.rect_ymin.pop_back();
							
						PrintInfo("Deleting last minimum point");
					}
					if (Par.is_min[0]) drawdm(true);
					Par.is_min[0] = true;
					gPad->Update();
				}
				else PrintInfo("Cannot delete last point since the current number of points is 0");
			}
			else if (CutMode.current_cut_mode == 2)
			{
				if (CutMode.line_xmin.size() != 0)
				{
					if (Par.is_min[1])
					{
						CutMode.line_xmin.pop_back();
						CutMode.line_xmax.pop_back();

						PrintInfo("Deleting last pair of lines");
					}
					else
					{
						CutMode.line_xmin.pop_back();
						PrintInfo("Deleting last line");
					}
					if (Par.is_min[1]) drawdm(true);
					Par.is_min[1] = true;
					gPad->Update();
				}
				else PrintInfo("Cannot delete last line/lines since the current number of lines is 0");
			}
			else if (CutMode.current_cut_mode == 3)
			{
				if (CutMode.line_ymin.size() != 0)
				{
					if (Par.is_min[2])
					{
						CutMode.line_ymin.pop_back();
						CutMode.line_ymax.pop_back();

						PrintInfo("Deleting last pair of lines");
					}
					else
					{
						CutMode.line_ymin.pop_back();
						PrintInfo("Deleting last line");
					}
					if (Par.is_min[2]) drawdm(true);
					Par.is_min[2] = true;
					gPad->Update();
				}
				else PrintInfo("Cannot delete last line/lines since the current number of lines is 0");
			}
			if (CutMode.current_cut_mode == 4)
			{
				if (CutMode.inv_rect_xmin.size() != 0)
				{
					if (Par.is_min[3])
					{
						CutMode.inv_rect_xmin.pop_back();
						CutMode.inv_rect_ymin.pop_back();
						CutMode.inv_rect_xmax.pop_back();
						CutMode.inv_rect_ymax.pop_back();

						PrintInfo("Deleting last minimum and maximum points");
					}
					else
					{
						CutMode.inv_rect_xmin.pop_back();
						CutMode.inv_rect_ymin.pop_back();
							
						PrintInfo("Deleting last minimum point");
					}
					if (Par.is_min[3]) drawdm(true);
					Par.is_min[3] = true;
					gPad->Update();
				}
				else PrintInfo("Cannot delete last point since the current number of points is 0");
			}
			break;
		case 'r':
			drawdm();
			PrintInfo("Resetting the range of the selected area");
			break;
		case 'p':
			for (int i = 0; i < CutMode.inv_rect_xmax.size(); i++)
			{
				std::cout <<
					"if (" + Par.xval_name + " < " << CutMode.inv_rect_xmin[i] <<
					" || " + Par.yval_name + " < " << CutMode.inv_rect_ymin[i] <<
					std::endl << "	|| " + Par.xval_name + " > " << CutMode.inv_rect_xmax[i] <<
					" || " + Par.yval_name + " > " << CutMode.inv_rect_ymax[i] <<
					") return true;" << std::endl;
			}
			for (int i = 0; i < CutMode.line_xmax.size(); i++)
			{
				std::cout <<
					"if (" + Par.xval_name + " > " << CutMode.line_xmin[i] <<
					" && " + Par.xval_name + " < " << CutMode.line_xmax[i] <<
					") return true;" << std::endl;
			}
			for (int i = 0; i < CutMode.line_ymax.size(); i++)
			{
				std::cout <<
					"if (" + Par.yval_name + " > " << CutMode.line_ymin[i] <<
					" && " + Par.yval_name + " < " << CutMode.line_ymax[i] <<
					") return true;" << std::endl;
			}
			for (int i = 0; i < CutMode.rect_xmax.size(); i++)
			{
				std::cout << "if (" + Par.xval_name + " > " << CutMode.rect_xmin[i] <<
					" && " + Par.yval_name + " > " << CutMode.rect_ymin[i] <<
					std::endl << "	&& " + Par.xval_name + " < " << CutMode.rect_xmax[i] <<
					" && " + Par.yval_name + " < " << CutMode.rect_ymax[i] <<
					") return true;" << std::endl;
			}
			break;
		case '0':
			Print("Deactivating cutting mode");
			CutMode.current_cut_mode = -1;
			break;
		case '1':
			Print("Activating rectangular cutting mode");
			CutMode.current_cut_mode = 1;
			break;
		case '2':
			Print("Activating linear cutting mode along x axis");
			CutMode.current_cut_mode = 2;
			break;
		case '3':
			Print("Activating linear cutting mode along y axis");
			CutMode.current_cut_mode = 3;
			break;
		case '4':
			Print("Activating inverse rectangular cutting mode");
			CutMode.current_cut_mode = 4;
			break;
	};
}

void GUIDM()
{
	gStyle->SetOptStat(0);
	Par.hist = (TH2F *) Par.file.Get("Heatmap: DCe0");
	Par.hist->SetTitle(Par.hist->GetName());

	Par.orig_integral = Par.hist->Integral();
	
	CutDeadAreas(Par.hist, &IsDeadDC, 1., 1.);
	//CutDeadAreas(Par.hist, &IsDeadPC1, 1.);
	//CutDeadAreas(Par.hist, &IsDeadPC2);
	//CutDeadAreas(Par.hist, &IsDeadPC3, 1.);
	//CutDeadAreas(Par.hist, &IsDeadEMCal, 2., 1., 3);
	//CutDeadAreas(Par.hist, &IsDeadTOFe, -1.);
	//CutDeadAreas(Par.hist, &IsDeadTOFw, -1.);

	TCanvas *canv = new TCanvas("canv");
   
	drawdm();
	
	gPad->AddExec("exec", "exec()");
}
