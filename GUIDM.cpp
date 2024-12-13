#include "IOTools.hpp"
#include "StrTools.hpp"

//#include "DataCutsSelector.hpp"
#include "DeadAreasCuts.h"

#include "MathTools.hpp"
#include "IOTools.hpp"

struct
{
	std::vector<std::string> name = {"Rectangle", "LineX", "LineY", "Inverse rectangle", "Angled line"};
	
	std::vector<double> rectXMin, rectXMax, rectYMin, rectYMax;
	std::vector<double> lineXMin, lineXMax;
	std::vector<double> lineYMin, lineYMax;
	std::vector<double> invRectXMin, invRectXMax, invRectYMin, invRectYMax;
	std::vector<double> angledLine1X1, angledLine1X2, angledLine1Y1, angledLine1Y2;
	std::vector<double> angledLine2X1, angledLine2X2, angledLine2Y1, angledLine2Y2;
	std::vector<double> tanAlpha1, tanAlpha2, shiftY1, shiftY2;
	
	int currentCutMode = -1;
} CutMode;

struct 
{
	double orig_integral;
	std::array<bool, 5> isMin = {true, true, true, true, true};
	TFile file = TFile("data/Real/Run14HeAu200/sum.root");
	//TFile file = TFile("../data/PHENIX_sim/Run7AuAu200/hotmaps.root");
	TH2F *hist;

	std::string xValName = "board";
	std::string yValName = "alpha";
} Par;

double Pol1(const double par0, const double par1, const double x)
{
   return par0 + x*par1;
}

void CutDeadAreas(TH2F *hist, bool (*cut_func)(const double, const double))
{
	for (int i = 1; i <= hist->GetXaxis()->GetNbins(); i++)
	{
		for (int j = 1; j <= hist->GetYaxis()->GetNbins(); j++)
		{
			const double valX = hist->GetXaxis()->GetBinCenter(i);
			const double valY = hist->GetYaxis()->GetBinCenter(j);
			
			if (cut_func(valX, valY)) 
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
			const double valX = hist->GetXaxis()->GetBinCenter(i);
			const double valY = hist->GetYaxis()->GetBinCenter(j);
			
			if (cut_func(aux_val, valX, valY)) 
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
			const double valX = hist->GetXaxis()->GetBinCenter(i);
			const double valY = hist->GetYaxis()->GetBinCenter(j);
			
			if (cut_func(aux_val1, aux_val2, valX, valY)) 
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
			const double valX = hist->GetXaxis()->GetBinCenter(i);
			const double valY = hist->GetYaxis()->GetBinCenter(j);
			
			if (cut_func(aux_val1, aux_val2, aux_val3, valX, valY)) 
			{
				hist->SetBinContent(i, j, 0.);
			}
		}
	}
}

void DrawDM(bool is_fixed_range = false)
{
	TH2F *histDM = (TH2F *) Par.hist->Clone();
	
	for (int i = 0; i < histDM->GetXaxis()->GetNbins(); i++)
	{
		for (int j = 0; j < histDM->GetYaxis()->GetNbins(); j++)
		{
			double xval = histDM->GetXaxis()->GetBinCenter(i);
			double yval = histDM->GetYaxis()->GetBinCenter(j);
			
			for (int k = 0; k < CutMode.rectXMax.size(); k++)
			{
				if (xval > CutMode.rectXMin[k] && xval < CutMode.rectXMax[k] && 
					yval > CutMode.rectYMin[k] && yval < CutMode.rectYMax[k]) 
					histDM->SetBinContent(histDM->GetBin(i, j), 0.);
			}
			for (int k = 0; k < CutMode.lineXMax.size(); k++)
			{
				if (xval > CutMode.lineXMin[k] && xval < CutMode.lineXMax[k]) 
					histDM->SetBinContent(histDM->GetBin(i, j), 0.);
			}
			for (int k = 0; k < CutMode.lineYMin.size(); k++)
			{
				if (yval > CutMode.lineYMin[k] && yval < CutMode.lineYMax[k]) 
					histDM->SetBinContent(histDM->GetBin(i, j), 0.);
			}
			for (int k = 0; k < CutMode.invRectXMax.size(); k++)
			{
				if (xval < CutMode.invRectXMin[k] || xval > CutMode.invRectXMax[k] ||
					yval < CutMode.invRectYMin[k] || yval > CutMode.invRectYMax[k]) 
					histDM->SetBinContent(histDM->GetBin(i, j), 0.);
			}
			for (int k = 0; k < CutMode.shiftY2.size(); k++)
			{ 
            const double y1Cut = Pol1(CutMode.shiftY1[k], CutMode.tanAlpha1[k], xval);
            const double y2Cut = Pol1(CutMode.shiftY2[k], CutMode.tanAlpha2[k], xval);
            
				if (yval > Minimum(y1Cut, y2Cut) && yval < Maximum(y1Cut, y2Cut)) 
            {
					histDM->SetBinContent(histDM->GetBin(i, j), 0.);
            }
			}
		}
	}

	if (is_fixed_range)
	{
		const int xMin = histDM->GetXaxis()->FindBin(gPad->GetUxmin());
		const int xMax = histDM->GetXaxis()->FindBin(gPad->GetUxmax())-1;
		const int yMin = histDM->GetYaxis()->FindBin(gPad->GetUymin())+1;
		const int yMax = histDM->GetYaxis()->FindBin(gPad->GetUymax())-1;
		
		histDM->GetXaxis()->SetRange(xMin, xMax);
		histDM->GetYaxis()->SetRange(yMin, yMax);
	}
	
	histDM->Draw("COLZ");
	gPad->Modified();
	gPad->Update();
   gPad->GetFrame()->SetBit(TBox::kCannotMove);

	PrintInfo("Data lost " + to_string((1.-histDM->Integral()/Par.orig_integral)*100.) + "\%");
}

void SetTLineStyle(TLine *line, const EColor color=kRed)
{
	line->SetLineColor(color);
	line->SetLineStyle(2);
	line->SetLineWidth(2);
}

void DrawCutLines(double xMax, double yMax, const int currentCutMode)
{
   if (currentCutMode < 0) return;
   
	switch (CutMode.currentCutMode)
	{
		case 0:
      {
         if (Par.isMin[currentCutMode]) return;
         
         TLine line1 = TLine(CutMode.rectXMin.back(), CutMode.rectYMin.back(), 
            xMax, CutMode.rectYMin.back());
         TLine line2 = TLine(CutMode.rectXMin.back(), CutMode.rectYMin.back(), 
            CutMode.rectXMin.back(), yMax);
         TLine line3 = TLine(xMax, CutMode.rectYMin.back(), xMax, yMax);
         TLine line4 = TLine(CutMode.rectXMin.back(), yMax, xMax, yMax);

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
         gPad->GetFrame()->SetBit(TBox::kCannotMove);
         
			break;
      }
		case 1:
      {
         if (Par.isMin[currentCutMode]) return;
         
         const double yMin = Par.hist->GetYaxis()->GetBinLowEdge(1);
         yMax = Par.hist->GetYaxis()->GetBinUpEdge(Par.hist->GetYaxis()->GetNbins());
         
         TLine line1 = TLine(CutMode.lineXMin.back(), yMin, CutMode.lineXMin.back(), yMax);
         TLine line2 = TLine(xMax, yMin, xMax, yMax);

         SetTLineStyle(&line1);
         SetTLineStyle(&line2);

         line1.Draw();
         line2.Draw();

         gPad->Modified();
         gPad->Update();
         gPad->GetFrame()->SetBit(TBox::kCannotMove);
         
			break;
      }
		case 2:
      {
         if (Par.isMin[currentCutMode]) return;
         
         const double xMin = Par.hist->GetXaxis()->GetBinLowEdge(1);
         xMax = Par.hist->GetXaxis()->GetBinUpEdge(Par.hist->GetXaxis()->GetNbins());
         
         TLine line1 = TLine(xMin, CutMode.lineYMin.back(), xMax, CutMode.lineYMin.back());
         TLine line2 = TLine(xMin, yMax, xMax, yMax);
         
         SetTLineStyle(&line1);
         SetTLineStyle(&line2);
         
         line1.Draw();
         line2.Draw();
         
         gPad->Modified();
         gPad->Update();
         gPad->GetFrame()->SetBit(TBox::kCannotMove);
         
			break;
      }
		case 3:
      {
         if (Par.isMin[currentCutMode]) return;
         
         TLine line1 = TLine(CutMode.invRectXMin.back(), CutMode.invRectYMin.back(), 
            xMax, CutMode.invRectYMin.back());
         TLine line2 = TLine(CutMode.invRectXMin.back(), CutMode.invRectYMin.back(), 
            CutMode.invRectXMin.back(), yMax);
         TLine line3 = TLine(xMax, CutMode.invRectYMin.back(), xMax, yMax);
         TLine line4 = TLine(CutMode.invRectXMin.back(), yMax, xMax, yMax);

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
         gPad->GetFrame()->SetBit(TBox::kCannotMove);
         
			break;
      }
      case 4:
         const double lineXMin = Par.hist->GetXaxis()->GetBinLowEdge(1);
         const double lineXMax = 
            Par.hist->GetXaxis()->GetBinUpEdge(Par.hist->GetXaxis()->GetNbins());
         
         if (CutMode.angledLine1X1.size() > CutMode.angledLine2X1.size())
         {
            if (CutMode.angledLine1X1.size() > CutMode.angledLine1X2.size())
            {
               if (Par.isMin[currentCutMode]) return;
               
               const double tanAlpha = (CutMode.angledLine1Y1.back() - yMax)/
                                 (CutMode.angledLine1X1.back() - xMax);
               
               TLine line = TLine(lineXMin, Pol1(yMax - xMax*tanAlpha, tanAlpha, lineXMin), 
                                   lineXMax, Pol1(yMax - xMax*tanAlpha, tanAlpha, lineXMax));

               SetTLineStyle(&line);
               
               line.Draw();
               
               gPad->Modified();
               gPad->Update();
               gPad->GetFrame()->SetBit(TBox::kCannotMove);
            }
            else
            {
               TLine line = TLine(lineXMin, 
                                  Pol1(CutMode.shiftY1.back(), CutMode.tanAlpha1.back(), lineXMin), 
                                  lineXMax, 
                                  Pol1(CutMode.shiftY1.back(), CutMode.tanAlpha1.back(), lineXMax));
               
               SetTLineStyle(&line, kGray);
               
               line.Draw();
               
               gPad->Modified();
               gPad->Update();
               gPad->GetFrame()->SetBit(TBox::kCannotMove);
            }
         }
         else if (CutMode.angledLine2X1.size() > 0 && !Par.isMin[currentCutMode])
         {
            const double tanAlpha = (CutMode.angledLine2Y1.back() - yMax)/
                              (CutMode.angledLine2X1.back() - xMax);
            
            TLine line1 = TLine(lineXMin, Pol1(yMax - xMax*tanAlpha, tanAlpha, lineXMin), 
                                lineXMax, Pol1(yMax - xMax*tanAlpha, tanAlpha, lineXMax));
            
            TLine line2 = TLine(lineXMin, 
                                Pol1(CutMode.shiftY1.back(), CutMode.tanAlpha1.back(), lineXMin), 
                                lineXMax, 
                                Pol1(CutMode.shiftY1.back(), CutMode.tanAlpha1.back(), lineXMax));
            
            SetTLineStyle(&line1);
            SetTLineStyle(&line2, kGray);
            
            line1.Draw();
            line2.Draw();
            
            gPad->Modified();
            gPad->Update();
            gPad->GetFrame()->SetBit(TBox::kCannotMove);
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
   
   if (event == kMouseMotion) DrawCutLines(x, y, CutMode.currentCutMode);
   
	if (event == kButton1Down)
	{
		switch (CutMode.currentCutMode)
		{
			case 0:
				if (Par.isMin[0])
				{
					CutMode.rectXMin.push_back(x);
					CutMode.rectYMin.push_back(y);
					
					Par.isMin[0] = false;

					PrintInfo("Setting the first point");	
				}
				else
				{
					if (x >= CutMode.rectXMin.back()) 
					{
						CutMode.rectXMax.push_back(Par.hist->GetXaxis()->GetBinUpEdge(
							Par.hist->GetXaxis()->FindBin(x)));
						CutMode.rectXMin.back() = Par.hist->GetXaxis()->GetBinLowEdge(
							Par.hist->GetXaxis()->FindBin(CutMode.rectXMin.back()));
					}
					else 
					{
						CutMode.rectXMax.push_back(Par.hist->GetXaxis()->GetBinUpEdge(
							Par.hist->GetXaxis()->FindBin(CutMode.rectXMin.back())));
						CutMode.rectXMin.back() = Par.hist->GetXaxis()->GetBinLowEdge(
							Par.hist->GetXaxis()->FindBin(x));
					}
					
					if (y >= CutMode.rectYMin.back()) 
					{
						CutMode.rectYMax.push_back(Par.hist->GetYaxis()->GetBinUpEdge(
							Par.hist->GetYaxis()->FindBin(y)));
						CutMode.rectYMin.back() = Par.hist->GetYaxis()->GetBinLowEdge(
							Par.hist->GetYaxis()->FindBin(CutMode.rectYMin.back()));
					}
					else 
					{
						CutMode.rectYMax.push_back(Par.hist->GetYaxis()->GetBinUpEdge(
							Par.hist->GetYaxis()->FindBin(CutMode.rectYMin.back())));
						CutMode.rectYMin.back() = Par.hist->GetYaxis()->GetBinLowEdge(
							Par.hist->GetYaxis()->FindBin(y));
					}

					Par.isMin[0] = true;
					PrintInfo("Setting the second point");
					if (Par.isMin[0]) DrawDM(true);
				}
				break;

			case 1:
				if (Par.isMin[1])
				{
					CutMode.lineXMin.push_back(x);
					
					Par.isMin[1] = false;

					PrintInfo("Setting the first line");	
				}
				else
				{
					if (x >= CutMode.lineXMin.back()) 
					{
						CutMode.lineXMax.push_back(Par.hist->GetXaxis()->GetBinUpEdge(
							Par.hist->GetXaxis()->FindBin(x)));
						CutMode.lineXMin.back() = Par.hist->GetXaxis()->GetBinLowEdge(
							Par.hist->GetXaxis()->FindBin(CutMode.lineXMin.back()));
					}
					else 
					{
						CutMode.lineXMax.push_back(Par.hist->GetXaxis()->GetBinUpEdge(
							Par.hist->GetXaxis()->FindBin(CutMode.lineXMin.back())));
						CutMode.lineXMin.back() = Par.hist->GetXaxis()->GetBinLowEdge(
							Par.hist->GetXaxis()->FindBin(x));
					}

					Par.isMin[1] = true;
					PrintInfo("Setting the second point");
					if (Par.isMin[1]) DrawDM(true);
				}
				break;

			case 2:
				if (Par.isMin[2])
				{
					CutMode.lineYMin.push_back(y);
					
					Par.isMin[2] = false;

					PrintInfo("Setting the first line");	
				}
				else
				{
					if (y >= CutMode.lineYMin.back()) 
					{
						CutMode.lineYMax.push_back(Par.hist->GetYaxis()->GetBinUpEdge(
							Par.hist->GetYaxis()->FindBin(y)));
						CutMode.lineYMin.back() = Par.hist->GetYaxis()->GetBinLowEdge(
							Par.hist->GetYaxis()->FindBin(CutMode.lineYMin.back()));
					}
					else 
					{
						CutMode.lineYMax.push_back(Par.hist->GetYaxis()->GetBinUpEdge(
							Par.hist->GetYaxis()->FindBin(CutMode.lineYMin.back())));
						CutMode.lineYMin.back() = Par.hist->GetYaxis()->GetBinLowEdge(
							Par.hist->GetYaxis()->FindBin(y));
					}

					Par.isMin[2] = true;
					PrintInfo("Setting the second point");
					if (Par.isMin[2]) DrawDM(true);
				}
				break;

			case 3:
				if (Par.isMin[3])
				{
					CutMode.invRectXMin.push_back(x);
					CutMode.invRectYMin.push_back(y);
					
					Par.isMin[3] = false;

					PrintInfo("Setting the first point");	
				}
				else
				{
					if (x >= CutMode.invRectXMin.back()) 
					{
						CutMode.invRectXMax.push_back(Par.hist->GetXaxis()->GetBinLowEdge(
							Par.hist->GetXaxis()->FindBin(x)));
						CutMode.invRectXMin.back() = Par.hist->GetXaxis()->GetBinUpEdge(
							Par.hist->GetXaxis()->FindBin(CutMode.invRectXMin.back()));
					}
					else 
					{
						CutMode.invRectXMax.push_back(Par.hist->GetXaxis()->GetBinLowEdge(
							Par.hist->GetXaxis()->FindBin(CutMode.invRectXMin.back())));
						CutMode.invRectXMin.back() = Par.hist->GetXaxis()->GetBinUpEdge(
							Par.hist->GetXaxis()->FindBin(x));
					}
					
					if (y >= CutMode.invRectYMin.back()) 
					{
						CutMode.invRectYMax.push_back(Par.hist->GetYaxis()->GetBinLowEdge(
							Par.hist->GetYaxis()->FindBin(y)));
						CutMode.invRectYMin.back() = Par.hist->GetYaxis()->GetBinUpEdge(
							Par.hist->GetYaxis()->FindBin(CutMode.invRectYMin.back()));
					}
					else 
					{
						CutMode.invRectYMax.push_back(Par.hist->GetYaxis()->GetBinLowEdge(
							Par.hist->GetYaxis()->FindBin(CutMode.invRectYMin.back())));
						CutMode.invRectYMin.back() = Par.hist->GetYaxis()->GetBinUpEdge(
							Par.hist->GetYaxis()->FindBin(y));
					}

					Par.isMin[3] = true;
					PrintInfo("Setting the second point");
					if (Par.isMin[3]) DrawDM(true);
				}
				break;
         case 4:
				if (Par.isMin[4])
				{
               if (CutMode.angledLine1X1.size() == CutMode.angledLine2X2.size())
               {
                  CutMode.angledLine1X1.push_back(Par.hist->GetXaxis()->
                     GetBinLowEdge(Par.hist->GetXaxis()->FindBin(x)));
                  CutMode.angledLine1Y1.push_back(Par.hist->GetYaxis()->
                     GetBinLowEdge(Par.hist->GetYaxis()->FindBin(y)));
                  
                  Par.isMin[4] = false;
                  PrintInfo("Setting the first point of the first line");	
               }
               else
               {
                  CutMode.angledLine2X1.push_back(Par.hist->GetXaxis()->
                     GetBinLowEdge(Par.hist->GetXaxis()->FindBin(x)));
                  CutMode.angledLine2Y1.push_back(Par.hist->GetYaxis()->
                     GetBinLowEdge(Par.hist->GetYaxis()->FindBin(y)));

                  Par.isMin[4] = false;
                  PrintInfo("Setting the second point of the first line");	
               }
				}
				else
				{
               if (CutMode.angledLine1X2.size() == CutMode.angledLine2X2.size())
               {
                  CutMode.angledLine1X2.push_back(Par.hist->GetXaxis()->GetBinLowEdge(
                     Par.hist->GetXaxis()->FindBin(x)));
                  CutMode.angledLine1Y2.push_back(Par.hist->GetYaxis()->GetBinLowEdge(
                     Par.hist->GetYaxis()->FindBin(y)));

                  CutMode.tanAlpha1.push_back((CutMode.angledLine1Y1.back() - 
                                               CutMode.angledLine1Y2.back())/
                                              (CutMode.angledLine1X1.back() - 
                                               CutMode.angledLine1X2.back()));
                  CutMode.shiftY1.push_back(CutMode.angledLine1Y2.back() - 
                                            CutMode.angledLine1X2.back()*CutMode.tanAlpha1.back());

                  Par.isMin[4] = true;
                  PrintInfo("Setting the second point of the first line");
                  if (Par.isMin[4]) DrawDM(true);
               }
               else
               {
                  CutMode.angledLine2X2.push_back(Par.hist->GetXaxis()->GetBinLowEdge(
                     Par.hist->GetXaxis()->FindBin(x)));
                  CutMode.angledLine2Y2.push_back(Par.hist->GetYaxis()->GetBinLowEdge(
                     Par.hist->GetYaxis()->FindBin(y)));

                  CutMode.tanAlpha2.push_back((CutMode.angledLine2Y1.back() - 
                                               CutMode.angledLine2Y2.back())/
                                              (CutMode.angledLine2X1.back() - 
                                               CutMode.angledLine2X2.back()));
                  CutMode.shiftY2.push_back(CutMode.angledLine2Y2.back() - 
                                            CutMode.angledLine2X2.back()*CutMode.tanAlpha2.back());
                   
                  Par.isMin[4] = true;
                  PrintInfo("Setting the second point of the second line");
                  if (Par.isMin[4]) DrawDM(true);
               }
				}
            break;
		}
	}
   else if (event == kKeyPress)
   {
      switch(px)
      {
         case 'u':
            switch (CutMode.currentCutMode)
            {
               case 0:
                  if (CutMode.rectXMin.size() != 0)
                  {
                     if (Par.isMin[0])
                     {
                        CutMode.rectXMin.pop_back();
                        CutMode.rectYMin.pop_back();
                        CutMode.rectXMax.pop_back();
                        CutMode.rectYMax.pop_back();

                        PrintInfo("Deleting last minimum and maximum points");
                     }
                     else
                     {
                        CutMode.rectXMin.pop_back();
                        CutMode.rectYMin.pop_back();
                           
                        PrintInfo("Deleting last minimum point");
                     }
                     if (Par.isMin[0]) DrawDM(true);
                     Par.isMin[0] = true;
                     gPad->Update();
                  }
                  else PrintInfo("Cannot delete last point since the current number of points is 0");
                  break;
               
               case 1:
                  if (CutMode.lineXMin.size() != 0)
                  {
                     if (Par.isMin[1])
                     {
                        CutMode.lineXMin.pop_back();
                        CutMode.lineXMax.pop_back();

                        PrintInfo("Deleting last pair of lines");
                     }
                     else
                     {
                        CutMode.lineXMin.pop_back();
                        PrintInfo("Deleting last line");
                     }
                     if (Par.isMin[1]) DrawDM(true);
                     Par.isMin[1] = true;
                     gPad->Update();
                  }
                  else PrintInfo("Cannot delete last line/lines since the current number of lines is 0");
                  break;
               case 2:
                  if (CutMode.lineYMin.size() != 0)
                  {
                     if (Par.isMin[2])
                     {
                        CutMode.lineYMin.pop_back();
                        CutMode.lineYMax.pop_back();

                        PrintInfo("Deleting last pair of lines");
                     }
                     else
                     {
                        CutMode.lineYMin.pop_back();
                        PrintInfo("Deleting last line");
                     }
                     if (Par.isMin[2]) DrawDM(true);
                     Par.isMin[2] = true;
                     gPad->Update();
                  }
                  else PrintInfo("Cannot delete last line/lines since the current number of lines is 0");
                  break;
               case 3:
                  if (CutMode.invRectXMin.size() != 0)
                  {
                     if (Par.isMin[3])
                     {
                        CutMode.invRectXMin.pop_back();
                        CutMode.invRectYMin.pop_back();
                        CutMode.invRectXMax.pop_back();
                        CutMode.invRectYMax.pop_back();

                        PrintInfo("Deleting last minimum and maximum points");
                     }
                     else
                     {
                        CutMode.invRectXMin.pop_back();
                        CutMode.invRectYMin.pop_back();
                           
                        PrintInfo("Deleting last minimum point");
                     }
                     if (Par.isMin[3]) DrawDM(true);
                     Par.isMin[3] = true;
                     gPad->Update();
                  }
                  else PrintInfo("Cannot delete last point since the current number of points is 0");
                  break;
               case 4:
                  if (CutMode.angledLine1X1.size() != 0)
                  {
                     if (CutMode.angledLine1X2.size() > CutMode.angledLine2X1.size())
                     {
                        if (Par.isMin[4])
                        {
                           CutMode.angledLine1X1.pop_back();
                           CutMode.angledLine1Y1.pop_back();
                           CutMode.angledLine1X2.pop_back();
                           CutMode.angledLine1Y2.pop_back();
                           CutMode.tanAlpha1.pop_back();
                           CutMode.shiftY1.pop_back();
                           
                           PrintInfo("Deleting last minimum and maximum points");
                        }
                        else
                        {
                           CutMode.angledLine1X1.pop_back();
                           CutMode.angledLine1Y1.pop_back();
                              
                           PrintInfo("Deleting last minimum point");
                        }
                        if (Par.isMin[4]) DrawDM(true);
                        Par.isMin[4] = true;
                        gPad->Update();
                     }
                     else
                     {
                        if (Par.isMin[4])
                        {
                           CutMode.angledLine2X1.pop_back();
                           CutMode.angledLine2Y1.pop_back();
                           CutMode.angledLine2X2.pop_back();
                           CutMode.angledLine2Y2.pop_back();
                           CutMode.tanAlpha2.pop_back();
                           CutMode.shiftY2.pop_back();

                           PrintInfo("Deleting last minimum and maximum points");
                        }
                        else
                        {
                           CutMode.angledLine2X1.pop_back();
                           CutMode.angledLine2Y1.pop_back();
                              
                           PrintInfo("Deleting last minimum point");
                        }
                     }
                     if (Par.isMin[4]) DrawDM(true);
                     Par.isMin[4] = true;
                     gPad->Update();
                  }
                  else PrintInfo("Cannot delete last point since the current number of points is 0");
                  break;
            }
            break;
            
         case 'r':
            DrawDM();
            PrintInfo("Resetting the range of the selected area");
            break;
            
         case 'p':
            for (int i = 0; i < CutMode.invRectXMax.size(); i++)
            {
               std::cout <<
                  "if (" + Par.xValName + " < " << CutMode.invRectXMin[i] <<
                  " || " + Par.yValName + " < " << CutMode.invRectYMin[i] <<
                  " ||" << std::endl << "    " + Par.xValName + " > " << CutMode.invRectXMax[i] <<
                  " || " + Par.yValName + " > " << CutMode.invRectYMax[i] <<
                  ") return true;" << std::endl;
            }
            for (int i = 0; i < CutMode.lineXMax.size(); i++)
            {
               std::cout <<
                  "if (" + Par.xValName + " > " << CutMode.lineXMin[i] <<
                  " && " + Par.xValName + " < " << CutMode.lineXMax[i] <<
                  ") return true;" << std::endl;
            }
            for (int i = 0; i < CutMode.lineYMin.size(); i++)
            {
               std::cout <<
                  "if (" + Par.yValName + " > " << CutMode.lineYMin[i] <<
                  " && " + Par.yValName + " < " << CutMode.lineYMin[i] <<
                  ") return true;" << std::endl;
            }
            for (int i = 0; i < CutMode.rectXMax.size(); i++)
            {
               std::cout << "if (" + Par.xValName + " > " << CutMode.rectXMin[i] <<
                  " && " + Par.yValName + " > " << CutMode.rectYMin[i] <<
                  " &&" << std::endl << "    " + Par.xValName + " < " << CutMode.rectXMax[i] <<
                  " && " + Par.yValName + " < " << CutMode.rectYMax[i] <<
                  ") return true;" << std::endl;
            }
            for (int i = 0; i < CutMode.shiftY2.size(); i++)
            {
               if ((CutMode.tanAlpha1[i] > 0 && CutMode.tanAlpha2[i] > 0) ||
                   (CutMode.tanAlpha1[i] < 0 && CutMode.tanAlpha2[i] < 0))
               {
                  std::cout << "if (" + Par.yValName + " > " << CutMode.shiftY1[i] << 
                               " + " + Par.xValName + "*" << CutMode.tanAlpha1[i] <<
                               " && " + Par.yValName + " < " << CutMode.shiftY2[i] << 
                               " + " + Par.xValName + "*" << CutMode.tanAlpha2[i] <<
                               ") return true;" <<std::endl;
               }
            }
            break;
            
         case '0':
            Print("Deactivating cutting mode");
            CutMode.currentCutMode = -1;
            break;
            
         case '1':
            Print("Activating rectangular cutting mode");
            CutMode.currentCutMode = 0;
            break;
            
         case '2':
            Print("Activating linear cutting mode along x axis");
            CutMode.currentCutMode = 1;
            break;
            
         case '3':
            Print("Activating linear cutting mode along y axis");
            CutMode.currentCutMode = 2;
            break;
            
         case '4':
            Print("Activating inverse rectangular cutting mode");
            CutMode.currentCutMode = 3;
            break;
            
         case '5':
            Print("Activating angled linear cutting mode");
            CutMode.currentCutMode = 4;
            break;
      }
   }
}

void GUIDM()
{
	gStyle->SetOptStat(0);
	Par.hist = (TH2F *) Par.file.Get("Heatmap: DCe1");
	Par.hist->SetTitle(Par.hist->GetName());

	Par.orig_integral = Par.hist->Integral();

   //DataCutsSelector dSel("Run14HeAu200MB");
	
	CutDeadAreas(Par.hist, &Run14HeAu200MBCuts::IsDeadDC, 1., -1.);
	//CutDeadAreas(Par.hist, &IsDeadPC1, 1.);
	//CutDeadAreas(Par.hist, &IsDeadPC2);
	//CutDeadAreas(Par.hist, &IsDeadPC3, 1.);
	//CutDeadAreas(Par.hist, &IsDeadEMCal, 2., 1., 3);
	//CutDeadAreas(Par.hist, &IsDeadTOFe, -1.);
	//CutDeadAreas(Par.hist, &IsDeadTOFw, -1.);

	TCanvas *canv = new TCanvas("canv");
   
	DrawDM();
	
	gPad->AddExec("exec", "exec()");
}
