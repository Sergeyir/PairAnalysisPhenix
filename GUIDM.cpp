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
	double origIntegral;
	std::array<bool, 5> isMin = {true, true, true, true, true};
	TFile realFile = TFile("data/Real/Run14HeAu200/sum.root");
	TFile simFile = TFile("data/PostSim/Run14HeAu200/Heatmaps/all.root");
   
	//TFile file = TFile("data/Real/Run15pp200/sum.root");
	TH2F *realHist;
   TH2F *realHistDM;
	TH2F *simHist;
   TH2F *simHistDM;

   bool useSimHist = true;
   bool isCurrentSim = false;

	std::string xValName;
	std::string yValName;
   const std::string tab = "   ";
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
	Par.realHistDM = (TH2F *) Par.realHist->Clone();
	if (Par.useSimHist) Par.simHistDM = (TH2F *) Par.simHist->Clone();
	
	for (int i = 0; i < Par.realHistDM->GetXaxis()->GetNbins(); i++)
	{
		for (int j = 0; j < Par.realHistDM->GetYaxis()->GetNbins(); j++)
		{
			double valX = Par.realHistDM->GetXaxis()->GetBinCenter(i);
			double valY = Par.realHistDM->GetYaxis()->GetBinCenter(j);
			
			for (int k = 0; k < CutMode.rectXMax.size(); k++)
			{
				if (valX > CutMode.rectXMin[k] && valX < CutMode.rectXMax[k] && 
					valY > CutMode.rectYMin[k] && valY < CutMode.rectYMax[k]) 
            {
					Par.realHistDM->SetBinContent(Par.realHistDM->GetBin(i, j), 0.);
					if (Par.useSimHist) Par.simHistDM->SetBinContent(Par.simHistDM->GetBin(i, j), 0.);
            }
			}
			for (int k = 0; k < CutMode.lineXMax.size(); k++)
			{
				if (valX > CutMode.lineXMin[k] && valX < CutMode.lineXMax[k]) 
            {
					Par.realHistDM->SetBinContent(Par.realHistDM->GetBin(i, j), 0.);
					if (Par.useSimHist) Par.simHistDM->SetBinContent(Par.simHistDM->GetBin(i, j), 0.);
            }
			}
			for (int k = 0; k < CutMode.lineYMin.size(); k++)
			{
				if (valY > CutMode.lineYMin[k] && valY < CutMode.lineYMax[k]) 
            {
					Par.realHistDM->SetBinContent(Par.realHistDM->GetBin(i, j), 0.);
					if (Par.useSimHist) Par.simHistDM->SetBinContent(Par.simHistDM->GetBin(i, j), 0.);
            }
			}
			for (int k = 0; k < CutMode.invRectXMax.size(); k++)
			{
				if (valX < CutMode.invRectXMin[k] || valX > CutMode.invRectXMax[k] ||
					valY < CutMode.invRectYMin[k] || valY > CutMode.invRectYMax[k]) 
            {
					Par.realHistDM->SetBinContent(Par.realHistDM->GetBin(i, j), 0.);
					if (Par.useSimHist) Par.simHistDM->SetBinContent(Par.simHistDM->GetBin(i, j), 0.);
            }
			}
			for (int k = 0; k < CutMode.shiftY2.size(); k++)
			{ 
            const double y1Cut = Pol1(CutMode.shiftY1[k], CutMode.tanAlpha1[k], valX);
            const double y2Cut = Pol1(CutMode.shiftY2[k], CutMode.tanAlpha2[k], valX);
            
				if (valY > Minimum(y1Cut, y2Cut) && valY < Maximum(y1Cut, y2Cut)) 
            {
					Par.realHistDM->SetBinContent(Par.realHistDM->GetBin(i, j), 0.);
					if (Par.useSimHist) Par.simHistDM->SetBinContent(Par.simHistDM->GetBin(i, j), 0.);
            }
			}
		}
	}

	if (is_fixed_range)
	{
		const int xMin = Par.realHistDM->GetXaxis()->FindBin(gPad->GetUxmin());
		const int xMax = Par.realHistDM->GetXaxis()->FindBin(gPad->GetUxmax())-1;
		const int yMin = Par.realHistDM->GetYaxis()->FindBin(gPad->GetUymin())+1;
		const int yMax = Par.realHistDM->GetYaxis()->FindBin(gPad->GetUymax())-1;
		
		Par.realHistDM->GetXaxis()->SetRange(xMin, xMax);
		Par.realHistDM->GetYaxis()->SetRange(yMin, yMax);
		if (Par.useSimHist) 
      {
         Par.simHistDM->GetXaxis()->SetRange(xMin, xMax);
		   Par.simHistDM->GetYaxis()->SetRange(yMin, yMax);
      }
	}
	
	if (!Par.isCurrentSim)
   {
      Par.realHistDM->Draw("COLZ");
   }
   else
   {
      Par.simHistDM->Draw("COLZ");
   }
	gPad->Modified();
	gPad->Update();
   gPad->GetFrame()->SetBit(TBox::kCannotMove);

	PrintInfo("Real data lost " + 
             to_string((1. - Par.realHistDM->Integral()/Par.origIntegral)*100.) + "\%");
   
   if (Par.useSimHist)
   {
      PrintInfo("Sim data lost " + 
                to_string((1. - Par.realHistDM->Integral()/Par.origIntegral)*100.) + "\%");
   }
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
         
         const double yMin = Par.realHist->GetYaxis()->GetBinLowEdge(1);
         yMax = Par.realHist->GetYaxis()->GetBinUpEdge(Par.realHist->GetYaxis()->GetNbins());
         
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
         
         const double xMin = Par.realHist->GetXaxis()->GetBinLowEdge(1);
         xMax = Par.realHist->GetXaxis()->GetBinUpEdge(Par.realHist->GetXaxis()->GetNbins());
         
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
         const double lineXMin = Par.realHist->GetXaxis()->GetBinLowEdge(1);
         const double lineXMax = 
            Par.realHist->GetXaxis()->GetBinUpEdge(Par.realHist->GetXaxis()->GetNbins());
         
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
						CutMode.rectXMax.push_back(Par.realHist->GetXaxis()->GetBinUpEdge(
							Par.realHist->GetXaxis()->FindBin(x)));
						CutMode.rectXMin.back() = Par.realHist->GetXaxis()->GetBinLowEdge(
							Par.realHist->GetXaxis()->FindBin(CutMode.rectXMin.back()));
					}
					else 
					{
						CutMode.rectXMax.push_back(Par.realHist->GetXaxis()->GetBinUpEdge(
							Par.realHist->GetXaxis()->FindBin(CutMode.rectXMin.back())));
						CutMode.rectXMin.back() = Par.realHist->GetXaxis()->GetBinLowEdge(
							Par.realHist->GetXaxis()->FindBin(x));
					}
					
					if (y >= CutMode.rectYMin.back()) 
					{
						CutMode.rectYMax.push_back(Par.realHist->GetYaxis()->GetBinUpEdge(
							Par.realHist->GetYaxis()->FindBin(y)));
						CutMode.rectYMin.back() = Par.realHist->GetYaxis()->GetBinLowEdge(
							Par.realHist->GetYaxis()->FindBin(CutMode.rectYMin.back()));
					}
					else 
					{
						CutMode.rectYMax.push_back(Par.realHist->GetYaxis()->GetBinUpEdge(
							Par.realHist->GetYaxis()->FindBin(CutMode.rectYMin.back())));
						CutMode.rectYMin.back() = Par.realHist->GetYaxis()->GetBinLowEdge(
							Par.realHist->GetYaxis()->FindBin(y));
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
						CutMode.lineXMax.push_back(Par.realHist->GetXaxis()->GetBinUpEdge(
							Par.realHist->GetXaxis()->FindBin(x)));
						CutMode.lineXMin.back() = Par.realHist->GetXaxis()->GetBinLowEdge(
							Par.realHist->GetXaxis()->FindBin(CutMode.lineXMin.back()));
					}
					else 
					{
						CutMode.lineXMax.push_back(Par.realHist->GetXaxis()->GetBinUpEdge(
							Par.realHist->GetXaxis()->FindBin(CutMode.lineXMin.back())));
						CutMode.lineXMin.back() = Par.realHist->GetXaxis()->GetBinLowEdge(
							Par.realHist->GetXaxis()->FindBin(x));
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
						CutMode.lineYMax.push_back(Par.realHist->GetYaxis()->GetBinUpEdge(
							Par.realHist->GetYaxis()->FindBin(y)));
						CutMode.lineYMin.back() = Par.realHist->GetYaxis()->GetBinLowEdge(
							Par.realHist->GetYaxis()->FindBin(CutMode.lineYMin.back()));
					}
					else 
					{
						CutMode.lineYMax.push_back(Par.realHist->GetYaxis()->GetBinUpEdge(
							Par.realHist->GetYaxis()->FindBin(CutMode.lineYMin.back())));
						CutMode.lineYMin.back() = Par.realHist->GetYaxis()->GetBinLowEdge(
							Par.realHist->GetYaxis()->FindBin(y));
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
						CutMode.invRectXMax.push_back(Par.realHist->GetXaxis()->GetBinLowEdge(
							Par.realHist->GetXaxis()->FindBin(x)));
						CutMode.invRectXMin.back() = Par.realHist->GetXaxis()->GetBinUpEdge(
							Par.realHist->GetXaxis()->FindBin(CutMode.invRectXMin.back()));
					}
					else 
					{
						CutMode.invRectXMax.push_back(Par.realHist->GetXaxis()->GetBinLowEdge(
							Par.realHist->GetXaxis()->FindBin(CutMode.invRectXMin.back())));
						CutMode.invRectXMin.back() = Par.realHist->GetXaxis()->GetBinUpEdge(
							Par.realHist->GetXaxis()->FindBin(x));
					}
					
					if (y >= CutMode.invRectYMin.back()) 
					{
						CutMode.invRectYMax.push_back(Par.realHist->GetYaxis()->GetBinLowEdge(
							Par.realHist->GetYaxis()->FindBin(y)));
						CutMode.invRectYMin.back() = Par.realHist->GetYaxis()->GetBinUpEdge(
							Par.realHist->GetYaxis()->FindBin(CutMode.invRectYMin.back()));
					}
					else 
					{
						CutMode.invRectYMax.push_back(Par.realHist->GetYaxis()->GetBinLowEdge(
							Par.realHist->GetYaxis()->FindBin(CutMode.invRectYMin.back())));
						CutMode.invRectYMin.back() = Par.realHist->GetYaxis()->GetBinUpEdge(
							Par.realHist->GetYaxis()->FindBin(y));
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
                  CutMode.angledLine1X1.push_back(Par.realHist->GetXaxis()->
                     GetBinLowEdge(Par.realHist->GetXaxis()->FindBin(x)));
                  CutMode.angledLine1Y1.push_back(Par.realHist->GetYaxis()->
                     GetBinLowEdge(Par.realHist->GetYaxis()->FindBin(y)));
                  
                  Par.isMin[4] = false;
                  PrintInfo("Setting the first point of the first line");	
               }
               else
               {
                  CutMode.angledLine2X1.push_back(Par.realHist->GetXaxis()->
                     GetBinLowEdge(Par.realHist->GetXaxis()->FindBin(x)));
                  CutMode.angledLine2Y1.push_back(Par.realHist->GetYaxis()->
                     GetBinLowEdge(Par.realHist->GetYaxis()->FindBin(y)));

                  Par.isMin[4] = false;
                  PrintInfo("Setting the second point of the first line");	
               }
				}
				else
				{
               if (CutMode.angledLine1X2.size() == CutMode.angledLine2X2.size())
               {
                  CutMode.angledLine1X2.push_back(Par.realHist->GetXaxis()->GetBinLowEdge(
                     Par.realHist->GetXaxis()->FindBin(x)));
                  CutMode.angledLine1Y2.push_back(Par.realHist->GetYaxis()->GetBinLowEdge(
                     Par.realHist->GetYaxis()->FindBin(y)));

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
                  CutMode.angledLine2X2.push_back(Par.realHist->GetXaxis()->GetBinLowEdge(
                     Par.realHist->GetXaxis()->FindBin(x)));
                  CutMode.angledLine2Y2.push_back(Par.realHist->GetYaxis()->GetBinLowEdge(
                     Par.realHist->GetYaxis()->FindBin(y)));

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
         {
            // Printing inverse box cut first
            for (int i = 0; i < CutMode.invRectXMax.size(); i++)
            {
               std::cout << Par.tab <<
                  "if (" + Par.xValName + " < " << CutMode.invRectXMin[i] <<
                  " || " + Par.yValName + " < " << CutMode.invRectYMin[i] <<
                  " ||" << std::endl << Par.tab << "    " + Par.xValName + " > " << 
                  CutMode.invRectXMax[i] << " || " + Par.yValName + 
                  " > " << CutMode.invRectYMax[i] << ") return true;" << std::endl;
            }
            
            // Determining bin ranges along X and Y axis
            int minXBin = 1;
            int minYBin = 1;
            int maxXBin = Par.realHistDM->GetXaxis()->GetNbins();
            int maxYBin = Par.realHistDM->GetYaxis()->GetNbins();

            if (CutMode.invRectXMax.size() != 0)
            {
               for (unsigned long i = 0; i < CutMode.invRectXMax.size(); i++)
               {
                  minXBin = Maximum(Par.realHist->GetXaxis()->FindBin(CutMode.invRectXMin[i]), 
                                    minXBin);
                  minYBin = Maximum(Par.realHist->GetYaxis()->FindBin(CutMode.invRectYMin[i]), 
                                    minYBin);
                  maxXBin = Minimum(Par.realHist->GetXaxis()->FindBin(CutMode.invRectXMax[i]), 
                                    maxXBin);
                  maxYBin = Minimum(Par.realHist->GetYaxis()->FindBin(CutMode.invRectYMax[i]), 
                                    maxYBin);
               }
               minXBin++;
               minYBin++;
               maxXBin--;
               maxYBin--;
            }
            else
            {
               for (int i = 1; i <= Par.realHistDM->GetXaxis()->GetNbins(); i++)
               {
                  bool isSelected = false;
                  for (int j = 1; j <= Par.realHistDM->GetYaxis()->GetNbins(); j++)
                  {
                     if (Par.realHist->GetBinContent(i, j) > 1e-7)
                     {
                        minXBin = i;
                        isSelected = true;
                        break;
                     }
                  }
                  if (isSelected) break;
               }
               for (int i = 1; i <= Par.realHistDM->GetYaxis()->GetNbins(); i++)
               {
                  bool isSelected = false;
                  for (int j = 1; j <= Par.realHistDM->GetXaxis()->GetNbins(); j++)
                  {
                     if (Par.realHist->GetBinContent(j, i) > 1e-7)
                     {
                        minYBin = i;
                        isSelected = true;
                        break;
                     }
                  }
                  if (isSelected) break;
               }
               for (int i = Par.realHistDM->GetXaxis()->GetNbins(); i >= minXBin; i--)
               {
                  bool isSelected = false;
                  for (int j = maxYBin; j >= minYBin; j--)
                  {
                     if (Par.realHist->GetBinContent(i, j) > 1e-7)
                     {
                        maxXBin = i;
                        isSelected = true;
                        break;
                     }
                  }
                  if (isSelected) break;
               }
               for (int i = Par.realHistDM->GetYaxis()->GetNbins(); i >= minYBin; i--)
               {
                  bool isSelected = false;
                  for (int j = maxXBin; j >= minXBin; j--)
                  {
                     if (Par.realHist->GetBinContent(j, i) > 1e-7)
                     {
                        maxYBin = i;
                        isSelected = true;
                        break;
                     }
                  }
                  if (isSelected) break;
               }
            }

            // Printing all other cuts
            std::cout << "switch(" << Par.yValName << "Bin)" << std::endl;
            std::cout << "{" << std::endl;
            
            for (int i = minYBin; i <= maxYBin; i++)
            {
               unsigned int numberOfEmptyBins = 0;
               unsigned int numberOfCuts = 0;
               
               for (int j = minXBin; j <= maxXBin; j++)
               {
                  if (Par.realHistDM->GetBinContent(j, i) < 1e-7) numberOfCuts++;
                  if (Par.realHist->GetBinContent(j, i) < 1e-7) numberOfEmptyBins++;
               }
               
               //if (numberOfEmptyBins == numberOfCuts) continue;
               
               if (numberOfCuts != maxXBin - minXBin + 1)
               {
                  std::cout << Par.tab << "case " << i << ": switch(" << Par.xValName << "Bin) {";

                  for (int j = minXBin; j <= maxXBin; j++)
                  {
                     if (Par.realHistDM->GetBinContent(j, i) < 1e-7)
                     {
                        std::cout << "case " << j << ": ";
                     }
                  }

                  std::cout << "return true; break;} ";
               }
               else 
               {
                  std::cout << Par.tab << "case " << i << ": return true; ";
               }
               std::cout << "break;" << std::endl;
            }
            
            std::cout << "}" << std::endl;
            break;
         }

         case 's':
         {
            if (Par.useSimHist)
            {
               Par.isCurrentSim = !Par.isCurrentSim;
               if (Par.isCurrentSim) PrintInfo("Showing sim histogram");
               else PrintInfo("Showing real data histogram");
            }
            else
            {
               PrintInfo("Cannot show sim histogram since the option for it was not specified");
            }
            DrawDM();
            break;
         }
            
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
   using namespace Run14HeAu200Cuts;
   
	gStyle->SetOptStat(0);
	Par.realHist = (TH2F *) Par.realFile.Get("Heatmap: DCe, zed>=0");
	if (Par.useSimHist) Par.simHist = (TH2F *) Par.simFile.Get("Heatmap: DCe, zed>=0");
	//Par.realHist->SetTitle(Par.realHist->GetName());

	Par.origIntegral = Par.realHist->Integral();

   //DataCutsSelector dSel("Run14HeAu200MB");

   //DC
	Par.xValName = "board"; Par.yValName = "alpha";
	CutDeadAreas(Par.realHist, &IsDeadDC, 2., 1.);
	if (Par.useSimHist) CutDeadAreas(Par.simHist, &IsDeadDC, 2., 1.);
   //PC1
	//Par.xValName = "pc1z"; Par.yValName = "pc1phi";
	//CutDeadAreas(Par.realHist, &IsDeadPC1, 1.);
   //PC2
	//Par.xValName = "pc2z"; Par.yValName = "pc2phi";
	//CutDeadAreas(Par.realHist, &IsDeadPC2);
   //PC3
	//Par.xValName = "pc3z"; Par.yValName = "pc3phi";
	//CutDeadAreas(Par.realHist, &IsDeadPC3, 1.);
   //EMCal
	//CutDeadAreas(Par.realHist, &IsDeadEMCal, 2., 1., 3);
   //TOFe
	//CutDeadAreas(Par.realHist, &IsDeadTOFe, -1.);
   //TOFw
	//CutDeadAreas(Par.realHist, &IsDeadTOFw, -1.);

	TCanvas *canv = new TCanvas("canv");
   
	DrawDM();
	
	gPad->AddExec("exec", "exec()");
}
