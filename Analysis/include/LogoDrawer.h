#pragma once

namespace LogoDrawer
{
   int uniqueId = 0;
}

void DrawPHENIXLogoPreliminary(const double x, const double y, const double size = 1.)
{
   gPad->Update();
   
   const double xMax = x + (gPad->PixeltoX(static_cast<int>(100.*size)) - gPad->PixeltoX(0))/(gPad->GetX2() - gPad->GetX1());
   const double yMax = y + (gPad->PixeltoY(0) - gPad->PixeltoY(static_cast<int>(45.*size)))/(gPad->GetY2() - gPad->GetY1());

   std::cout << xMax << " " << yMax << std::endl;

   TPad *logo_pad = new TPad(("LogoDrawer_" + to_string(LogoDrawer::uniqueId)).c_str(), "", 
                             0., 0., 1., 1.);
   
   logo_pad->SetFillStyle(4000);
   LogoDrawer::uniqueId++;
   
   logo_pad->SetPad(x, y, xMax, yMax);
   
   logo_pad->Draw();
   logo_pad->cd();
   logo_pad->Update();
   
   TText normalText;
   TText boldText; 
   
   normalText.SetTextFont(43);
   boldText.SetTextFont(63);
   
   normalText.SetTextSize(static_cast<int>(20.*size));
   boldText.SetTextSize(static_cast<int>(20.*size));

   normalText.DrawText(0.05, 0.1, "preliminary");
   boldText.DrawText(0.066, 0.45, "PH");
   boldText.DrawText(0.51, 0.45, "ENIX");
   
   TLine line;
   if (size > 1.) line.SetLineWidth(static_cast<int>(size));
   line.SetLineColor(kGray+1);
   
   const double centerX = 0.4195;
   const double centerY = 0.615;
   
   const double lineDistanceFromCenterX = 0.018;
   const double lineDistanceFromCenterY = 0.04;
   
   const double lineLengthX = 0.095;
   const double lineLengthY = 0.211;
   
   //Drawing O
   for (double phi = M_PI/12; phi < M_PI*2.; phi += M_PI/6.)
   {
      line.DrawLineNDC(centerX+lineDistanceFromCenterX*cos(phi), 
                       centerY+lineDistanceFromCenterY*sin(phi), 
                       centerX+lineLengthX*cos(phi), 
                       centerY+lineLengthY*sin(phi));
   }
   
   //Drawing red thing
   TF1 leftPartOfRedThingFunction("leftPartOfRedThing", "pol2");
   TF1 rightPartOfRedThingFunction("rightPartOfRedThing", "pol2");
   
   TGraph grShade;
   grShade.SetFillStyle(1001);
   grShade.SetFillColor(kRed);
   
   for (double grX = 0.17; grX < centerX; grX += 0.0025)
   {
      grShade.AddPoint(grX, -3.23*grX*grX + 1.411*grX + 0.8447);
   }
   
   for (double grX = centerX; grX < 0.69; grX += 0.0026)
   {
      grShade.AddPoint(grX, sqrt(grX - 0.405)/1.7 - pow(grX + 0.6, 2)/4.3 + 1.0392);
   }

   for (double grX = 0.69; grX > centerX; grX -= 0.0026)
   {
      grShade.AddPoint(grX, sqrt(grX - 0.405)/1.3 - pow(grX + 0.6, 2)/3.51 + 1.0203);
   }
   
   for (double grX = centerX; grX > 0.17; grX -= 0.0025)
   {
      grShade.AddPoint(grX, -3.53*grX*grX + 1.411*grX + 0.8447);
   }
   
   grShade.DrawClone("F");
}
