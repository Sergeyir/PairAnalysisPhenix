// $SOURCE$
//------------------------------------------------------------------------------------------------
//                DeadAreasCuts functions declarations and realisations
//------------------------------------------------------------------------------------------------
// DeadAreasCuts
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic set of functions for rejecting data in bad/dead areas of different detectors
 **/
//------------------------------------------------------------------------------------------------

#ifndef DEAD_AREAS_CUTS_H
#define DEAD_AREAS_CUTS_H

namespace Run7AuAu200Cuts
{
   bool IsDeadDC(const double phi, const double zed, const double board,
                 const double alpha) // Run7AuAu200
   {
      if (alpha > (0.276 * board - (0.2))) 	return true;
      if (alpha < (-0.276 * board + (0.2))) 	return true;
      if (alpha > (-0.242 * board + (19.2)))	return true;
      if (alpha < (0.242 * board + (-19.3)))	return true;	
      
      if (phi > M_PI/2. && zed >= 0)
      {
         if (!( (alpha<(-0.36*board+0.36*3.5)) || 
            (alpha>(-0.36*board+0.36*4)) ) && alpha > 0) return true;
         if (!( (alpha>(0.36*board-0.36*5.7)) || 
            (alpha<(0.36*board-0.36*6.2)) )) return true;
         if (!( (alpha>(0.36*board-0.36*20.5)) || 
            (alpha<(0.36*board-0.36*22.)) )) return true;
         if (!( (alpha<(-0.36*board+0.36*20.8)) || 
            (alpha>(-0.36*board+0.36*21.5)) )) return true;
         if (!( (alpha<(-10.*board+10.*20.6)) || 
            (alpha>(-10.*board+10.*21.7)) )) return true;
         if (!( (alpha>(0.355*board-0.36*24.7)) || 
            (alpha<(0.355*board-0.36*25.6)) )) return true;
         if (!( (alpha>(0.36*board-0.36*31.8)) || 
            (alpha<(0.36*board-0.36*36.)) )) return true;
         if (!( (alpha<(-0.4*board+0.36*40.6)) || 
            (alpha>(-0.4*board+0.36*45.2)) )) return true;
         if (!( (alpha<(-0.4*board+0.36*75.3)) || 
            (alpha>(-0.4*board+0.36*81.)) )) return true;
         if (!( (alpha<(-0.4*board+0.36*57.3)) || 
            (alpha>(-0.4*board+0.36*59.5)) ) && alpha < 0) return true;
      }		    
          
      if (phi > M_PI/2. && zed < 0) 
      {
         if (!( (alpha<(-0.36*board+0.36*3)) || 
            (alpha>(-0.36*board+0.36*3.5)) )) return true;
         if (!( (alpha>(0.4*board-0.36*12.)) || 
            (alpha<(0.4*board-0.36*13.2)) )) return true;
         if (!( (alpha<(-0.4*board+0.36*19.8)) || 
            (alpha>(-0.4*board+0.36*21)) )) return true;
         if (!( (alpha>(0.4*board-0.4*27.8)) || 
            (alpha<(0.4*board-0.4*28.7)) )) return true;
         if (!( (alpha>(0.36*board-0.36*31.8)) || 
            (alpha<(0.36*board-0.36*36)) )) return true;
         if (!( (alpha<(-0.4*board+0.36*43.1)) || 
            (alpha>(-0.4*board+0.36*44.8)) )) return true;
         if (!( (alpha<(-0.4*board+0.36*40.5)) || 
            (alpha>(-0.4*board+0.36*41.2)) )) return true;
         if (!( (alpha<(-0.45*board+0.45*49.5)) || 
            (alpha>(-0.45*board+0.45*50.3)) )) return true;
         if (!( (alpha>(0.45*board-0.45*49.3)) || 
            (alpha<(0.45*board-0.45*50.9)) )) return true;
         if (!( (alpha>(50.*board-50.*49.3)) || 
            (alpha<(50.*board-50.*50.9)) )) return true;
         if (!( (alpha<(-0.45*board+0.36*68.4)) || 
            (alpha>(-0.45*board+0.36*69.8)) )) return true;
         if (!( (alpha>(0.45*board-0.45*58.8)) || 
            (alpha<(0.45*board-0.45*59.5)) )) return true;
         if (!( (alpha<(-0.4*board+0.36*75.2)) || 
            (alpha>(-0.4*board+0.36*80.8)) )) return true;
      }
            
      if (phi < 1.5 && zed >= 0)
      {
         if (!( (alpha>(0.5*board-0.36*13)) || 
            (alpha<(0.5*board-0.36*15.8)) )) return true;
         if (!( (alpha>(0.38*board-0.36*17.4)) || 
            (alpha<(0.38*board-0.36*19.3)) )) return true;
         if (!( (alpha>(0.38*board-0.38*32.)) || 
            (alpha<(0.38*board-0.38*43.)) )) return true;
         if (!( (alpha<(-0.5*board+0.5*55.5)) || 
            (alpha>(-0.5*board+0.5*56)) )) return true;
         if (!( (alpha<(-0.5*board+0.5*57.5)) || 
            (alpha>(-0.5*board+0.5*58)) )) return true;
         if (!( (alpha<(-0.38*board+0.38*62.)) || 
            (alpha>(-0.38*board+0.38*63.3)) )) return true;
         if (!( (alpha<(-50.*board+49.5*62.)) || 
            (alpha>(-50.*board+50.5*63.3)) )) return true;
         
         if ( !(alpha > 10000*(board-76.5)) ) return true;
      }
         
      if (phi < 1.5 && zed < 0)
      {
         if (!( (alpha>(0.38*board-0.36*4)) || 
            (alpha<(0.38*board-0.36*4.5)) )) return true;
         if (!( (alpha>(0.38*board-0.36*19.5)) || 
            (alpha<(0.38*board-0.36*20.5)) )) return true;
         if (!( (alpha<(-0.38*board+0.36*19.5)) || 
            (alpha>(-0.38*board+0.36*20.5)) )) return true;
         if (!( (alpha<(-50.*board+50.*18.5)) || 
            (alpha>(-50.*board+50.*19.5)) )) return true;
         if (!( (alpha>(0.5*board-0.36*32.5)) || 
            (alpha<(0.5*board-0.36*33.8)) )) return true;
         if (!( (alpha<(-0.5*board+0.36*32)) || 
            (alpha>(-0.5*board+0.36*34)) )) return true;
         if (!( (alpha>(10000*(board-23.2))) || 
            (alpha<(10000*(board-24.5))) )) return true;
         if (!( (alpha>(0.38*board-0.38*31.2)) || 
            (alpha<(0.38*board-0.38*43.)) )) return true;
         if (!( (alpha>(0.38*board-0.38*62.)) || 
            (alpha<(0.38*board-0.38*63.3)) )) return true;
         if (!( (alpha>(0.5*board-0.5*53.5)) || 
            (alpha<(0.5*board-0.5*54)) )) return true;
         if (!( (alpha<(-0.38*board+0.38*62.)) || 
            (alpha>(-0.38*board+0.38*63.3)) )) return true;
         if (!( (alpha<(-50.*board+50.*62.)) || 
            (alpha>(-50.*board+50.*63.3)) )) return true;

         if (!( (alpha>(0.15*board-0.15*71.5)) || 
            (alpha<(0.15*board-0.15*73.5)) )) return true;

         if (alpha < (0.242 * board + (-19))) return true;	
      }
      return false;
   }

   bool IsDeadPC1(const double phi, const double pc1z, const double pc1phi) // Run7AuAu200
   {
      if (phi > M_PI/2.)
      {
         if (pc1z < -83 || pc1phi < 2.165
            || pc1z > 84 || pc1phi > 3.7) return true;
         if (pc1z > -5 && pc1z < 4) return true;
         if (pc1phi > 3.03 && pc1phi < 3.075) return true;
         if (pc1phi > 2.995 && pc1phi < 3.015) return true;
         if (pc1phi > 2.93 && pc1phi < 2.955) return true;
         if (pc1phi > 2.28 && pc1phi < 2.4) return true;
         if (pc1phi > 2.74 && pc1phi < 2.745) return true;
         if (pc1phi > 2.54 && pc1phi < 2.545) return true;
         if (pc1phi > 3.13 && pc1phi < 3.14) return true;
         if (pc1phi > 3.33 && pc1phi < 3.335) return true;
         if (pc1phi > 3.53 && pc1phi < 3.535) return true;
         if (pc1z < 0)
         {
            if (pc1phi > 2.95)
            {
               if (pc1phi > 3.64) return true;
               if (pc1phi > 3.345 && pc1phi < 3.39) return true;
               if (pc1z > -82 && pc1phi > 3.535
                  && pc1z < -74 && pc1phi < 3.57) return true;
               if (pc1z > -10 && pc1phi > 3.535
                  && pc1z < -4 && pc1phi < 3.57) return true;
               if (pc1z > -81 && pc1phi > 3.435
                  && pc1z < -78 && pc1phi < 3.45) return true;
               if (pc1z > -48 && pc1phi > 3.405
                  && pc1z < -43 && pc1phi < 3.45) return true;
               if (pc1z > -51 && pc1phi > 3.405
                  && pc1z < -44 && pc1phi < 3.425) return true;
               if (pc1z > -84 && pc1phi > 3.29
                  && pc1z < -79 && pc1phi < 3.33) return true;
               if (pc1z > -67 && pc1phi > 3.29
                  && pc1z < -65 && pc1phi < 3.3) return true;
               if (pc1z > -28 && pc1phi > 3.33
                  && pc1z < -24 && pc1phi < 3.35) return true;
               if (pc1z > -47 && pc1phi > 3.29
                  && pc1z < -41 && pc1phi < 3.335) return true;
               if (pc1z > -42 && pc1phi > 3.29
                  && pc1z < -39 && pc1phi < 3.305) return true;
               if (pc1z > -42 && pc1phi > 3.32
                  && pc1z < -39 && pc1phi < 3.335) return true;
               if (pc1z > -26 && pc1phi > 3.305
                  && pc1z < -23 && pc1phi < 3.32) return true;
               if (pc1z > -20 && pc1phi > 3.25
                  && pc1z < -15 && pc1phi < 3.295) return true;
               if (pc1z > -10 && pc1phi > 3.25
                  && pc1z < -5 && pc1phi < 3.295) return true;
               if (pc1z > -16 && pc1phi > 3.25
                  && pc1z < -9 && pc1phi < 3.265) return true;
               if (pc1z > -84 && pc1phi > 3.14
                  && pc1z < -80 && pc1phi < 3.175) return true;
               if (pc1z > -74 && pc1phi > 3.14
                  && pc1z < -70 && pc1phi < 3.175) return true;
               if (pc1z > -51 && pc1phi > 3.18
                  && pc1z < -46 && pc1phi < 3.205) return true;
               if (pc1z > -39 && pc1phi > 3.145
                  && pc1z < -36 && pc1phi < 3.155) return true;
               if (pc1z > -81 && pc1phi > 3.075
                  && pc1z < -70 && pc1phi < 3.1) return true;
               if (pc1z > -43 && pc1phi > 3.085
                  && pc1z < -40 && pc1phi < 3.1) return true;
               if (pc1z > -31 && pc1phi > 3.075
                  && pc1z < -19 && pc1phi < 3.1) return true;
               if (pc1z > -37 && pc1phi > 3.095
                  && pc1z < -34 && pc1phi < 3.105) return true;
               if (pc1z > -31 && pc1phi > 3.11
                  && pc1z < -30 && pc1phi < 3.12) return true;
            }
            else
            {
               if (pc1phi > 2.615 && pc1phi < 2.655) return true;
               if (pc1z > -65 && pc1phi < 2.2) return true;
               if (pc1z > -60 && pc1phi > 2.415
                  && pc1z < -49 && pc1phi < 2.47) return true;
               if (pc1z > -59 && pc1phi > 2.81
                  && pc1z < -54 && pc1phi < 2.83) return true;
               if (pc1z > -57 && pc1phi > 2.83
                  && pc1z < -54 && pc1phi < 2.855) return true;
               if (pc1z > -15 && pc1phi > 2.815
                  && pc1z < -12 && pc1phi < 2.86) return true;
               if (pc1z > -14 && pc1phi > 2.815
                  && pc1z < -10 && pc1phi < 2.84) return true;
               if (pc1z > -38 && pc1phi > 2.78
                  && pc1z < -36 && pc1phi < 2.795) return true;
               if (pc1z > -84 && pc1phi > 2.645
                  && pc1z < -81 && pc1phi < 2.665) return true;
               if (pc1z > -29 && pc1phi > 2.71
                  && pc1z < -26 && pc1phi < 2.72) return true;
               if (pc1z > -39 && pc1phi > 2.51
                  && pc1z < -36 && pc1phi < 2.525) return true;
               if (pc1z > -84 && pc1phi > 2.165
                  && pc1z < -79 && pc1phi < 2.195) return true;
               if (pc1z > -77 && pc1phi > 2.185
                  && pc1z < -74 && pc1phi < 2.225) return true;
               if (pc1z > -58 && pc1phi > 2.24
                  && pc1z < -56 && pc1phi < 2.25) return true;
               if (pc1z > -45 && pc1phi > 2.27
                  && pc1z < -43 && pc1phi < 2.285) return true;
               if (pc1z > -33 && pc1phi > 2.195
                  && pc1z < -29 && pc1phi < 2.235) return true;
            }
         }
         else
         {
            if (pc1phi > 2.95)
            {
               if (pc1phi < 3.005) return true;
               if (pc1phi > 3.635 && pc1phi < 3.645) return true;
               if (pc1z > 3 && pc1phi > 3.68
                  && pc1z < 7 && pc1phi < 3.705) return true;
               if (pc1z > 18 && pc1phi > 3.685
                  && pc1z < 23 && pc1phi < 3.7) return true;
               if (pc1z > 38 && pc1phi > 3.605
                  && pc1z < 51 && pc1phi < 3.655) return true;
               if (pc1z > 53 && pc1phi > 3.345
                  && pc1z < 56 && pc1phi < 3.37) return true;
               if (pc1z > 38 && pc1phi > 3.25
                  && pc1z < 43 && pc1phi < 3.29) return true;
               if (pc1z > 46 && pc1phi > 3.25
                  && pc1z < 50 && pc1phi < 3.285) return true;
               if (pc1z > 63 && pc1phi > 3.255
                  && pc1z < 68 && pc1phi < 3.3) return true;
               if (pc1z > 51 && pc1phi > 3.205
                  && pc1z < 54 && pc1phi < 3.22) return true;
               if (pc1z > 53 && pc1phi > 3.17
                  && pc1z < 56 && pc1phi < 3.19) return true;
               if (pc1z > 76 && pc1phi > 3.265
                  && pc1z < 78 && pc1phi < 3.28) return true;
               if (pc1z > 79 && pc1phi > 3.23
                  && pc1z < 83 && pc1phi < 3.245) return true;
               if (pc1z > 59 && pc1phi > 3.07
                  && pc1z < 71 && pc1phi < 3.1) return true;
               if (pc1z > 76 && pc1phi > 3.19
                  && pc1z < 78 && pc1phi < 3.2) return true;
            }
            else
            {
               if (pc1phi > 2.91) return true;
               if (pc1z < 38 && pc1phi < 2.195) return true;
               if (pc1phi > 2.675 && pc1phi < 2.715) return true;
               if (pc1z > 22 && pc1phi > 2.895
                  && pc1z < 30 && pc1phi < 2.92) return true;
               if (pc1z > 3 && pc1phi > 2.74
                  && pc1z < 10 && pc1phi < 2.78) return true;
               if (pc1z > 29 && pc1phi > 2.74
                  && pc1z < 33 && pc1phi < 2.775) return true;
               if (pc1z > 44 && pc1phi > 2.775
                  && pc1z < 51 && pc1phi < 2.815) return true;
               if (pc1z > 53 && pc1phi > 2.71
                  && pc1z < 56 && pc1phi < 2.745) return true;
               if (pc1z > 63 && pc1phi > 2.75
                  && pc1z < 69 && pc1phi < 2.765) return true;
               if (pc1z > 8 && pc1phi > 2.655
                  && pc1z < 12 && pc1phi < 2.68) return true;
               if (pc1z > 15 && pc1phi > 2.67
                  && pc1z < 18 && pc1phi < 2.68) return true;
               if (pc1z > 3 && pc1phi > 2.54
                  && pc1z < 10 && pc1phi < 2.58) return true;
               if (pc1z > 30 && pc1phi > 2.6
                  && pc1z < 33 && pc1phi < 2.615) return true;
               if (pc1z > 19 && pc1phi > 2.495
                  && pc1z < 31 && pc1phi < 2.545) return true;
               if (pc1z > 8 && pc1phi > 2.42
                  && pc1z < 20 && pc1phi < 2.47) return true;
               if (pc1z > 48 && pc1phi > 2.505
                  && pc1z < 54 && pc1phi < 2.545) return true;
               if (pc1z > 53 && pc1phi > 2.395
                  && pc1z < 58 && pc1phi < 2.42) return true;
               if (pc1z > 57 && pc1phi > 2.415
                  && pc1z < 60 && pc1phi < 2.425) return true;
               if (pc1z > 74 && pc1phi > 2.415
                  && pc1z < 81 && pc1phi < 2.47) return true;
               if (pc1z > 79 && pc1phi > 2.395
                  && pc1z < 85 && pc1phi < 2.44) return true;
               if (pc1z > 75 && pc1phi > 2.475
                  && pc1z < 84 && pc1phi < 2.495) return true;
               if (pc1z > 79 && pc1phi > 2.16
                  && pc1z < 85 && pc1phi < 2.195) return true;
               if (pc1z > 20 && pc1phi > 2.225
                  && pc1z < 23 && pc1phi < 2.235) return true;
               if (pc1z > 59 && pc1phi > 2.205
                  && pc1z < 63 && pc1phi < 2.22) return true;
               if (pc1z > 60 && pc1phi > 2.185
                  && pc1z < 64 && pc1phi < 2.2) return true;
               if (pc1z > 40 && pc1phi > 2.225
                  && pc1z < 43 && pc1phi < 2.24) return true;
            }
         }
      }
      else
      {
         if (pc1z < -84 || pc1phi < -0.555
            || pc1z > 84 || pc1phi > 0.97) return true;
         if (pc1z > -4 && pc1z < 4) return true;
         if (pc1phi > 0.045 && pc1phi < 0.285) return true;
         if (pc1phi > 0.63 && pc1phi < 0.69) return true;
         if (pc1phi > 0.79 && pc1phi < 0.795) return true;
         if (pc1phi > 0.59 && pc1phi < 0.595) return true;
         if (pc1phi > 0.395 && pc1phi < 0.4) return true;
         if (pc1phi > 0.005 && pc1phi < 0.01) return true;
         
         if (pc1z < 0)
         {
            if (pc1phi > 0.2)
            {
               if (pc1z > -33 && pc1phi > 0.95) return true;
               if (pc1phi > 0.835 && pc1phi < 0.875) return true;
               if (pc1z > -71 && pc1phi > 0.79
                  && pc1z < -59 && pc1phi < 0.84) return true;
               if (pc1z > -54 && pc1phi > 0.82
                  && pc1z < -49 && pc1phi < 0.84) return true;
               if (pc1z > -40 && pc1phi > 0.825
                  && pc1z < -28 && pc1phi < 0.84) return true;
               if (pc1z > -81 && pc1phi > 0.69
                  && pc1z < -75 && pc1phi < 0.705) return true;
               if (pc1z > -78 && pc1phi > 0.695
                  && pc1z < -74 && pc1phi < 0.715) return true;
               if (pc1z > -36 && pc1phi > 0.595
                  && pc1z < -34 && pc1phi < 0.61) return true;
               if (pc1z > -30 && pc1phi > 0.555
                  && pc1z < -19 && pc1phi < 0.595) return true;
               if (pc1z > -71 && pc1phi > 0.4
                  && pc1z < -59 && pc1phi < 0.44) return true;
               if (pc1z > -30 && pc1phi > 0.36
                  && pc1z < -18 && pc1phi < 0.44) return true;
               if (pc1z > -15 && pc1phi > 0.44
                  && pc1z < -13 && pc1phi < 0.46) return true;
               if (pc1z > -17 && pc1phi > 0.555
                  && pc1z < -14 && pc1phi < 0.57) return true;
               if (pc1z > -18 && pc1phi > 0.525
                  && pc1z < -15 && pc1phi < 0.54) return true;
               if (pc1z > -43 && pc1phi > 0.5
                  && pc1z < -38 && pc1phi < 0.51) return true;
               if (pc1z > -38 && pc1phi > 0.51
                  && pc1z < -35 && pc1phi < 0.52) return true;
            }
            else
            {
               if (pc1phi > -0.11 && pc1phi < -0.085) return true;
               if (pc1phi > -0.21 && pc1phi < -0.18) return true;
               if (pc1phi > -0.505 && pc1phi < -0.495) return true;
               if (pc1phi > -0.485 && pc1phi < -0.48) return true;
               if (pc1z > -10 && pc1phi > 0.005
                  && pc1z < -3 && pc1phi < 0.05) return true;
               if (pc1z > -30 && pc1phi > -0.39
                  && pc1z < -19 && pc1phi < -0.345) return true;
               if (pc1z > -46 && pc1phi > -0.125
                  && pc1z < -39 && pc1phi < -0.1) return true;
               if (pc1z > -44 && pc1phi > -0.135
                  && pc1z < -42 && pc1phi < -0.125) return true;
               if (pc1z > -41 && pc1phi > -0.155
                  && pc1z < -39 && pc1phi < -0.145) return true;
               if (pc1z > -63 && pc1phi > -0.02
                  && pc1z < -60 && pc1phi < -0.01) return true;
               if (pc1z > -84 && pc1phi > -0.155
                  && pc1z < -81 && pc1phi < -0.145) return true;
               if (pc1z > -83 && pc1phi > -0.18
                  && pc1z < -80 && pc1phi < -0.165) return true;
               if (pc1z > -66 && pc1phi > -0.53
                  && pc1z < -63 && pc1phi < -0.515) return true;
               if (pc1z > -25 && pc1phi > -0.54
                  && pc1z < -23 && pc1phi < -0.53) return true;
               if (pc1z > -16 && pc1phi > -0.52
                  && pc1z < -14 && pc1phi < -0.505) return true;
               if (pc1z > -46 && pc1phi > -0.405
                  && pc1z < -44 && pc1phi < -0.395) return true;
            }
         }
         else
         {
            if (pc1phi > 0.2)
            {
               if (pc1phi > 0.93) return true;
               if (pc1z > 60 && pc1phi > 0.905) return true;
               if (pc1phi > 0.535 && pc1phi < 0.555) return true;
               if (pc1phi > 0.87 && pc1z < 10) return true;
               if (pc1z > 13 && pc1phi > 0.515
                  && pc1z < 15 && pc1phi < 0.56) return true;
               if (pc1z > 69 && pc1phi > 0.355
                  && pc1z < 81 && pc1phi < 0.44) return true;
               if (pc1z > 54 && pc1phi > 0.295
                  && pc1z < 57 && pc1phi < 0.31) return true;
               if (pc1z > 55 && pc1phi > 0.305
                  && pc1z < 62 && pc1phi < 0.325) return true;
               if (pc1z > 44 && pc1phi > 0.28
                  && pc1z < 47 && pc1phi < 0.295) return true;
               if (pc1z > 14 && pc1phi > 0.735
                  && pc1z < 19 && pc1phi < 0.76) return true;
               if (pc1z > 14 && pc1phi > 0.715
                  && pc1z < 16 && pc1phi < 0.745) return true;
               if (pc1z > 73 && pc1phi > 0.69
                  && pc1z < 75 && pc1phi < 0.7) return true;
               if (pc1z > 75 && pc1phi > 0.75
                  && pc1z < 78 && pc1phi < 0.765) return true;
            }
            else
            {
               if (pc1phi < -0.54) return true;
               if (pc1z > 7 && pc1phi > -0.265
                  && pc1z < 8 && pc1phi < -0.235) return true;
               if (pc1phi > -0.395 && pc1phi < -0.34) return true;
               if (pc1phi > -0.25 && pc1phi < -0.205) return true;
               if (pc1z > 49 && pc1phi > 0.005
                  && pc1z < 61 && pc1phi < 0.05) return true;
               if (pc1z > 77 && pc1phi > -0.03
                  && pc1z < 79 && pc1phi < -0.015) return true;
               if (pc1z > 24 && pc1phi > -0.52
                  && pc1z < 26 && pc1phi < -0.5) return true;
               if (pc1z > 27 && pc1phi > -0.545
                  && pc1z < 30 && pc1phi < -0.53) return true;
               if (pc1z > 2 && pc1phi > -0.275
                  && pc1z < 7 && pc1phi < -0.215) return true;
               if (pc1z > 63 && pc1phi > -0.3
                  && pc1z < 66 && pc1phi < -0.28) return true;
            }
         }
      }
      return false;
   }

   bool IsDeadPC2(const double pc2z, const double pc2phi) // Run7AuAu200
   {
      if (pc2z < -136 || pc2phi < -0.555
         || pc2z > 140 || pc2phi > 0.965) return true;
      if (pc2z > -12 && pc2z < 14) return true;
      if (pc2phi > 0.57 && pc2phi < 0.625) return true;
      if (pc2phi > 0.04 && pc2phi < 0.275) return true;
      if (pc2phi > -0.22 && pc2phi < -0.165) return true;
      if (pc2z < 0)
      {
         if (pc2phi > 0.615)
         {
            if (pc2z > -88 && pc2phi > 0.81
               && pc2z < -76 && pc2phi < 0.87) return true;
            if (pc2z > -82 && pc2z < -79) return true;
            if (pc2z < -82)
            {
               if (pc2z > -124 && pc2phi > 0.87
                  && pc2z < -84 && pc2phi < 0.95) return true;
               if (pc2z > -124 && pc2phi > 0.785
                  && pc2z < -110 && pc2phi < 0.835) return true;
               if (pc2z > -118 && pc2phi > 0.73
                  && pc2z < -104 && pc2phi < 0.76) return true;
               if (pc2z > -126 && pc2phi > 0.67
                  && pc2z < -102 && pc2phi < 0.73) return true;
               if (pc2z > -92 && pc2phi > 0.65
                  && pc2z < -88 && pc2phi < 0.66) return true;
            }
            else
            {
               if (pc2z > -58 && pc2phi > 0.885) return true;
               if (pc2z > -56 && pc2z < -36 && pc2phi < 0.68) return true;
               if (pc2z > -72 && pc2phi > 0.875
                  && pc2z < -58 && pc2phi < 0.955) return true;
               if (pc2z > -56 && pc2phi > 0.82
                  && pc2z < -34 && pc2phi < 0.87) return true;
               if (pc2z > -74 && pc2phi > 0.76
                  && pc2z < -62 && pc2phi < 0.795) return true;
               if (pc2z > -40 && pc2phi > 0.765
                  && pc2z < -24 && pc2phi < 0.795) return true;
               if (pc2z > -74 && pc2phi > 0.635
                  && pc2z < -52 && pc2phi < 0.685) return true;
               if (pc2z > -56 && pc2phi > 0.705
                  && pc2z < -42 && pc2phi < 0.755) return true;
            }
         }
         else if (pc2phi > 0.275)
         {
            if (pc2z > -84 && pc2z < -82) return true;
            if (pc2z > -108 && pc2phi > 0.47
               && pc2z < -76 && pc2phi < 0.52) return true;
            if (pc2z > -100 && pc2phi > 0.385
               && pc2z < -72 && pc2phi < 0.45) return true;
            if (pc2z < -83)
            {
               if (pc2phi > 0.35 && pc2z < -128 && pc2phi < 0.375) return true;
               if (pc2z < -126 && pc2phi < 0.34) return true;
            }
            else
            {
               if (pc2phi < 0.315) return true;
               if (pc2phi > 0.405 && pc2phi < 0.44) return true;
               if (pc2z > -40 && pc2phi > 0.5 && pc2z < -24) return true;
               if (pc2z > -22 && pc2phi > 0.475 && pc2phi < 0.52) return true;
               if (pc2z > -54 && pc2phi > 0.39 && pc2phi < 0.455) return true;
               if (pc2z > -60 && pc2phi < 0.37) return true;
               if (pc2z > -74 && pc2phi > 0.505
                  && pc2z < -52 && pc2phi < 0.555) return true;
            }
         }
         else if (pc2phi > -0.165)
         {
            if (pc2z > -82 && pc2z < -80) return true;
            if (pc2z > -92 && pc2phi > -0.065
               && pc2z < -70 && pc2phi < -0.015) return true;
            if (pc2z < -81)
            {
               if (pc2z > -114 && pc2phi > -0.07
                  && pc2z < -96 && pc2phi < -0.025) return true;
               if (pc2z > -100 && pc2phi > -0.145
                  && pc2z < -86 && pc2phi < -0.095) return true;
            }
            else
            {
               if (pc2z > -76 && pc2phi > 0 && pc2z < -54) return true;
               if (pc2z > -30 && pc2phi > -0.06
                  && pc2z < -18 && pc2phi < -0.025) return true;
            }
         }
         else
         {
            if (pc2z > -82 && pc2z < -80) return true;
            if (pc2z < -81)
            {
               if (pc2phi > -0.28 && pc2z < -114) return true;
               if (pc2phi > -0.38 && pc2z < -124) return true;
               if (pc2z < -112 && pc2phi < -0.435) return true;
               if (pc2phi > -0.34 && pc2z < -108 && pc2phi < -0.275) return true;
               if (pc2z > -94 && pc2phi < -0.445) return true;
               if (pc2z > -102 && pc2phi > -0.35
                  && pc2z < -94 && pc2phi < -0.3) return true;
               if (pc2z > -108 && pc2phi > -0.46
                  && pc2z < -98 && pc2phi < -0.415) return true;
               if (pc2z > -96 && pc2phi > -0.265
                  && pc2z < -88 && pc2phi < -0.245) return true;
               if (pc2z > -122 && pc2phi > -0.41
                  && pc2z < -112 && pc2phi < -0.385) return true;
            }
            else
            {
               if (pc2phi > -0.27 && pc2z < -68 && pc2phi < -0.23) return true;
               if (pc2z > -76 && pc2phi > -0.315
                  && pc2z < -56 && pc2phi < -0.26) return true;
               if (pc2z > -66 && pc2phi > -0.345
                  && pc2z < -60 && pc2phi < -0.305) return true;
               if (pc2z > -74 && pc2phi > -0.495
                  && pc2z < -60 && pc2phi < -0.44) return true;
               if (pc2z > -46 && pc2phi > -0.32
                  && pc2z < -34 && pc2phi < -0.275) return true;
               if (pc2z > -34 && pc2phi > -0.395
                  && pc2z < -26 && pc2phi < -0.38) return true;
               
               if (pc2phi > -0.535 && pc2z < -70 && pc2phi < -0.49) return true;
               if (pc2z > -48 && pc2z < -36 && pc2phi < -0.515) return true;
               if (pc2z > -24 && pc2phi < -0.49) return true;
               if (pc2z > -22 && pc2phi > -0.455 && pc2phi < -0.375) return true;
               if (pc2z > -16 && pc2phi > -0.27) return true;
            }
         }
      }
      else	
      {
         if (pc2phi > 0.6)
         {
            if (pc2z > 83 && pc2z < 86) return true;
            if (pc2phi > 0.925 && pc2z < 22) return true;
            if (pc2z > 76 && pc2phi > 0.68
               && pc2z < 90 && pc2phi < 0.835) return true;
            if (pc2z < 85)
            {
               if (pc2z > 12 && pc2phi > 0.64
                  && pc2z < 26 && pc2phi < 0.68) return true;
               if (pc2z > 22 && pc2phi > 0.75
                  && pc2z < 32 && pc2phi < 0.795) return true;
               if (pc2z > 38 && pc2phi > 0.79
                  && pc2z < 48 && pc2phi < 0.82) return true;
               if (pc2z > 38 && pc2phi > 0.635
                  && pc2z < 50 && pc2phi < 0.675) return true;
               if (pc2z > 56 && pc2phi > 0.675
                  && pc2z < 76 && pc2phi < 0.73) return true;
               if (pc2z > 30 && pc2phi > 0.83
                  && pc2z < 42 && pc2phi < 0.87) return true;
               if (pc2z > 64 && pc2phi > 0.94
                  && pc2z < 78 && pc2phi < 0.975) return true;
            }
            else
            {
               if (pc2z > 102 && pc2phi > 0.865
                  && pc2z < 114 && pc2phi < 0.91) return true;
               if (pc2z > 92 && pc2phi > 0.75
                  && pc2z < 112 && pc2phi < 0.795) return true;
               if (pc2z > 96 && pc2phi > 0.79
                  && pc2z < 108 && pc2phi < 0.83) return true;
            }
         }
         else if (pc2phi > 0.2)
         {
            if (pc2z > 82 && pc2z < 84) return true;
            if (pc2z > 70 && pc2phi > 0.36
               && pc2z < 96 && pc2phi < 0.405) return true;
            if (pc2z < 83)
            {
               if (pc2z > 22 && pc2phi > 0.28
                  && pc2z < 38 && pc2phi < 0.325) return true;
               if (pc2z > 22 && pc2phi > 0.47
                  && pc2z < 42 && pc2phi < 0.52) return true;
               if (pc2z > 56 && pc2phi > 0.46
                  && pc2z < 74 && pc2phi < 0.525) return true;
               if (pc2z > 60 && pc2phi > 0.4
                  && pc2z < 68 && pc2phi < 0.435) return true;
               if (pc2z > 56 && pc2phi > 0.325
                  && pc2z < 76 && pc2phi < 0.37) return true;
               if (pc2z > 40 && pc2z < 50 && pc2phi < 0.295) return true;
               if (pc2z > 16 && pc2phi > 0.4
                  && pc2z < 22 && pc2phi < 0.42) return true;
               if (pc2z > 28 && pc2phi > 0.39
                  && pc2z < 34 && pc2phi < 0.405) return true;
            }
            else
            {
               if (pc2z > 98 && pc2phi < 0.335) return true;
               if (pc2phi > 0.545 && pc2z < 94) return true;
               if (pc2z > 132 && pc2phi > 0.55) return true;
            }
         }
         else if (pc2phi > -0.2)
         {
            if (pc2z > 82 && pc2z < 84) return true;
            if (pc2z > 70 && pc2phi > -0.07
               && pc2z < 94 && pc2phi < -0.02) return true;
            if (pc2z < 83)
            {
               if (pc2z > 36 && pc2phi > -0.155
                  && pc2z < 56 && pc2phi < -0.11) return true;
               if (pc2z > 54 && pc2z < 74 && pc2phi < -0.145) return true;
               if (pc2z > 20 && pc2z < 50 && pc2phi < -0.135) return true;
               if (pc2phi > 0 && pc2z < 20 && pc2phi < 0.01) return true;
               if (pc2z > 56 && pc2phi > 0.03 && pc2z < 66) return true;
               if (pc2z > 40 && pc2phi > -0.105
                  && pc2z < 46 && pc2phi < -0.075) return true;
               if (pc2z > 72 && pc2phi > -0.125
                  && pc2z < 80 && pc2phi < -0.11) return true;
            }
            else
            {
               if (pc2phi > -0.025 && pc2z < 106 && pc2phi < 0.01) return true;
               if (pc2z > 116 && pc2phi > 0) return true;
               if (pc2z > 130 && pc2phi < -0.135) return true;
            }
         }
         else
         {
            if (pc2z > 82 && pc2z < 84) return true;
            if (pc2phi > -0.33 && pc2phi < -0.26) return true;
            else if (pc2z < 83)
            {
               if (pc2z > 8 && pc2phi > -0.27
                  && pc2z < 26 && pc2phi < -0.235) return true;
               if (pc2z > 8 && pc2phi > -0.27
                  && pc2z < 26 && pc2phi < -0.235) return true;
               if (pc2z > 12 && pc2phi > -0.535
                  && pc2z < 26 && pc2phi < -0.49) return true;
               if (pc2z > 24 && pc2phi > -0.415
                  && pc2z < 36 && pc2phi < -0.38) return true;
               if (pc2z > 32 && pc2phi > -0.39
                  && pc2z < 42 && pc2phi < -0.36) return true;
               if (pc2z > 40 && pc2phi > -0.565
                  && pc2z < 50 && pc2phi < -0.515) return true;
               if (pc2z > 44 && pc2phi > -0.49
                  && pc2z < 52 && pc2phi < -0.465) return true;
               if (pc2z > 70 && pc2phi > -0.36) return true;
            }
            else
            {
               if (pc2z > 120 && pc2phi < -0.52) return true;
               if (pc2z > 90 && pc2phi > -0.5
                  && pc2z < 110 && pc2phi < -0.45) return true;
            }
         }
      }
      
      return false;
   }

   bool IsDeadPC3(const double phi, const double pc3z, const double pc3phi) // Run7AuAu200
   {
      if (phi > M_PI/2.)
      {
         if (pc3z < -166 || pc3phi < 2.17
            || pc3z > 162 || pc3phi > 3.705) return true;
         if (pc3z > -12 && pc3z < 12) return true;
         if (pc3phi > 3.315 && pc3phi < 3.345) return true;
         if (pc3phi > 2.925 && pc3phi < 2.95) return true;
         if (pc3phi > 2.53 && pc3phi < 2.56) return true;
         if (pc3phi > 2.195 && pc3phi < 2.27) return true;
         if (pc3z < 0)
         {
            if (pc3phi > 3.32)
            {
               if (pc3z > -96 && pc3z < -94) return true;
               if (pc3z > -108 && pc3phi > 3.515
                  && pc3z < -80 && pc3phi < 3.57) return true;
               if (pc3z < -95)
               {
                  if (pc3z > -162 && pc3phi > 3.44
                     && pc3z < -154 && pc3phi < 3.46) return true;
               }
               else
               {
                  if (pc3z > -18 && pc3phi > 3.67) return true;
                  if (pc3z > -26 && pc3phi < 3.38) return true;
                  if (pc3z > -26 && pc3phi > 3.595 && pc3phi < 3.645) return true;
                  if (pc3z > -26 && pc3phi > 3.445 && pc3phi < 3.495) return true;
                  if (pc3z > -66 && pc3z < -42 && pc3phi < 3.38) return true;
                  if (pc3z > -68 && pc3phi > 3.41
                     && pc3z < -42 && pc3phi < 3.46) return true;
                  if (pc3z > -42 && pc3phi > 3.415
                     && pc3z < -22 && pc3phi < 3.46) return true;
               }
            }
            else if (pc3phi > 2.947)
            {
               if (pc3z > -106 && pc3z < -102) return true;
               if (pc3z > -132 && pc3phi > 3.13
                  && pc3z < -96 && pc3phi < 3.175) return true;
               if (pc3z < -105)
               {
                  if (pc3z > -156 && pc3z < -146 && pc3phi > 3.285) return true;
                  if (pc3phi > 3.135 && pc3z < -142 && pc3phi < 3.175) return true;
                  if (pc3phi > 3.02 && pc3z < -162 && pc3phi < 3.06) return true;
                  if (pc3z > -128 && pc3phi > 2.97 && pc3phi < 3.035) return true;
                  if (pc3z > -130 && pc3phi > 3.065
                     && pc3z < -116 && pc3phi < 3.08) return true;
               }
               else
               {
                  if (pc3z > -86 && pc3phi > 3.135
                     && pc3z < -62 && pc3phi < 3.175) return true;
                  if (pc3z > -34 && pc3phi > 3.18
                     && pc3z < -18 && pc3phi < 3.205) return true;
               }
            }
            else if (pc3phi > 2.557)
            {
               if (pc3z > -98 && pc3z < -94) return true;
               if (pc3phi > 2.65 && pc3phi < 2.71) return true;
               if (pc3z < -97)
               {
                  if (pc3phi > 2.89 && pc3z < -162) return true;
                  if (pc3z > -146 && pc3phi > 2.885 && pc3z < -102) return true;
                  if (pc3phi > 2.625 && pc3z < -140 && pc3phi < 2.655) return true;
                  if (pc3z > -126 && pc3phi > 2.585
                     && pc3z < -104 && pc3phi < 2.655) return true;
               }
               else
               {
                  if (pc3z > -88 && pc3phi > 2.885 && pc3z < -42) return true;
                  if (pc3z > -30 && pc3phi > 2.89) return true;
                  if (pc3z > -48 && pc3phi > 2.705
                     && pc3z < -22 && pc3phi < 2.75) return true;
                  if (pc3z > -80 && pc3phi > 2.79
                     && pc3z < -64 && pc3phi < 2.815) return true;
               }
            }
            else
            {
               if (pc3z > -90 && pc3z < -86) return true;
               if (pc3z > -108 && pc3phi > 2.42
                  && pc3z < -60 && pc3phi < 2.47) return true;
               if (pc3z > -108 && pc3phi > 2.455
                  && pc3z < -80 && pc3phi < 2.51) return true;
               if (pc3z < -89)
               {
                  if (pc3z > -148 && pc3z < -120 && pc3phi > 2.27) return true;
                  if (pc3phi > 2.31 && pc3z < -160 && pc3phi < 2.51) return true;
                  if (pc3phi > 2.345 && pc3z < -144 && pc3phi < 2.435) return true;
                  if (pc3z > -116 && pc3phi > 2.31
                     && pc3z < -106 && pc3phi < 2.35) return true;
               }
               else
               {
                  if (pc3z > -46 && pc3phi > 2.495 && pc3z < -34) return true;
                  if (pc3z > -26 && pc3phi > 2.42 && pc3phi < 2.47) return true;
                  if (pc3z > -66 && pc3phi > 2.205
                     && pc3z < -48 && pc3phi < 2.245) return true;
                  if (pc3z > -78 && pc3phi > 2.17
                     && pc3z < -64 && pc3phi < 2.2) return true;
               }
            }
         }
         else
         {
            if (pc3phi > 3.32)
            {
               if (pc3z > 96 && pc3z < 100) return true;
               if (pc3z > 82 && pc3phi > 3.525
                  && pc3z < 148 && pc3phi < 3.57) return true;
               if (pc3z > 80 && pc3phi > 3.37
                  && pc3z < 106 && pc3phi < 3.46) return true;
               if (pc3z < 97)
               {
                  if (pc3phi > 3.41 && pc3phi < 3.46) return true;
                  if (pc3z > 22 && pc3z < 46 && pc3phi < 3.42) return true;
                  if (pc3z > 22 && pc3phi > 3.635
                     && pc3z < 48 && pc3phi < 3.68) return true;
                  if (pc3z > 30 && pc3phi > 3.565
                     && pc3z < 38 && pc3phi < 3.58) return true;
                  if (pc3z > 38 && pc3phi > 3.525
                     && pc3z < 66 && pc3phi < 3.57) return true;
               }
               else
               {
                  if (pc3z > 124 && pc3phi > 3.67 && pc3z < 146) return true;
                  if (pc3phi > 3.6 && pc3z < 118 && pc3phi < 3.65) return true;
                  if (pc3z > 142 && pc3phi > 3.595 && pc3phi < 3.65) return true;
                  if (pc3z > 122 && pc3phi > 3.485 && pc3phi < 3.53) return true;
                  if (pc3z > 140 && pc3phi > 3.37 && pc3phi < 3.46) return true;
                  if (pc3z > 122 && pc3phi > 3.37 && pc3phi < 3.42) return true;
                  if (pc3z > 136 && pc3phi > 3.575
                     && pc3z < 146 && pc3phi < 3.595) return true;
               }
            }
            else if (pc3phi > 2.947)
            {
               if (pc3phi > 3 && pc3phi < 3.035) return true;
               if (pc3z > 102 && pc3z < 106) return true;
               if (pc3z > 82 && pc3phi > 3.16 && pc3phi < 3.215) return true;
               if (pc3z > 82 && pc3phi > 3.24
                  && pc3z < 108 && pc3phi < 3.29) return true;
               if (pc3z < 103)
               {
                  if (pc3phi > 3.2 && pc3z < 32 && pc3phi < 3.255) return true;
                  if (pc3phi > 3.135 && pc3z < 26 && pc3phi < 3.175) return true;
                  if (pc3phi > 2.97 && pc3z < 68 && pc3phi < 3.015) return true;
                  if (pc3z > 42 && pc3z < 68 && pc3phi < 2.98) return true;
                  if (pc3z > 42 && pc3phi > 3.205
                     && pc3z < 52 && pc3phi < 3.24) return true;
                  if (pc3z > 40 && pc3phi > 3.095
                     && pc3z < 68 && pc3phi < 3.14) return true;
                  if (pc3z > 68 && pc3phi > 3.105
                     && pc3z < 76 && pc3phi < 3.12) return true;
               }
               else
               {
                  if (pc3phi > 3.13 && pc3phi < 3.18) return true;
                  if (pc3z > 142 && pc3phi > 3.09 && pc3phi < 3.145) return true;
                  if (pc3z > 122 && pc3phi > 2.975 && pc3phi < 3.035) return true;
                  if (pc3z > 116 && pc3z < 126 && pc3phi < 2.97) return true;
                  if (pc3z > 120 && pc3phi > 3.03
                     && pc3z < 146 && pc3phi < 3.105) return true;
                  if (pc3z > 122 && pc3phi > 3.215
                     && pc3z < 150 && pc3phi < 3.29) return true;
               }
            }
            else if (pc3phi > 2.557)
            {
               if (pc3phi > 2.86 && pc3phi < 2.87) return true;
               if (pc3z > 96 && pc3z < 98) return true;
               if (pc3z > 22 && pc3z < 128 && pc3phi < 2.6) return true;
               if (pc3z > 82 && pc3phi > 2.62
                  && pc3z < 108 && pc3phi < 2.675) return true;
               if (pc3z > 90 && pc3phi > 2.705
                  && pc3z < 108 && pc3phi < 2.75) return true;
               if (pc3z < 97)
               {
                  if (pc3z > 24 && pc3phi > 2.89
                     && pc3z < 32 && pc3phi < 2.91) return true;
                  if (pc3z > 34 && pc3phi > 2.91
                     && pc3z < 44 && pc3phi < 2.93) return true;
                  if (pc3z > 54 && pc3phi > 2.7
                     && pc3z < 68 && pc3phi < 2.74) return true;
                  if (pc3z > 72 && pc3phi > 2.77
                     && pc3z < 78 && pc3phi < 2.785) return true;
               }
               else
               {
                  if (pc3z > 118 && pc3phi > 2.65
                     && pc3z < 122 && pc3phi < 2.66) return true;
                  if (pc3z > 126 && pc3phi > 2.76
                     && pc3z < 132 && pc3phi < 2.775) return true;
                  if (pc3z > 142 && pc3phi > 2.81
                     && pc3z < 164 && pc3phi < 2.865) return true;
                  if (pc3z > 146 && pc3phi > 2.77
                     && pc3z < 168 && pc3phi < 2.815) return true;
               }
            }
            else
            {
               if (pc3z > 96 && pc3z < 98) return true;
               if (pc3z < 97)
               {
                  if (pc3z < 42) return true;
                  if (pc3phi > 2.495 && pc3z < 46) return true;
                  if (pc3phi > 2.42 && pc3z < 12 && pc3phi < 2.46) return true;
                  if (pc3phi < 2.245) return true;
                  if (pc3phi > 2.42
                     && pc3z < 66 && pc3phi < 2.47) return true;
                  if (pc3z > 40 && pc3phi > 2.305
                     && pc3z < 68 && pc3phi < 2.355) return true;
                  if (pc3z > 64 && pc3phi > 2.495 && pc3z < 78) return true;
                  if (pc3z > 82 && pc3phi > 2.38 && pc3phi < 2.425) return true;
               }
               else
               {
                  if (pc3z > 108 && pc3phi > 2.31
                     && pc3z < 128 && pc3phi < 2.35) return true;
                  if (pc3z > 142 && pc3phi > 2.35
                     && pc3z < 152 && pc3phi < 2.385) return true;
                  if (pc3z > -42 && pc3phi > 2.385
                     && pc3z < -32 && pc3phi < 2.425) return true;
                  if (pc3z > -80 && pc3phi > 2.265
                     && pc3z < -64 && pc3phi < 2.29) return true;
               }
            }
         }
      }
      else
      {
         if (pc3z < -168 || pc3phi < -0.565
            || pc3z > 164 || pc3phi > 0.97) return true;
         if (pc3z > -8 && pc3z < 10) return true;
         if (pc3z > 96 && pc3z < 100) return true;
         if (pc3z > -98 && pc3z < -94) return true;
         if (pc3phi > 0.58 && pc3phi < 0.615) return true;
         if (pc3phi > 0.08 && pc3phi < 0.23) return true;
         if (pc3phi > -0.205 && pc3phi < -0.17) return true;
         if (pc3z < 0)
         {
            if (pc3phi > 0.6)
            {
               if (pc3z < -96)
               {
                  if (pc3z > -126 && pc3phi > 0.86
                     && pc3z < -102 && pc3phi < 0.915) return true;
                  if (pc3z > -122 && pc3phi > 0.84
                     && pc3z < -114 && pc3phi < 0.865) return true;
                  if (pc3z > -126 && pc3phi > 0.67
                     && pc3z < -102 && pc3phi < 0.72) return true;
               }
               else
               {
                  if (pc3phi > 0.64 && pc3z < -86 && pc3phi < 0.675) return true;
                  if (pc3z > -84 && pc3z < -64 && pc3phi < 0.64) return true;
                  if (pc3z > -84 && pc3phi > 0.86
                     && pc3z < -70 && pc3phi < 0.955) return true;
                  if (pc3z > -44 && pc3phi > 0.9
                     && pc3z < -18 && pc3phi < 0.95) return true;
                  if (pc3z > -40 && pc3phi > 0.63
                     && pc3z < -22 && pc3phi < 0.675) return true;
               }
            }
            else if (pc3phi > 0.2)
            {
               if (pc3z > -106 && pc3phi > 0.395
                  && pc3z < -80 && pc3phi < 0.44) return true;
               if (pc3z > -126 && pc3phi > 0.47
                  && pc3z < -84 && pc3phi < 0.52) return true;
               if (pc3z < -96)
               {
                  if (pc3phi > 0.47 && pc3z < -160 && pc3phi < 0.51) return true;
                  if (pc3z > -136 && pc3z < -122 && pc3phi < 0.255) return true;
                  if (pc3z > -104 && pc3phi > 0.325 && pc3phi < 0.34) return true;
                  if (pc3z > -142 && pc3phi > 0.4
                     && pc3z < -130 && pc3phi < 0.43) return true;
                  if (pc3z > -114 && pc3phi > 0.38
                     && pc3z < -110 && pc3phi < 0.39) return true;
               }
               else
               {
                  if (pc3z > -74 && pc3phi > 0.545 && pc3z < -60) return true;
                  if (pc3z > -46 && pc3phi > 0.51 && pc3z < -30) return true;
                  if (pc3z > -66 && pc3phi > 0.35 && pc3phi < 0.445) return true;
                  if (pc3z > -20 && pc3phi < 0.265) return true;
                  if (pc3z > -76 && pc3phi > 0.37
                     && pc3z < -68 && pc3phi < 0.39) return true;
                  if (pc3z > -72 && pc3phi > 0.32
                     && pc3z < -58 && pc3phi < 0.37) return true;
                  if (pc3z > -64 && pc3phi > 0.245
                     && pc3z < -42 && pc3phi < 0.295) return true;
                  if (pc3z > -66 && pc3phi > 0.515
                     && pc3z < -60 && pc3phi < 0.53) return true;
               }
            }
            else if (pc3phi > -0.2)
            {
               if (pc3z < -96)
               {
                  if (pc3z > -122 && pc3z < -116 && pc3phi < -0.16) return true;
                  if (pc3z > -108 && pc3phi > -0.105
                     && pc3z < -104 && pc3phi < -0.095) return true;
               }
               else
               {
                  if (pc3z > -88 && pc3phi > 0.005 && pc3z < -74) return true;
                  if (pc3z > -78 && pc3phi > 0.04 && pc3z < -62) return true;
                  if (pc3z > -16 && pc3phi > -0.095 && pc3phi < -0.06) return true;
                  if (pc3z > -18 && pc3z < -12 && pc3phi < -0.14) return true;
                  if (pc3z > -84 && pc3phi > -0.075
                     && pc3z < -62 && pc3phi < -0.03) return true;
               }
            }
            else
            {
               if (pc3phi > -0.495 && pc3phi < -0.485) return true;
               if (pc3z < -96)
               {
                  if (pc3phi > -0.29 && pc3z < -144 && pc3phi < -0.24) return true;
                  if (pc3phi > -0.535 && pc3z < -142 && pc3phi < -0.485) return true;
                  if (pc3phi > -0.34 && pc3z < -164 && pc3phi < -0.31) return true;
                  if (pc3phi > -0.49 && pc3z < -158 && pc3phi < -0.38) return true;
                  if (pc3z > -102 && pc3phi > -0.395 && pc3phi < -0.38) return true;
                  if (pc3z > -158 && pc3phi > -0.31
                     && pc3z < -144 && pc3phi < -0.265) return true;
                  if (pc3z > -144 && pc3phi > -0.35
                     && pc3z < -126 && pc3phi < -0.305) return true;
                  if (pc3z > -120 && pc3phi > -0.32
                     && pc3z < -112 && pc3phi < -0.305) return true;
                  if (pc3z > -144 && pc3phi > -0.42
                     && pc3z < -130 && pc3phi < -0.38) return true;
                  if (pc3z > -112 && pc3phi > -0.425
                     && pc3z < -106 && pc3phi < -0.41) return true;
                  if (pc3z > -148 && pc3phi > -0.495
                     && pc3z < -130 && pc3phi < -0.45) return true;
               }
               else
               {
                  if (pc3phi > -0.315 && pc3z < -84 && pc3phi < -0.28) return true;
                  if (pc3z > -56 && pc3z < -42 && pc3phi < -0.525) return true;
                  if (pc3z > -86 && pc3phi > -0.49
                     && pc3z < -72 && pc3phi < -0.45) return true;
                  if (pc3z > -52 && pc3phi > -0.315
                     && pc3z < -42 && pc3phi < -0.275) return true;
               }
            }
         }
         else
         {
            if (pc3phi > 0.6)
            {
               if (pc3z < 96)
               {
                  if (pc3z > 14 && pc3phi > 0.805
                     && pc3z < 22 && pc3phi < 0.83) return true;
                  if (pc3z > 36 && pc3phi > 0.835
                     && pc3z < 50 && pc3phi < 0.865) return true;
                  if (pc3z > 74 && pc3phi > 0.83
                     && pc3z < 78 && pc3phi < 0.84) return true;
                  if (pc3z > 32 && pc3phi > 0.685
                     && pc3z < 44 && pc3phi < 0.715) return true;
                  if (pc3z > 86 && pc3phi > 0.755 && pc3phi < 0.795) return true;
                  if (pc3z > 76 && pc3phi > 0.945
                     && pc3z < 90 && pc3phi < 0.985) return true;
               }
               else
               {
                  if (pc3z > 106 && pc3phi > 0.94 && pc3z < 124) return true;
                  if (pc3z > 130 && pc3phi > 0.905
                     && pc3z < 140 && pc3phi < 0.95) return true;
                  if (pc3z > 104 && pc3phi > 0.89
                     && pc3z < 108 && pc3phi < 0.9) return true;
                  if (pc3z > 104 && pc3phi > 0.79
                     && pc3z < 118 && pc3phi < 0.835) return true;
                  if (pc3z > 136 && pc3phi > 0.855
                     && pc3z < 142 && pc3phi < 0.865) return true;
                  if (pc3z > 136 && pc3phi > 0.825
                     && pc3z < 140 && pc3phi < 0.835) return true;
                  if (pc3z > 130 && pc3phi > 0.755
                     && pc3z < 146 && pc3phi < 0.795) return true;
                  if (pc3z > 104 && pc3phi > 0.675
                     && pc3z < 110 && pc3phi < 0.69) return true;
               }
               
            }
            else if (pc3phi > 0.2)
            {
               if (pc3phi > 0.49 && pc3phi < 0.5) return true;
               if (pc3z < 98)
               {
                  if (pc3z > 14 && pc3phi > 0.48
                     && pc3z < 28 && pc3phi < 0.52) return true;
                  if (pc3z > 62 && pc3phi > 0.47
                     && pc3z < 88 && pc3phi < 0.52) return true;
                  if (pc3z > 26 && pc3phi > 0.445
                     && pc3z < 30 && pc3phi < 0.455) return true;
                  if (pc3z > 54 && pc3phi > 0.34
                     && pc3z < 60 && pc3phi < 0.35) return true;
                  if (pc3z > 72 && pc3z < 82 && pc3phi < 0.25) return true;
                  if (pc3z > 86 && pc3phi > 0.55 && pc3phi < 0.57) return true;
               }
               else
               {
                  if (pc3z > 98 && pc3phi > 0.475
                     && pc3z < 108 && pc3phi < 0.52) return true;
                  if (pc3z > 134 && pc3phi > 0.285
                     && pc3z < 150 && pc3phi < 0.335) return true;
                  if (pc3z > 148 && pc3phi > 0.34
                     && pc3z < 154 && pc3phi < 0.36) return true;
               }
            }
            else if (pc3phi > -0.2)
            {
               if (pc3z < 98)
               {
                  if (pc3phi > 0.045 && pc3z < 24) return true;
                  if (pc3phi > -0.07 && pc3z < 18 && pc3phi < -0.02) return true;
                  if (pc3phi > -0.125 && pc3z < 18 && pc3phi < -0.1) return true;
                  if (pc3z > 46 && pc3z < 60 && pc3phi < -0.14) return true;
                  if (pc3z > 66 && pc3phi > 0.035 && pc3z < 80) return true;
                  if (pc3z > 92 && pc3phi > -0.15
                     && pc3z < 98 && pc3phi < -0.14) return true;
               }
               else
               {
                  if (pc3z > 154 && pc3phi > -0.115 && pc3phi < -0.065) return true;
                  if (pc3z > 104 && pc3phi > -0.02
                     && pc3z < 114 && pc3phi < 0.01) return true;
                  if (pc3z > 144 && pc3phi > 0.01
                     && pc3z < 156 && pc3phi < 0.05) return true;
                  if (pc3z > 114 && pc3phi > -0.11
                     && pc3z < 126 && pc3phi < -0.075) return true;
               }
            }
            else
            {
               if (pc3z < 98)
               {
                  if (pc3z > 24 && pc3phi > -0.235 && pc3z < 46) return true;
                  if (pc3phi > -0.345 && pc3z < 26 && pc3phi < -0.3) return true;
                  if (pc3z > 46 && pc3z < 56 && pc3phi < -0.53) return true;
                  if (pc3z > 82 && pc3phi > -0.35 && pc3phi < -0.265) return true;
                  if (pc3z > 60 && pc3phi > -0.315
                     && pc3z < 74 && pc3phi < -0.275) return true;
                  if (pc3z > 52 && pc3phi > -0.505
                     && pc3z < 62 && pc3phi < -0.455) return true;
                  if (pc3z > 44 && pc3phi > -0.275
                     && pc3z < 50 && pc3phi < -0.255) return true;
                  if (pc3z > 52 && pc3phi > -0.385
                     && pc3z < 56 && pc3phi < -0.375) return true;
                  if (pc3z > 72 && pc3phi > -0.35
                     && pc3z < 78 && pc3phi < -0.335) return true;
               }
               else
               {
                  if (pc3z > 142 && pc3phi < -0.525) return true;
                  if (pc3phi > -0.35 && pc3z < 108 && pc3phi < -0.3) return true;
                  if (pc3z > 106 && pc3phi > -0.27
                     && pc3z < 114 && pc3phi < -0.25) return true;
                  if (pc3z > 136 && pc3phi > -0.375
                     && pc3z < 144 && pc3phi < -0.36) return true;
                  if (pc3z > 134 && pc3phi > -0.465
                     && pc3z < 140 && pc3phi < -0.45) return true;
                  if (pc3z > 146 && pc3phi > -0.3
                     && pc3z < 154 && pc3phi < -0.275) return true;
               }
            }
         }
      }
      return false;
   }

   bool IsDeadEMCal(const double phi, const double zed, const int sect, 
                    const double pemcy, const double pemcz) // Run7AuAu200
   {
      if (phi >= 1.5 && zed >= 0)
      {
         switch (sect)
         {
            case 0:
               if (pemcy < -274.8 || pemcz < 6.57143
                  || pemcy > -113.1 || pemcz > 190.571) return true;
               if (pemcy > -150.9 && pemcy < -127.8) return true;
               if (pemcy > -176.1 && pemcy < -148.8 && pemcz < 51.1429) return true;
               if (pemcy > -243.3 && pemcz > 80.8571
                  && pemcy < -226.5 && pemcz < 86.5714) return true;
               if (pemcy > -226.5 && pemcz > 140.286
                  && pemcy < -218.1 && pemcz < 148.286) return true;
               if (pemcy > -224.4 && pemcz > 42
                  && pemcy < -213.9 && pemcz < 47.7143) return true;
               break;
            case 1:
               if (pemcy < -70.4 || pemcz < 6.57143
                  || pemcy > 46.2 || pemcz > 181.429) return true;
               if (pemcz > 34 && pemcy < -41.8 && pemcz < 42) return true;
               if (pemcy > 8.8 && pemcz > 88.8571
                  && pemcy < 15.4 && pemcz < 95.7143) return true;
               if (pemcy > 33 && pemcz > 64.8571
                  && pemcy < 37.4 && pemcz < 70.5714) return true;
               break;
            case 2:
               if (pemcy < 108.9 || pemcz < 3.14286
                  || pemcy > 287.4 || pemcz > 176.857) return true;
               if (pemcy > 274.8 && pemcz > 31.7143 && pemcz < 38.5714) return true;
               if (pemcy > 237 && pemcz > 130 && pemcy < 251.7) return true;
               if (pemcy > 247.5 && pemcz > 96.8571
                  && pemcy < 260.1 && pemcz < 132.286) return true;
               if (pemcy > 146.7 && pemcz > 59.1429
                  && pemcy < 155.1 && pemcz < 67.1429) return true;
               if (pemcy > 138.3 && pemcz > 11.1429
                  && pemcy < 144.6 && pemcz < 18) return true;
               if (pemcy > 243.3 && pemcz > 53.4286
                  && pemcy < 253.8 && pemcz < 60.2857) return true;
               if (pemcy > 232.8 && pemcz > 31.7143
                  && pemcy < 245.4 && pemcz < 38.5714) return true;
               if (pemcy > 262.2 && pemcz > 27.1429
                  && pemcy < 272.7 && pemcz < 39.7143) return true;
               break;
            case 3:
               if (pemcy < 292.8 || pemcz < 11.1429
                  || pemcy > 428.8 || pemcz > 167.714) return true;
               if (pemcy < 300.8 && pemcz < 67.1429) return true;
               if (pemcy > 412.8 && pemcz > 66 && pemcz < 78.5714) return true;
               if (pemcy > 412.8 && pemcz > 44.2857 && pemcz < 56.8571) return true;
               if (pemcy > 321.6 && pemcz > 66
                  && pemcy < 331.2 && pemcz < 131.143) return true;
               if (pemcy > 328 && pemcz > 66
                  && pemcy < 337.6 && pemcz < 77.4286) return true;
               if (pemcy > 307.2 && pemcz > 44.2857
                  && pemcy < 312 && pemcz < 51.1429) return true;
               if (pemcy > 372.8 && pemcz > 152.857
                  && pemcy < 377.6 && pemcz < 160.857) return true;
               break;
         }
      }
      else if (phi >= 1.5 && zed < 0)
      {
         switch (sect)
         {
            case 0:
               if (pemcy < -289.5 || pemcz < -189.429
                  || pemcy > -113.1 || pemcz > 1.42857) return true;
               if (pemcy > -150.9 && pemcy < -127.8) return true;
               if (pemcy > -197.1 && pemcz > -60.2857
                  && pemcy < -184.5 && pemcz < -53.4286) return true;
               if (pemcy > -209.7 && pemcz > -54.5714
                  && pemcy < -199.2 && pemcz < -44.2857) return true;
               if (pemcy > -205.5 && pemcz > -95.7143
                  && pemcy < -199.2 && pemcz < -88.8571) return true;
               if (pemcy > -174 && pemcz > -83.1429
                  && pemcy < -163.5 && pemcz < -76.2857) return true;
               if (pemcy > -241.2 && pemcz > -63.7143
                  && pemcy < -224.4 && pemcz < -56.8571) return true;
               if (pemcy > -237 && pemcz > -140.286
                  && pemcy < -226.5 && pemcz < -133.429) return true;
               if (pemcy > -264.3 && pemcz > -60.2857
                  && pemcy < -253.8 && pemcz < -53.4286) return true;
               if (pemcy > -283.2 && pemcz > -94.5714
                  && pemcy < -276.9 && pemcz < -87.7143) return true;
               break;
            case 1:
               if (pemcy < -92.4 || pemcz < -189.429
                  || pemcy > 101.2 || pemcz > -3.14286) return true;
               if (pemcy > 74.8 && pemcy < 92.4) return true;
               if (pemcy > 2.2 && pemcy < 19.8) return true;
               if (pemcy > -55 && pemcz > -43.1429
                  && pemcy < -44 && pemcz < -36.2857) return true;
               if (pemcy > -35.2 && pemcz > -62.5714
                  && pemcy < -28.6 && pemcz < -58) return true;
               if (pemcy > -77 && pemcz > -117.429
                  && pemcy < -59.4 && pemcz < -109.429) return true;
               if (pemcy > -77 && pemcz > -109.429
                  && pemcy < -63.8 && pemcz < -95.7143) return true;
               if (pemcy > -55 && pemcz > -51.1429
                  && pemcy < -46.2 && pemcz < -46.5714) return true;
               break;
            case 2:
               if (pemcy < 108.9 || pemcz < -159.714
                  || pemcy > 287.4 || pemcz > 0.285714) return true;
               if (pemcy > 218.1 && pemcy < 234.9) return true;
               if (pemcy > 123.6 && pemcz > -59.1429
                  && pemcy < 129.9 && pemcz < -53.4286) return true;
               if (pemcy > 127.8 && pemcz > -136.857
                  && pemcy < 134.1 && pemcz < -128.857) return true;
               if (pemcy > 161.4 && pemcz > -163.143
                  && pemcy < 167.7 && pemcz < -155.143) return true;
               if (pemcy > 211.8 && pemcz > -22.5714
                  && pemcy < 220.2 && pemcz < -15.7143) return true;
               if (pemcy > 234.9 && pemcz > -92.2857
                  && pemcy < 243.3 && pemcz < -85.4286) return true;
               if (pemcy > 262.2 && pemcz > -103.714
                  && pemcy < 268.5 && pemcz < -98) return true;
               if (pemcy > 171.9 && pemcz > -132.286
                  && pemcy < 182.4 && pemcz < -119.714) return true;
               break;
            case 3:
               if (pemcy < 292.8 || pemcz < -182.571
                  || pemcy > 427.2 || pemcz > -2) return true;
               if (pemcz > -59.1429 && pemcy < 299.2 && pemcz < -52.2857) return true;
               if (pemcy > 337.6 && pemcz > -163.143
                  && pemcy < 342.4 && pemcz < -152.857) return true;
               if (pemcy > 377.6 && pemcz > -21.4286
                  && pemcy < 382.4 && pemcz < -14.5714) return true;
               break;
         }
      }
      else if (phi < 1.5 && zed >= 0)
      {
         switch (sect)
         {
            case 0:
               if (pemcy < -228.6 || pemcz < 5.42857
                  || pemcy > -102.6 || pemcz > 182.571) return true;
               if (pemcy > -186.6 && pemcz > 70.5714
                  && pemcy < -178.2 && pemcz < 77.4286) return true;
               if (pemcy > -117.3 && pemcz > 44.2857
                  && pemcy < -108.9 && pemcz < 50) return true;
               if (pemcy > -121.5 && pemcz > 38.5714
                  && pemcy < -115.2 && pemcz < 45.4286) return true;
               break;
            case 1:
               if (pemcy < -92.4 || pemcz < 4.28571
                  || pemcy > 52.8 || pemcz > 192.857) return true;
               if (pemcy > -74.8 && pemcz > 102.571
                  && pemcy < -63.8 && pemcz < 110.571) return true;
               break;
            case 2:
               if (pemcy < 113.1 || pemcz < 5.42857
                  || pemcy > 287.4 || pemcz > 187.143) return true;
               if (pemcy > 180.3 && pemcz > 27.1429
                  && pemcy < 188.7 && pemcz < 35.1429) return true;
               break;
            case 3:
               if (pemcy < 292.8 || pemcz < 8.85714
                  || pemcy > 425.6 || pemcz > 180.286) return true;
               if (pemcy > 406.4 && pemcz > 134.571 && pemcz < 158.571) return true;
               if (pemcy > 409.6 && pemcz > 114
                  && pemcy < 425.6 && pemcz < 135.714) return true;
               if (pemcy > 344 && pemcz > 151.714
                  && pemcy < 353.6 && pemcz < 175.714) return true;
               if (pemcy > 352 && pemcz > 21.4286
                  && pemcy < 363.2 && pemcz < 55.7143) return true;
               if (pemcy > 371.2 && pemcz > 86.5714
                  && pemcy < 377.6 && pemcz < 94.5714) return true;
               break;
         }
      }
      else
      {
         switch (sect)
         {
            case 0:
               if (pemcy < -274.8 || pemcz < -187.143
                  || pemcy > -104.7 || pemcz > 0.285714) return true;
               if (pemcy > -115.2 && pemcz > -76.2857 && pemcz < -63.7143) return true;
               if (pemcy > -132 && pemcz > -151.714
                  && pemcy < -119.4 && pemcz < -128.857) return true;
               if (pemcy > -121.5 && pemcz > -108.286
                  && pemcy < -115.2 && pemcz < -101.429) return true;
               break;
            case 1:
               if (pemcy < -92.4 || pemcz < -182.571
                  || pemcy > 37.4 || pemcz > -0.857143) return true;
               if (pemcy > -19.8 && pemcz > -174.571
                  && pemcy < -6.6 && pemcz < -150.571) return true;
               break;
            case 2:
               if (pemcy < 127.8 || pemcz < -183.714
                  || pemcy > 285.3 || pemcz > -2) return true;
               if (pemcy > 195 && pemcy < 203.4 && pemcz < -175.714) return true;
               if (pemcy > 276.9 && pemcz > -140.286 && pemcz < -131.143) return true;
               if (pemcy > 266.4 && pemcz > -167.714 && pemcz < -159.714) return true;
               if (pemcy > 274.8 && pemcz < -144.857) return true;
               if (pemcy > 268.5 && pemcz < -174.571) return true;
               if (pemcy > 136.2 && pemcz > -108.286
                  && pemcy < 148.8 && pemcz < -94.5714) return true;
               if (pemcy > 138.3 && pemcz > -131.143
                  && pemcy < 148.8 && pemcz < -117.429) return true;
               if (pemcy > 171.9 && pemcz > -179.143
                  && pemcy < 178.2 && pemcz < -170) return true;
               if (pemcy > 205.5 && pemcz > -10
                  && pemcy < 211.8 && pemcz < -4.28571) return true;
               if (pemcy > 266.4 && pemcz > -141.429
                  && pemcy < 279 && pemcz < -126.571) return true;
               break;
            case 3:
               if (pemcy < 292.8 || pemcz < -174.571
                  || pemcy > 427.2 || pemcz > -0.857143) return true;
               if (pemcz > -80.8571 && pemcy < 297.6 && pemcz < -74) return true;
               if (pemcz > -96.8571 && pemcy < 297.6 && pemcz < -87.7143) return true;
               if (pemcz > -131.143 && pemcy < 296 && pemcz < -120.857) return true;
               if (pemcy > 404.8 && pemcz > -12.2857 && pemcy < 411.2) return true;
               if (pemcy > 347.2 && pemcz < -140.286) return true;
               if (pemcy > 409.6 && pemcy < 424 && pemcz < -127.714) return true;
               if (pemcy > 313.6 && pemcz > -88.8571
                  && pemcy < 321.6 && pemcz < -78.5714) return true;
               if (pemcy > 316.8 && pemcz > -54.5714
                  && pemcy < 324.8 && pemcz < -40.8571) return true;
               if (pemcy > 329.6 && pemcz > -98
                  && pemcy < 339.2 && pemcz < -74) return true;
               if (pemcy > 336 && pemcz > -166.571
                  && pemcy < 385.6 && pemcz < -133.429) return true;
               if (pemcy > 380.8 && pemcz > -131.143
                  && pemcy < 420.8 && pemcz < -115.143) return true;
               if (pemcy > 401.6 && pemcz > -59.1429
                  && pemcy < 414.4 && pemcz < -51.1429) return true;
               if (pemcy > 344 && pemcz > -134.571
                  && pemcy < 371.2 && pemcz < -120.857) return true;
               break;
         }
      }
      
      return false;
   }

   bool IsBadSlat(const int slat) // Run7AuAu200
   {
      switch (slat)
      { 
         //fallthrough for bad slat values
         case 0: 
         case 1: 
         case 2: 
         case 3: 
         case 4: 
         case 5: 
         case 8: 
         case 16: 
         case 17: 
         case 18: 
         case 19: 
         case 20: 
         case 32: 
         case 33: 
         case 34: 
         case 35: 
         case 36: 
         case 48: 
         case 49: 
         case 50: 
         case 51: 
         case 52: 
         case 64: 
         case 65: 
         case 66: 
         case 67: 
         case 68: 
         case 69: 
         case 70: 
         case 71: 
         case 80: 
         case 81: 
         case 82: 
         case 83: 
         case 84: 
         case 85: 
         case 86: 
         case 87: 
         case 88: 
         case 89: 
         case 90: 
         case 91: 
         case 125: 
         case 142: 
         case 151: 
         case 179: 
         case 180: 
         case 214: 
         case 248: 
         case 251: 
         case 253: 
         case 255: 
         case 292: 
         case 304: 
         case 305: 
         case 306: 
         case 309: 
         case 310: 
         case 311: 
         case 324: 
         case 335: 
         case 339: 
         case 354: 
         case 362: 
         case 367: 
         case 376: 
         case 377: 
         case 378: 
         case 379: 
         case 383: 
         case 401: 
         case 452: 
         case 464: 
         case 465: 
         case 469: 
         case 501: 
         case 539: 
         case 552: 
         case 553: 
         case 554: 
         case 555: 
         case 556: 
         case 595: 
         case 602: 
         case 628: 
         case 681: 
         case 682: 
         case 683: 
         case 684: 
         case 685: 
         case 686: 
         case 687: 
         case 699: 
         case 700: 
         case 701: 
         case 702: 
         case 703: 
         case 714: 
         case 715: 
         case 716: 
         case 717: 
         case 718: 
         case 719: 
         case 728: 
         case 729: 
         case 730: 
         case 731: 
         case 732: 
         case 733: 
         case 734: 
         case 735: 
         case 739: 
         case 742: 
         case 745: 
         case 746: 
         case 747: 
         case 748: 
         case 749: 
         case 750: 
         case 751: 
         case 757: 
         case 758: 
         case 759: 
         case 760: 
         case 761: 
         case 762: 
         case 763: 
         case 764: 
         case 765: 
         case 766: 
         case 767: 
         case 783: 
         case 785: 
         case 788: 
         case 791: 
         case 800: 
         case 802: 
         case 803: 
         case 820: 
         case 835: 
         case 848: 
         case 863: 
         case 874: 
         case 904: 
         case 905: 
         case 906: 
         case 907: 
         case 909: 
         case 910: 
         case 911: 
         case 919: 
         case 930: 
         case 958: 
            return true; 
            break; 
      } 
   return false; 
   }

   bool IsBadStripTOFw(const int strip) // Run7AuAu200
   {
      switch (strip) 
      { 
         //fallthrough for bad strip values
         case 64: 
         case 65: 
         case 66: 
         case 67: 
         case 68: 
         case 69: 
         case 70: 
         case 72: 
         case 184: 
         case 187: 
         case 188: 
         case 189: 
         case 190: 
         case 191: 
         case 248: 
         case 249: 
         case 250: 
         case 251: 
         case 252: 
         case 253: 
         case 254: 
         case 255: 
         case 256: 
         case 257: 
         case 258: 
         case 320: 
         case 321: 
         case 322: 
         case 323: 
         case 324: 
         case 440: 
         case 441: 
         case 442: 
         case 443: 
         case 444: 
         case 445: 
         case 446: 
         case 447: 
         case 504: 
         case 507: 
         case 508: 
         case 509: 
         case 510: 
         case 511: 
            return true; 
            break; 
      } 
      return false;
   }

   bool IsDeadTOFe(const double zed, const double tofy, const double tofz) // Run7AuAu200
   {
      if (tofy < -270 || tofy > 82) return true;
      if (zed >= 0)
      {
         if (tofz < 3.75 || tofz > 172.5) return true;
         if (tofy > -120 && tofy < -78) return true;
         if (tofz > 50 && tofz < 52.5) return true;
         
         if (tofy < -100)
         {
            if (tofy > -212 && tofy < -210) return true;
            if (tofy > -232 && tofy < -230) return true;
            if (tofy > -150 && tofy < -148) return true;
            if (tofy > -150 && tofz > 45) return true;
            if (tofy > -208 && tofz > 47.5 && tofy < -168) return true;
            if (tofz > 35 && tofy < -232 && tofz < 40) return true;
            if (tofy > -168 && tofz > 10 && tofz < 15) return true;
            if (tofy > -172 && tofy < -168 && tofz < 27.5) return true;
            if (tofy > -212 && tofz > 27.5
               && tofy < -168 && tofz < 30) return true;
            if (tofy > -212 && tofz > 30
               && tofy < -208 && tofz < 51.25) return true;
         }
         else
         {
            if (tofy > -20 && tofy < -16) return true;
            if (tofy > 24 && tofy < 28) return true;
            if (tofy > 48 && tofy < 50) return true;
            if (tofy > -42 && tofy < -40) return true;
            if (tofz > 148.75 && tofz < 150) return true;
            if (tofz > 100 && tofz < 101.25) return true;
            if (tofz > 66.25 && tofy < -20 && tofz < 72.5) return true;
            if (tofz > 110 && tofy < -20 && tofz < 115) return true;
            if (tofz > 131.25 && tofy < -20 && tofz < 136.25) return true;
            if (tofy > -82 && tofz > 2.5
               && tofy < -16 && tofz < 11.25) return true;
            if (tofy > 44 && tofz < 10) return true;
            if (tofy > 48 && tofz > 162.5) return true;
            if (tofy > -18 && tofz > 112.5
               && tofy < 48 && tofz < 118.75) return true;
            if (tofy > -18 && tofz > 83.75
               && tofy < 48 && tofz < 91.25) return true;
            if (tofy > 48 && tofz > 83.75 && tofz < 100) return true;
         }
      }
      else
      {
         if (tofy < -268 || tofz < -173.75
            || tofy > 82 || tofz > 0) return true;
         if (tofy > -122 && tofy < -78) return true;
         if (tofz > -48.75 && tofz < -46.25) return true;
         if (tofy < -100)
         {
            if (tofy > -234 && tofy < -230) return true;
            if (tofy > -212 && tofy < -210) return true;
            if (tofy > -150 && tofy < -146) return true;
            if (tofy > -170 && tofz > -31.25 && tofy < -168) return true;
            if (tofz > -36.25 && tofy < -150 && tofz < -30) return true;
            if (tofy > -230 && tofy < -208 && tofz < -20) return true;
            if (tofy > -212 && tofy < -168 && tofz < -41.25) return true;
            if (tofy > -146 && tofz > -36.25 && tofz < -32.5) return true;
         }
         else
         {
            if (tofy > -44 && tofy < -40) return true;
            if (tofy > -20 && tofy < -18) return true;
            if (tofy > 24 && tofy < 28) return true;
            if (tofy > 48 && tofy < 50) return true;
            if (tofz > -147.5 && tofz < -145) return true;
            if (tofz > -97.5 && tofz < -96.25) return true;
            if (tofy > 48 && tofz < -167.5) return true;
            if (tofy > 48 && tofz > -23.75) return true;
            if (tofz > -33.75 && tofy < -18 && tofz < -25) return true;
            if (tofz > -78.75 && tofy < -18 && tofz < -73.75) return true;
            if (tofz > -107.5 && tofy < -18 && tofz < -101.25) return true;
            if (tofz > -93.75 && tofy < -42 && tofz < -90) return true;
            if (tofz > -135 && tofy < -20 && tofz < -132.5) return true;
            if (tofy > 28 && tofz > -41.25 && tofz < -35) return true;
            if (tofy > 28 && tofz > -17.5 && tofz < -11.25) return true;
            if (tofy > 50 && tofz > -138.75 && tofz < -132.5) return true;
            if (tofy > 28 && tofz > -97.5 && tofz < -93.75) return true;
            if (tofy > -20 && tofz > -38.75
               && tofy < 24 && tofz < -28.75) return true;
            if (tofy > -18 && tofz > -63.75
               && tofy < 50 && tofz < -47.5) return true;
            if (tofy > -18 && tofz > -73.75
               && tofy < 50 && tofz < -67.5) return true;
            if (tofy > -42 && tofz > -102.5
               && tofy < 24 && tofz < -96.25) return true;
            if (tofy > -18 && tofz > -126.25
               && tofy < 48 && tofz < -121.25) return true;
         }
      }
      return false;
   }

   bool IsDeadTOFw(const double zed, const double board, const double alpha) // Run7AuAu200
   {
      if (board > 60.845 || board < 19.505) return true;
      if (alpha > (board-24.)*0.033 || alpha < (board-58.3)*0.034) return true;
      if (alpha > (board-28.2)*0.034 && alpha < (board-27.8)*0.034) return true;
      if (alpha > (board-50.3)*0.034 && alpha < (board-32.2)*0.034) return true;
      if (alpha > (board-54.5)*0.034 && alpha < (board-54.2)*0.034) return true;
      if (alpha > (board-29.6)*0.17 && alpha < (board-29.3)*0.17) return true;
      if (alpha > (board-49.3)*0.17 && alpha < (board-49.)*0.17) return true;
      if (alpha > (board-59.4)*0.17 && alpha < (board-59.1)*0.17) return true;
      if (zed >= 0)
      {
         if (board > 56.87 && board < 57.93) return true;
         if (board > 55.015 && board < 56.075) return true;
      }
      else
      {
         if (board < 24.54) return true;
         if (board > 30.635 && board < 33.285) return true;
         if (alpha > (board-54)*0.4 && alpha < (board-53.5)*0.4) return true;
      }
      return false;
   }
}

namespace Run14HeAu200MBCuts
{
   /*
   short runGroup;
   
   void SetRunGroup(const int runNumber)
   {
      if (runNumber < 415835) runGroup = 0;
      else if (runNumber < 416442) runGroup = 3;
      else if (runNumber < 416379) runGroup = 2;
      else runGroup = 1;
   }
   */
   
   bool IsDeadDC(const double phi, const double zed, 
                 const double board, const double alpha) // Run7AuAu200
   {
      if (phi > M_PI/2.)
      {
         if (alpha > -0.102667 + board*0.266667) return true;
         if (alpha < 0.0957352 + board*-0.206272)  return true;
         if (alpha > 14.3202 + board*-0.181185) return true;
         if (alpha < -87.2597 + board*1.1122) return true;
         
         if (zed >= 0)
         {
            if (alpha > -6.6834 + board*0.310027 && alpha < -7.09342 + board*0.347317) return true;
            if (alpha < 12.4846 + board*-0.581463 && alpha > 7.2617 + board*-0.358537) return true;
            if (alpha > -7.77282 + board*0.303833 && alpha < -11.0491 + board*0.476098) return true;
            if (alpha < -10.1516 + board*0.323345 && alpha > -10.6781 + board*0.301045) return true;
            if (alpha < 24.9954 + board*-0.62439 && alpha > 16.329 + board*-0.426016) return true;
            if (alpha < 30.7301 + board*-0.483902 && alpha > 24.6596 + board*-0.41626) return true;
            if (alpha < 41.6035 + board*-0.636098 && alpha > 36.7232 + board*-0.565854) return true;
         }
         else
         {
            if (alpha < -9.19242 + board*0.393496 && alpha > -6.07446 + board*0.239024) return true;
            if (alpha > 13.9444 + board*-0.363415 && alpha < 32.132 + board*-0.8) return true;
            if (alpha > 23.077 + board*-0.423693 && alpha < 33.5927 + board*-0.604878) return true;
            if (alpha > 21.1754 + board*-0.357724 && alpha < 26.5033 + board*-0.417561) return true;
            if (alpha < -16.2884 + board*0.517073 && alpha > -19.709 + board*0.55935) return true;
            if (alpha < -15.0165 + board*0.546341 && alpha > -10.0039 + board*0.356098) return true;
            if (alpha > 3.66712 + board*-0.456585 && alpha < 6.08358 + board*-0.539837) return true;
            if (alpha > 0. && alpha > 15.4388 + board*-0.44878 && 
                alpha < 22.8834 + board*-0.62439) return true;
            if (alpha > 0. && alpha > 6.45502 + board*-0.107317 && 
                alpha < 6.06966 + board*-0.097561) return true;
         }
      }
      else
      {
         if (alpha > -0.71322 + board*0.165854) return true;
         if (alpha < 1.62663 + board*-0.390244) return true;
         if (alpha > 65.1375 + board*-0.839024) return true;
         if (alpha < -11.9492 + board*0.15122) return true;
         
         if (zed >= 0)
         {
            if (alpha < -4.67151 + board*0.520325 && alpha > -5.55649 + board*0.546341) return true;
            if (alpha > 6.31318 + board*-0.304878 && alpha < 6.15108 + board*-0.284878) return true;
            if (alpha < -19.712 + board*0.8 && alpha > -11.9702 + board*0.458537) return true;
            if (alpha < -15.0711 + board*0.456585 && alpha > -16.6123 + board*0.487805) return true;
            if (alpha < -25.7918 + board*0.741463 && alpha > -25.2103 + board*0.702439) return true;
            if (alpha < -26.2051 + board*0.536585 && alpha > -32.087 + board*0.630894) return true;
            if (alpha > 8.47415 + board*-0.243902 && alpha < 12.593 + board*-0.292683) return true;
            if (alpha < -24.8661 + board*0.643902 && alpha > -19.9127 + board*0.461789) return true;
            if (alpha > 9.26263 + board*-0.390244 && alpha < 16.5577 + board*-0.630894) return true;
            if (alpha > 15.3738 + board*-0.279675 && alpha < 32.2756 + board*-0.539837) return true;
            if (alpha < -19.1533 + board*0.342547 && alpha > -32.7472 + board*0.565854) return true;
            if (alpha < -27.4379 + board*0.442276 && alpha > -25.4605 + board*0.386341) return true;
            if (alpha < -28.9814 + board*0.409756 && alpha > -35.8149 + board*0.478049) return true;
            if (alpha > 0. && alpha > 20.345 + board*-0.292683 && 
                alpha < 19.1742 + board*-0.273171) return true;
         }
         else
         {
            if (alpha < -9.028 + board*0.533333 && alpha > -9.7874 + board*0.538537) return true;
            if (alpha > 4.94498 + board*-0.292683 && alpha < 6.04644 + board*-0.331707) return true;
            if (alpha > 6.25234 + board*-0.302439 && alpha < 7.60376 + board*-0.32892) return true;
            if (alpha > 6.53424 + board*-0.273171 && alpha < 8.39234 + board*-0.323035) return true;
            if (alpha < -10.347 + board*0.414634 && alpha > -10.9841 + board*0.421951) return true;
            if (alpha < -16.4177 + board*0.497561 && alpha > -22.3426 + board*0.656911) return true;
            if (alpha > 12.0949 + board*-0.344715 && alpha < 15.4352 + board*-0.358537) return true;
            if (alpha < -28.6937 + board*0.589268 && alpha > -32.5245 + board*0.639024) return true;
            if (alpha > 14.3998 + board*-0.261463 && alpha < 21.3267 + board*-0.347967) return true;
            if (alpha < -22.1929 + board*0.396748 && alpha > -20.6637 + board*0.323902) return true;
            if (alpha < -42.3668 + board*0.634146 && alpha > -32.7348 + board*0.437073) return true;
         }
      }
      return false;
   }

   bool IsDeadPC1(const double phi, const double pc1z, const double pc1phi) // Run7AuAu200
   {
      if (phi > M_PI/2.)
      {
         if (pc1z < -83.5 || pc1phi < 2.17 ||
             pc1z > 84. || pc1phi > 3.71) return true;
         
         const short pc1zBin = static_cast<short>((pc1z + 90.)*2.) + 1;
         const short pc1phiBin = static_cast<short>((pc1phi - 2.05)*100.) + 1;
         
         switch(pc1phiBin)
         {
            case 13: switch(pc1zBin) {case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 21: case 22: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 340: case 341: case 342: case 343: case 344: case 345: case 346: case 347: case 348: return true; break;} break;
            case 14: switch(pc1zBin) {case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 21: case 22: case 27: case 28: case 29: case 30: case 31: case 32: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 340: case 341: case 342: case 343: case 344: case 345: case 346: case 347: case 348: return true; break;} break;
            case 15: switch(pc1zBin) {case 27: case 28: case 29: case 30: case 31: case 32: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 16: switch(pc1zBin) {case 27: case 28: case 29: case 30: case 31: case 32: case 115: case 116: case 117: case 118: case 119: case 120: case 121: case 122: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 17: switch(pc1zBin) {case 27: case 28: case 29: case 30: case 31: case 32: case 115: case 116: case 117: case 118: case 119: case 120: case 121: case 122: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 18: switch(pc1zBin) {case 27: case 28: case 29: case 30: case 31: case 32: case 115: case 116: case 117: case 118: case 119: case 120: case 121: case 122: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 19: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 20: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 21: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 22: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 23: switch(pc1zBin) {case 90: case 91: case 92: case 93: case 94: case 95: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 24: switch(pc1zBin) {case 90: case 91: case 92: case 93: case 94: case 95: case 115: case 116: case 117: case 118: case 119: case 120: case 121: case 122: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 25: switch(pc1zBin) {case 115: case 116: case 117: case 118: case 119: case 120: case 121: case 122: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 26: switch(pc1zBin) {case 115: case 116: case 117: case 118: case 119: case 120: case 121: case 122: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 27: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 28: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 29: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 30: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 31: switch(pc1zBin) {case 149: case 150: case 151: case 152: case 153: case 154: case 155: case 156: case 157: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 32: switch(pc1zBin) {case 149: case 150: case 151: case 152: case 153: case 154: case 155: case 156: case 157: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 33: switch(pc1zBin) {case 149: case 150: case 151: case 152: case 153: case 154: case 155: case 156: case 157: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 34: switch(pc1zBin) {case 149: case 150: case 151: case 152: case 153: case 154: case 155: case 156: case 157: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 35: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 36: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 37: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 38: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 39: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 40: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 41: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 42: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 43: return true; break;
            case 44: return true; break;
            case 45: return true; break;
            case 46: return true; break;
            case 47: return true; break;
            case 48: return true; break;
            case 49: return true; break;
            case 50: return true; break;
            case 51: return true; break;
            case 52: return true; break;
            case 53: switch(pc1zBin) {case 170: case 171: case 172: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 54: switch(pc1zBin) {case 170: case 171: case 172: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 55: switch(pc1zBin) {case 170: case 171: case 172: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 56: switch(pc1zBin) {case 170: case 171: case 172: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 57: switch(pc1zBin) {case 170: case 171: case 172: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 58: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 59: switch(pc1zBin) {case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 21: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 60: switch(pc1zBin) {case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 21: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 61: switch(pc1zBin) {case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 21: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 62: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 197: case 198: case 199: case 200: case 201: case 202: case 203: case 204: return true; break;} break;
            case 63: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 197: case 198: case 199: case 200: case 201: case 202: case 203: case 204: case 211: case 212: case 213: case 214: case 215: case 216: case 217: return true; break;} break;
            case 64: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 197: case 198: case 199: case 200: case 201: case 202: case 203: case 204: case 211: case 212: case 213: case 214: case 215: case 216: case 217: return true; break;} break;
            case 65: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 197: case 198: case 199: case 200: case 201: case 202: case 203: case 204: return true; break;} break;
            case 66: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 287: case 288: case 289: case 290: case 291: case 292: return true; break;} break;
            case 67: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 287: case 288: case 289: case 290: case 291: case 292: return true; break;} break;
            case 68: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 287: case 288: case 289: case 290: case 291: case 292: return true; break;} break;
            case 69: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 287: case 288: case 289: case 290: case 291: case 292: return true; break;} break;
            case 70: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 189: case 190: case 191: case 192: case 193: case 194: case 195: case 196: case 197: case 198: case 199: case 200: case 240: case 241: case 242: case 243: case 244: case 245: case 246: return true; break;} break;
            case 71: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 200: case 240: case 241: case 242: case 243: case 244: case 245: case 246: return true; break;} break;
            case 72: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 197: case 198: case 199: case 200: case 240: case 241: case 242: case 243: case 244: case 245: case 246: return true; break;} break;
            case 73: switch(pc1zBin) {case 99: case 100: case 101: case 102: case 103: case 104: case 105: case 106: case 107: case 108: case 109: case 110: case 111: case 112: case 113: case 114: case 115: case 116: case 117: case 118: case 119: case 120: case 121: case 122: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 189: case 190: case 191: case 192: case 193: case 194: case 195: case 196: case 197: case 198: case 199: case 200: case 201: case 202: case 203: case 204: case 205: case 206: case 207: case 208: case 209: case 210: case 211: case 212: case 213: case 214: case 215: case 216: case 217: case 218: case 219: case 220: case 221: case 273: case 274: case 275: case 276: case 277: case 278: case 279: case 280: case 281: case 282: return true; break;} break;
            case 74: switch(pc1zBin) {case 99: case 100: case 101: case 102: case 103: case 104: case 105: case 106: case 107: case 108: case 109: case 110: case 111: case 112: case 113: case 114: case 115: case 116: case 117: case 118: case 119: case 120: case 121: case 122: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 189: case 190: case 191: case 192: case 193: case 194: case 195: case 196: case 197: case 198: case 199: case 200: case 201: case 202: case 203: case 204: case 205: case 206: case 207: case 208: case 209: case 210: case 211: case 212: case 213: case 214: case 215: case 216: case 217: case 218: case 219: case 220: case 221: case 273: case 274: case 275: case 276: case 277: case 278: case 279: case 280: case 281: case 282: return true; break;} break;
            case 75: switch(pc1zBin) {case 99: case 100: case 101: case 102: case 103: case 105: case 121: case 122: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 197: case 198: case 199: case 200: case 201: case 202: case 209: case 210: case 211: case 212: case 213: case 214: case 215: case 216: case 217: case 218: case 219: case 220: case 221: case 273: case 274: case 275: case 276: case 277: case 278: case 279: case 280: case 281: case 282: return true; break;} break;
            case 76: switch(pc1zBin) {case 99: case 100: case 101: case 102: case 103: case 104: case 105: case 107: case 108: case 109: case 110: case 111: case 112: case 113: case 114: case 115: case 116: case 119: case 120: case 121: case 122: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 197: case 198: case 201: case 208: case 209: case 210: case 211: case 212: case 213: case 214: case 215: case 216: case 217: case 218: case 219: case 220: case 221: case 273: case 274: case 275: case 276: case 277: case 278: case 279: case 280: case 281: case 282: return true; break;} break;
            case 77: switch(pc1zBin) {case 67: case 68: case 69: case 70: case 71: case 72: case 99: case 100: case 101: case 102: case 103: case 104: case 105: case 106: case 107: case 108: case 109: case 110: case 111: case 112: case 113: case 114: case 115: case 116: case 117: case 118: case 119: case 120: case 121: case 122: case 150: case 151: case 152: case 153: case 154: case 155: case 156: case 157: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 197: case 198: case 199: case 200: case 201: case 202: case 203: case 204: case 205: case 206: case 207: case 208: case 209: case 210: case 211: case 212: case 213: case 214: case 215: case 216: case 217: case 218: case 219: case 220: case 221: case 273: case 274: case 275: case 276: case 277: case 278: case 279: case 280: case 281: case 282: return true; break;} break;
            case 78: switch(pc1zBin) {case 67: case 68: case 69: case 70: case 71: case 72: case 150: case 151: case 152: case 153: case 154: case 155: case 156: case 157: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 79: switch(pc1zBin) {case 67: case 68: case 69: case 70: case 71: case 72: case 150: case 151: case 152: case 153: case 154: case 155: case 156: case 157: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 80: switch(pc1zBin) {case 67: case 68: case 69: case 70: case 71: case 72: case 150: case 151: case 152: case 153: case 154: case 155: case 156: case 157: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 81: switch(pc1zBin) {case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 21: case 67: case 68: case 69: case 70: case 71: case 72: case 150: case 151: case 152: case 153: case 154: case 155: case 156: case 157: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 82: switch(pc1zBin) {case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 21: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 83: switch(pc1zBin) {case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 21: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 84: switch(pc1zBin) {case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 21: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 85: switch(pc1zBin) {case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 21: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 86: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 224: case 225: case 226: case 227: case 228: case 229: case 230: case 231: case 232: case 233: case 234: case 235: case 236: case 237: case 238: case 239: case 240: case 241: return true; break;} break;
            case 87: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 224: case 225: case 226: case 227: case 228: case 229: case 230: case 231: case 232: case 233: case 234: case 235: case 236: case 237: case 238: case 239: case 240: case 241: return true; break;} break;
            case 88: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 224: case 225: case 226: case 227: case 228: case 229: case 230: case 231: case 232: case 233: case 234: case 235: case 236: case 237: case 238: case 239: case 240: case 241: return true; break;} break;
            case 89: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 224: case 225: case 226: case 227: case 228: case 229: case 230: case 231: case 232: case 233: case 234: case 235: case 236: case 237: case 238: case 239: case 240: case 241: return true; break;} break;
            case 90: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 91: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 92: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 223: case 224: case 225: case 226: case 227: case 228: case 229: case 230: case 231: return true; break;} break;
            case 93: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 199: case 200: case 201: case 202: case 203: case 204: case 205: case 223: case 224: case 225: case 226: case 227: case 228: case 229: case 230: case 231: return true; break;} break;
            case 94: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 199: case 200: case 201: case 202: case 203: case 204: case 205: case 223: case 224: case 225: case 226: case 227: case 228: case 229: case 230: case 231: return true; break;} break;
            case 95: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 199: case 200: case 201: case 202: case 203: case 204: case 205: case 223: case 224: case 225: case 226: case 227: case 228: case 229: case 230: case 231: return true; break;} break;
            case 96: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 199: case 200: case 201: case 202: case 203: case 204: case 205: case 223: case 224: case 225: case 226: case 227: case 228: case 229: case 230: case 231: return true; break;} break;
            case 97: switch(pc1zBin) {case 140: case 141: case 142: case 143: case 144: case 145: case 146: case 147: case 148: case 149: case 150: case 151: case 152: case 153: case 154: case 155: case 156: case 157: case 158: case 159: case 160: case 161: case 162: case 163: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 223: case 224: case 225: case 226: case 227: case 228: case 229: case 230: case 231: return true; break;} break;
            case 98: switch(pc1zBin) {case 140: case 141: case 142: case 143: case 144: case 145: case 146: case 147: case 148: case 149: case 150: case 151: case 152: case 153: case 155: case 156: case 157: case 158: case 159: case 160: case 161: case 162: case 163: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 99: switch(pc1zBin) {case 140: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 100: switch(pc1zBin) {case 119: case 120: case 121: case 122: case 123: case 124: case 125: case 126: case 127: case 128: case 129: case 130: case 131: case 132: case 133: case 134: case 135: case 136: case 137: case 138: case 139: case 140: case 160: case 161: case 162: case 163: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 101: switch(pc1zBin) {case 18: case 19: case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27: case 28: case 29: case 30: case 31: case 32: case 33: case 34: case 35: case 36: case 37: case 38: case 39: case 40: case 41: case 119: case 120: case 121: case 122: case 123: case 124: case 125: case 126: case 127: case 128: case 129: case 130: case 131: case 132: case 133: case 134: case 135: case 136: case 137: case 138: case 139: case 140: case 143: case 144: case 145: case 146: case 147: case 148: case 149: case 150: case 151: case 152: case 153: case 154: case 155: case 156: case 157: case 158: case 159: case 160: case 161: case 162: case 163: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 299: case 300: case 301: case 302: case 303: case 304: case 305: case 306: case 307: case 308: case 309: case 310: case 311: case 312: case 313: case 314: case 315: case 316: case 317: case 318: case 319: case 320: case 321: case 322: return true; break;} break;
            case 102: switch(pc1zBin) {case 18: case 19: case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27: case 28: case 29: case 30: case 31: case 32: case 33: case 34: case 35: case 36: case 37: case 38: case 39: case 40: case 41: case 119: case 120: case 121: case 122: case 123: case 124: case 125: case 126: case 127: case 128: case 129: case 130: case 131: case 132: case 133: case 134: case 135: case 136: case 137: case 138: case 139: case 140: case 143: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 299: case 300: case 301: case 302: case 303: case 304: case 305: case 306: case 307: case 308: case 309: case 310: case 311: case 312: case 313: case 314: case 315: case 316: case 317: case 318: case 319: case 320: case 321: case 322: return true; break;} break;
            case 103: switch(pc1zBin) {case 18: case 19: case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27: case 28: case 29: case 30: case 31: case 32: case 33: case 34: case 35: case 36: case 37: case 38: case 39: case 40: case 41: case 119: case 120: case 121: case 122: case 142: case 143: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 299: case 300: case 301: case 302: case 303: case 304: case 305: case 306: case 307: case 308: case 309: case 310: case 311: case 312: case 313: case 314: case 315: case 316: case 317: case 318: case 319: case 320: case 321: case 322: return true; break;} break;
            case 104: switch(pc1zBin) {case 18: case 19: case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27: case 28: case 29: case 30: case 31: case 32: case 33: case 34: case 35: case 36: case 37: case 38: case 39: case 40: case 41: case 119: case 120: case 121: case 125: case 126: case 128: case 133: case 137: case 138: case 139: case 140: case 141: case 142: case 143: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 299: case 300: case 301: case 302: case 303: case 304: case 305: case 306: case 307: case 308: case 309: case 310: case 311: case 312: case 313: case 314: case 315: case 316: case 317: case 318: case 319: case 320: case 321: case 322: return true; break;} break;
            case 105: switch(pc1zBin) {case 18: case 19: case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27: case 28: case 29: case 30: case 31: case 32: case 33: case 34: case 35: case 36: case 37: case 38: case 39: case 40: case 41: case 119: case 120: case 121: case 122: case 123: case 124: case 125: case 126: case 127: case 128: case 129: case 130: case 131: case 132: case 133: case 134: case 135: case 136: case 137: case 138: case 139: case 140: case 141: case 142: case 143: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 299: case 300: case 301: case 302: case 303: case 304: case 305: case 306: case 307: case 308: case 309: case 310: case 311: case 312: case 313: case 314: case 315: case 316: case 317: case 318: case 319: case 320: case 321: case 322: return true; break;} break;
            case 106: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 107: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 108: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 109: switch(pc1zBin) {case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 21: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 110: switch(pc1zBin) {case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 21: case 33: case 34: case 35: case 36: case 37: case 38: case 39: case 40: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 111: switch(pc1zBin) {case 20: case 21: case 33: case 34: case 35: case 36: case 37: case 38: case 39: case 40: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 112: switch(pc1zBin) {case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 21: case 33: case 34: case 35: case 36: case 37: case 38: case 39: case 40: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 113: switch(pc1zBin) {case 33: case 34: case 35: case 36: case 37: case 38: case 39: case 40: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 287: case 288: case 289: case 290: case 291: case 292: return true; break;} break;
            case 114: switch(pc1zBin) {case 33: case 34: case 35: case 36: case 37: case 38: case 39: case 40: case 54: case 55: case 56: case 57: case 58: case 59: case 60: case 81: case 82: case 83: case 84: case 85: case 86: case 87: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 287: case 288: case 289: case 290: case 291: case 292: return true; break;} break;
            case 115: switch(pc1zBin) {case 54: case 55: case 56: case 57: case 58: case 59: case 60: case 81: case 82: case 83: case 84: case 85: case 86: case 87: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 116: switch(pc1zBin) {case 54: case 55: case 56: case 57: case 58: case 59: case 60: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 117: switch(pc1zBin) {case 54: case 55: case 56: case 57: case 58: case 59: case 60: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 118: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 119: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 120: switch(pc1zBin) {case 140: case 141: case 142: case 143: case 144: case 145: case 146: case 147: case 148: case 149: case 150: case 151: case 152: case 153: case 154: case 155: case 156: case 157: case 158: case 159: case 160: case 161: case 162: case 163: case 164: case 165: case 166: case 167: case 168: case 169: case 170: case 171: case 172: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 121: switch(pc1zBin) {case 140: case 141: case 142: case 143: case 144: case 145: case 146: case 147: case 148: case 149: case 150: case 151: case 152: case 153: case 154: case 155: case 156: case 157: case 158: case 159: case 160: case 161: case 162: case 163: case 164: case 165: case 166: case 167: case 168: case 169: case 170: case 171: case 172: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 258: case 259: case 260: case 261: case 262: case 263: case 264: case 265: case 266: case 267: case 268: case 269: case 270: case 271: case 272: case 273: case 274: case 275: case 276: case 277: case 278: case 279: case 280: case 281: return true; break;} break;
            case 122: switch(pc1zBin) {case 140: case 141: case 142: case 143: case 144: case 145: case 146: case 147: case 148: case 149: case 150: case 151: case 152: case 153: case 154: case 155: case 156: case 157: case 158: case 159: case 160: case 161: case 162: case 163: case 164: case 165: case 166: case 167: case 168: case 169: case 170: case 171: case 172: case 173: case 174: case 175: case 176: case 177: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 258: case 259: case 260: case 261: case 262: case 263: case 264: case 265: case 266: case 267: case 268: case 269: case 270: case 271: case 272: case 273: case 274: case 275: case 276: case 277: case 278: case 279: case 280: case 281: case 309: case 310: case 311: case 312: case 313: case 314: return true; break;} break;
            case 123: switch(pc1zBin) {case 140: case 141: case 142: case 143: case 144: case 145: case 146: case 147: case 148: case 149: case 150: case 151: case 152: case 153: case 154: case 155: case 156: case 157: case 158: case 159: case 160: case 161: case 162: case 163: case 164: case 165: case 166: case 167: case 168: case 169: case 170: case 171: case 172: case 173: case 174: case 175: case 176: case 177: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 258: case 259: case 260: case 263: case 267: case 275: case 276: case 278: case 279: case 280: case 281: case 309: case 310: case 311: case 312: case 313: case 314: return true; break;} break;
            case 124: switch(pc1zBin) {case 140: case 141: case 142: case 143: case 144: case 145: case 146: case 147: case 148: case 149: case 150: case 151: case 152: case 153: case 154: case 155: case 156: case 157: case 158: case 159: case 160: case 161: case 162: case 163: case 164: case 165: case 166: case 167: case 168: case 169: case 170: case 171: case 172: case 173: case 174: case 175: case 176: case 177: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 258: case 259: case 260: case 261: case 262: case 263: case 264: case 265: case 266: case 267: case 268: case 269: case 270: case 271: case 272: case 273: case 274: case 275: case 276: case 277: case 278: case 279: case 280: case 281: case 309: case 310: case 311: case 312: case 313: case 314: return true; break;} break;
            case 125: switch(pc1zBin) {case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 21: case 91: case 92: case 93: case 94: case 95: case 140: case 141: case 142: case 143: case 144: case 145: case 146: case 147: case 148: case 149: case 150: case 151: case 152: case 153: case 154: case 155: case 156: case 157: case 158: case 159: case 160: case 161: case 162: case 163: case 164: case 165: case 166: case 167: case 168: case 169: case 170: case 171: case 172: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 258: case 259: case 260: case 261: case 262: case 263: case 264: case 265: case 266: case 267: case 268: case 269: case 270: case 271: case 272: case 273: case 274: case 275: case 276: case 277: case 278: case 279: case 280: case 281: case 309: case 310: case 311: case 312: case 313: case 314: return true; break;} break;
            case 126: switch(pc1zBin) {case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 21: case 91: case 92: case 93: case 94: case 95: case 140: case 141: case 142: case 143: case 144: case 145: case 146: case 147: case 148: case 149: case 150: case 151: case 152: case 153: case 154: case 155: case 156: case 157: case 158: case 159: case 160: case 161: case 162: case 163: case 164: case 165: case 166: case 167: case 168: case 169: case 170: case 171: case 172: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 127: switch(pc1zBin) {case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 21: case 91: case 92: case 93: case 94: case 95: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 128: switch(pc1zBin) {case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 21: case 91: case 92: case 93: case 94: case 95: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 129: switch(pc1zBin) {case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 21: case 126: case 127: case 128: case 129: case 130: case 131: case 132: case 133: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 130: switch(pc1zBin) {case 126: case 127: case 128: case 129: case 130: case 131: case 132: case 133: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 131: switch(pc1zBin) {case 126: case 127: case 128: case 129: case 130: case 131: case 132: case 133: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 132: switch(pc1zBin) {case 126: case 127: case 128: case 129: case 130: case 131: case 132: case 133: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 133: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 134: switch(pc1zBin) {case 74: case 75: case 76: case 77: case 78: case 79: case 80: case 81: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 135: switch(pc1zBin) {case 74: case 75: case 76: case 77: case 78: case 79: case 80: case 81: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 136: switch(pc1zBin) {case 74: case 75: case 76: case 77: case 78: case 79: case 80: case 81: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 137: switch(pc1zBin) {case 86: case 87: case 88: case 89: case 90: case 91: case 92: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 138: switch(pc1zBin) {case 86: case 87: case 88: case 89: case 90: case 91: case 92: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 139: switch(pc1zBin) {case 86: case 87: case 88: case 89: case 90: case 91: case 92: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 140: switch(pc1zBin) {case 86: case 87: case 88: case 89: case 90: case 91: case 92: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 141: switch(pc1zBin) {case 86: case 87: case 88: case 89: case 90: case 91: case 92: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 142: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 143: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 144: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 145: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 146: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 147: switch(pc1zBin) {case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27: case 28: case 29: case 30: case 31: case 32: case 33: case 34: case 35: case 36: case 37: case 38: case 39: case 40: case 41: case 42: case 43: case 44: case 45: case 46: case 47: case 48: case 49: case 50: case 51: case 52: case 53: case 54: case 55: case 56: case 57: case 58: case 59: case 60: case 61: case 62: case 63: case 64: case 65: case 66: case 67: case 68: case 69: case 70: case 71: case 72: case 73: case 74: case 75: case 76: case 77: case 78: case 79: case 80: case 81: case 82: case 83: case 84: case 85: case 86: case 87: case 88: case 89: case 90: case 91: case 92: case 93: case 94: case 95: case 96: case 97: case 98: case 99: case 100: case 101: case 102: case 103: case 104: case 105: case 106: case 107: case 108: case 109: case 110: case 111: case 112: case 113: case 114: case 115: case 116: case 117: case 118: case 119: case 120: case 121: case 122: case 123: case 124: case 125: case 126: case 127: case 128: case 129: case 130: case 131: case 132: case 133: case 134: case 135: case 136: case 137: case 138: case 139: case 140: case 141: case 142: case 143: case 144: case 145: case 146: case 147: case 148: case 149: case 150: case 151: case 152: case 153: case 154: case 155: case 156: case 157: case 158: case 159: case 160: case 161: case 162: case 163: case 164: case 165: case 166: case 167: case 168: case 169: case 170: case 171: case 172: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 148: switch(pc1zBin) {case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27: case 28: case 29: case 30: case 31: case 32: case 33: case 34: case 35: case 36: case 37: case 38: case 39: case 40: case 41: case 42: case 43: case 44: case 45: case 46: case 47: case 48: case 49: case 50: case 51: case 52: case 53: case 54: case 55: case 56: case 57: case 58: case 59: case 60: case 61: case 62: case 63: case 64: case 65: case 66: case 67: case 68: case 69: case 70: case 71: case 72: case 73: case 74: case 75: case 76: case 77: case 78: case 79: case 80: case 81: case 82: case 83: case 84: case 85: case 86: case 87: case 88: case 89: case 90: case 91: case 92: case 93: case 94: case 95: case 96: case 97: case 98: case 99: case 100: case 101: case 102: case 103: case 104: case 105: case 106: case 107: case 108: case 109: case 110: case 111: case 112: case 113: case 114: case 115: case 116: case 117: case 118: case 119: case 120: case 121: case 122: case 123: case 124: case 125: case 126: case 127: case 128: case 129: case 130: case 131: case 132: case 133: case 134: case 135: case 136: case 137: case 138: case 139: case 140: case 141: case 142: case 143: case 144: case 145: case 146: case 147: case 148: case 149: case 150: case 151: case 152: case 153: case 154: case 155: case 156: case 157: case 158: case 159: case 160: case 161: case 162: case 163: case 164: case 165: case 166: case 167: case 168: case 169: case 170: case 171: case 172: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 149: switch(pc1zBin) {case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27: case 28: case 29: case 30: case 31: case 32: case 33: case 34: case 35: case 36: case 37: case 38: case 39: case 40: case 41: case 42: case 43: case 44: case 45: case 46: case 47: case 48: case 49: case 50: case 51: case 52: case 53: case 54: case 55: case 56: case 57: case 58: case 59: case 60: case 61: case 62: case 63: case 64: case 65: case 66: case 67: case 68: case 69: case 70: case 71: case 72: case 73: case 74: case 75: case 76: case 77: case 78: case 79: case 80: case 81: case 82: case 83: case 84: case 85: case 86: case 87: case 88: case 89: case 90: case 91: case 92: case 93: case 94: case 95: case 96: case 97: case 98: case 99: case 100: case 101: case 102: case 103: case 104: case 105: case 106: case 107: case 108: case 109: case 110: case 111: case 112: case 113: case 114: case 115: case 116: case 117: case 118: case 119: case 120: case 121: case 122: case 123: case 124: case 125: case 126: case 127: case 128: case 129: case 130: case 131: case 132: case 133: case 134: case 135: case 136: case 137: case 138: case 139: case 140: case 141: case 142: case 143: case 144: case 145: case 146: case 147: case 148: case 149: case 150: case 151: case 152: case 153: case 154: case 155: case 156: case 157: case 158: case 159: case 160: case 161: case 163: case 164: case 165: case 166: case 167: case 168: case 169: case 170: case 171: case 172: case 173: case 174: case 175: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 150: switch(pc1zBin) {case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27: case 28: case 29: case 30: case 31: case 32: case 33: case 34: case 35: case 36: case 37: case 38: case 39: case 40: case 41: case 42: case 43: case 44: case 45: case 46: case 47: case 48: case 49: case 50: case 51: case 52: case 53: case 54: case 55: case 56: case 57: case 58: case 59: case 60: case 61: case 62: case 63: case 64: case 65: case 66: case 67: case 68: case 69: case 70: case 71: case 72: case 73: case 74: case 75: case 76: case 77: case 78: case 79: case 80: case 81: case 82: case 83: case 84: case 85: case 86: case 87: case 88: case 89: case 90: case 91: case 92: case 93: case 94: case 95: case 96: case 97: case 98: case 99: case 100: case 101: case 102: case 103: case 104: case 105: case 106: case 107: case 108: case 109: case 110: case 111: case 112: case 113: case 114: case 115: case 116: case 117: case 118: case 119: case 120: case 121: case 122: case 123: case 124: case 125: case 126: case 127: case 128: case 129: case 130: case 131: case 132: case 133: case 134: case 135: case 136: case 137: case 138: case 139: case 140: case 141: case 142: case 143: case 144: case 145: case 146: case 147: case 148: case 149: case 150: case 151: case 152: case 153: case 154: case 155: case 156: case 157: case 158: case 159: case 160: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 151: switch(pc1zBin) {case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27: case 28: case 29: case 30: case 31: case 32: case 33: case 34: case 35: case 36: case 37: case 38: case 39: case 40: case 41: case 42: case 43: case 44: case 45: case 46: case 47: case 48: case 49: case 50: case 51: case 52: case 53: case 54: case 55: case 56: case 57: case 58: case 59: case 60: case 61: case 62: case 63: case 64: case 65: case 66: case 67: case 68: case 69: case 70: case 71: case 72: case 73: case 74: case 75: case 76: case 77: case 78: case 79: case 80: case 81: case 82: case 83: case 84: case 85: case 86: case 87: case 88: case 89: case 90: case 91: case 92: case 93: case 94: case 95: case 96: case 97: case 98: case 99: case 100: case 101: case 102: case 103: case 104: case 105: case 106: case 107: case 108: case 109: case 110: case 111: case 112: case 113: case 114: case 115: case 116: case 117: case 118: case 119: case 120: case 121: case 122: case 123: case 124: case 125: case 126: case 127: case 128: case 129: case 130: case 131: case 132: case 133: case 134: case 135: case 136: case 137: case 138: case 139: case 140: case 141: case 142: case 143: case 144: case 145: case 146: case 147: case 148: case 149: case 150: case 151: case 152: case 153: case 154: case 155: case 156: case 157: case 158: case 159: case 160: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 152: switch(pc1zBin) {case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27: case 28: case 29: case 30: case 31: case 32: case 33: case 34: case 35: case 36: case 37: case 38: case 39: case 40: case 41: case 42: case 43: case 44: case 45: case 46: case 47: case 48: case 49: case 50: case 51: case 52: case 53: case 54: case 55: case 56: case 57: case 58: case 59: case 60: case 61: case 62: case 63: case 64: case 65: case 66: case 67: case 68: case 69: case 70: case 71: case 72: case 73: case 74: case 75: case 76: case 77: case 78: case 79: case 80: case 81: case 82: case 83: case 84: case 85: case 86: case 87: case 88: case 89: case 90: case 91: case 92: case 93: case 94: case 95: case 96: case 97: case 98: case 99: case 100: case 101: case 102: case 103: case 104: case 105: case 106: case 107: case 108: case 109: case 110: case 111: case 112: case 113: case 114: case 115: case 116: case 117: case 118: case 119: case 120: case 121: case 122: case 123: case 124: case 125: case 126: case 127: case 128: case 129: case 130: case 131: case 132: case 133: case 134: case 135: case 136: case 137: case 138: case 139: case 140: case 141: case 142: case 143: case 144: case 145: case 146: case 147: case 148: case 149: case 150: case 151: case 152: case 153: case 154: case 155: case 156: case 157: case 158: case 159: case 160: case 161: case 162: case 163: case 164: case 165: case 166: case 167: case 168: case 169: case 170: case 171: case 172: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 153: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 154: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 155: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 156: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 258: case 259: case 260: case 261: case 262: case 263: case 264: case 265: case 266: case 267: case 268: case 269: case 270: case 271: case 272: case 273: case 274: case 275: case 276: case 277: case 278: case 279: case 280: case 281: case 282: return true; break;} break;
            case 157: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 258: case 259: case 260: case 261: case 262: case 263: case 264: case 265: case 266: case 267: case 268: case 269: case 270: case 271: case 272: case 273: case 274: case 275: case 276: case 277: case 278: case 279: case 280: case 281: case 282: return true; break;} break;
            case 158: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 258: case 259: case 260: case 261: case 262: case 263: case 264: case 265: case 266: case 267: case 268: case 269: case 270: case 271: case 272: case 273: case 274: case 275: case 276: case 277: case 278: case 279: case 280: case 281: case 282: return true; break;} break;
            case 159: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 258: case 259: case 260: case 261: case 262: case 263: case 264: case 265: case 266: case 267: case 268: case 269: case 270: case 271: case 272: case 273: case 274: case 275: case 276: case 277: case 278: case 279: case 280: case 281: case 282: return true; break;} break;
            case 160: switch(pc1zBin) {case 79: case 80: case 81: case 82: case 83: case 84: case 85: case 86: case 87: case 88: case 89: case 90: case 91: case 92: case 93: case 94: case 95: case 96: case 97: case 98: case 99: case 100: case 101: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 258: case 259: case 260: case 261: case 262: case 263: case 264: case 265: case 266: case 267: case 268: case 269: case 270: case 271: case 272: case 273: case 274: case 275: case 276: case 277: case 278: case 279: case 280: case 281: case 282: return true; break;} break;
            case 161: switch(pc1zBin) {case 79: case 80: case 81: case 82: case 83: case 84: case 85: case 86: case 87: case 88: case 89: case 90: case 91: case 92: case 93: case 94: case 95: case 96: case 97: case 98: case 99: case 100: case 101: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 162: switch(pc1zBin) {case 79: case 80: case 81: case 82: case 83: case 84: case 85: case 86: case 87: case 88: case 89: case 90: case 91: case 92: case 93: case 94: case 95: case 96: case 97: case 98: case 101: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 163: switch(pc1zBin) {case 79: case 80: case 81: case 82: case 100: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 164: switch(pc1zBin) {case 79: case 80: case 81: case 82: case 83: case 84: case 85: case 86: case 87: case 88: case 89: case 90: case 91: case 92: case 93: case 94: case 95: case 96: case 97: case 98: case 99: case 100: case 101: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: return true; break;} break;
            case 165: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 189: case 190: case 191: case 192: case 193: case 218: case 219: case 220: case 221: case 222: case 223: case 224: case 225: case 226: return true; break;} break;
            case 166: switch(pc1zBin) {case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 189: case 190: case 191: case 192: case 193: case 218: case 219: case 220: case 221: case 222: case 223: case 224: case 225: case 226: return true; break;} break;
         }
      }
      else
      {
         if (pc1z < -85 || pc1phi < -0.48 ||
             pc1z > 83.5 || pc1phi > 0.98) return true; 
         
         const short pc1zBin = static_cast<short>((pc1z + 90.)*2.) + 1;
         const short pc1phiBin = static_cast<short>((pc1phi + 0.55)*100.) + 1;
         
         switch(pc1phiBin)
         {
            case 8: switch(pc1zBin) {case 59: case 60: case 61: case 62: case 63: case 64: case 65: case 66: case 67: case 68: case 69: case 70: case 71: case 72: case 73: case 74: case 75: case 76: case 77: case 78: case 79: case 80: case 81: case 82: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 9: switch(pc1zBin) {case 59: case 60: case 61: case 62: case 63: case 64: case 65: case 66: case 67: case 68: case 69: case 70: case 71: case 72: case 73: case 74: case 75: case 76: case 77: case 78: case 79: case 80: case 81: case 82: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 10: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 11: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 12: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 13: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 14: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 15: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 16: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 17: switch(pc1zBin) {case 120: case 121: case 122: case 123: case 124: case 125: case 126: case 127: case 128: case 129: case 130: case 131: case 132: case 133: case 134: case 135: case 136: case 137: case 138: case 139: case 140: case 141: case 142: case 143: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 18: switch(pc1zBin) {case 120: case 121: case 122: case 123: case 124: case 125: case 126: case 127: case 128: case 129: case 130: case 131: case 132: case 133: case 134: case 135: case 136: case 137: case 139: case 140: case 142: case 143: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 19: switch(pc1zBin) {case 120: case 142: case 143: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 20: switch(pc1zBin) {case 120: case 121: case 122: case 123: case 124: case 125: case 126: case 127: case 128: case 129: case 130: case 131: case 132: case 133: case 134: case 135: case 136: case 137: case 138: case 139: case 140: case 141: case 142: case 143: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 21: switch(pc1zBin) {case 120: case 121: case 122: case 123: case 124: case 125: case 126: case 127: case 128: case 129: case 130: case 131: case 132: case 133: case 134: case 135: case 136: case 137: case 138: case 139: case 140: case 141: case 142: case 143: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 22: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 23: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 24: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 25: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 26: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 304: case 305: case 306: case 307: case 308: case 309: case 310: case 311: case 312: return true; break;} break;
            case 27: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 304: case 305: case 306: case 307: case 308: case 309: case 310: case 311: case 312: return true; break;} break;
            case 28: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 29: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 189: case 190: case 191: case 192: case 193: case 194: case 195: case 196: case 197: case 198: case 199: case 200: case 201: return true; break;} break;
            case 30: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 200: case 201: return true; break;} break;
            case 31: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 201: return true; break;} break;
            case 32: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 183: case 184: case 185: case 186: case 187: case 188: case 189: case 190: case 191: case 192: case 193: case 194: case 195: case 196: case 197: case 198: case 199: case 200: case 201: return true; break;} break;
            case 33: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 189: case 190: case 191: case 192: case 193: case 194: case 195: case 196: case 197: case 198: case 199: case 200: case 201: return true; break;} break;
            case 34: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 35: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 36: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 37: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 38: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 39: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 40: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 41: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 299: case 300: case 301: case 302: case 303: case 304: case 305: case 306: case 307: case 308: case 309: case 310: case 311: case 312: case 313: case 314: case 315: case 316: case 317: case 318: case 319: case 320: case 321: case 322: return true; break;} break;
            case 42: switch(pc1zBin) {case 93: case 94: case 95: case 96: case 97: case 98: case 99: case 100: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 299: case 300: case 301: case 302: case 303: case 304: case 305: case 306: case 307: case 308: case 309: case 310: case 311: case 312: case 313: case 314: case 315: case 316: case 317: case 318: case 319: case 320: case 321: case 322: return true; break;} break;
            case 43: switch(pc1zBin) {case 93: case 94: case 95: case 96: case 97: case 98: case 99: case 100: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 299: case 300: case 301: case 302: case 303: case 304: case 305: case 306: case 307: case 308: case 309: case 310: case 311: case 312: case 313: case 314: case 315: case 316: case 317: case 318: case 319: case 320: case 321: case 322: return true; break;} break;
            case 44: switch(pc1zBin) {case 93: case 94: case 95: case 96: case 97: case 98: case 99: case 100: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 299: case 300: case 301: case 302: case 303: case 304: case 305: case 306: case 307: case 308: case 309: case 310: case 311: case 312: case 313: case 314: case 315: case 316: case 317: case 318: case 319: case 320: case 321: case 322: return true; break;} break;
            case 45: switch(pc1zBin) {case 93: case 94: case 95: case 96: case 97: case 98: case 99: case 100: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 46: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 47: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 48: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 49: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 50: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 51: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 52: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 53: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 54: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 55: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 56: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 57: switch(pc1zBin) {case 11: case 12: case 13: case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 279: case 280: case 281: case 282: case 283: case 284: case 285: case 286: case 287: case 288: case 289: case 290: case 291: case 292: case 293: case 294: case 295: case 296: case 297: case 298: case 299: case 300: return true; break;} break;
            case 58: switch(pc1zBin) {case 11: case 15: case 16: case 17: case 18: case 19: case 20: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 279: case 280: case 281: case 282: case 283: case 284: case 285: case 286: case 287: case 288: case 289: case 290: case 291: case 292: case 293: case 294: case 295: case 296: case 297: case 298: case 299: case 300: return true; break;} break;
            case 59: switch(pc1zBin) {case 11: case 12: case 13: case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 279: case 280: case 281: case 282: case 283: case 284: case 285: case 286: case 287: case 288: case 289: case 290: case 291: case 292: case 293: case 294: case 295: case 296: case 297: case 298: case 299: case 300: return true; break;} break;
            case 60: switch(pc1zBin) {case 11: case 12: case 13: case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 279: case 280: case 281: case 282: case 283: case 284: case 285: case 286: case 287: case 288: case 289: case 290: case 291: case 292: case 293: case 294: case 295: case 296: case 297: case 298: case 299: case 300: return true; break;} break;
            case 61: switch(pc1zBin) {case 11: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 62: switch(pc1zBin) {case 11: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 63: switch(pc1zBin) {case 11: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 64: switch(pc1zBin) {case 11: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 65: return true; break;
            case 66: return true; break;
            case 67: return true; break;
            case 68: switch(pc1zBin) {case 79: case 80: case 81: case 82: case 83: case 84: case 85: case 86: case 87: case 88: case 89: case 90: case 91: case 92: case 93: case 94: case 95: case 96: case 97: case 98: case 99: case 100: case 101: case 102: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 69: switch(pc1zBin) {case 79: case 80: case 81: case 82: case 83: case 84: case 85: case 86: case 87: case 89: case 90: case 91: case 92: case 93: case 94: case 95: case 96: case 97: case 98: case 99: case 100: case 101: case 102: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 70: switch(pc1zBin) {case 79: case 80: case 81: case 82: case 102: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 189: case 190: case 191: case 192: case 193: case 194: case 195: case 196: case 197: case 198: case 199: case 200: case 201: case 202: case 203: case 204: case 205: case 206: case 207: case 208: case 209: case 210: case 211: case 212: case 213: case 214: case 215: case 216: case 217: case 218: case 219: case 220: case 221: case 222: case 223: case 224: case 225: case 226: case 227: case 228: case 229: case 230: case 231: case 232: case 233: case 234: case 235: case 236: case 237: case 238: case 239: case 240: case 241: case 242: case 243: case 244: case 245: case 246: case 247: case 248: case 249: case 250: case 251: case 252: case 253: case 254: case 255: case 256: case 257: case 258: case 259: case 260: case 261: case 262: case 263: case 264: case 265: case 266: case 267: case 268: case 269: case 270: case 271: case 272: case 273: case 274: case 275: case 276: case 277: case 278: case 279: case 280: case 281: case 282: case 283: case 284: case 285: case 286: case 287: case 288: case 289: case 290: case 291: case 292: case 293: case 294: case 295: case 296: case 297: case 298: case 299: case 300: case 301: case 302: case 303: case 304: case 305: case 306: case 307: case 308: case 309: case 310: case 311: case 312: case 313: case 314: case 315: case 316: case 317: case 318: case 319: case 320: case 321: case 322: case 323: case 324: case 325: case 326: case 327: case 328: case 329: case 330: case 331: case 332: case 333: case 334: case 335: case 336: case 337: case 338: case 339: case 340: case 341: case 342: case 343: case 344: case 345: case 346: case 347: return true; break;} break;
            case 71: switch(pc1zBin) {case 79: case 80: case 81: case 85: case 95: case 96: case 99: case 100: case 101: case 102: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 189: case 190: case 191: case 192: case 193: case 194: case 195: case 196: case 197: case 198: case 199: case 200: case 201: case 202: case 203: case 204: case 205: case 206: case 207: case 208: case 209: case 210: case 211: case 212: case 213: case 214: case 215: case 216: case 217: case 218: case 219: case 220: case 221: case 222: case 223: case 224: case 225: case 226: case 227: case 228: case 229: case 230: case 231: case 232: case 233: case 234: case 235: case 236: case 237: case 238: case 239: case 240: case 241: case 242: case 243: case 244: case 245: case 246: case 247: case 248: case 249: case 250: case 251: case 252: case 253: case 254: case 255: case 256: case 257: case 258: case 259: case 260: case 261: case 262: case 263: case 264: case 265: case 266: case 267: case 268: case 269: case 270: case 271: case 272: case 273: case 274: case 275: case 276: case 277: case 278: case 279: case 280: case 281: case 282: case 283: case 284: case 285: case 286: case 287: case 288: case 289: case 290: case 291: case 292: case 293: case 294: case 295: case 296: case 297: case 298: case 299: case 300: case 301: case 302: case 303: case 304: case 305: case 306: case 307: case 308: case 309: case 310: case 311: case 312: case 313: case 314: case 315: case 316: case 317: case 318: case 319: case 320: case 321: case 322: case 323: case 324: case 325: case 326: case 327: case 328: case 329: case 330: case 331: case 332: case 333: case 334: case 335: case 336: case 337: case 338: case 339: case 340: case 341: case 342: case 343: case 344: case 345: case 346: case 347: return true; break;} break;
            case 72: return true; break;
            case 73: return true; break;
            case 74: return true; break;
            case 75: return true; break;
            case 76: return true; break;
            case 77: return true; break;
            case 78: return true; break;
            case 79: return true; break;
            case 80: return true; break;
            case 81: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 189: case 190: case 191: case 192: case 193: case 194: case 195: case 196: case 197: case 198: case 199: case 200: case 201: case 202: case 203: case 204: case 205: case 206: case 207: case 208: case 209: case 210: case 211: case 212: case 213: case 214: case 215: case 216: case 217: case 218: case 219: case 220: case 221: case 222: case 223: case 224: case 225: case 226: case 227: case 228: case 229: case 230: case 231: case 232: case 233: case 234: case 235: case 236: case 237: case 238: case 239: case 240: case 241: case 242: case 243: case 244: case 245: case 246: case 247: case 248: case 249: case 250: case 251: case 252: case 253: case 254: case 255: return true; break;} break;
            case 82: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 83: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 84: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 85: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 86: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 87: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 88: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 89: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 90: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 91: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 92: switch(pc1zBin) {case 120: case 121: case 122: case 123: case 124: case 125: case 126: case 127: case 128: case 129: case 130: case 131: case 132: case 133: case 134: case 135: case 136: case 137: case 138: case 139: case 140: case 141: case 142: case 143: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 320: case 321: case 322: case 323: case 324: case 325: case 326: case 327: case 328: case 329: case 330: case 331: case 332: case 333: case 334: case 335: case 336: case 337: case 338: case 339: case 340: case 341: case 342: case 343: return true; break;} break;
            case 93: switch(pc1zBin) {case 120: case 121: case 122: case 123: case 124: case 131: case 134: case 143: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 320: case 321: case 322: case 323: case 324: case 325: case 326: case 327: case 328: case 329: case 330: case 331: case 332: case 333: case 334: case 335: case 336: case 337: case 338: case 339: case 340: case 341: case 342: case 343: return true; break;} break;
            case 94: switch(pc1zBin) {case 120: case 143: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 320: case 321: case 322: case 323: case 324: case 325: case 326: case 327: case 328: case 329: case 330: case 331: case 332: case 333: case 334: case 335: case 336: case 337: case 338: case 339: case 340: case 341: case 342: case 343: return true; break;} break;
            case 95: switch(pc1zBin) {case 120: case 121: case 122: case 123: case 124: case 125: case 126: case 127: case 128: case 129: case 130: case 131: case 132: case 133: case 134: case 135: case 136: case 137: case 138: case 139: case 140: case 141: case 142: case 143: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 320: case 321: case 322: case 323: case 324: case 325: case 326: case 327: case 328: case 329: case 330: case 331: case 332: case 333: case 334: case 335: case 336: case 337: case 338: case 339: case 340: case 341: case 342: case 343: return true; break;} break;
            case 96: switch(pc1zBin) {case 120: case 121: case 122: case 123: case 124: case 125: case 126: case 127: case 128: case 129: case 130: case 131: case 132: case 133: case 134: case 135: case 136: case 137: case 138: case 139: case 140: case 141: case 142: case 143: case 174: case 175: case 176: case 177: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 97: switch(pc1zBin) {case 120: case 121: case 122: case 123: case 124: case 125: case 126: case 127: case 128: case 129: case 130: case 131: case 132: case 133: case 134: case 135: case 136: case 137: case 138: case 139: case 140: case 141: case 142: case 143: case 174: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 98: switch(pc1zBin) {case 120: case 121: case 122: case 123: case 124: case 125: case 126: case 127: case 128: case 129: case 130: case 131: case 132: case 133: case 134: case 135: case 136: case 137: case 138: case 139: case 140: case 141: case 142: case 143: case 174: case 175: case 176: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 99: switch(pc1zBin) {case 120: case 121: case 122: case 123: case 124: case 125: case 126: case 127: case 128: case 129: case 130: case 131: case 132: case 133: case 134: case 135: case 136: case 137: case 138: case 139: case 140: case 141: case 142: case 143: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 218: case 219: case 220: case 221: case 222: case 223: case 224: case 225: case 226: case 227: case 228: case 229: case 230: case 231: case 232: case 233: case 234: case 235: case 236: case 237: case 238: case 239: case 240: case 241: case 242: return true; break;} break;
            case 100: switch(pc1zBin) {case 150: case 151: case 152: case 153: case 154: case 155: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 218: case 219: case 220: case 221: case 222: case 223: case 229: case 237: case 238: case 239: case 240: case 241: case 242: return true; break;} break;
            case 101: switch(pc1zBin) {case 150: case 151: case 152: case 153: case 154: case 155: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 218: case 219: case 220: case 239: case 242: return true; break;} break;
            case 102: switch(pc1zBin) {case 150: case 151: case 152: case 153: case 154: case 155: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 218: case 219: case 220: case 228: case 238: case 239: case 240: case 241: case 242: return true; break;} break;
            case 103: switch(pc1zBin) {case 150: case 151: case 152: case 153: case 154: case 155: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 218: case 219: case 220: case 221: case 222: case 223: case 224: case 225: case 226: case 227: case 228: case 229: case 230: case 231: case 232: case 233: case 234: case 235: case 236: case 237: case 238: case 239: case 240: case 241: case 242: return true; break;} break;
            case 104: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 105: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 106: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 107: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 108: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 109: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 110: return true; break;
            case 111: return true; break;
            case 112: return true; break;
            case 113: return true; break;
            case 114: switch(pc1zBin) {case 31: case 35: case 36: case 37: case 38: case 39: case 40: case 41: case 42: case 43: case 44: case 45: case 46: case 47: case 48: case 49: case 50: case 51: case 52: case 53: case 54: case 55: case 57: case 60: case 62: case 64: case 66: case 70: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 115: switch(pc1zBin) {case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 116: switch(pc1zBin) {case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 117: switch(pc1zBin) {case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 118: switch(pc1zBin) {case 45: case 48: case 180: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 119: switch(pc1zBin) {case 29: case 31: case 32: case 33: case 34: case 35: case 36: case 37: case 38: case 39: case 40: case 41: case 42: case 43: case 44: case 45: case 46: case 47: case 48: case 49: case 50: case 51: case 52: case 53: case 54: case 55: case 56: case 57: case 58: case 59: case 60: case 61: case 62: case 64: case 65: case 66: case 70: case 71: case 77: case 98: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 120: switch(pc1zBin) {case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27: case 28: case 29: case 30: case 31: case 32: case 33: case 34: case 35: case 36: case 37: case 38: case 39: case 40: case 41: case 42: case 43: case 44: case 45: case 46: case 47: case 48: case 49: case 50: case 51: case 52: case 53: case 54: case 55: case 56: case 57: case 58: case 59: case 60: case 61: case 62: case 63: case 64: case 65: case 66: case 67: case 68: case 69: case 70: case 71: case 72: case 73: case 74: case 75: case 76: case 77: case 78: case 79: case 80: case 81: case 82: case 83: case 84: case 85: case 86: case 87: case 88: case 89: case 90: case 91: case 92: case 93: case 94: case 95: case 96: case 97: case 98: case 99: case 100: case 101: case 102: case 103: case 106: case 109: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 121: switch(pc1zBin) {case 11: case 12: case 13: case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27: case 28: case 29: case 30: case 31: case 32: case 33: case 34: case 35: case 36: case 37: case 38: case 39: case 40: case 41: case 42: case 43: case 44: case 45: case 46: case 47: case 48: case 49: case 50: case 51: case 52: case 53: case 54: case 55: case 56: case 57: case 58: case 59: case 60: case 61: case 62: case 63: case 64: case 65: case 66: case 67: case 68: case 69: case 70: case 71: case 72: case 73: case 74: case 75: case 76: case 77: case 78: case 79: case 80: case 81: case 82: case 83: case 84: case 85: case 86: case 87: case 88: case 89: case 90: case 91: case 92: case 93: case 94: case 95: case 96: case 97: case 98: case 99: case 100: case 101: case 102: case 103: case 104: case 105: case 106: case 107: case 108: case 109: case 110: case 111: case 112: case 113: case 114: case 115: case 116: case 117: case 118: case 119: case 120: case 121: case 122: case 123: case 124: case 125: case 126: case 127: case 128: case 129: case 130: case 131: case 132: case 133: case 134: case 135: case 136: case 137: case 138: case 139: case 140: case 141: case 143: case 144: case 145: case 147: case 148: case 149: case 150: case 151: case 153: case 155: case 156: case 160: case 165: case 179: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 122: return true; break;
            case 123: return true; break;
            case 124: return true; break;
            case 125: switch(pc1zBin) {case 19: case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27: case 28: case 29: case 30: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 184: case 186: case 187: case 191: case 192: case 193: case 194: case 195: case 196: case 197: case 198: case 199: case 200: case 201: case 202: case 203: case 204: case 205: case 206: case 207: case 208: case 209: case 210: case 211: case 212: case 213: case 214: case 215: case 216: case 217: case 218: case 219: case 220: case 221: case 222: case 223: case 224: case 225: case 226: case 227: case 228: case 229: case 230: case 231: case 232: case 233: case 234: case 235: case 236: case 237: case 238: case 239: case 240: case 241: case 242: case 243: case 244: case 245: case 246: case 247: case 248: case 249: case 250: case 251: case 252: case 253: case 254: case 255: case 256: case 257: case 258: case 259: case 260: case 261: case 262: case 263: case 264: case 265: case 266: case 267: case 268: case 269: case 270: case 271: case 272: case 273: case 274: case 275: case 276: case 277: case 278: case 279: case 280: case 281: case 282: case 283: case 284: case 285: case 286: case 287: case 288: case 289: case 290: case 291: case 292: case 293: case 294: case 295: case 296: case 297: case 298: case 299: case 300: case 301: case 302: case 303: case 304: case 305: case 306: case 307: case 308: case 309: case 310: case 311: case 312: case 313: case 314: case 315: case 316: case 317: case 318: case 319: case 320: case 321: case 322: case 323: case 324: case 325: case 326: case 327: case 328: case 329: case 330: case 331: case 332: case 333: case 334: case 335: case 336: case 337: case 338: case 339: case 340: case 341: case 342: case 343: case 344: case 345: case 346: case 347: return true; break;} break;
            case 126: switch(pc1zBin) {case 19: case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27: case 28: case 29: case 30: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 189: case 190: case 191: case 192: case 193: case 194: case 195: case 196: case 197: case 198: case 199: case 200: case 201: case 202: case 203: case 204: case 205: case 206: case 207: case 208: case 209: case 210: case 211: case 212: case 213: case 214: case 215: case 216: case 217: case 218: case 219: case 220: case 221: case 222: case 223: case 224: case 225: case 226: case 227: case 228: case 229: case 230: case 231: case 232: case 233: case 234: case 235: case 236: case 237: case 238: case 239: case 240: case 241: case 242: case 243: case 244: case 245: case 246: case 247: case 248: case 249: case 250: case 251: case 252: case 253: case 254: case 255: case 256: case 257: case 258: case 259: case 260: case 261: case 262: case 263: case 264: case 265: case 266: case 267: case 268: case 269: case 270: case 271: case 272: case 273: case 274: case 275: case 276: case 277: case 278: case 279: case 280: case 281: case 282: case 283: case 284: case 285: case 286: case 287: case 288: case 289: case 290: case 291: case 292: case 293: case 294: case 295: case 296: case 297: case 298: case 299: case 300: case 301: case 302: case 303: case 304: case 305: case 306: case 307: case 308: case 309: case 310: case 311: case 312: case 313: case 314: case 315: case 316: case 317: case 318: case 319: case 320: case 321: case 322: case 323: case 324: case 325: case 326: case 327: case 328: case 329: case 330: case 331: case 332: case 333: case 334: case 335: case 336: case 337: case 338: case 339: case 340: case 341: case 342: case 343: case 344: case 345: case 346: case 347: return true; break;} break;
            case 127: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 189: case 190: case 191: case 192: case 193: case 194: case 195: case 196: case 197: case 198: case 199: case 200: case 201: case 202: case 203: case 204: case 205: case 206: case 207: case 208: case 209: case 210: case 211: case 212: case 213: case 214: case 215: case 216: case 217: case 218: case 219: case 220: case 221: case 222: case 223: case 224: case 225: case 226: case 227: case 228: case 229: case 230: case 231: case 232: case 233: case 234: case 235: case 236: case 237: case 238: case 239: case 240: case 241: case 242: case 243: case 244: case 245: case 246: case 247: case 248: case 249: case 250: case 251: case 252: case 253: case 254: case 255: case 256: case 257: case 258: case 259: case 260: case 261: case 262: case 263: case 264: case 265: case 266: case 267: case 268: case 269: case 270: case 271: case 272: case 273: case 274: case 275: case 276: case 277: case 278: case 279: case 280: case 281: case 282: case 283: case 284: case 285: case 286: case 287: case 288: case 289: case 290: case 291: case 292: case 293: case 294: case 295: case 296: case 297: case 298: case 299: case 300: case 301: case 302: case 303: case 304: case 305: case 306: case 307: case 308: case 309: case 310: case 311: case 312: case 313: case 314: case 315: case 316: case 317: case 318: case 319: case 320: case 321: case 322: case 323: case 324: case 325: case 326: case 327: case 328: case 329: case 330: case 331: case 332: case 333: case 334: case 335: case 336: case 337: case 338: case 339: case 340: case 341: case 342: case 343: case 344: case 345: case 346: case 347: return true; break;} break;
            case 128: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 189: case 190: case 191: case 192: case 193: case 194: case 195: case 196: case 197: case 198: case 199: case 200: case 201: case 202: case 203: case 204: case 205: case 206: case 207: case 208: case 209: case 210: case 211: case 212: case 213: case 214: case 215: case 216: case 217: case 218: case 219: case 220: case 221: case 222: case 223: case 224: case 225: case 226: case 227: case 228: case 229: case 230: case 231: case 232: case 233: case 234: case 235: case 236: case 237: case 238: case 239: case 240: case 241: case 242: case 243: case 244: case 245: case 246: case 247: case 248: case 249: case 250: case 251: case 252: case 253: case 254: case 255: case 256: case 257: case 258: case 259: case 260: case 261: case 262: case 263: case 264: case 265: case 266: case 267: case 268: case 269: case 270: case 271: case 272: case 273: case 274: case 275: case 276: case 277: case 278: case 279: case 280: case 281: case 282: case 283: case 284: case 285: case 286: case 287: case 288: case 289: case 290: case 291: case 292: case 293: case 294: case 295: case 296: case 297: case 298: case 299: case 300: case 301: case 302: case 303: case 304: case 305: case 306: case 307: case 308: case 309: case 310: case 311: case 312: case 313: case 314: case 315: case 316: case 317: case 318: case 319: case 320: case 321: case 322: case 323: case 324: case 325: case 326: case 327: case 328: case 329: case 330: case 331: case 332: case 333: case 334: case 335: case 336: case 337: case 338: case 339: case 340: case 341: case 342: case 343: case 344: case 345: case 346: case 347: return true; break;} break;
            case 129: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 189: case 190: case 191: case 192: case 193: case 194: case 195: case 196: case 197: case 198: case 199: case 200: case 201: case 202: case 203: case 204: case 205: case 206: case 207: case 208: case 209: case 210: case 211: case 212: case 213: case 214: case 215: case 216: case 217: case 218: case 219: case 220: case 221: case 222: case 223: case 224: case 225: case 226: case 227: case 228: case 229: case 230: case 231: case 232: case 233: case 234: case 235: case 236: case 237: case 238: case 239: case 240: case 241: case 242: case 243: case 244: case 245: case 246: case 247: case 248: case 249: case 250: case 251: case 252: case 253: case 254: case 255: case 256: case 257: case 258: case 259: case 260: case 261: case 262: case 263: case 264: case 265: case 266: case 267: case 268: case 269: case 270: case 271: case 272: case 273: case 274: case 275: case 276: case 277: case 278: case 279: case 280: case 281: case 282: case 283: case 284: case 285: case 286: case 287: case 288: case 289: case 290: case 291: case 292: case 293: case 294: case 295: case 296: case 297: case 298: case 299: case 300: case 301: case 302: case 303: case 304: case 305: case 306: case 307: case 308: case 309: case 310: case 311: case 312: case 313: case 314: case 315: case 316: case 317: case 318: case 319: case 320: case 321: case 322: case 323: case 324: case 325: case 326: case 327: case 328: case 329: case 330: case 331: case 332: case 333: case 334: case 335: case 336: case 337: case 338: case 339: case 340: case 341: case 342: case 343: case 344: case 345: case 346: case 347: return true; break;} break;
            case 130: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 208: case 209: case 210: case 211: case 212: case 213: case 214: case 215: case 216: return true; break;} break;
            case 131: switch(pc1zBin) {case 11: case 12: case 13: case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27: case 28: case 29: case 30: case 31: case 32: case 33: case 34: case 35: case 36: case 37: case 38: case 39: case 40: case 41: case 42: case 43: case 44: case 45: case 46: case 47: case 48: case 49: case 50: case 51: case 52: case 53: case 54: case 55: case 56: case 57: case 58: case 59: case 60: case 61: case 62: case 63: case 64: case 65: case 66: case 67: case 68: case 69: case 70: case 71: case 72: case 73: case 74: case 75: case 76: case 77: case 78: case 79: case 80: case 81: case 82: case 83: case 84: case 85: case 86: case 87: case 88: case 89: case 90: case 91: case 92: case 93: case 94: case 95: case 96: case 97: case 98: case 99: case 100: case 101: case 102: case 103: case 104: case 105: case 106: case 107: case 108: case 109: case 110: case 111: case 112: case 113: case 114: case 115: case 116: case 117: case 118: case 119: case 120: case 121: case 122: case 123: case 124: case 125: case 126: case 127: case 128: case 129: case 130: case 131: case 132: case 133: case 134: case 135: case 136: case 137: case 138: case 139: case 140: case 141: case 142: case 143: case 144: case 145: case 146: case 147: case 148: case 149: case 150: case 151: case 152: case 153: case 154: case 155: case 156: case 157: case 158: case 159: case 160: case 161: case 162: case 163: case 164: case 165: case 166: case 167: case 168: case 169: case 170: case 171: case 172: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 132: switch(pc1zBin) {case 11: case 12: case 13: case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27: case 28: case 29: case 30: case 31: case 32: case 33: case 34: case 35: case 36: case 37: case 38: case 39: case 40: case 41: case 42: case 43: case 44: case 45: case 46: case 47: case 48: case 49: case 50: case 51: case 52: case 53: case 54: case 55: case 56: case 57: case 58: case 59: case 60: case 61: case 62: case 63: case 64: case 65: case 66: case 67: case 68: case 69: case 70: case 71: case 72: case 73: case 74: case 75: case 76: case 77: case 78: case 79: case 80: case 81: case 82: case 83: case 84: case 85: case 86: case 87: case 88: case 89: case 90: case 91: case 92: case 93: case 94: case 95: case 96: case 97: case 98: case 99: case 100: case 101: case 102: case 103: case 104: case 105: case 106: case 107: case 108: case 109: case 110: case 111: case 112: case 113: case 114: case 115: case 116: case 117: case 118: case 119: case 120: case 121: case 122: case 123: case 124: case 125: case 126: case 127: case 128: case 129: case 130: case 131: case 132: case 133: case 134: case 135: case 136: case 137: case 138: case 139: case 140: case 141: case 142: case 143: case 144: case 145: case 146: case 147: case 148: case 149: case 150: case 151: case 152: case 153: case 154: case 155: case 156: case 157: case 158: case 159: case 160: case 161: case 162: case 163: case 164: case 165: case 166: case 167: case 168: case 169: case 170: case 171: case 172: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 133: switch(pc1zBin) {case 11: case 12: case 13: case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27: case 28: case 29: case 30: case 31: case 32: case 33: case 34: case 35: case 36: case 37: case 38: case 39: case 40: case 41: case 42: case 43: case 44: case 45: case 46: case 47: case 48: case 49: case 50: case 51: case 52: case 53: case 54: case 55: case 56: case 57: case 58: case 59: case 60: case 61: case 62: case 63: case 64: case 65: case 66: case 67: case 68: case 69: case 70: case 71: case 72: case 73: case 74: case 75: case 76: case 77: case 78: case 79: case 80: case 81: case 82: case 83: case 84: case 85: case 86: case 87: case 88: case 89: case 90: case 91: case 92: case 93: case 94: case 95: case 96: case 97: case 98: case 99: case 100: case 101: case 102: case 103: case 104: case 105: case 106: case 107: case 108: case 109: case 110: case 111: case 112: case 113: case 114: case 115: case 116: case 117: case 118: case 119: case 120: case 121: case 122: case 123: case 124: case 125: case 126: case 127: case 128: case 129: case 130: case 131: case 132: case 133: case 134: case 135: case 136: case 137: case 138: case 139: case 140: case 141: case 142: case 143: case 144: case 145: case 146: case 147: case 148: case 149: case 150: case 151: case 152: case 153: case 154: case 155: case 156: case 157: case 158: case 159: case 160: case 161: case 162: case 163: case 164: case 165: case 166: case 167: case 168: case 169: case 170: case 171: case 172: case 173: case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 134: switch(pc1zBin) {case 12: case 13: case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27: case 28: case 29: case 30: case 31: case 32: case 33: case 34: case 35: case 36: case 37: case 38: case 39: case 40: case 41: case 42: case 43: case 44: case 45: case 46: case 47: case 48: case 49: case 50: case 51: case 52: case 53: case 54: case 55: case 56: case 57: case 58: case 59: case 60: case 61: case 62: case 63: case 64: case 65: case 66: case 67: case 68: case 69: case 70: case 71: case 72: case 73: case 74: case 75: case 76: case 77: case 78: case 79: case 80: case 81: case 82: case 83: case 84: case 85: case 86: case 87: case 88: case 89: case 90: case 91: case 92: case 93: case 94: case 95: case 96: case 97: case 98: case 99: case 100: case 101: case 102: case 103: case 104: case 105: case 106: case 107: case 108: case 109: case 110: case 111: case 112: case 113: case 114: case 115: case 116: case 117: case 118: case 119: case 120: case 121: case 122: case 123: case 125: case 126: case 127: case 131: case 133: case 134: case 135: case 136: case 137: case 139: case 147: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 135: switch(pc1zBin) {case 15: case 16: case 17: case 25: case 26: case 27: case 28: case 29: case 30: case 31: case 32: case 33: case 34: case 35: case 36: case 37: case 38: case 39: case 40: case 41: case 42: case 43: case 44: case 45: case 46: case 47: case 48: case 49: case 50: case 51: case 52: case 53: case 54: case 55: case 56: case 57: case 58: case 59: case 60: case 61: case 62: case 63: case 64: case 65: case 66: case 67: case 68: case 69: case 70: case 71: case 72: case 73: case 74: case 75: case 76: case 77: case 78: case 81: case 84: case 85: case 87: case 154: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 136: switch(pc1zBin) {case 32: case 33: case 34: case 35: case 36: case 37: case 38: case 39: case 40: case 41: case 42: case 43: case 44: case 45: case 46: case 47: case 48: case 49: case 50: case 52: case 54: case 56: case 57: case 58: case 59: case 61: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 137: switch(pc1zBin) {case 174: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 138: switch(pc1zBin) {case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 139: return true; break;
            case 140: return true; break;
            case 141: return true; break;
            case 142: return true; break;
            case 143: return true; break;
            case 144: return true; break;
            case 145: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 189: case 190: case 191: case 192: case 193: case 194: case 195: case 196: case 197: case 198: case 199: case 200: case 201: return true; break;} break;
            case 146: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 189: case 190: case 191: case 192: case 193: case 194: case 195: case 196: case 197: case 198: case 199: case 200: case 201: case 303: case 304: case 305: case 306: case 307: case 308: case 309: case 310: return true; break;} break;
            case 147: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 184: case 185: case 186: case 187: case 188: case 189: case 190: case 191: case 192: case 193: case 194: case 195: case 196: case 197: case 198: case 199: case 200: case 201: case 303: case 304: case 305: case 306: case 307: case 308: case 309: case 310: return true; break;} break;
            case 148: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 200: case 201: case 303: case 304: case 305: case 306: case 307: case 308: case 309: case 310: return true; break;} break;
            case 149: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 199: case 200: case 201: case 303: case 304: case 305: case 306: case 307: case 308: case 309: case 310: return true; break;} break;
            case 150: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: case 187: case 188: case 189: case 190: case 191: case 192: case 193: case 194: case 195: case 196: case 197: case 198: case 199: case 200: case 201: return true; break;} break;
            case 151: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 152: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
            case 153: switch(pc1zBin) {case 174: case 175: case 176: case 177: case 178: case 179: case 180: case 181: case 182: case 183: case 184: case 185: case 186: return true; break;} break;
         }
      }
      return false;
   }

   bool IsDeadPC2(const double pc2z, const double pc2phi) // Run7AuAu200
   {
      return false;
   }

   bool IsDeadPC3(const double phi, const double pc3z, const double pc3phi) // Run7AuAu200
   {
      return false;
   }

   bool IsDeadEMCal(const double phi, const double zed, const int sect, 
                    const double pemcy, const double pemcz) // Run7AuAu200
   {
        return false;
   }

   bool IsBadSlat(const int slat) // Run7AuAu200
   {
      return false; 
   }

   bool IsBadStripTOFw(const int strip) // Run7AuAu200
   {
      return false;
   }

   bool IsDeadTOFe(const double zed, const double tofy, const double tofz) // Run7AuAu200
   {
      return false;
   }

   bool IsDeadTOFw(const double zed, const double board, const double alpha) // Run7AuAu200
   {
        return false;
   }
}

#endif /* DEAD_AREAS_CUTS_H */
