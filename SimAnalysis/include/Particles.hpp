// $HEADER$
//------------------------------------------------------------------------------------------------
//                                 Particles global containers
//------------------------------------------------------------------------------------------------
// Particles
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic global containers for storing constant particle properties
 **/
//------------------------------------------------------------------------------------------------

#ifndef PARTICLES_HPP
#define PARTICLES_HPP

#include <string>
#include <map>
#include <vector>
#include <array>

struct 
{
   const int pion = 211;
   const int kaon = 321;
   const int proton = 2212;
   const int electron = 11;
   const int noPID = 1999;
   const int junk = -9999;
} PartId;

struct
{
   std::map<std::string, int> iterMap = 
   {
      {"electron", 0},
      {"positron", 1},
      {"pion", 2},
      {"apion", 3},
      {"kaon", 4},
      {"akaon", 5},
      {"proton", 6},
      {"aproton", 7}
   };

   std::vector<int> id = {11, -11, 
                          211, -211, 
                          321, -321, 
                          2212, -2212};

   std::vector<std::string> name = {"electron", "positron", 
                                    "pion", "apion", 
                                    "kaon", "akaon", 
                                    "proton", "aproton"};
   //short name
   std::vector<std::string> sname = {"e", "e", 
                                     "pi", "pi", 
                                     "K", "K", 
                                     "P", "P"};
   
   //for embedding
   std::vector<std::string> abs_name = {"electron", "electron", 
                                        "pion", "pion", 
                                        "kaon", "kaon", 
                                        "proton", "proton"};
   
   std::vector<double> mass = {0.5110999e-3, 0.5110999e-3, 
                               0.139570, 0.139570, 
                               0.493677, 0.493677, 
                               0.938272, 0.938272};
   
   std::vector<int> geant_id = {1, 2, 
                                8, 9, 
                                11, 12, 
                                14, 15};
   
   std::vector<int> charge = {-1, 1, 
                              1, -1, 
                              1, -1, 
                              1, -1};
} ParticleProperties;

struct Lambda1520
{
   const std::string name = "lambda1520";
   const float mass = 1.5195;
   const float gamma = 0.0165;
};

struct KStar892
{
   const std::string name = "KStar892";
   const float mass = 0.89166;
   const float gamma = 0.0508;
};

#endif
