/** 
 *  @file   Constants.hpp
 *  @brief Contains constants used for code in PairAnalysisPhenix
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysisPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <unordered_map>
#include <string>

/// speed of light in vacuum [cm/ns]
const double SPEED_OF_LIGHT = 29.9792458;
/// mass of a charged pion [GeV/c^2]
const double MASS_PION = 139.57039e-3;

/// @namespace ParticleMap
/// maps id of particles to their properties (mass, GEANT id, name, etc.)
namespace ParticleMap
{
   /// massses [GeV/c]
   std::unordered_map<int, double> mass = 
   {
      {11, 0.510998e-3},
      {211, 139.570e-3},
      {321, 493.677e-3},
      {2212, 938.272e-3},
      {-11, 0.510998e-3},
      {-211, 139.570e-3},
      {-321, 493.677e-3},
      {-2212, 938.272e-3}
   };
   /// names
   std::unordered_map<int, std::string> name =
   {
      {11, "e"},
      {211, "pi+"},
      {321, "k+"},
      {2212, "p"},
      {-11, "ebar"},
      {-211, "pi-"},
      {-321, "k-"},
      {-2212, "pbar"}
   };
   /// ids by GEANT numbering scheme (https://www.star.bnl.gov/public/comp/simu/newsite/gstar/Manual/particle_id.html)
   std::unordered_map<int, int> idGEANT =
   {
      {11, 3},
      {211, 8},
      {321, 11},
      {2212, 14},
      {-11, 2},
      {-211, 9},
      {-321, 12},
      {-2212, 15}
   };
   /* template to add some property later
   std::unordered_map<int, >  =
   {
      {11, },
      {211, },
      {321, },
      {2212, },
      {-11, },
      {-211, },
      {-321, },
      {-2212, }
   };
   */
}


#endif /* CONSTANTS_HPP */
