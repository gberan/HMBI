#ifndef _constants_h
#define _constants_h

/*
  A namespace for physical constants and conversion factors.  To access
these, put this in your C file:

     #include "constants.h"
     using namespace hmbi_constants;

Then the variables will be accessible.


  GJB 11/09
 */

namespace hmbi_constants {

  const double pi = 3.14159265359;
  const double ec = 1.60217733e-19; // e- charge, C per au
  const double epsilon = 8.85418781762e-12; // permittivity,  C^2/m  
  const double Na = 6.0221367e23; // Avogadro's number
  const double h_planck = 6.62617363e-34; // Planck's constant, J*s
  const double c_light = 2.99792458e8; // Speed of light, m/s


  // Conversion factors
  const double HartreesToKJpermole = 2625.500;
  const double HartreesToKcalpermole = 627.509;
  const double BohrToAng = 0.5291772085;
  const double AngToBohr = 1.0/BohrToAng;
  const double MtoAng = 1.0e10; // meters -> Angstroms
  const double DegreesToRadians = pi/180.0;
  const double RadiansToDegrees = 180.0/pi;
  
}


#endif
