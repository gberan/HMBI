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
  const double k_boltz = 1.3806488e-23; //Boltzmann constant, J/K
  const double R_gas_constant = 8.3144621; // in J/(mol*K)

  // Conversion factors
  const double HartreesToKJpermole = 2625.500;
  const double HartreesToKcalpermole = 627.509;
  const double HartreesToJ = HartreesToKJpermole/Na*1000;
  const double HartreesToEV = 27.212;
  const double EVToHartrees = 1.0/HartreesToEV;
  const double BohrToAng = 0.5291772085;
  const double AngToBohr = 1.0/BohrToAng;
  const double MtoAng = 1.0e10; // meters -> Angstroms
  const double DegreesToRadians = pi/180.0;
  const double RadiansToDegrees = 180.0/pi;
  const double WaveNumberToFreq = 100*c_light;//cm^-1 to Hertz
  const double FreqtoWaveNumber = 1/WaveNumberToFreq;


  
  const double GigaPascalsToHartreesperBohrcube = 1/HartreesToKJpermole*Na/AngToBohr/AngToBohr/AngToBohr/1.0e24;
}


#endif
