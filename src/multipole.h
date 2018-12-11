#ifndef _multipole_h
#define _multipole_h

#include <iostream>
#include <math.h>
using std::istringstream; // part of iostream
#include <stdio.h>
#include <fstream>
#include <cassert>
#include <iomanip>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <unistd.h> //for getcwd
using std::string;
#include "params.h"
#include "vector.h"
/* 
   A class for multipole expansions.  Currently treats up to hexadecapoles 
   (rank 4).  Use spherical tensor moments throughout.

   GJB 09/09
*/

class Multipole {
  int Rank;
  int Nmom; // length of the multipole moment vector of rank Rank.
  Vector Moments; // Vector of Nmom 1, 4, 9, etc.

 public:
  // Constructors & Destructor
  Multipole(); // blank constructor
  Multipole(int rank, double* moments); 
  Multipole(const Multipole& other, bool copy_data=true);
  ~Multipole();
  
  void Initialize(int rank);
  void Initialize(const Multipole& other, bool copy_data=true);

  void Set() {Moments.Set();};
 
  // For accessing elements of the multipole vector matrix
  double& operator()(int i) const
    {
      if (i >= Moments.GetLength()) {
	printf("ERROR: Multipole Moment vector element %d is out of bounds.\n",i);
	exit(1);
      }
      return Moments[i];
    } 


  Multipole& Scale(double factor);
  Multipole operator+(const Multipole& other);
  Multipole& operator=(const Multipole& other);

  // Some basic functions
  int GetRank() const {return Rank;};
  int GetLength() const {return Nmom;};
  Vector GetMoments() const {return Moments;};
  Vector GetMomentsVector() {return Moments;}; // JDH get rid of this?

  // Print out the multipoles
  void Print(string title);
  void VectorPrint(string title) {Moments.Print(title);};

  // JDH NMR code
  Vector Spherical_to_Cartesian(Multipole Spherical);

};

#endif
