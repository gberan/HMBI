#ifndef _polarizability_h
#define _polarizability_h

#include <iostream>
#include <math.h>
using std::istringstream; // part of iostream
#include <stdio.h>
#include <fstream.h>
#include <cassert>
#include <iomanip.h>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <unistd.h> //for getcwd
using std::string;
#include "params.h"
#include "vector.h"
#include "matrix.h"
/* 
   A class for storing polarizabilities.  Currently treats up to quadrupole 
   polarizabilities (rank 2).  Use spherical tensor moments throughout.

   GJB 09/09
*/

class Polarizability {
  int Rank;
  int Npol; // dimension of the polarizability matrix
  Matrix Pols; // size (Npol x Npol)



 public:
  // Constructors & Destructor
  Polarizability(); // blank constructor
  Polarizability(int rank, double* pols); 
  Polarizability(const Polarizability& other, bool copy_data=true);
  ~Polarizability();
  
  Polarizability& operator=(const Polarizability& other);

  // Some basic functions
  int GetRank() const {return Rank;};
  int GetLength() const {return Npol;};
  Matrix GetPolarizabilities() const {return Pols;};


  // For accessing elements of the polarizability matrix
  double& operator()(int i, int j) const
    {
      if (i >= Pols.GetRows() || j >= Pols.GetCols()) {
	printf("ERROR: Polarizability matrix element (%d,%d) is out of bounds.\n",i,j);
	exit(1);
      }
      return Pols(i,j);
    } 

  // Print out the multipoles
  void Print(string title);
  void PrintAll(string title);

};

#endif
