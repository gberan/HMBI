#include <string>
#include <algorithm>
using std::string;
#include "multipole.h"


Multipole::Multipole(): Rank(-1), Nmom(0) {

}

// primary constructor
Multipole::Multipole(int rank, double* moments) : Rank(rank), Nmom(0) {

  // Figure out how many spherical tensor moments this rank corresponds to.
  Nmom = 0;
  for (int i=0;i<=rank;i++) 
    Nmom += 2*i + 1;

  // Note, dipole moments are read as 10 (z), 11c (x), 11s (y).  We will store them
  // as x,y,z, or 11c, 11s, 10, so we need to reorder them.
  double z = moments[1];
  double x = moments[2];
  double y = moments[3];
  moments[1] = x;
  moments[2] = y;
  moments[3] = z;

  // Initialize the moments vector
  Vector tmp(moments, Nmom);
  

  Moments.Initialize(Nmom);
  Moments = tmp;

}

// Copy constructor
Multipole::Multipole(const Multipole& other, bool copy_data) : Rank(other.Rank), Nmom(other.Nmom) {

  // Read the ranks
  //Rank = other.GetRank();
  //Nmom = other.GetLength();

  // Initialize the vectors
  Moments.Initialize(Nmom);

  if (copy_data) 
    Moments = other.GetMoments();
  else 
    Moments.Set();

}

Multipole::~Multipole() {

}

void Multipole::Initialize(int rank) {

  Rank = rank;
  // Figure out how many spherical tensor moments this rank corresponds to.
  Nmom = 0;
  for (int i=0;i<=rank;i++) 
    Nmom += 2*i + 1;
  
  Moments.Initialize(Nmom);

}

void Multipole::Initialize(const Multipole& other, bool copy_data) {

   // Read the ranks
  Rank = other.GetRank();
  Nmom = other.GetLength();

  // Initialize the vectors
  Moments.Initialize(Nmom);

  if (copy_data) 
    Moments = other.Moments;
  else 
    Moments.Set();

}


// Multiply the multipole moments by a double
Multipole& Multipole::Scale(double factor) {
  Moments.Scale(factor);
  return *this;
}

// Add two multipole vectors of equivalent rank or different rank.
Multipole Multipole::operator+(const Multipole& other) {

  // Create an empty multipole object with the proper rank.  If the
  // rank of the input multipoles differ, the final multipole adopts
  // the larger of the two ranks.
  Multipole result;
  int max_rank;

  (Rank >= other.Rank) ? max_rank = Rank : max_rank = other.Rank;
  //if (Rank >= other.Rank) max_rank = Rank;
  //else max_rank = other.Rank;
  result.Initialize(max_rank);

  result.Moments =  Moments.Add(other.Moments,true);
  //result.Moments += other.Moments;

  /*
  if (other.Rank == Rank) {
    result.Moments = Moments;
    result.Moments += other.Moments;
  }
  else { // ranks differ, just add element-by-element.

    for (int i=0;i<Nmom;i++)
      result.GetMoments()(i) = Moments.GetMoments()(i);
    for (i=0;i<other.Nmom;i++)
      result.GetMoments()(i) += other.Moments.GetMoments()(i);

  }
    */

  return result;

}

// Add another multipole to the current one.
Multipole& Multipole::operator+=(const Multipole& other) {

  // First test that the ranks of both are equal. of other is less than or equal to rank of this.


  if (other.Rank != Rank) {
    printf("ERROR: Multipole += operator requires both multipoles have the same rank.\n");
    printf("       The multipole + operator can handle mixed-rank addition.\n");
    exit(1);
  }


  Moments += other.Moments;

  return *this;

}


// Overload the '=' operator
Multipole& Multipole::operator=(const Multipole& other) {
  if (this!=&other) {
    // Read the ranks
    Rank = other.Rank;
    Nmom = other.Nmom;
    
    // Initialize the moments
    if (Nmom > 0) {
      Moments = other.Moments;
    }
  }
  return *this;
}

// Print out the multipole moments with pretty formatting
void Multipole::Print(string title) {  

  if (Rank < 0) {
    printf("Error: Multipole::Print(): Unknown Rank = %d\n",Rank); 
    exit(1);
  }

  printf("%s:  Rank %d\n",title.c_str(),Rank);
  printf("  Charge:       (00)  %10.6f\n",Moments[0]);
  if (Rank >=1) {
    printf("  Dipole:       (10)  %10.6f  (11c) %10.6f  (11s) %10.6f\n",
	   Moments[3],Moments[1],Moments[2]);
  }
  if (Rank >=2) {
    printf("  Quadrupole:   (20)  %10.6f  (21c) %10.6f  (21s) %10.6f  (22c) %10.6f  (22s) %10.6f\n",
	   Moments[4],Moments[5],Moments[6],Moments[7],Moments[8]);
  }
  if (Rank >=3) {
    printf("  Octopole:     (30)  %10.6f  (31c) %10.6f  (31s) %10.6f  (32c) %10.6f  (32s) %10.6f\n",
	   Moments[9],Moments[10],Moments[11],Moments[12],Moments[13]);
    printf("                (33c) %10.6f  (33s) %10.6f\n",
	   Moments[14],Moments[15]);
  }
  if (Rank >=4) {
    printf("  Hexadecapole: (40)  %10.6f  (41c) %10.6f  (41s) %10.6f  (42c) %10.6f  (42s) %10.6f\n",
	   Moments[16],Moments[17],Moments[18],Moments[19],Moments[20]);
    printf("                (43c) %10.6f  (43s) %10.6f  (44c) %10.6f  (44s) %10.6f\n",
	   Moments[21],Moments[22],Moments[23],Moments[24]);
  }

  printf("\n");
}

// debugging routine: use to zero out all moments except dipole & quadrupole.  Helpful with
// testing some induction routines, but serves no general purpose.
void Multipole::MaskElements() {

  printf("Masking multipole moments\n");
  Moments[0] = 0.0;
  if (Nmom > 9) {
    for (int i=9;i<Nmom;i++)
      Moments[i]=0.0;
  }

  //Print("Masked multipole moments");

}
